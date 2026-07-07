"""
nexrad_to_zarr.py
=================
Convert a list of NEXRAD Level II volume scans into a single, consolidated,
time-appended Zarr store using xradar.

The problem this solves
-----------------------
NEXRAD WSR-88D radars run adaptive volume-coverage patterns (VCPs). Within a
single VCP the *number and position* of elevation sweeps changes volume to
volume because of:

* **AVSET** - the radar truncates high tilts when there is no echo aloft, so a
  quiet volume may have 14 sweeps and an active storm volume 20.
* **SAILS/MRLE** - supplemental low-level base tilts are inserted at variable
  index positions mid-volume.
* **Split cuts** - the lowest tilts are scanned twice, once in a long-range
  "surveillance" mode carrying the polarimetric moments (DBZH, ZDR, PHIDP,
  RHOHV) and once in a Doppler mode carrying velocity/spectrum-width
  (VRADH, WRADH).

Because of this, sweep *index* (sweep_0, sweep_1, ...) does not correspond to a
fixed elevation angle across volumes, so naive concatenation is wrong. This
module normalizes every volume onto a fixed **canonical** set of elevation
angles keyed by fixed angle, merges the split cuts into one sweep per angle,
drops SAILS duplicates, and fills AVSET-truncated tilts with NaN. The result is
a rectangular hypercube per sweep that concatenates cleanly along a new
``volume_time`` dimension.

Store layout
------------
``<out>.zarr``
  ``/sweep_0`` ... ``/sweep_14``    (one group per canonical elevation)
      dims:  (volume_time, azimuth, range)
      vars:  DBZH VRADH WRADH ZDR PHIDP RHOHV CCORH  (float32)
      coords: volume_time, azimuth, range, sweep_fixed_angle,
              latitude, longitude, altitude

Open the whole store as a DataTree:
    import datatree ; dt = datatree.open_datatree("kiwa.zarr", engine="zarr")
or a single sweep:
    import xarray as xr ; ds = xr.open_zarr("kiwa.zarr", group="sweep_0")
"""
import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import xradar
import xradar.util as xu

# --- Canonical VCP-215 geometry (KIWA) -------------------------------------
# Fixed elevation angles (deg) after merging split cuts and SAILS duplicates.
CANONICAL_ANGLES = [0.5, 0.9, 1.3, 1.8, 2.4, 3.1, 4.0, 5.1,
                    6.4, 8.0, 10.0, 12.0, 14.0, 16.7, 19.5]
# The 7 dual-pol / Doppler moments we retain.
MOMENTS = ["DBZH", "VRADH", "WRADH", "ZDR", "PHIDP", "RHOHV", "CCORH"]
# Target (n_azimuth, n_range) per canonical angle, taken from a full 20-sweep
# volume. Low tilts are 0.5-deg super-resolution (720 rays); higher tilts 1 deg.
TARGET_GEOM = {0.5: (720, 1832), 0.9: (720, 1832), 1.3: (720, 1712),
               1.8: (360, 1536), 2.4: (360, 1336), 3.1: (360, 1160),
               4.0: (360, 984), 5.1: (360, 820), 6.4: (360, 676),
               8.0: (360, 540), 10.0: (360, 452), 12.0: (360, 384),
               14.0: (360, 328), 16.7: (360, 276), 19.5: (360, 240)}
RANGE_START, RANGE_RES = 125.0, 250.0   # gate-center start (m) & spacing (m)


def _nearest_canon(a):
    return min(CANONICAL_ANGLES, key=lambda c: abs(c - a))


def _fixed_azimuth(naz):
    """Canonical ray-center azimuths for a sweep with `naz` rays."""
    res = 360.0 / naz
    return (np.arange(naz) * res + res / 2.0).astype("float64")


def _fixed_range(nrng):
    return (RANGE_START + RANGE_RES * np.arange(nrng)).astype("float64")


def normalize_sweep(ds, naz, nrng):
    """Snap one source sweep onto the fixed (azimuth, range) grid.

    We do NOT rely on ``reindex_angle`` to emit exactly ``naz`` rays - a sweep
    whose azimuth coverage is slightly short (e.g. 357 or 717 rays) would then
    fail to align on append. Instead we sort by azimuth, make it a real index,
    and ``reindex`` onto the canonical azimuth grid, guaranteeing exactly naz
    rays. Range is handled the same way.
    """
    # sort rays by azimuth and index them
    if "azimuth" in ds.coords:
        ds = ds.sortby("azimuth")
        if "azimuth" not in ds.indexes:
            ds = ds.set_xindex("azimuth")
    tgt_az = _fixed_azimuth(naz)
    ds = ds.reindex(azimuth=tgt_az, method="nearest",
                    tolerance=(360.0 / naz))
    if "range" not in ds.indexes and "range" in ds.coords:
        ds = ds.set_xindex("range")
    ds = ds.reindex(range=_fixed_range(nrng), method="nearest", tolerance=1.0)
    return ds


# pyart field name -> our moment name
_PYART_FIELD_MAP = {
    "reflectivity": "DBZH", "velocity": "VRADH", "spectrum_width": "WRADH",
    "differential_reflectivity": "ZDR", "differential_phase": "PHIDP",
    "cross_correlation_ratio": "RHOHV", "clutter_filter_power_removed": "CCORH"}


def _iter_source_sweeps(reader_obj):
    """Yield (fixed_angle, sweep_Dataset) from either an xradar DataTree or a
    Py-ART Radar. Each sweep_Dataset has dims (azimuth, range), our moment
    variable names, and azimuth/range coordinates."""
    if hasattr(reader_obj, "groups"):                      # xradar DataTree
        for g in reader_obj.groups:
            if not g.startswith("/sweep"):
                continue
            sds = reader_obj[g].ds
            yield float(sds["sweep_fixed_angle"].values), sds
    else:                                                  # Py-ART Radar
        radar = reader_obj
        rng = radar.range["data"].astype("float64")
        for s in range(radar.nsweeps):
            sl = radar.get_slice(s)
            az = radar.azimuth["data"][sl].astype("float64")
            ang = float(radar.fixed_angle["data"][s])
            dvars = {}
            for pf, mom in _PYART_FIELD_MAP.items():
                if pf in radar.fields:
                    arr = np.ma.filled(
                        radar.fields[pf]["data"][sl].astype("float32"), np.nan)
                    dvars[mom] = (("azimuth", "range"), arr)
            yield ang, xr.Dataset(dvars, coords={"azimuth": az, "range": rng})


def build_canonical_volume(reader_obj):
    """dict: canonical angle -> merged Dataset of the retained moments.

    Accepts an xradar DataTree or a Py-ART Radar (fallback reader)."""
    groups = {}
    for ang, sds in _iter_source_sweeps(reader_obj):
        groups.setdefault(_nearest_canon(ang), []).append(sds)
    out = {}
    for canon, sweeps in groups.items():
        naz, nrng = TARGET_GEOM[canon]
        # surveillance sweep (most range gates) first -> its DBZH/dual-pol win
        sweeps = sorted(sweeps, key=lambda s: -s.sizes["range"])
        merged = None
        for sds in sweeps:
            nds = normalize_sweep(sds, naz, nrng)
            part = xr.Dataset({m: nds[m] for m in MOMENTS if m in nds.data_vars})
            if merged is None:
                merged = part
            else:
                for m in part.data_vars:          # add Doppler moments only
                    if m not in merged.data_vars:
                        merged[m] = part[m]
        out[canon] = merged
    return out


def volume_to_sweep_datasets(dt, volume_time):
    """Return {sweep_i: Dataset(volume_time, azimuth, range)} for one volume,
    with all 15 canonical sweeps present (AVSET-truncated ones all-NaN)."""
    vol = build_canonical_volume(dt)
    lat, lon, alt = _site_geometry(dt)
    out = {}
    for i, ang in enumerate(CANONICAL_ANGLES):
        naz, nrng = TARGET_GEOM[ang]
        if ang in vol:
            ds = vol[ang][[m for m in MOMENTS if m in vol[ang].data_vars]]
            for m in MOMENTS:                      # backfill any missing moment
                if m not in ds.data_vars:
                    ds[m] = (("azimuth", "range"),
                             np.full((naz, nrng), np.nan, "float32"))
        else:                                       # AVSET-truncated tilt
            ds = xr.Dataset(
                {m: (("azimuth", "range"),
                     np.full((naz, nrng), np.nan, "float32")) for m in MOMENTS},
                coords={"azimuth": _fixed_azimuth(naz),
                        "range": _fixed_range(nrng)})
        ds = ds[MOMENTS]                            # fixed variable order
        for m in MOMENTS:
            ds[m] = ds[m].astype("float32")
        ds = ds.drop_vars([c for c in ("elevation", "time")
                           if c in ds.coords], errors="ignore")
        ds = ds.expand_dims(volume_time=[np.datetime64(volume_time)])
        ds = ds.assign_coords(sweep_fixed_angle=float(ang),
                              latitude=lat, longitude=lon, altitude=alt)
        ds.attrs = {}
        out[f"sweep_{i}"] = ds
    return out


def _site_geometry(reader_obj):
    """(lat, lon, alt) from an xradar DataTree or Py-ART Radar."""
    if hasattr(reader_obj, "groups"):
        root = reader_obj["/"].ds
        return (float(root.latitude), float(root.longitude),
                float(root.altitude))
    r = reader_obj
    return (float(r.latitude["data"][0]), float(r.longitude["data"][0]),
            float(r.altitude["data"][0]))


def _read_volume(local_path):
    """Open a Level II volume, preferring xradar; fall back to Py-ART for the
    small fraction of volumes whose split-cut layout xradar cannot reconcile.
    Returns (reader_obj, reader_name)."""
    try:
        return xradar.io.open_nexradlevel2_datatree(local_path), "xradar"
    except Exception:
        import pyart
        return pyart.io.read_nexrad_archive(local_path), "pyart"


def build_zarr_store(scan_list, out_zarr, s3_client=None, bucket="unidata-nexrad-level2",
                     scratch="scratch", download_workers=8, overwrite=True,
                     progress=True):
    """Convert a list of NEXRAD Level II scans into one time-appended Zarr store.

    Parameters
    ----------
    scan_list : list of dict
        Each item must have ``key`` (S3 key) and ``timestamp`` (pandas Timestamp).
        Use ``list_scan_keys(radar, start, end)`` to build this.
    out_zarr : str
        Path of the output ``.zarr`` store (created/overwritten).
    s3_client : boto3 client, optional
        UNSIGNED client for the Unidata bucket. Defaults to ``nexrad_s3_client``
        from the nexrad-site-rainfall skill if available.
    download_workers : int
        Threads used to prefetch the next scans while the current one is written.

    Returns
    -------
    dict  summary: {n_scans, n_written, volume_times, out_zarr, failures}
    """
    from concurrent.futures import ThreadPoolExecutor
    if s3_client is None:
        s3_client = nexrad_s3_client()          # noqa: F821  (skill helper)
    os.makedirs(scratch, exist_ok=True)
    if overwrite and os.path.exists(out_zarr):
        shutil.rmtree(out_zarr)

    scans = sorted(scan_list, key=lambda d: d["timestamp"])

    def _fetch(item):
        lp = os.path.join(scratch, os.path.basename(item["key"]))
        if not os.path.exists(lp):
            s3_client.download_file(bucket, item["key"], lp)
        return item, lp

    written, failures, vtimes = 0, [], []
    readers = {}
    first = True
    # prefetch pool: download ahead while we read/write
    with ThreadPoolExecutor(max_workers=download_workers) as pool:
        for item, lp in pool.map(_fetch, scans):
            vt = pd.Timestamp(item["timestamp"])
            try:
                dt, reader = _read_volume(lp)
                readers[reader] = readers.get(reader, 0) + 1
                sweeps = volume_to_sweep_datasets(dt, vt)
                for name, ds in sweeps.items():
                    if first:
                        enc = {m: {"chunks": (1,) + ds[m].shape[1:]}
                               for m in MOMENTS}
                        # Pin CF time units to a FIXED epoch. Without this,
                        # xarray re-derives "since <first value>" on every
                        # append and the appended dates drift.
                        enc["volume_time"] = {
                            "units": "seconds since 1970-01-01T00:00:00",
                            "calendar": "proleptic_gregorian",
                            "dtype": "int64"}
                        ds.to_zarr(out_zarr, group=name, mode="w", encoding=enc)
                    else:
                        ds.to_zarr(out_zarr, group=name, append_dim="volume_time")
                written += 1
                vtimes.append(vt)
                first = False
                if progress and written % 20 == 0:
                    print(f"  {written}/{len(scans)} volumes written")
            except Exception as e:                # keep going on a bad volume
                failures.append({"key": item["key"], "error": repr(e)})
                if progress:
                    print(f"  SKIP {item['key']}: {e!r}")
            finally:
                if os.path.exists(lp):
                    os.remove(lp)

    # consolidate metadata for fast opening
    try:
        import zarr
        zarr.consolidate_metadata(out_zarr)
    except Exception as e:
        print("consolidate_metadata warning:", repr(e))

    return {"n_scans": len(scans), "n_written": written,
            "volume_times": vtimes, "out_zarr": out_zarr,
            "readers": readers, "failures": failures}


# --- Read-back / QC / plotting helpers -------------------------------------
def open_sweep(store, sweep=0):
    """Open one canonical sweep group as an xarray Dataset."""
    return xr.open_zarr(store, group=f"sweep_{sweep}")


def qc_report(store):
    """Return a DataFrame summarizing per-sweep volume coverage in the store."""
    rows = []
    for i in range(len(CANONICAL_ANGLES)):
        ds = xr.open_zarr(store, group=f"sweep_{i}")
        cov = float(np.isfinite(ds.DBZH).any(dim=("azimuth", "range")).mean())
        rows.append({"sweep": f"sweep_{i}",
                     "fixed_angle": float(ds.sweep_fixed_angle),
                     "n_vol": ds.sizes["volume_time"],
                     "naz": ds.sizes["azimuth"], "nrng": ds.sizes["range"],
                     "vol_with_DBZH_frac": round(cov, 3)})
        ds.close()
    return pd.DataFrame(rows)


def ppi_from_store(store, volume_time=None, sweep=0, moment="DBZH"):
    """Read one PPI from the Zarr store and return (X_km, Y_km, field, meta).

    X is east-west, Y is north-south, both in km from the radar. If
    ``volume_time`` is None the volume with the highest max DBZH on ``sweep`` is
    used (the storm peak)."""
    ds = xr.open_zarr(store, group=f"sweep_{sweep}")
    if volume_time is None:
        vmax = ds.DBZH.max(dim=("azimuth", "range")).values
        volume_time = ds.volume_time.values[int(np.nanargmax(vmax))]
    sel = ds.sel(volume_time=np.datetime64(volume_time), method="nearest")
    field = sel[moment].values
    az = np.deg2rad(ds.azimuth.values)
    rng = ds.range.values / 1000.0
    AZ, RNG = np.meshgrid(az, rng, indexing="ij")
    X, Y = RNG * np.sin(AZ), RNG * np.cos(AZ)
    meta = {"volume_time": str(sel.volume_time.values),
            "fixed_angle": float(ds.sweep_fixed_angle),
            "lat": float(ds.latitude), "lon": float(ds.longitude),
            "max": float(np.nanmax(field))}
    ds.close()
    return X, Y, field, meta


def plot_ppi(store, out_png, volume_time=None, sweep=0, moment="DBZH",
             vmin=-10, vmax=70, cmap="ChaseSpectral", extent_km=150, floor=5.0):
    """Render a Cartesian-km PPI from the store to ``out_png``.

    Uses a plain matplotlib axis (not cartopy GeoAxes) for robustness with the
    large super-resolution pcolormesh."""
    import matplotlib.pyplot as plt
    try:
        import cmweather  # noqa: F401  registers ChaseSpectral
    except Exception:
        cmap = "viridis"
    X, Y, field, meta = ppi_from_store(store, volume_time, sweep, moment)
    fig, ax = plt.subplots(figsize=(8, 7))
    pm = ax.pcolormesh(X, Y, np.ma.masked_less(field, floor), cmap=cmap,
                       vmin=vmin, vmax=vmax, shading="auto")
    th = np.linspace(0, 2 * np.pi, 200)
    for r in (50, 100, 150):
        if r <= extent_km:
            ax.plot(r * np.cos(th), r * np.sin(th), color="0.5", lw=0.6, ls="--")
            ax.text(0, r, f"{r} km", color="0.4", fontsize=7, ha="center",
                    va="bottom")
    ax.plot(0, 0, "k^", ms=9)
    ax.set_xlabel("East-west distance from radar (km)")
    ax.set_ylabel("North-south distance (km)")
    ax.set_aspect("equal")
    ax.set_xlim(-extent_km, extent_km)
    ax.set_ylim(-extent_km, extent_km)
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(f"{moment} (dBZ)" if moment == "DBZH" else moment)
    ax.set_title(f"{meta['fixed_angle']:.1f}deg {moment} - {meta['volume_time']}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    return meta

