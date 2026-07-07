"""
arm_cfradial_to_zarr.py
=======================
Convert a day of ARM CfRadial-1 radar volumes into a single, consolidated,
time-appended Zarr store using xradar.

This is the ARM sibling of ``nexrad_to_zarr.py``. Both modules solve the same
core problem -- turn a directory/list of single-volume radar files into one
Zarr store with a ``volume_time`` dimension you can slice, subset, and plot
from -- but the *normalization* challenge is different for each radar.

The problem this solves
-----------------------
A single ARM radar day is not necessarily a single scan strategy. The busiest
CACTI C-SAPR2 day (2018-11-13, datastream ``corcsapr2cfrppiqcM1.b1``,
Cordoba, Argentina) **interleaves two incompatible strategies** in one stream:

* **PPI surveillance** - 15 sweeps, fixed tilts 0.5 -> 32.7 deg, full 360 deg
  azimuth, 1100 gates to ~110 km. One volume every ~15 min, all day.
* **Sector rapid-scan** - 1 sweep at 3.6 deg, a *moving* azimuth wedge, 275
  gates to ~27 km. Confined to 13:14 -> 15:22 UTC (this is what makes the day
  ~2.4x the normal volume count).

The two strategies have different sweep counts, range depths, and azimuth
coverage, so they cannot share arrays. Rather than force one grid, this module
writes **two parallel group trees** in one store -- ``ppi/`` and ``sector/`` --
each internally rectangular and cleanly appendable along ``volume_time``. A
homogeneous PPI-only day simply collapses to just the ``ppi/`` tree.

Where NEXRAD's problem was adaptive VCP sweep counts (AVSET / SAILS / split
cuts) normalized onto canonical *elevation angles* in a single tree, ARM's
problem here is *distinct scan strategies* separated into distinct *trees*. The
shared machinery -- reindex each sweep onto a fixed azimuth grid, pin the
``volume_time`` CF encoding so appends don't drift the dates, stream one volume
at a time, consolidate, and read/QC/plot straight from the store -- is the same.

Store layout
------------
``<out>.zarr``
  ``/ppi/sweep_0`` ... ``/ppi/sweep_14``    (one group per PPI elevation)
      dims:  (volume_time, azimuth, range)   e.g. (85, 360, 1100)
  ``/sector/sweep_0``                        (the single sector tilt)
      dims:  (volume_time, azimuth, range)   e.g. (148, 360, 275)
  vars (both):  the retained QC polarimetric moments (float32) + int masks
  coords:       volume_time, azimuth (fixed 360-bin 1deg grid), range,
                sweep_fixed_angle, latitude, longitude, altitude

Open one sweep:
    import xarray as xr
    ds = xr.open_zarr("cacti.zarr", group="ppi/sweep_0")
"""
import os
import re
import glob
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import xradar
import xradar.util as xu

# --- Retained fields -------------------------------------------------------
# The 12 QC moments kept from the ~35 in the b1 product. Uncorrected / lag-1 /
# V-channel duplicates are dropped to hold store size down.
FIELDS = [
    "attenuation_corrected_reflectivity_h",
    "attenuation_corrected_differential_reflectivity",
    "reflectivity",
    "differential_reflectivity",
    "copol_correlation_coeff",
    "specific_differential_phase",
    "differential_phase",
    "mean_doppler_velocity",
    "spectral_width",
    "signal_to_noise_ratio_copolar_h",
    "censor_mask",
    "classification_mask",
]

# Fixed azimuth grid: 360 bins, ray centers 0.5 .. 359.5 deg.
AZ_NBIN = 360
# Sweep count that identifies a full PPI surveillance volume; anything else is
# treated as a sector rapid-scan (the CACTI day has 15-sweep PPI + 1-sweep
# sector, but the split is by "is it the multi-tilt survey" not a magic number).
PPI_NSWEEP = 15

# CfRadial-1 encodes the volume time in the filename: ....YYYYMMDD.HHMMSS.nc
_TIME_RE = re.compile(r"\.(\d{8})\.(\d{6})\.nc$")


def volume_time_from_name(fname):
    """Parse the volume start time from an ARM CfRadial-1 filename."""
    m = _TIME_RE.search(os.path.basename(fname))
    if not m:
        raise ValueError(f"no YYYYMMDD.HHMMSS timestamp in {fname!r}")
    return pd.to_datetime(m.group(1) + m.group(2), format="%Y%m%d%H%M%S")


def _fixed_azimuth(nbin=AZ_NBIN):
    """Canonical ray-center azimuths (deg) for an `nbin`-ray sweep."""
    res = 360.0 / nbin
    return (np.arange(nbin) * res + res / 2.0).astype("float64")


def normalize_sweep(sw, root, volume_time, fields=FIELDS, nbin=AZ_NBIN):
    """Snap one source sweep onto the fixed azimuth grid and attach coords.

    Two source quirks are handled here (the same two that bite any CfRadial ->
    Zarr append):

    1. **Duplicate azimuths.** Raw sweeps oscillate between 360 and 361 rays
       with a fractional start angle; the wrapping ray duplicates an azimuth.
       A duplicated index makes ``reindex`` raise, so we de-duplicate (keep
       first) before snapping onto the fixed grid.
    2. **Fixed grid, not just reindex_angle.** We reindex onto exactly `nbin`
       canonical azimuths (nearest, tolerance half a bin) so every volume --
       including the moving sector wedge -- aligns exactly on append. Bins the
       sweep did not illuminate become NaN.

    ``latitude``/``longitude``/``altitude`` live at the DataTree *root* in
    CfRadial-1, not on the sweep, so they are pulled from ``root``.
    """
    keep = [f for f in fields if f in sw.data_vars]
    ds = sw[keep]
    if "azimuth" in ds.coords:
        ds = ds.sortby("azimuth")
        _, uidx = np.unique(ds["azimuth"].values, return_index=True)
        if len(uidx) != ds.sizes["azimuth"]:
            ds = ds.isel(azimuth=np.sort(uidx))
        if "azimuth" not in ds.indexes:
            ds = ds.set_xindex("azimuth")
    ds = ds.reindex(azimuth=_fixed_azimuth(nbin), method="nearest",
                    tolerance=(360.0 / nbin) / 2.0)
    # drop per-ray coords that would otherwise vary volume to volume
    ds = ds.drop_vars([c for c in ("time", "elevation")
                       if c in ds.coords], errors="ignore")
    ds = ds.assign_coords(
        sweep_fixed_angle=float(sw["sweep_fixed_angle"].values),
        latitude=float(root["latitude"].values),
        longitude=float(root["longitude"].values),
        altitude=float(root["altitude"].values),
        volume_time=volume_time,
    ).expand_dims("volume_time")
    ds.attrs = {}
    return ds


def classify_file(path):
    """Return ('ppi'|'sector', n_sweep) for one CfRadial-1 volume.

    Uses the ``sweep`` dimension size, read cheaply without a full xradar open.
    """
    import netCDF4
    with netCDF4.Dataset(path) as d:
        n = d.dimensions["sweep"].size
    return ("ppi" if n == PPI_NSWEEP else "sector"), n


def survey_scan_strategy(files):
    """Classify every file in `files`; return a DataFrame + a split dict.

    Columns: file, volume_time, tree, n_sweep. The returned ``(df, groups)``
    has ``groups = {"ppi": [paths...], "sector": [paths...]}`` sorted by time.
    """
    rows, groups = [], {"ppi": [], "sector": []}
    for f in files:
        tree, n = classify_file(f)
        vt = volume_time_from_name(f)
        rows.append({"file": os.path.basename(f), "volume_time": vt,
                     "tree": tree, "n_sweep": n})
        groups[tree].append(f)
    for k in groups:
        groups[k] = sorted(groups[k], key=volume_time_from_name)
    df = pd.DataFrame(rows).sort_values("volume_time").reset_index(drop=True)
    return df, groups


# CF time units pinned to a FIXED epoch. Without this, xarray re-derives
# "seconds since <first value>" on every append and the appended dates drift
# (a 15-minute-later volume can be written 15 *days* later).
_TIME_ENC = {"units": "seconds since 1970-01-01T00:00:00",
             "calendar": "proleptic_gregorian", "dtype": "int64"}


def _encoding_for(ds):
    enc = {v: {"chunks": (1,) + ds[v].shape[1:]} for v in ds.data_vars}
    enc["volume_time"] = dict(_TIME_ENC)
    return enc


def build_zarr_store(files, out_zarr, fields=FIELDS, overwrite=True,
                     progress=True):
    """Convert a directory or list of ARM CfRadial-1 volumes into one Zarr store.

    Parameters
    ----------
    files : str or list of str
        A glob directory (``"raw_20181113"`` -> ``raw_20181113/*.nc``) or an
        explicit list of file paths.
    out_zarr : str
        Output ``.zarr`` path (created / overwritten).
    fields : list of str
        Source variable names to retain (default: the 12 QC moments).

    Returns
    -------
    dict  summary: {n_files, ppi_vols, sector_vols, trees, out_zarr}
    """
    if isinstance(files, str):
        files = sorted(glob.glob(os.path.join(files, "*.nc")))
    if overwrite and os.path.exists(out_zarr):
        shutil.rmtree(out_zarr)

    _, groups = survey_scan_strategy(files)
    written = {"ppi": 0, "sector": 0}
    first_write = True

    for tree in ("ppi", "sector"):
        flist = groups[tree]
        for vi, f in enumerate(flist):
            vt = volume_time_from_name(f)
            dt = xradar.io.open_cfradial1_datatree(f)
            root = dt.ds
            nsw = sum(1 for g in dt.groups if g.startswith("/sweep_"))
            for si in range(nsw):
                ds = normalize_sweep(dt[f"sweep_{si}"].ds, root, vt, fields)
                gpath = f"{tree}/sweep_{si}"
                if vi == 0:
                    mode = "w" if first_write else "a"
                    ds.to_zarr(out_zarr, group=gpath, mode=mode,
                               encoding=_encoding_for(ds), consolidated=False)
                    first_write = False
                else:
                    ds.to_zarr(out_zarr, group=gpath,
                               append_dim="volume_time", consolidated=False)
            dt.close()
            written[tree] += 1
            if progress and written[tree] % 20 == 0:
                print(f"  {tree}: {written[tree]}/{len(flist)} volumes")
        if progress:
            print(f"{tree} done: {written[tree]} volumes")

    try:
        import zarr
        for tree, ntree in _tree_sweeps(out_zarr).items():
            for si in range(ntree):
                zarr.consolidate_metadata(f"{out_zarr}/{tree}/sweep_{si}")
    except Exception as e:  # noqa: BLE001
        print("consolidate_metadata warning:", repr(e))

    return {"n_files": len(files), "ppi_vols": written["ppi"],
            "sector_vols": written["sector"],
            "trees": _tree_sweeps(out_zarr), "out_zarr": out_zarr}


# --- Read-back / QC / plotting helpers -------------------------------------
def _tree_sweeps(store):
    """{tree: n_sweep} present in the store (e.g. {'ppi': 15, 'sector': 1})."""
    out = {}
    for tree in ("ppi", "sector"):
        base = os.path.join(store, tree)
        if not os.path.isdir(base):
            continue
        n = 0
        while os.path.isdir(os.path.join(base, f"sweep_{n}")):
            n += 1
        if n:
            out[tree] = n
    return out


def open_sweep(store, tree="ppi", sweep=0):
    """Open one sweep group (e.g. ``ppi/sweep_0``) as an xarray Dataset."""
    return xr.open_zarr(store, group=f"{tree}/sweep_{sweep}", consolidated=False)


def qc_report(store, field="attenuation_corrected_reflectivity_h"):
    """Per-(tree, sweep) coverage table read back from the store."""
    rows = []
    for tree, nsw in _tree_sweeps(store).items():
        for si in range(nsw):
            ds = open_sweep(store, tree, si)
            fld = field if field in ds.data_vars else "reflectivity"
            cov = float(np.isfinite(ds[fld]).any(dim=("azimuth", "range")).mean())
            rows.append({"tree": tree, "sweep": f"sweep_{si}",
                         "fixed_angle": round(float(ds.sweep_fixed_angle), 2),
                         "n_vol": ds.sizes["volume_time"],
                         "naz": ds.sizes["azimuth"], "nrng": ds.sizes["range"],
                         "vol_with_echo_frac": round(cov, 3)})
            ds.close()
    return pd.DataFrame(rows)


def ppi_from_store(store, tree="ppi", sweep=0, volume_time=None,
                   field="attenuation_corrected_reflectivity_h"):
    """Read one PPI from the store; return (X_km, Y_km, field, meta).

    X is east-west, Y is north-south, both km from the radar. ``volume_time``
    None picks the volume with the largest 99th-percentile of ``field`` on that
    sweep (the storm peak)."""
    ds = open_sweep(store, tree, sweep)
    if field not in ds.data_vars:
        field = "reflectivity"
    if volume_time is None:
        score = np.array([np.nanpercentile(ds[field].isel(volume_time=i).values, 99)
                          for i in range(ds.sizes["volume_time"])])
        volume_time = ds.volume_time.values[int(np.nanargmax(score))]
    sel = ds.sel(volume_time=np.datetime64(volume_time), method="nearest")
    arr = sel[field].values
    el = np.deg2rad(float(ds.sweep_fixed_angle))
    az = np.deg2rad(ds.azimuth.values)
    rng = ds.range.values / 1000.0
    AZ, RNG = np.meshgrid(az, rng, indexing="ij")
    X = RNG * np.cos(el) * np.sin(AZ)
    Y = RNG * np.cos(el) * np.cos(AZ)
    meta = {"volume_time": str(sel.volume_time.values)[:19],
            "fixed_angle": float(ds.sweep_fixed_angle),
            "lat": float(ds.latitude), "lon": float(ds.longitude),
            "max": float(np.nanmax(arr))}
    ds.close()
    return X, Y, arr, meta


def plot_ppi(store, out_png, tree="ppi", sweep=0, volume_time=None,
             field="attenuation_corrected_reflectivity_h",
             vmin=-8, vmax=64, cmap="ChaseSpectral", extent_km=110, floor=None):
    """Render a Cartesian-km PPI straight from the store to ``out_png``."""
    import matplotlib.pyplot as plt
    try:
        import cmweather  # noqa: F401  registers ChaseSpectral
    except Exception:  # noqa: BLE001
        cmap = "viridis"
    X, Y, arr, meta = ppi_from_store(store, tree, sweep, volume_time, field)
    if floor is not None:
        arr = np.ma.masked_less(arr, floor)
    fig, ax = plt.subplots(figsize=(8, 7))
    pm = ax.pcolormesh(X, Y, arr, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto")
    th = np.linspace(0, 2 * np.pi, 200)
    for r in (50, 100):
        if r <= extent_km:
            ax.plot(r * np.cos(th), r * np.sin(th), color="0.5", lw=0.6, ls="--")
    ax.plot(0, 0, "k^", ms=9)
    ax.set_aspect("equal")
    ax.set_xlim(-extent_km, extent_km)
    ax.set_ylim(-extent_km, extent_km)
    ax.set_xlabel("East-west distance from radar (km)")
    ax.set_ylabel("North-south distance (km)")
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(field)
    ax.set_title(f"{meta['fixed_angle']:.1f}deg {tree} - {meta['volume_time']}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    return meta
