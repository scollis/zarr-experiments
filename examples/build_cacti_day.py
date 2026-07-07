#!/usr/bin/env python
"""Build the CACTI C-SAPR2 2018-11-13 Zarr store from ARM CfRadial-1 volumes.

This is the ARM counterpart to ``build_kiwa_storm.py``. It works from a
directory of already-downloaded ``corcsapr2cfrppiqcM1.b1`` CfRadial-1 files and
writes a single consolidated Zarr store with two group trees (``ppi/`` +
``sector/``), then renders a peak-volume PPI straight from the store.

Getting the raw files
---------------------
The b1 volumes stream from the ARM Live data service. **Do not hardcode your
token** -- read it from the environment. Request a token at
https://adc.arm.gov/armlive/ and export it:

    export ARM_LIVE_USER="YourName"
    export ARM_LIVE_TOKEN="xxxxxxxxxxxxxxxx"

Then fetch a day with ``download_arm_day`` below (curl-based; the ARM Live
endpoint is reached with ``user=<USER>:<TOKEN>`` as a query parameter). The QC'd
**b1** product is the one that streams; the raw **a1** product 404s from the
live cache.
"""
import os
import re
import json
import subprocess
import arm_cfradial_to_zarr as a2z

DATASTREAM = "corcsapr2cfrppiqcM1.b1"
DATE = "2018-11-13"          # b1 available 2018-09-23 .. 2018-12-10
RAW_DIR = "raw_20181113"
OUT = "cacti_20181113.zarr"

ARM_BASE = "https://adc.arm.gov/armlive"


def _arm_credential():
    """Return the ``user=<USER>:<TOKEN>`` query value from the environment.

    Never hardcode the token. Set ARM_LIVE_USER and ARM_LIVE_TOKEN (or a single
    ARM_LIVE_CREDENTIAL of the form ``User:token``).
    """
    cred = os.environ.get("ARM_LIVE_CREDENTIAL")
    if cred:
        return cred
    user = os.environ.get("ARM_LIVE_USER")
    token = os.environ.get("ARM_LIVE_TOKEN")
    if not (user and token):
        raise SystemExit(
            "Set ARM_LIVE_USER and ARM_LIVE_TOKEN (or ARM_LIVE_CREDENTIAL) in "
            "your environment; request a token at https://adc.arm.gov/armlive/")
    return f"{user}:{token}"


def download_arm_day(datastream=DATASTREAM, date=DATE, out_dir=RAW_DIR):
    """Download one UTC day of an ARM datastream via the ARM Live service.

    Uses ``curl`` so the sandbox proxy is honored. The token is taken from the
    environment and is NEVER written to disk or logged.
    """
    os.makedirs(out_dir, exist_ok=True)
    cred = _arm_credential()
    d = date.replace("-", "")
    # 1) query the file list for the day
    q = subprocess.run(
        ["curl", "-s",
         f"{ARM_BASE}/data/query?user={cred}"
         f"&ds={datastream}&start={d}&end={d}"],
        capture_output=True, text=True, check=True)
    files = json.loads(q.stdout).get("files", [])
    print(f"{len(files)} files for {datastream} {date}")
    # 2) fetch each file (resumable: skip files already on disk)
    for i, fn in enumerate(files):
        lp = os.path.join(out_dir, fn)
        if os.path.exists(lp) and os.path.getsize(lp) > 0:
            continue
        subprocess.run(
            ["curl", "-s", "-o", lp,
             f"{ARM_BASE}/data/saveData?user={cred}&file={fn}"],
            check=True)
        if (i + 1) % 20 == 0:
            print(f"  {i + 1}/{len(files)} downloaded")
    return out_dir


def main():
    if not os.path.isdir(RAW_DIR) or not os.listdir(RAW_DIR):
        download_arm_day()

    summary = a2z.build_zarr_store(RAW_DIR, OUT)
    print(json.dumps(summary, default=str, indent=2))

    meta = a2z.plot_ppi(OUT, "cacti_ppi_peak.png", tree="ppi", sweep=0)
    print("peak PPI:", meta)

    qc = a2z.qc_report(OUT)
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
