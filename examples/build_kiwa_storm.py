#!/usr/bin/env python
"""Build the KIWA 2025-11-18 storm-window Zarr store and a peak-volume PPI.

This is the end-to-end example the repository README describes. It streams
NEXRAD Level II volumes from the public ``unidata-nexrad-level2`` S3 bucket,
normalizes the adaptive VCP geometry onto canonical elevation angles, and
writes a single consolidated Zarr store.

Requires the ``list_scan_keys`` helper from the nexrad-site-rainfall workflow
(or supply your own list of ``{"key", "timestamp"}`` dicts).
"""
import nexrad_to_zarr as n2z

RADAR = "KIWA"
START = "2025-11-18 00:00"
END   = "2025-11-19 00:00"
OUT   = "kiwa_20251118.zarr"


def main():
    # Enumerate the Level II volume keys for the window.
    # from your scan-listing helper:
    #   keys = list_scan_keys(RADAR, START, END)
    # keys must be a list of {"key": <s3 key>, "timestamp": <pandas Timestamp>}
    raise SystemExit(
        "Provide `keys` via your scan-listing helper, then call "
        "n2z.build_zarr_store(keys, OUT, download_workers=10)."
    )


if __name__ == "__main__":
    main()
