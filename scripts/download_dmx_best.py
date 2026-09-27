#!/usr/bin/env python3
"""Download the small GEO demuxlet BEST members needed for donor rescue.

The complete GSE313736 RAW archive is intentionally not downloaded.  This
script uses the existing byte-range index and downloads only *_dmx.best.gz
members, which contain barcode-level donor-assignment evidence.
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
RAW = ROOT / "data" / "raw" / "GSE313736" / "demux_best"
INDEX = META / "GSE313736_remote_tar_index.json"
LOG = ROOT / "logs" / "phase0R_dmx_download.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(member: dict) -> tuple[str, str, int, str]:
    name = member["name"]
    out = RAW / name
    part = out.with_suffix(out.suffix + ".part")
    expected = int(member["size"])
    if out.exists() and out.stat().st_size == expected:
        return name, "EXISTS", expected, sha256(out)
    RAW.mkdir(parents=True, exist_ok=True)
    start = part.stat().st_size if part.exists() else 0
    if start > expected:
        part.unlink()
        start = 0
    url = json.loads(INDEX.read_text())["url"]
    chunk_size = 16 * 1024 * 1024
    try:
        with part.open("ab") as fh:
            while start < expected:
                end = min(expected - 1, start + chunk_size - 1)
                headers = {"Range": f"bytes={member['data_offset'] + start}-{member['data_offset'] + end}"}
                with requests.get(url, headers=headers, stream=True, timeout=(30, 600)) as response:
                    response.raise_for_status()
                    if response.status_code != 206:
                        raise RuntimeError(f"expected HTTP 206, got {response.status_code}")
                    received = 0
                    for payload in response.iter_content(1024 * 1024):
                        if payload:
                            fh.write(payload)
                            received += len(payload)
                if received <= 0:
                    raise RuntimeError("empty range response")
                start += received
        if part.stat().st_size != expected:
            raise RuntimeError(f"size {part.stat().st_size} != expected {expected}")
        part.replace(out)
        return name, "DOWNLOADED", expected, sha256(out)
    except Exception as exc:  # pragma: no cover - surfaced in the log
        return name, f"ERROR:{exc!r}", 0, ""


def main() -> None:
    members = [m for m in json.loads(INDEX.read_text())["members"] if m["name"].endswith("_dmx.best.gz")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(fetch, members))
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("w") as fh:
        fh.write("filename\tstatus\tbytes\tsha256\n")
        for row in results:
            fh.write("\t".join(map(str, row)) + "\n")
    failures = [row for row in results if row[1].startswith("ERROR")]
    print(f"selected_members={len(results)} completed={len(results) - len(failures)} errors={len(failures)}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
