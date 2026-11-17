"""Download the public NHANES modules used in the two-period analysis.

The files are retrieved directly from CDC/NCHS. Existing files are retained
unless ``--force`` is supplied. Participant-level files remain local and are
excluded from version control.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PERIODS = {
    "2017-March 2020": {
        "directory": ROOT / "data" / "raw_2017_2020",
        "base_url": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles",
        "modules": (
            "P_DEMO",
            "P_LUX",
            "P_DIQ",
            "P_GHB",
            "P_ALQ",
            "P_BMX",
            "P_HEPA",
            "P_HEPB_S",
            "P_HEPBD",
            "P_IMQ",
        ),
    },
    "August 2021-August 2023": {
        "directory": ROOT / "data" / "raw",
        "base_url": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles",
        "modules": (
            "DEMO_L",
            "LUX_L",
            "HEPA_L",
            "HEPB_S_L",
            "HEPBD_L",
            "HEPC_L",
            "IMQ_L",
            "DIQ_L",
            "GHB_L",
            "ALQ_L",
            "BMX_L",
            "HIQ_L",
            "HUQ_L",
            "SMQ_L",
            "BIOPRO_L",
            "ALB_CR_L",
            "BPXO_L",
            "BPQ_L",
            "HDL_L",
        ),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_period(period: str, specification: dict[str, object], force: bool) -> None:
    destination_dir = specification["directory"]
    base_url = specification["base_url"]
    modules = specification["modules"]
    assert isinstance(destination_dir, Path)
    assert isinstance(base_url, str)
    assert isinstance(modules, tuple)

    destination_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for module in modules:
        destination = destination_dir / f"{module}.xpt"
        url = f"{base_url}/{module}.xpt"
        if force or not destination.exists():
            print(f"Downloading {period}: {module} ...")
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "nhanes-temporal-reproducibility-code/2.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                destination.write_bytes(response.read())
        manifest.append(
            {
                "period": period,
                "module": module,
                "file": destination.name,
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
                "source_url": url,
            }
        )

    manifest_path = destination_dir / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Verified {len(manifest)} modules for {period}: {manifest_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="overwrite existing files")
    args = parser.parse_args()

    for period, specification in PERIODS.items():
        download_period(period, specification, args.force)


if __name__ == "__main__":
    main()
