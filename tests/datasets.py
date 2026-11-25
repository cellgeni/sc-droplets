"""Public 10x datasets the integration tests check Cell Ranger and STARsolo against.

Each dataset is a raw matrix and the cells a known Cell Ranger release called
on it. Files are downloaded once into `$SC_DROPLETS_DATA` (default
`~/.cache/sc-droplets/datasets`) and checked against their SHA-256. They are
10x Genomics' public datasets (CC BY 4.0); nothing is redistributed.

    uv run python tests/datasets.py           # download every dataset
    uv run python tests/datasets.py pbmc4k_2_1   # one
"""

from __future__ import annotations

import gzip
import hashlib
import io
import os
import shutil
import sys
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

BASE = "https://cf.10xgenomics.com/samples/"


@dataclass(frozen=True)
class Remote:
    url: str
    sha256: str


@dataclass(frozen=True)
class Dataset:
    name: str
    cellranger: str  # the release that made the filtered matrix
    chemistry: str
    raw: Remote
    filtered: Remote
    note: str = ""


DATASETS: Dict[str, Dataset] = {
    d.name: d
    for d in [
        Dataset(
            "pbmc4k_2_1",
            "2.1.0",
            "SC3Pv2",
            Remote(BASE + "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices_h5.h5", ""),
            Remote(BASE + "cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", ""),
        ),
        Dataset(
            "pbmc_1k_v2_3_0",
            "3.0.0",
            "SC3Pv2",
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_raw_feature_bc_matrix.h5", ""),
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5", ""),
        ),
        Dataset(
            "pbmc_1k_v3_3_0",
            "3.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_raw_feature_bc_matrix.h5", ""),
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5", ""),
        ),
        Dataset(
            "pbmc5k_nextgem_8_0",
            "8.0.0",
            "SC3Pv3",
            Remote(
                BASE + "cell-exp/8.0.0/5k_Human_PBMC_TotalSeqB_3p_nextgem_Multiplex/"
                "5k_Human_PBMC_TotalSeqB_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
                "",
            ),
            Remote(
                BASE + "cell-exp/8.0.0/5k_Human_PBMC_TotalSeqB_3p_nextgem_5k_Human_PBMC_TotalSeqB_3p_nextgem/"
                "5k_Human_PBMC_TotalSeqB_3p_nextgem_5k_Human_PBMC_TotalSeqB_3p_nextgem_count_sample_filtered_feature_bc_matrix.h5",
                "",
            ),
            note="Gene Expression + Antibody Capture",
        ),
        Dataset(
            "pbmc5k_gemx_9_0",
            "9.0.0",
            "SC3Pv4",
            Remote(
                BASE + "cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_Multiplex/"
                "5k_Human_Donor1_PBMC_3p_gem-x_Multiplex_count_raw_feature_bc_matrix.h5",
                "",
            ),
            Remote(
                BASE + "cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x/"
                "5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5",
                "",
            ),
        ),
        Dataset(
            "spinal_cord_10_0",
            "10.0.0",
            "SC3Pv4",
            Remote(
                BASE + "cell-exp/10.0.0/WT_mouse_spinal_cord_P112_specimen_1_Multiplex/"
                "WT_mouse_spinal_cord_P112_specimen_1_Multiplex_raw_feature_bc_matrix.h5",
                "",
            ),
            Remote(
                BASE + "cell-exp/10.0.0/WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1/"
                "WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5",
                "",
            ),
            note="nuclei",
        ),
    ]
}


def cache_dir() -> Path:
    return Path(os.environ.get("SC_DROPLETS_DATA", Path.home() / ".cache" / "sc-droplets" / "datasets"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(remote: Remote) -> Path:
    """The local copy of a remote file, downloaded and checked if needed."""
    target = cache_dir() / remote.url[len(BASE) :]
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + ".part")
        # cf.10xgenomics.com refuses urllib's default agent.
        request = urllib.request.Request(remote.url, headers={"User-Agent": "sc-droplets-tests"})
        with urllib.request.urlopen(request) as response, open(partial, "wb") as out:
            shutil.copyfileobj(response, out)
        partial.rename(target)
    if remote.sha256 and _sha256(target) != remote.sha256:
        raise RuntimeError(f"{target} does not match its recorded SHA-256; delete it to download again")
    return target


def raw_path(name: str) -> Path:
    return fetch(DATASETS[name].raw)


def filtered_barcodes(name: str, genome: Optional[str] = None) -> List[str]:
    """The barcodes of the filtered matrix (one genome's, for a Cell Ranger 2 MEX)."""
    path = fetch(DATASETS[name].filtered)
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path) as archive:
            members = [m for m in archive.getmembers() if m.name.endswith("barcodes.tsv")]
            if genome is not None:
                members = [m for m in members if f"/{genome}/" in m.name]
            if len(members) != 1:
                raise ValueError(f"{path}: expected one barcodes.tsv, found {[m.name for m in members]}")
            text = archive.extractfile(members[0]).read().decode()
        return text.split()
    import h5py

    with h5py.File(path, "r") as f:
        group = f["matrix"] if "matrix" in f else f[genome or next(iter(f))]
        return [b.decode() for b in group["barcodes"][()]]


def main(argv: List[str]) -> int:
    names = argv or list(DATASETS)
    for name in names:
        dataset = DATASETS[name]
        for remote in (dataset.raw, dataset.filtered):
            path = fetch(remote)
            print(f"{name}: {path} {path.stat().st_size / 2**20:.1f} MiB sha256={_sha256(path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
