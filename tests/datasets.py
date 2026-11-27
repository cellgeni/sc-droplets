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
            Remote(BASE + "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices_h5.h5", "4972237f2e8ee11efea5c37850a2e26b0103f317bd656cbe81bfc9c6eb041225"),
            Remote(BASE + "cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", "097ff873bae645bec5586898ac002ccebbce952d52ac19ea0a35c27c4820df26"),
        ),
        Dataset(
            "pbmc_1k_v2_3_0",
            "3.0.0",
            "SC3Pv2",
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_raw_feature_bc_matrix.h5", "fd258b397faf0e38eed00f8f06f6068cd380963d15e90e25a67c5cd9e19e4f71"),
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5", "344673d4eca5af4b4f8343ef85c039fa5b24fea75a7d714a54427548b7b60213"),
        ),
        Dataset(
            "pbmc_1k_v3_3_0",
            "3.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_raw_feature_bc_matrix.h5", "ed17bec5896a8208acfd6202e035fc36c65da14b3605aa577e93496b5d1ebd81"),
            Remote(BASE + "cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5", "8191f576550c1b449d03441b9eb098ee9f73fa82513d171ba87d31d551e3ffda"),
        ),
        Dataset(
            "pbmc5k_nextgem_8_0",
            "8.0.0",
            "SC3Pv3",
            Remote(
                BASE + "cell-exp/8.0.0/5k_Human_PBMC_TotalSeqB_3p_nextgem_Multiplex/"
                "5k_Human_PBMC_TotalSeqB_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
                "1861c0b6f433a980b39a18f0b6fe0df0ad16e57f94d255d51bf2cee749175b76",
            ),
            Remote(
                BASE + "cell-exp/8.0.0/5k_Human_PBMC_TotalSeqB_3p_nextgem_5k_Human_PBMC_TotalSeqB_3p_nextgem/"
                "5k_Human_PBMC_TotalSeqB_3p_nextgem_5k_Human_PBMC_TotalSeqB_3p_nextgem_count_sample_filtered_feature_bc_matrix.h5",
                "20eec804fdc58e247c773da4a777d3ff147292612e0cadcf1656fcd0ecffcf3f",
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
                "a93a27fa3c48ff5011b625702fcba95943d80d053d19e3ac46f9842fb7de488d",
            ),
            Remote(
                BASE + "cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x/"
                "5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5",
                "6c46fc939ed9d39d8cc228c29f2e6e6b0af0771a1a2c2bf1abc71e6788b122d5",
            ),
        ),
        Dataset(
            "spinal_cord_10_0",
            "10.0.0",
            "SC3Pv4",
            Remote(
                BASE + "cell-exp/10.0.0/WT_mouse_spinal_cord_P112_specimen_1_Multiplex/"
                "WT_mouse_spinal_cord_P112_specimen_1_Multiplex_raw_feature_bc_matrix.h5",
                "70c934c2b25fb7352b24894829b00bd8af8aefec93d0338b6761131f39175bf0",
            ),
            Remote(
                BASE + "cell-exp/10.0.0/WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1/"
                "WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5",
                "9eb052d76f5dc167023e509977f4d22d0e81b5e0b4f42bee71045b99ac9e17c9",
            ),
            note="nuclei",
        ),
        Dataset(
            "flex_lymphnode1_7_2",
            "7.2.0",
            "MFRP",
            Remote(
                BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode1_BC1_AB1/"
                "4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode1_BC1_AB1_count_sample_raw_feature_bc_matrix.h5",
                "ecc062ac14fdb050a08446c9173a243f272ecb6bcc8df73d25c20796e27eb9c3",
            ),
            Remote(
                BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode1_BC1_AB1/"
                "4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode1_BC1_AB1_count_sample_filtered_feature_bc_matrix.h5",
                "d7f93c736efce7b625bea14bd808148836cd76d25d1afe9425c4ee549364aa85",
            ),
            note="multiplexed Flex, one sample's own raw matrix; Gene Expression + Antibody Capture",
        ),
        Dataset(
            "flex_lymphnode2_7_2",
            "7.2.0",
            "MFRP",
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode2_BC2_AB2/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode2_BC2_AB2_count_sample_raw_feature_bc_matrix.h5", "fcf97ab27c1797830f90f533d4348ebdd139244dea254e2be2c1c25ef332a84a"),
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode2_BC2_AB2/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_LymphNode2_BC2_AB2_count_sample_filtered_feature_bc_matrix.h5", "dd920e8d33a955a9c9cb854c99b9ac7071cfd8f0df3b9c67fe76c8939cd1f970"),
            note="the same 4-plex Flex run as flex_lymphnode1_7_2",
        ),
        Dataset(
            "flex_spleen1_7_2",
            "7.2.0",
            "MFRP",
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen1_BC3_AB3/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen1_BC3_AB3_count_sample_raw_feature_bc_matrix.h5", "ba0b1e56c551a43b21878ee00986f8306faccae3d1c4f85995688a20cbb5c6eb"),
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen1_BC3_AB3/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen1_BC3_AB3_count_sample_filtered_feature_bc_matrix.h5", "3437ffb1069cba94e9d228b9ddb5de2d71dd149321fe13bd1d9e8412d810defb"),
            note="the same 4-plex Flex run as flex_lymphnode1_7_2",
        ),
        Dataset(
            "flex_spleen2_7_2",
            "7.2.0",
            "MFRP",
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen2_BC4_AB4/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen2_BC4_AB4_count_sample_raw_feature_bc_matrix.h5", "8484a7af9cd76586f9421e265faab4f3a7c0d28fd3a041ef08cc3d5137c205b1"),
            Remote(BASE + "cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen2_BC4_AB4/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Spleen2_BC4_AB4_count_sample_filtered_feature_bc_matrix.h5", "e633bce55426bca7f37df2f2d7f4ffb2a14258686f6ad9100954bcb3d8f4cc1a"),
            note="the same 4-plex Flex run as flex_lymphnode1_7_2",
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
