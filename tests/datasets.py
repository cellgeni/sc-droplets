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
        # Cell Ranger 4.0 to 7.1 (no public source): what these releases called, to see which rule
        # they follow. `expect` in the note is the --expect-cells 10x's dataset page gives.
        Dataset(
            "pbmc10k_nextgem_4_0",
            "4.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_raw_feature_bc_matrix.h5", "a0a48ffdb9e8d92bc2627f857367b4b626b3a6c6d66f4766db7247e8b001758d"),
            Remote(BASE + "cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.h5", "161df2ff3b5457138d33a53450214ec4a5b8c8790b8bd2aeec58952e621f1490"),
            note="expect 10,000",
        ),
        Dataset(
            "neuron10k_nextgem_4_0",
            "4.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/4.0.0/SC3_v3_NextGem_SI_Neuron_10K/SC3_v3_NextGem_SI_Neuron_10K_raw_feature_bc_matrix.h5", "84049467d3c5605a513b8fe2558a37151d3ecc3570dc4ba058f0a7ab97f65098"),
            Remote(BASE + "cell-exp/4.0.0/SC3_v3_NextGem_SI_Neuron_10K/SC3_v3_NextGem_SI_Neuron_10K_filtered_feature_bc_matrix.h5", "67d6f4bc3f3c11c8df9104ccb7862efb858a704186997fce51db307b20c31633"),
            note="mouse; expect 10,000",
        ),
        Dataset(
            "pbmc_parent_4_0",
            "4.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_raw_feature_bc_matrix.h5", "af88e73c76c8afdd33ae6945207b726987616c7469159163c8d31c2dfd083524"),
            Remote(BASE + "cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5", "bc4b56daf63d276f2f47ffdabd746e549c8099788e4ad887e0dfe45110588763"),
            note="the whole-transcriptome parent of the targeted PBMC datasets; expect not published",
        ),
        Dataset(
            "pbmc10k_5p_v2_5_0",
            "5.0.0",
            "SC5P-R2",
            Remote(BASE + "cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_raw_feature_bc_matrix.h5", "28f774d91e35549a0e97180761fb9ddadd7f33b752940e579424566bff8bc6cc"),
            Remote(BASE + "cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5", "d82e92f73fbe01bc2b41b1aa96220617cdbad53ef5387edead67daf5f012cc96"),
            note="5′ v2 gene expression of a GEX + VDJ + antibody run; expect 10,000",
        ),
        Dataset(
            "brain_tumor_6_0",
            "6.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_raw_feature_bc_matrix.h5", "12cbb355cc299d107a94f4381a291ab8e14903d231469d61f477bd921cb14db7"),
            Remote(BASE + "cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_filtered_feature_bc_matrix.h5", "0b4aade7c11650b7f7ee50b2383dc4f6e1d985e969560275ccb4f492144652d5"),
            note="expect 2,000",
        ),
        Dataset(
            "breast_cancer_6_0",
            "6.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.0.0/Breast_Cancer_3p/Breast_Cancer_3p_raw_feature_bc_matrix.h5", "f0ee5f4ca584e44c8fd82e46b377710eb9b50ff773dc5ea69ba19a6543d1c1af"),
            Remote(BASE + "cell-exp/6.0.0/Breast_Cancer_3p/Breast_Cancer_3p_filtered_feature_bc_matrix.h5", "87edbb067a025c899b541a1d79fd77ebc43e998b395f7c29bc23c890598ffa97"),
            note="expect 7,500",
        ),
        Dataset(
            "hgmm1k_6_0",
            "6.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.0.0/1k_hgmm_3p/1k_hgmm_3p_raw_feature_bc_matrix.h5", "8586d1e98b76a95da7c2d84588f6f4793abf2a96b77f2bcbe586ad8cc7b2145d"),
            Remote(BASE + "cell-exp/6.0.0/1k_hgmm_3p/1k_hgmm_3p_filtered_feature_bc_matrix.h5", "fcdc9e3c03536d728ffb07b7b6352d073b6a49067a973c3a1a582260cc70ad62"),
            note="human and mouse; expect 1,000",
        ),
        Dataset(
            "hgmm1k_lt_6_0",
            "6.0.0",
            "SC3Pv3LT",
            Remote(BASE + "cell-exp/6.0.0/1k_hgmm_3p_LT/1k_hgmm_3p_LT_raw_feature_bc_matrix.h5", "a53d4e54725f9f41e3307bc5c6b813961e61a7e10a95697421ed38ea1bfb1fc2"),
            Remote(BASE + "cell-exp/6.0.0/1k_hgmm_3p_LT/1k_hgmm_3p_LT_filtered_feature_bc_matrix.h5", "788620555cbc5aa086eee6dc90d72702bc2f5314c5633c806d1392245f135a2a"),
            note="human and mouse, 3′ LT; expect 1,000",
        ),
        Dataset(
            "hgmm10k_controller_6_1",
            "6.1.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.1.0/10k_hgmm_3p_nextgem_Chromium_Controller/10k_hgmm_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.h5", "3a4f400bd4b5d80319bc034054e0be79d48ff6a0db8d92a021c81b6bd576f737"),
            Remote(BASE + "cell-exp/6.1.0/10k_hgmm_3p_nextgem_Chromium_Controller/10k_hgmm_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5", "6e92c3bebc45e8fa43b9b293f9ba2a117467c608bd569a928ed493efaa109f35"),
            note="human and mouse; expect 10,000",
        ),
        Dataset(
            "pbmc500_lt_6_1",
            "6.1.0",
            "SC3Pv3LT",
            Remote(BASE + "cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_Controller/500_PBMC_3p_LT_Chromium_Controller_raw_feature_bc_matrix.h5", "73d043ce7e552dcc9168d1d2493b5bb3e18c1dca3c02fcd9e315557894e923f8"),
            Remote(BASE + "cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_Controller/500_PBMC_3p_LT_Chromium_Controller_filtered_feature_bc_matrix.h5", "1a0ccf0b34299aa4dfd025a05c8e0a92f9fd16443c5fbbf893467631b0a62995"),
            note="3′ LT; expect 500",
        ),
        Dataset(
            "pbmc10k_controller_6_1",
            "6.1.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.h5", "3941e9e2ce765590ef13765df29b33d98c612866c51c15d36ccdd6eec78e0fc0"),
            Remote(BASE + "cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5", "f4077cec270dd11e857bc388f9e5e77d4824c43aca70461222f5960593d8f09d"),
            note="expect 10,000",
        ),
        Dataset(
            "pbmc10k_chromium_x_6_1",
            "6.1.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5", "c2572841abc0db0a59ea086b038d5bd9dc272d9c1b43dbcd78d1811e7e1c3b71"),
            Remote(BASE + "cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", "42a789ff3f92f0e117d1ff4ffcb2fa8edee87dd2a03d83a81502709344a38bf0"),
            note="expect 10,000",
        ),
        Dataset(
            "pbmc10k_chromium_x_6_1_2",
            "6.1.2",
            "SC3Pv3",
            Remote(BASE + "cell-exp/6.1.2/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5", "c2572841abc0db0a59ea086b038d5bd9dc272d9c1b43dbcd78d1811e7e1c3b71"),
            Remote(BASE + "cell-exp/6.1.2/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", "3a56663124a31aac49148f972514329c8064eb618a63566c3e933b86d3bc4eb3"),
            note="the raw matrix of pbmc10k_chromium_x_6_1, filtered again by 6.1.2; expect 10,000",
        ),
        Dataset(
            "mouse_brain_nuclei_7_0",
            "7.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/7.0.0/5k_mouse_brain_CNIK_3pv3/5k_mouse_brain_CNIK_3pv3_raw_feature_bc_matrix.h5", "ce224cfa07551e2091f907e0fb699d3acd171a77526fd30d0cbcde777a7e310b"),
            Remote(BASE + "cell-exp/7.0.0/5k_mouse_brain_CNIK_3pv3/5k_mouse_brain_CNIK_3pv3_filtered_feature_bc_matrix.h5", "169acce12b039401f36f14893e4b683ca8c8fe20eb10a0545f8f4f6b212a0748"),
            note="nuclei; expect not published",
        ),
        Dataset(
            "jejunum_nuclei_7_0",
            "7.0.0",
            "SC3Pv3",
            Remote(BASE + "cell-exp/7.0.0/5k_human_jejunum_CNIK_3pv3/5k_human_jejunum_CNIK_3pv3_raw_feature_bc_matrix.h5", "3f33def80c79d326f8f69d558f82d65c0c02ffb976068e09b6f84399945d4210"),
            Remote(BASE + "cell-exp/7.0.0/5k_human_jejunum_CNIK_3pv3/5k_human_jejunum_CNIK_3pv3_filtered_feature_bc_matrix.h5", "3ad88d316f2828fc98ded3934f53c80e91282a80804d85b3d6959fa7ba146fc9"),
            note="nuclei; expect not published",
        ),
        Dataset(
            "pbmc5k_7_0_1",
            "7.0.1",
            "SC3Pv3",
            Remote(BASE + "cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_raw_feature_bc_matrix.h5", "fd357aa54f63b794d451552a023cd5d215d28bd9948dc2a03f697ff74f51a6e7"),
            Remote(BASE + "cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5", "eddc0645d9adab0fe0b3a7bef5c3aa2b993c13d96d9e92afb2f3011e1bf02497"),
            note="expect not published",
        ),
        Dataset(
            "flex_spleen_ffpe_7_1",
            "7.1.0",
            "SFRP",
            Remote(
                BASE + "cell-exp/7.1.0/10k_mouse_spleen_scFFPE_singleplex_Multiplex/"
                "10k_mouse_spleen_scFFPE_singleplex_Multiplex_count_raw_feature_bc_matrix.h5",
                "a47cde71cbfdb0948256e2c3b56715fd07abc778619cff6967e0f782f99e26a6",
            ),
            Remote(
                BASE + "cell-exp/7.1.0/10k_mouse_spleen_scFFPE_singleplex_10k_mouse_spleen_scFFPE_singleplex/"
                "10k_mouse_spleen_scFFPE_singleplex_10k_mouse_spleen_scFFPE_singleplex_count_sample_filtered_feature_bc_matrix.h5",
                "f42ec2cae1c478966fa31fcb3acf22db1cb6a0316d3a65e6496b4a880462224d",
            ),
            note="singleplex Flex (cellranger multi); expect not published",
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
