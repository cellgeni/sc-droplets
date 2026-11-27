"""velocyto's loom files."""
from __future__ import annotations


import numpy as np
import pytest
from conftest import write_loom

import sc_droplets as scd
from sc_droplets.io.loom import is_loom


def test_a_loom_reads_as_velocyto_writes_it(spliced, tmp_path):
    path = write_loom(spliced, tmp_path / "sample.loom")
    assert is_loom(path)
    adata = scd.read(path)
    assert adata.obs_names[0] == "sample1:" + spliced.obs_names[0].split("-")[0] + "x"
    assert (adata.layers["spliced"] != spliced.layers["spliced"]).nnz == 0
    assert (adata.layers["unspliced"] != spliced.layers["unspliced"]).nnz == 0
    assert (adata.X != spliced.layers["spliced"] + spliced.layers["unspliced"]).nnz == 0


def test_a_loom_needs_spliced_and_unspliced_layers(spliced, tmp_path):
    import h5py

    path = write_loom(spliced, tmp_path / "sample.loom")
    with h5py.File(path, "a") as handle:
        del handle["layers/unspliced"]
    with pytest.raises(ValueError, match="layers/unspliced"):
        scd.read(path)
