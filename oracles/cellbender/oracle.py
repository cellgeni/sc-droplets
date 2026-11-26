"""Build one CellBender release's own SingleCellRNACountsDataset on a raw matrix and write its priors.

    python oracle.py RELEASE RAW.h5 OUT_DIR EXPECTED_CELLS

CellBender's source for RELEASE must be on PYTHONPATH (run.sh checks out the tag). The dataset
object is what `cellbender remove-background` builds before training: it loads the file, works
out the priors and chooses the droplets to analyse. Writes OUT_DIR/priors.tsv
(`case<TAB>name<TAB>value`, or `case<TAB>error<TAB>message`), OUT_DIR/<case>.classes.txt.gz
(`barcode<TAB>class` for the probable cells, the other analysed droplets and the surely-empty
droplets, in CellBender's order) and OUT_DIR/versions.txt.

Cases: default; expected (--expected-cells); total (--total-droplets-included 15000); low30
(--low-count-threshold 30); from 0.3.0 also forced (--force-cell-umi-prior 3000
--force-empty-umi-prior 30) and expected_total (both of the first two).
"""

import gzip
import logging
import os
import platform
import sys

import numpy as np
import scipy
import torch

from cellbender.remove_background.data.dataset import SingleCellRNACountsDataset

logging.disable(logging.CRITICAL)


def version_tuple(release):
    return tuple(int(p) for p in release.split("."))


def cases(release, expected):
    out = [
        ("default", {}),
        ("expected", {"expected_cell_count": expected}),
        ("total", {"total_droplet_barcodes": 15000}),
        ("low30", {"low_count_threshold": 30}),
    ]
    if version_tuple(release) >= (0, 3, 0):
        out += [
            ("forced", {"force_cell_umi_prior": 3000.0, "force_empty_umi_prior": 30.0}),
            ("expected_total", {"expected_cell_count": expected, "total_droplet_barcodes": 15000}),
        ]
    return out


def build(release, path, options):
    """The dataset object as each release's command line builds it, with its defaults."""
    v = version_tuple(release)
    if v < (0, 3, 0):
        kwargs = dict(input_file=path, expected_cell_count=None, total_droplet_barcodes=25000,
                      fraction_empties=0.5, model_name="full", gene_blacklist=[], low_count_threshold=15)
        if v >= (0, 2, 0):
            kwargs.update(exclude_antibodies=False, fpr=[0.01])
    else:
        kwargs = dict(input_file=path, model_name="full", exclude_features=[], low_count_threshold=5,
                      fpr=[0.01], expected_cell_count=None, total_droplet_barcodes=None,
                      force_cell_umi_prior=None, force_empty_umi_prior=None, fraction_empties=0.2,
                      ambient_counts_in_cells_low_limit=0.1, gene_blacklist=[])
    kwargs.update(options)
    return SingleCellRNACountsDataset(**kwargs)


def values_of(release, dataset):
    priors = dataset.priors
    if version_tuple(release) < (0, 3, 0):
        n_cells = min(int(priors["n_cells"]), int(np.sum(np.asarray(dataset.data["matrix"].sum(axis=1)).ravel() > 0)))
        names = ["n_cells", "cell_counts", "empty_counts", "log_counts_crossover", "d_std", "cell_prob", "cell_logit"]
        values = {name: priors[name] for name in names}
        if hasattr(dataset, "empty_UMI_threshold"):
            values["empty_UMI_threshold"] = dataset.empty_UMI_threshold
    else:
        n_cells = int(priors["expected_cells"])
        names = ["cell_counts", "empty_counts", "empty_count_upper_limit", "expected_cells", "total_droplets",
                 "transition_point", "log_counts_crossover", "surely_empty_counts", "d_std", "cell_logit"]
        values = {name: priors[name] for name in names}
        values["low_count_cutoff"] = dataset.low_count_cutoff
        values["empty_UMI_threshold"] = dataset.empty_UMI_threshold
        values["n_features"] = len(dataset.analyzed_gene_inds)
    values["n_analysed"] = len(dataset.analyzed_barcode_inds)
    values["n_empty"] = len(dataset.empty_barcode_inds)
    return n_cells, values


def main(release, path, out, expected):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "versions.txt"), "w") as handle:
        handle.write("cellbender %s\npython %s\nnumpy %s\nscipy %s\ntorch %s\nplatform %s %s\n" % (
            release, platform.python_version(), np.__version__, scipy.__version__, torch.__version__,
            sys.platform, platform.machine()))
    with open(os.path.join(out, "priors.tsv"), "w") as handle:
        for name, options in cases(release, expected):
            try:
                dataset = build(release, path, options)
            except (AssertionError, IndexError, ValueError) as error:
                message = str(error).splitlines()[0] if str(error) else type(error).__name__
                handle.write("%s\terror\t%s\n" % (name, message))
                print("%s %s: error %s" % (release, name, message))
                continue
            n_cells, values = values_of(release, dataset)
            for key, value in values.items():
                handle.write("%s\t%s\t%r\n" % (name, key, float(value)))
            barcodes = dataset.data["barcodes"]
            analysed = np.asarray(dataset.analyzed_barcode_inds, dtype=int)
            empty = np.asarray(dataset.empty_barcode_inds, dtype=int)
            with gzip.open(os.path.join(out, "%s.classes.txt.gz" % name), "wt") as classes:
                for rank, index in enumerate(analysed):
                    classes.write("%s\t%s\n" % (_text(barcodes[index]), "cell" if rank < n_cells else "transition"))
                for index in empty:
                    classes.write("%s\tempty\n" % _text(barcodes[index]))
            print("%s %s: %d cells, %d analysed, %d empty" % (release, name, n_cells, len(analysed), len(empty)))


def _text(barcode):
    return barcode.decode() if isinstance(barcode, bytes) else str(barcode)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
