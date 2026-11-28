---
title: Reproducing a pipeline's calls
parent: Guides
nav_order: 2
---

# Reproducing a pipeline's calls
{: .no_toc }

To get the cells a pipeline called, from its raw matrix: find the release that ran, find the
options it ran with, run the matching method with both, and check the result against the
pipeline's own filtered barcodes.

1. TOC
{:toc}

## Find the release

| Tool | Where the release is written |
|---|---|
| Cell Ranger | the `software_version` attribute of the HDF5 matrices of recent releases (Cell Ranger 3.0's have none); otherwise the web summary or the pipeline's log |
| STARsolo | the first line of `Log.out` (`STAR version=2.7.10a`); the `soloCellFilter` line shows the rule and its numbers |
| DropletUtils | `packageVersion("DropletUtils")` in the R session that ran it, and the seed set before the call |
| salmon / alevin, alevin-fry, kb, bustools, UMI-tools, CellBender, cellsweep | the tool's `--version`, or its log |

```bash
python -c "import h5py; print(h5py.File('raw_feature_bc_matrix.h5').attrs.get('software_version'))"
```

The release matters: the same matrix gives different cells under different releases. Each
method's *Versions* section lists what changed and the `--tool-version` values it accepts. A
release that cannot be reproduced, such as Cell Ranger 6.0 to 7.1, is refused with a pointer to
its page.

## Map the pipeline to a method

| The pipeline ran | Run |
|---|---|
| `cellranger count` | `cellranger`, with `--expect-cells` or `--force-cells` if they were given |
| `cellranger multi` (Flex, OCM) | `cellranger` on each sample: see [Multiplexed and Flex runs](multiplexed.md) |
| STARsolo `--soloCellFilter CellRanger2.2 N p r` | `starsolo-cellranger2.2 --expected-cells N --max-percentile p --max-min-ratio r` |
| STARsolo `--soloCellFilter EmptyDrops_CR ...` | `starsolo-emptydrops-cr`, its ten numbers as `--expected-cells` to `--sim-n` |
| STARsolo `--soloCellFilter TopCells N` | `starsolo-topcells --n-cells N` |
| DropletUtils `emptyDrops`, `barcodeRanks`, `defaultDrops` | `dropletutils-emptydrops`, `dropletutils-barcoderanks`, `dropletutils-defaultdrops`, with `--seed` for emptyDrops |
| alevin-fry `generate-permit-list` | `alevinfry-permit-list` (`--expect-cells`, `--force-cells`, `--unfiltered` as given) |
| salmon alevin's initial whitelist | `alevin-initial-whitelist` |
| kb `count --filter bustools` | `bustools-allowlist` |
| UMI-tools `whitelist --knee-method` | `umitools-distance` or `umitools-density` |
| CellBender `remove-background` | `cellbender-prior` reproduces its priors; its own cells come from its [model](../methods/cellbender/model.md) |
| cellsweep, zUMIs, Drop-seq (dropbead), dropEst | `cellsweep-empty`, `zumis-cell-barcodes`, `dropseq-dropbead`, `dropest-bounds` |

```bash
sc-droplets call cellranger raw_feature_bc_matrix.h5 --tool-version 7.2.0 --barcodes cells.tsv.gz
sc-droplets call starsolo-emptydrops-cr Solo.out --tool-version 2.7.10a --barcodes cells.tsv.gz
sc-droplets call dropletutils-emptydrops raw_feature_bc_matrix.h5 --tool-version 1.26 --seed 100 --table emptydrops.csv
```

## Match the inputs too

The original's call depends on more than its options:

- **Counts.** UMI-tools, alevin, alevin-fry, bustools, zUMIs and Drop-seq rank reads by default. A
  count matrix holds UMIs; give the reads with `--reads` and count them with `--unit reads`
  ([Read-based methods](reads.md)).
- **Features.** Most tools count Gene Expression only, as the package does by default. A
  pipeline that counted everything needs `--feature-type all`.
- **Genome.** In a barnyard, Cell Ranger calls each genome on its own features; `--genome`
  counts one.
- **Matrix.** STARsolo's `Gene` and `GeneFull` give different totals; read the one the pipeline
  filtered with `--solo-feature`.

## Check the result

Compare the cells with the pipeline's own filtered barcodes:

```python
import gzip
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "cellranger@7.2.0")
ours = set(result.barcodes[result.is_cell])
theirs = set(gzip.open("filtered_feature_bc_matrix/barcodes.tsv.gz", "rt").read().split())
len(ours & theirs), len(ours - theirs), len(theirs - ours)
```

For a deterministic method the sets should be the same. Where they are not, the usual causes are,
in order: another release or another option than the one assumed; other counts (reads, features,
genome); and ties at the cut. Cell Ranger picks some barcodes by numpy's default sort, which on
recent numpy with AVX-512 leaves tied barcodes in another order. For a Monte Carlo test, a few
barcodes at the FDR cut can differ when the original's random numbers came from another build,
such as STAR compiled with clang. [Reads, UMIs, ties and rounding](../overview/arithmetic.md)
lists these details tool by tool. Each method's *Validation* section shows how close the package
came on public data.
