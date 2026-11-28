---
title: dropbead
parent: Drop-seq
grand_parent: Methods
nav_order: 1
method_ids: [dropseq-dropbead]
status: implemented
---

# Drop-seq cell number (dropbead)
{: .no_toc }

dropbead's `estimateCellNumber`: on the cumulative fraction of reads of the top barcodes, the
barcode farthest from the diagonal. It returns a number of cells; the cells are that many barcodes
from the top.

Why the rule lands where it does, drawn on real data:
[dropbead's diagonal](../../explained/knee/dropbead.md), in *Explained*.

1. TOC
{:toc}

## At a glance

| Id | `dropseq-dropbead` |
|---|---|
| Tool | dropbead 0.3.1 (`estimateCellNumber`, and the knee drawn by `plotCumulativeFractionOfReads`) |
| Category | knee |
| Counts | reads per barcode (Drop-seq's `BAMTagHistogram`) in dropbead; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | GPL-3.0 |
| Validated against | dropbead 0.3.1 (commit `73fd121`), its own function in R (exact) |

## How it works

The counts sorted in decreasing order, $$c_1 \ge c_2 \ge \dots$$, and $$N$$ = `max.cells`
(15,000 by default):

1. $$x_i = i / N$$ and the cumulative fraction $$s_i = \sum_{k \le i} c_k / \sum_{k \le N} c_k$$
   for $$i = 1, \dots, N$$, as R's `cumsum` and `sum` compute them (in long double).
2. $$m = (s_N - s_1) / (N - x_1)$$ (dropbead's slope, which mixes the two scales), and the
   distances $$d_i = \sin(\arctan m)\,|s_i - x_i|$$.
3. The number of cells is the first $$i$$ with the largest $$d_i$$.

With fewer than $$N$$ barcodes dropbead's sums take missing values and it returns no number;
sc-droplets stops with an error.

## Versions

| dropbead | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.3.1 | `default` | | the rule above |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--max-cells` | 15,000 | `max.cells` | the top barcodes the curve is drawn over (dropbead's plot uses 10,000) |
| `--unit` | `umis` | reads per barcode | count UMIs (the matrix) or reads (`--reads`) |

## Run it

```bash
sc-droplets call dropseq-dropbead raw_feature_bc_matrix.h5 --json
sc-droplets call dropseq-dropbead raw_feature_bc_matrix.h5 --max-cells 10000 --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "dropseq-dropbead")
result.thresholds["cells"]
```

The thresholds hold the number of cells (`cells`) and the count of the last one (`umis` or
`reads`). `sc-droplets plot diagnostics` draws the rank curve with the cut.

## Differences from the original

- Drop-seq counts reads per barcode; sc-droplets counts UMIs by default, reads with `--reads`.
- dropbead takes the counts as sorted in its input file, so among equal counts the cells are the
  first in that file. Here equal counts keep the matrix's order.

## Validation

`oracles/dropseq/run.sh` sources dropbead's `R/Tools.R` at commit `73fd121` in R from conda-forge
and calls `estimateCellNumber` on the sorted UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw
matrices, with `max.cells` 15,000 and 10,000. Identical numbers of cells.

## References

- Alles, Karaiskos, Praktiknjo et al. 2017, Cell fixation and preservation for droplet-based
  single-cell transcriptomics, *BMC Biology* 15:44,
  [doi:10.1186/s12915-017-0383-5](https://doi.org/10.1186/s12915-017-0383-5).
- Source: [`R/Tools.R`](https://github.com/rajewsky-lab/dropbead/blob/73fd121/R/Tools.R)
  (`estimateCellNumber`).
