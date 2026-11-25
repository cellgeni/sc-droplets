---
title: OrdMag
parent: Cell Ranger
grand_parent: Methods
nav_order: 2
method_ids: [cellranger-ordmag]
status: implemented
---

# Cell Ranger OrdMag
{: .no_toc }

"Order of magnitude": the cells are the barcodes within a factor of 10 of a robust maximum, the
99th percentile of the top expected cells. Cell Ranger 2 called cells with it alone; from 3.0 it
is the first step of [Cell Ranger cell calling](cellranger.md), and it is what Cell Ranger runs on
antibody-only libraries.

1. TOC
{:toc}

## At a glance

| Id | `cellranger-ordmag` |
|---|---|
| Tool | Cell Ranger 2.0 – 3.1 and 7.2 – 10.1 |
| Category | rank-threshold |
| Counts | UMIs, Gene Expression features, per genome and GEM group |
| Needs expected cells | 2.x/3.x: 3,000 per GEM group unless given; 7.2+: estimated unless given |
| Deterministic | yes: the bootstraps use Cell Ranger's seed |
| Original licence | MIT (3.0.2 source); 10x licence (7.2+) |
| Validated against | 10x's filtered matrices: pbmc4k (2.1.0); as the first step of [Cell Ranger cell calling](cellranger.md) for 3.0 – 10.0 |

## How it works

For one genome's features and one GEM group, with $$N$$ expected cells (`--expect-cells` divided
by the number of GEM groups, or 3,000):

1. $$b = \operatorname{round}(N \cdot (1 - 0.99))$$, halves rounded up (Python 2), capped at the
   number of nonzero barcodes − 1. For 3,000 expected cells $$b = 30$$.
2. Draw 100 bootstrap samples of the nonzero totals (with replacement, as many as there are),
   from numpy's legacy random stream seeded 0. In each sample, the baseline is its $$b$$-th
   largest value, the cut is $$\max(1, \operatorname{round}(0.1 \times \text{baseline}))$$, and
   $$n_i$$ counts the values $$\ge$$ the cut.
3. $$n = \operatorname{round}(\overline{n_i})$$.
4. The cells are the $$n$$ barcodes with the largest totals, as numpy's (pre-SIMD) quicksort
   orders them; ties at the boundary fall where that sort leaves them.

With $$b = 0$$ (fewer than 50 expected cells) the baseline is the *smallest* value, as Cell
Ranger indexes `[-0]`; 7.1 fixed this.

### From 7.2

- **Baseline**: the $$(b+1)$$-th largest, with $$b = \operatorname{round}(N \cdot q)$$, where
  $$q = 1 - 0.99$$ computed in floating point (slightly above 0.01), and numpy's rounding (halves
  to even). The cut is $$\max(1, \operatorname{round}(0.1 \times \text{baseline}))$$.
- **Random stream**: a fresh `RandomState(0)` for each genome and GEM group.
- **Expected cells**: given, they become $$\max(N, 50)$$; not given, they are estimated. For each
  of 100 bootstraps and each candidate $$r$$ on a grid of 2,000 log-spaced values from 2 to a cap
  (the low end of the [ambient range](cellranger.md#3-the-ambient-range-and-the-fdr-by-release-and-chemistry),
  at most 262,144), count the barcodes $$f(r)$$ within the order of magnitude of the
  $$(\operatorname{round}(r q) + 1)$$-th largest; keep the $$r$$ with the smallest
  $$(f(r) - r)^2 / r$$. The estimate is the mean of the 100 picks, rounded, at least 50.
- **Cells**: $$n_0$$, the mean count over 100 more bootstraps (continuing the stream), then a
  tie extension that, as written, always takes one barcode past the run of barcodes tied with the
  $$n_0$$-th: $$n_0 + 1$$ cells when there are no ties. If the extension would add more than 20 %
  of $$n_0$$, 7.2 keeps $$n_0$$ and 8.0+ keeps the extension up to 20 %.
- **Selection**: the top $$n$$ by a stable sort, ties to the barcode later in the matrix.

The genomes are taken in the order their features first appear; each genome's cells are called
on its own features, and the cells are their union.

## Versions

| Cell Ranger | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 1.0 – 1.3 | — | refused | no source; the paper's rule is "within an order of magnitude of the top 1 % of expected cells", without saying whether it bootstrapped |
| 2.0 – 3.1 | `v2` | `2.1.0`, `3.0.2` … | as above |
| 4.0 – 7.1 | — | refused | no source |
| 7.2 | `v7.2` | `7.2.0` | as "From 7.2"; tie extension reverted beyond 20 % |
| 8.0 – 10.1 | `v8.0`, `v8.0.1`, `v9.0`, `v10.0`, `v10.1` | `8.0.0` … `10.1.0` | tie extension truncated at 20 %; the estimate's cap follows each release's ambient ranges |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expect-cells` | 3,000 per GEM group (2.x/3.x); estimated (7.2+) | `--expect-cells` | expected number of cells |
| `--chemistry` | the input file's | `--chemistry` | chemistry, for the estimate's cap (7.2+) |
| `--n-probe-barcodes` | none | `probe_barcode_ids` | probe barcodes / OCM tags of this sample, for the cap (7.2+) |

## Run it

```bash
sc-droplets call cellranger-ordmag raw_gene_bc_matrices_h5.h5 --tool-version 2.1.0 --expect-cells 5000 --barcodes cells.tsv
sc-droplets plot diagnostics raw_gene_bc_matrices_h5.h5 -m cellranger-ordmag@2.1.0:expect_cells=5000 -o ordmag.png
```

```python
import sc_droplets as scd

adata = scd.read("raw_gene_bc_matrices_h5.h5")
scd.call(adata, "cellranger-ordmag@2.1.0", expect_cells=5000)
```

The diagnostics plot shows the cut, the robust maximum and the range of the 100 bootstraps' cell
counts.

## Differences from the original

- Ties at the boundary are ordered by a reimplementation of numpy's classic quicksort
  (`sc_droplets.primitives.legacysort`), because recent numpy sorts with SIMD code on some
  processors and leaves ties in a different order. On processors without it, the two agree.
- Cell Ranger first removes antibody aggregates; that belongs to
  [Cell Ranger cell calling](cellranger.md), not to OrdMag alone.

## Validation

On 10x's pbmc4k (Cell Ranger 2.1.0, raw `raw_gene_bc_matrices_h5.h5`), `--expect-cells 5000`
gives exactly the 4,340 barcodes of 10x's filtered matrix. 10x does not publish the option it ran
with; 5,000 is the value that reproduces it (1,000 to 4,000 call 19 to 236 fewer, 6,000 and above
call more).

## References

- Zheng et al. 2017, Massively parallel digital transcriptional profiling of single cells,
  *Nat Commun* 8:14049, [doi:10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049).
- Source (MIT): [`stats.py`](https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/stats.py)
  and [`filter_barcodes`](https://github.com/lfperales/cellranger/blob/3.0.2/mro/stages/counter/filter_barcodes/__init__.py),
  Cell Ranger 3.0.2.
