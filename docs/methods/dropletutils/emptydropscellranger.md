---
title: emptyDropsCellRanger
parent: DropletUtils
grand_parent: Methods
nav_order: 3
method_ids: [dropletutils-emptydropscellranger]
status: implemented
---

# DropletUtils emptyDropsCellRanger
{: .no_toc }

An R rewrite of [STARsolo's EmptyDrops_CR](../starsolo/emptydrops-cr.md), itself STAR's
reading of Cell Ranger 3: a robust-maximum cut for the initial cells, an ambient profile from the
droplets ranked 45,000 to 90,000, and the EmptyDrops test of the droplets just below the cut. It
runs on DropletUtils' own machinery (edgeR's Good-Turing, R's indexing and rounding, DropletUtils'
random numbers), so it does not give STAR's cells.

1. TOC
{:toc}

## At a glance

| Id | `dropletutils-emptydropscellranger` |
|---|---|
| Tool | DropletUtils 1.14 and later (`emptyDropsCellRanger()`) |
| Category | statistical-test |
| Counts | the matrix as given (UMIs), every feature with a count |
| Needs expected cells | uses 3,000 unless told otherwise |
| Deterministic | yes, given R's seed (default 100) |
| Original licence | GPL-3 |
| Validated against | `DropletUtils::emptyDropsCellRanger` in R with the same seed (exact) |

## How it works

With totals $$t$$ (`as.integer(round(colSums(m)))`, after rounding the matrix) sorted in
decreasing order, ties in matrix order, $$t_{(1)} \ge t_{(2)} \ge \dots$$ (1-based, every column
counted):

1. $$m = \text{round}(N_{\text{exp}} (1 - q))$$ with R's half-to-even `round`: 30 for 3,000 cells
   and $$q = 0.99$$.
2. **Initial cut** `retain` $$= \max(\text{round}(t_{(m)} / 10),\ 1)$$; the $$n_s$$ droplets with
   $$t \ge$$ `retain` are initial cells.
3. **Candidates' minimum** $$u = \max(500,\ \text{round}(0.01 \cdot t_{(\lfloor n_s/2 \rfloor)}))$$.
   In 1.14 – 1.18.0 the total at that position is taken from the matrix in its own order, not
   from the sorted totals.
4. **Ambient droplets** are those at sorted positions $$\min(45000, n)$$ to $$\min(90000, n)$$,
   both ends included; their profile is edgeR's Good-Turing estimate with DropletUtils' pseudo-count
   for zero probabilities (as [emptyDrops](emptydrops.md) step 3).
5. **Tested droplets** are those at sorted positions $$n_s + 1$$ to
   $$\min(n_s + 20000,\ \#\{t > u\})$$, written as R's `a:b`, which counts down when $$b < a$$.
6. **Test**: emptyDrops' multinomial Monte Carlo ($$\alpha = \infty$$, 10,000 iterations, the same
   random numbers).
7. Droplets with $$t \ge$$ `retain` get $$p = 0$$; Benjamini-Hochberg over the tested droplets.

The cells are the droplets with FDR $$\le$$ `--fdr` (0.01, STARsolo's and Cell Ranger's value).

### Where it differs from STARsolo

| | STARsolo 2.7.11b | emptyDropsCellRanger |
|---|---|---|
| Robust maximum | the 31st largest total (0-based index 30) | the 30th largest |
| Rounding | $$\lfloor x + 0.5 \rfloor$$ | half to even |
| Barcodes ranked | detected barcodes only | every column |
| Ambient ranks | 45,000 to 89,999 (0-based) | 45,000 to 90,000 (1-based, inclusive) |
| Candidates | stop at the first total below $$u$$ | every total above $$u$$ |
| Good-Turing | STAR's `sgt.h` | edgeR's, plus a pseudo-count |
| Random numbers | `std::mt19937` per simulation | PCG32 per iteration, seeded from R |
| p-values | $$\#\{\text{sim} < \text{obs}\}$$ | $$\#\{\text{sim} \le \text{obs}\}$$ |
| Adjustment | not capped at 1 | `p.adjust` |

## Versions

| DropletUtils | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| before 1.14 | — | — | no `emptyDropsCellRanger` |
| 1.14 – 1.18.0 | `v1.14` | `1.14` | candidates' minimum from the unsorted totals |
| 1.18.1 – 1.34 | `v1.18.1` | `1.18.1` | as above |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expected-cells` | 3000 | `n.expected.cells` | expected number of cells |
| `--max-percentile` | 0.99 | `max.percentile` | robust maximum percentile |
| `--max-min-ratio` | 10 | `max.min.ratio` | robust maximum ÷ this is the initial cut |
| `--umi-min` | 500 | `umi.min` | fewest UMIs a candidate needs |
| `--umi-min-frac-median` | 0.01 | `umi.min.frac.median` | fraction of the middle initial cell's total a candidate needs |
| `--cand-max-n` | 20000 | `cand.max.n` | most candidates |
| `--ind-min` | 45000 | `ind.min` | first ambient position (1-based) |
| `--ind-max` | 90000 | `ind.max` | last ambient position |
| `--round-counts` / `--no-round-counts` | rounds | `round` | round non-integer counts |
| `--niters` | 10000 | `niters` | Monte Carlo iterations |
| `--fdr` | 0.01 | | largest FDR called a cell |

## Run it

```bash
sc-droplets call dropletutils-emptydropscellranger raw_feature_bc_matrix.h5 -o calls.h5ad
sc-droplets call dropletutils-emptydropscellranger raw_feature_bc_matrix.h5 --tool-version 1.14 --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
scd.call(adata, "dropletutils-emptydropscellranger", expected_cells=5000)
```

Classes are `initial`, `rescued`, `candidate` (tested, not called) and `background`; tested
droplets get `pvalue`, `fdr`, `logprob` and `limited` columns.

## Differences from the original

None known; see [emptyDrops](emptydrops.md#differences-from-the-original) for what reproducing
R's arithmetic takes.

## Validation

Compared with `DropletUtils::emptyDropsCellRanger` 1.14.2, 1.18.0, 1.20.0, 1.26.0 and 1.30.0
from bioconda on 10x's pbmc_1k_v3 raw matrix with `set.seed(100)`: identical p-values, FDR and
calls. On small matrices where edgeR's Good-Turing has too few distinct counts, R stops with an
error, and so does this.

## References

- Kaminow, Yunusov & Dobin 2021, STARsolo, bioRxiv,
  [doi:10.1101/2021.05.05.442755](https://doi.org/10.1101/2021.05.05.442755).
- Lun et al. 2019, EmptyDrops, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- Source: [`R/emptyDropsCellRanger.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/emptyDropsCellRanger.R)
  (1.30); the bug fix is commit `49071509` (1.18.1).
