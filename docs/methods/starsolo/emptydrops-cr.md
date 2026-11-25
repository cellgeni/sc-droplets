---
title: EmptyDrops_CR
parent: STARsolo
grand_parent: Methods
nav_order: 3
method_ids: [starsolo-emptydrops-cr]
status: implemented
---

# STARsolo EmptyDrops_CR
{: .no_toc }

STAR's reimplementation of Cell Ranger 3's two-step calling: the CellRanger2.2 cut, then a
multinomial test of the barcodes just below it against the ambient profile of the droplets ranked
45,000 to 90,000. It follows Cell Ranger's description, not its code, and differs in the details
that move borderline calls: the Good-Turing estimate, the candidates' minimum, the random numbers,
and the order of ties.

1. TOC
{:toc}

## At a glance

| Id | `starsolo-emptydrops-cr` |
|---|---|
| Tool | STAR 2.7.8a and later (`--soloCellFilter EmptyDrops_CR`) |
| Category | statistical-test |
| Counts | UMIs of the detected barcodes, every feature |
| Needs expected cells | uses 3,000 unless told otherwise |
| Deterministic | yes: each simulation has a fixed seed |
| Original licence | MIT |
| Validated against | `STAR --runMode soloCellFiltering` 2.7.10a and 2.7.11b (exact) |

## How it works

1. **Initial cells.** The [CellRanger2.2](cellranger2-2.md) cut (v3); the barcodes above it are
   nCellsSimple cells.
2. **Order.** Detected barcodes by decreasing total, ties by increasing position in the barcode
   list.
3. **Ambient profile.** Sum the counts of the barcodes at ranks `indMin` to `indMax − 1` (45,000
   to 89,999, 0-based) per feature. With $$n_r$$ features summing to $$r$$ (features detected in
   no barcode at all are not counted among the zeros), STAR's Simple Good-Turing gives each
   feature with sum $$r > 0$$ the probability of one item seen $$r$$ times, and shares
   $$P_0 = n_1/N$$ equally among the detected features with sum 0. The profile is normalised to sum to 1.
   STAR's Good-Turing switches from the Turing estimate to the smoothed one at the first
   frequency whose successor is missing, or where the two agree within 1.96 standard deviations;
   it needs at least 5 distinct frequencies.
4. **Candidates.** The barcodes following the initial cells in that order, up to `candMaxN`
   (20,000), stopping at the first with fewer than
   $$\max(\text{umiMin},\ \lfloor \text{umiMinFracMedian} \cdot u_{(\lfloor nCellsSimple/2 \rfloor)} \rfloor)$$
   UMIs: 500, or 1 % of the middle initial cell's total, truncated.
5. **Observed likelihood.** For a candidate with counts $$x$$ and total $$n$$,
   $$\log L = \log n! - \sum_g \log x_g! + \sum_g x_g \log p_g,$$
   with log-factorials from a table of cumulative sums of $$\log k$$.
6. **Simulations.** For $$s = 0 \dots \text{simN} - 1$$, a `std::mt19937` seeded
   $$19760110\,(s+1)$$ draws one feature at a time from `std::discrete_distribution` over the
   profile's nonzero entries, up to the largest candidate's total; after each draw the simulated
   log-likelihood grows by $$\log p_g + \log n - \log c_g$$ ($$c_g$$: how often $$g$$ has been
   drawn).
7. **p-values.** $$p = (1 + \#\{s : L_s(n) < \log L\}) / (1 + \text{simN})$$, at the candidate's
   own total.
8. **Adjustment.** Benjamini-Hochberg without capping at 1; candidates with adjusted
   $$p \le \text{FDR}$$ (0.01) are cells.

STAR stops early, adding no cells, when there are no more than `indMin` detected barcodes, when
the ambient droplets hold no counts, or when there is no candidate.

## Versions

| STAR | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| before 2.7.8a | — | — | no EmptyDrops_CR |
| 2.7.8a – 2.7.11b | `v1` | `2.7.8a` | as above; 2.7.10a only changed a log message |

DropletUtils' `emptyDropsCellRanger` (1.14 onwards) is a port of this STAR code to R.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expected-cells` | 3000 | `nExpectedCells` | expected number of cells |
| `--max-percentile` | 0.99 | `maxPercentile` | robust maximum percentile |
| `--max-min-ratio` | 10 | `maxMinRatio` | robust maximum ÷ this is the initial cut |
| `--ind-min` | 45000 | `indMin` | first ambient rank (0-based) |
| `--ind-max` | 90000 | `indMax` | rank after the last ambient droplet |
| `--umi-min` | 500 | `umiMin` | fewest UMIs a candidate may have |
| `--umi-min-frac-median` | 0.01 | `umiMinFracMedian` | fraction of the middle initial cell's total a candidate needs |
| `--cand-max-n` | 20000 | `candMaxN` | most candidates |
| `--fdr` | 0.01 | `FDR` | largest adjusted p-value called a cell |
| `--sim-n` | 10000 | `simN` | number of simulations |

The parameters are in STAR's order: `--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500
0.01 20000 0.01 10000`.

## Run it

```bash
sc-droplets call starsolo-emptydrops-cr Solo.out -o calls.h5ad
sc-droplets call starsolo-emptydrops-cr raw_feature_bc_matrix.h5 --fdr 0.001 --table calls.tsv
```

```python
import sc_droplets as scd

adata = scd.read("Solo.out")
result = scd.call(adata, "starsolo-emptydrops-cr")
adata.obs["starsolo_emptydrops_cr_class"].value_counts()
```

Each barcode's class is `initial` (above the cut), `rescued`, `candidate` (tested, not called) or
`background`; candidates also get `pvalue`, `fdr` and `loglik` columns.

## Differences from the original

- The random draws are libstdc++'s (GCC), which STAR's Linux binaries use. A STAR built with
  libc++ (clang, macOS) draws differently and so calls slightly different borderline cells.
- Where STAR's Good-Turing fails (fewer than 5 distinct ambient counts), STAR's result is
  undefined; this method adds no cells and says so.
- STAR could read past its arrays when there are fewer candidates than `candMaxN` at the end of
  the list; this stops at the last barcode.

## Validation

Compared with `STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR` (2.7.10a and
2.7.11b) on two public 10x raw matrices; the cell sets are identical. The building blocks are
also checked one by one: the Good-Turing port against STAR's `sgt.h` compiled with GCC, and the
random draws against libstdc++'s `discrete_distribution` with `mt19937`.

## References

- Kaminow, Yunusov & Dobin 2021, STARsolo, bioRxiv,
  [doi:10.1101/2021.05.05.442755](https://doi.org/10.1101/2021.05.05.442755).
- Lun et al. 2019, EmptyDrops, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- Gale & Sampson 1995, Good-Turing frequency estimation without tears, *J Quant Linguist* 2:217,
  [doi:10.1080/09296179508590051](https://doi.org/10.1080/09296179508590051).
- Source: [`SoloFeature_emptyDrops_CR.cpp`](https://github.com/alexdobin/STAR/blob/2.7.11b/source/SoloFeature_emptyDrops_CR.cpp),
  [`sgt.h`](https://github.com/alexdobin/STAR/blob/2.7.11b/source/SimpleGoodTuring/sgt.h) (2.7.11b).
