---
title: BAP bead knee
parent: Planned
grand_parent: Methods
nav_order: 14
method_ids: [bap-knee]
status: planned
---

# BAP bead knee: a density valley on log bead counts
{: .no_toc }

The first step of BAP (bead-based ATAC processing), written for droplet ATAC with several beads
per droplet (Bio-Rad's dscATAC-seq). The counts per bead are smoothed on the log scale with R's
`density`, a valley chosen by fixed rules sets the threshold, and a guard falls back to 500 when
the threshold looks implausible. It follows the density valley of
[UMI-tools](../umitools/density.md) with other constants and rules. BAP then merges beads into
droplets, which is not part of this call. Planned for 0.3.

1. TOC
{:toc}

## At a glance

| Id | `bap-knee` (planned) |
|---|---|
| Tool | BAP, tag `biorad-paper1` (2018) to 0.7.1 (`master`, 2026), entry points `bap` and `bap2`; the knee unchanged since 2018-07-26 |
| Category | density |
| Counts | per bead: `bap`, distinct read starts (MAPQ above 30, nuclear chromosomes); `bap2`, distinct fragments (mitochondrial included) |
| Needs expected cells | no |
| Deterministic | yes; R 4.4.0 moved `density`'s grid, so the threshold can depend on the R version |
| Original licence | MIT |
| Extra inputs | a BAM with a bead tag, for BAP's own counts; the knee needs only one count per bead |

## How it works

The input is one count $$c_b$$ per bead.

1. **Floor.** $$f = 0.001 \times \operatorname{mode}(c)$$, the most frequent count (ties: the
   value seen first in file order, which for `bap2` is awk's hash order). Beads with $$c > f$$
   are kept; with a mode below 1,000, that is every bead with a count.
2. **Log.** $$z = \log_{10} c$$.
3. **Density.** R's `density` with a Gaussian kernel and an absolute bandwidth of 0.1 (in
   $$\log_{10}$$ units), from $$\min z$$ to $$\max z$$ at 10,000 points: R's binned FFT estimate,
   interpolated to $$x_1 < \dots < x_{10000}$$ with values $$y_i$$. R 4.4.0 changed the grid
   coordinates, moving the values by about 0.1 % (`old.coords = TRUE` restores the old ones).
4. **Candidates.** Each step is "down" when $$y_i - y_{i-1} < 0$$ (with $$y_0 = +\infty$$, so the
   first step is down; flat steps count as up). The candidates are the ends of the runs of down
   and up steps. With more than 3 runs, every other end is kept, starting with the first: the
   ends of the falling runs, the local minima. With 3 or fewer runs all ends are kept, so a local
   maximum can be a candidate. If $$y_1 = y_2$$ the first candidate is dropped; so are candidates
   at index 1 or 10,000.
5. **Choice.** Scanning the candidates from the right (largest count first), the first index
   $$i$$ with $$i \ge 2{,}000$$ (20 % of the grid; the code's comment says "of barcodes") and
   either $$\max z - x_i > 0.5$$ or $$x_i < \max z / 2$$.
6. **Threshold.** $$T = 10^{x_i}$$, or 0 if no candidate qualifies.
7. **Guard.** If $$T > 100{,}000$$ or $$T < 100$$, the threshold is 500. Both values are written
   (`bead_threshold`, `bead_threshold_nosafety`).
8. **Call.** The beads with $$c_b \ge$$ the guarded threshold.

No random numbers. `--minimum-barcode-fragments N` with $$N > 0$$ skips the knee and keeps
$$c_b \ge N$$; `--barcode-whitelist` takes the listed beads.

Against UMI-tools' density valley: the floor is 0.001 × the mode, not the maximum; the bandwidth
is absolute, not 0.1 × the standard deviation; minima come from run lengths, not from comparing
neighbours; beads at the threshold are called ($$\ge T$$, not $$> T$$); a threshold outside
[100, 100,000], including the 0 of a failed search, becomes 500; there is no expected-cells mode.

After this step BAP computes the Jaccard index of Tn5 insertions between called beads, applies
the same function to the largest $$10^6$$ values (guard outside $$[10^{-6}, 0.5]$$ → 0.005)
unless `--minimum-jaccard-index` is given, and merges beads into droplet barcodes. That is a
different question (which beads share a droplet) and is not part of this method.

## Versions

The knee (`00_knee_CL.R`) and its caller (`10b_knee_execute.R`) last changed in commit `f9d30a4`
(2018-07-26) and are byte-identical at `biorad-paper1`, `v0.6.0` and `master`. What changed is
the count fed to them:

| Entry point | Releases | Count per bead |
|---|---|---|
| `bap` | 2018 – today | reads with MAPQ above 30 on chromosomes of at least 1,000,000 bp (the nuclear ones), once per read start and bead |
| `bap2` | 0.6.0 (2019-07) – today | distinct fragments (chromosome, start and end shifted for Tn5 by +4 and −5) from proper pairs, insert at most 2,000 bp, pair score at least 30 (taken to be the lower mate's MAPQ, not checked), off the blacklist; summed over every chromosome, the mitochondrial one included |

Tags: `bap-barcode` (2018-06-21), `biorad-paper1` (`7872481`, 2018-08-06, the paper's code),
`v0.6.0` (`8b35816`, 2019-07-26). PyPI's `bap-atac` runs from 0.6.4 (2019-10-29) to 0.7.0
(2021-08-29); 0.7.1 is only in the repository. Since the knee is a function of the count vector,
one method can take the count as input and document both definitions. R before 4.4.0 and from
4.4.0 would be two Variants of the dependency.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--minimum-barcode-fragments` | 0 (the knee) | above 0, the threshold itself |
| `--barcode-whitelist` | none | call the listed beads instead |
| `--mapq` | 30 | read (`bap`) or pair (`bap2`) quality |
| `--min_chromosome_size` | 1,000,000 bp | `bap`: the chromosomes counted |
| `--max-insert` | 2,000 | `bap2`: the longest fragment |
| `--all-pairs` | off | `bap2`: keep pairs that are not proper |
| `--blacklist-file`, `--reference-genome` | the genome's | `bap2`: regions left out |
| `--bead-tag` | `XB` | the BAM tag holding the bead |
| fixed | 0.001; bandwidth 0.1; 10,000 points; index 2,000; 0.5; [100, 100,000] → 500 | steps 1–7 |

## Status

Planned for 0.3. Implementing it needs:

- a count per bead. The knee runs on any per-barcode total, so the package's UMI totals, or
  the unique fragments per barcode in 10x's `singlecell.csv` (a stand-in for `bap2`'s count),
  serve; BAP's own counts need a BAM reader with bead tags, which the package does not have;
- a port of R's `density.default` (linear binning, FFT, interpolation) in both grid conventions;
  the rest is integer arithmetic and exact comparisons;
- the input order of the counts, which settles ties for the mode.

The oracle runs the function itself, as the UMI-tools oracle does: `Rscript
bap/bin/R/10b_knee_execute.R counts.csv 1 V2` from tag `v0.6.0` (the same as `master`) in a pixi
environment (conda-forge R, r-data.table, r-dplyr), in R 4.3 and in R 4.4 or later, on count
vectors: `passed_filters` from 10x ATAC `singlecell.csv` files and per-barcode totals of RNA
matrices. A pipeline-level check would need a BAM with bead tags (the paper's data, GEO
GSE123581, not checked) and is not needed. No blockers are known beyond the `density` port.

## References

- Lareau et al. 2019, Droplet-based combinatorial indexing for massive-scale single-cell
  chromatin accessibility, *Nat Biotechnol* 37:916–924,
  [doi:10.1038/s41587-019-0147-6](https://doi.org/10.1038/s41587-019-0147-6).
- Source (MIT): [caleblareau/bap](https://github.com/caleblareau/bap/tree/v0.6.0) at tag
  `v0.6.0`, `bap/bin/R/00_knee_CL.R` and `bap/bin/R/10b_knee_execute.R`; the paper's code at tag
  [`biorad-paper1`](https://github.com/caleblareau/bap/tree/biorad-paper1) (the paper cites the
  fork `buenrostrolab/bap`).
- Research notes: docs/_research/planned-atac.md §4.
