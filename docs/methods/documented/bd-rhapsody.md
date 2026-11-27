---
title: BD Rhapsody
parent: Documented only
grand_parent: Methods
nav_order: 3
method_ids: [bd-rhapsody]
status: docs-only
---

# BD Rhapsody: the second derivative of cumulative reads
{: .no_toc }

The putative-cell step of BD's Rhapsody Sequence Analysis Pipeline. Cell labels are sorted by
reads, and the cut is a deep minimum of the second derivative of the cumulative reads on a log
scale. An optional refined algorithm then removes and adds labels using highly variable genes
(mRNA) or Gaussian mixtures (ATAC). The code is closed and is known only from BD's documentation,
which leaves out constants. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `bd-rhapsody` (documented only) |
|---|---|
| Tool | BD Rhapsody Sequence Analysis Pipeline 1.6 – 4.0 (Docker image `bdgenomics/rhapsody`, CWL workflows, Seven Bridges) |
| Category | knee |
| Counts | reads per cell label, of the RSEC-adjusted molecules of the chosen bioproduct type (mRNA by default); for ATAC, transposase sites in peaks |
| Needs expected cells | no (`Expected_Cell_Count` and `Exact_Cell_Count` optional) |
| Deterministic | not stated; no randomness is described, though the Gaussian mixtures of the ATAC refined step may start at random |
| Original licence | closed (no software licence stated; the guide is "All rights reserved") |
| Extra inputs | reads per cell label, which a molecule matrix does not carry |

## How it works

The input is, per cell label, the number of reads of all RSEC-adjusted molecules (DBEC-adjusted in
the 1.12 handbook) of the selected bioproduct type: mRNA by default, AbSeq optionally. For ATAC it
is the number of transposase sites in peaks. The unit is reads, not molecules.

### Basic algorithm, 4.0 (the default since 2.0)

1. Sort the cell labels by reads, decreasing; take the cumulative reads on a $$\log_{10}$$ scale
   and its second derivative. No smoothing is mentioned.
2. Shift the second derivative so that it is negative everywhere and take its mean $$\bar m$$
   over the whole curve. Test the local minima against depths of $$2.0\,\bar m$$, then
   $$1.75\,\bar m$$, then $$1.5\,\bar m$$; the first level reached by at least one minimum ends
   the search. If none reaches $$1.5\,\bar m$$, the deepest minimum is taken, flagged as a sign of
   poor sample quality.
3. With one qualifying minimum, it is the cut. With several, the one giving the most cells such
   that the number of cells is below $$\max(1000,\ 0.6 \times \#\{\text{labels with} \ge 10 \text{ reads}\})$$
   and the label at the inflection point has more than 50 reads (75 for VDJ calling).
4. With `Expected_Cell_Count`, the inflection point closest to it is taken instead of step 3.
5. With `Exact_Cell_Count` $$N$$, the top $$N$$ labels by reads (by transposase sites in peaks for
   ATAC). For joint mRNA–ATAC calling, labels are ranked in each modality and the $$N$$ with the
   best mean rank are taken.

### Basic algorithm, 1.12 (handbook)

1. The same curve and second derivative.
2. An inflection point is a minimum of the second derivative. It is valid if it is at least half
   as deep as the global minimum and at most $$-0.3$$.
3. The smoothing window of the second derivative grows until there are two valid inflection
   points; the one giving more cells is taken, to include low-content populations (lymphocytes
   beside plasma cells, in BD's example). Neither the window sizes nor the smoother are given.

### Refined algorithm, mRNA and AbSeq (the default before 2.0)

1. The basic set $$A$$.
2. False positives: highly variable genes as in Macosko et al. ($$\log_{10}(\text{count} + 1)$$;
   per gene, the mean and the dispersion, variance over mean; 20 bins by mean; normalised
   dispersion = (dispersion − bin mean) / bin sd). For each of several cutoffs on the normalised
   dispersion (three in BD's example, values not given), the second-derivative analysis on those
   genes gives $$B_i$$; the noise is $$N_i = A \setminus B_i$$, and $$C = A \setminus \bigcap_i N_i$$.
3. False negatives: per gene, reads over all labels against reads in $$C$$, with a line of best
   fit. Genes whose residuals lie at least one standard deviation from the median residual are
   "enriched in noise"; the analysis on those genes gives $$D$$.
4. $$E = C \cup D$$, less the labels below a "dynamically determined minimum threshold" of
   molecules (not given).

### Refined algorithm, ATAC

Labels with less than 0.1 of their transposase sites in peaks are dropped. A count filter scans
thresholds on $$\log_{10}$$ site counts, from the 1st to the 99th percentile (the lower at least
1, that is 10 sites), in steps of 0.1, fitting a two-component Gaussian mixture at each and
keeping the most stable threshold. A two-component mixture on site counts then sets the putative cells (above the largest
count of the non-cell component), and one on site count and fraction in peaks may redefine the
cell cluster when its 90 % ellipses do not overlap. Last, one 2-D Gaussian is fitted to the
putative cells and one to the non-cells; each label's log-likelihood ratio, cut by the basic
algorithm without its half-depth rule, removes or recovers labels. That the half-depth rule is
named here suggests it still exists in 4.0 (an inference).

For multiome, the default joint call is the intersection of the mRNA and ATAC putative cells;
alternatively the call is made on one assay and cells absent from the other are dropped. Calls on
AbSeq flag suspected protein aggregates in `<sample>_Protein_Aggregates_Experimental.csv`.

## Versions

| Pipeline (date) | Change to the call, as the release notes give it |
|---|---|
| 1.6 (2019-06-10) | `Exact_Cell_Count` added; refined calling can be disabled |
| 1.9.1 (2020-10-06) | changed "to reduce overcalling of putative cells in high cell input experiments" |
| 1.10 (2022-01-24) | calling on AbSeq read counts, "for troubleshooting only" |
| 1.11 (2022-08-18) | intronic reads counted in molecules (with a parameter); protein-aggregate flag for AbSeq calls |
| 1.11.1 (2022-12-15) | "improved basic and refined putative cell calling algorithms" |
| 2.0 (2023-06-14) | the basic algorithm becomes the default, refined optional; `Expected_Cell_Count` added |
| 2.1 (2023-11-10) | basic calling changed to "prevent most inappropriately high cell calls derived from noise signatures" |
| 2.3 (2025-03-10) | `Exact_Cell_Count` fixed for ATAC-only and joint mRNA–ATAC calling |
| 3.0 (2025-10-29) | the basic algorithm updated "to fix several edge cases" (about 1 % more cells); VDJ-only calling |
| 4.0 (2026-09-15) | fix for calling that failed "when the last point of the curve was the only minimum" |

None of these changes is described beyond those lines. The guides for 2.0 – 3.0 could not be read,
so how the 1.12 rule became the 4.0 rule is known only from these lines.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `Cell_Calling_Data` | `mRNA` | what to call on: `mRNA`, `AbSeq`, `ATAC`, `mRNA_and_ATAC` or `VDJ` |
| `Cell_Calling_Bioproduct_Algorithm` | `Basic` | `Basic` or `Refined`, for mRNA and AbSeq (`Enable_Refined_Cell_Call` in 2.0; before 2.0 refined unless disabled) |
| `Cell_Calling_ATAC_Algorithm` | `Basic` | `Basic` or `Refined` |
| `Expected_Cell_Count` | none | take the inflection point nearest it |
| `Exact_Cell_Count` | none | call exactly the top $$N$$ labels ($$N \ge 1$$) |
| `Exclude_Intronic_Reads` | `False` | changes the molecules counted, and perhaps the reads per label (not stated) |

## Status

Documented only. The code is closed, and the documentation leaves out what a reimplementation
would need: the smoothing windows (1.12), how the derivative is taken, the dispersion cutoffs and
the minimum-molecule threshold of the refined step. No version of it could be validated. The rule
also ranks reads per cell label, which a molecule matrix does not carry.

BD's own call can be set beside the package's. `<sample>_RSEC_MolsPerCell_Unfiltered_MEX.zip`
holds every cell label with at least 10 reads; `<sample>_RSEC_MolsPerCell_MEX.zip` (`_DBEC_` for
targeted or AbSeq panels) holds the putative cells; in 1.x, `<sample>_Putative_Cells_Origin.csv`
says whether each came from the basic or the refined step. Run the package's methods on the
unfiltered matrix, mark the putative cells in a boolean `obs` column, and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls. Labels with fewer than 10 reads are absent from that matrix, so it holds
fewer empty droplets than a full raw matrix.

## References

- [BD Rhapsody Sequence Analysis Pipeline 4.0 documentation](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/):
  [putative cell calling](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_putative_cells.html),
  [parameters](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/setup/input/parameters.html),
  [release notes](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/release_notes.html),
  [outputs](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/outputs/top_outputs.html).
- BD Biosciences 2023, *BD Single-Cell Multiomics Bioinformatics Handbook*, 23-21713(04)
  (pipeline 1.12).
- Macosko et al. 2015, Highly parallel genome-wide expression profiling of individual cells using
  nanoliter droplets, *Cell* 161:1202–1214,
  [doi:10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002).
- Research notes: docs/_research/docs-only.md §4.
