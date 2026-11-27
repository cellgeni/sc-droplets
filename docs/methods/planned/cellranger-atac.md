---
title: Cell Ranger ATAC
parent: Planned
grand_parent: Methods
nav_order: 11
method_ids: [cellranger-atac]
status: planned
---

# Cell Ranger ATAC: a negative-binomial mixture on fragments in peaks
{: .no_toc }

`cellranger-atac count`'s cell caller. Barcodes that look like low-targeting noise, gel-bead
doublets or barcode multiplets are set aside. A two-component negative-binomial mixture is fitted
by EM to the fragments in peaks of the rest, and a barcode is a cell when the signal component is
at least 100,000 times as likely as the noise component at its count. No paper describes the
caller; this page follows 10x's algorithm pages, a reading of the one public release (1.2) and
10x's published outputs, which pin down the other releases. Planned for 0.3.

1. TOC
{:toc}

## At a glance

| Id | `cellranger-atac` (planned) |
|---|---|
| Tool | Cell Ranger ATAC 1.0.0 – 2.2.0 (`cellranger-atac count`) |
| Category | mixture |
| Counts | fragments per barcode that overlap a peak (`peak_region_fragments`), not cut sites |
| Needs expected cells | no (`--force-cells` optional) |
| Deterministic | yes: no step draws random numbers |
| Original licence | 1.2 source: MIT wording plus a clause limiting use to 10x products (not a free licence); binaries under 10x's EULA. Never copied |
| Extra inputs | the fragments file, the peaks and the reference's contig lengths (for the exclusions), or 10x's `singlecell.csv` |

## How it works

For each species $$s$$ of the reference (one, or two in a barnyard) and each barcode $$b$$, over
the species' primary contigs:

1. **Counts.** $$F_s(b)$$ is the number of fragments of $$b$$ (lines of the fragments file, which
   holds deduplicated nuclear fragments with both mates at MAPQ above 30; read-pair counts
   ignored). $$P_s(b)$$ counts those that overlap a peak, once however many they touch: fragment
   $$[x, y)$$ overlaps peak $$[a, e)$$ when $$a < y$$ and $$e \ge x$$. $$D$$ is the number of
   fragments over all barcodes. The peak-barcode matrix counts cut sites: its column sums are
   not $$P$$.
2. **Low-targeting barcodes (1.1+).** $$g_s$$ is the fraction of the species' primary genome
   within 2,000 bp of a peak. A barcode with $$F_s > 0$$ and $$P_s / F_s < g_s$$ is excluded for
   $$s$$.
3. **Gel-bead doublets (1.2+).** Per GEM group, barcodes with more than 250 fragments are
   candidates; above 15,000 of them, the 15,000 alphabetically last barcode sequences are kept,
   not the 15,000 largest as documented. A fragment of $$b_1$$ that starts exactly where a
   fragment of $$b_2$$ ended (a shared Tn5 cut) adds one link between $$b_1$$ and $$b_2$$.
   Positions where the barcodes ending there plus the fragments starting there exceed 20 are
   skipped. Each candidate's neighbour is the other candidate with most links (ties: the
   alphabetically first). Mutual neighbours are a doublet, whatever their link count; the one with
   fewer fragments (equal: the alphabetically later) is excluded for every species. 10x's page
   speaks of more links with each other than with themselves; self-links are never compared.
4. **Barcode multiplets (1.2+).** A barcode is read as parts A (bases 1–7), C (8–9) and B
   (10–16). Per GEM group, links $$L$$ are counted as in step 3 between barcodes that differ only
   in A, and between barcodes that differ only in B. For a major barcode $$M$$ with at least 100
   fragments, a barcode $$m$$ that differs from it only in A, or only in B, is excluded for every
   species when $$L(m, M) > 0.4\,L(m, m)$$. A minor matched by several majors keeps the largest
   link; if $$m$$ and $$M$$ each match the other, the one with more fragments is kept (ties follow
   Python 2 dictionary order).
5. **Merge.** A barcode keeps the first reason in the order multiplet, doublet, low-targeting.
   2.x writes it as `excluded_reason`: 0 none, 1 gel-bead doublet, 2 low-targeting, 3 multiplet.
6. **Shift.** $$c = \max(5, \lfloor 0.02\,D / W \rfloor)$$, with $$W$$ the size of the barcode
   list (737,280): the fragments a barcode expects from other GEMs. The data are
   $$y = P_s - c$$ for the barcodes not excluded for $$s$$ with $$P_s \ge c$$, zeros kept. With
   fewer than 10 distinct values of $$y$$ the species has no cells, forced or not.
7. **Mixture.** $$\mathrm{NB}(k; \mu, \alpha)$$ has mean $$\mu$$ and variance
   $$\mu + \alpha\mu^2$$. The fit runs on the distinct values $$v$$ with their multiplicities.
   - *Start*: $$t_0$$ is the largest $$v$$ whose cumulative share of fragments is at most 0.01,
     clamped between the second-smallest and second-largest values; $$v \le t_0$$ is noise, the
     rest signal.
   - *M-step*: each component's $$\mu$$ is its weighted mean, its $$\alpha$$ the estimate below,
     $$\pi_n$$ the noise share of the weight; if $$\mu_n > \mu_s$$ the components swap.
   - *Dispersion*: with weighted mean $$m$$ and population variance $$s^2$$,
     $$\alpha_0 = \max(10^{-6}, (s^2 - m)/m^2)$$ (all zeros: $$10^{-6}$$), kept when
     $$s^2 \le m$$ or a value exceeds $$10^8$$; otherwise the fixed point of the NB score equation
     at $$\mu = m$$, $$\alpha \leftarrow \log(1 + \alpha m) / (m - \overline{S_\alpha})$$ with
     $$S_\alpha(k) = \sum_{i=1}^{k-1} \alpha i / (1 + \alpha i)$$, to a relative change of
     $$10^{-6}$$ or 10,000 steps.
   - *E-step*: responsibilities in proportion to $$\pi_n \mathrm{NB}(v; \mu_n, \alpha_n)$$ and
     $$(1 - \pi_n) \mathrm{NB}(v; \mu_s, \alpha_s)$$; a value where both underflow drops out.
   - E and M alternate until all five parameters change by a relative $$10^{-6}$$ or less, or 250
     iterations.
8. **Threshold.** Scan $$k$$ upwards from $$k_{min}$$ to $$\lfloor \mu_s \rfloor$$ for the first
   $$k$$ with $$(1 - \pi_n)\,\mathrm{NB}(k; \mu_s, \alpha_s) \ge R\,\pi_n\,\mathrm{NB}(k; \mu_n, \alpha_n)$$,
   $$R = 100{,}000$$; if none passes, $$\lfloor \mu_s \rfloor$$. In 1.2 (and 1.1, whose outputs
   match it), $$k_{min} = \max(\lfloor \mu_n \rfloor, q)$$ with $$q$$ the signal component's 0.001
   quantile; 2.x starts at $$\lfloor \mu_n \rfloor$$. With $$T = k + c$$, the cells are the
   barcodes not excluded with $$P_s \ge T$$, united over species.
9. **Cap and force.** Up to 2.1, more than 20,000 cells for a species forces 20,000. Forcing $$N$$
   (`--force-cells`): $$T$$ is the $$N$$-th largest positive $$y$$ plus $$c$$, so ties at it can
   give more than $$N$$ cells.
10. **Barnyards.** 10x documents a refit on the cell barcodes. In 1.2 the refit's noise
    component is empty and the first threshold stands: the pass changes nothing (the reported
    thresholds of atac_hgmm_1k_nextgem are the first pass's).

The reported `goodness_of_fit` and `estimated_cells_present` do not enter the call.

## Versions

| Cell Ranger ATAC | Source | Change to the call |
|---|---|---|
| 1.0.0, 1.0.1 | none; pinned from 10x's outputs | no exclusions; $$R = 1{,}000$$ (100,000 gives 287 instead of 224 on atac_v1_pbmc_5k) |
| 1.1.0 | none; outputs | low-targeting exclusion (±2,000 bp); $$R = 100{,}000$$; the fit matches the 1.2 rule to every printed digit |
| 1.2.0 | public, one commit | gel-bead doublet and barcode-multiplet exclusions; the scan floored at the signal's 0.001 quantile |
| 2.0.0 | none; outputs | new peak caller and barcode-aware duplicate marking (change $$P$$); primary contigs are those with genes; `excluded_reason`; the scan starts at $$\lfloor \mu_n \rfloor$$, without 1.2's quantile floor |
| 2.1.0 | none; outputs | low-targeting padding ±250 bp; "improved calling on multi-species samples", perhaps the barnyard refit (inferred) |
| 2.2.0 | none; release notes only | 20,000-cell cap removed; `--force-cells` takes any positive integer; fragment strand used by the doublet and multiplet filters |

From `singlecell.csv` alone the call rule is exact on 10x's 1.0.1, 1.1.0, 2.0.0 and 2.1.0
datasets; 1.2.0 needs the fragment-level exclusions. A refit of 2.x data agrees with 10x's fitted
parameters to 4–5 significant digits and gives thresholds 802 (10x: 803, no barcode between) and
304; whether 2.x keeps some other floor is not known. No 2.2 dataset was found.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--force-cells` | unset | call this many cells per species; 1 – 20,000 up to 2.1, any positive integer in 2.2 |
| `reanalyze --barcodes` | unset | take the cells from a given `singlecell.csv` (1.1+) |
| `--peaks` | the pipeline's | custom peaks, which change $$P$$ and $$g_s$$ (`reanalyze` 1.1+, `count` 2.0+) |
| reference | — | primary contigs, species prefixes and contig lengths |
| internal | $$R$$: 100,000 (1,000 in 1.0) | the odds a cell needs (step 8) |
| internal | none (1.0), 2,000 bp (1.1 – 2.0), 250 bp (2.1+) | low-targeting padding (step 2) |
| internal | 0.02 and 5; 250, 15,000; 100, 0.4; 20 | the shift; doublet and multiplet limits; pileup cap |
| internal | 20,000 (1.0 – 2.1) | cell cap per species |

## Status

Planned for 0.3. Implementing it needs:

- inputs the package does not read: per barcode $$P$$, $$F$$ and the exclusions, plus $$D$$ and
  $$W$$. A count matrix is not enough (it holds cut sites per peak). A reader for 10x's
  `singlecell.csv` gives $$P$$, $$F$$ and, from 2.0, `excluded_reason`; a reader for the
  fragments file, the peaks and the contig lengths gives everything from scratch;
- the mixture of steps 6–8, with SciPy's negative-binomial pmf and quantile;
- Variants for 1.0, 1.1, 1.2, 2.0 and 2.1. 2.2 has neither source nor a published dataset, so
  it would be a Gap.

The oracle is 10x's published outputs, as for [Cell Ranger](../cellranger/cellranger.md):
`singlecell.csv` and `summary.json` of atac_v1_pbmc_5k (1.0.1), atac_pbmc_5k_v1 (1.1.0),
atac_pbmc_5k_nextgem and atac_hgmm_1k_nextgem (1.2.0), atac_pbmc_10k_nextgem (2.0.0) and
10k_pbmc_ATACv2_nextgem_Chromium_Controller (2.1.0). The check is exact from `singlecell.csv` for
1.0, 1.1 and 2.x; 1.2's exclusions need its fragments and peaks. Running the original is not
needed (1.2 builds only against 10x's dependency bundle; 2.x is binary only). Blockers to
exactness: the SciPy 10x bundled; the remaining 2.x fit difference (count definition, EM details,
or some 40 low-count barcodes whose low-targeting status differs); the multiplet ties; and 2.x's
exclusion rules, which have no source (use `excluded_reason`).

## References

- 10x Genomics, Cell Ranger ATAC, *Algorithms overview*, one page per release, e.g.
  [1.2](https://support.10xgenomics.com/single-cell-atac/software/pipelines/1.2/algorithms/overview),
  and the [release notes](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/release-notes/release-notes).
- Satpathy et al. 2019, Massively parallel single-cell chromatin landscapes of human immune cell
  development and intratumoral T cell exhaustion, *Nat Biotechnol* 37:925–936,
  [doi:10.1038/s41587-019-0206-z](https://doi.org/10.1038/s41587-019-0206-z): the assay, with its
  own cell rule, not the pipeline's.
- Source (10x licence, read, never copied): [10XGenomics/cellranger-atac](https://github.com/10XGenomics/cellranger-atac/tree/5753003c286de55c193f02e1d66c847cb635aa60)
  at commit `5753003` (1.2), `mro/atac/stages/processing/cell_calling/`.
- Research notes: docs/_research/planned-atac.md §1.
