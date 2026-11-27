# Planned for 0.3: ATAC and multiome cell calling

Source reading for the 0.3 method pages (November 2026): Cell Ranger ATAC, Cell Ranger ARC,
EmptyDropsMultiome, BAP's bead knee and SADE. Each section gives identity and licence, releases,
the algorithm in our own words, parameters, per-barcode output and what an oracle run needs.
Tags as in [README.md](README.md): **[SRC]** read in the source at the tag or commit named,
**[DOC]** documentation or release notes only, **[OUT]** inferred from published outputs,
**UNVERIFIED** inferred with nothing to check it against.

Licence rule: 10x Genomics code (cellranger-atac, cellranger-arc) is under 10x's own licence and
is never copied or transliterated. Its sections describe the algorithm in prose and formulas and
quote at most an identifier or a constant. The method pages are written from these descriptions.

The Cell Ranger ATAC summary in [cellranger-starsolo.md](cellranger-starsolo.md) §4 is the
starting point; corrections to it are marked **Correction**.

Findings that shape the 0.3 pages:

- Cell Ranger ATAC's mixture is pinned down by 10x's own outputs: from `singlecell.csv` and
  `summary.json` the call rule is exact on every release checked, 1.1's fit matches to every
  printed digit, 1.0 used an odds ratio of 1000 (not 100000), and 2.x starts the threshold scan
  at the noise mean without 1.2's floor at the signal component's 0.1 % quantile.
- Cell Ranger ARC's calls are reproduced exactly (three datasets, 1.0 and 2.0) by K = 2 k-means
  on de-duplicated $$(\log A, \log G)$$ points, with $$A$$ = ATAC **cut sites** in peaks (not
  fragments), $$G$$ = exonic + intronic UMIs, and > 1 count in each library.
- EmptyDropsMultiome has two call rules: the README's `FDR_multi < 0.001` (large droplets
  accepted) and the paper's `FDR < 0.001` (not accepted). Its results depend on the DropletUtils
  version (DM before 1.28) and its RNA mixture fit starts from R's unseeded RNG.
- BAP's knee function is unchanged since 2018; its result depends on the per-bead count (bap1
  vs bap2) and on the R version (`density` changed in R 4.4.0).
- SADE has no licence and no release, and its entropy formula changed many times after the
  preprint; the preprint describes k = 2 and a depth normalisation closer to the current code
  than to the commit of its posting date.

## 1. Cell Ranger ATAC: negative-binomial mixture on fragments in peaks

### 1.1 Identity

- **Paper**: none describes the cell caller. Satpathy et al. 2019, *Nat Biotechnol* 37:925–936,
  doi:10.1038/s41587-019-0206-z (checked: resolves to "Massively parallel single-cell chromatin
  landscapes…") introduced the 10x scATAC assay and processed it with Cell Ranger ATAC 1.0, but
  its methods call cells with their own rule (≥ 1,000 unique fragments, TSS enrichment ≥ 8,
  ≤ 45,000 fragments), not the pipeline's mixture [DOC]. The specification is 10x's
  "Algorithms overview" page (one per release) plus the 1.2 source.
- **Repository**: `https://github.com/10XGenomics/cellranger-atac`, one branch (`main`), no tags,
  one commit `5753003c286de55c193f02e1d66c847cb635aa60` (2020-10-06, message "Cell Ranger ATAC
  1.2"). Python 2.7 stages under the Martian pipeline framework, plus Go and C [SRC].
- **Licence**: `LICENSE` at that commit, "Copyright (c) 2020 10x Genomics". MIT wording with an
  extra clause 2: the rights may be exercised only with a 10x Genomics product or data generated
  with one. Redistribution and modification are permitted, but the field-of-use clause makes it
  non-free and incompatible with GPL-3 [SRC]. Bundled third-party code has its own licences
  (`lib/plsa/src/LICENSE`, `tenkit/.../docopt/LICENSE`). **Never copied**: the method page is
  written from the description below.
- **Distribution**: binary tarballs from 10x (EULA); the 1.2 source builds only against 10x's
  "Ranger 1.2.0" dependency bundle (README) [SRC].

**Correction** to `cellranger-starsolo.md` §1: that table gives cellranger-atac the licence of
Cell Ranger 7.2+ ("internal use … no redistribution"). The cellranger-atac 1.2 `LICENSE` is the
MIT text plus the 10x-product clause; the conclusion (never copy) is unchanged.

### 1.2 Releases

From 10x's release notes (Wayback copy of
`https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/release-notes/release-notes`,
2026-07-11) and the per-release "Algorithms overview" pages
(`https://support.10xgenomics.com/single-cell-atac/software/pipelines/<v>/algorithms/overview`,
Wayback 2025). The 10x site itself now answers automated requests with a bot checkpoint, so
these notes use the archived copies.

| Release | Date | Cell-calling change | Source | Tag |
|---|---|---|---|---|
| 1.0.0 | 2018-10-29 | Initial: NB mixture on fragments in peaks, odds ratio **1000**, no exclusions, < 20k cells per species, `--force-cells` | No | [DOC] 1.0 page; [OUT] 1.0.1 |
| 1.0.1 | 2018-12-17 | none listed | No | [DOC] |
| 1.1.0 | 2019-04-17 | Low-targeting exclusion added (peaks padded ±2000 bp); odds ratio **100000** | No | [DOC] 1.1 page; [OUT] |
| 1.2.0 | 2019-11-21 | Gel-bead doublet and barcode-multiplet exclusions added | **Yes** | [SRC] |
| 2.0.0 | 2021-05-03 | New wavelet peak caller and barcode-aware duplicate marking (both change the counts); primary contigs = gene-containing contigs; `excluded_reason` and `nonprimary` columns; threshold search no longer floored at the signal component's 0.1 % quantile | No | [DOC]; threshold rule [OUT] |
| 2.1.0 | 2022-04-04 | "Less stringent" low-targeting filter (peaks padded ±**250** bp); "improved calling on multi-species samples" | No | [DOC]; 250 bp [OUT] (`fraction_of_genome_within_250bp_of_peaks`) |
| (2.1) | 2025-02-20 | References 2024-A only | — | [DOC] |
| 2.2.0 | 2025-04-30 | 20,000-cell cap removed; `--force-cells` takes any positive integer; strand in fragments "improves … doublet and inclusion list contamination detection" | No | [DOC] |

Only 1.2 has public source. 1.0, 1.1 and 2.x are Gaps for code, but published outputs pin down
most of their behaviour (§1.6).

### 1.3 The algorithm (1.2) [SRC]

Stages under `mro/atac/stages/processing/cell_calling/` at commit `5753003`:
`remove_low_targeting_barcodes`, `remove_gel_bead_doublet_barcodes`, `remove_barcode_multiplets`,
`merge_excluded_barcodes`, `detect_cell_barcodes`; EM helpers in `lib/python/analysis/em_fitting.py`;
constants in `lib/python/constants.py`; the pipeline in `mro/atac/_produce_cell_barcodes.mro`.
No step uses random numbers.

**Inputs.** The fragments file (one line per deduplicated fragment: contig, Tn5-adjusted start
and end, corrected barcode, read-pair count), the peaks BED, the reference's contig definitions
(primary contigs including sex chromosomes; species prefixes for barnyard references, `""` for a
single species) and contig lengths, and the barcode whitelist (`737K-cratac-v1`, 737,280 barcodes).
The fragments file holds only fragments with both mates MAPQ > 30, not mitochondrial, not chimeric,
on a primary contig, with a barcode (`mark_duplicates`).

**Step 0: counts.** For each species s and barcode b, over the primary contigs of s:

- $$F_s(b)$$ = number of fragment lines of b (the read-pair count is ignored);
- $$P_s(b)$$ = number of those fragments whose span overlaps at least one peak. A fragment
  $$[x, y)$$ overlaps a peak $$[a, e)$$ when $$a < y$$ and $$e \ge x$$: a peak ending exactly at
  the fragment's start counts, one base more permissive than half-open overlap. Peaks are taken
  as sorted and non-overlapping. A fragment counts once however many peaks it touches.
- $$D = \sum_{b,s} F_s(b)$$, the total number of fragments.

These are the `passed_filters` and `peak_region_fragments` columns of `singlecell.csv` (barnyard:
`passed_filters_<sp>`, `peak_region_fragments_<sp>`) [SRC], [OUT]. They are fragments, not cut
sites: the peak-barcode matrix counts cut sites, so its column sums are not $$P$$ [DOC].

**Step 1: low-targeting barcodes (1.1+).** For species s, let $$g_s$$ = (bases covered by the
union of all peaks each padded by 2,000 bp on both sides, clipped to the contig) / (total length
of the species' primary contigs). Exclude b for s when $$F_s(b) > 0$$ and
$$P_s(b)/F_s(b) < g_s$$. Reported as `fraction_of_genome_within_2000bp_of_peaks` and
`number_of_low_targeting_barcodes`.

**Step 2: gel-bead doublets (1.2).** Per GEM group, the candidates are barcodes with more than
250 fragments (all species together). If there are more than 15,000 the code keeps the 15,000
**alphabetically last** barcode strings, not the 15,000 largest (the comment says "top N"; the
sort is over the names): a finding, and a difference from the docs. A linkage count
$$L(b_1, b_2)$$ is built per primary contig: walking fragments by start position, a fragment of
$$b_1$$ starting at $$x$$ is
linked to every earlier candidate fragment that ended exactly at $$x$$ (a shared Tn5 event); each
such earlier fragment of $$b_2$$ adds to $$L(b_1,b_2)$$ and $$L(b_2,b_1)$$ (once if
$$b_1 = b_2$$). A position is skipped when (distinct barcodes with a fragment ending at this
fragment's end) + (fragments starting at $$x$$) exceeds 20. Each candidate's nearest neighbour is
the other candidate with the largest $$L$$ (ties: alphabetically first). Mutual nearest
neighbours form a putative doublet, with no minimum link count; the member with fewer fragments
is excluded for every species (equal counts: the alphabetically later one). 10x's page describes
the test as "more linked fragments with each other than with themselves"; the code never compares
with the self-links.

**Step 3: barcode multiplets (1.2).** A barcode splits into part A (bases 1–7), part C (8–9) and
part B (10–16). Per GEM group, two linkage tables like step 2 (adding 1 per linked pair, same
pileup cap 20): between barcodes that share C and B (differ in A), and between barcodes that share
C and A (differ in B). For every "major" barcode M with at least 100 fragments and every other
observed part-A (resp. part-B) sequence with the same C, the "minor" barcode m that differs from M
only there is excluded when $$L(m, M) > 0.4 \cdot L(m, m)$$. A minor matched by several majors
keeps the one with the largest link (first seen on ties). Then, if m → M and M → m both hold, the
barcode with more fragments is reinstated (ties depend on Python 2 dict order); a one-step chain
collapse follows, which does not change the excluded set. Excluded for every species, label
`whitelist_contam`.

**Merging.** Exclusion sets are merged in the order multiplets, gel-bead doublets, low-targeting;
a barcode keeps the first reason. 2.x writes it as `excluded_reason` (0 none, 1 gel-bead doublet,
2 low-targeting, 3 barcode multiplet) [DOC]; on 2.0 and 2.1 outputs, barcodes that are
low-targeting by the step 1 rule but labelled 1 or 3 (16,441 on pbmc_10k_nextgem 2.0) show the
same precedence [OUT].

**Step 4: shift.** $$c = \max\left(5, \left\lfloor 0.02\,D / W \right\rfloor\right)$$ with
$$W$$ = whitelist size (constants `MINIMUM_COUNT = 5`, `WHITELIST_CONTAM_RATE = 0.02`), the
expected fragments per barcode from other GEMs. For each species, the data are
$$y_b = P_s(b) - c$$ for the barcodes not excluded for s with $$P_s(b) \ge c$$ (zeros kept).
If the $$y$$ take fewer than 10 distinct values, the species gets no cells (and `--force-cells`
is ignored for it).

**Step 5: two-component NB mixture by EM.** Write the distinct values $$v_j$$ with multiplicities
$$w_j$$, and $$\mathrm{NB}(k;\mu,\alpha)$$ for the negative binomial with mean $$\mu$$ and
variance $$\mu + \alpha\mu^2$$ (scipy `nbinom` with $$n = 1/\alpha$$, $$p = 1/(1+\alpha\mu)$$).

- *Start.* $$t_0$$ = the largest $$v_j$$ whose cumulative share of fragments
  $$\sum_{i \le j} v_i w_i / \sum_i v_i w_i$$ is at most 0.01 (if none, the second-smallest
  value), clamped between the second-smallest and second-largest distinct values. Hard
  responsibilities: noise $$r_n(v) = [v \le t_0]$$, signal $$r_s(v) = [v > t_0]$$.
- *M-step*, per component, over the values with $$r(v_j) > 0$$ and weights $$u_j = w_j r(v_j)$$:
  $$\mu = \sum u_j v_j / \sum u_j$$; $$\alpha$$ from the dispersion estimate below;
  $$\pi_n = \sum_j w_j r_n(v_j) / \sum_j w_j (r_n(v_j) + r_s(v_j))$$. If $$\mu_n > \mu_s$$ the
  components are swapped ($$\pi_n \to 1 - \pi_n$$).
- *Dispersion.* If every value is 0: $$10^{-6}$$. Weighted mean $$m$$ and population variance
  $$s^2$$; $$\alpha_0 = \max(10^{-6}, (s^2 - m)/m^2)$$. Return $$\alpha_0$$ if the largest value
  exceeds $$10^8$$ or $$s^2 \le m$$. Otherwise iterate
  $$\alpha \leftarrow \log(1 + \alpha m) \big/ \left(m - \overline{S_\alpha}\right)$$, where
  $$S_\alpha(k) = \sum_{i=1}^{k-1} \alpha i / (1 + \alpha i)$$ ($$S_\alpha(0) = 0$$) and the bar is
  the $$u$$-weighted mean, until $$2|\alpha' - \alpha| / (\alpha' + \alpha) < 10^{-6}$$ (return
  the new value) or 10,000 iterations (return the last). This is the fixed-point form of the NB
  score equation for $$\alpha$$ at $$\mu = m$$; the code credits RiboPy.
- *E-step.* $$a_j = \pi_n \mathrm{NB}(v_j;\mu_n,\alpha_n)$$,
  $$b_j = (1-\pi_n)\mathrm{NB}(v_j;\mu_s,\alpha_s)$$; responsibilities $$a/(a+b)$$, $$b/(a+b)$$;
  where $$a + b = 0$$ (both pmfs underflow) both are set to 0, so the value drops out of the next
  M-step.
- Iterate E then M until every one of the five parameters changes by a relative amount
  $$|\theta_{old} - \theta|/\theta \le 10^{-6}$$, at most 250 iterations.

**Step 6: threshold.** $$k_{min} = \max(\lfloor \mu_n \rfloor, q)$$, $$q$$ = the 0.001 quantile
of the fitted signal NB. Scan $$k = k_{min}, k_{min}+1, \dots, \lfloor \mu_s \rfloor$$ for the
first $$k$$ with
$$\dfrac{(1-\pi_n)\,\mathrm{NB}(k;\mu_s,\alpha_s)}{\pi_n\,\mathrm{NB}(k;\mu_n,\alpha_n)} \ge R$$,
$$R = 100000$$ (`CELL_CALLING_THRESHOLD`); if none, $$\lfloor \mu_s \rfloor$$. A zero noise pmf
gives an infinite ratio, which passes. If $$k_{min} > \lfloor \mu_s \rfloor$$ the original fails
(no value assigned). The cell threshold is $$T = k + c$$; b is a cell of species s when it is not
excluded for s and $$P_s(b) \ge T$$. The cell barcodes are the union over species.

**Step 7: cap and force.** If more than 20,000 barcodes pass, the species is forced to 20,000.
Forcing N cells (from `--force-cells`, a per-species map or `"default"`, capped at 20,000): among
$$y > 0$$, $$T$$ = the N-th largest $$y$$ plus $$c$$ (the smallest positive $$y$$ plus $$c$$ if
N is at least their number); cells are $$P_s \ge T$$, so ties at the N-th value are all called
and the count can exceed N. N ≤ 0 is an error.

**Step 8: barnyard second pass.** Documented as a refit on the cell barcodes only. In 1.2 the
masking condition is inverted: every non-excluded barcode is passed in with a start threshold of
−1, the noise component is empty, its mean undefined (NaN), the refit's noise fraction NaN, and
the first-pass threshold is kept. The pass changes nothing [SRC]; on atac_hgmm_1k_nextgem (1.2)
the reported `cell_threshold_hg19` 1016 and `cell_threshold_mm10` 1758 equal the first-pass
thresholds recomputed from the reported fitted parameters [OUT]. Skipped when `--force-cells` is
given. 2.1's "improved calling on multi-species samples" may be this fix (UNVERIFIED).

`goodness_of_fit` (the KS statistic of $$y$$ against the fitted mixture CDF) and
`estimated_cells_present` $$= \lfloor (1-\pi_n)\,|y| \rfloor$$ are reported, not used.

### 1.4 Parameters that affect the call

| Parameter | Default | Releases |
|---|---|---|
| `count --force-cells N` | unset; 1 ≤ N ≤ 20,000 (wrapper check) until 2.1; any positive integer in 2.2 | all |
| `reanalyze --force-cells N`, `reanalyze --barcodes singlecell.csv` (use given cells; exclusive with force) | unset | 1.1+ |
| `--peaks` (custom peaks; changes $$P$$ and $$g$$) | pipeline peaks | `reanalyze` 1.1+, `count` 2.0+ |
| Reference (primary contigs, species prefixes, contig lengths) | — | all |
| Internal: odds ratio $$R$$ | 1000 (1.0), 100000 (1.1+) | — |
| Internal: low-targeting padding | none (1.0), 2000 bp (1.1–2.0), 250 bp (2.1+) | — |
| Internal: contamination rate 0.02, minimum shift 5, cap 20,000 (until 2.1), doublet minimum 250 fragments and 15,000 candidates, multiplet minimum 100 fragments and factor 0.4, pileup 20 | — | 1.2 [SRC] |

### 1.5 Output per barcode

`singlecell.csv`: `is_<sp>_cell_barcode` (single species `is__cell_barcode`, the prefix being
empty), `cell_id` (1.x, `<sp>_cell_<i>` numbered over the sorted cell barcodes, else `None`),
`passed_filters`, `peak_region_fragments` (and per-species columns for barnyards), 2.x
`excluded_reason`. `summary.json`: `cell_threshold[_<sp>]`, `fitted_mean_noise`,
`fitted_dispersion_noise`, `fitted_mean_signal`, `fitted_dispersion_signal`,
`fraction_cell_calling_noise` (1.1+), `goodness_of_fit`, `estimated_cells_present`,
`annotated_cells`, `cells_detected`, `number_of_low_targeting_barcodes`,
`fraction_of_genome_within_{2000,250}bp_of_peaks`, `putative_gelbead_doublets_found`,
`putative_barcode_multiplets_found` [SRC], [OUT]. Raw and filtered peak-barcode matrices
(cut sites per peak).

### 1.6 Checks against 10x's published outputs [OUT]

Run on `singlecell.csv` and `summary.json` from 10x's CDN with a scratch implementation written
from §1.3 (not kept):

| Dataset (release) | Call rule | Low-targeting from `g` | Refit of $$P$$ | Threshold |
|---|---|---|---|---|
| atac_v1_pbmc_5k (1.0.1) | cells = $$P \ge 224$$ exactly, no exclusions | — | parameters to 7 digits | 224 only with $$R = 1000$$ (100000 gives 287) |
| atac_pbmc_5k_v1 (1.1.0) | cells = not low-targeting and $$P \ge 146$$ exactly | 187,855 = reported | every printed digit equal | 146 |
| atac_pbmc_5k_nextgem (1.2.0) | all cells $$P \ge 1441$$; 155 non-cells above it are not low-targeting (doublets or multiplets) | 192,273 = reported | needs the fragment-level exclusions | 1441 = signal 0.1 % quantile 1436 + 5 (the floor binds; the odds rule alone gives 387) |
| atac_pbmc_10k_nextgem (2.0.0) | cells = (`excluded_reason` 0 and $$P \ge 803$$) exactly | 227,308 vs 227,257 reported | 4–5 significant digits (μ_s 11,342.9 vs 11,344.8) | 803 = first $$k \ge \lfloor\mu_n\rfloor$$ with odds ≥ 100000 from the reported parameters; the 1.2 floor would give 1125 |
| 10k_pbmc_ATACv2_nextgem_Chromium_Controller (2.1.0) | same rule, $$P \ge 304$$, exact | 102,823 vs 102,782 | 4–5 digits | 304 likewise (floor: 1842); shift 6 ($$D$$ = 257.6 M) |

So: 1.0 used $$R = 1000$$; 1.1 matches the 1.2 code apart from the 1.2 exclusions; 2.x keeps the
model and $$R$$ but drops the signal-quantile floor, starting the scan at $$\lfloor\mu_n\rfloor$$
(the exact form of any remaining floor is UNVERIFIED). A refit of 2.x data agrees with 10x to 4–5
digits and gives thresholds 802 (vs 803; no barcode lies between, so identical calls) and 304; the
remaining difference (count definition, EM details, or the ~40 low-count barcodes whose
low-targeting status differs) is UNVERIFIED. The 2.x low-targeting discrepancies are barcodes
with < 40 fragments, mostly $$P = 0$$; use `excluded_reason` rather than recomputing.

### 1.7 Oracle

- Data on `https://cf.10xgenomics.com/samples/cell-atac/`, each with `_singlecell.csv`,
  `_summary.json`, `_fragments.tsv.gz`, `_peaks.bed`, `_raw_peak_bc_matrix.h5`,
  `_filtered_peak_bc_matrix.h5` (all six answer HTTP 200 for the first six below):
  `1.0.1/atac_v1_pbmc_5k/atac_v1_pbmc_5k_*`, `1.1.0/atac_pbmc_5k_v1/atac_pbmc_5k_v1_*`,
  `1.2.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_*`,
  `1.2.0/atac_hgmm_1k_nextgem/atac_hgmm_1k_nextgem_*` (barnyard),
  `2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_*`,
  `2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_*`;
  also `1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_*` and
  `2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X_*` (only
  `singlecell.csv` checked). No 2.2 dataset found.
- Route A (exact, cheap): from `singlecell.csv` alone ($$P$$, $$F$$, `excluded_reason`; for 1.1
  the low-targeting set from the reported genome fraction) refit and compare threshold, fitted
  parameters and cell set. Works for 1.0, 1.1, 2.0, 2.1. 1.2 needs Route B for the exclusions.
- Route B (full): from `fragments.tsv.gz` + `peaks.bed` + the reference's contig lengths
  (primary contigs, e.g. hg19 1.2 or GRCh38 2020-A) recompute steps 0–3. Only 1.2 has source for
  steps 2–3; 2.x exclusions can be checked against `excluded_reason` (2.2 changed them, per the
  release notes).
- Running the original is not needed: 10x's published calls are the oracle, as for Cell Ranger.
  The 1.2 source needs 10x's Ranger bundle (EULA) and Python 2.7; 2.x is binary-only.
- Blocks to exactness: scipy's `nbinom` pmf and quantile in the bundled (unknown) scipy version;
  the 2.x count and fit differences above; Python 2 dict order in the multiplet ties.
- What sc-droplets needs as input: per-barcode $$P$$ (fragments overlapping peaks), $$F$$, the
  exclusions, $$D$$ and $$W$$. A count matrix is not enough (it holds cut sites per peak).
  Readers for `singlecell.csv` (and ARC's `per_barcode_metrics.csv`) or a fragments file are new
  input work.

## 2. Cell Ranger ARC: joint ATAC + GEX calling (de-duplicated points, OrdMag start, k-means K=2)

### 2.1 Identity

- **Paper**: none. Specification: 10x's "Calling Cell Barcodes" page, one per release
  (`https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/{1.0,2.0}/algorithms/cell-calling`,
  Wayback 2025-09) [DOC]. It cites OrdMag to Zheng et al. 2017 (doi:10.1038/ncomms14049).
- **Repository**: none for the pipeline. The only public code is the command-line wrapper inside
  Cell Ranger: `https://github.com/10XGenomics/cellranger/tree/cellranger-10.1.0/lib/rust/cr_wrap/src/arc/`
  (`count.rs`, `reanalyze.rs`, `types.rs`; Rust), which shows the options and how they reach the
  pipeline, not the caller [SRC].
- **Licence**: that tree's `LICENSE` (tag `cellranger-10.1.0`, "Copyright (c) 2026 10x
  Genomics"): use, copy and modify for internal use only, only with 10x products or data from
  them, no redistribution or sublicensing, derivatives bound by the same terms [SRC]. The ARC
  binaries are under 10x's EULA. **Never copied.**
- **Language and distribution**: binary tarball (Martian pipeline; ATAC processing rewritten in
  Rust in 2.0 [DOC]).

### 2.2 Releases

From the Wayback copy (2026-08-13) of
`https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/release-notes/release-notes` [DOC]:

| Release | Date | Change that affects the call |
|---|---|---|
| 1.0.0 | 2020-09-09 | Initial joint caller. Override `--min-atac-count X --min-gex-count Y` calls barcodes meeting **either** threshold |
| 1.0.1 | 2020-12-10 | none (reference contig limits) |
| 2.0.0 | 2021-05-03 | Wavelet peak caller and barcode-aware duplicate marking (change ATAC counts); barcode correction at Hamming distance 1 (was 2); override now requires **both** ($$A \ge X$$ and $$G \ge Y$$); raw matrix keeps only barcodes with a count. `gex_umis_count` in `per_barcode_metrics.csv` is all zeros (output bug) |
| 2.0.1 | 2022-02-03 | `gex_umis_count` fixed |
| 2.0.2 | 2022-08-18 | none (mkfastq); reference note 2025-02-20 |
| 2.1.0 | 2025-10-06 | "The ATAC filtering thresholds have been adjusted to avoid removing real cells" (may change calls substantially on low-quality ATAC); 20,000-cell cap removed; fragment strand used by the doublet and inclusion-list contamination filters |
| 2.2.0 | 2026-04-20 | GEM-X Epi Multiome chemistry; no cell-calling change listed |

No release has public source; every release is a Gap for code. What 2.1 changed in the ATAC
filters is not documented beyond that sentence (UNVERIFIED; by analogy with Cell Ranger ATAC 2.1
it may be the 250 bp low-targeting padding). The 2.1 algorithm page was not archived in readable
form.

### 2.3 The algorithm [DOC], with what the outputs show [OUT]

Per barcode: $$A$$ = ATAC **transposition events (cut sites) in peaks**
(`atac_peak_region_cutsites`, two per fragment lying inside a peak) and $$G$$ = GEX UMIs
(`gex_umis_count` = exonic + intronic UMIs; introns are counted by default). The ATAC-only
pipeline counts fragments; ARC counts cut sites (the docs say "transposition events in peaks", and
the [OUT] test below separates only on cut sites).

1. **Exclusions** (any one masks the barcode): ATAC low-targeting (fraction of fragments in peaks
   below the fraction of the genome in peaks padded by 2,000 bp; as §1.3 step 1), ATAC gel-bead
   doublets (minor barcode of a linked pair; as §1.3 step 2), and **more than 1 count in each
   library** ($$A > 1$$ and $$G > 1$$). The first two are recorded in `excluded_reason`
   (1 doublet, 2 low-targeting, 3 multiplet); the third is not recorded. ARC 1.0 and 2.0 outputs
   contain no reason 3 (three datasets) [OUT]: no barcode-multiplet filter there.
2. **De-duplication**: the eligible barcodes become points $$(A, G)$$; identical points are
   collapsed to one, so the many low-count background barcodes weigh no more than one point each.
3. **OrdMag start**: on the de-duplicated points, a threshold per axis "10-fold less than the
   maximum value after removing outliers" (Cell Ranger's OrdMag); points above both thresholds
   start as cells, the rest as non-cells. N, quantile and any bootstrap are not documented.
4. **k-means, K = 2**, from the two groups' centroids, one initialisation, on the de-duplicated
   points in **log space** [DOC for k-means; OUT for the space]: the final cells are the points
   nearer the cell centroid in $$(\log A, \log G)$$, Euclidean.
5. **Map back**: each eligible barcode takes the label of its point.
6. **Cap**: < 20,000 cells (mechanism undocumented); removed in 2.1.

Override: with `--min-atac-count X --min-gex-count Y` the caller is skipped; eligible barcodes
(the filters of step 1 still apply) with $$A \ge X$$ or $$G \ge Y$$ (1.0), and (2.0+), are cells.
The wrapper refuses one without the other and passes them as a `"default"` map entry
[SRC wrapper].

No randomness is documented; k-means has one deterministic start. If OrdMag kept Cell Ranger 2's
bootstrap it would be seeded, and the outcome below does not depend on the start.

**[OUT] reproduction** (scratch code, `per_barcode_metrics.csv` only): eligible = `excluded_reason`
0 and $$A \ge 2$$ and $$G \ge 2$$; de-duplicate; $$\log_{10}$$; Lloyd iterations to a fixed
labelling from any of three starts (per-axis 99th percentile ÷ 10 over all points, over the top
3,000, or max ÷ 10). Every start converges in 3–4 iterations to **exactly** 10x's cell set:

| Dataset (release) | Eligible points | 10x cells | Mismatches |
|---|---|---|---|
| pbmc_granulocyte_sorted_10k (1.0.0) | 25,420 | 11,909 | 0 |
| pbmc_granulocyte_sorted_10k (2.0.0) | 29,162 | 11,898 (11,895 points) | 0 |
| pbmc_granulocyte_sorted_3k (2.0.0) | 6,520 | 2,711 | 0 |

The specification is that sharp: $$\log(1+x)$$ instead of $$\log x$$ gives 2–6 mismatches; a
minimum of 1 count gives 12, a minimum of 3 gives 6, keeping the excluded barcodes gives 107;
with fragments in peaks or all fragments instead of cut sites, cells and non-cells are not
separable by any line in log or linear space. The exact OrdMag rule could not be pinned down
because it does not change the fixed point on these data (it could on data with two local
optima).

### 2.4 Parameters

`--min-atac-count`, `--min-gex-count` (unset), `--peaks` (custom peaks change $$A$$ and the
low-targeting fraction), `--gex-exclude-introns` (2.0+; presumably $$G$$ = exonic UMIs only,
UNVERIFIED), `reanalyze --barcodes` (given cells). Internal: at least 2 counts per library,
K = 2, the 20,000 cap until 2.0.2.

### 2.5 Output per barcode

`per_barcode_metrics.csv`: `is_cell`, `excluded_reason`, `atac_peak_region_cutsites`,
`atac_peak_region_fragments`, `atac_fragments`, `gex_umis_count` (zeros in 2.0.0: use
`gex_exonic_umis + gex_intronic_umis`), plus QC columns; `summary.csv` "Estimated number of
cells"; the filtered feature-barcode matrix [DOC], [OUT].

### 2.6 Oracle

- Data (on `https://cf.10xgenomics.com/samples/cell-arc/`; `_per_barcode_metrics.csv` checked for
  each; `_summary.csv`, `_atac_fragments.tsv.gz`, `_atac_peaks.bed`, `_raw_feature_bc_matrix.h5`
  and `_filtered_feature_bc_matrix.h5` checked for both pbmc_granulocyte_sorted_10k releases,
  the others by name pattern, UNVERIFIED):
  `1.0.0/pbmc_granulocyte_sorted_10k/`, `2.0.0/pbmc_granulocyte_sorted_10k/`,
  `1.0.0/pbmc_granulocyte_sorted_3k/`, `2.0.0/pbmc_granulocyte_sorted_3k/`,
  `1.0.0/e18_mouse_brain_fresh_5k/`, `2.0.0/e18_mouse_brain_fresh_5k/`,
  `2.0.0/lymph_node_lymphoma_14k/`, `2.0.0/human_brain_3k/`. No 2.1 or 2.2 dataset found.
- Route: exact comparison from `per_barcode_metrics.csv` alone ($$A$$, $$G$$, `excluded_reason`),
  as above. Recomputing the exclusions needs the fragments (§1.3 steps 1–2, with ARC's ATAC
  barcode list, believed to be `737K-arc-v1`, UNVERIFIED).
- Blocks: the OrdMag start (undocumented; matters only if k-means has several fixed points); the
  2.1 filter change and the cap mechanism; no source.
- A count matrix holds $$G$$ per gene and ATAC cut sites per peak; if the ATAC rows of
  `raw_feature_bc_matrix.h5` sum to $$A$$ (UNVERIFIED), the method can run from the raw matrix
  given the exclusions.

**Correction** to `cellranger-starsolo.md` §4: "OrdMag N, log scaling and k-means metric
undocumented" can be narrowed. 10x's own calls are reproduced exactly with log scaling of the raw
counts and Euclidean K = 2 k-means on unweighted de-duplicated points, ATAC = cut sites in peaks,
GEX = exonic + intronic UMIs, eligibility > 1 in each library [OUT]. ARC 1.0/2.0 has no
barcode-multiplet exclusion [OUT].

## 3. EmptyDropsMultiome (Megas, Lorenzi, Marioni 2024)

### 3.1 Identity

- **Paper**: Megas S, Lorenzi V, Marioni JC. EmptyDropsMultiome discriminates real cells from
  background in single-cell multiomics assays. *Genome Biol* 2024;25:121.
  doi:10.1186/s13059-024-03259-x (PubMed 38741206, PMC11520057; checked) [DOC].
- **Repository**: `https://github.com/MarioniLab/EmptyDropsMultiome`. One tag, `first_release`
  = commit `a1b1bbbd55a19c4306ef6e8c438a25e979b3cedd` (2024-02-21; GitHub release 2024-04-23;
  Zenodo doi:10.5281/zenodo.11051741). `main` = `5fd8e2f154` (2024-06-03) differs from the tag
  only in `README.md` and a figure [SRC]. Reproducibility code:
  `https://github.com/MarioniLab/EmptyDropsMultiome2023_final` (tag `first_release_reprod_repo`,
  Zenodo record 11051756), which vendors mixtools 2.0.0 with one added `print` of the
  log-likelihood (no arithmetic change) [SRC].
- **Licence**: `DESCRIPTION` "GPL (>= 3) + file LICENSE"; `LICENSE` (YEAR 2022, "EmptyDropsMultiome
  authors"), `LICENSE.md` (GPL-3 text). MIT until commit `cd3de95` (2023-05-24, "change license
  from MIT to GPL3") [SRC]. May be ported closely with its notice.
- **Language and distribution**: R package `EmptyDropsMultiome` 0.0.0.9000, installed with
  `devtools::install_github`; not on CRAN, Bioconductor or conda (checked 2026-11-27). Imports
  DropletUtils, mixtools, S4Vectors, Matrix, ggplot2, stats.

### 3.2 Releases

One release. Development history before it, for reference [SRC]: Fisher's method then the
arithmetic mean for combining p-values (`d2f5d15` 2023-05-21, `4bd1b8d` 2023-05-24); the
higher line moved from 1/2 to 2/3 of the way to the tip (`1a46dc3`, 2023-07-03); `seed`
argument added (`4b798b7`, 2023-08-01). Pre-release states are not worth a Variant.

The result depends on the installed **DropletUtils** (from our notes: `emptyDrops`' default is a
Dirichlet-multinomial with estimated α in 1.2–1.26 and the plain multinomial from 1.28;
EmptyDropsMultiome passes no `alpha`) and on **mixtools** (`normalmixEM`). The README installs
R 4.3.0 with the matching Bioconductor DropletUtils (1.20 or 1.22, DM) and describes the model as
Dirichlet-multinomial, so the paper's results are the DM ones (UNVERIFIED which minor version).
A method page should treat DropletUtils < 1.28 and ≥ 1.28 as two Variants of the dependency.

### 3.3 The algorithm (`first_release`) [SRC]

Inputs: a gene × barcode RNA count matrix and a peak × barcode ATAC matrix with barcode names.
The README and the paper's scripts take both from Cell Ranger ARC's `raw_feature_bc_matrix.h5`
(`Read10X_h5`, assays "Gene Expression" and "Peaks"), so the ATAC counts are cut sites per peak.
Main function `emptydrops_multiome()` (`R/emptydrops_multiome.R`).

**Step 1: per-modality bounds** (`R/fit_3_normals.R`), unless the user gives both `lower` and
`barhop` for that modality. On the raw totals of every barcode in the matrix (not logged), fit a
three-component normal mixture with mixtools `normalmixEM` (`epsilon = 1e-9`, `maxit = 1000`,
ECM off):

| | RNA (`fit_3_normals`) | ATAC (`fit_3_normals_atac`) |
|---|---|---|
| Start means | 5, 500, 1200 | 5, 100, 10000 |
| Start SDs | 5, 50, 800 | 5, 50, 800 |
| Start weights | **random** (mixtools draws $$U(0,1)^3$$, normalised, redrawn until all ≥ 0.05) | 6/7, 0.8/7, 0.1/7 (mixtools renormalises them to sum 1) |
| Sanity check (else defaults) | $$\mu_1 > 100$$ or $$\mu_2 > 700$$ or $$\mu_1 > \mu_2 - 15$$ | $$\mu_1 > 60$$ or $$\mu_2 > 600$$ or $$\mu_1 > \mu_2 - 15$$ |
| Defaults (also on fit error) | μ (25, 200, 1000), σ (1, 60, 300), λ (1, 0.10, 0.05)/1.15 | μ (1, 40, 800), σ (1, 60, 600), λ (1, 0.10, 0.05) |
| `lower` $$l$$ | $$\mu_2 + 1.5\sigma_2$$ | $$\mu_2 + 2\sigma_2$$ |

Components are ordered by mean (1 = barcode-swap "barhop" noise, 2 = ambient, 3 = cells).
`barhop` $$b$$ is where the weighted densities of components 1 and 2 are equal: the root
$$x = (-B + \sqrt{B^2 - 4AC})/(2A)$$ of $$Ax^2 + Bx + C = 0$$ with
$$A = 1/\sigma_1^2 - 1/\sigma_2^2$$, $$B = -2(\mu_1/\sigma_1^2 - \mu_2/\sigma_2^2)$$,
$$C = \mu_1^2/\sigma_1^2 - \mu_2^2/\sigma_2^2 + 2\log(\lambda_2\sigma_1/(\lambda_1\sigma_2))$$;
$$-C/B$$ if $$A = 0$$. If the discriminant is ≤ 0 nothing is returned and the later steps fail.
mixtools restarts from fully random values when a variance collapses (< 1e-8) and errors after
20 restarts (EmptyDropsMultiome then uses the defaults). EmptyDropsMultiome does not seed the RNG
before this step, so the RNA fit depends on the caller's RNG state; the paper's script calls
`set.seed(0)` first.

**Step 2: emptyDrops per modality** (`R/emptydrops_with_barhop.R`): drop barcodes with total
$$< b$$; `set.seed(seed)` (default 42); DropletUtils `emptyDrops(m, lower = l, niters,
test.ambient = TRUE)` with every other argument at its default; `niters` 10,000 for RNA and
25,000 for ATAC. So the ambient profile comes from barcodes with $$b \le$$ total $$\le l$$, and
p-values are computed for every remaining barcode with a positive total, ambient ones included.
Only the `PValue` column is used, so emptyDrops' `retain` (the knee) and its own FDR do not
affect the call. Dropped barcodes get NA.

**Step 3: combine** (`R/aggr_p_values.R`): join the two tables by barcode (union, so a barcode
missing from one matrix gets NA), $$p = (p_{RNA} + p_{ATAC})/2$$ (arithmetic mean; NA if either is
NA; the code's Fisher option is commented out), and BH over all non-NA $$p$$ (ambient barcodes
included). The paper cites Vovk & Wang 2020 for the mean [DOC].

**Step 4: k-means line** (`R/cellranger_k_means_v2.R`, `call_cells`; the authors' own
approximation of Cell Ranger ARC). Totals NA → 0. Sort by (ATAC, RNA); keep one barcode per
distinct (ATAC, RNA) pair with both > 0 (ARC requires > 1 and removes the ATAC exclusions; this
does neither). Per axis, on the kept raw totals $$x$$: $$q_1$$ = the 0.99 quantile (R type 7);
$$S = \{x > q_1/10\}$$ (if empty, threshold 0); $$q_2$$ = the 0.99 quantile of $$x_S$$;
threshold = 0.99 quantile of $$\{x : x < q_2\}$$ (over all kept values, not just $$S$$) ÷ 10,
or $$q_1/10$$ if no value is below $$q_2$$.
Initial cells: ATAC and RNA both above their thresholds. If that is everyone **or no one**, every
kept point is called a cell and the line is (0, 0) (a bug for "no one"). Otherwise transform to
$$\log_{10}(x + 0.1)$$, start centres at the means of the two initial groups, and run R's
`stats::kmeans` (Hartigan–Wong, `iter.max = 10`, one start; `set.seed(123)` is called but
unused). The cell cluster is the one with the larger RNA centre. Labels pass to duplicates; a
barcode with a zero in either modality is not a k-means cell. The k-means line is the
perpendicular bisector of the two final centres in $$(x, y) = (\log_{10}(A+0.1),
\log_{10}(R+0.1))$$: $$y = \beta x + \alpha$$. Hartigan–Wong's final partition need not equal the
nearest-centre split, so the `k_means` column and the line can disagree near the line.

**Step 5: two parallel lines** (`R/calculate_perp_bisector.R`, `R/find_tip.R`,
`R/accept_k_means.R`):

- *Tip*: $$s = \log_{10}(R+0.1) - \beta\log_{10}(A+0.1)$$ for every barcode; the tip is
  $$s >$$ its 0.99 quantile; its centre is (mean $$R$$, mean $$A$$) on the raw scale, then
  $$\log_{10}(\cdot + 0.1)$$.
- *Higher line*: intercept $$\alpha_{hi} = \alpha + \tfrac{2}{3} d \sqrt{\beta^2+1}$$, where
  $$d$$ is meant to be the distance from the tip centre to the k-means line. **Bug**: the point is
  passed as (log RNA, log ATAC) while the line is in (log ATAC, log RNA), so
  $$d = |\beta \log_{10}(\bar R + 0.1) - \log_{10}(\bar A + 0.1) + \alpha| / \sqrt{\beta^2+1}$$,
  the distance of the mirrored point.
- *Lower line*: through $$(\log_{10}(l_A + 0.1), \log_{10}(l_R + 0.1))$$ with slope $$\beta$$;
  intercept capped at $$\alpha$$.

**Step 6: final FDR.**

1. `FDR`: the BH value of step 3, set to 1 for barhop barcodes of either modality, for
   barcodes with $$\log_{10}(R+0.1) < \beta\log_{10}(A+0.1) + \alpha_{lo}$$, and where NA.
2. `FDR_multi`: `FDR`, but 0 for barcodes with
   $$\log_{10}(R+0.1) - \beta\log_{10}(A+0.1) > \alpha_{hi}$$ (`above_ambiguous`).

Call: the README counts `FDR_multi < 0.001`. The paper says the above-the-line acceptance is
optional and was **not** used ("FDR ≤ 0.001 without setting FDR to 0 for very large droplets"),
and its script counts `FDR < 0.001` [DOC], [SRC reproducibility repository,
`real/tutorial_w_mixture_fit_and_umap_qc.R`]. The package has no switch: both columns are
returned, and the two rules are two call modes.

### 3.4 Parameters

`lower_rna`, `barhop_rna`, `lower_atac`, `barhop_atac` (NULL: estimated, per modality, unless
both of the pair are given), `niter_rna = 10000`, `niter_atac = 25000`, `seed = 42`, `verbose =
TRUE` (prints and draws plots; needs a graphics device). Not exposed: FDR 0.001, the mixture
starts, sanity bounds and defaults, 1.5σ / 2σ, the pseudocount 0.1, the tip quantile 0.99, the
factor 2/3, `iter.max = 10`.

### 3.5 Output per barcode

An `S4Vectors::DataFrame`, one row per barcode in either matrix: `Total_RNA`, `LogProb_RNA`,
`PValue_RNA`, `Limited_RNA`, `FDR_RNA`, `barhop_amb_tent_RNA` (2 barhop, 1 ambient, 3 above
`lower`), the same for `_chromatin`, `PValue_multi`, `FDR_multi`, `k_means`, `above_ambiguous`,
`FDR`; metadata `lower_*`, `barhop_*`, `k_means_slope`, `k_means_intercept`,
`higher_intercept`, `lower_intercept` [SRC].

### 3.6 Oracle

- Original: R in a pixi environment (conda-forge R, bioconda DropletUtils, conda-forge
  r-mixtools), EmptyDropsMultiome installed from commit `a1b1bbbd` with `R CMD INSTALL`, a
  `pdf()` device for the plots, `set.seed(0)` before the call (as the paper). Two environments:
  DropletUtils 1.22 (DM, the paper's) and ≥ 1.28 (multinomial). Compare `lower_*`, `barhop_*`,
  every p-value, `FDR`, `FDR_multi`, `k_means`, the three intercepts and both call sets.
- Data: Cell Ranger ARC `raw_feature_bc_matrix.h5` of pbmc_granulocyte_sorted_10k (ARC 2.0.0, the
  paper's PBMC sample) and pbmc_granulocyte_sorted_3k; the paper's gonad samples are
  Garcia-Alonso et al. 2022 (not checked).
- Exactness: the Monte Carlo is DropletUtils' (already reproduced in this project, seeded by
  `set.seed(42)`); the RNA mixture start uses R's `runif` after the caller's seed (reproducible
  with `rrandom`); `normalmixEM` and R's Hartigan–Wong `kmeans` (Fortran AS 136, GPL via R) need
  ports. With those, an exact match is plausible; the validation rule for Monte Carlo methods
  applies otherwise.

## 4. BAP: the bead knee (Lareau et al. 2019)

### 4.1 Identity

- **Paper**: Lareau CA, Duarte FM, Chew JG, et al. Droplet-based combinatorial indexing for
  massive-scale single-cell chromatin accessibility. *Nat Biotechnol* 2019;37:916–924.
  doi:10.1038/s41587-019-0147-6 (PubMed 31235917, PMC10299900; checked). Methods: "highly
  abundant barcodes were detected … by quantifying each unique barcode sequence among
  nuclear-mapping reads. Our knee-calling algorithm then established a per-experiment bead
  threshold" [DOC].
- **Repository**: `https://github.com/caleblareau/bap` (the paper cites
  `github.com/buenrostrolab/bap`, a fork last pushed 2018-10-30). Tags `bap-barcode`
  (2018-06-21), `biorad-paper1` (`7872481`, 2018-08-06, the paper's code), `v0.6.0` (`8b35816`,
  2019-07-26). `master` = `763fb13` (2026-07-06, packaging modernised, version 0.7.1) [SRC].
- **Licence**: `LICENSE`, MIT, "Copyright (c) 2018 Caleb Lareau" [SRC]. May be ported with the
  notice.
- **Language and distribution**: Python (click CLIs, Snakemake) calling R scripts; PyPI
  `bap-atac` 0.6.4 (2019-10-29) to 0.7.0 (2021-08-29); 0.7.1 only in the repository.

### 4.2 Releases

The knee function `get_density_threshold_CL` (`bap/bin/R/00_knee_CL.R`) and its caller
`bap/bin/R/10b_knee_execute.R` last changed in `f9d30a4` (2018-07-26) and are byte-identical at
`biorad-paper1`, `v0.6.0` and `master` [SRC]. What changed is the per-bead count fed to it:

| Entry point | Releases | Count per bead (input to the knee) |
|---|---|---|
| `bap` (bap1; `bap/cli.py` at `biorad-paper1`, later `cli_old_dontEdit.py`; `bin/python/10_quantBarcode_Filt.py`) | 2018 – today | reads with MAPQ **> 30** (`--mapq`) on chromosomes of at least 1,000,000 bp (`--min_chromosome_size`; the nuclear ones), counted once per (read start, bead) [SRC] |
| `bap2` (`bap/cli_bap2.py`, `bin/snake/Snakefile.bap2.chr`, `bin/R/21_annotate_fragments.R`) | 0.6.0 (2019-07) – today | distinct fragments (chromosome, Tn5-shifted start +4/−5, end) built from proper pairs (default; `--all-pairs` off), insert ≤ 2000 (`--max-insert`), BEDPE score ≥ 30 (`--mapq`; bedtools' score for the pair, taken to be the lower mate MAPQ, UNVERIFIED), not overlapping the blacklist; summed over every chromosome split **including the mitochondrial one** (the Snakefile appends the mito sample to the list it sums) [SRC] |

NEWS 0.6.0: "Alpha version of bap2 that can be deployed for all 10x and BioRad scATAC-seq
datasets". The two are one method with two count definitions; a page can take the count as input
and document both.

### 4.3 The algorithm [SRC, `00_knee_CL.R` at `f9d30a4`]

Input: one count $$c_b$$ per bead (the file's second column, `V2`), `logTransform = TRUE`,
`type = "bead"`.

1. **Floor**: $$f = 0.001 \times \mathrm{mode}(c)$$ (the most frequent value; ties go to the
   value seen first in file order, which for bap2 is awk's hash order). Keep $$c > f$$. With a
   mode below 1000 this keeps every bead with a count.
2. $$z = \log_{10} c$$.
3. **Density**: R `density(z, bw = 0.1, kernel = "gaussian", n = 10000, from = min z, to = max z)`:
   an absolute bandwidth of 0.1 in log10 units, R's FFT-binned estimate on its internal grid
   interpolated to 10,000 points $$x_1 < \dots < x_{10000}$$ with values $$y_i$$. R 4.4.0
   changed this function's grid coordinates ("correcting density values typically by a factor
   of about 0.999"; `old.coords = TRUE` restores the old ones) [DOC R NEWS], so the result can
   depend on the R version.
4. **Local minima**: mark each step "down" when $$y_i - y_{i-1} < 0$$ (with $$y_0 = +\infty$$,
   so the first step is down; flat steps count as "up"); take the run-length ends of the
   down/up sequence; if there are more than 3 runs keep every other end starting with the first
   (the ends of descending runs); if 3 or fewer keep all ends, so a local maximum (the end of a
   rising run) can then be a candidate; if $$y_1 = y_2$$ drop the first. Drop candidates at
   index 1 or 10,000.
5. **Choice**: scanning candidates from the right (highest count first), take the first index
   $$i$$ with $$i \ge 2000$$ (20 % of the grid; the comment says "20% of barcodes") and
   ($$\max z - x_i > 0.5$$ or $$x_i < \max z / 2$$).
6. **Threshold** $$T = 10^{x_i}$$; if no candidate qualifies, $$T = 0$$.
7. **Guard** ("bead"): if $$T > 100000$$ or $$T < 100$$, use 500. Both are written
   (`bead_threshold`, `bead_threshold_nosafety`).
8. **Call**: beads with $$c_b \ge$$ the guarded threshold, bead name not `NA`.

No random numbers. Overrides: `--minimum-barcode-fragments N` (N > 0 skips the knee and keeps
$$c_b \ge N$$), `--barcode-whitelist FILE` (the listed beads).

Against UMI-tools' density method (our `umitools.py`, notes B2), which this follows: the floor is
0.001 × the **mode**, not the maximum; the bandwidth is absolute (0.1 log10 units, R's FFT
`density`), not scipy's `gaussian_kde` factor 0.1 × the SD; minima come from run lengths, not
`argrelextrema`; beads are $$\ge T$$, not $$> T$$; the 100–100,000 guard falls back to 500; there
is no expected-cells mode. Notes B10 are right; to add: the 500 fallback applies to any
threshold outside [100, 100000] (including the 0 of a failed search), and the bead counts differ
between `bap` and `bap2`.

After this step (not part of the bead call): pairwise Jaccard of Tn5 insertions among the called
beads, a second knee with the same function (`type = "jaccard"`, on the top 10^6 Jaccard values,
guard outside [1e-6, 0.5] → 0.005) unless `--minimum-jaccard-index`, and merging of beads into
droplet barcodes.

### 4.4 Parameters

`--minimum-barcode-fragments` (0 = knee), `--barcode-whitelist`, `--mapq` (30), `--max-insert`
(2000, bap2), `--all-pairs` (off), `--blacklist-file` / `--reference-genome`, `--bead-tag` (`XB`).
Constants: 0.001, bw 0.1, n 10,000, 0.2 × n, 0.5, guard [100, 100000] → 500.

### 4.5 Output per bead

`knee/<name>.barcodeQuantSimple.csv` (bead, count), `knee/<name>.bapParams.csv`
(`bead_threshold`, `bead_threshold_nosafety`), `final/<name>.HQbeads.tsv` (called beads); then
`final/<name>.barcodeTranslate.tsv` (bead → droplet) [SRC].

### 4.6 Oracle

- Function level, as the UMI-tools oracle: run `Rscript bap/bin/R/10b_knee_execute.R counts.csv 1
  V2` from tag `v0.6.0` (identical to `master`) in a pixi environment (conda-forge R,
  r-data.table, r-dplyr) on count vectors, in two R versions (4.3 and ≥ 4.4) because of the
  `density` change. Inputs: per-barcode `passed_filters` from 10x ATAC `singlecell.csv`
  (unique fragments, a stand-in for bap2's count) and per-barcode totals from RNA matrices.
- Pipeline level would need a BAM with bead tags (the paper's dscATAC-seq data, GEO GSE123581,
  accession from the paper; contents not checked); not needed, since the step is a function of
  the count vector.
- Exactness: port R's `density.default` (linear binning, FFT, `approx`) in both coordinate
  conventions; everything else is integer or exact comparisons.

## 5. SADE: Shannon entropy for scATAC-seq denoising (2026 preprint)

### 5.1 Identity

- **Preprint**: Yang S, An K, Lee JJ, Gehring J, Booeshaghi AS, Huang H, Streets AM. Denoising
  single-cell ATAC-seq data with Shannon entropy. *bioRxiv* 2026, posted 2026-10-23, v1, CC-BY.
  doi:10.64898/2026.10.20.760589 (resolves to bioRxiv; Crossref and the bioRxiv API agree;
  Europe PMC PPR1350460). The full text was read once through a summarising fetch (bioRxiv
  rate-limited every other request; its equations are images), so its statements below are
  [DOC] at second hand.
- **Repository**: `https://github.com/streetslab/SADE` (the preprint's code availability; it
  also names `https://github.com/streetslab/SADE_manuscript`, analysis notebooks and run
  scripts, no licence, last pushed 2026-10-21, which call a local clone without a commit).
  SADE: created 2025-06-10, no tags or releases, 359 commits, branches `main`, `norm`,
  `Np1div2`, `feature/macs-peak-caller`; very active, 92 commits in November 2026 alone [SRC].
- **Licence**: **none**. No `LICENSE` file and no licence in the GitHub metadata, so all rights
  are reserved by default. Treat like the 10x code: read, describe in our own words, never copy.
  `filter_bam_with_barcodes.sh` adapts a 10x knowledge-base snippet for BAM filtering
  (irrelevant here).
- **Language and distribution**: Bash driver (`sade.sh`) calling Python 3 scripts (numpy,
  pandas, scipy, matplotlib, tqdm; unpinned, `EnvironmentSetup/pip_requirements.txt`) and
  ripgrep. Not on PyPI or conda; run from a clone.

### 5.2 Versions

No releases, so the method page must name commits. Two states matter [SRC]:

| State | Commit | Date | Differences |
|---|---|---|---|
| Preprint | `eef42ee` (last default-branch commit before the 2026-10-20 submission) | 2026-09-30 | no blacklist; right insertion at the fragment end; adjustment $$\lfloor (n+1)/2 \rfloor$$; entropy without depth terms (below); spline $$s = 1$$; $$k = 2$$; no fallbacks (no candidate → crash); debris flag "$$P_{closed} < 1 - 0.5$$" |
| Current | `f83ca20` | 2026-11-24 | blacklist masking (default on); right insertion at end − 1; $$\lfloor (n+2)/2 \rfloor$$; depth-adjusted entropy; $$k = 3$$ in `sade.sh`; shoulder and single-population fallbacks; `-n` targeted cells; `-p` threads; `autosomes` option; debris flag on the capped open-state fraction |

Between them the entropy normalisation alone changed about ten times (commits of 2026-11-18 to
2026-11-21: "Adjust read-depth for state entropy", "Use (insertion +2)/2", "Scale mixture state
entropy…", "Come back to the og normalization", "Modify normalization"). Any implementation
should pin a commit and expect it to go stale.

### 5.3 The algorithm at `f83ca20` [SRC]

Input: a 10x-style fragments file (5 tab-separated columns: contig, start, end, barcode, count;
rows with another number of fields are dropped), the genome's contig sizes (`ref/<genome>/`),
the ENCODE blacklist for that genome, chromosome(s) (default `chr1` only), window size $$w$$
(default 3000 bp).

**Step 1: insertion counts per window.** For each fragment on the chosen chromosome, two Tn5
insertions: at `start` and at `end − 1`; each falls in window $$\lfloor pos / w \rfloor$$ of
$$M = \lfloor L / w \rfloor + 1$$ windows ($$L$$ = contig length). An insertion in a window that
overlaps a blacklist interval is set aside; the blacklisted windows still count in $$M$$. The
count column is ignored (one fragment line = two insertions). Several chromosomes are
concatenated into one vector.

**Step 2: entropy per barcode.** Let $$n_j > 0$$ be the non-empty windows' counts, $$M_+$$ their
number, $$a_j = \lfloor (n_j + 2)/2 \rfloor$$ (described as a paired-end adjustment), $$f_v$$ the
number of windows with $$a_j = v$$, and $$\bar a = \sum_v v f_v / M_+$$.

- If $$\bar a \le 1$$ the barcode gets no entropy and is dropped.
- Zero-truncated Poisson MLE: $$\lambda = W_0(-\bar a e^{-\bar a}) + \bar a$$ (principal
  Lambert W); $$P_0 = e^{-\lambda}$$.
- Expected empty windows inside the open region: $$M_0 = M_+ P_0 / (1 - P_0)$$.
- Open-state probability, depth-adjusted: $$p = (M_0 + M_+) / (M \bar a)$$; if $$p > 1$$ the
  barcode is dropped. Reported `P_open_state` $$= \min(p\,\bar a, 1)$$.
- $$H_{states} = -p\log_2 p - (1-p)\log_2(1-p)$$.
- With $$\pi_v = \mathrm{Pois}(v;\lambda)(1 - P_0)$$:
  $$H_{open} = -\sum_v f_v \pi_v \log_2 \pi_v - M_0 P_0 \log_2 P_0$$.
- $$H = H_{states} + p \cdot H_{open} (1 - P_0) / (M_+ \bar a)$$ (column `Entropy`).

Written per barcode: `Entropy`, `Entropy_open_region`, `Mp` ($$M_+$$), `mle_lambda`,
`P0_open_region`, `P_open_state`.

**Step 3: cutoff on the ranked curve** (`find_threshold.py`, `autothreshold.py`). Sort entropies
in decreasing order, $$y_r = \log_{10} H_{(r)}$$, $$r = 0, 1, \dots$$; keep the first 20,000.
Fit a cubic smoothing spline with SciPy `make_splrep(x=r, y, k=3, s=1)` (SciPy ≥ 1.15). Compute
the first derivative $$d_1$$ and the second derivative $$d_2$$ scaled by $$\max|d_2|$$ at every
rank.

- Candidates, scanning ranks upwards until $$k$$ are found ($$k$$ = `-k`, default 3 in `sade.sh`;
  the README and the tutorial's table say 2): a region opens at the first rank with $$d_2 < 0$$
  and closes at the next rank with $$d_2 > 0$$; its candidate is the rank of minimal $$d_1$$
  (steepest fall) within it.
- $$k = 1$$: the candidate is the "shoulder", the rank where the spline (both axes scaled to
  [0, 1]) lies farthest above the chord from the first to the last point. For $$k > 1$$ the
  shoulder is used only if a region opened and never closed and no region was found.
- No candidate: median and MAD of $$y$$; the cutoff is the number of ranks with
  $$y - \mathrm{median} > -3\,\mathrm{MAD}$$, capped at the last rank.
- Among the candidates, the cutoff rank $$R$$ is the one with the smallest $$d_1$$.
- With `-n N` (targeted cells, default 0 = off): if $$R \notin [\lfloor N/5 \rfloor, 20000]$$, take
  the next-steepest candidates in turn until one falls inside (or none is left).

**Step 4: call.** The $$R$$ barcodes with the highest entropy (`sort -g -r` on the column, then
the first $$R$$; the barcode at rank $$R$$, whose value is reported as the cutoff, is excluded).
Of these, barcodes with `P_open_state` > `-s` (default 0.5; values outside [0.1, 1] revert to 0.5)
are flagged as debris and removed. Output: `entropy_filtered_bc_df.tsv`, `bc_pass_entropy.tsv`,
`Entropy_filtered_bc_CBZ.txt`, `filtered_fragments.tsv` (ripgrep fixed-word match of the
barcodes), `entropy_cutoff.csv` (cutoff, log10 cutoff, rank), `all_candidate_entropy_cutoffs.csv`.

No randomness. Multiprocessing (`-p`, default 4) changes only the row order of the table; the
ties at rank $$R$$ are broken by `sort`'s whole-line comparison. The result depends on SciPy's
`make_splrep` (new in 1.15), on the blacklist file shipped in the repository and on the
chromosome choice.

### 5.4 The preprint-date state `eef42ee` [SRC]

Same pipeline with: insertions at `start` and `end` (no − 1; an `end` equal to the contig length
would raise), no blacklist, $$a_j = \lfloor (n_j + 1)/2 \rfloor$$,
$$p = (M_0 + M_+)/M$$ (barcode dropped if $$p > 1$$), $$\pi_v = (f_v/M_+)(1 - P_0)$$ (the
empirical zero-truncated frequencies, not a Poisson pmf),
$$H_{open} = -\sum_v \pi_v\log_2\pi_v - P_0\log_2 P_0$$,
$$H = H_{states} + p\,H_{open}(1-P_0)/M_+$$; spline $$s = 1$$, $$k = 2$$, candidates as above
but without the shoulder or MAD fallbacks; debris when $$1 - p < 1 - 0.5$$ (i.e. $$p > 0.5$$).

**What the preprint says** [DOC, via the fetch summary]: a per-barcode mixture of a Poisson
for accessible windows and a point mass at zero for inaccessible ones; the Poisson rate from a
zero-truncated MLE through the principal Lambert W branch; imputed zero windows $$M_0$$; an
entropy decomposition in which the conditional entropy of the accessible regions is normalised
"by the total number of observed Tn5 insertions"; invalid barcodes get no entropy; a cubic
spline on the ranked entropies, candidates from the second derivative, **k = 2** by default, the
candidate with the lowest first derivative chosen; debris = accessible genome fraction > 0.5;
defaults 3 kb windows and chr1; no paired-end adjustment mentioned; no versions. Normalising by
the total insertions matches the current code's $$M_+\bar a$$ (the sum of the adjusted counts)
better than `eef42ee`'s $$M_+$$, so the preprint may describe a later state than the
default-branch commit of its posting date, or a branch (UNVERIFIED). Data: the multiome PBMC
"MO V2" of De Rop et al. (ATAC SRR24422953, RNA SRR24446480) and 10x's PBMC 10k ATAC v1.1, human
brain 3k, fresh and cryopreserved E18 mouse brain 5k.

### 5.5 Parameters

`-c` chromosome(s) (`chr1`; a space-separated list, or `autosomes`), `-w` window (3000), `-k`
candidates (3), `-s` saturation cutoff (0.5), `-b` blacklist Y/N (Y), `-n` targeted cells (0),
`-p` threads (4), `-g` genome (hg38, mm10, hg19 references shipped). Fixed: 20,000-rank limit,
spline degree 3 and $$s = 1$$, 3 MAD.

### 5.6 Oracle

- Original: a clone at a pinned commit, a Python venv with the SciPy of the day, `-p 1`.
- Data: the repository's tutorial runs `f83ca20` (or a commit shortly before) on 10x
  `https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz`
  with `-g hg38 -k 3`: 54,046 scored barcodes, candidates at ranks 157, 534, 995, cutoff rank
  534, entropy cutoff 0.01763 [OUT, the notebook's saved output]. Any 10x ATAC or ARC
  `fragments.tsv.gz` (§1.7, §2.6) works.
- Exactness: deterministic; needs exact ports of the Lambert W step (SciPy's `lambertw`),
  Poisson pmf and `make_splrep`. The input is a fragments file, not a count matrix, so this
  method needs the same fragment reader as Cell Ranger ATAC's route B.
- Status advice: "planned" until the authors tag a release and add a licence; the page can
  describe the method meanwhile.

## Summary

| Tool | Category | Licence | Public source? | Deterministic? | Inputs | Oracle route |
|---|---|---|---|---|---|---|
| Cell Ranger ATAC 1.0–2.2 | Two-component NB mixture (EM) on fragments in peaks, odds-ratio threshold; fragment-level exclusions first | 10x: 1.2 `LICENSE` = MIT + 10x-product clause; binaries EULA. Never copied | 1.2 only (one commit `5753003`); 1.0, 1.1, 2.x none | Yes (no RNG; Python 2 dict order only in multiplet ties) | Per barcode: fragments overlapping peaks, total fragments, exclusions (need fragments + peaks + contig lengths); total fragments; barcode-list size | 10x `singlecell.csv` + `summary.json`: rule exact on 1.0.1, 1.1, 2.0, 2.1; fit exact on 1.1; 1.0 needs odds 1000; 2.x drops the quantile floor; 1.2 exclusions need the fragments |
| Cell Ranger ARC 1.0–2.2 | 2-D: de-duplicated (ATAC, GEX) points, OrdMag start, K = 2 k-means in log space | 10x EULA; wrapper under Cell Ranger's 10x licence. Never copied | No (CLI wrapper only) | Yes | Per barcode: ATAC cut sites in peaks, GEX UMIs (exonic + intronic), ATAC exclusions | `per_barcode_metrics.csv`: calls reproduced exactly on 3 datasets (1.0.0, 2.0.0); 2.1+ unchecked |
| EmptyDropsMultiome (`first_release`, `a1b1bbbd`) | Two emptyDrops Monte Carlo tests, mean of p-values, BH; k-means line with two parallel lines | GPL ≥ 3 | Yes | Seeded (`set.seed(42)` per modality); RNA mixture start random unless the caller seeds | RNA gene × barcode and ATAC peak × barcode matrices (ARC raw h5) | Run the R original in pixi, DropletUtils 1.22 (DM) and ≥ 1.28 (multinomial), `set.seed(0)`, ARC PBMC 10k 2.0.0 |
| BAP bead knee (`biorad-paper1` … 0.7.1) | KDE valley on log10 per-bead counts (R `density`, bw 0.1), guard → 500 | MIT | Yes (knee unchanged since 2018-07-26) | Yes; R ≥ 4.4 `density` grid change | Per-bead counts: bap1 unique read starts (nuclear, MAPQ > 30), bap2 unique fragments (filters, blacklist, incl. mito) | Function level: `Rscript 10b_knee_execute.R` from `v0.6.0` on count vectors, R 4.3 and ≥ 4.4 |
| SADE (`eef42ee` preprint date; `f83ca20` current) | Per-barcode Shannon entropy of Tn5 insertions in genome windows (zero-truncated Poisson), spline inflection on the ranked log entropy, debris flag | None (all rights reserved). Describe only | Yes, no tags; changes weekly | Yes | Fragments file, contig sizes, blacklist; chr1 by default | Run a pinned commit; the repository's tutorial (atac_pbmc_500_nextgem, cutoff rank 534) |

Inputs beyond a count matrix: Cell Ranger ATAC needs fragments in peaks per barcode (from
`singlecell.csv` or a fragments file), ARC needs `per_barcode_metrics.csv` (or the raw h5 plus
the ATAC exclusions), BAP needs per-bead counts, SADE needs the fragments file. Only
EmptyDropsMultiome runs from the raw matrices alone.
