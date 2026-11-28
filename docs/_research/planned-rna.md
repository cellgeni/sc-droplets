# Planned RNA methods (0.2)

Source-reading notes, November 2026, for the methods planned for release 0.2 on RNA: CB2 /
scCB2, DIEM, DropletQC, DropSift, valiDrops, SiftCell, miQC, dropEst's quality classifier and
salmon alevin's naive-Bayes whitelist. They go deeper than section C of
[other-methods-and-reviews.md](other-methods-and-reviews.md); where it was wrong, the
corrections are listed before the summary table at the end. Tags as in the
[README](README.md): **[SRC]** read in the source at the tag or commit named, **[DOC]**
documentation or paper only, **UNVERIFIED** inferred. Source was read from clones of each
repository (Bioconductor git for scCB2 and miQC) and, for dependencies whose behaviour matters
(ks, hdrcde, segmented, mixtools, flexmix, e1071, Seurat, xgboost, pastecs, pcaPP, squat), from
their CRAN mirrors or repositories at the versions named.

Each section: identity, releases, the algorithm, parameters, output, and what an oracle run
would need.

## 1. CB2 / scCB2

### Identity

- Paper: Ni Z, Chen S, Brown J, Kendziorski C. CB2 improves power of cell detection in
  droplet-based single-cell RNA sequencing data. *Genome Biology* 21:137 (2020),
  doi:10.1186/s13059-020-02054-8 (PMID 32513247; resolves to this paper).
- Code: Bioconductor `scCB2` (<https://git.bioconductor.org/packages/scCB2>), mirrored at
  <https://github.com/zijianni/scCB2>. R only (the DESCRIPTION says `NeedsCompilation: yes` but
  there is no `src/`). Imports DropletUtils (`emptyDrops`), edgeR (`goodTuringProportions`),
  doParallel/foreach, Seurat.
- Licence: `License: GPL-3` in DESCRIPTION; no LICENSE file in the repository. [SRC]
- Read at Bioconductor branch `RELEASE_3_24` (scCB2 1.24.0, commit 69d9ccd), `R/CB2FindCell.R`,
  `R/FilterGB.R`, `R/GetCellMat.R`, `R/QuickCB2.R`, `R/CheckBackgroundCutoff.R`. [SRC]

### Releases

| scCB2 | Bioconductor | Change to the call |
|---|---|---|
| 0.99.x (GitHub, 2019 – Aug 2020) | — | the paper's code; 0.99.29 changed the sparse correlation, 0.99.30 added `TopNGene`, 0.99.31 the barnyard filter (after the paper) |
| 1.0.0 | 3.12 (Oct 2020) | first release; `CB2FindCell` identical from here to 1.24.0 |
| 1.2.0, 1.4.0 | 3.13, 3.14 | none |
| 1.6.0 | 3.15 (Apr 2022, commit a42167c) | `GetCellMat(MTfilter=)` default 1 → **0.25** (cells with > 25 % `MT-`/`mt-` UMIs dropped); new `CheckBackgroundCutoff`. `CB2FindCell` unchanged; `QuickCB2` still passes `MTfilter = 1` |
| 1.8.0 – 1.24.0 | 3.16 – 3.24 | none (version bumps only) |

`git diff RELEASE_3_15 RELEASE_3_24 -- R` is empty. [SRC] bioconda has
`bioconductor-sccb2` 1.0.0 – 1.20.0. A pipeline today runs 1.20.0 – 1.24.0; for the call they
are all one variant. The call does depend on the installed **DropletUtils** (its `emptyDrops`
decides every barcode outside the clusters), so the variant is really (scCB2 any, DropletUtils
version). [SRC]

### Algorithm (`CB2FindCell`, 1.0.0 – 1.24.0) [SRC]

Input: a feature × barcode UMI matrix (`dgCMatrix` or dense), raw.

1. **Features.** With `GeneExpressionOnly = TRUE`, rows whose name matches
   `TotalSeqB|Cell Multiplexing Oligo` are dropped (`Read10xRaw` names rows by symbol by default,
   so this only works when the antibody names carry that text). `FilterGB(x, 0, 0)` then drops
   barcodes with total 0 and genes with total 0.
2. **Upper threshold** (when `upper = NULL`). `Calc_upper(dat, lower)` is a modified
   `barcodeRanks`: run-length encode the sorted totals, mid-rank per run
   $$r = \text{cumsum}(len) - (len - 1)/2$$; fit `smooth.spline(log10 r, log10 total, df = 20)` on
   **all** runs (not only those above `lower`, unlike DropletUtils), keep the runs with total >
   `lower`; on the fitted $$y$$: $$d_1 = \Delta y / \Delta x$$, right edge = argmin $$d_1$$, left
   edge = argmax $$d_1$$ over `1:right`; curvature $$y''/(1 + y'^2)^{3/2}$$ from the spline's
   derivatives on [left, right]; knee = `round(10^y[argmin curvature])`, inflection =
   `round(10^y[right])` (fitted values, R's half-to-even `round`). Then iterate: `step =
   (knee - inflection)/10` from the **first** fit only; repeat { `upper_temp = knee`; `new_lower
   = inflection + step` (inflection of the latest fit); stop if at most 3 barcodes have total >
   `new_lower`; refit with `lower = new_lower`; stop if the knee is unchanged }. `upper =
   upper_temp`. No iteration cap. `upper ≤ lower` is an error.
3. **Retained cells**: barcodes with total ≥ `upper` (totals over the gene rows). Genes now zero in
   the rest are dropped again.
4. **Background** $$B_0$$: total ≤ `lower`; candidates $$B_1$$: `lower` < total < `upper`. Ambient
   $$p$$ = `edgeR::goodTuringProportions(rowSums(B0))`.
5. **Barnyard**: if **every** gene name contains `_`, candidates are kept only when one species
   prefix holds ≥ 0.75 of their UMIs; the others are never tested (not cells).
6. **Top genes**: keep genes with `rank(-p) ≤ TopNGene` (30000; `rank` averages ties, so the
   unseen genes, all tied at Good-Turing's $$p_0/N_0$$, are kept or dropped together).
7. **Cluster-threshold profile** $$c$$: $$c = p$$, unless ≥ 10 barcodes were retained and the
   entropy of $$p$$ is ≤ the entropy of the retained cells' Good-Turing profile, in which case
   $$c$$ = that profile. Entropy $$-\sum p\log p$$ with zeros set to 1.
8. **Statistic** per candidate: Pearson correlation of its count vector with $$p$$.
9. **Clusters.** Sort candidates by total, ascending (`sort`, radix, stable: ties in column order);
   chunks of 1000 consecutive barcodes, the last merged into the one before if < 500. In a chunk
   of $$m > 20$$ barcodes: Pearson correlation between barcodes (population covariance via
   `crossprod`), `hclust(as.dist(1 - cor))` (complete linkage), `cutree(k = m/20)`
   (`as.integer`, so ⌊m/20⌋ clusters). Each cluster gets: statistic = median of its members'
   correlations with $$p$$; size label = the total of the member whose correlation is nearest
   that median (first on ties); tightness = mean off-diagonal within-cluster correlation (0 for a
   singleton). Chunks of ≤ 20 barcodes give no clusters.
10. **Cluster filter.** For $$b = 1..10$$: $$\tau_b$$ = mean off-diagonal correlation of
    `rmultinom(100, 100 b, c)` (master R RNG). Index $$= \lfloor \text{size}/100 \rfloor$$, set to
    9 if > 10 and to 2 if < 3 (so sizes 1000–1099 use $$\tau_{10}$$ but ≥ 1100 use $$\tau_9$$).
    Keep clusters with tightness > $$\tau_{\text{index}}$$; the rest go back to step 12 as single
    barcodes.
11. **Cluster test.** For each distinct size label $$s$$: 1000 draws `rmultinom(1000, s, p)`,
    statistic of each draw; $$P = (\#\{\text{sim} \le \text{obs}\} + 1)/1001$$ (low correlation =
    cell-like). BH over all kept clusters; clusters with adjusted $$P \le$$ `FDR_threshold` are
    cells, all their members.
12. **Single barcodes.** `DropletUtils::emptyDrops(GE rows, columns = remaining candidates + B0 +
    retained, lower = lower, retain = upper)` with every other argument at its default (10,000
    iterations, multinomial). Cells: `FDR ≤ FDR_threshold` (NA → not a cell); the retained
    barcodes come back here with FDR 0.
13. Cells = emptyDrops cells (incl. retained) ∪ cluster cells. If no single candidates remain,
    retained ∪ cluster cells.

Randomness: steps 10 and 12 use the master R RNG (`set.seed` reproduces them). Step 11 runs in
`foreach %dopar%` on a PSOCK cluster made inside the function (`makeCluster(Ncores)`), so the
simulations use the workers' own unseeded RNG: **not reproducible**, whatever the caller seeds,
even with `Ncores = 1`. P-values have resolution 1/1001. [SRC]

### Parameters

`CB2FindCell(RawDat, FDR_threshold = 0.01, lower = 100, upper = NULL, GeneExpressionOnly =
TRUE, Ncores = 2, TopNGene = 30000, verbose = TRUE)`. `GetCellMat(CBout, MTfilter = 0.25,
MTgene = NULL)` (1.6.0+; 1 before) drops cells with mito fraction > `MTfilter` (genes matching
`\<MT-` or `\<mt-`). `QuickCB2(dir | h5file, FDR_threshold = 0.01, MTfilter = 1, ...)`.
`CheckBackgroundCutoff` (1.6.0+) suggests the smallest of 100, 150, …, 500 that puts ≥ 90 % of
barcodes or ≥ 10 % of UMIs at or below it (advice only). [SRC]

### Output

A `SummarizedExperiment`: assay `cell_matrix` (all rows, cell columns: emptyDrops cells then
cluster cells); metadata `ClusterStat` (per kept cluster: `corr`, `count` = size label, `pval`,
`padj`), `Cluster` (member barcodes per kept cluster, named by tightness), `BarcodeStat` (per
single candidate: `logLH` = emptyDrops LogProb, `count`, `pval`, `padj`), `background` (ambient
counts). Classes per barcode: retained / cluster cell / emptyDrops cell / not a cell. [SRC]

### Oracle

R ≥ 4.0 with Bioconductor scCB2 (bioconda `bioconductor-sccb2` 1.20.0, or Bioc 3.24's 1.24.0)
and a pinned DropletUtils. The knee/upper and the cluster structure are deterministic (and
`hclust`'s tie handling must be matched exactly); the cluster test is not reproducible as
shipped. Options: compare within a noise band (the repo's Monte Carlo rule), or make the workers
deterministic for the oracle only, e.g. `Ncores = 1` with a worker profile
(`R_PROFILE_USER`) that calls `set.seed`, or `assignInNamespace` a `makeCluster` that runs
`clusterSetRNGStream`. Exact reproduction of step 12 also needs R's `rmultinom` (sequential
`rbinom`, BTPE) for step 10, because those draws advance the master stream before `emptyDrops`
takes its seeds; up to 10 × 100 × (genes − 1) binomial draws (R's `rmultinom` stops once the
size is used up). UNVERIFIED whether the `upper` iteration can cycle.

## 2. DIEM

### Identity

- Paper: Alvarez M, Rahmani E, Jew B, et al. Enhancing droplet-based single-nucleus RNA-seq
  resolution using the semi-supervised machine learning classifier DIEM. *Scientific Reports*
  10:11019 (2020), doi:10.1038/s41598-020-67513-5 (PMID 32620816; resolves to this paper).
  Preprint doi:10.1101/786285.
- Code: <https://github.com/marcalva/diem>, R package `diem` with C++ (Rcpp, RcppEigen,
  RcppProgress, optional OpenMP); imports irlba, Matrix. Not on CRAN, Bioconductor or conda:
  `devtools::install_github("marcalva/diem")`.
- Licence: `License: GPL-3` in DESCRIPTION; no LICENSE file. [SRC]
- Read at tag `v2.4.1` (commit a9a4d9d, 2022-12-19; the repository's last commit):
  `R/filter.R`, `R/var.R`, `R/normalize.R`, `R/init_clust.R`, `R/em.R`, `R/debris_id.R`, `R/post.R`,
  `src/MultEM.cpp`, `src/utils.cpp`, `src/var.cpp`. [SRC]

### Releases

| diem | Tag / commit | Change [SRC unless noted] |
|---|---|---|
| 1.0 | `v1.0` (Nov 2019) | bioRxiv v1: multinomial mixture; `set_debris_test_set(top_n = 1e4, min_counts = 100, min_genes = 100)`, EM `eps = 1e-8`, call on posterior `CleanProb ≥ 0.95` and ≥ 200 genes |
| 2.0.0 – 2.2.0 | `v2.0.0`, `v2.0.1`, `v2.1.0`, `v2.2.0` (Feb 2020) | Dirichlet-multinomial mixtures (`LlkDirMultSparsePar`); `init(k_init = 30, nstart_init = 30)`; call on posterior `pp_thresh = 0.95` |
| 2.3.0 | untagged, ca31a55 – f1c5a39 (Apr – Jun 2020) | debris score and `call_targets(thresh_score = 0.5)`; multinomial default; k-means `k_init = 20, nstart = 30`, clusters < 10 droplets merged; EM `max_iter = 300`. The published paper's method (UNVERIFIED which commit) |
| 2.4.0 | `v2.4.0` (Feb 2022) | EM rewritten in C++ (`mult_em`, `dirmult_em`); `nstart_init` and `min_size_init` no longer used (k-means `nstart = 1`, no merging); `max_iter = 1000` |
| 2.4.1 | `v2.4.1` (Dec 2022) | fixes an out-of-bounds read in `vec_cmplmnt` (the free/fixed droplet split) and `divide_by_colsum` dropping dimnames (with recent Matrix, 2.4.0's debris score can fail) |

A pipeline today installs `master` = 2.4.1. The paper ([DOC]) describes: debris set < 100
counts, multinomial mixture, k-means on 30 PCs of 2000 variable genes over droplets with ≥ 200
genes, k = 20, EM to ε = 1e-4 with pseudocount 1e-10, debris score from Welch-t debris genes
scaled so the lowest cluster mean is 0 and debris clusters' mean 1, keep score < 0.5.

### Algorithm (2.4.1, the README workflow) [SRC]

The README runs: `create_SCE` → `set_debris_test_set` → `filter_genes` → `get_pcs` → `init` →
`run_em` → `assign_clusters` → `estimate_dbr_score` → `call_targets(thresh_score = 0.5)`
(single-nucleus), or `call_targets(clusters = "debris", thresh_score = NULL)` (single-cell).

1. `create_SCE`: drop barcodes with total 0; `make.unique` on gene and barcode names. Counts are
   UMIs over all rows given.
2. `set_debris_test_set(min_counts = 100, min_genes = 0, top_n = NULL)`: test set = barcodes
   with total ≥ `min_counts` and genes ≥ `min_genes`, after raising `min_counts` to the total of
   the `top_n`-th barcode when `top_n` is given; background ("debris set") = every other barcode
   with a count.
3. `filter_genes(cpm_thresh = 0)`: genes with CPM > 0 over all barcodes (the per-group check in
   the function is computed and then overwritten).
4. `get_pcs(min_genes = 200, n_var_genes = 2000, lss = 0.3, n_pcs = 30, seedn = 1)` on the test
   droplets with ≥ 200 genes ("cluster set"): per gene, $$\log_{10}(\text{mean} + 1)$$ and
   $$\log_{10}(\text{var} + 1)$$ of raw counts (variance with divisor $$n - 1$$); `loess(log_var ~
   log_mean, span = 0.3)` (R defaults otherwise); the 2000 genes with the largest residual.
   Normalise: counts of all expressed genes divided by the droplet total, times the median total,
   `log1p`. PCA of the variable genes, centred and scaled: `irlba::prcomp_irlba(n = 30)` after
   `set.seed(1, "Mersenne-Twister")` (or `prcomp` when 30 ≥ half the smaller dimension).
5. `init(k_init = 20, seedn = 1, model = "mltn")`: `scale()` the PCs (each PC to unit variance),
   `set.seed(1)`, `kmeans(pcs, centers = 20)` with R's defaults (Hartigan-Wong, `iter.max = 10`,
   **`nstart = 1`**; the documented `nstart_init = 30` is not passed on), warnings suppressed.
   Labels: background = 1, k-means clusters 2..21; test droplets with < 200 genes start
   unassigned (a zero row of $$Z$$). Initial parameters from the hard labels:
   $$\alpha_{gk} = (c + \sum_i Z_{ik} x_{gi}) / \sum_g(\cdot)$$ and $$\pi_k \propto \sum_i Z_{ik}$$.
   Because `params_frm_clust(x, ac, genes.use, psc, ...)` passes `psc` positionally into
   `alpha_prior`, $$c = 2\cdot10^{-10}$$ here (prior 1e-10 + pseudocount 1e-10).
6. `run_em(model = "mltn", max_iter = 1000, eps = 1e-4, psc = 1e-10, alpha_prior = 0,
   pi_prior = 0, threads = 1)`: multinomial mixture over all droplets with a count, expressed
   genes. Background droplets are fixed in component 1 throughout; test droplets are free. First
   an E-step with the initial parameters, then repeat {M-step: $$\pi$$ from all rows of $$Z$$
   (fixed ones included), $$\alpha_{gk} \propto 10^{-10} + \sum_i Z_{ik} x_{gi}$$; E-step:
   $$Z_{ik} = \text{softmax}_k(\sum_g x_{gi}\log\alpha_{gk} + \log\pi_k)$$} until the mean over
   free droplets of $$\sum_k |Z^{new}_{ik} - Z^{old}_{ik}|$$ is < 1e-4 or 1000 iterations.
   `model = "DM"` runs a Dirichlet-multinomial EM instead (not the README's default).
7. `assign_clusters`: each test droplet to its most probable component (first on ties).
8. `estimate_dbr_score(dbr_clust = 1, thresh_genes = 200, thresh_p = 0.05, thresh_logfc = 0,
   max_genes = 5000, p_method = "fdr", sf = 1)`: debris clusters = 1 plus the clusters whose
   test droplets average < 200 genes. Expression $$\log(1 + x_{gi}/n_i)$$ (scale factor 1).
   Welch t-test per gene, test droplets in debris clusters against the rest (genes with mean > 0
   in both groups); BH; debris genes = adjusted p < 0.05 and $$\log_2$$ fold change > 0, the first
   5000 by mean difference. Score = sum of the droplet's expression over debris genes, minus the
   smallest cluster mean, divided by the mean over droplets in debris clusters.
   **Index bug:** the low-gene clusters are found with `which(tapply(n_genes, Cluster, mean) <
   200)`, which returns positions in the alphabetical order of the labels ("1", "10", "11", …,
   "19", "2", "20", …), and those positions are then matched against the labels. With ≥ 10
   clusters the wrong clusters are flagged (checked in R: a low-gene cluster "10" at position 2
   flags cluster "2"). Present in 2.3.0 too.
9. `call_targets(thresh_score = 0.5, clusters = NULL, min_genes = 200)`: a test droplet is
   `Clean` if score < 0.5 and genes ≥ 200 (and not in `clusters`); otherwise `Debris`. With
   `clusters = "debris", thresh_score = NULL`: `Clean` if not in a debris cluster (same flagged
   set) and genes ≥ 200. Background droplets get no call (not cells).

Randomness and determinism: seeded (`set.seed(1)` before `prcomp_irlba` and before `kmeans`), so
deterministic for a given R, irlba and BLAS. With `threads > 1`, `mult_alpha_ml` updates the
shared $$\alpha$$ from an OpenMP loop over droplets without atomics (a data race), so results are
then not reproducible; the default is 1 thread. [SRC]

### Parameters

`set_debris_test_set(min_counts = 100, min_genes = 0, top_n = NULL, debris_ids = NULL)`,
`filter_genes(cpm_thresh = 0)`, `get_pcs(min_genes = 200, n_var_genes = 2000, lss = 0.3,
n_pcs = 30, seedn = 1)`, `init(k_init = 20, model = "mltn", seedn = 1, n_sample = NULL)`,
`run_em(model = "mltn", max_iter = 1000, eps = 1e-4, psc = 1e-10)`, `estimate_dbr_score(...)`
as above, `call_targets(thresh_score = 0.5, clusters = NULL, min_genes = 200)`. [SRC]

### Output

Per test droplet in `sce@test_data`: `total_counts`, `n_genes`, `Cluster` (1 = debris
component), `ClusterProb`, `score.debris`, `Call` (`Clean` / `Debris`); `get_clean_ids()`
returns the cells. `summarize_clusters()` gives per-cluster type, size, mean genes, mean score
and top genes. [SRC]

### Oracle

R with `diem` built from GitHub at `v2.4.1` (needs a C++ compiler with Eigen; conda-forge
`r-rcppeigen`, `r-rcppprogress`, `r-irlba`), `threads = 1`. Deterministic given seed 1, but
exact reproduction needs R's `kmeans` (Hartigan-Wong, Fortran) from the same start
(`sample.int` under Mersenne-Twister), `irlba`'s random start vector and its iterations, and
R's `loess`; the PCs are sign- and BLAS-sensitive, and k-means with `nstart = 1` is sensitive
to them. A practical oracle compares the final calls with a tolerance, or reproduces exactly
by taking the original's k-means labels as an input fixture and checking the EM, score and call
exactly.

## 3. DropletQC

### Identity

- Paper: Muskovic W, Powell JE. DropletQC: improved identification of empty droplets and damaged
  cells in single-cell RNA-seq data. *Genome Biology* 22:329 (2021),
  doi:10.1186/s13059-021-02547-0 (PMID 34857027; resolves to this paper). Preprint
  doi:10.1101/2021.08.02.454717.
- Code: <https://github.com/powellgenomicslab/DropletQC>, R package `DropletQC` (pure R; imports
  Rsamtools, GenomicAlignments, GenomicFeatures, rtracklayer, ks, mclust, future/furrr). Not on
  CRAN, Bioconductor or conda: `devtools::install_github("powellgenomicslab/DropletQC")`.
- Licence: `License: MIT + file LICENSE`; `LICENSE` (YEAR 2021, Walter Muskovic) and
  `LICENSE.md` (MIT). [SRC]
- Read at tag `v1.0.0` (commit 5d7dadc, 2021-11-18; the last commit):
  `R/nuclear_fraction_tags.R`, `R/nuclear_fraction_annotation.R`, `R/identify_empty_drops.R`,
  `R/identify_damaged_cells.R`. [SRC]

### Releases

One tag, `v1.0.0` (Nov 2021), whose DESCRIPTION still says `Version: 0.0.0.9000`. The functions
reached their final form in Feb 2021 (commit 57f6f44 fixed `nuclear_fraction_tags` returning
scores in alphabetical rather than input order; 2a9ebb6 added the plot option); the package was
renamed from `dropletQC` in May 2021. One variant. [SRC]

### Algorithm [SRC unless noted]

DropletQC flags barcodes it is given; the paper ([DOC]) applies it **after** another caller
(EmptyDrops with `lower` raised to 500 and FDR 1 %, or Cell Ranger, CellBender, EmptyNN), and
`nuclear_fraction_tags(outs = ...)` reads `filtered_feature_bc_matrix/barcodes.tsv.gz` by
default. On a whole raw matrix the density would be dominated by empty droplets; whether the
rule still works there is untested (UNVERIFIED).

**Nuclear fraction**, `nuclear_fraction_tags(bam, barcodes, tiles = 100, cell_barcode_tag =
"CB", region_type_tag = "RE", exon_tag = "E", intron_tag = "N")`: the genome (sequences with
mapped reads, from `idxstatsBam`) is cut into `tiles` ranges; for each, `scanBam` reads the
`CB` and `RE` tags of every alignment overlapping the range (default `ScanBamParam` flags: so
PCR duplicates, secondary alignments and multimappers are all counted; an alignment that crosses
a tile boundary is counted in both tiles). Records without `CB` or `RE`, or whose `CB` is not
in `barcodes`, are dropped. Per barcode: $$\text{NF} = N/(N + E)$$ over alignment records with
`RE = N` (intronic) and `RE = E` (exonic); `RE = I` (intergenic) is ignored; 0/0 → 0. These are
**reads, not UMIs**. `nuclear_fraction_annotation(annotation_path, bam, barcodes, tiles =
1000)` is the alternative for BAMs without `RE`: exons (merged, strand ignored) and introns
minus exons from the GTF/GFF, `findOverlaps(minoverlap = 5)`, a read overlapping both counted
in both.

**Empty droplets**, `identify_empty_drops(nf_umi, nf_rescue = 0.05, umi_rescue = 1000)`, input
a data frame (NF, total UMIs):

1. `ks::kdde(nf, deriv.order = 1)`: binned kernel estimate of the density's first derivative,
   Gaussian kernel, bandwidth `hpi(nf, nstage = 2, binned = TRUE, deriv.order = 1)` (the
   plug-in selector for the derivative), 401 grid points from $$\min - (3.7 + 1)h$$ to
   $$\max + (3.7 + 1)h$$ (ks 1.15.3 `binning`). (A density estimate with its own bandwidth is
   also computed, for the plot only.)
2. Run-length encode `derivative > 0` along the grid; cutoff = the grid point at index
   (length of run 1 + length of run 2), i.e. the last point of the second run. When the
   derivative starts positive this is the first local minimum after the first peak; if it
   starts negative or zero it is the first peak instead (no check). Fewer than 4 runs gives a
   warning only.
3. `empty_droplet` if NF < cutoff, else `cell`; then `cell` again if NF > `nf_rescue` and UMI >
   `umi_rescue` (strict inequalities).

**Damaged cells** (optional), `identify_damaged_cells(nf_umi_ed_ct, nf_sep = 0.15,
umi_sep_perc = 50)`: input (NF, UMI, the empty-droplet call, a **cell type** per barcode). For
the barcodes called `cell`, per cell type: `mclust::Mclust(cbind(NF, log10 UMI), G = 1:2,
modelNames = "EEI")` (BIC picks G). Damaged only if G = 2, the component with the higher NF mean
has the lower UMI mean, the NF means differ by more than 0.15, and $$10^{\mu_{\text{umi,dam}}}
< 10^{\mu_{\text{umi,cell}}}(1 - 0.5)$$; then members of the high-NF component (mclust's MAP
classification) are `damaged_cell`.

Determinism: no random numbers anywhere (ks's plug-in bandwidth and binned estimate, mclust's
hierarchical initialisation are deterministic). [SRC]

### Parameters

`nuclear_fraction_tags(outs | bam + barcodes, tiles = 100, cores, cell_barcode_tag = "CB",
region_type_tag = "RE", exon_tag = "E", intron_tag = "N")`; `identify_empty_drops(nf_umi,
nf_rescue = 0.05, umi_rescue = 1000)`; `identify_damaged_cells(nf_umi_ed_ct, nf_sep = 0.15,
umi_sep_perc = 50)`. [SRC]

### Output

`nuclear_fraction_*`: a data frame of NF per barcode. `identify_empty_drops`: the input with
`cell_status` ∈ {`cell`, `empty_droplet`}. `identify_damaged_cells`: `df` with `cell_status` ∈
{`cell`, `empty_droplet`, `damaged_cell`} (+ plots). Flags only; nothing is removed. [SRC]

### Oracle

R with Bioconductor (Rsamtools, GenomicAlignments, GenomicFeatures, rtracklayer), CRAN ks and
mclust, and DropletQC from GitHub at `v1.0.0`. The rule itself needs only (NF, UMI) per barcode,
so the oracle can feed `identify_empty_drops` an NF table we compute, and check it exactly; that
splits the test in two: (a) the NF from a Cell Ranger BAM (needs a BAM: none of the repo's
public matrices has one; a small public 10x BAM, or the maintainer's run), and (b) the KDE rule
on (NF, UMI). Exactness in (b) needs ks's `hpi` (two-stage plug-in, binned) and binned `kdde`
reproduced to the grid point; pin the ks version (1.15.3 current; check older ones for changes
to `hpi`). An NF from velocyto/STARsolo/kb **UMIs** (the repo's `unspliced_fraction`) is a
different quantity from DropletQC's read-based NF: usable as an approximation, not as the
original.

## 4. DropSift

### Identity

- No paper. Broad Institute (McCarroll lab; authors James Nemesh, Alec Wysoker); used on Drop-seq
  and 10x single-nucleus data, inputs from Drop-seq tools or the Optimus pipeline.
- Code: <https://github.com/broadinstitute/DropSift>, R package `DropSift` (pure R; imports
  e1071, Seurat, hdrcde, pastecs, cluster, Matrix, data.table, rhdf5). Not on CRAN, Bioconductor
  (git.bioconductor.org has no such package; NEWS's "1.0.0 (Bioconductor Release)" heading is
  aspirational) or conda; `remotes::install_github("broadinstitute/DropSift")`, and a Docker
  image built by the repository's workflow.
- Licence: `License: MIT + file LICENSE`; `LICENSE` (YEAR 2024, Broad Institute) and
  `LICENSE.md` (MIT). [SRC]
- Read at commit ccaefb7 (2026-10-29, `main`, DESCRIPTION `Version: 0.99.2`; there are **no
  tags**): `R/SvmNucleusCaller.R`, `R/svmIntronicCellSelection.R`, `R/findTrainingDataBounds.R`,
  `R/CellCallingUtils.R`, `R/geneModuleScore.R`, `R/parseOptimus.R`. [SRC]

### Releases

No tags or releases; DESCRIPTION went 1.0 (Dec 2024) → 0.99.0 (Feb 2025) → 0.99.1 (May 2026) →
0.99.2 (Jun 2026) and has not changed since, while the calling changed several times. The
variants have to be named by commit [SRC]:

| Commits | Behaviour |
|---|---|
| f502f47 (2024-12-17) – be3decb (2026-04-23) | one SVM, nucleus vs empty; call = the SVM's predicted class (or probability ≥ `cellProbabilityThreshold` if given). 4b1bef3 (Apr 2025) added `useCBRBInitialization`; 28d30d0 (Jun 2025) `random.seed`; 66f4aaf (Nov 2025) a debris fix |
| c0794aa (#35, 2026-05-27), ed93bff (#36, 06-01) | initialisation revised (debris-filtered threshold learning, candidate stability score); experimental 2-D HDR refinement, off by default |
| 6cbee01 (#38, 2026-06-06) – 38bfe17 | a third class, **debris** (high UMI, low intronic), with its own gene module and two more SVMs; class from probabilities ≥ 0.5 |
| f059e3d (#43, 2026-09-18) – ccaefb7 | automatic fallback to the two-peak solution (`twoClusterFallbackRatio = 2`); current |

A pipeline installing from GitHub today gets ccaefb7. The notes below describe it.

### Inputs [SRC]

A cell-features table, one row per barcode: `cell_barcode`, `num_transcripts` (UMIs; computed
from the matrix if absent), `pct_intronic`, `pct_mt` (fractions, despite the names), and with
CellBender `frac_contamination` (+ `num_retained_transcripts`, else approximated as
$$\text{round}((1 - f)\,n)$$). From Optimus h5ad: `pct_intronic = (reads_mapped_intronic +
reads_mapped_intronic_as) / reads_mapped_uniquely`, `pct_mt = reads_mapped_mitochondrial /
reads_mapped_uniquely`, barcodes with ≥ 20 UMIs. From Drop-seq tools: the per-barcode read
metrics. So `pct_intronic` is a **read** fraction from the aligner, not from UMIs. Plus the UMI
matrix (genes × barcodes) for the gene-module scores; zero-count barcodes are dropped.

### Algorithm (ccaefb7) [SRC]

`SvmNucleusCaller(cellFeatures, dgeMatrix, cellProbabilityThreshold = NULL, maxUmisEmpty = 50,
useCBRBFeatures = TRUE, useCBRBInitialization = useCBRBFeatures, forceTwoClusterSolution =
FALSE, use2DTrainingRefinement = FALSE, twoClusterFallbackRatio = 2)`; `runIntronicSVM(...,
random.seed = NA)` wraps it with file input/output.

**Features**: `num_transcripts`, `pct_intronic`, `pct_mt`, + `frac_contamination` when
`useCBRBFeatures`, + `empty_gene_module_score` when a matrix is given. Before each SVM:
`num_transcripts` → $$\log_{10}$$, then every feature z-scored over **all** barcodes.

**Exemplars, CellBender path** (`useCBRBInitialization = TRUE`, the default when CellBender
features are present): empty = barcodes with `frac_contamination == 1` and > 20 UMIs, bounded by
the 75 % highest-density region (HDR) of $$\log_{10}$$ UMI and of `pct_intronic`. Nuclei =
barcodes with `frac_contamination < 1` at or below its 25th percentile, 95 % HDR bounds; then
extended: barcodes with $$\log_{10}$$ UMI > empty lower bound + 1, `pct_intronic` inside the
nucleus bounds and `frac_contamination ≤ 0.1`, 95 % HDR, and the union of the two boxes. Labels:
empty if inside the empty box and `frac_contamination == 1`; nucleus if inside the nucleus box
and `frac_contamination ≤ 0.1`. No debris class on this path.

**Exemplars, density path** (no CellBender, or `useCBRBInitialization = FALSE`):

1. Two runs, with UMI floors 50 (`maxUmisEmpty`) and 20; the default rule is the pit after the
   highest peak. `x` = $$\log_{10}(\text{UMI} + 1)$$ of barcodes ≥ the floor.
2. Threshold rule: `stats::density(x, adjust = a, n = 200)` (Gaussian, `bw.nrd0` × a) for
   $$a = b^{m}$$, $$b \in \{0.25, 0.5, \dots, 3\}$$ (each value raised to the power of the
   smoothing multiplier $$m$$ of step 3); turning points with `pastecs::turnpoints`; the first
   pit after the highest peak. Across the 12 values of $$a$$: for each, the median and MAD of
   the thresholds in a window of ±2 neighbours; the threshold = the local median with the
   smallest relative MAD. (Two-peak variant: the lowest pit between the two highest peaks,
   median over the grid.)
3. For smoothing multipliers $$m = 1.1^{0..14}$$ (15 iterations): threshold from all `x` →
   provisional empty bounds → drop barcodes above the empty UMI bound with `pct_intronic < 0.25`
   (likely debris) → threshold again. Empty box: among barcodes below the threshold, the 75 % HDR
   of $$\log_{10}$$ UMI (`hdrcde::hdr`, bandwidth `hdrbw` ×1, ×2, ×8, … until ≤ 2 intervals),
   the interval containing the density maximum, then the 95 % HDR of `pct_intronic` inside it;
   UMI upper bound capped at the threshold. Nucleus box: barcodes above the threshold with
   `pct_intronic` ≥ a floor (two candidates: 0.9 × the empty box's lower intronic bound capped at
   0.25, or 0.25), 75 % HDRs of `pct_intronic` and then of UMI; UMI upper bound = the 99.5th
   percentile; lower bound ≥ threshold.
4. Each candidate (multiplier × floor) is scored by the mean silhouette (`cluster::silhouette`,
   Euclidean, on min-max-rescaled $$\log_{10}$$ UMI and `pct_intronic`) of a **random 10 %**
   (`sample`) of its exemplars; selection score = median silhouette over ±2 neighbouring
   candidates − 0.02 × (1 − stability), stability = the smaller of the mean overlaps
   (min(containment, Jaccard)) of the empty and nucleus sets with neighbours. Best score wins.
5. Between the 50- and 20-UMI runs: the higher selection score, among those with
   empty/nucleus exemplar ratio > 0.2 (if any). If that score's silhouette ≤ 0.5, the two-peak
   solution is also computed and adopted when its silhouette is ≥ 2× (or when the default's is
   ≤ 0 and the two-peak's > 0).
6. Debris box: barcodes above the threshold with `pct_intronic` < the floor; those at or above
   their median UMI; 95 % HDRs of `pct_intronic` then UMI, clipped to UMI ≥ threshold and
   intronic ≤ floor; UMI upper bound the 99.5th percentile; ≥ 10 barcodes needed.

Labels are given by box membership, testing $$\log_{10}(\text{UMI})$$ (without the +1 used to
build the boxes); debris overrides empty overrides nucleus where boxes overlap.

**Gene-module scores** (for `negative_class` = empty, then debris; each needs ≥ 100 exemplars
per class): the negative exemplars' counts are pooled in **random** order (`sample`) into
pseudo-barcodes whose sizes match the nuclei's (shuffled) UMI totals; Seurat `LogNormalize`
(10,000); `FindMarkers(nuclei vs pseudobulk, test.use = "wilcox", min.pct = 0.25,
logfc.threshold = 1)`, adjusted p ≤ 0.05; the 100 genes most **down** in nuclei; the score is
Seurat `AddModuleScore` on all barcodes (which calls `set.seed(1)` itself). A failed module is
dropped with a warning.

**SVMs**: `e1071::svm(type = "C-classification", kernel = "radial", probability = TRUE)`, other
arguments default (cost 1, gamma = 1/#features), trained on exemplars only, predicting every
barcode: nucleus vs empty (features above); nucleus vs debris (empty module replaced by the
debris module); debris vs empty (both modules). `is_cell_prob` = min(P(nucleus|vs empty),
P(nucleus|vs debris)); `is_debris_prob` = min(P(debris|vs empty), 1 − P(nucleus|vs debris)).

**Call**: `nucleus` if `is_cell_prob` ≥ 0.5 (or `cellProbabilityThreshold`), `debris` if
`is_debris_prob` ≥ 0.5, both → the larger (ties nucleus), else `empty_or_other`; barcodes with
$$\log_{10}(\text{UMI} + 1)$$ below the empty box's upper UMI bound get NA.

Randomness: R's RNG in the silhouette subsample, the pseudobulk shuffles, and libsvm's
probability cross-validation (e1071 uses `unif_rand`); Seurat's `AddModuleScore` re-seeds with 1
in the middle. `random.seed = NA` by default: **not reproducible run to run** unless the caller
seeds. Versions of Seurat (FindMarkers' fold change and Wilcoxon backend changed in Seurat 5;
presto if installed), e1071, hdrcde and R's `density` all enter the result.

### Output

`cell_features` with `training_label_class`, `empty_gene_module_score`,
`debris_gene_module_score`, `p_nucleus_vs_empty`, `p_nucleus_vs_debris`, `p_debris_vs_empty`,
`is_cell_prob`, `is_debris_prob`, `barcode_class` ∈ {`nucleus`, `debris`, `empty_or_other`,
NA}; plus the bounds, exemplar barcodes, DE gene tables, the SVM models, and
`getCBRBArgs()` (CellBender `expected_cells` and `total_droplets_included` estimated from the
calls). [SRC]

### Oracle

R ≥ 4.3 with DropSift at a pinned commit, Seurat 5 and SeuratObject (CRAN), e1071, hdrcde 3.5,
pastecs; `random.seed = 1`. Inputs need per-barcode intronic and mitochondrial **read** fractions
(an Optimus h5ad, Drop-seq metrics, or our own from a BAM or STARsolo's per-barcode
`CellReads.stats`, UNVERIFIED that its columns suffice), and optionally CellBender's
`frac_contamination`. Exact reproduction would require R's RNG stream through `sample`, Seurat's
Wilcoxon and module score, `density`, hdrcde's root-finding, and libsvm's Platt scaling:
practically, an oracle compares calls with a tolerance, or checks the stages separately (exemplar
boxes for a fixed seed; SVM probabilities given fixed exemplars and features).

## 5. valiDrops

### Identity

- Paper: Kavaliauskaite G, Madsen JGS. Automatic quality control of single-cell and
  single-nucleus RNA-seq using valiDrops. *NAR Genomics and Bioinformatics* 5(4):lqad101 (2023),
  doi:10.1093/nargab/lqad101 (PMID 38025048; resolves to this paper).
- Code: <https://github.com/madsen-lab/valiDrops>, R package `valiDrops` (pure R; imports
  segmented, zoo, mixtools, inflection, robustbase, irlba, Seurat, scry, presto, glmnet). Not on
  CRAN, Bioconductor or conda: `devtools::install_github("madsen-lab/valiDrops")`.
- Licence: `License: GPL-3` in DESCRIPTION; a GPL-3 `LICENSE` file was added on 2025-03-18
  (f67d751, the only commit after the tag). [SRC]
- Read at tag `v1.0.0` (commit d2da5be, 2023-09-25; DESCRIPTION still `Version: 0.1.0`):
  `R/valiDrops.R`, `R/rank_barcodes.R`, `R/quality_metrics.R`, `R/quality_filter.R`,
  `R/expression_metrics.R`, `R/expression_filter.R`, `R/label_dead.R`. [SRC]

### Releases

One tag, `v1.0.0` (Sep 2023); `main` differs only by the LICENSE file. Before the tag: the
three-stage pipeline settled in Nov 2022; 0834284 (Sep 2023) made stage 3 optional
(`stageThree`); 905cd3a and the next two commits (Sep 2023) accept `DelayedMatrix`. One variant
(v1.0.0). Behaviour does depend on the installed **segmented**: up to 2.0-2 its bootstrap
restarts draw from the global RNG; from 2.0-3 (Feb 2024) with `seed = NULL` it seeds itself from
the digits of the mean of y and restores the caller's RNG state, so stage 1 became
deterministic. [SRC, segmented `R/seg.lm.fit.boot.r` across its CRAN versions]

### Algorithm (`valiDrops(counts, rank_barcodes = TRUE, stageThree = TRUE, ...)`) [SRC]

Input: the raw UMI matrix (genes × barcodes). Gene identity matters: species and annotation are
auto-detected against bundled tables (human, mouse, rat, zebrafish, worm, fly; symbol, Entrez,
Ensembl, HGNC, MGI; version suffixes stripped) to find mitochondrial (chromosome MT), ribosomal
(`^RPL`/`^RPS`) and protein-coding genes.

**Stage 1, `rank_barcodes(type = "UMI", psi.min = 2, psi.max = 5, alpha = 0.001, alpha.max =
0.05, boot = 10, factor = 1.5)`**: totals > 0, `rank(-total)` (average ranks), one point per
distinct total (its average rank); $$x = \ln \text{rank}$$, $$y = \ln \text{total}$$; centred
moving averages of width $$\lceil 2N^{1/3}\rceil$$ ($$N$$ distinct totals, `zoo::rollmean`). For
$$k = 2..5$$ breakpoints: `segmented(lm(y ~ x), psi = k values evenly spaced between the
0.001 and 0.999 quantiles of x, seg.control(alpha = 0.001 * 0.999, n.boot = 10))`, on a "psi
values too close" error retrying with alpha + 0.001 up to 0.05. Model: the first (fewest
breakpoints) whose RMSE is at most 1.5 × the smallest. Angles between consecutive segments
$$\theta_i = \arctan\frac{s_i - s_{i+1}}{1 + s_i s_{i+1}}$$ (degrees); breakpoint = argmin
over $$i \ge 2$$ (the first breakpoint is never chosen); threshold = the distinct total whose
$$\ln$$ rank is nearest that breakpoint, returned as `exp(log(total))`. Barcodes with total ≥
that value go on; since `exp(log(x))` need not equal $$x$$ in floating point, the barcodes at
exactly the threshold total can fall on either side (a port must compute the same `exp(log)`).

**Stage 2, `quality_metrics` + `quality_filter(mito = TRUE, distance = TRUE, coding = TRUE,
mito.nreps = 10, mito.max = 0.3, npsi = 3, dist.threshold = 5, coding.threshold = 3)`**, per
barcode: $$\ln$$ UMIs, $$\ln$$ genes, mitochondrial, ribosomal and protein-coding UMI fractions.
Applied in sequence, each on the survivors of the last:

1. Mitochondrial: 10 times, `mixtools::normalmixEM(ln genes)` (2 components, **random**
   initialisation: `normalmix.init` draws `rexp`/`rnorm`), take the barcodes in the higher-mean
   component; over $$t$$ from their median mito fraction to 1 by 0.001, $$c(t)$$ = the number
   with mito ≤ $$t$$; threshold `inflection::uik(t, c)`; median of the 10. If > 0.3: 10 times on
   a random subsample of min(5000, 0.8n) barcodes, `segmented(lm(ln genes ~ mito), npsi = 1..5)`
   until the smallest breakpoint is ≤ 0.3; median. Keep mito ≤ threshold.
2. Distance: `segmented(lm(ln genes ~ ln UMIs), npsi = 3)` (fewer if it fails); keep residuals
   within median ± 5 × `robustbase::Sn`.
3. Coding: keep the coding fraction within median ± 3 × Sn.

**Stage 3, `expression_metrics(nfeats = 5000, npcs = 10, k.min = 5, res.shallow = 0.1, top.n =
10)` + `expression_filter(mito = 3, ribo = 3, min.significant = 1, min.target.pct = 0.3,
max.background.pct = 0.7, min.diff.pct = 0.2, min.de.frac = 0.01)`**, on protein-coding genes of
the stage-2 barcodes: normalise to 10,000 and `log1p`; 5000 genes by `scry`'s binomial deviance;
z-score; `irlba` (10 PCs, random start); Seurat `FindNeighbors` (SNN) and `FindClusters` (Louvain)
at resolution 0.1 ("shallow") and at a "deep" resolution chosen so the smallest cluster has
about 5 barcodes (coarse grid 1..20, then ±0.9 by 0.1). For each deep cluster against the
barcodes outside its majority shallow cluster: genes detected in ≥ 10 % of either group with
$$\log_2$$ fold change ≥ 0.25; Wilcoxon (presto if installed, else `wilcox.test`), Bonferroni over
all expressed genes. From the top markers (at most 10 with FDR ≤ 0.05, by p; the first gene alone
when none is significant): mean $$(p_1 - p_2)/p_1$$, mean $$p_1$$, mean $$p_2$$, the number with
$$(p_1 - p_2)/p_1 < -0.01$$, the number significant, the smallest FDR, the DE fraction, and the
cluster's median mito and ribo fractions. A cluster is valid if: no negative markers, mean
$$(p_1-p_2)/p_1 \ge 0.2$$, mean $$p_1 \ge 0.3$$, mean $$p_2 \le 0.7$$, at least one significant,
$$-\log_{10}$$ of the smallest FDR at or above a data-driven level (the smaller of median + 3 Sn
over clusters with pct.diff ≤ 0.4, and a one-breakpoint `segmented` fit of pct.diff on
$$-\log_{10}$$ FDR), DE fraction > 0.01, and mito and ribo medians at most median + 3 Sn across
clusters. Cells = the barcodes of valid clusters.

`label_dead = TRUE` (off by default) adds a dead/live/uncertain label from a ridge `cv.glmnet`
with resampling; not part of the call. [SRC]

Randomness: `normalmixEM` (and the mito fallback's subsample), `irlba`'s start vector, and with
old segmented the bootstrap; Seurat's `FindClusters` uses its own `random.seed = 0`. The package
sets no seed: **not reproducible without `set.seed`** by the caller. [SRC]

### Parameters

As listed above; `valiDrops(counts, rank_barcodes = TRUE, mitochondrial_clusters = 3,
ribosomal_clusters = 3, label_dead = FALSE, stageThree = TRUE, ...)` passes the rest to each
stage (`species`, `annotation`, and the `mito`, `ribo`, `coding` gene lists can be given). [SRC]

### Output

A data frame per barcode that passed stage 1: `barcode`, `logUMIs`, `logFeatures`,
`mitochondrial_fraction`, `ribosomal_fraction`, `coding_fraction`, `qc.pass` (`pass`/`fail`),
and `label` (`live`/`dead`/`uncertain`) with `label_dead`. Barcodes below the stage-1 threshold
are not listed (not cells). [SRC]

### Oracle

R ≥ 4.0.2 with valiDrops from GitHub at `v1.0.0`, Seurat, scry (Bioconductor), presto (GitHub
`immunogenomics/presto`), segmented ≥ 2.0-3, inflection, mixtools, robustbase; `set.seed` before
the call. Stage 1 alone (`rank_barcodes`) is deterministic with segmented ≥ 2.0-3 and a good exact
oracle target; segmented's fitting (bootstrap restarts from a seed derived from the decimal
digits of the mean of y) would have to be ported exactly. Stages 2–3 involve EM with random
starts, irlba, Louvain and Wilcoxon with presto: an oracle can only compare with a tolerance, or
check stages given the original's intermediate outputs (e.g. its clusters).

## 6. SiftCell

### Identity

- Paper: Xi J, Park SR, Lee JH, Kang HM. SiftCell: a robust framework to detect and isolate
  cell-containing droplets from single-cell RNA sequence reads. *Cell Systems* 14(7):620–628.e3
  (2023), doi:10.1016/j.cels.2023.06.002 (PMID 37473732; resolves to this paper).
- Code: <https://github.com/jyxi7676/SiftCell>, R package `SiftCell` (R + Rcpp/RcppArmadillo) that
  shells out to a Python script (`inst/runXGboost.py`: pandas, scikit-learn, xgboost). Depends on
  `squat` (<https://github.com/hyunminkang/squat>, v1.11, Apache-2.0) for its binomial
  overdispersion test. Not on CRAN, Bioconductor, PyPI or conda.
- Licence: `License: GPL (>= 2)` in DESCRIPTION; no LICENSE file. [SRC]
- Read at `master` commit adaaec1 (2023-07-23), compared with tag `v1.0.0` (298dab1,
  2023-01-04): `R/SiftCellShuffle.R`, `R/SiftCellBoost.R`, `R/SiftCellMix.R`,
  `src/dgeShuffle.cpp`, `src/getAmbientProp.cpp`, `inst/runXGboost.py`; squat at 3da0ca4
  (`R/squat.R`, `src/squat.cpp`). [SRC]

### Releases

One tag, `v1.0.0` (Jan 2023); DESCRIPTION `Version: 1.0`. Between the tag and `master`
(Mar–Jul 2023): `SiftCellShuffle` reads the MTX with `readMM` instead of parsing lines;
`SiftCellBoost` only frees memory; the XGBoost script is unchanged in substance;
`SiftCellMix` was rewritten and `SiftCellBayes` removed. For the call (Shuffle + Boost) one
variant. The repository also commits compiled `.o`/`.so` files. [SRC]

### Algorithm (SiftCell-Shuffle + SiftCell-Boost) [SRC]

Input: a directory `workingdir/DGE/` with uncompressed `matrix.mtx`, `barcodes.tsv`,
`genes.tsv` (10x v2 names; gene names, column 2, made unique), the raw UMI matrix.

**SiftCell-Shuffle** (`SiftCellShuffle(workingdir)`): expand the matrix into one (gene, barcode)
pair per UMI, then shuffle the gene labels with a Fisher–Yates-style pass where the partner of
position $$i$$ is $$i + \text{round}(U(0, T - i - 1))$$ (R's `runif` through Rcpp; the rounding
gives the two end positions half weight, so the permutation is not exactly uniform). Barcode
totals and gene totals are preserved; each barcode's composition becomes a draw from the pooled
gene profile, i.e. a cell-free null. Written to `workingdir/shuffleDGE/`.

**SiftCell-Boost** (`SiftCellBoost(workingdir, threshold = 100, expectedN = 1000, dataName =
"Sample", mitoGenes = TRUE, otherFlagGenes = NULL, N.PCs = 100, geneUMIthreshold = 50, seed = 0,
quantileReg = FALSE, labels = FALSE)`):

1. Barcodes with > 100 UMIs (`threshold`), in the original and the shuffled matrix (the same
   barcodes, since totals are kept).
2. Flagged barcodes: for each `^MT-|^mt-` gene, its fraction per barcode; cutoff = median + 3 SD
   over barcodes; flagged if any mito gene exceeds its cutoff.
3. Positive exemplars ("Nonsoup"): with $$p_0$$ = pooled gene proportions over all barcodes with a
   count and weights $$w_g = \sqrt{p_{0g}/(1-p_{0g})}$$, each barcode's squat score
   `squat_binom_overdisp_test(x, size = total, p = p0, w)`: over its detected genes, a per-gene
   z for how extreme $$x_g$$ is under Binomial(total, $$p_{0g}$$) conditioned on $$x_g > 0$$
   (expectation-based, variance-adjusted, capped at ±10), combined as
   $$\sum w z / \sqrt{\sum w^2}$$. Positives = barcodes ranked ≤ 1000 (`expectedN`) both by UMIs
   and by score (`rank(..., ties.method = "random")`).
4. Variable genes: genes with > 50 UMIs; per gene a squat test across barcodes (sizes = barcode
   totals, weights $$\sqrt{\text{size}}$$); the top 1000 by $$-\log_{10} p$$ (or, with
   `quantileReg`, those above a logistic quantile-regression curve at τ = 0.8).
5. `log(x / total × 10^4 + 1)` on those genes for original and shuffled barcodes together; `prcomp`
   (centred, scaled, 100 PCs) of the stacked $$2N$$ barcodes.
6. Training set: positives labelled `Nonsoup`; flagged barcodes and **all** shuffled barcodes
   labelled `Soup`. Prediction set: every other original barcode. (`labels =` a file adds manual
   positives/negatives.)
7. Python: `xgboost.XGBClassifier(objective = "binary:logistic")` with every other setting at
   the installed version's default, fit on the PCs, `predict` (class, i.e. probability > 0.5).
8. Cells = the positives + predicted `Nonsoup`; written to `CellContaining.txt`. Flagged
   barcodes and barcodes ≤ 100 UMIs are never cells.

Defects that affect a run [SRC]: `seed = 0` (default) sets no seed (only `seed = NULL` seeds from
the clock), so the shuffle and the random tie-breaking are unseeded; the Python script is run
twice and its failure is silenced (`try(..., silent = TRUE)`); string class labels need
**xgboost < 1.6** (1.5.2 encodes them; 1.6.0 and later raise "Invalid classes inferred",
checked in xgboost's `sklearn.py` at v1.5.2 / v1.6.0), so with a current xgboost no output file is
written; `s.strip(dataName + "Org_")` strips characters, not a prefix (harmless for "Sample",
corrupts barcodes when `dataName` contains A, C, G or T); `workingdir` must end in `/`.

**SiftCell-Mix** (not a caller; needs a cell-type label for every barcode above the threshold):
ambient profile = mean over barcodes with a count of their gene proportions weighted by
$$\ln(\min(\text{total}, 100) + 1)$$, normalised; type profiles = $$0.99 \times$$ the type's pooled
proportions + $$0.01 \times$$ ambient; per barcode, `Rsolnp::solnp` maximises
$$\log\sum_k \pi_k \exp(\sum_g x_g \log q_{kg})$$ over the simplex (parallel workers). As coded
this is a mixture over whole-droplet likelihoods, so the optimum tends to a single type
(UNVERIFIED against the paper's description of an ambient fraction per droplet).

### Parameters

As above; the call is driven by `threshold` (100), `expectedN` (1000), `geneUMIthreshold` (50),
`N.PCs` (100), `mitoGenes` (TRUE), and the xgboost defaults of the installed version. [SRC]

### Output

`CellContaining.txt` (one barcode per line) in `workingdir`; `labeledData.csv` and
`predData.csv` (PCs) beside it. No per-barcode probability is kept (the script calls `predict`,
not `predict_proba`). [SRC]

### Oracle

R with SiftCell built from GitHub (Rcpp, RcppArmadillo; squat from GitHub; quantreg, Rsolnp,
reticulate), Python with pandas, scikit-learn and **xgboost 1.5.x** (conda-forge `py-xgboost` has
1.5.0 and 1.5.1; PyPI has 1.5.2); `set.seed()` called by the harness before `SiftCellShuffle` and
`SiftCellBoost`. Deterministic then for fixed versions (xgboost's default `tree_method` and threads
aside). Exact reproduction needs R's `runif` stream for the shuffle (the repo has R's
Mersenne-Twister), squat's test ported exactly, `prcomp`, and xgboost itself (callable from Python,
so the port could call the same xgboost on the same PCs). The 1000-barcode positive set and the
shuffled null are the parts worth validating exactly.

## 7. miQC

### Identity

- Paper: Hippen AA, Falco MM, Weber LM, et al. miQC: an adaptive probabilistic framework for
  quality control of single-cell RNA-sequencing data. *PLoS Computational Biology*
  17(8):e1009290 (2021), doi:10.1371/journal.pcbi.1009290 (PMID 34428202; resolves to this
  paper).
- Code: Bioconductor `miQC` (<https://git.bioconductor.org/packages/miQC>), mirrored at
  <https://github.com/greenelab/miQC>. Pure R; imports flexmix (CRAN, GPL ≥ 2), splines,
  SingleCellExperiment, ggplot2. bioconda `bioconductor-miqc` 1.0.0 – 1.18.0.
- Licence: `License: BSD_3_clause + file LICENSE`; `LICENSE` is BSD 3-Clause, "Copyright (c)
  2021, Greene Laboratory". [SRC]
- Read at Bioconductor branch `RELEASE_3_24` (miQC 1.22.0): `R/mixtureModel.R`,
  `R/filterCells.R`, `R/get1DCutoff.R`; earlier releases by `git diff` between `RELEASE_3_*`
  branches. [SRC]

### What it is

Not an empty-droplet caller: miQC removes **compromised cells** (high mitochondrial fraction for
their number of detected genes) from barcodes already called as cells. Its input is scater's
per-cell QC (`addPerCellQC(sce, subsets = list(mito = ...))`) on the called cells. On a whole
raw matrix the two-component mixture would model empty droplets versus cells instead (untested,
UNVERIFIED). In sc-droplets it is a "miQC-style mixture regression" QC step, applied after a
caller, as STATUS.md plans.

### Releases [SRC]

| miQC | Bioconductor | Change |
|---|---|---|
| 1.0.0 | 3.13 (May 2021) | `mixtureModel(model_type = "linear" / "spline" / "polynomial")`; `filterCells` keeps cells with posterior of the compromised component ≤ 0.75 |
| 1.2.0 | 3.14 (Oct 2021) | `model_type = "one_dimensional"`; `filterCells(keep_all_below_boundary = TRUE, enforce_left_cutoff = TRUE)`, both on by default (changes the default result) |
| 1.4.0 | 3.15 | none |
| 1.6.0 | 3.16 (Nov 2022) | `filterCells` stops with an error when no cell fails the posterior cut |
| 1.8.0 | 3.17 (Apr 2023) | `get1DCutoff` added; filtering unchanged |
| 1.10.0 – 1.22.0 | 3.18 – 3.24 | none |

Variants for the filter: 1.0.0 (posterior only), 1.2.0 – 1.22.0 (with the boundary and left
cutoff; 1.6.0+ errors where 1.2–1.4 returned every cell). A pipeline today runs 1.18 – 1.22.

### Algorithm (1.2.0 – 1.22.0) [SRC]

Inputs per cell: `detected` (genes with count > 0) and `subsets_mito_percent` (mito UMIs / total
× 100), read from `colData` under exactly those names (the `detected` and `subsets_mito_percent`
arguments of `mixtureModel` exist but are not used in the formula).

1. `mixtureModel(sce, model_type = "linear")`: `flexmix(subsets_mito_percent ~ detected, k = 2)`,
   a two-component mixture of Gaussian linear regressions fitted by EM (flexmix defaults:
   **random** initial hard assignment with `sample`, `iter.max = 200`, `minprior = 0.05` (a
   component whose prior falls below is dropped; a warning if only one remains), `tolerance =
   1e-6`, one start). `"spline"`: `~ bs(detected)` (cubic B-spline, df 3); `"polynomial"`:
   `~ poly(detected, 2)`; `"one_dimensional"`: `~ 1` (a two-component Gaussian mixture on mito %).
2. Compromised component = the one with the larger intercept (`parameters(model)[1]`).
3. `filterCells(sce, model, posterior_cutoff = 0.75)`: keep if the posterior of the compromised
   component ≤ 0.75. Error (1.6.0+) if every cell is kept.
4. `keep_all_below_boundary = TRUE`: also keep every cell whose mito % is below the intact
   component's fitted line at its `detected`.
5. `enforce_left_cutoff = TRUE`: let $$m$$ = the lowest mito % among discarded cells and $$d$$ = the
   `detected` of the first cell (column order) with mito % equal to $$m$$; discard every cell with
   `detected` ≤ $$d$$ and mito % ≥ $$m$$.
6. Returns the kept cells. `get1DCutoff` (1.8.0+) returns $$m$$ for a one-dimensional model, or
   errors if kept and discarded mito ranges overlap.

Randomness: flexmix's initialisation uses R's RNG; the package sets no seed and its error message
suggests "setting a different random seed". Deterministic under `set.seed` for a fixed flexmix.

### Parameters

`mixtureModel(sce, model_type = "linear")`; `filterCells(sce, model = NULL, posterior_cutoff =
0.75, keep_all_below_boundary = TRUE, enforce_left_cutoff = TRUE)` (with `model = NULL` it fits
the linear model itself, with a warning). [SRC]

### Output

The filtered `SingleCellExperiment`; inside `filterCells` a per-cell `prob_compromised` and
`keep` (the returned object carries only the kept cells; `plotFiltering` shows both). [SRC]

### Oracle

R with Bioconductor miQC (bioconda `bioconductor-miqc` 1.18.0, or 1.22.0 from Bioc 3.24), flexmix
and scater for the QC columns; `set.seed()` before `mixtureModel`. Exact reproduction needs
flexmix's EM (Gaussian regression M-step with weighted least squares, the `minprior` pruning)
and R's `sample` for the initial labels; the repo has R's Mersenne-Twister and `sample.int`
would need R's rejection sampling (R ≥ 3.6). Feasible: the model is small. Alternatively give
the original's fitted parameters as a fixture and check steps 2–5 exactly.

## 8. dropEst quality classifier (dropestr)

### Identity

- Paper: Petukhov V, Guo J, Baryawno N, et al. dropEst: pipeline for accurate estimation of
  molecular counts in droplet-based single-cell RNA-seq experiments. *Genome Biology* 19:78
  (2018), doi:10.1186/s13059-018-1449-6 (PMID 29921301; resolves to this paper).
- Code: <https://github.com/kharchenkolab/dropEst> (formerly `hms-dbmi/dropEst`), R package
  `dropestr` in `dropestr/` (R + Rcpp; imports ks ≥ 1.10.4, pcaPP ≥ 1.9, Matrix, dplyr). bioconda
  `dropest` 0.8.5 and 0.8.6 (Linux only) builds the C++ pipeline and installs `dropestr` with
  `R CMD INSTALL`.
- Licence: `LICENSE.md` is GPL-3 for the repository; `dropestr/DESCRIPTION` says `License:
  GPL-3`. [SRC]
- Read at tag `v0.8.6` (commit 8f59b44, 2019-08-01; dropestr 0.7.9):
  `dropestr/R/low_quality_cells.R`, `dropestr/R/kde_classifier.R`,
  `dropestr/R/quality_control.R`, `scripts/report.Rmd`,
  `dropestr/vignettes/low-quality-cells.Rmd`, `configs/*.xml`. `master` (abcc42c, Sep 2020) only
  changes documentation after the tag. [SRC]

### Releases [SRC]

dropEst tags `v0.8.5` (Apr 2019) and `v0.8.6` (Aug 2019) both ship dropestr 0.7.9, with the same
classifier. Earlier changes, untagged: 72e6428 (Jan 2018) the current scoring scheme; 0.7.7
(eaf0411, Mar 2018) clamps negative KDE values to 0 and skips the intergenic filter for
pseudoaligners; d994488 (Feb 2019) makes `bandwidth.mult` take effect (default 1, so no change
at the default). One variant (dropestr 0.7.9). The bounds that label the training sets are the
repo's `dropest-bounds` method (same tag).

### Algorithm (`ScorePipelineCells`, dropestr 0.7.9) [SRC]

Input: dropEst's `cell.counts.rds` list: `cm_raw` (UMI counts for barcodes with ≥
`min_genes_before_merge` genes, 20 in the shipped configs), `cm` (after barcode merging, ≥
`min_genes_after_merge` = 100 genes), `aligned_reads_per_cell`, `aligned_umis_per_cell`,
`reads_per_chr_per_cells$Exon` / `$Intergenic` (reads per chromosome and region type),
`merge_targets`; optionally droptag's `reads_per_cb` (all reads, aligned or not).

1. Training labels: barcodes of `cm_raw` sorted by UMIs; `EstimateCellsNumber` bounds give
   `High` (rank ≤ min), `Low` (rank ≥ max), `Unknown` (between).
2. Features per barcode (`PrepareLqCellsData`): `ReadsPerUmi` = aligned reads / UMIs;
   `UmiPerGene` = UMIs / genes; `LowExpressedGenesFrac` = genes with exactly 1 UMI / genes;
   `IntergenicFrac` = intergenic reads / (intergenic + aligned reads); with `reads_per_cb`,
   `NotAlignedUmisFrac` = max(total − aligned reads, 0) / total reads / ReadsPerUmi (totals
   summed over merged barcodes); with mito genes or the mito chromosome name,
   `MitochondrionFraction` (UMI fraction from genes, or exonic-read fraction from the
   chromosome). Features that are ~0 everywhere are dropped; each is scaled to [0, 1] by
   $$(x - \min)/\max(x - \min)$$.
3. Robust sparse PCA (`GetOptimalPcs`, `max.pcs = 3`): `pcaPP::sPCAgrid(scaled, k = #features)`
   (pcaPP defaults: `method = "mad"`, `lambda = 1`, `maxiter = 10`, `splitcircle = 25`, L1-median
   centring; a deterministic grid search); keep the fewest PCs reaching 98 % of explained variance,
   at most 3; "used features" = those with normalised |loading| > 0.075 on a kept PC; scores
   rescaled to [0, 1].
4. Mitochondrial (if given) and intergenic (if present) outliers: fraction >
   `mean(x, trim = 0.2) + 4 · mad(x)` (or a given threshold). If the feature is among the used
   features, those barcodes are relabelled `Low` before training; otherwise their final score
   is set to the minimum score.
5. Training (`TrainClassifier`): `High` = class 1, `Low` = class 0; if there are more than
   1.5 × #High Low barcodes, keep the round(1.5 × #High) Low barcodes with the most UMIs. Priors
   0.5/0.5. Per class, the normal-scale bandwidth `ks::Hns` $$= (4/(n(d+2)))^{2/(d+4)}\,\Sigma$$.
6. Score (`PredictKDE`): for the barcodes of `cm` (or all of `cm_raw` with `predict.all =
   TRUE`), class densities `ks::kde(class data, H = H × bandwidth.mult, eval.points = x)`
   (binned on ks's default grid when n > 500 in 2-D or > 1000 in 3-D, then interpolated; exact
   otherwise), clamped at 0; score $$= f_1/(f_1 + f_0)$$; 0.5 where $$f_0 + f_1 < 10^{-10}$$.
7. Call: none inside the function. dropEst's report (`scripts/report.Rmd`) and the vignette count
   **score > 0.9** as cells ("Cell scores (N cells > 0.9)").

`ScoreQualityData(umi.counts, scoring.data, ...)` is the same classifier on a user-supplied
feature table (scaled, PCA only if `max.pcs` < #features), without the outlier filters.

Determinism: no random numbers (sPCAgrid's grid search, Hns and KDE are deterministic). Result
depends on the ks version (binning grid, `predict` interpolation). [SRC]

### Parameters

`ScorePipelineCells(pipeline.data, mitochondrion.genes = NULL, mit.chromosome.name = NULL,
tags.data = NULL, filter.mitochondrial = NULL (TRUE if mito info given), filter.intergenic =
TRUE, mit.fraction.threshold = NULL, intergenic.fraction.threshold = NULL, max.pcs.number = 3,
predict.all = FALSE, kde.bandwidth.mult = 1, cell.number = NULL)`; the 0.9 cut is a convention
of the report. [SRC]

### Output

A named vector of scores in [0, 1] (probability of the `High` class) per barcode of `cm` (or
`cm_raw`), in decreasing-UMI order. [SRC]

### Oracle

R with `dropestr` from bioconda `dropest` 0.8.6 (Linux), or the R files sourced from the
`v0.8.6` tag as the repo's `dropest-bounds` oracle does, with conda-forge `r-ks` and `r-pcapp`.
Inputs beyond the UMI matrix: aligned reads per barcode (Cell Ranger's `molecule_info.h5` read
counts per molecule summed per barcode give this, or STARsolo's per-barcode read statistics) and,
for `IntergenicFrac`, intergenic reads per barcode (from a BAM, or STARsolo `CellReads.stats`);
without reads only `UmiPerGene` and `LowExpressedGenesFrac` can be computed (passing UMIs as reads
gives a constant `ReadsPerUmi`, which is not dropped, as only all-zero features are, and makes the
[0, 1] scaling divide by zero), so such inputs go through `ScoreQualityData` with those two
features. Exact reproduction needs `sPCAgrid` (pcaPP's C++ grid search) and ks's
binned KDE ported; a staged oracle (fixed PCs as input, then the KDE scores exactly) is the
practical route.

## 9. salmon alevin naive-Bayes whitelist

### Identity

- Paper: Srivastava A, Malik L, Smith T, Sudbery I, Patro R. Alevin efficiently estimates
  accurate gene abundances from dscRNA-seq data. *Genome Biology* 20:65 (2019),
  doi:10.1186/s13059-019-1670-y (PMID 30917859; resolves to this paper).
- Code: <https://github.com/COMBINE-lab/salmon>, C++: `src/WhiteList.cpp` (the classifier),
  `src/CollapsedCellOptimizer.cpp` (per-barcode features, the call into the classifier),
  `src/Alevin.cpp` (the initial whitelist and the low-confidence band, already the repo's
  `alevin-initial-whitelist`). bioconda `salmon` 0.11.0 – 1.10.3 with alevin.
- Licence: GPL-3 (`LICENSE` at `v1.10.3`, "GNU GENERAL PUBLIC LICENSE Version 3"). [SRC]
- Read at tags `v0.11.0`, `v0.11.3`, `v0.12.0`, `v0.13.1`, `v0.14.0`, `v0.14.1`, `v1.0.0`,
  `v1.1.0`, `v1.4.0`, `v1.10.3` (`src/WhiteList.cpp`), and `v0.14.0`, `v1.0.0`, `v1.1.0`,
  `v1.10.3` (`src/CollapsedCellOptimizer.cpp`, `include/SalmonDefaults.hpp`). [SRC]

**Correction** to [other-methods-and-reviews.md](other-methods-and-reviews.md): alevin was not
removed by salmon 2.0. It was removed in **salmon 1.11.0** (March 2026): from `v1.11.0` the tree
has `src/commands/AlevinRemoved.cpp` (the command prints that it was removed and points to
alevin-fry, or "Salmon v1.10.2" for the legacy implementation) and no `src/WhiteList.cpp`.
1.11.x (March 2026) and 1.12.x (12 and 22 June 2026) are bulk-only C++ releases; 2.0.0 (13 June
2026) is the Rust rewrite. The last release with alevin is **1.10.3** (March 2024). [SRC]

### Releases [SRC]

`WhiteList.cpp` has three behaviours:

| salmon | Features | Notes |
|---|---|---|
| 0.11.0 – 0.13.1 | computed in memory from the quantified matrix: mapping rate, mean/max gene count with the mean over **all** genes, deduplication rate, genes over that mean (4), + mito, ribo fraction (if `--mrna`/`--rrna` files exist), + maximum correlation with the true set (`--useCorrelation`, off by default) | 0.11.0 → 0.11.3 only rewrites the prediction loop (same result) |
| 0.14.0 – 1.0.x | read back from `featureDump.txt` (written with C++'s default 6 significant digits): mapping rate, deduplication rate, mean/max with the mean over **expressed** genes, expressed genes, genes over that mean, mean molecules per UMI-graph arborescence (6), + mito, + ribo | |
| 1.1.0 – 1.10.3 | as 0.14 without the arborescence feature (5), + mito, + ribo | `numFeatures {6}` → `{5}` in 1.1.0; nothing else changes through 1.10.3 |

The classifier itself (`naive_bayes_fit`, `get_log_likelihood`, `classifyBarcodes`) is the same
code from 0.11.0 to 1.10.3. A pipeline still running alevin today most likely runs 1.10.x
(bioconda's last is 1.10.3; the removal message points to 1.10.2).

### Algorithm (1.1.0 – 1.10.3) [SRC]

Input: the initial whitelist from `Alevin.cpp` (knee set, or expected/forced cells, top $$k$$ by
corrected read count) followed by the low-confidence band below it ($$\lfloor 0.5k\rfloor$$
clamped to [`--lowRegionMinNumBarcodes` = 200, 1000]), in decreasing count order (`trueBarcodes`,
$$N$$ barcodes, $$L$$ of them in the band). Barcodes in this list with no mapped read are dropped
first (and $$L$$ is reduced for those at index $$>$$ $$N - L$$; one at exactly $$N - L$$ is not
counted, an off-by-one).

Skipped (no classifier; every barcode in the list is kept) when `--whitelist` is given,
`--noWhitelist`, hash/naive-eqclass modes, $$L$$ < 200, or $$N - L$$ < 90.

1. Features per barcode, from the quantification (`CollapsedCellOptimizer`): with raw reads
   $$r$$ (corrected-barcode read count), mapped reads $$m$$, deduplicated UMIs $$u = \sum_g
   \alpha_g$$ (alevin's per-gene estimates $$\alpha_g$$), expressed genes $$G = \#\{\alpha_g >
   0\}$$:
   mapping rate $$m/r$$; deduplication rate $$1 - u/m$$; mean-by-max $$(u/G)/\max_g \alpha_g$$;
   $$G$$; genes with $$\alpha_g > u/G$$; and, if the `--mrna` / `--rrna` gene lists give more than
   one gene, mito and ribo fractions $$\sum_{g \in M} \alpha_g / u$$. Each is printed with 6
   significant digits to `featureDump.txt` and parsed back with `strtod`. (With exactly one
   listed gene the classifier expects the column but the writer omits it: alevin exits with
   code 84.)
2. Classes by position: true = the first $$T = \lfloor (N - L)/2 \rfloor$$ barcodes (top half of
   the initial whitelist); ambiguous = the next $$N - L - T$$; false = the last $$L$$ (the band).
3. Gaussian naive Bayes fit, with two off-by-ones: the true class uses rows $$0..T-2$$ and the
   false class rows $$N-L..N-2$$ (each range drops its last barcode), but both are divided by the
   full class size; per feature $$\mu = \sum x / n_c$$ and $$\sigma^2 = \sum (x - \mu)^2 / n_c$$
   over the rows used. Priors `classCount / norm` are integer divisions, so both are **0**, and a
   zero prior contributes $$\log$$ prior = 0: the classes are compared without priors.
4. Score per class: $$-\tfrac12\left(0.5 + \sum_{j:\,\sigma_j^2 \ne 0} \left[\log(2\pi\sigma_j^2)
   + (x_j - \mu_j)^2/\sigma_j^2\right]\right)$$ (a feature with zero variance in that class is
   skipped in that class only).
5. Whitelist = the $$T$$ true barcodes + every ambiguous barcode whose true score ≥ false score
   (ties go to true). The band is never kept. Written to `whitelist.txt` in list order;
   `--dumpFeatures` also writes `predictions.txt` (the two scores of the kept ambiguous
   barcodes).

0.11 – 0.13 differ only in the features (step 1): mean-by-max and genes-over-mean use the mean
over all genes of the reference, the features are not rounded, there is no expressed-genes
feature, and `--useCorrelation` adds each barcode's maximum Pearson correlation with the true
set. 0.14 – 1.0 add the arborescence feature (mean number of molecules per UMI-graph
arborescence).

Determinism: the classifier has no randomness. Upstream, the hash-map order of tied barcodes (which
sets the list order and so the class boundaries) and alevin's EM determine the features. One more
hazard, read in the source but not observed in a run (UNVERIFIED in effect): `featureDump.txt` is
written by the worker threads in completion order (`GZipWriter::writeSparseAbundances`, under a
mutex), while `performWhitelisting` ignores its CB column and takes row $$i$$ as `trueBarcodes[i]`.
With one worker (`-p` ≤ 2) the orders agree; with more, neighbouring rows can swap, so a barcode's
features can be classified and named as its neighbour's, and the calls near the class boundaries can
change from run to run.

### Parameters

`--mrna FILE`, `--rrna FILE` (gene lists that add features), `--lowRegionMinNumBarcodes` (200),
`--noWhitelist`, `--whitelist FILE`, `--dumpFeatures`; plus everything that shapes the initial
whitelist (`--expectCells`, `--forceCells`, `--keepCBFraction`, `--maxNumBarcodes`,
`--freqThreshold`). Defaults from `include/SalmonDefaults.hpp` at v1.10.3. [SRC]

### Output

`whitelist.txt` (the final cells); `featureDump.txt` (the features of every barcode in the list,
plus raw, mapped and deduplicated counts); `predictions.txt` with `--dumpFeatures`. The quantified
matrix (`quants_mat`) holds every barcode in the list; the whitelist marks the cells. [SRC]

### Oracle

Two levels. (a) The classifier alone: build a program from `WhiteList.cpp`'s functions copied out
of each tag (as `oracles/alevin/harness.py` does for `Alevin.cpp`), feed it a feature table, and
compare the whitelist exactly; our port must reproduce the 6-significant-digit round trip and the
off-by-ones. (b) End to end: alevin's features need raw reads per barcode (including unmapped),
mapped reads and alevin's own per-gene estimates, i.e. a run of salmon alevin 1.10.3 (bioconda)
on FASTQs; from a count matrix alone only mean-by-max, expressed genes, genes over the mean and
the mito/ribo fractions can be computed, and mapping and deduplication rates need read counts
(STARsolo's `CellReads.stats` or Cell Ranger's `molecule_info.h5` approximate them, not exactly).

## Corrections to other-methods-and-reviews.md

The eight DOIs it gives for these tools (DropSift has no paper) resolve to the right papers
(checked against PubMed); the licences in its section C and summary table match the
repositories. What needs correcting:

- Introduction: salmon **1.11.0** (March 2026) removed `salmon alevin`, not salmon 2.0 (13 June
  2026); the last release with alevin is 1.10.3 (March 2024); the last C++ release is 1.12.1
  (22 June 2026), not 1.12.0. Naive-Bayes whitelisting spans 0.11.0 – 1.10.3, not 0.14 – 1.x.
  (Section 9.)
- C6 DropletQC: the nuclear fraction comes from the BAM's `RE` tags or from read overlaps with a
  GTF/GFF (`nuclear_fraction_annotation`), not from velocyto; it counts alignment records
  (reads, duplicates included), not UMIs. The cut-off is taken on the KDE of the density's
  **derivative** (its own plug-in bandwidth): the end of the second sign run, which is the first
  minimum after the first peak only when the derivative starts positive. (Section 3.)
- C1 CB2: the knee iteration's step $$(knee - inflection)/10$$ is fixed from the first fit; the
  cluster p-values come from unseeded PSOCK workers, so CB2 is not reproducible under
  `set.seed`; clusters must also pass a tightness filter before testing. (Section 1.)
- C3 DIEM: the low-gene debris clusters are selected by position in alphabetical label order
  (a bug that flags the wrong clusters once there are ≥ 10); k-means runs with `nstart = 1` since
  2.4.0. (Section 2.)
- C8 DropSift: since June 2026 (#38) there is a third, debris class with two more SVMs and a
  debris gene module, and since September 2026 an automatic two-peak fallback; the repository has
  no tags. (Section 4.)
- C13: miQC is post-calling QC (compromised cells), as listed; its licence is BSD-3-Clause.

## Summary

| Tool | Category | Licence | Language | Deterministic? | Extra inputs beyond the raw UMI matrix | Oracle difficulty |
|---|---|---|---|---|---|---|
| CB2 / scCB2 1.0–1.24 | knee + clustered Monte Carlo test + EmptyDrops | GPL-3 | R | no: cluster test in unseeded PSOCK workers; the rest under `set.seed` | none | medium–high: noise band; exact needs R's `rbinom` and seeded workers |
| DIEM 2.4.1 | multinomial-mixture EM + debris score | GPL-3 | R + C++ | yes (seeds 1 itself; 1 thread) | none (mito/MALAT1 only for plots) | medium–high: `irlba`, Hartigan–Wong k-means, `loess`; staged oracle with fixed k-means labels |
| DropletQC 1.0.0 | KDE-derivative valley on nuclear fraction + rescue; damaged cells by bivariate GMM | MIT | R | yes | BAM with `RE` (or GTF + BAM); cell types for damaged cells | medium: port ks `hpi` + binned `kdde`; NF from a BAM |
| DropSift (ccaefb7) | KDE/HDR exemplar boxes + RBF SVMs with gene-module scores | MIT | R | no, unless the caller seeds | intronic and mito **read** fractions per barcode; optional CellBender `frac_contamination` | high: Seurat, e1071, hdrcde, `sample`; compare with tolerance or by stage |
| valiDrops 1.0.0 | segmented-regression knee + adaptive QC + cluster marker filter | GPL-3 | R | stage 1 yes (segmented ≥ 2.0-3); stages 2–3 no unless seeded | gene annotation (species auto-detected) | stage 1 medium (port `segmented`); stages 2–3 high |
| SiftCell 1.0 (adaaec1) | shuffled-UMI null + XGBoost on PCs | GPL ≥ 2 | R + C++ + Python | no by default (`seed = 0` sets nothing); yes under `set.seed` with fixed xgboost | none (uncompressed MTX directory) | medium–high: needs xgboost < 1.6; R `runif` stream, squat test |
| miQC 1.2–1.22 | two-component mixture of regressions (mito % ~ genes); post-calling QC | BSD-3-Clause | R | no, unless seeded (flexmix random start) | mito gene list; called cells | medium: small flexmix EM, R `sample` |
| dropEst classifier (dropestr 0.7.9) | KDE Bayes classifier on robust sparse PCs of QC features | GPL-3 | R | yes | aligned reads per barcode; intergenic reads; mito | medium–high: `sPCAgrid`, ks binned KDE; staged oracle |
| alevin naive Bayes (salmon 0.11–1.10.3) | Gaussian naive Bayes on per-barcode features | GPL-3 | C++ | the classifier yes; but with `-p` > 2 the feature rows follow thread completion order (section 9) | raw and mapped reads per barcode, alevin's gene estimates | low for the classifier (harness from `WhiteList.cpp`); high end to end (needs salmon on FASTQs) |
