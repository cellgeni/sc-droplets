# DropletUtils, CellBender and cellsweep: cell and empty-droplet calling by version

Read from tagged source, release branches and commits (GitHub and code.bioconductor.org,
`RELEASE_3_22`–`RELEASE_3_24`). UNVERIFIED marks what could not be confirmed.

Three findings that shape the package:

- DropletUtils `emptyDrops` used a Dirichlet-multinomial with an estimated `alpha` by default for
  1.2–1.26; only 1.0 and 1.28+ use the plain multinomial. Its knee definition changed three times
  (1.10, 1.28, 1.30), and the knee sets `retain`, so `emptyDrops` calls drift across versions with
  identical arguments.
- CellBender v0.3.0, v0.3.2 and v0.4.0 share byte-identical prior heuristics (diffs are type
  annotations only). v0.3.1 was never tagged (withdrawn). v0.3.0_rc, the paper's code, uses a
  pyro GMM instead.
- cellsweep is an ambient-RNA denoiser, not a cell caller. Its only knee-type component marks
  droplets as empty; that step changed in 1.0.2 (the `kneed` library), and its result also depends
  on the installed `kneed` version.

## 1. CellBender `remove-background`

Licence BSD-3-Clause (Broad Institute), every release.

| Tag | GitHub release | PyPI | Python | Key deps |
|---|---|---|---|---|
| v0.1.0 | 2020-09-25 | 0.1.0 (uploaded 2023-03-23) | 3.7 | `pyro-ppl>=0.3.2`; torch unpinned (conda) |
| v0.2.0 | 2022-12-09 | 0.2.0 | 3.7 | same |
| v0.2.1 | 2022-12-09 | 0.2.1 | 3.7 | + `anndata>=0.7` |
| v0.2.2 | 2022-12-09 | not on PyPI | 3.7 | same |
| v0.3.0_rc | tag only (2022-06-27); code of the Nature Methods paper | — | — | — |
| v0.3.0 | 2023-08-08 | 0.3.0 | 3.7 | `pyro-ppl>=1.8.4`, torch unpinned |
| v0.3.1 | never tagged; withdrawn (integer-overflow bug, PR #347) | — | — | — |
| v0.3.2 | 2024-07-17 | 0.3.2 (uploaded 2026-05-05) | classifier 3.7 | same |
| v0.4.0 | 2026-08-15 | 0.4.0 | `>=3.10`; CI 3.10–3.12 | `pyro-ppl>=1.8.4`, torch, `anndata>=0.9`, `dill` |

v0.4.0 changes no algorithm or output (changelog, release notes, and the diffs of `priors.py`,
`dataset.py`, `consts.py`).

Python API: no in-memory API in any version. `SingleCellRNACountsDataset(input_file=...)` always
loads a file; `run_remove_background(args)` reads and writes files. AnnData only as `.h5ad`
(v0.2.1+: `layers['counts']`, then `layers['spliced']` if `X` is None, else `X`). In v0.3+,
`cellbender.remove_background.data.priors.get_priors(umi_counts, low_count_threshold)` is plain
numpy/scipy but the module imports torch at the top.

### 1a. `cellbender_v0.1_prior` (v0.1.0)

Source: `cellbender/remove_background/data/dataset.py` at v0.1.0:
`_trim_dataset_for_analysis`, `_estimate_priors`, `get_d_priors_from_dataset`,
`estimate_cell_count_from_dataset`, `estimate_chi_ambient_from_dataset`; `consts.py`, `argparse.py`.

| Constant | Value |
|---|---|
| `EMPIRICAL_LOW_UMI_TO_EMPTY_DROPLET_THRESHOLD` | 0.8 |
| `LOW_UMI_CUTOFF` (CLI `--low-count-threshold`) | 15 (the `__init__` default 30 is overridden) |
| `TOTAL_DROPLET_DEFAULT` | 25000 |
| `CELL_PROB_CUTOFF` | 0.5 |
| `FRACTION_EMPTIES` | 0.5 |
| `SIMPLE_MODEL_D_STD_PRIOR` | 0.2 |

```
G   = genes with column sum > 0, minus --blacklist-genes
c   = row sums of X[:, G]                      # size priors
cA  = row sums of X (all genes), int           # n_cells estimate and trimming order
empty_counts = int(expm1( scipy.stats.mode(round(log1p(c[c > low_count_threshold]), 1))[0] ))   # ties -> smallest
cell_counts  = int(expm1( median(log1p(c[c > 5*empty_counts])) ))
if expected_cells is None:
    cut = expm1( mean(log1p(cell_counts), log1p(empty_counts)) )
    n_cells = sum(cA > cut)
n_cells = min(n_cells, sum(cA > 0))
order   = argsort(cA)[::-1]
low_cut = max(low_count_threshold, int(empty_counts * 0.8))
n_above = sum(cA > low_cut)                          # asserts n_above > n_cells
n_trans = min(total_droplets - n_cells, n_above - n_cells)
analyzed = order[:n_cells + n_trans]
empties  = order[n_cells+n_trans : n_above] if n_trans == total_droplets - n_cells else []
log_crossover = mean(log1p(cell_counts), log1p(empty_counts))
d_std     = (log1p(cell_counts) - log_crossover) / 5
cell_logit = logit((1 - fraction_empties) * n_cells / len(analyzed))
chi_ambient: among analyzed barcodes with np.log(total) < log_crossover  (log, not log1p; fixed in v0.2)
```

Model call (needs CellBender): `p = sigmoid(p_y logit)` for analyzed barcodes; cells `p > 0.5`.

### 1b. `cellbender_v0.2_prior` (v0.2.0 = v0.2.1 = v0.2.2)

Same functions as v0.1.0. Changes: `EMPIRICAL_LOW_UMI_TO_EMPTY_DROPLET_THRESHOLD = 0.5`; with
`--expected-cells` given, `cell_counts = int(median(top n_cells of c))`, else the v0.1 log-median
rule; chi_ambient empties use `log1p` (fix); new `empty_UMI_threshold =
sorted_desc(cA)[len(analyzed)]`; `TRAINING_FRACTION = 0.9`; `--fpr` default `[0.01]`;
`--exclude-antibody-capture`; `D_STD_PRIOR = 0.02`; learning rate 1e-4, `z_dim` 100; v0.2.1 h5ad
input; v0.2.2 retry options. Default `--low-count-threshold` 15.

Model call: training masks "surely empty" `counts < empty_UMI_threshold`, "surely cell"
`counts >= exp(log1p(cell_counts))`; prior `Bernoulli(logits=cell_logit - 100)`; cells `p > 0.5`.

### 1c. `cellbender_v0.3rc_gmm` (tag v0.3.0_rc) — docs only

```
c = c[c > low_count_threshold]
log_mode = mode(round(5*log(c))/5)
sd = distance in log space from the mode to where the histogram (100 bins) drops below 0.6 x peak
surely_empty = exp(log_mode + 2*sd); total_droplets = sum(c > surely_empty) + 5000
cut_low = exp(log_mode - 3*sd)
pyro GMM on log(c[c > cut_low]): 10 components, alpha prior 0.1, 500 epochs
crossover = mean of first 2 component locs; n_cells = #(c > exp(crossover)); cell_counts = their median
```

Needs torch/pyro, stochastic, never on PyPI.

### 1d. `cellbender_v0.3_prior` (v0.3.0 = v0.3.2 = v0.4.0)

Source: `cellbender/remove_background/data/priors.py` (`_threshold_otsu`, `_create_histogram`,
`_peak_density_given_cutoff`, `get_cell_count_empty_count`, `get_expected_cells_and_total_droplets`,
`get_priors`, `compute_crossover_surely_empty_and_stds`, `get_cell_count_given_expected_cells`,
`get_empty_count_given_expected_cells_and_total_droplets`) and `dataset.py`
(`SingleCellRNACountsDataset.__init__`, `_exclude_feature_types`, `_clean_features`,
`_trim_noiseless_features`, `_trim_droplets`, `_estimate_chi_ambient`).

| Constant | Value |
|---|---|
| `LOW_UMI_CUTOFF` (`--low-count-threshold`) | 5 |
| `EMPIRICAL_LOW_UMI_TO_EMPTY_DROPLET_THRESHOLD` | 0.5 |
| `MAX_TOTAL_DROPLETS_GUESSED` | 70000 |
| `MAX_EMPTIES_TO_INCLUDE` | 20000 |
| `NUM_EMPTIES_INCREMENT` | 20000 |
| `AMBIENT_COUNTS_IN_CELLS_LOW_LIMIT` | 0.1 |
| `FRACTION_EMPTIES` | 0.2 |
| `CELL_PROB_CUTOFF` | 0.5 |
| `--expected-cells` / `--total-droplets-included` / `--force-cell-umi-prior` / `--force-empty-umi-prior` | None |

```
features = drop --exclude-feature-types, zero-total features, --ignore-features
u = row sums over these features (all barcodes, zeros included)

get_cell_count_empty_count(u, low_count_threshold=L):
  s = sort(u)[::-1]
  otsu_input = s[: sum(u > L) // 4]
  cell_low = exp(Otsu(log(otsu_input)))   # np.histogram 256 bins; threshold = bin_centers[argmax(w1[:-1]*w2[1:]*(m1[:-1]-m2[1:])**2)]
  cell_prior = mean(u[u > cell_low])
  k = s[s > L]
  cutoff = max(0.1*cell_low, k[-100]);  empty_prior = -100
  while |delta| > 10 (at most 5 iterations):
      nonc = k[k <= cutoff];  n_c = sum(k > cell_low)
      idx = logspace(log10(n_c), log10(n_c + len(nonc)), 1000);  idx = clip(int(idx - n_c), 0, len(nonc)-1)
      v = log(sort(nonc)[::-1][idx])
      grid = arange(floor(min v) - .01, ceil(max v) + .01, 0.1)
      dens = scipy.stats.gaussian_kde(v).evaluate(grid)    # Scott bandwidth
      new_empty = exp(grid[argmax dens])
      upper = exp(first grid point at/after the peak with dens < 0.6*max), else last grid point
      cutoff = 0.75*sqrt(cell_prior*new_empty); delta = new_empty - empty_prior; empty_prior = new_empty
  cell_counts = median(u[u > sqrt(cell_prior*empty_prior)])
  return cell_counts, empty_counts = empty_prior, empty_count_upper_limit = upper

get_expected_cells_and_total_droplets:
  expected_cells = sum(u >= cell_counts)                   # deliberately conservative
  total_droplets = sum(u >= sqrt(empty_counts * empty_count_upper_limit))
  transition     = sum(u >= sqrt(cell_counts * empty_counts))
  total_droplets = min(total_droplets, transition + 20000)

get_priors: run both; while total_droplets > 70000 (at most 6 retries): rerun with L = empty_count_upper_limit

User overrides (SingleCellRNACountsDataset.__init__):
  expected_cells given: cell_counts = exp(mean(log(top-expected u)))
      if expected + 20000 > total_droplets: total_droplets = expected + 20000
          and, if total_droplets was not given:
          empty_counts = median(u_sorted_desc[max(expected, total-500) : total])
          empty_count_upper_limit = min(sqrt(cell_counts*empty_counts), 1.5*empty_counts)
  total_droplets given: set it and recompute empty_counts / upper limit as above
  force_cell_umi_prior -> cell_counts
  force_empty_umi_prior -> empty_counts; upper = min(sqrt(cell*empty), 2*empty)

compute_crossover_surely_empty_and_stds:
  surely_empty_counts = sort(u)[::-1][total_droplets]
  log_counts_crossover = (log(surely_empty_counts) + log(cell_counts)) / 2
  d_std = std(log u[u>0] above crossover) / 5;  d_empty_std = 0.01

_trim_noiseless_features: empties = (u < empty_count_upper_limit) & (u > L)
  keep feature g if mean_empties[g] * max(expected_cells, 5000) > 0.1
_trim_droplets (u on the final feature set):
  low_cut  = max(L, int(0.5 * empty_counts)); analyzed = top total_droplets barcodes
  empties  = ranks [total_droplets, #(u > low_cut)); empty_UMI_threshold = u_sorted_desc[total_droplets]
  cell_logit = logit((expected_cells / total_droplets) * (1 - 0.2))
```

Model part (requires CellBender): per-droplet prior logit −100 if
`lc <= (log(empty_UMI_threshold) + log1p(cell_counts)) / 2`, +10 if `lc >= log1p(cell_counts)`,
`cell_logit` otherwise; encoder blends toward ±10 outside the band; final
`cell_probability = sigmoid(p_y)` for the top `total_droplets`; cells `p > 0.5` (`run.py` uses `>`
for outputs, `>=` in `collect_output_metrics`).

Outputs: `out.h5` (`/matrix` CR v3 layout; `/droplet_latents`: `barcode_indices_for_latents`,
`gene_expression_encoding`, `cell_size`, `cell_probability`, `droplet_efficiency`,
`background_fraction`; `/global_latents`; `/metadata`), `out_filtered.h5`, `out_cell_barcodes.csv`,
`out_posterior.h5`, `out_report.html`, `out.pdf`, `out.log`, `out_metrics.csv`.

Reimplementable: everything in 1a, 1b, 1d, plus the deterministic hard-prior bands. Requires
CellBender: v0.3.0_rc GMM priors, `cell_probability`, denoised counts. None of the heuristics is
CellBender's cell call; `expected_cells` is intentionally conservative.

## 2. cellsweep (`github.com/pachterlab/cellsweep`)

An EM fit of a multinomial mixture for ambient plus "bulk" contamination (a denoiser): cell type
`(1-β)(1-α_n)·p_{k,g}`, ambient `(1-β)·α_n·a_g`, bulk `β·m_g`. Requires `adata.obs[celltype_key]`.
Empty droplets fixed at α=1 anchor the ambient profile. Licence BSD-2-Clause. Paper: Caskey, Rich,
Weber, Mortazavi, Pachter, Hallgrímsdóttir, "Single-Cell Genomics Decontamination with CellSweep",
bioRxiv 2026-03-06, doi:10.64898/2026.03.04.709349.

API: `cellsweep.denoise(adata_or_path, ..., empty_droplet_method, umi_cutoff, expected_cells,
random_state=42, ...)` (named `denoise_count_matrix` before 1.0.0); CLI `cellsweep denoise`;
`kb sweep` pins `cellsweep==1.0.2`.

No git tags; PyPI versions → commits: 0.1.0 `ce249907` (2026-03-06), 0.1.1 `41243ec0`/`2e0ee6ad`,
0.1.2 `faa3a6ad`, 0.1.3 `c33474fe`, 1.0.0 `3e12c260` (2026-10-07), 1.0.1 `ed8083a8`, 1.0.2
`c81ab360` (2026-11-10). 1.0.2 adds `kneed>=0.8.5`.

Empty-droplet rules (`cellsweep/utils/data_utils.py`, `infer_empty_droplets`; empty = total < cutoff):

- **autoknee** (0.1.0 → 1.0.1; method `"threshold"`): `umi_cutoff` if given; else
  `sort(totals)[::-1][expected_cells-1]`; else `automatic_umi_cutoff_detection(min_counts=10)`:
  ```
  c = sort(totals)[::-1]; c = c[c > 10]; x = arange(1, n+1)
  ys = scipy.ndimage.gaussian_filter1d(log10(c), sigma=2)
  d2 = np.gradient(np.gradient(ys, x), x); cutoff = c[argmin(d2)]
  ```
- **kneed** (1.0.2; default `"kneed"`; `umi_cutoff`/`expected_cells` force `manual`):
  ```
  t = totals[totals > 10]; umis, n = unique(t, return_counts=True)
  umis = umis[::-1]; ranks = cumsum(n[::-1])
  cutoff = KneeLocator(ranks, umis, curve="convex", direction="decreasing", online=True).knee_y
  ```
  kneed 0.8.6 (2026-03-20) added `detection_active` gating in `find_knee`, which can change
  `online=True` results versus 0.8.5; kneed 0.9.0 mirrors by x value.
- **mx_filter** (all versions; Python API only before 1.0.2), from the cellatlas/mx project:
  ```
  x = sort(totals)[::-1]; v = log1p(x)
  GaussianMixture(n_components=2, random_state=42).fit_predict(v)
  cutoff = x[argmax(entropy(predict_proba(v), axis=1))]
  ```

## 3. Bioconductor DropletUtils

Licence GPL-3. Paper: Lun et al., Genome Biology 2019, doi:10.1186/s13059-019-1662-y.

Versions ↔ Bioconductor: 1.0 (3.7), 1.2 (3.8), 1.4 (3.9), 1.6 (3.10), 1.8 (3.11), 1.10 (3.12),
1.12 (3.13), 1.14 (3.14), 1.16 (3.15), 1.18 (3.16), 1.20 (3.17), 1.22 (3.18), 1.24 (3.19), 1.26
(3.20), 1.28 (3.21), 1.30 (3.22), 1.32 (3.23), 1.34 (3.24, current release); devel 1.35.

Preamble of every `barcodeRanks` version (`R/barcodeRanks.R`):

```
totals = colSums(m)          # as.integer(round()) from 1.12 via .intColSums
o = order(totals, decreasing); rle over sorted totals: run.totals, lengths
mid-rank = cumsum(len) - (len-1)/2
y = log10(run.totals[run.totals > lower]); x = log10(mid-rank[...])     # lower = 100
```

- **barcoderanks v1** (1.0–1.8; `df=20`, `fit.bounds=NULL`):
  ```
  d1n = diff(y)/diff(x); right = which.min(d1n); left = which.max(d1n[1:right])
  inflection = 10^y[right]
  idx = left:right   (or y within log10(fit.bounds))
  fit = smooth.spline(x[idx], y[idx], df=20); curv = d2/(1+d1^2)^1.5
  knee = 10^y[which.min(curv)]   # BUG: indexes the full y, not y[idx]
  ```
  The knee lands `left-1` points too early. Fixed in commit `584320df` (1.10).
- **v2** (1.10–1.26): `exclude.from=50`; `skip = min(len(d1n)-1, #(x <= log10(exclude.from)))`;
  drop the first `skip` derivatives, offset indices back; `knee = 10^y[idx][which.min(curv)]`; if
  `len(idx) < 4`, `knee = 10^y[idx[1]]`.
- **v3** (1.28; commit `407cf2a6`): same edges and inflection; `df` ignored;
  ```
  chord from (x,y)[left] to (x,y)[right]; above = points with y >= chord
  knee = 10^y of the above-point with max |g·x - y + b|/sqrt(g²+1)
  ```
- **v4** (1.30–1.34; commit `8f18cce6`; `window=1`, `gradient.threshold=-1`):
  ```
  keep = run.totals > lower & mid-rank > exclude.from(50)     # need >= 2 points
  s = cumulative arc length along (x,y); windows start at each point with s+1 <= s_max
  right end and midpoint (s+0.5) linearly interpolated (findInterval left.open=TRUE)
  per window: gap = chord length, grad = chord slope, mid_above = mid.y > chord at mid.x
  elbow windows = !mid_above & grad < -1
  if none: infl = argmin grad; knee candidates = all mid_above
  else:    e = first elbow window; knee candidates = mid_above among windows[1:e]
           infl = argmin grad over windows[1 : findInterval(mid.x[e], left.x)]
  knee = 10^mid.y[candidate with min gap]   (fallback: the inflection window)
  inflection = 10^mid.y[infl]; if no window fits: both = 10^y[last]
  ```

`emptyDrops` / `testEmptyDrops` (`R/emptyDrops.R`, `R/ambientProfileEmpty.R`,
`src/montecarlo_pval.cpp`), every version:

```
drop genes with rowSums(m)==0 over ALL barcodes
totals = as.integer(round(colSums(m))); ambient = totals <= lower(100)
prof = rowSums(m[, ambient]);  a = edgeR::goodTuringProportions(prof)
test set = !ambient (test.ambient=TRUE: totals > 0) & totals > ignore
LogProb = sum_g [x log a_g - lfactorial(x)] + lfactorial(T)                          # multinomial
     or   sum_g [lgamma(αa_g+x) - lgamma(αa_g) - lfactorial(x)] + lfactorial(T) + lgamma(α) - lgamma(T+α)   # DM
Monte Carlo (niters=10000): each iteration draws ONE sequential sample path from a (DM: p_g ~ Gamma(α·a_g,1)),
  reused for all observed totals in increasing order; n.above = #iterations with sim_logp <= obs_logp
  PValue = (n.above+1)/(niters+1); Limited = (n.above == 0)
emptyDrops: retain = knee of barcodeRanks unless given; p'[Total >= retain] = 0; FDR = p.adjust(p', "BH")
```

| DropletUtils | emptyDrops behaviour | Default α | Knee for `retain` |
|---|---|---|---|
| 1.0 | as above; ambient p-values enter BH when test.ambient | multinomial | v1 |
| 1.2 | `alpha=NULL` → DM MLE via `optimize(c(0.01,10000))` over ambient columns | estimated | v1 |
| 1.4 | pcg32 seeds per iteration from the R RNG (`dqrng`); stop on empty ambient profile | estimated | v1 |
| 1.6 | Boost `discrete_distribution`; auto-round non-integer matrices | estimated | v1 |
| 1.8 | `round=`; `.safe_good_turing` (pseudo = 1/sum(prof) over zeros, rest × (1−pseudo)) | estimated | v1 |
| 1.10 | `by.rank=` (lower = totals[order(desc)[by.rank+1]]); original lower still passed to barcodeRanks (bug) | estimated | v2 |
| 1.12 | `test.ambient=TRUE`: Total ≤ lower get FDR NA (excluded from BH); `NA` restores old; by.rank lower passed to barcodeRanks | estimated | v2 |
| 1.14–1.18 | refactor; no behaviour change found | estimated | v2 |
| 1.20 | `known.empty=`; by.rank: `rank(-totals, ties="first") > by.rank` ambient, lower = max(totals[ambient]) | estimated | v2 |
| 1.28 | `alpha=Inf` default | multinomial | v3 |
| 1.30–1.34 | unchanged | multinomial | v4 |

Other functions:

- **emptyDropsCellRanger** (1.14+, port of STARsolo 2.7.9a `EmptyDrops_CR`): `n.expected.cells=3000`,
  `max.percentile=0.99`, `max.min.ratio=10`, `umi.min=500`, `umi.min.frac.median=0.01`,
  `cand.max.n=20000`, `ind.min=45000`, `ind.max=90000`, `niters=10000`, `alpha=Inf`.
  ```
  retain    = max(round(totals_desc[round(3000*(1-0.99))] / 10), 1)
  ncs       = #(totals >= retain)
  min.umi   = max(500, round(0.01 * totals_desc[ncs/2]))       # R truncates the index
  cand_last = min(ncs + 20000, #(totals > min.umi))
  ambient   = rank positions 45000..90000 (descending); tested = (ncs+1)..cand_last
  FDR via BH with Total >= retain set to 0
  ```
  Bug in 1.14.x, 1.16.x, 1.18.0: used the unsorted `totals[ncs/2]`; fixed in 1.18.1 (`49071509`).
- **defaultDrops** (1.0 → present): `thr = quantile(libs[top 3000], 0.99, type=7) * 0.1`; `cell = libs > thr`.
- **ambientProfileEmpty** (renamed `estimateAmbience`, 1.12); **ambientProfileBimodal** (renamed
  `inferAmbience`, 1.10/1.12: per feature `kmeans(log1p(x), centers=quantile(c(0.05,0.95)))`);
  **cleanTagCounts** (1.12).

Porting traps for exact equality: R `smooth.spline(df=20)` (df→λ search, knot rule; GCV fallback
if df > unique x); edgeR's SGT C++; R `optimize` (Brent, tol ≈ 1.22e-4); `quantile` type 7;
half-to-even `round`; `p.adjust` NA handling; the Monte Carlo RNG (pcg32 + Boost) cannot be
bit-matched — validate p-values distributionally.

## 4. Existing Python ports

| Package | What it implements | Fidelity | Licence |
|---|---|---|---|
| `cellbender` 0.4.0 | the real thing; `priors.get_priors` callable on arrays | exact for v0.3/0.4 priors | BSD-3 |
| `emptydrops` 0.0.5 (nh3/emptydrops) | Cell Ranger 3.0.2 `find_nonambient_barcodes` + SGT | Cell Ranger, not DropletUtils | MIT (port of the MIT-era 3.0.2 code) |
| `cell-filter` 0.1.2 (ArcInstitute) | STARsolo/CR-style EmptyDrops, own SGT | not faithful by default (n_expected 20000, max_percentile 0.95, ratio 5, FDR 0.01) | MIT |
| `cytearc` 1.4.0 | `barcode_ranks` / `empty_drops` rebuilt from the paper | matches no single DropletUtils version | BSD-3 |
| `kneed` 0.9.0 | Kneedle | changed in 0.8.6 and 0.9.0 | BSD-3 |
| `dropkick` 1.2.8 | multi-Otsu labels + glmnet logistic regression | different algorithm | MIT |
| scanpy, pegasus, scvi-tools | no EmptyDrops or knee caller | — | — |

## Timeline

| Date | Tool / version | Change |
|---|---|---|
| 2018-05 | DropletUtils 1.0 | barcodeRanks v1; emptyDrops multinomial; defaultDrops |
| 2018-10 | DU 1.2 | `alpha`; DM MLE default |
| 2019-05 | DU 1.4 | pcg32 seeding |
| 2019-10 | DU 1.6 | Boost sampler; rounding |
| 2020-04 | DU 1.8 | safe Good-Turing |
| 2020-09 | CellBender v0.1.0 | mode/median priors; factor 0.8; total 25000 |
| 2020-10 | DU 1.10 | barcodeRanks v2; `by.rank` |
| 2021-05 | DU 1.12 | test.ambient FDR change |
| 2021-10 | DU 1.14 | emptyDropsCellRanger (median bug until 1.18.0) |
| 2022-12 | CellBender v0.2.x | factor 0.5; median prior with expected cells |
| 2023-04 | DU 1.20 | known.empty; by.rank semantics |
| 2023-08 | CellBender v0.3.0 | Otsu + KDE priors; LOW_UMI 5 |
| 2024-07 | CellBender v0.3.2 | same heuristics |
| 2025-04 | DU 1.28 | barcodeRanks v3 (chord); α = Inf |
| 2025-10 | DU 1.30 | barcodeRanks v4 (curve tracing) |
| 2026-03 | cellsweep 0.1.x | gaussian_filter1d knee; mx_filter |
| 2026-08 | CellBender v0.4.0 | Python 3.10+; outputs unchanged |
| 2026-10/11 | cellsweep 1.0.0–1.0.2 | `denoise` rename; `kneed` default (1.0.2) |
