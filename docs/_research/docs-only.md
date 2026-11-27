# Docs-only methods

Source-reading notes, November 2026, for the methods that will be documented but not implemented:
CellBender's model as a cell caller (what comes after the priors on
[`cellbender/prior.md`](../methods/cellbender/prior.md)), EmptyNN, dropkick, BD Rhapsody's cell
calling, ScaleBio's ScaleRna, Singleron's CeleScope, Ambimux and AmbientMapper / scifi-demux. They
go deeper than B11 and C of [other-methods-and-reviews.md](other-methods-and-reviews.md) and the
CellBender part of [dropletutils-cellbender-cellsweep.md](dropletutils-cellbender-cellsweep.md),
and correct them where they were wrong (each correction is marked **Correction**). Tags as in the
[README](README.md): **[SRC]** read in the source at the tag or commit named, **[DOC]**
documentation or paper only, **[OUT]** inferred from published outputs, **UNVERIFIED** inferred.

Each section: identity, releases that changed the call, the method, parameters, output per
barcode, and why it is documented rather than implemented (and whether the pipeline's own output
can still be compared with other methods).

## 1. CellBender `remove-background`: the model's cell call

What happens after the priors on [`cellbender/prior.md`](../methods/cellbender/prior.md). Read in
`broadinstitute/CellBender` at tags v0.1.0, v0.2.2, v0.3.0_rc, v0.3.0, v0.3.2 and v0.4.0
(`model.py`, `vae/encoder.py`, `run.py`, `posterior.py`, `train.py`, `data/dataprep.py`,
`consts.py`, `argparser.py`; `infer.py` and `data/dataset.py` before 0.3.0).

### Identity

- Paper: Fleming, Chaffin, Arduini, …, Marioni, Philippakis, Ellinor, Babadi. "Unsupervised removal
  of systematic background noise from droplet-based single-cell experiments using CellBender."
  *Nat Methods* 20, 1323–1335 (2023), doi:10.1038/s41592-023-01943-7 (PMID 37550580; checked).
  The paper's code is the tag v0.3.0_rc (changelog: "The code run in the paper is technically
  v0.3.0_rc"). [DOC]
- Repository: <https://github.com/broadinstitute/CellBender>. Licence: BSD-3-Clause, file
  `LICENSE` ("Copyright (c) 2019, Broad Institute, Inc."). [SRC]
- Language: Python on PyTorch and Pyro. Distribution: PyPI `cellbender` (0.1.0, 0.2.0, 0.2.1,
  0.3.0, 0.3.2, 0.4.0), a GPU Docker image `us.gcr.io/broad-dsde-methods/cellbender`, and a WDL
  workflow for Terra (`wdl/`). [SRC] [DOC]

### Releases that changed the call

Tags and dates are in [dropletutils-cellbender-cellsweep.md](dropletutils-cellbender-cellsweep.md)
§1. What changed in the part after the priors:

| Release | Cell-probability network | Prior on the cell indicator $$y$$ | Supervision of the logit | Seed | Training defaults |
|---|---|---|---|---|---|
| 0.1.0 | own MLP for $$p$$ (`--p-layers 100 10`): input log total, normalised counts and counts minus ambient; initialised as a sigmoid of the log total around the crossover | `Bernoulli(logits=cell_logit)`, the same for every droplet | none | 0 | 150 epochs, z 20, `--z-layers 500`, Adam lr 1e-3 constant, empties 0.5 of each batch, NB likelihood |
| 0.2.0–0.2.2 | one MLP (`[100, 50]`) for $$p$$, $$d$$, $$\epsilon$$; input log total, log genes detected, a scaled Poisson log-likelihood under the ambient profile, scaled normalised counts | `Bernoulli(logits=cell_logit − 100)` for every droplet (marked `# TODO`) | logit $$\sim N(-5, 1)$$ for surely empty (count < `empty_UMI_threshold`), weight 1; $$N(5, 1)$$ for surely cell (count ≥ cell prior), weight 10 | 0 | 150 epochs, z 100, `--z-layers 500`, lr 1e-4 with OneCycleLR to 1e-3, empties 0.5, NB–Poisson convolution |
| 0.3.0_rc | as 0.2 (output scale 2; offset at the 25th percentile of the cell-like outputs) | as 0.2 | $$N(\mp 5, 1)$$ weight 0.1 for below / above the crossover; $$N(-5, 1)$$ weight 10 for surely empty; $$N(\mp 5, 0.2)$$ weight 10 / 1 when a surely-empty droplet has logit ≥ 0 or a surely-cell droplet logit ≤ 0 | 1234 | z 100, lr 1e-4 OneCycleLR, empties 0.2; priors from a pyro GMM (§1c of the older notes) |
| 0.3.0 = 0.3.2 = 0.4.0 | two 512-unit layers with batch norm and 50 % dropout; four extra features; a signed square-root term in the log total; hard blends to ±10 outside the band (below) | per droplet: −100 below the band, `cell_logit` inside, +10 above | $$N(\mp 10, 1)$$ weight 10 for below / at or above the crossover | 1234 | 150 epochs, z 64, `--z-layers 512`, lr 1e-4 OneCycleLR to 1e-3, empties 0.2, batch $$\min(256, 0.9T/2)$$ |

All [SRC]. v0.3.0 → v0.3.2 changes only `train.py`'s ELBO-failure checks (used only with
`--epoch-elbo-fail-fraction` / `--final-elbo-fail-fraction`); v0.3.2 → v0.4.0 is formatting and
types in the model files (checked with `git diff -w`). The changelog says the same ("does not
introduce any changes in the algorithm or the output"). [SRC] [DOC]

### The method (0.3.0 and later)

Inputs: the raw matrix over every feature not excluded (`--exclude-feature-types`,
`--ignore-features`), minus features whose projected ambient count in cells is below
`--projected-ambient-count-threshold` 0.1; the barcodes analysed are the top $$T$$ by total
(`total_droplets`, from the priors). Notation: $$n$$ is a droplet's total over the analysed
features, $$s = \log(1 + n)$$, $$C$$ the cell count prior, $$t$$ = `empty_UMI_threshold` (the total at
rank $$T$$, the largest surely-empty droplet), $$c = (\log t' + \log C)/2$$ the log crossover, with
$$t'$$ = `surely_empty_counts`, the same rank-$$T$$ total on the prior's feature set. [SRC]

Generative model (`model.py`, `--model full`, the default): per droplet a latent expression code
$$z \sim N(0, I_{64})$$ decoded to a profile $$\chi$$; cell size $$d_c \sim \mathrm{LogNormal}(\log(1+C), \sigma_d)$$;
empty size $$d_e \sim \mathrm{LogNormal}(\log(1+E), 0.01)$$ ($$E$$ the empty prior); efficiency
$$\epsilon \sim \mathrm{Gamma}(50, 50)$$; swapping fraction $$\rho \sim \mathrm{Beta}(1.5, 50)$$; cell
indicator $$y \sim \mathrm{Bernoulli}$$; counts from a negative-binomial (cell) plus Poisson
(background) convolution with mean $$(1-\rho)\,y\,\epsilon\,d_c\,\chi$$ and noise rate
$$\epsilon[(1-\rho)\chi_a d_e + \rho\,\bar\chi\,(y d_c + d_e)]$$, $$\chi_a$$ the ambient profile and
$$\bar\chi$$ the average profile. [SRC]

Prior on $$y$$ (`get_p_logit_prior`): logit −100 if $$\log n \le (\log t + \log(1+C))/2$$, +10 if
$$\log n \ge \log(1+C)$$, otherwise `cell_logit` $$= \mathrm{logit}\big((E_c/T)(1 - f)\big)$$ with
$$E_c$$ the expected cells and $$f$$ `--empty-drop-training-fraction` (0.2). The lower edge is
essentially the crossover $$c$$ (it uses $$\log(1+C)$$ where $$c$$ uses $$\log C$$). [SRC]

Variational posterior (`EncodeNonZLatents.forward`), the source of the cell probability:

1. Features: $$s$$; $$\log(1 + \text{genes detected})$$; an "overlap"
   $$-\tfrac12\big(\max(0, s - \mu_e)/0.1\big)^2$$ with $$\mu_e$$ the learnt `d_empty_loc`; and
   $$\lVert z_\text{loc} \rVert_2$$ from the $$z$$ encoder. Counts are normalised then
   `log1p`-transformed, batch-normalised and passed through `Dropout1d(0.5)`; two
   Linear(512)–BatchNorm–Softplus layers, each fed the four features again; a linear output $$o$$. [SRC]
2. On the first minibatch (training mode) the offset $$o_0$$ is the mean of $$o$$ over that batch
   (and it asserts more than 4 droplets have $$s > c$$). [SRC]
3. Logit: $$\ell = (o - o_0) + 10\,\mathrm{sign}(s - c)\sqrt{|s - c|}$$. [SRC]
4. Known cells: $$a = \sigma(50(s - \log(1+C)))$$, $$\ell \leftarrow (1-a)\ell + 10a$$. Known empties:
   $$b = \sigma(50(s - \log(1+t)))$$, $$\ell \leftarrow b\ell - 10(1-b)$$. [SRC]
5. Training (`TraceEnum_ELBO`, $$y$$ enumerated): the ELBO plus, for every droplet, a
   $$N(\ell; \text{cell\_logit}, 2)$$ term ("balance"), a $$N(\ell; -10, 1)$$ term weighted 10 when
   $$n <$$ the crossover count $$e^c$$ and $$N(\ell; +10, 1)$$ weighted 10 when $$n \ge e^c$$; the
   ambient profile is also fitted to the droplets the guide samples as empty (Poisson, weight 0.01);
   the medians of $$\epsilon$$ among cells and among empties are pulled to 1 ($$N(\cdot; 1, 0.01)$$).
   ClippedAdam (clip norm 10) under OneCycleLR, max lr 10 × `--learning-rate`. [SRC]
6. Cell probability (`Posterior._get_latents_map`, after `model.eval()`, so batch norm uses running
   statistics and dropout is off): $$p = \sigma(\ell)$$ for each of the $$T$$ analysed barcodes, in
   batches of 500. [SRC]
7. Call: cell if $$p > 0.5$$ (`consts.CELL_PROB_CUTOFF`; `run.py` writes `_cell_barcodes.csv` and
   `_filtered.h5` with `>`; `collect_output_metrics` counts cells with `>=`). Barcodes outside the
   top $$T$$ are never cells. The call does not depend on `--fpr`, which only sets how many counts
   are removed. [SRC]

Reading of the above (UNVERIFIED, inferred from the terms, not measured): with the −100 prior below
the crossover, the $$\pm 10$$ supervision on either side of it and the signed square-root term, the
call is anchored to "total above the crossover $$e^c$$"; the network mainly decides droplets above
the crossover that look like ambient (calling them empty) and droplets near it. Above the cell
prior $$C$$ the blend forces $$p \to \sigma(10)$$. The docs say the same in words: "cell
probability" means "the probability that a given droplet is not empty", and many low-quality
droplets are expected among the calls (`docs/source/troubleshooting`, "It seems like CellBender
called too many cells"). [DOC]

Earlier releases:

- 0.1.0: $$p$$ from its own MLP, the guide samples $$z$$ only where the logit is ≥ 0, prior
  `cell_logit` = $$\mathrm{logit}((1-0.5)\,n_\text{cells}/T)$$ for every droplet, no supervision;
  $$p = \sigma(\ell)$$ for analysed barcodes (NaN set to 0); cell if $$p > 0.5$$. [SRC]
- 0.2.x: the soft supervision above (−5 / +5) and a prior that, at `cell_logit − 100`, says
  "empty" for every droplet, so cells come only from the likelihood and the +5 supervision of
  droplets at or above the cell prior. [SRC]
- 0.3.0_rc: as 0.2 for the network, with the extra crossover supervision listed in the table, and
  the GMM priors: crossover $$c$$ = mean of the locations of the two heaviest components (after
  dropping components with weight < 0.01 and those within 3 sd of the lowest component). [SRC]

Randomness: `pyro.set_rng_seed` (Python, NumPy and torch) with 0 (0.1, 0.2) or 1234 (0.3+), so
a rerun on the same machine, versions and thread count should repeat [SRC for the seeding;
repeatability UNVERIFIED]. Sources of variation: the train/test split and minibatch order (NumPy,
global state), empties drawn with replacement for each batch, dropout, the Monte Carlo ELBO,
torch's CPU thread count (`--cpu-threads`) and GPU kernels. A user reports about 30 of 11,001
cells differing between two runs of v0.3.0 three months apart in an unchanged environment (issue
#401 comment). [DOC] `Dropout1d` applied to a 2-D (batch × genes) input is treated by recent torch
as one unbatched sample of channels, so it may drop whole droplets rather than genes, and torch's
handling of that case has changed between releases (UNVERIFIED; torch warns about 2-D input to
`dropout1d`).

### Parameters that affect the call (0.3.0+; defaults)

`--expected-cells` (None), `--total-droplets-included` (None), `--force-cell-umi-prior`,
`--force-empty-umi-prior` (None), `--low-count-threshold` 5, `--model full` (`naive`, `simple`,
`ambient`, `swapping`), `--epochs` 150, `--learning-rate` 1e-4, `--constant-learning-rate` off,
`--z-dim` 64, `--z-layers 512`, `--training-fraction` 0.9, `--empty-drop-training-fraction` 0.2,
`--exclude-feature-types` none, `--ignore-features` none, `--projected-ambient-count-threshold` 0.1,
`--num-training-tries` 1 with `--learning-rate-retry-mult` 0.2 (0.5 in 0.2.x), `--cuda`. The cutoff
0.5 is a constant. With `--model naive` no network is trained and there is no cell probability
(UNVERIFIED what the code then does: `latents_map['p']` is None and `run.py` compares it with 0.5).
[SRC]

### Output per barcode

- 0.3.0+: `out.h5` `/droplet_latents/cell_probability` (float, one per analysed barcode, aligned
  with `/droplet_latents/barcode_indices_for_latents`, indices into the raw barcode list), plus
  `cell_size`, `droplet_efficiency`, `background_fraction`, `gene_expression_encoding`;
  `out_cell_barcodes.csv` (the $$p > 0.5$$ barcodes); `out_filtered.h5`; `out_metrics.csv`.
  `cellbender.remove_background.downstream.anndata_from_h5` puts `cell_probability` in `obs`. [SRC]
- 0.1–0.2: `latent_cell_probability` with `barcode_indices_for_latents` in the output HDF5 group,
  `latent_scale`, `latent_RT_efficiency` (0.2), `_cell_barcodes.csv`. [SRC]

### Why docs only, and what can still be compared

Reimplementing it would mean reimplementing a Pyro VAE trained by stochastic variational
inference: hours on a CPU for a full sample (the docs recommend a GPU), dependent on torch / Pyro
versions and hardware, with no claim of bit-for-bit repeatability, so the oracle rule (≥ 99.5 %
agreement within a noise band) has no noise band to check against. The deterministic part is
already the `cellbender-prior` method. CellBender's own output can still be compared: reading
`cell_probability` and `barcode_indices_for_latents` from `out.h5` (or `_cell_barcodes.csv`) gives
a per-barcode call that `compare` and the plots could take as a stored call. The version is not
in the HDF5 (v0.4.0's `/metadata` holds the learning curve, analysed barcodes and features, FPR
and estimator settings); `out.log` starts with "CellBender <version>" (`cli.py`). [SRC]

## 2. EmptyNN

Read in `lkmklsmn/empty_nn` at commit e1e6226 (2024-01-04, the last; `R/emptynn.R`,
`R/neg_create_model.R`, `DESCRIPTION`, and the bundled `EmptyNN_1.0.tar.gz`, whose R files are
identical), and the paper (PMC8369248).

### Identity

- Paper: Yan, Zhao, Simon. "EmptyNN: A neural network based on positive and unlabeled learning to
  remove cell-free droplets and recover lost cells in scRNA-seq data." *Patterns* 2(8), 100311
  (2021), doi:10.1016/j.patter.2021.100311 (PMID 34430929; checked). [DOC]
- Repository: <https://github.com/lkmklsmn/empty_nn>. No licence: no licence file, and
  `DESCRIPTION` reads `License: What license is it under?` (the template placeholder). The
  article's CC BY-NC-ND 4.0 covers the text, not the code. [SRC] [DOC]
- Language: R calling Keras / TensorFlow (`DESCRIPTION` imports `keras (== 2.8.0)`,
  `tensorflow (== 2.8.0)`, `Matrix`; the code also `require`s `caret`, not declared).
  Distribution: GitHub only (`install_github("lkmklsmn/empty_nn")` or the bundled tarball); not on
  CRAN or Bioconductor. No tags or releases; version 1.0 throughout. [SRC]

### Releases that changed the call

None. `R/emptynn.R` and `R/neg_create_model.R` are unchanged since commit 880bd51 (2020-11-01,
"add package file"); commit 660e9fa (2022-03-23) only pinned keras and tensorflow to 2.8.0 (from
`keras (>= 2.3.0)`) and added a copy of the package under `EmptyNN/`. [SRC]

### The method (`emptynn(counts, threshold = 100, k = 5, iteration = 5, batch_size = 16, epoch = 10)`)

Input: a raw count matrix with barcodes as rows (it stops if there are fewer rows than columns).
[SRC]

1. $$n$$ = row sums over all genes; keep barcodes with $$n \ge 10$$. [SRC]
2. P ("negative", empty) = barcodes with $$n \le$$ `threshold`; if there are 10,000 or more, a
   random 10,000 of them (`sample`). U ("unlabelled") = barcodes with $$n >$$ `threshold`. [SRC]
3. Features: the 2,000 genes with the largest column sums over P (`tail(sort(colSums))`, so ties
   follow R's sort order), each barcode's counts divided by its total $$n$$ over all genes. No log.
   [SRC]
4. For each of `iteration` rounds: `caret::createFolds` splits U into `k` folds. For each fold
   $$j$$: train on P labelled 0 plus fold $$j$$ labelled 1 (shuffled with `sample`); predict the other
   $$k-1$$ folds; store the softmax probability of class 1 in that round's column. [SRC]
5. Network (`neg_create_model`): Dense(256, ReLU) → Dropout(0.4) → Dense(128, ReLU) → Dropout(0.3)
   → Dense(2, softmax); binary cross-entropy on one-hot labels, RMSprop (Keras default rate),
   `batch_size` 16, `epoch` 10, `validation_split = 0.2` (Keras holds out the last 20 % of the
   shuffled training rows; the validation loss is not used for anything). A new model per fold.
   [SRC]
6. Score: the mean over rounds of the stored probabilities; cell if the mean $$> 0.5$$. Barcodes
   in P or below 10 counts are never cells. [SRC]

Randomness: `sample` (P subsample, shuffling), `caret::createFolds`, Keras weight initialisation,
dropout and minibatch order. No seed is set anywhere (neither `set.seed` nor
`tensorflow::set_random_seed`), so two runs give different calls unless the caller seeds both R and
TensorFlow. The paper mentions no seed. [SRC] [DOC]

Code and paper disagree [SRC vs DOC]:

- Each round stores one prediction per barcode, not $$k-1$$: `cv_p[test, i] <- y_test_pred[, 2]`
  overwrites the column at every fold, so a barcode keeps the prediction from the last fold in
  which it was a test sample (fold $$k$$'s model, or fold $$k-1$$'s for barcodes in fold $$k$$). The
  paper and README say each barcode is predicted $$(k-1) \times N$$ times and averaged.
- The paper gives defaults $$k = 10$$, $$N = 10$$ and layers of 128, 64 and 2 neurons; the code has
  $$k = 5$$, `iteration` 5, layers 256, 128, 2 with dropout. The README says $$N$$ defaults to 10.
- The paper keeps barcodes with counts "greater than 10"; the code keeps $$\ge 10$$. The paper does
  not mention the 10,000 cap on P.

### Parameters that affect the call

`threshold` 100 (the P / U boundary, in total counts), `k` 5, `iteration` 5, `batch_size` 16,
`epoch` 10; fixed: the 10-count floor, the 10,000 cap on P, 2,000 genes, the architecture,
RMSprop, validation split 0.2, the 0.5 cut. [SRC]

### Output per barcode

A list: `nn.keep`, a logical vector over `rownames(counts)` (every input barcode); `prediction`, a
data frame over U with one column per round (`iteration1`, …) and `mean.crossval`, the score.
[SRC]

### Why docs only, and what can still be compared

No licence (the code cannot be ported into a GPL-3 package), a neural network trained per run on
Keras / TensorFlow 2.8 from R, and no seed: there is no fixed output to validate against, and the
original does not reproduce its own call. A user's own result (`nn.keep` or `mean.crossval`, saved
with the barcodes) could be read as a stored call for `compare`.

## 3. dropkick

Read in `KenLauLab/dropkick` at every tag v1.0.0–v1.2.5 and `master` bd2647a (2023-11-13):
`dropkick/api.py` (`recipe_dropkick`, `auto_thresh_obs`, `filter_thresh_obs`, `dropkick`),
`logistic.py` (`LogitNet`), `util.py` (`_score_lambda_path`, `_fit_and_score`), `scorer.py`,
`setup.py`, `LICENSE`, `LICENSES/`; and the paper (PMC8494217).

### Identity

- Paper: Heiser, Wang, Chen, Hughey, Lau. "Automated quality control and cell identification of
  droplet-based single-cell data using dropkick." *Genome Res* 31(10), 1742–1752 (2021),
  doi:10.1101/gr.271908.120 (PMID 33837131; checked). [DOC]
- Repository: <https://github.com/KenLauLab/dropkick>. Licence: `LICENSE` is MIT ("Copyright (c)
  2020 Cody N Heiser") for dropkick's own code, but the package compiles and ships glmnet's Fortran
  (`dropkick/src/glmnet`, via python-glmnet) under GPL-2 (`LICENSES/GLMNET_LICENSE`,
  `LICENSES/PYTHON_GLMNET_LICENSE`: "under the terms of GNU General Public License Version 2"),
  and adapted scikit-learn scorer code under BSD (`LICENSES/SCIKIT_LEARN_LICENSE`). Whether the
  GPL-2 part is "or later" is not stated in those files (UNVERIFIED). [SRC]
- Language: Python on scanpy, scikit-image and scikit-learn, with Fortran glmnet built by
  `numpy.distutils` (`setup.py`), so building from source needs gfortran and a NumPy that still has
  `numpy.distutils` (removed in NumPy 2.0 / Python 3.12; UNVERIFIED whether the PyPI wheels avoid
  this). Distribution: PyPI `dropkick` (1.0.0–1.2.8; 1.1.3, 1.2.6, 1.2.7 and 1.2.8 have no git
  tag; 1.2.8 uploaded 2023-11-14, the day after the last commit), and a command line
  (`python -m dropkick` / `dropkick run`). [SRC]

### Releases that changed the call

| Release (tag date) | Change |
|---|---|
| 1.0.0 (2020-03-11) | labels from Otsu on two metrics: `arcsinh_n_genes_by_counts` (keep above) and `pct_counts_ambient` (keep below); barcodes with < 10 genes dropped; `LogitNet` defaults `standardize=True`, CV scored by accuracy; `max_iter` 1000, `n_jobs` −1 |
| 1.0.1 (2020-04-03) | `min_genes` 50 replaces the fixed 10 |
| 1.0.2 (2020-04-17) | `n_ambient` (10): ambient genes are the `n_ambient` with the lowest dropout rate (ties at the boundary included) |
| 1.1.0 (2020-04-25) | multi-Otsu on genes, Otsu on ambient percentage; barcodes below the first multi-Otsu threshold dropped from training |
| 1.1.2 (2020-05-13) | that dropping removed (all barcodes ≥ `min_genes` trained on) |
| 1.1.4 (2020-06-05) | default labels from genes only (multi-Otsu on `arcsinh_n_genes_by_counts`, keep above); `max_iter` 2000, `n_jobs` 2 |
| 1.1.5 (2020-06-16) | dropping below the first threshold restored: those barcodes are neither trained on nor scored (score 0) |
| 1.2.0 (2020-07-28) | `LogitNet(standardize=False, scoring="log_loss")` (was standardised features, accuracy) |
| 1.2.1–1.2.5, master | output handling only (overwriting existing columns; label stored as the strings `"True"`/`"False"` from 1.2.3; master initialises the `.obs` columns to NaN before filling) |

All [SRC], from diffs of the four call functions and `logistic.py` / `util.py` between consecutive
tags. 1.1.3 and 1.2.6–1.2.8 exist only on PyPI; the commits after v1.2.5 change no call function,
so 1.2.6–1.2.8 should call as 1.2.5 (UNVERIFIED against the PyPI files).

### The method (1.2.0 and later; `dropkick(adata, ...)`)

Input: an AnnData of raw counts, all barcodes. [SRC]

1. Preprocessing (`recipe_dropkick` on a copy): drop barcodes with fewer than `min_genes` (50)
   genes detected and genes with zero counts; flag mitochondrial genes (`^mt-|^MT-`) and the
   `n_ambient` (10) genes with the lowest dropout rate (fraction of barcodes with zero) as
   "ambient"; `sc.pp.calculate_qc_metrics`; `arcsinh_n_genes_by_counts` $$= \operatorname{arcsinh}(\text{genes detected})$$.
   `sc.pp.normalize_total` (target the median total), `log1p` and
   `sc.pp.highly_variable_genes(n_top_genes=2000, flavor="seurat", n_bins=20)`; the model input
   is $$\operatorname{arcsinh}$$ of the normalised counts, centred and scaled per gene
   (`sc.pp.scale`) over the barcodes kept. [SRC]
2. Training labels: `skimage.filters.threshold_multiotsu` (3 classes, 256 bins) on
   `arcsinh_n_genes_by_counts` gives $$t_1 < t_2$$. Barcodes with value $$\le t_1$$ are dropped;
   label 0 if $$t_1 <$$ value $$\le t_2$$; label 1 if value $$> t_2$$. [SRC]
3. Model: glmnet elastic-net logistic regression (`LogitNet`) on the 2,000 HVG columns, `alpha`
   0.1 (the L1 share), 100 lambdas down to $$10^{-4}$$ of the largest, intercept, no further
   standardisation, `tol` 1e-7, `max_iter` 2000. [SRC]
4. Lambda: 5-fold `StratifiedKFold(shuffle=True, random_state=seed)` with `seed` 18; in each fold
   the HVGs are chosen again on the training part (from the `log1p` layer) and the path refitted;
   each lambda is scored by negative log loss. $$\lambda_\text{best}$$ is the largest lambda whose
   mean score is at least the best mean score minus one standard error (`scipy.stats.sem`), the
   one-standard-error rule. The final coefficients are the full-data fit's at
   $$\lambda_\text{best}$$. (With several `alphas`, the alpha whose model has the best training
   accuracy at its own $$\lambda_\text{best}$$ wins.) [SRC]
5. Score and label: `dropkick_score` = predicted probability of class 1 for the barcodes kept in
   step 2; `dropkick_label` = argmax of the two class probabilities, so cell if score $$> 0.5$$ (at
   exactly 0.5 the argmax takes class 0). Barcodes dropped in steps 1 or 2 get score 0 and label
   False. The paper says "≥ 0.5"; the code's tie goes the other way. [SRC] [DOC]

Randomness: only the CV fold split, seeded (18), and glmnet's coordinate descent is
deterministic; the CV runs in threads (`n_jobs` 2) but each fold is independent. So a run is
deterministic for fixed versions [SRC]. The call depends on the versions of scanpy
(`normalize_total`, `highly_variable_genes`, `scale`), scikit-image (`threshold_multiotsu`; its
results changed when 0.20 moved to 32-bit floats, see the generic Otsu notes) and scikit-learn
(`StratifiedKFold`) (UNVERIFIED how much).

### Parameters that affect the call (1.2.x defaults)

`min_genes` 50, `mito_names` `^mt-|^MT-` (only for QC columns), `n_ambient` 10 (only for QC
columns in the default), `n_hvgs` 2000, `metrics ["arcsinh_n_genes_by_counts"]`,
`thresh_methods ["multiotsu"]` (also `otsu`, `li`, `mean`), `directions ["above"]`, `alphas [0.1]`,
`max_iter` 2000, `seed` 18. `n_jobs` does not change the result. The CLI exposes the same. [SRC]

### Output per barcode

In the caller's `adata.obs`: `dropkick_score` (float, 0 for barcodes not scored),
`dropkick_label` (categorical `"True"` / `"False"` from 1.2.3; bool categorical in 1.2.1–1.2.2),
the metric columns used (`arcsinh_n_genes_by_counts`); in `.var`: `dropkick_coef`; in `.uns`:
`dropkick_thresholds`, `dropkick_args` (chosen alpha and lambda, lambda path, CV scores, seed).
Returns the fitted `LogitNet`. Note that dropkick's label is "good cell" (it is a QC filter as well
as a cell caller), so it removes low-quality cells EmptyDrops would keep. [SRC] [DOC]

### Why docs only, and what can still be compared

It is deterministic and could in principle be reproduced, but its arithmetic is glmnet's Fortran
coordinate descent (GPL-2, possibly incompatible with GPL-3), on top of scanpy preprocessing whose
results move between scanpy versions; an exact port would mean porting glmnet. Its label is a
quality call, not a cell / empty call. dropkick is itself pip-installable and writes into the same
AnnData, so its `dropkick_label` / `dropkick_score` can be compared with other methods directly
as a stored call. The labelling heuristic (multi-Otsu on arcsinh genes detected) is already
expressible with the package's `generic` multi-Otsu, but that is not dropkick's call.

## 4. BD Rhapsody Sequence Analysis Pipeline: putative cells

Closed source. Read from BD's documentation only: the online guide
<https://bd-rhapsody-bioinfo-docs.genomics.bd.com/> ("BD Rhapsody Sequence Analysis Pipeline
4.0"; pages `steps/steps_putative_cells.html`, `setup/input/parameters.html`,
`release_notes.html`, `outputs/top_outputs.html`, read November 2026) and the *BD Single-Cell
Multiomics Bioinformatics Handbook* 23-21713(04), 2023-01, which covers pipeline v1.12
(bdbiosciences.com PDF). Older online pages (2.0–3.0) could not be read (web.archive.org was not
reachable), so the basic rule between 1.12 and 4.0 is known only through the release notes.

### Identity

- No method paper; the documentation is the only description. [DOC]
- Source: not public. Distributed as the Docker image `bdgenomics/rhapsody` (tags per release),
  CWL workflows and example YAML from `bitbucket.org/CRSwDev/cwl`, Seven Bridges, and from 2.1 an
  experimental Docker-free install bundle. The documentation states no software licence; its
  copyright page covers the guide ("© 2023 BD. All rights reserved", "For Research Use Only").
  Whether the image holds readable source was not checked. [DOC]
- Language: not stated (UNVERIFIED). [DOC]

### Releases that changed the call (release notes)

| Release | Date | Change |
|---|---|---|
| 1.6 | 2019-06-10 | `Exact_Cell_Count` added; option to disable refined calling (basic only) |
| 1.9.1 | 2020-10-06 | algorithm changed "to reduce overcalling of putative cells in high cell input experiments" |
| 1.10 | 2022-01-24 | option to call cells on AbSeq read counts ("for troubleshooting only") |
| 1.11 | 2022-08-18 | intronic reads counted in molecules (with a parameter); protein-aggregate flag for AbSeq calls |
| 1.11.1 | 2022-12-15 | "Improved basic and refined putative cell calling algorithms" |
| 2.0 | 2023-06-14 | basic algorithm becomes the default (refined optional); `Expected_Cell_Count` added |
| 2.1 | 2023-11-10 | basic calling changed to "prevent most inappropriately high cell calls derived from noise signatures"; VDJ high-quality cell designation |
| 2.3 | 2025-03-10 | `Exact_Cell_Count` fixed for ATAC-only and joint mRNA–ATAC calling |
| 3.0 | 2025-10-29 | "Basic putative cell calling algorithm updated to fix several edge cases and get more precise cell calls" (about 1 % more cells); VDJ-only cell calling |
| 4.0 | 2026-09-15 | fix: calling failed "when the last point of the curve was the only minimum" |

All [DOC]. Before 2.0 the refined algorithm was the default for mRNA. None of the changes is
described beyond the quoted lines.

### The method

Input: per cell label, the number of reads of all RSEC-adjusted molecules (DBEC-adjusted in the
v1.12 handbook) of the selected bioproduct type (mRNA by default; AbSeq optional); for ATAC, the
number of transposase sites in peaks. Note the unit: reads, not molecules. [DOC]

Basic algorithm, v1.12 (handbook) [DOC]:

1. Sort cell labels by reads, descending; take the cumulative reads curve on a log10 scale; compute
   its second derivative.
2. An inflection point is a minimum of the second derivative. It is valid if it is at least half as
   deep as the global minimum and $$\le -0.3$$.
3. The smoothing window of the second derivative grows until there are two valid inflection
   points; the one giving the larger cell number is taken (to include low-content populations such
   as lymphocytes beside plasma cells). The window sizes and the smoother are not given.

Basic algorithm, 4.0 (online guide) [DOC]:

1. Same curve and second derivative. No smoothing is mentioned.
2. The second-derivative curve is shifted so that it is entirely negative and its mean over the
   whole curve taken. Local minima are checked against depth multiples of that mean: first
   $$\ge 2.0\times$$, then $$1.75\times$$, then $$1.5\times$$; the first level with at least one
   minimum ends the search. If none reaches $$1.5\times$$, the deepest minimum is used (flagged as a
   sign of poor sample quality).
3. One qualifying minimum: use it. Several: take the one with the most cells such that the cell
   count is below a dynamic maximum, $$\max(1000,\ 0.6 \times \#\{\text{labels with} \ge 10 \text{ reads}\})$$,
   and the read count of the cell at the inflection point is above 50 (75 for VDJ cell calling).
4. `Expected_Cell_Count` given: the inflection point closest to it is chosen, bypassing step 3.
5. `Exact_Cell_Count` $$N$$ given: the top $$N$$ labels by reads (by transposase sites in peaks for
   ATAC). For joint mRNA–ATAC calling, labels are ranked in each modality and the $$N$$ with the best
   mean rank are taken.

Refined algorithm, mRNA / AbSeq (optional from 2.0; the default before) [DOC]:

1. Basic set A.
2. False positives: highly variable genes as in Macosko et al. (log10(count + 1); mean and
   dispersion = variance / mean per gene; 20 bins by mean; normalised dispersion = (dispersion −
   bin mean) / bin sd). For each of several cutoffs on normalised dispersion (three in the example;
   values not given), rerun the second-derivative analysis on those genes, giving $$B_i$$; noise
   $$N_i = A - B_i$$; $$C = A - \bigcap_i N_i$$.
3. False negatives: per gene, total reads over all labels against reads in C; a line of best fit;
   genes with residuals at least one standard deviation from the median residual are "enriched in
   noise"; rerun the second-derivative analysis on those genes, giving D.
4. $$E = C \cup D$$, minus labels below a "dynamically determined minimum threshold" of molecules
   (not given). The handbook says the refined step applies only to mRNA-based calling.

Refined algorithm, ATAC (`Cell_Calling_ATAC_Algorithm Refined`) [DOC]:

1. Drop labels whose fraction of transposase sites in peaks is below 0.1.
2. Count filter: bounds at the bottom and top 1 % of log10 transposase-site counts (lower bound at
   least 1, i.e. 10 sites); for each threshold from lower to upper bound in steps of 0.1, fit a
   two-component Gaussian mixture and record the weight difference (non-cell minus cell; the cell
   component is the one with the higher 90 % confidence interval); a threshold is invalid when the
   two 90 % intervals overlap completely; a stability score counts how many following thresholds
   keep decreasing the weight difference; take the valid threshold with the highest stability,
   ties to the larger weight difference, then the larger threshold.
3. A two-component mixture on site counts; the largest count in the non-cell component is the
   threshold for putative cells.
4. A two-component mixture on (site count, fraction in peaks); ignored if the 90 % ellipses
   overlap, otherwise the component with the higher mean count is the redefined cell cluster.
5. Fit one 2-D Gaussian to putative cells and one to non-cells (excluding low-quality labels);
   per-label log-likelihood ratio; threshold it with the basic algorithm (without the half-depth
   criterion, which suggests that criterion still exists in 4.0 [UNVERIFIED]); remove labels below
   it (if step 4 was ignored) or recover labels at or above it (if step 4 applied).

Multiome (`Cell_Calling_Data mRNA_and_ATAC`): the default joint call is the intersection of the
mRNA and ATAC putative cells; alternatively call on one assay and drop cells absent from the other.
AbSeq-based calls flag suspected protein aggregates (high expression across most AbSeq targets) in
`<sample>_Protein_Aggregates_Experimental.csv`. [DOC]

Randomness: none described; the Gaussian-mixture fits of the ATAC refined algorithm may involve a
random initialisation (UNVERIFIED).

### Parameters that affect the call (4.0)

`Cell_Calling_Data` (`mRNA` default; `AbSeq`, `ATAC`, `mRNA_and_ATAC`, `VDJ`),
`Cell_Calling_Bioproduct_Algorithm` (`Basic` default, `Refined`), `Cell_Calling_ATAC_Algorithm`
(`Basic` default, `Refined`), `Expected_Cell_Count` (none), `Exact_Cell_Count` (none, ≥ 1),
`Exclude_Intronic_Reads` (`False`; changes molecule counts, and so perhaps reads per label
[UNVERIFIED]). In 2.0 the refined switch was `Enable_Refined_Cell_Call`; before 2.0 the option
was "Disable Refined Putative Cell Calling". [DOC]

### Output per barcode

Putative-cell tables `<sample>_RSEC_MolsPerCell_MEX.zip` (and `_DBEC_` for targeted or AbSeq);
`<sample>_RSEC_MolsPerCell_Unfiltered_MEX.zip`, "all cell labels with ≥ 10 reads"; Seurat
`.rds`; in v1.x `<sample>_Putative_Cells_Origin.csv` (basic / refined per putative cell). So a
per-barcode call is "in the putative-cell table or not", over the labels with at least 10 reads.
[DOC] The 4.0 output index lists no per-label reads or second-derivative data. [DOC]

### Why docs only, and what can still be compared

The code is closed and the documentation leaves out constants (smoothing windows, derivative
method, dispersion cutoffs, the minimum-molecule threshold), so no reimplementation could be
validated; the rule also works on reads per cell label, which a molecule matrix does not carry.
Comparable: yes. The unfiltered MEX (labels with ≥ 10 reads) gives the raw matrix and the
putative-cell MEX gives BD's call, so other methods can run on the former and be compared with the
latter, keeping in mind that labels with < 10 reads are absent.

## 5. ScaleBio ScaleRna: cell calling

Source-available, not open source. Read in `ScaleBio/ScaleRna` at tags v1.4.1 (2024-02-09),
v1.5.0 (2024-03-19), v1.6.0, v2.0.1 (2025-05-06), v2.1.0 (2025-07-01) and v2.1.2 (2026-08-21, the
last; commit 32fac72): `bin/call_cells.py`, `bin/cell_finder.py` (`bin/callCells.py` and
`bin/cellFinder.py` in 1.5–1.6, `bin/getSampleMetrics.py` before 1.5), `bin/filter_beads.py`, `modules/create_mtx.nf`,
`modules/create_mtx.config`, `nextflow.config`, `docs/cellCalling.md`, `changelog.md`,
`LICENSE.md`. The "ScaleRna workflow handbook" that `docs/cellCalling.md` cites was not read.

### Identity

- No method paper; CellFinder's docstring cites Lun et al. 2019 (EmptyDrops). [SRC]
- Repository: <https://github.com/ScaleBio/ScaleRna>. Licence: `LICENSE.md`, Scale
  Biosciences' end-user licence agreement: a "limited, non-exclusive, non-transferable,
  non-sublicensable license to use the Licensed Software … solely in connection with the use of
  ScaleBio Product(s)"; §3.1 forbids modifying it, creating "any derivative work", or combining
  it with other software, and limits use to "data generated using ScaleBio Products". So the
  source can be read but not ported or reused. [SRC]
- Language: Nextflow (DSL2) with Python (pandas, NumPy, SciPy, DuckDB) and R; counting by STARsolo
  (`GeneFull_Ex50pAS`, multimappers `PropUnique`, so counts can be fractional). Distribution:
  GitHub, containers on AWS ECR (changelog 1.5). [SRC]

### Releases that changed the call

| Release | Change |
|---|---|
| ≤ 1.4.1 | TopCell only (`getSampleMetrics.py`): threshold $$\max(P_{99}/\texttt{minCellRatio},\ \texttt{minReads})$$ on `umis`, no rounding; `minReads` 100; or a fixed number of cells (`useSTARthreshold` takes STARsolo's estimate) |
| 1.5.0 | `--cellFinder` (EmptyDrops-like) added, off by default; `minReads` renamed `minUTC`; `UTC`, `fixedCells`, `expectedCells` reorganised; threshold rounded |
| 2.0.0/2.0.1 | `cellFinder` on by default (QuantumScale RNA becomes the default kit); barnyard samples called per species (CellFinder off for barnyard) |
| 2.1.0 | cell calling per sub-library; `medianFraction` 0.05: barcodes below 5 % of the median count of the threshold-passing cells are never tested; BH now includes zero p-values for "already called" barcodes (see the note below); beads with an ambient-like RT distribution filtered by default (`filterAmbientBeads`) |
| 2.1.1–2.1.2 | no change to calling (metrics, sample sheets) |

[SRC] (code diffs and `changelog.md`). **Correction** to B11 in
[other-methods-and-reviews.md](other-methods-and-reviews.md): the EmptyDrops-like mode is the
default only from 2.0; in 1.5–1.6 the default is TopCell.

### The method (2.1.x)

Per sample (and per sub-library), on `counts`, the unique transcript count per barcode (UTC). [SRC]

1. Ambient barcodes: count $$<$$ `minUTC` (100).
2. Threshold: `UTC` if set (> 0); otherwise TopCell: the preliminary set is the top
   `expectedCells` barcodes if $$0 <$$ `expectedCells` $$<$$ number of barcodes, else all
   non-ambient barcodes; threshold $$= \max\big(\mathrm{round}(P_{\texttt{topCellPercent}} / \texttt{minCellRatio}),\ \texttt{minUTC}\big)$$
   with NumPy's linear `percentile` (99) and ratio 10. [SRC]
3. Without CellFinder: cell if count $$\ge$$ threshold. With `fixedCells`: the top `expectedCells`
   by count (`argsort`, so ties follow NumPy's sort). [SRC]
4. With CellFinder (default, not for barnyard samples): cell if count $$\ge$$ threshold;
   "ambiguous" if `minUTC` $$\le$$ count $$<$$ threshold and count $$\ge$$ `medianFraction` ×
   median count of the threshold cells. Each ambiguous barcode is tested: [SRC]
   - Genes with zero total over the sample are dropped. Ambient profile: simple Good-Turing
     (Gale–Sampson) on the per-gene sums over ambient barcodes, with a pseudo-gene of count 1
     appended and removed; genes unseen among ambient barcodes share the unseen mass in proportion
     to their total over all barcodes. The switch from Turing to smoothed estimates uses
     $$0.95 \times$$ the standard deviation, where Gale and Sampson and edgeR's `goodTuring`
     (`conf = 1.96`, used by DropletUtils) use 1.96. [SRC]
   - Overdispersion $$\alpha$$ (if `--alpha` is 0, the default): 1,000 ambient barcodes with
     count $$\ge 10$$ drawn without replacement; maximise the 0.1 %-trimmed mean of their
     Dirichlet-multinomial log-likelihoods under $$\alpha \times$$ profile with SciPy's bounded
     scalar search on $$[1, 10000]$$, `xatol` 0.001; 3000 if the search fails. [SRC]
   - Counts are rounded to integers; each ambiguous barcode's DM log-likelihood under
     $$\alpha \times$$ profile is computed. [SRC]
   - Monte Carlo: 10,000 probability vectors from Dirichlet($$\alpha \times$$ profile), each drawn
     multinomially to the smallest ambiguous total, then extended one count at a time (gene drawn
     from that vector's probabilities) up to the largest ambiguous total, the log-likelihood
     updated in closed form ($$+\log\frac{n}{n+\alpha-1} + \log\frac{x_g+\alpha_g-1}{x_g}$$ after
     the increment). For a barcode with total $$n$$, $$p = (R + 1)/(10001)$$, $$R$$ = number of
     simulated vectors at total $$n$$ with log-likelihood $$\le$$ the barcode's. [SRC]
   - Benjamini–Hochberg (`scipy.stats.false_discovery_control`) on the p-values with $$k$$ zeros
     prepended; rescued if adjusted $$p \le$$ `cellFinderFdr` (0.001); flagged `cellFinder`. [SRC]
5. Ambient beads (QuantumScale, `filterAmbientBeads` true): per bead, the KL divergence (base 2)
   between its distribution of counts over RT barcodes and the library's, divided by
   $$\log_2(\#\text{RT barcodes})$$, over beads with at least two RT barcodes at $$\ge$$ `minUTC`
   and libraries with at least 8 RT barcodes; passing barcodes on beads with normalised divergence
   $$<$$ `minBeadDivergence` (0.05) are failed and flagged `ambient_bead`. [SRC]
6. MAD outliers (log reads ±5 MAD, passing-read fraction −8 MAD, mito fraction +5 MAD) are flagged
   but not removed unless `filterOutliers`. [SRC]

Note on step 4's BH: `call_cell_barcodes` passes `len(cell_bcs)`, the length of the Boolean array
over every barcode in the sample, as the number of "already called (highly confident) cells", so
$$k$$ is the number of barcodes, not the number above the threshold. With $$k$$ much larger than the
number of ambiguous barcodes the adjustment is close to none. [SRC; the effect on calls is
UNVERIFIED] In 1.5–2.0, BH ran on the ambiguous barcodes alone. [SRC]

Randomness: a module-level `numpy.random.default_rng(42)` (PCG64) drives the ambient sample, the
Dirichlet and multinomial draws and the increments, so the call is deterministic for the same input
and NumPy version; in barnyard runs the state carries over between species. [SRC]

Differences from DropletUtils' `emptyDrops` (from the code, not tested): the upper bound is the
TopCell threshold, not the knee; the lower bound is `minUTC` 100 with "ambient" meaning $$<$$ 100
(DropletUtils: $$\le$$ `lower`); $$\alpha$$ is always estimated, on a 1,000-barcode sample with a
trimmed likelihood (DropletUtils 1.28+ defaults to the plain multinomial); the extra
`medianFraction` floor; the BH set; fractional STARsolo counts rounded; NumPy's generator rather
than R's. [SRC] UNVERIFIED how much these change the calls.

### Parameters that affect the call (2.1.x defaults)

`minUTC` 100, `UTC` 0 (off), `cellFinder` true, `cellFinderFdr` 0.001, `medianFraction` 0.05
(`--medianFraction` defaults to 0 in the script, but the workflow passes the config value), the
script-only `--alpha` 0 (estimate), `topCellPercent` 99, `minCellRatio` 10, `expectedCells` 0 (per
sample in `samples.csv`), `fixedCells` false, `filterAmbientBeads` true with `minBeadDivergence`
0.05, `filterOutliers` false, `starFeature` `GeneFull_Ex50pAS`, `starMulti` `PropUnique`
(`Unique` for barnyard), `roundCounts` false. [SRC]

### Output per barcode

`<sample>_allBarcodes.parquet` (every barcode: `counts`, `pass`, `flags` such as `cellFinder`,
`ambient_bead`, `low_reads`, `high_mito`; barnyard `species`), `<sample>_allCells.csv` (passing
barcodes), `<sample>_cellcalling_stats.csv` (numbers of ambiguous and ambient barcodes, the
ambiguous range, the estimated $$\alpha$$, the smallest p-value), the filtered STARsolo matrix, and
optionally an AnnData (`--annData`). [SRC]

### Why docs only, and what can still be compared

The licence forbids derivative works, so nothing of it may be ported, and the documentation is too
thin to write an independent version that could be checked. The workflow's own outputs are good for
comparison: STARsolo's raw matrix plus `allBarcodes.parquet` give every barcode's `pass` and flags,
which can be loaded as a stored call beside other methods run on the same raw matrix. TopCell
alone ($$\max(\mathrm{round}(P_{99}/10), 100)$$) is a generic rule the `generic` module could
express, but it would not be ScaleBio's code.

## 6. Singleron CeleScope (`celescope rna`): cell calling

Read in `singleron-RD/CeleScope` at tags v1.0.1 (2020-08-07), v1.3.0, v1.5.2, v1.8.0, v1.9.0
(2022-04-01), v1.17.0 (2023-08-11), v2.0.0 (2023-10-26), v2.0.9, v2.15.6 (2026-11-11) and `master`
2679b7e (2026-11-20): `tools/count.py` (`celescope/tools/count.py` from 1.x), `celescope/tools/
emptydrop_cr/cell_calling_3.py` (+ `sgt.py`, `stats.py`), `celescope/tools/starsolo.py`,
`celescope/tools/cells.py`, `celescope/rna/__init__.py`, `celescope/tools/multi.py`,
`doc/CHANGELOG.md`, `LICENSE`, `conda_pkgs.txt`.

### Identity

- No method paper for the pipeline. [DOC]
- Repository: <https://github.com/singleron-RD/CeleScope> (mirror on gitee). Licence: MIT, file
  `LICENSE` ("Copyright (c) 2020 SingleronBio"). The 1.x EmptyDrops_CR module is marked as a
  modified copy of Cell Ranger's `sgt.py` and `stats.py` (Cell Ranger 3.0.2 was MIT; whether the
  copied files keep 10x's notice was not checked). [SRC]
- Language: Python, with STAR (bioconda `star=2.7.11a` from 2.0.0; 2.6.1 before), R for some
  steps. Distribution: PyPI / conda environment from `conda_pkgs.txt`, Docker (`Dockerfile`). [SRC]

### Releases that changed the call

| Releases | Default rna cell call |
|---|---|
| 1.0.1 – 1.8.x (2020-08 – 2022-03) | `auto` (CeleScope's own rule, below), `--expected_cell_num` 3000 (`--cells 3000` in 1.0.1); options `cellranger3` (1.x port of Cell Ranger 3.0.2's EmptyDrops) and, until 1.5.x, `inflection` (an R script, `rescue.R`, removed in 1.6.0) |
| 1.9.0 – 1.17.x (2022-04 – 2023-08) | `EmptyDrops_CR`, the same port renamed (`cellranger3` → `EmptyDrops_CR`), expected 3000; `auto` still available |
| 2.0.0 – 2.15.6, master (2023-10 – ) | STARsolo 2.7.11a's own filter: `--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000` on the `GeneFull_Ex50pAS` matrix (the reported feature); `tools/count.py` (with `auto`, `EmptyDrops_CR`, `--expected_cell_num` now 8000) is no longer in any assay's step list |

[SRC] (argument defaults at each tag; the rna `STEPS`; changelog 2.0.0 "Use Starsolo to quantify
the gene expression"). **Correction** to B11 in
[other-methods-and-reviews.md](other-methods-and-reviews.md): "auto = Cell Ranger v2 ordmag
(expected 8000)" describes dead code in 2.x; `auto` was the default only in 1.0–1.8, with expected
3000, and it is not Cell Ranger's ordmag (no percentile, no bootstrap, see below); the default
since 2.0 is STARsolo's EmptyDrops_CR, not CeleScope's own.

### The method

**2.0 and later (default).** The `starsolo` step runs STAR 2.7.11a with
`--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000`, i.e.
nExpectedCells 3000, maxPercentile 0.99, maxMinRatio 10, indMin 45000, indMax 90000, umiMin 500,
umiMinFracMedian 0.01, candMaxN 20000, **FDR 0.001** (STAR's own default is 0.01), simN 10000;
`--soloFeatures GeneFull_Ex50pAS Gene`, `--report_soloFeature GeneFull_Ex50pAS`. The filtered
matrix of the reported feature is the call; the algorithm is STARsolo's (see
[cellranger-starsolo.md](cellranger-starsolo.md), `starsolo_emptydrops_cr`). [SRC] The 2.0.0
changelog also adds a filter of UMIs with fewer supporting reads than `read_threshold`, which
changes counts before calling (not read further). [DOC]

The optional `celescope rna cells` step (not run by `multi_rna`) re-filters the raw matrix: a
user barcode list (`--barcode`), or `--force_cells N` (the top N rows of the counts file, by UMI),
or `--soloCellFilter <string>` (runs `STAR --runMode soloCellFiltering`), then optionally
`--max_mito` (fraction of mitochondrial UMIs) and `--min_gene`; the original call is kept as
`default_filtered`. [SRC]

**1.9–1.17 `EmptyDrops_CR` (and 1.0–1.8 `cellranger3`)** (`cell_calling_3.py`, a port of Cell
Ranger 3.0.2's): on CeleScope's own UMI matrix (featureCounts-based in 1.x); ambient profile by
simple Good-Turing on barcodes at descending ranks [45000, 90000) with nonzero counts, zero-count
features sharing $$p_0$$ equally; the initial call is Cell Ranger's bootstrapped ordmag
(`filter_cellular_barcodes_ordmag`, `recovered_cells` = `--expected_cell_num` 3000, baseline index
$$\mathrm{round}(0.01 \times 3000)$$, 100 bootstrap resamples); candidates are the next up to
20,000 barcodes by count with count $$\ge \max(500,\ \lceil 0.01 \times$$ median of the initial
cells$$\rceil)$$; multinomial log-likelihoods, 10,000 simulations, Monte Carlo p-values,
Benjamini–Hochberg, cell if adjusted $$p \le 0.01$$. Python `random` and NumPy's global generator
are seeded 0 at import. [SRC] UNVERIFIED how the "modified" `sgt.py` and `stats.py` differ from
Cell Ranger 3.0.2's.

**1.0–1.8 `auto`**: sort barcodes by UMI count, descending; $$i = \lfloor 0.01 \times E \rfloor$$
($$E$$ = expected cells, 3000); threshold $$= \max(1, \lfloor 0.1 \times u_{(i)} \rfloor)$$ with
$$u_{(i)}$$ the count of the $$i$$-th barcode (1-based); cells are barcodes with count $$\ge$$
threshold. So with defaults: one tenth of the 30th-largest count, truncated. Deterministic; ties
in the sort do not matter because the threshold is a count. [SRC] (Cell Ranger 2's ordmag takes
the 99th percentile of the top $$E$$ counts instead; CeleScope's rule is an order statistic.)

`--force_cell_num N` (1.x): from 1.7.2, the count threshold whose cumulative number of barcodes
(at or above it) is closest to N, looking only at the two counts around the first that reaches N;
before 1.7.2, among the counts whose cumulative numbers span N ± 10 %, the count with the most
barcodes (the largest step). Cells are barcodes at or above that count, so ties make the number
differ from N. [SRC]

### Parameters that affect the call

2.x: `--soloCellFilter` (default above), `--soloFeatures`, `--report_soloFeature`
(`GeneFull_Ex50pAS`; `Gene` to drop intronic reads), and the `cells` step's `--force_cells`,
`--barcode`, `--max_mito`, `--min_gene`. 1.x: `--cell_calling_method` (`auto` ≤ 1.8,
`EmptyDrops_CR` ≥ 1.9), `--expected_cell_num` 3000, `--force_cell_num` (or column 4 of the
`multi` mapfile). [SRC]

### Output per barcode

2.x: `outs/raw` and `outs/filtered` (MEX, gzip), the counts file with a `mark` column (`CB` for
called barcodes, `UB` otherwise; written by the `cells` step), h5ad, the HTML report. 1.x: the
marked counts file (`*_counts.txt`, `mark` CB/UB), `*_matrix_10X` raw and filtered. [SRC]

### Why docs only, and what can still be compared

The 2.x default is not a method of its own: it is `starsolo-emptydrops-cr` at STAR 2.7.11a with
`fdr=0.001`, which sc-droplets already implements, run on the `GeneFull_Ex50pAS` raw matrix
CeleScope writes (UNVERIFIED end to end against a CeleScope run; the matrix's UMI filtering by
`read_threshold` happens before). So the page can point to that spec and say how to reproduce it
rather than add an id. The 1.x rules (the `auto` order-statistic rule and the Cell Ranger 3.0.2
port) are historical; `auto` is simple enough to implement if wanted. CeleScope's raw and filtered
matrices compare directly with other methods.

## 7. Ambimux

Read in `marcalva/ambimux` at tag v0.5.0 (commit 11f55af, 2025-09-19, the last; `mod.c`, `mod.h`,
`math_util.c`, `main.c`, `README.md`, `LICENSE.txt`) and the preprint's abstract.

### Identity

- Preprint: Alvarez, Li, Lee, …, Kaikkonen, Zaitlen, Pajukanta. "Integrated ambient modeling and
  genetic demultiplexing of single-cell RNA+ATAC multiome experiments with Ambimux." bioRxiv
  (2025), doi:10.1101/2025.08.21.671671 (posted 2025-08-26, v1; checked against bioRxiv; not
  peer-reviewed). [DOC]
- Repository: <https://github.com/marcalva/ambimux>. Licence: MIT, `LICENSE.txt` ("Copyright (c)
  2021 Marcus Alvarez"); it builds htslib 1.17 (MIT/BSD) locally. [SRC]
- Language: C (C99, pthreads). Distribution: source only (`make hts && make`), Linux; no package.
  A simulator, `marcalva/ambisim`, is separate. [SRC]

### Releases that changed the call

One git tag, v0.5.0, at 11f55af (2025-09-19). The version string in `main.c` went 0.2.0
(2023-11-16), 0.4.0 (2024-01), 0.5.0 (commit 832de39, 2024-09-17); later commits kept "v0.5.0"
while changing the model: ambient-sample-proportion estimation refactored (2024-10-28, "refactor
ambient"; 2025-08-20, "Add robust pi_as estimate", "scale vareps by num mol"), 376 changed lines in
`mod.c`. So "ambimux v0.5.0" names two different models; only the tag is reproducible. [SRC]

### The method (tag v0.5.0)

Inputs: coordinate-sorted, indexed RNA and/or ATAC BAMs with barcode (`CB`) and UMI (`UB`) tags;
a VCF with `GT` genotypes of the pooled donors (bi-allelic SNPs only); a GTF and a peak BED
(reads outside genes / peaks are dropped from the model by default, `--mdl-reads intra`); an
optional exclusion BED. It needs several genotyped donors in one channel: the call rests on
alleles, not counts. [SRC] [DOC]

1. Counting: RNA molecules deduplicated by barcode and UMI, ATAC fragments deduplicated; per
   barcode, the molecules that overlap a variant ("informative"), with each base's Phred quality
   (≥ `--min-phredq` 30), MAPQ ≥ 30. [SRC]
2. Fixed empties: barcodes with fewer than `--out-min` (100) RNA UMIs **and** fewer than 100 ATAC
   fragments are fixed as empty (their ambient fraction set to 1) and are not tested; they inform
   the ambient pool. [SRC] [DOC]
3. Model: for droplet $$d$$, $$H_d \in \{0, 1, 2\}$$ (empty, singlet, doublet) with prior
   $$\lambda$$, donors $$S_d$$ with proportions $$\pi$$ (a doublet's pair has probability
   $$\pi_{s_1}\pi_{s_2}/\sum_{i \ne j}\pi_i\pi_j$$). Each informative molecule comes from the
   ambient pool with probability $$\alpha_d$$ (per modality, per assignment) or from the nucleus
   (split equally between the two nuclei of a doublet); an empty droplet's molecules are all
   ambient. A base's probability given its source is $$(1-\varepsilon)\,p(\text{allele}) +
   \varepsilon/2$$, $$\varepsilon$$ from its base quality (0.01 if missing); $$p(\text{allele})$$ is
   the donor's genotype ($$\{0, 0.5, 1\}$$ from GT) or, for the ambient source, the
   $$\pi_\text{amb}$$-weighted mean of the donors' genotypes. Molecules without variants contribute
   nothing, so the likelihood is allelic only. [SRC]
4. Fitting: $$\pi_\text{amb}$$ first, from the fixed-empty droplets (iterated to convergence, at
   most 10,000 iterations); then EM: per droplet and assignment, $$\alpha$$ by a Newton-type
   maximisation of the marginal likelihood with a Beta-like prior centred on the mean singlet
   $$\alpha$$ (clipped to $$[10^{-4}, 1 - 10^{-4}]$$) and weight `--amb-prior-w` $$10^{-8}$$;
   inner loop updates $$\lambda$$ and $$\pi$$ (pseudocount $$10^{-8}$$) and the ambient genotype
   until the relative parameter change is below `--eps` $$10^{-6}$$ (at most `--max-iter` 100);
   outer loop until the mean singlet $$\alpha$$ moves by less than $$10^{-6}$$ (at most 100).
   Initial $$\alpha$$ 0.5, $$\lambda$$ = (fixed-empty share, 0.9 × and 0.1 × the tested share). [SRC]
5. Output per droplet: $$\log P(H, S, X)$$ for every assignment; posteriors PP0, PP1, PP2 of
   empty / singlet / doublet; `best_type` = the single assignment $$(H, S)$$ with the highest
   joint probability. [SRC]
6. Call: the program does not threshold. The README recommends calling singlets with PP1 above
   0.9, 0.95 or 0.99; a cell-versus-empty call is "not Empty" (or $$1 - \text{PP0}$$ above a
   threshold). [DOC] [SRC]

Reading of the above (UNVERIFIED): because $$\alpha = 1$$ reproduces the empty likelihood, a
droplet whose informative molecules match the ambient genotype mixture is decided by $$\lambda$$;
a droplet is called non-empty when its alleles lean towards one (or two) donors. Droplets with few
informative molecules carry little evidence either way, whatever their total count.

Randomness: none; droplets are processed in an order from a fixed hash-based permutation
(`permute` in `math_util.c`, no seed). With several threads, sums are accumulated under a lock in
thread order, so the last bits may differ between runs (UNVERIFIED). [SRC]

### Parameters that affect the call

`--out-min` 100, `--mdl-reads intra` (`inter`, `all`), `--max-nh` 0 (off), `--min-phredq` 30,
`--rna-mapq` 30, `--atac-mapq` 30, `--eps` 1e-6, `--amb-eps` 1e-6, `--max-iter` 100,
`--amb-prior-w` 1e-8, `--tx-basic`, `--exclude`, the VCF and the donor list (`--samples`), and the
user's threshold on PP0 / PP1. [SRC]

### Output per barcode

`summary.txt` for tested droplets: counts (`n_rna_molecules`, `n_atac_molecules`, informative
counts, `FRIG`, `FRIP`, mitochondrial percentages), `best_type`, `best_sample`,
`best_rna_ambient`, `best_atac_ambient` and their information, `PP0`–`PP2`, `LLK0`–`LLK2`;
`sample_llk.txt.gz`, `alpha_rna.txt.gz`, `alpha_atac.txt.gz`, `lambda.txt.gz`, `pi.txt.gz`; gene,
peak and allele count matrices (MEX) over all barcodes. Fixed empties are absent from
`summary.txt`. [DOC] [SRC]

### Why docs only, and what can still be compared

It needs read-level BAMs and donor genotypes, which a count matrix does not have, and it answers a
different question (which donor, and how much ambient, per droplet) for pooled multi-donor
multiome channels. Its EM is deterministic in principle and could be ported (MIT), but it would
sit outside the package's matrix-in, call-out model. Comparison is possible where a pool has
genotypes: `summary.txt` gives, per barcode, PP0 and `best_type`, which can be joined to the raw
matrix as a stored call (fixed empties count as empty).

## 8. AmbientMapper and scifi-demux

Read in `gomezcan/ambientmapper` at tag v0.1.0 (commit 7e77faf, 2026-09-20; the version the
preprint used, per its README): `src/ambientmapper/genotyping.py`, `assign_streaming.py`,
`pipeline.py`, `cli.py`, `README.md`, `CHANGELOG.md`, `LICENSE`; and `gomezcan/scifi-demux` at
v0.1.3 (2026-09-20, the version the preprint cites) and v0.1.4 (2026-11-23): `src/scifi_demux/`,
`README.md`, `LICENSE`. The preprint's full text was not read, only its abstract.

### Identity

- Preprint: Gomez-Cano, Jiang, Welch, Marand. "Empirical Estimation of Ambient Contamination in
  Combinatorial Single-Cell Methods Using Multi-Reference Mapping." bioRxiv (2026),
  doi:10.64898/2026.09.11.750809 (posted 2026-09-16, v1; checked against bioRxiv; not
  peer-reviewed). [DOC]
- Repositories: <https://github.com/gomezcan/ambientmapper>,
  <https://github.com/gomezcan/scifi-demux>, and the analysis scripts in
  `gomezcan/ambientmapper_manuscript`. Licence: MIT for both (`LICENSE`, "Copyright (c) 2025
  gomezcan"). [SRC]
- Language: Python 3.9–3.12 (pandas, NumPy, DuckDB, pyarrow, pysam, typer). Distribution: GitHub
  (`pip install git+…`); not on PyPI or bioconda as far as the READMEs say. scifi-demux wraps
  UMI-tools, cutadapt, bwa, samtools, Picard and GNU parallel. [SRC]

### What each tool does

- **scifi-demux** is preprocessing only for scifi-ATAC: UMI extraction, adapter trimming,
  demultiplexing by the combinatorial and 10x barcodes, mapping (bwa), duplicate marking, a MAPQ /
  proper-pair / multi-mapping (XA, NM) filter and barcode tagging. It has no cell-calling step: no
  knee, threshold or model on counts per barcode (searched the whole `src/`). [SRC] So it is not a
  cell caller and needs no page beyond a mention.
- **AmbientMapper** calls each barcode empty, singlet, doublet or ambiguous from reads mapped
  competitively to several reference genomes (species or genotypes pooled in one library). [SRC]

### Releases that changed the call

Only v0.1.0 is tagged. The unreleased 0.2 series (CHANGELOG) moves the intermediate files from TSV
to Parquet with "no behavior change" claimed for the calls. [SRC] [DOC]

### The method (AmbientMapper v0.1.0; `ambientmapper run`: extract → filter → chunks → assign → genotyping)

Inputs: one coordinate-sorted BAM per reference genome, same reads, with barcode (`CB`/`BC`),
`AS`, `NM`, MAPQ and optionally `XA`. [DOC]

1. `extract`, `filter`: per genome, a per-read table; keep barcodes with at least
   `min_barcode_freq` (10) reads in that genome's table. [SRC]
2. `assign`: per read, the best and second-best genome by (lower NM, higher AS, higher MAPQ);
   empirical distributions of $$\Delta$$AS and $$\Delta$$MAPQ (best minus second) by decile of the
   winner's score, learnt on 50,000 barcodes sampled with `random.Random(42)`; per read, p-values
   of its $$\Delta$$s; the read is a "winner" for its best genome or "ambiguous" at `alpha` 0.05
   with an NM rule and a dominance filter across genomes. [SRC] [DOC]
3. `genotyping` pass 1: per read and genome, a score $$w_\text{AS}\mathrm{AS} +
   w_\text{MAPQ}\mathrm{MAPQ} - w_\text{NM}\mathrm{NM}$$ (weights 1; the code at v0.1.0 uses the
   raw values, the README says robust z-scores); per read a softmax over genomes with
   $$\beta = 1$$ and an extra ambient constant $$10^{-3}$$ in the denominator gives $$L(r, g)$$.
   With `winner_only` (on by default in `MergeConfig`, off by default in the `run` pipeline unless
   `--genotyping-winner-only`), $$L(r, g) = 1$$ for the read's best genome only. Expected counts
   $$C(b, g) = \sum_r L(r, g)$$. [SRC]
4. Ambient profile $$\eta(g)$$: seed = the 2 % of barcodes with the smallest effective count
   (`eta_seed_quantile`); $$\eta$$ = their summed $$L$$ by genome, normalised; then twice
   (`eta_iters` 2): candidates with top-genome share ≤ 0.6 and top1/top2 ratio ≤ 1.5, keep those
   whose Jensen–Shannon divergence from $$\eta$$ is at or below its 95th percentile, recompute
   $$\eta$$. (The README's "barcodes with < 200 reads" does not match this code.) [SRC]
5. Per barcode (pass 2), with $$E_r = \sum_g \eta(g) L(r, g)$$ and the top 3 genomes by
   $$C(b, g)$$ as candidates:
   - empty: $$\ell_0 = \sum_r \log E_r$$, 0 parameters;
   - singlet $$g_1$$: $$\ell_1 = \max_\alpha \sum_r \log[(1-\alpha) L(r, g_1) + \alpha E_r]$$, $$\alpha$$ on
     a grid 0, 0.02, …, 0.5, 1 parameter;
   - doublet $$g_1, g_2$$: $$\max_{\alpha, \rho} \sum_r \log[(1-\alpha)(\rho L_{g_1} + (1-\rho) L_{g_2}) +
     \alpha E_r]$$, $$\rho$$ on 0.1, 0.15, …, 0.9, 2 parameters;
   - $$\mathrm{BIC} = -2\ell + k \log n$$, $$n$$ the number of reads. [SRC]
6. Call: fewer than `min_reads` (5) reads → `low_reads`; **empty** if
   $$\mathrm{BIC}_0 + 10 < \min(\mathrm{BIC}_1, \mathrm{BIC}_2)$$ (`empty_bic_margin`) **and** the
   top genome holds at most 60 % of the mass **and** top1/top2 ≤ 1.5 (optional read-count and JSD
   gates off by default); otherwise singlet or doublet by a $$\Delta$$BIC of 6 (`bic_margin`),
   else ambiguous, with sublabels (`single_clean`, `dirty_singlet`, `weak_doublet`,
   `indistinguishable`). An optional cross-mapping correction (`xmap`) replaces $$L(r, g_1)$$ by a
   learnt mixture when a genome's cross-mapping share exceeds 0.05. [SRC]

Randomness: only the seeded barcode subsamples (42); grids and thresholds are deterministic. [SRC]

### Parameters that affect the call (v0.1.0 defaults)

`min_barcode_freq` 10, assign `alpha` 0.05 (and the 50,000-barcode subsamples), `beta` 1,
`w_as`/`w_mapq`/`w_nm` 1, `ambient_const` 1e-3, `winner_only`, `min_reads` 5, `empty_bic_margin`
10, `bic_margin` 6, `empty_top1_max` 0.6, `empty_ratio12_max` 1.5, `empty_reads_max` and
`empty_jsd_max` (off), `eta_seed_quantile` 0.02, `eta_iters` 2, `empty_tau_quantile` 0.95,
`alpha_grid` 0.02, `max_alpha` 0.5, `rho_grid` 0.05, `topk_genomes` 3, the `xmap_*` options. [SRC]

### Output per barcode

`final/<sample>_cells_calls.tsv.gz`: `call` (`single_clean`, `dirty_singlet`, `doublet`,
`weak_doublet`, `indistinguishable`, `ambiguous`, `empty`, `low_reads`), `genome_1`, `genome_2`,
`alpha`, `rho`, BICs, `n_reads`, `n_effective`, `p_top1`–`p_top3`, `entropy`, `status_flag`;
`final/<sample>_expected_counts_by_genome.csv` ($$C(b, g)$$). `decontam` then writes a read drop
list and a barcode policy. Barcodes removed by `min_barcode_freq` do not appear. [DOC] [SRC]

### Why docs only, and what can still be compared

AmbientMapper works on reads mapped to several references (multi-species or multi-genotype pools,
designed for scifi-ATAC), not on a count matrix, and the call needs those per-read genome
likelihoods. It is deterministic and MIT, so it could be reproduced, but only from BAMs. Its
`cells_calls.tsv.gz` can be joined to a raw barcode × peak (or gene) matrix of the same library as
a stored call ("not empty and not low_reads"); scifi-demux has nothing to compare.

## Corrections to the existing notes

- [other-methods-and-reviews.md](other-methods-and-reviews.md) B11, ScaleBio: the EmptyDrops-like
  CellFinder is the default only from 2.0; 1.5–1.6 default to TopCell; ≤ 1.4 has TopCell only. The
  TopCell percentile is taken over the barcodes at or above `minUTC` (or the top `expectedCells`)
  and the threshold is rounded and floored at `minUTC`. (§5)
- B11, CeleScope: "auto = Cell Ranger v2 ordmag (expected 8000)" describes code no assay runs in
  2.x. 1.0–1.8 default `auto` with expected 3000 (an order statistic, not ordmag); 1.9–1.17
  default a port of Cell Ranger 3.0.2's EmptyDrops (FDR 0.01); 2.0+ default STARsolo 2.7.11a's
  `EmptyDrops_CR` with FDR 0.001. (§6)
- B11, BD Rhapsody: the 2.0× / 1.75× / 1.5× rule is the 4.0 documentation; the v1.12 handbook
  describes valid minima as at least half the global minimum's depth and ≤ −0.3, with a growing
  smoothing window. (§4)
- C2, EmptyNN: the code stores one prediction per barcode per round (the last fold that tested it),
  not $$k-1$$; the paper's defaults ($$k = 10$$, $$N = 10$$, layers 128/64) differ from the code's
  (5, 5, 256/128 with dropout); barcodes below 10 counts are dropped first. (§2)
- C5, dropkick: the repository is MIT but ships glmnet's Fortran under GPL-2; multi-Otsu on
  arcsinh genes alone is the default only from 1.1.4 (Otsu on genes and ambient percentage
  before). (§3)
- C12, scifi-demux: preprocessing only, no cell calling; AmbientMapper is the caller. (§8)
- [dropletutils-cellbender-cellsweep.md](dropletutils-cellbender-cellsweep.md) §1c: at v0.3.0_rc
  `LOW_UMI_CUTOFF` is 15, and "the first 2 component locs" are the two heaviest components (sorted
  by weight, after dropping weights < 0.01 and components within 3 sd of the lowest). (§1)

## Summary

| Tool | Category | Licence | Open source? | Deterministic? | Inputs | Why docs only |
|---|---|---|---|---|---|---|
| CellBender `remove-background` model (0.1–0.4) | deep generative model (VAE, SVI); cell = $$\sigma(\text{logit}) > 0.5$$ | BSD-3-Clause | yes | seeded (0 / 1234) but not reproducible across hardware or versions | raw matrix, all features | trained Pyro VAE, GPU in practice; no fixed output to validate; priors already implemented as `cellbender-prior` |
| EmptyNN (1.0) | classifier (PU-learning neural network) | none (`License: What license is it under?`) | source visible, no licence | no (no seed) | raw matrix | no licence; Keras/TensorFlow model trained per run; unseeded |
| dropkick (1.0–1.2.8) | classifier (elastic-net logistic regression on self-made labels) | MIT + GPL-2 (glmnet) | yes | yes for fixed versions (seed 18) | raw AnnData | glmnet Fortran (GPL-2) and scanpy-version-dependent preprocessing; a quality call, not cell/empty |
| BD Rhapsody pipeline (1.x–4.0) | knee (second derivative of log cumulative reads); refined HVG and GMM steps | proprietary, none stated | no | not stated | reads per cell label (BAM-level); ATAC: sites in peaks | closed code; constants not documented; works on reads |
| ScaleBio ScaleRna (≤ 1.4–2.1.2) | threshold (TopCell) + EmptyDrops-like test (CellFinder) | ScaleBio EULA (no derivative works) | source-available, not open | yes (NumPy seed 42) | STARsolo raw matrix (fractional counts), RT/bead barcodes | licence forbids porting; outputs comparable |
| Singleron CeleScope (1.0–2.15.6) | 1.x: order-statistic threshold / Cell Ranger 3.0.2 EmptyDrops port; 2.x: STARsolo EmptyDrops_CR | MIT | yes | yes (seeded) | FASTQ → own (1.x) or STARsolo (2.x) matrix | 2.x is `starsolo-emptydrops-cr@2.7.11a` with `fdr=0.001`, already implemented; 1.x historical |
| Ambimux (v0.5.0) | genotype mixture model (EM), empty/singlet/doublet posteriors | MIT | yes | yes in principle (thread order may change last bits) | RNA/ATAC BAMs + donor VCF | needs reads and genotypes of pooled donors; answers a demultiplexing question |
| AmbientMapper (v0.1.0) | multi-reference read mixture, BIC model selection | MIT | yes | yes (seed 42) | one BAM per reference genome | needs reads mapped to several references; scifi-ATAC design |
| scifi-demux (v0.1.3–0.1.4) | preprocessing only | MIT | yes | n/a | scifi-ATAC FASTQ | not a cell caller |
