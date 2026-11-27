---
title: dropkick
parent: Documented only
grand_parent: Methods
nav_order: 2
method_ids: [dropkick]
status: docs-only
---

# dropkick: a logistic regression trained on its own labels
{: .no_toc }

dropkick labels the barcodes with the most genes detected as cells and those on the ambient
plateau as empty, by multi-Otsu thresholds on $$\operatorname{arcsinh}$$ of genes detected. It then
trains an elastic-net logistic regression on 2,000 highly variable genes and scores every barcode
it kept. Its label means "good cell": it is a quality filter as well as a cell caller, and it drops
low-quality cells that a test against the ambient profile would keep. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `dropkick` (documented only) |
|---|---|
| Tool | dropkick 1.0.0 – 1.2.8 (PyPI; git tags v1.0.0 – v1.2.5), `dropkick.dropkick` and `dropkick run` |
| Category | classifier |
| Counts | genes detected per barcode for the labels; normalised counts of 2,000 highly variable genes for the model |
| Needs expected cells | no |
| Deterministic | yes for fixed releases of dropkick, scanpy, scikit-image and scikit-learn: the folds are seeded (18) and glmnet's coordinate descent draws nothing |
| Original licence | MIT (dropkick's code); GPL-2 (the glmnet Fortran it compiles and ships); BSD (adapted scikit-learn scorers) |
| Extra inputs | none; the raw AnnData |

## How it works

1.2.0 and later; the earlier rules are under [Versions](#versions).

1. **Preprocessing** (`recipe_dropkick`, on a copy). Drop barcodes with fewer than `min_genes`
   (50) genes detected, and genes with no count. Flag the mitochondrial genes, and flag as ambient
   the `n_ambient` (10) genes with the lowest dropout rate (the share of barcodes where they are
   zero); by default both only fill QC columns. $$g = \operatorname{arcsinh}(\text{genes detected})$$
   (`arcsinh_n_genes_by_counts`). Normalise each barcode to the median total
   (`sc.pp.normalize_total`), take `log1p`, and pick 2,000 highly variable genes
   (`flavor="seurat"`, 20 bins). The model's input is $$\operatorname{arcsinh}$$ of the normalised
   counts, centred and scaled per gene over the barcodes kept (`sc.pp.scale`).
2. **Labels.** scikit-image's `threshold_multiotsu` (3 classes, 256 bins) on $$g$$ gives
   $$t_1 < t_2$$. Barcodes with $$g \le t_1$$ are dropped; $$t_1 < g \le t_2$$ are labelled 0
   (empty); $$g > t_2$$ are labelled 1 (cell). This is the [multi-Otsu](../generic/multiotsu.md)
   rule.
3. **Model.** glmnet's elastic-net logistic regression (`LogitNet`) on the 2,000 genes: `alpha`
   0.1 (the L1 share), a path of 100 lambdas down to $$10^{-4}$$ of the largest, an intercept, no
   further standardisation, tolerance $$10^{-7}$$, at most 2,000 iterations.
4. **Lambda.** Five folds, `StratifiedKFold(shuffle=True, random_state=18)`. In each fold the
   highly variable genes are chosen again on the training part (from the `log1p` layer) and the
   path is refitted; each lambda is scored by its negative log loss. $$\lambda_\text{best}$$ is
   the largest lambda whose mean score is at least the best mean score less one standard error
   (`scipy.stats.sem`). The coefficients are the full-data fit's at $$\lambda_\text{best}$$. With
   several `alphas`, the one whose model has the best training accuracy at its own
   $$\lambda_\text{best}$$ wins.
5. **Score and label.** `dropkick_score` is the predicted probability of class 1 for the barcodes
   kept in step 2. `dropkick_label` is the argmax of the two class probabilities, so a cell when
   the score is above 0.5; at exactly 0.5 the argmax takes class 0, where the paper says
   $$\ge 0.5$$. Barcodes dropped in steps 1 and 2 get score 0 and label False.

The only randomness is the fold split, seeded with 18. The folds run in threads (`n_jobs` 2) but
independently, so a run is deterministic for fixed releases. The call still moves with scanpy
(`normalize_total`, `highly_variable_genes`, `scale`), scikit-image (`threshold_multiotsu` moved to
32-bit floats in 0.20) and scikit-learn (`StratifiedKFold`); by how much was not measured.

## Versions

| dropkick (tag date) | Change to the call |
|---|---|
| 1.0.0 (2020-03-11) | labels from Otsu on two metrics: $$\operatorname{arcsinh}$$ genes detected (cells above) and the percentage of counts in ambient genes (cells below); barcodes with fewer than 10 genes dropped; standardised features, folds scored by accuracy; 1,000 iterations, `n_jobs` −1 |
| 1.0.1 (2020-04-03) | `min_genes` 50 replaces the fixed 10 |
| 1.0.2 (2020-04-17) | `n_ambient` (10): the ambient genes are the 10 with the lowest dropout rate, ties at the boundary included |
| 1.1.0 (2020-04-25) | multi-Otsu on genes, Otsu on the ambient percentage; barcodes below the first multi-Otsu threshold left out of training |
| 1.1.2 (2020-05-13) | every barcode with at least `min_genes` genes trained on again |
| 1.1.4 (2020-06-05) | labels from genes alone by default (multi-Otsu on $$\operatorname{arcsinh}$$ genes detected, cells above); 2,000 iterations, `n_jobs` 2 |
| 1.1.5 (2020-06-16) | barcodes below the first threshold neither trained on nor scored (score 0) |
| 1.2.0 (2020-07-28) | `LogitNet(standardize=False, scoring="log_loss")`, where it standardised the features and scored by accuracy |
| 1.2.1 – 1.2.5 | output handling only; the label is the strings `"True"` / `"False"` from 1.2.3 |

1.1.3 and 1.2.6 – 1.2.8 are on PyPI without a git tag. The commits after 1.2.5 change none of the
functions of the call, so 1.2.6 – 1.2.8 should call as 1.2.5 (not checked against the PyPI files).

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `min_genes` | 50 | barcodes with fewer genes detected are dropped, never cells |
| `mito_names` | `^mt-` or `^MT-` | mitochondrial genes, for the QC columns only |
| `n_ambient` | 10 | ambient genes (lowest dropout rate), for the QC columns only by default |
| `n_hvgs` | 2000 | highly variable genes in the model |
| `metrics` | `["arcsinh_n_genes_by_counts"]` | the metrics the labels come from |
| `thresh_methods` | `["multiotsu"]` | per metric: `multiotsu`, `otsu`, `li` or `mean` |
| `directions` | `["above"]` | per metric: whether cells lie above or below the threshold |
| `alphas` | `[0.1]` | elastic-net mixing (the L1 share); several are compared if given |
| `max_iter` | 2000 | glmnet iterations |
| `seed` | 18 | the cross-validation split |
| `n_jobs` | 2 | threads; does not change the result |

The command line (`dropkick run`) exposes the same.

## Status

Documented only. dropkick is deterministic and could in principle be reproduced, but its
arithmetic is glmnet's Fortran coordinate descent. That code is GPL-2 (whether "or later" is not
stated, so it may not combine with GPL-3), and it runs on scanpy preprocessing whose results move
between scanpy releases: an exact port would mean porting glmnet. Its label also answers another
question, a good cell rather than a cell instead of an empty droplet.

Its own call is easy to set beside the package's. dropkick installs with pip and writes into the
AnnData it is given: `obs["dropkick_score"]` and `obs["dropkick_label"]`, with the coefficients in
`var["dropkick_coef"]` and the thresholds, the chosen alpha and lambda in `uns`. Run it on a copy
of the same raw AnnData, turn the label into a boolean `obs` column, and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

The labelling step alone is the package's [multi-Otsu](../generic/multiotsu.md) on
$$\operatorname{arcsinh}$$ genes detected, with a floor of 50 genes. Its top class should be
dropkick's training cells, up to ties at the threshold and the features counted (not checked
against dropkick); it is not dropkick's call:

```bash
sc-droplets call multiotsu raw_feature_bc_matrix.h5 --unit genes --transform arcsinh --min-counts 50
```

## References

- Heiser, Wang, Chen, Hughey & Lau 2021, Automated quality control and cell identification of
  droplet-based single-cell data using dropkick, *Genome Research* 31:1742–1752,
  [doi:10.1101/gr.271908.120](https://doi.org/10.1101/gr.271908.120).
- Friedman, Hastie & Tibshirani 2010, Regularization paths for generalized linear models via
  coordinate descent, *Journal of Statistical Software* 33:1–22,
  [doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01).
- Source: [KenLauLab/dropkick](https://github.com/KenLauLab/dropkick/tree/v1.2.5) at v1.2.5
  (`dropkick/api.py`, `dropkick/logistic.py`, `dropkick/util.py`, `LICENSE`, `LICENSES/`); PyPI
  [`dropkick`](https://pypi.org/project/dropkick/) 1.2.8.
- Research notes: docs/_research/docs-only.md §3.
