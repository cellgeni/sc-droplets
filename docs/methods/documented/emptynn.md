---
title: EmptyNN
parent: Documented only
grand_parent: Methods
nav_order: 1
method_ids: [emptynn]
status: docs-only
---

# EmptyNN: a neural network trained on positive and unlabelled barcodes
{: .no_toc }

EmptyNN takes the barcodes with few counts as surely empty and leaves the rest unlabelled. It
splits the unlabelled barcodes into folds, trains a small neural network to tell each fold from
the empty barcodes, and scores the other folds with it. A barcode whose mean score over several
rounds is above 0.5 is a cell. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `emptynn` (documented only) |
|---|---|
| Tool | EmptyNN 1.0 (R, GitHub `lkmklsmn/empty_nn` only; no tags or releases) |
| Category | classifier |
| Counts | every gene's counts, as proportions of the barcode's total |
| Needs expected cells | no |
| Deterministic | no: no seed is set for R's `sample`, `caret::createFolds` or TensorFlow |
| Original licence | none (`DESCRIPTION` keeps the template's "What license is it under?") |
| Extra inputs | none; the raw matrix with barcodes as rows. Keras and TensorFlow 2.8 |

## How it works

`emptynn(counts, threshold = 100, k = 5, iteration = 5, batch_size = 16, epoch = 10)`, on a raw
matrix with barcodes as rows (it stops if there are fewer rows than columns):

1. **Floor.** $$n$$ is each barcode's total over all genes. Barcodes with $$n < 10$$ are dropped.
2. **Empty and unlabelled.** The labelled set $$P$$, taken as empty: barcodes with
   $$n \le$$ `threshold` (100); a random 10,000 of them (`sample`) when there are 10,000 or more.
   The unlabelled set $$U$$: barcodes with $$n >$$ `threshold`.
3. **Features.** The 2,000 genes with the largest sums over $$P$$ (`tail(sort(colSums))`, so ties
   follow R's sort). Each barcode's counts are divided by its total $$n$$ over all genes; no log.
4. **Cross-validation.** In each of `iteration` (5) rounds, `caret::createFolds` splits $$U$$ into
   `k` (5) folds. For each fold $$j$$, a new network is trained on $$P$$ labelled 0 and fold $$j$$
   labelled 1, shuffled with `sample`, and predicts the other $$k - 1$$ folds. The round's column
   stores the softmax probability of class 1.
5. **Network.** Dense(256, ReLU), Dropout(0.4), Dense(128, ReLU), Dropout(0.3), Dense(2, softmax).
   Binary cross-entropy on one-hot labels; RMSprop at Keras's default rate; `batch_size` 16,
   `epoch` 10. Keras holds out the last 20 % of the shuffled rows (`validation_split = 0.2`); their
   loss is not used for anything.
6. **Call.** The score `mean.crossval` is the mean of the rounds' columns. A barcode of $$U$$ is a
   cell if its score is above 0.5. Barcodes at or below `threshold`, and those below 10 counts,
   are never cells.

Each round stores one prediction per barcode, not $$k - 1$$: the column is overwritten at every
fold, so a barcode keeps the prediction of the last model that tested it (fold $$k$$'s, or fold
$$k - 1$$'s for barcodes in fold $$k$$). The score is a mean over the rounds of those.

The code and the paper disagree:

- The paper and README say each barcode is predicted $$(k - 1) \times N$$ times and averaged; the
  code keeps one prediction per round, as above.
- The paper gives $$k = 10$$, $$N = 10$$ rounds and layers of 128, 64 and 2 neurons; the code has
  $$k = 5$$, 5 rounds, and layers of 256, 128 and 2 with dropout. The README says $$N$$ defaults
  to 10.
- The paper keeps barcodes with more than 10 counts; the code keeps those with at least 10. The
  paper does not mention the cap of 10,000 empty barcodes.

The randomness is R's `sample` (the subsample of $$P$$, the shuffling), `caret::createFolds`,
Keras's weight initialisation, dropout and the minibatch order. Nothing is seeded, neither with
`set.seed` nor with `tensorflow::set_random_seed`, and the paper mentions no seed. Two runs give
different calls unless the caller seeds both R and TensorFlow.

## Versions

| EmptyNN | Change to the call |
|---|---|
| 880bd51 (2020-11-01, "add package file") | `emptynn` as above; version 1.0 throughout |
| 660e9fa (2022-03-23) | `keras` and `tensorflow` pinned to 2.8.0 (`keras (>= 2.3.0)` before); a copy of the package under `EmptyNN/`; the R code unchanged |
| e1e6226 (2024-01-04, the last commit) | no change to `R/emptynn.R` or `R/neg_create_model.R` |

The repository also bundles `EmptyNN_1.0.tar.gz`, whose R files are the same.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `threshold` | 100 | totals at or below it are taken as empty; above it, unlabelled |
| `k` | 5 | folds of the unlabelled barcodes |
| `iteration` | 5 | rounds of cross-validation |
| `batch_size` | 16 | minibatch size |
| `epoch` | 10 | epochs per network |

Fixed in the code: the floor of 10 counts, the cap of 10,000 empty barcodes, the 2,000 genes, the
architecture, RMSprop, the validation split of 0.2 and the cut at 0.5.

## Status

Documented only. The code has no licence, so it cannot be ported into this GPL-3 package (the
article's CC BY-NC-ND 4.0 licence covers its text, not the code). The call is a set of neural
networks trained per run with Keras and TensorFlow 2.8 from R, and nothing is seeded: the original
does not repeat its own call, so there is no fixed output to validate against.

Its own result can still be set beside the package's calls. `emptynn` returns a list:
`nn.keep`, a logical vector over `rownames(counts)`, every input barcode; and `prediction`, a data
frame over the unlabelled barcodes with one column per round (`iteration1`, …) and the score
`mean.crossval`. Saved with the barcodes, `nn.keep` becomes a boolean `obs` column of the same
raw matrix read with `scd.read`, and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

## References

- Yan, Zhao & Simon 2021, EmptyNN: A neural network based on positive and unlabeled learning to
  remove cell-free droplets and recover lost cells in scRNA-seq data, *Patterns* 2:100311,
  [doi:10.1016/j.patter.2021.100311](https://doi.org/10.1016/j.patter.2021.100311).
- Source: [lkmklsmn/empty_nn](https://github.com/lkmklsmn/empty_nn/tree/e1e6226) at e1e6226
  (`R/emptynn.R`, `R/neg_create_model.R`, `DESCRIPTION`).
- Research notes: docs/_research/docs-only.md §2.
