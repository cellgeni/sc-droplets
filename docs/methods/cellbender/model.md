---
title: CellBender model
parent: CellBender
grand_parent: Methods
nav_order: 2
method_ids: [cellbender-model]
status: docs-only
---

# CellBender model: the posterior cell probability
{: .no_toc }

What `cellbender remove-background` itself calls a cell. After its [priors](prior.md), it trains
a variational autoencoder of the raw counts in which every droplet either holds a cell or not,
and it reports each analysed droplet's posterior probability of not being empty. Droplets above
0.5 are its cells. The priors are implemented; the model and its probability are documented, not
implemented.

1. TOC
{:toc}

## At a glance

| Id | `cellbender-model` (documented only) |
|---|---|
| Tool | CellBender 0.1.0 – 0.4.0 (`cellbender remove-background`); 0.3.0_rc, the paper's code, was never on PyPI |
| Category | classifier |
| Counts | the counts of every feature not excluded, in the droplets the priors analyse |
| Needs expected cells | no (optional; it feeds the priors) |
| Deterministic | no: seeded, but minibatches, dropout and the Monte Carlo ELBO depend on torch, Pyro, the thread count and the GPU |
| Original licence | BSD-3-Clause |
| Extra inputs | none; the raw matrix. A GPU in practice |

## How it works

This is 0.3.0 and later; the earlier releases are under [Versions](#versions). From the
[priors](prior.md): $$T$$ droplets are analysed, the top $$T$$ by total; $$C$$ is the cell prior,
$$E$$ the empty prior and $$E_c$$ the expected cells; $$t$$ is the largest surely-empty count (the
total at rank $$T$$); $$c$$ is the log crossover. For a droplet, $$n$$ is its total over the
analysed features and $$s = \log(1 + n)$$.

1. **Input.** The raw matrix over every feature not excluded (`--exclude-feature-types`,
   `--ignore-features`), less the features whose projected ambient count in cells is below 0.1.
   Barcodes outside the top $$T$$ are never cells.
2. **Model** (`--model full`). Per droplet: a latent code $$z \sim N(0, I_{64})$$ decoded to an
   expression profile $$\chi$$; a cell indicator $$y \in \{0, 1\}$$; a cell size
   $$d_c \sim \mathrm{LogNormal}(\log(1 + C), \sigma_d)$$ and an empty size
   $$d_e \sim \mathrm{LogNormal}(\log(1 + E), 0.01)$$; an efficiency
   $$\epsilon \sim \mathrm{Gamma}(50, 50)$$; a barcode-swapping fraction
   $$\rho \sim \mathrm{Beta}(1.5, 50)$$. The counts are a negative binomial (the cell) convolved
   with a Poisson (the background), of mean $$(1-\rho)\,y\,\epsilon\,d_c\,\chi$$ and rate
   $$\epsilon[(1-\rho)\chi_a d_e + \rho\,\bar\chi\,(y d_c + d_e)]$$, where $$\chi_a$$ is the
   ambient profile and $$\bar\chi$$ the average profile.
3. **Prior on $$y$$.** Logit $$-100$$ when $$\log n \le \tfrac12(\log t + \log(1 + C))$$, which is
   essentially the crossover; $$+10$$ when $$\log n \ge \log(1 + C)$$; in between, the prior cell
   logit $$\mathrm{logit}\big((E_c / T)(1 - f)\big)$$, with $$f$$ = `--empty-drop-training-fraction`
   (0.2).
4. **The logit of $$y$$.** The encoder takes the counts, normalised, `log1p`-transformed,
   batch-normalised and passed through `Dropout1d(0.5)`, and four features: $$s$$,
   $$\log(1 + \text{genes detected})$$, an overlap $$-\tfrac12\big(\max(0, s - \mu_e)/0.1\big)^2$$
   with $$\mu_e$$ the learnt location of the empty size, and the norm of the $$z$$ encoder's mean.
   Two layers of 512 units (linear, batch norm, softplus), each fed the four features again, and
   a linear output $$o$$. Then
   $$\ell = (o - o_0) + 10\,\mathrm{sign}(s - c)\sqrt{|s - c|},$$
   where $$o_0$$ is the mean of $$o$$ over the first training minibatch (which must hold more than
   4 droplets with $$s > c$$). Two blends follow: with $$a = \sigma(50(s - \log(1 + C)))$$,
   $$\ell \leftarrow (1 - a)\ell + 10a$$; with $$b = \sigma(50(s - \log(1 + t)))$$,
   $$\ell \leftarrow b\ell - 10(1 - b)$$. Droplets well above the cell prior get $$\ell \approx 10$$,
   droplets below the surely-empty count $$\ell \approx -10$$.
5. **Training.** Stochastic variational inference (Pyro's `TraceEnum_ELBO`, with $$y$$ summed
   out) for 150 epochs on 90 % of the droplets, in minibatches of $$\min(256, 0.9T/2)$$ of which
   a fifth are surely-empty droplets drawn with replacement. ClippedAdam (clip norm 10), learning
   rate $$10^{-4}$$ under a one-cycle schedule that peaks at 10 times it. Added to the ELBO, for
   every droplet: a normal term (sd 2, weight 1) holding $$\ell$$ near the prior cell logit, and
   one (sd 1, weight 10) pulling it to $$-10$$ when $$n < e^c$$ or to $$+10$$ when $$n \ge e^c$$.
   The ambient profile is also fitted to the droplets the guide samples as empty (Poisson, weight
   0.01), and the medians of $$\epsilon$$ among cells and among empty droplets are held near 1.
6. **Probability.** After training, in evaluation mode (batch norm on its running statistics, no
   dropout), $$p = \sigma(\ell)$$ for each of the $$T$$ analysed droplets, in batches of 500.
7. **Call.** A cell if $$p > 0.5$$, a constant (`CELL_PROB_CUTOFF`). The barcode list and the
   filtered matrix use $$>$$; the metrics file counts cells with $$\ge$$. `--fpr` does not change
   the call: it sets only how many counts are removed.

A reading of these terms, not measured: with $$-100$$ below the crossover, the supervision of
$$\pm 10$$ on either side of it and the signed square root, the call is anchored to "total above
the crossover $$e^c$$". The network mostly decides the droplets above it that look ambient, and
those near it; above the cell prior the blend forces $$p$$ to $$\sigma(10)$$. CellBender's
documentation says as much in words: the probability is that a droplet is not empty, and many
low-quality droplets are expected among the calls.

`pyro.set_rng_seed(1234)` seeds Python, NumPy and torch. Runs can still differ through the
train/test split and minibatch order, the empty droplets drawn into each batch, dropout, the
Monte Carlo ELBO, torch's CPU thread count (`--cpu-threads`) and the GPU's kernels. One user saw 30
of 11,001 cells differ between two runs of 0.3.0, three months apart, in an unchanged environment
(issue #401).

## Versions

| CellBender | Change to the call |
|---|---|
| 0.1.0 | its own network for $$p$$ (layers 100, 10) on the log total, the normalised counts and the counts less the ambient, started as a sigmoid of the log total around the crossover; prior logit $$\mathrm{logit}(0.5\,E_c/T)$$ for every droplet; no supervision; the guide samples $$z$$ only where the logit is $$\ge 0$$; negative-binomial counts; $$z$$ of 20, constant learning rate $$10^{-3}$$, half of each batch empty; seed 0; NaN probabilities set to 0 |
| 0.2.0 – 0.2.2 | one network (100, 50) for $$p$$, $$d$$ and $$\epsilon$$, on the log total, log genes detected, a scaled Poisson log-likelihood under the ambient profile and the scaled normalised counts; the prior logit is the prior cell logit less 100 for every droplet (marked `# TODO`), so cells come from the likelihood and the supervision: $$N(-5, 1)$$, weight 1, below the surely-empty count, and $$N(+5, 1)$$, weight 10, at or above the cell prior; $$z$$ of 100, one-cycle learning rate from $$10^{-4}$$; seed 0 |
| 0.3.0_rc (the paper's code) | as 0.2, with an output scale of 2 and the offset at the 25th percentile of the cell-like outputs; also $$N(\mp 5, 1)$$, weight 0.1, either side of the crossover, $$N(-5, 1)$$, weight 10, for surely-empty droplets, and $$N(\mp 5, 0.2)$$, weight 10 or 1, when a surely-empty droplet's logit is $$\ge 0$$ or a surely-cell droplet's $$\le 0$$; priors from a pyro Gaussian mixture; a fifth of each batch empty; seed 1234 |
| 0.3.0 | the model above |
| 0.3.2, 0.4.0 | no change to the call: 0.3.2 changes only the ELBO-failure checks (used only with `--epoch-elbo-fail-fraction` or `--final-elbo-fail-fraction`), 0.4.0 formatting and types |

In every release the cut is $$p > 0.5$$ and only analysed droplets can be cells, so the call also
moves with the [priors](prior.md#versions), which changed in 0.2 and 0.3.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--expected-cells` | estimated | the expected cells of the priors; sets the prior cell logit |
| `--total-droplets-included` | estimated | $$T$$, the droplets analysed; the others are never cells |
| `--force-cell-umi-prior`, `--force-empty-umi-prior` | estimated | the cell and empty priors |
| `--low-count-threshold` | 5 | barcodes at or below it are left out of the priors |
| `--model` | `full` | also `naive`, `simple`, `ambient`, `swapping`; `naive` trains no network and has no cell probability (what it then calls was not checked) |
| `--epochs` | 150 | training epochs |
| `--learning-rate` | $$10^{-4}$$ | the base rate of the one-cycle schedule (none with `--constant-learning-rate`) |
| `--z-dim`, `--z-layers` | 64, 512 | the latent code and the encoder's hidden layer |
| `--training-fraction` | 0.9 | droplets trained on; the rest are a test set |
| `--empty-drop-training-fraction` | 0.2 | share of surely-empty droplets in each minibatch; enters the prior cell logit |
| `--projected-ambient-count-threshold` | 0.1 | features with fewer projected ambient counts in cells are dropped |
| `--exclude-feature-types`, `--ignore-features` | none | features left out |
| `--num-training-tries` | 1 | retries after a failed ELBO check, each at the rate times `--learning-rate-retry-mult` (0.2; 0.5 in 0.2.x) |
| `--cuda` | off | train on the GPU |
| `--fpr` | 0.01 | counts removed; does not change the call |

## Status

Documented only. The call is the output of a Pyro variational autoencoder trained by stochastic
variational inference: hours on a CPU for a full sample (CellBender recommends a GPU), dependent
on the torch and Pyro releases and the hardware, with no claim to repeat bit for bit: there is no
noise band to hold a reimplementation to. The deterministic part, the priors, is
[`cellbender-prior`](prior.md); its `cell` class is not this call, but its analysed droplets bound
it.

CellBender's own call can still be set beside the package's. From 0.3.0, `out.h5` holds
`/droplet_latents/cell_probability`, one value per analysed barcode, and
`/droplet_latents/barcode_indices_for_latents`, each one's index into the raw barcode list;
`out_cell_barcodes.csv` lists the barcodes above 0.5, and
`cellbender.remove_background.downstream.anndata_from_h5` puts `cell_probability` in `obs`. In 0.1
and 0.2 the output HDF5 holds `latent_cell_probability` with `barcode_indices_for_latents`, and
`_cell_barcodes.csv` the cells. Mark those barcodes in a boolean `obs` column of the raw matrix
read with `scd.read`, and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls. The release is not in
the HDF5; `out.log` starts with "CellBender" and the version.

## References

- Fleming, Chaffin, Arduini, Akkad, Banks, Marioni, Philippakis, Ellinor & Babadi 2023,
  Unsupervised removal of systematic background noise from droplet-based single-cell experiments
  using CellBender, *Nature Methods* 20:1323–1335,
  [doi:10.1038/s41592-023-01943-7](https://doi.org/10.1038/s41592-023-01943-7).
- Source: [`cellbender/remove_background`](https://github.com/broadinstitute/CellBender/tree/v0.3.0/cellbender/remove_background)
  at v0.3.0 (`model.py`, `vae/encoder.py`, `train.py`, `posterior.py`, `run.py`, `consts.py`),
  the same in the call at v0.3.2 and v0.4.0; the tags
  [v0.3.0_rc](https://github.com/broadinstitute/CellBender/tree/v0.3.0_rc),
  [v0.2.2](https://github.com/broadinstitute/CellBender/tree/v0.2.2) and
  [v0.1.0](https://github.com/broadinstitute/CellBender/tree/v0.1.0).
- [CellBender's documentation](https://cellbender.readthedocs.io/), Troubleshooting; issue
  [#401](https://github.com/broadinstitute/CellBender/issues/401) on repeated runs.
- Research notes: docs/_research/docs-only.md §1.
