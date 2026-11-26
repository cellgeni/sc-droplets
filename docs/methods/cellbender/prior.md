---
title: priors
parent: CellBender
grand_parent: Methods
nav_order: 1
method_ids: [cellbender-prior]
status: implemented
---

# CellBender priors
{: .no_toc }

What `cellbender remove-background` decides before it trains its model: prior counts for cells and
for empty droplets, how many droplets are probably cells, how many it analyses, and which empty
droplets it learns the ambient profile from. Each barcode gets one of four classes: `cell` (the
probable cells), `transition` (the other analysed droplets, which the model decides), `empty`
(the surely-empty droplets it trains on) and `excluded` (never looked at). Only `cell` is called.

This is not CellBender's cell call. Its cells are those whose posterior probability after training
is above 0.5, and from 0.3.0 the expected cells are deliberately conservative (CellBender's own
comment: "expected cells does well when you give it a very conservative estimate").

1. TOC
{:toc}

## At a glance

| Id | `cellbender-prior` |
|---|---|
| Tool | CellBender 0.1.0 and later (`cellbender remove-background`, `SingleCellRNACountsDataset`) |
| Category | density |
| Counts | the totals over every feature (CellBender analyses every feature unless told to exclude some; `--feature-type` selects) |
| Needs expected cells | no (optional) |
| Deterministic | yes |
| Original licence | BSD-3-Clause |
| Validated against | CellBender 0.1.0, 0.2.2, 0.3.0, 0.3.2 and 0.4.0 (exact) |

## How it works

$$u$$ is each barcode's total count (an integer), $$L$$ the low-count threshold, $$E$$ the expected
cells and $$T$$ the total droplets analysed. "Sorted" means in decreasing order, and rank $$r$$
(from 0) is the $$(r+1)$$-th largest.

### 0.3.0 and later

**Cell and empty priors** (`get_cell_count_empty_count`):

1. Otsu's threshold (256-bin histogram, scikit-image's fast form: the centre of the bin that
   maximises the between-class variance) of $$\log u$$ over the first
   $$\lfloor \#\{u > L\} / 4 \rfloor$$ sorted totals gives the cells' lower limit $$c_\text{low}$$.
   The first cell prior is the mean of $$u$$ over $$u > c_\text{low}$$.
2. The empty prior comes from a kernel density, refined up to 5 times until it moves by no more
   than 10. The sorted totals above $$L$$ are $$k$$; the first cutoff is
   $$\max(0.1\,c_\text{low},\ k_{-100})$$, the 100th smallest of them. Each round:
   - Of the $$k$$ at or below the cutoff ($$m$$ of them), take 1000 at log-spaced ranks
     $$\lfloor 10^{t} - n_c \rfloor$$ clamped to $$[0, m-1]$$, for 1000 values of $$t$$ evenly
     spaced from $$\log_{10} n_c$$ to $$\log_{10}(n_c + m)$$, where $$n_c = \#\{k > c_\text{low}\}$$.
     This resampling makes the histogram look like the log-log rank plot.
   - Fit scipy's `gaussian_kde` (Scott's bandwidth) to their natural logs and evaluate it on a grid
     of step 0.1 from $$\lfloor \min \rfloor - 0.01$$ to $$\lceil \max \rceil + 0.01$$.
   - The empty prior is $$e^{x}$$ at the density's maximum. Its upper limit is $$e^{x}$$ at the
     first grid point at or after the peak where the density is below 0.6 of the peak (else the
     last grid point).
   - The next cutoff is $$0.75\sqrt{\text{cell prior} \times \text{empty prior}}$$.
3. The cell prior becomes the median of $$u$$ over
   $$u > \sqrt{\text{cell prior} \times \text{empty prior}}$$.

**Expected cells and total droplets** (`get_expected_cells_and_total_droplets`):
$$E = \#\{u \ge \text{cell prior}\}$$, and
$$T = \#\{u \ge \sqrt{\text{empty prior} \times \text{upper limit}}\}$$, but at most
$$\#\{u \ge \sqrt{\text{cell prior} \times \text{empty prior}}\} + 20{,}000$$. If $$T > 70{,}000$$,
steps 1 to 3 and this step run again with $$L$$ set to the empty prior's upper limit, at most six
times.

**Given values** replace the estimates, in this order:

- `--expected-cells` $$E$$: the cell prior becomes $$\exp(\text{mean of } \log u)$$ over the $$E$$
  largest. If $$E + 20{,}000 > T$$, then $$T = E + 20{,}000$$, and, unless $$T$$ is given too, the
  empty prior becomes the median of the sorted totals at ranks $$[\max(E, T - 500), T)$$, with
  upper limit $$\min(\sqrt{c' \times \text{empty}},\ 1.5 \times \text{empty})$$. Here $$c'$$ is
  the geometric mean of the $$E$$ largest totals.
- `--total-droplets-included` $$T$$: the empty prior and its upper limit are re-estimated the same
  way.
- `--force-cell-umi-prior` sets the cell prior. `--force-empty-umi-prior` sets the empty prior,
  with upper limit $$\min(\sqrt{\text{cell} \times \text{empty}},\ 2 \times \text{empty})$$.

Then the surely-empty count is the sorted total at rank $$T$$. The log crossover is the mean of
its log and the log of the cell prior. $$d_\text{std}$$ is a fifth of the standard deviation of
$$\log u$$ over the barcodes above the crossover.

**Features.** CellBender drops the features with almost no ambient counts. A feature stays if its
mean count over the droplets with $$L < u <$$ upper limit, times $$\max(E, 5000)$$, is above 0.1
(`--projected-ambient-count-threshold`). $$u'$$ is each barcode's total over the features that
stay.

**Droplets.** The barcodes are ranked by $$u'$$, decreasing, as numpy's
`argsort(u')[::-1]` orders them:

- ranks $$[0, E)$$ are `cell`;
- ranks $$[E, T)$$ are `transition`;
- the remaining barcodes with $$u' > \max(L, \lfloor 0.5 \times \text{empty prior} \rfloor)$$
  are `empty`;
- every other barcode is `excluded`.

The largest surely-empty count is $$u'$$ at rank $$T$$. The prior cell logit is
$$\operatorname{logit}\!\left(\tfrac{E}{T}(1 - 0.2)\right)$$ (0.2: `--empty-drop-training-fraction`).

The model then gives each analysed droplet a prior logit: $$-100$$ when
$$\log u' \le \tfrac{1}{2}(\log \text{surely-empty} + \log(1 + \text{cell prior}))$$, $$+10$$
when $$\log u' \ge \log(1 + \text{cell prior})$$, and the prior cell logit in between. Training
starts from there.

### 0.1.0 and 0.2.x

1. Empty prior: $$\lfloor e^{m} - 1 \rfloor$$, where $$m$$ is the most frequent value (the smallest,
   if several) of $$\log(1 + u)$$ rounded to one decimal, over $$u > L$$.
2. Cell prior: $$\lfloor \exp(\text{median of } \log(1 + u) \text{ over } u > 5 \times
   \text{empty}) - 1 \rfloor$$. From 0.2.0, with `--expected-cells` $$E$$: the median of the
   $$E$$ largest totals, truncated to an integer.
3. Expected cells: $$E$$ if given, else the number of barcodes above
   $$\exp(\tfrac{1}{2}(\log(1 + \text{cell}) + \log(1 + \text{empty}))) - 1$$; $$n$$ is that
   number, but at most the number of barcodes with $$u > 0$$.
4. Droplets, ranked by $$u$$ as numpy's `argsort(u)[::-1]` orders them:
   - ranks $$[0, n)$$ are `cell`;
   - the next $$\min(T - n,\ N - n)$$ barcodes are `transition`, where $$N$$ is the number of
     barcodes with $$u > \max(L, \lfloor f \times \text{empty} \rfloor)$$, $$f = 0.8$$ in 0.1.0
     and $$0.5$$ from 0.2.0;
   - the rest of those $$N$$ are `empty`, but only when the transition got all $$T - n$$
     barcodes;
   - every other barcode is `excluded`.
5. The log crossover is the mean of $$\log(1 + \text{cell})$$ and $$\log(1 + \text{empty})$$;
   $$d_\text{std} = (\log(1 + \text{cell}) - \text{crossover}) / 5$$. The prior cell logit is
   $$\operatorname{logit}(0.5\,E / \#\text{analysed})$$, where this $$E$$ is the given or
   estimated count before the cap at the barcodes with $$u > 0$$. From 0.2.0, the largest
   surely-empty count is $$u$$ at the first rank after the analysed droplets.

CellBender stops with an error when there are no cells, when no barcode above the low cut is left
for the empty droplets, when $$T$$ is not above the number of cells, or (before 0.3) when no
transition droplet is left; so does this.

## Versions

| CellBender | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.1.0 | `v0.1` | `0.1.0` | mode and log-median priors; low cut 0.8 × empty prior; $$T$$ = 25,000; $$L$$ = 15 |
| 0.2.0 – 0.2.2 | `v0.2` | `0.2.0` | low cut 0.5 × empty prior; median of the top $$E$$ when $$E$$ is given |
| 0.3.0 – 0.4.0 | `v0.3` | `0.3.0` | Otsu and kernel-density priors; $$E$$ and $$T$$ estimated; feature trimming; $$L$$ = 5 |

0.3.0_rc (the paper's code, a pyro Gaussian mixture) is not offered.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expected-cells` | estimated | `--expected-cells` | the number of probable cells |
| `--total-droplets-included` | 25,000 (0.1, 0.2), estimated (0.3) | `--total-droplets-included` | the number of droplets analysed, cells included |
| `--low-count-threshold` | 15 (0.1, 0.2), 5 (0.3) | `--low-count-threshold` | barcodes at or below this count are left out of the priors and the empty droplets |
| `--force-cell-umi-prior` | none | `--force-cell-umi-prior` | the cell prior itself (0.3) |
| `--force-empty-umi-prior` | none | `--force-empty-umi-prior` | the empty prior itself (0.3) |
| `--empty-drop-training-fraction` | 0.5 (0.1, 0.2), 0.2 (0.3) | `--empty-drop-training-fraction` | enters the prior cell logit only |
| `--projected-ambient-count-threshold` | 0.1 | `--projected-ambient-count-threshold` | the feature trimming's cut on projected ambient counts in cells (0.3) |

## Run it

```bash
sc-droplets call cellbender-prior raw_feature_bc_matrix.h5 --json
sc-droplets call cellbender-prior raw_feature_bc_matrix.h5 --tool-version 0.2.2 --expected-cells 5000
sc-droplets call cellbender-prior raw_feature_bc_matrix.h5 --feature-type "Gene Expression" -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "cellbender-prior")
result.thresholds["cell_counts"], result.thresholds["total_droplets"]
adata.obs["cellbender_prior_class"].value_counts()
```

The thresholds hold the priors under CellBender's names: `cell_counts`, `empty_counts`,
`expected_cells`, `total_droplets`, `low_count_cutoff`, `log_counts_crossover`, `d_std` and
`cell_logit`. From 0.2 they also hold `empty_umi_threshold`, the largest surely-empty count. From
0.3 they add `empty_count_upper_limit`, `surely_empty_counts`, `cell_count_low_limit` (Otsu's) and
`n_features`, the number of features kept. `sc-droplets plot diagnostics` draws the rank curve
with the four classes and the priors, and the histogram the priors came from, with the kernel
density from 0.3.

These numbers are what CellBender writes to its log ("Prior on counts for cells is …", "Using …
probable cell barcodes, plus an additional … barcodes, and … empty droplets").

## Differences from the original

- CellBender 0.2's `--exclude-antibody-capture` leaves antibodies out of the priors but not out of
  the ranking; here the selected features are used throughout. `--ignore-features`
  (`--blacklist-genes`) and `--model simple` are not offered.
- Ties are ranked as numpy's classic quicksort leaves them, which is what numpy before 1.25 does
  everywhere. Later numpy releases on some x86-64 processors (with AVX-512, or AVX2 in newer
  releases) sort with SIMD code that leaves equal totals in another order. CellBender run there can then put
  different barcodes at the class boundaries $$E$$ and $$T$$, all of them with the same total
  as the barcodes here.
- CellBender reads an `.h5ad` as it is stored, often float32, and then computes in float32. Here
  the totals are summed in float64, and treated as integers when they are whole.

## Validation

`oracles/cellbender/run.sh` builds each release's `SingleCellRNACountsDataset`, from its tagged
source in its own environment (torch, pyro and numpy before 1.25), on 10x's raw pbmc_1k_v3
(Cell Ranger 3.0) and pbmc4k (Cell Ranger 2.1) matrices. The cases are the defaults,
`--expected-cells`, `--total-droplets-included`, `--low-count-threshold 30`, and, from 0.3,
forced priors and $$E$$ with $$T$$. The barcodes of every class and every count match exactly.
The priors that come from `log` and `exp` match to the last bits. CellBender's own values
depend on the processor here: on Linux with AVX-512, numpy computes float64 `log` and `exp`
with Intel's SVML, which can differ from the C library by a unit in the last place. A Linux
runner and macOS gave pbmc4k's empty prior as 36.234075926476564 and 36.23407592647656.

## References

- Fleming, Chaffin, Arduini, Akkad, Banks, Marioni, Philippakis, Ellinor & Babadi 2023,
  Unsupervised removal of systematic background noise from droplet-based single-cell experiments
  using CellBender, *Nature Methods* 20:1323–1335,
  [doi:10.1038/s41592-023-01943-7](https://doi.org/10.1038/s41592-023-01943-7).
- Source: [`priors.py`](https://github.com/broadinstitute/CellBender/blob/v0.3.0/cellbender/remove_background/data/priors.py)
  and [`dataset.py`](https://github.com/broadinstitute/CellBender/blob/v0.3.0/cellbender/remove_background/data/dataset.py)
  (0.3.0, unchanged in 0.3.2 and 0.4.0);
  [`dataset.py` at 0.2.2](https://github.com/broadinstitute/CellBender/blob/v0.2.2/cellbender/remove_background/data/dataset.py)
  and [at 0.1.0](https://github.com/broadinstitute/CellBender/blob/v0.1.0/cellbender/remove_background/data/dataset.py).
