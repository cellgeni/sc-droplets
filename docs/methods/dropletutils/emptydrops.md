---
title: emptyDrops
parent: DropletUtils
grand_parent: Methods
nav_order: 2
method_ids: [dropletutils-emptydrops]
status: implemented
---

# DropletUtils emptyDrops
{: .no_toc }

EmptyDrops (Lun et al. 2019) estimates the ambient RNA profile from the droplets with the
fewest counts, then tests every other droplet against it: a droplet whose counts are unlikely to
have been drawn from that profile holds a cell. Droplets at or above the knee of the barcode-rank
curve are kept whatever the test says. The test was a multinomial in 1.0, a Dirichlet-multinomial
with an estimated dispersion from 1.2 to 1.26, and a multinomial again from 1.28; the knee
changed in 1.10, 1.28 and 1.30.

1. TOC
{:toc}

## At a glance

| Id | `dropletutils-emptydrops` |
|---|---|
| Tool | DropletUtils 1.0 and later (`emptyDrops()`) |
| Category | statistical-test |
| Counts | the matrix as given (UMIs), every feature with a count |
| Needs expected cells | no |
| Deterministic | yes, given R's seed (default 100) |
| Original licence | GPL-3 |
| Validated against | `DropletUtils::emptyDrops` in R with the same seed (exact) |

## How it works

1. **Counts.** From 1.6, a matrix whose row or column sums are not whole numbers (by
   `all.equal`) is rounded; from 1.12 it is always rounded (`round = TRUE`). Totals
   $$t_b$$ are `as.integer(round(colSums(m)))`. Features with no count anywhere are dropped.
2. **Ambient droplets.** Those with $$t_b \le$$ `lower` (100). With `by.rank = k`, `lower`
   becomes the $$(k+1)$$-th largest total (1.10 – 1.18); from 1.20, the ambient droplets are
   exactly those ranked below $$k$$ (ties ranked by position) and `lower` the largest of their
   totals.
3. **Ambient profile.** Sum the ambient droplets' counts per feature, then edgeR's
   `goodTuringProportions`: a feature seen $$r$$ times gets the Simple Good-Turing probability of
   an item seen $$r$$ times, and the features seen 0 times share $$P_0 = n_1 / N$$ equally. edgeR's
   estimator switches from the Turing to the smoothed estimate at the first count whose successor
   is missing or where the two agree within 1.96 standard deviations, and keeps the smoothed one
   after that. From 1.8, if any feature still has probability 0 (no singletons), they share
   $$1/\sum$$ (the reciprocal of the ambient total) and the rest are scaled by $$1 - 1/\sum$$.
4. **Dispersion.** With `alpha` estimated (the default from 1.2 to 1.26), $$\alpha$$ maximises the
   Dirichlet-multinomial likelihood of the ambient droplets,
   $$\ell(\alpha) = n\log\Gamma(\alpha) - \sum_b \log\Gamma(t_b + \alpha) + \sum_{b,g:\ x_{bg} > 0}
   [\log\Gamma(x_{bg} + \alpha a_g) - \log\Gamma(\alpha a_g)]$$, found by R's `optimize` on
   (0.01, 10000) (Brent's method, tolerance $$\varepsilon^{1/4}$$).
5. **Tested droplets.** Those with $$t_b >$$ `lower` (with `test.ambient`, every droplet with
   $$t_b > 0$$), and $$t_b >$$ `ignore` if given.
6. **Observed statistic.** The data-dependent part of the log-probability,
   $$L_b = \sum_{g:\ x_{bg} > 0} [x_{bg} \log a_g - \log x_{bg}!]$$ (multinomial), or
   $$\sum_g [\log\Gamma(x_{bg} + \alpha a_g) - \log x_{bg}! - \log\Gamma(\alpha a_g)]$$
   (Dirichlet-multinomial). The `LogProb` column adds $$\log t_b!$$ (and
   $$\log\Gamma(\alpha) - \log\Gamma(t_b + \alpha)$$).
7. **Monte Carlo.** Each of `niters` (10,000) iterations draws one path of features from the
   ambient profile (from a profile drawn from the Dirichlet, with finite $$\alpha$$), one feature
   at a time up to the largest tested total, accumulating the same statistic:
   $$+\log a_g$$ (or $$+\log(\alpha a_g + c_g)$$), then $$-\log(c_g + 1)$$, with $$c_g$$ the times
   $$g$$ was drawn before. For each tested droplet, $$n_b$$ counts the iterations whose value at
   $$t_b$$ draws is at most $$L_b$$.
8. **p-values.** $$p_b = (n_b + 1)/(\text{niters} + 1)$$; `Limited` marks $$n_b = 0$$.
9. **Retained droplets.** `retain` is the knee of [barcodeRanks](barcoderanks.md) of the same
   release, called with the same `lower`; droplets with $$t_b \ge$$ `retain` get $$p = 0$$.
10. **FDR.** Benjamini-Hochberg over the tested droplets (`p.adjust`). From 1.12, with
    `test.ambient = TRUE`, droplets with $$t_b \le$$ `lower` get FDR NA and do not count towards
    the adjustment (`test.ambient = NA` keeps them, as before 1.12).

DropletUtils returns the FDR and leaves the cut to the user; its documentation and the OSCA book
use FDR $$\le 0.001$$, which is the default `--fdr` here.

### Random numbers

What makes the p-values reproducible is R's random number generator, seeded by `set.seed` before
`emptyDrops` is called (sc-droplets' `--seed`, default 100):

| DropletUtils | Seeding | Draws | Dirichlet profile |
|---|---|---|---|
| 1.0 – 1.2 | R's Mersenne-Twister, one stream for every iteration | `unif_rand()` times the profile's sum, inverted through the cumulative profile | R's `rgamma` |
| 1.4 | iteration $$i$$: `pcg32(seed_i, i)`, `seed_i` from dqrng's `generateSeedVectors` (two 32-bit words from R's `R_unif_index`) | Boost `uniform_real_distribution(0, sum)`, inverted through the cumulative profile | Boost `gamma_distribution` |
| 1.6 on | as 1.4 | Boost `discrete_distribution` (alias table: `uniform_int` for the column, `uniform_01` for the test) | Boost `gamma_distribution` |

sc-droplets reproduces each piece bit for bit: R's seeding of the Mersenne-Twister, dqrng's seed
vectors, PCG32's XSH-RR output, Boost's integer buckets, alias table, ziggurat exponential and
gamma rejection samplers. With the seed R was given, the p-values are R's.

## Versions

| DropletUtils | Variant | `--tool-version` | Model by default | Knee for `retain` | Changes |
|---|---|---|---|---|---|
| 1.0 | `v1.0` | `1.0` | multinomial | v1.0 | |
| 1.2 | `v1.2` | `1.2` | Dirichlet-multinomial | v1.0 | `alpha` |
| 1.4 | `v1.4` | `1.4` | Dirichlet-multinomial | v1.0 | PCG32 per iteration; error if the ambient profile is empty |
| 1.6 | `v1.6` | `1.6` | Dirichlet-multinomial | v1.0 | alias-table draws; non-integer matrices rounded |
| 1.8 | `v1.8` | `1.8` | Dirichlet-multinomial | v1.0 | `round`; pseudo-count for zero probabilities |
| 1.10 | `v1.10` | `1.10` | Dirichlet-multinomial | v1.10 | `by.rank` (`retain` still uses the given `lower`) |
| 1.12 – 1.18 | `v1.12` | `1.12` | Dirichlet-multinomial | v1.12 | always rounds; `test.ambient` FDR rule; `retain` uses the `by.rank` lower |
| 1.20 – 1.26 | `v1.20` | `1.20` | Dirichlet-multinomial | v1.12 | `by.rank` by exact rank |
| 1.28 | `v1.28` | `1.28` | multinomial | v1.28 | `alpha = Inf` by default |
| 1.30 – 1.34 | `v1.30` | `1.30` | multinomial | v1.30 | |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--lower` | 100 | `lower` | droplets with this total or less are ambient |
| `--niters` | 10000 | `niters` | Monte Carlo iterations |
| `--test-ambient` | `false` | `test.ambient` | `true` tests the ambient droplets too; `na` (1.12 on) also keeps them in the FDR |
| `--ignore` | none | `ignore` | droplets with this total or less are not tested |
| `--alpha` | the release's | `alpha` | `estimate`, `inf` (multinomial) or a number; 1.0 has no `alpha` |
| `--retain` | the knee | `retain` | droplets with at least this total are kept untested |
| `--by-rank` | none (1.10 on) | `by.rank` | the droplets ranked below this many are ambient |
| `--round-counts` / `--no-round-counts` | rounds (1.8 on) | `round` | round non-integer counts |
| `--fdr` | 0.001 | | largest FDR called a cell |

## Run it

```bash
sc-droplets call dropletutils-emptydrops raw_feature_bc_matrix.h5 -o calls.h5ad
sc-droplets call dropletutils-emptydrops raw_feature_bc_matrix.h5 --tool-version 1.26 --seed 42 --table calls.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "dropletutils-emptydrops", tool_version="1.26", random_state=42)
adata.obs["dropletutils_emptydrops_1_26_0_fdr"].describe()
```

Each droplet's class is `retained` (at or above `retain`), `rescued` (tested and called below
`retain`), `tested` (tested, not called), `ignored` (above `lower`, at or below `ignore`) or
`ambient`; tested droplets also get `pvalue`, `fdr`, `logprob` and `limited` columns.

## Differences from the original

None in the results, on x86-64. What it takes, and where R itself differs between platforms:

- The dispersion's likelihood is so flat near its optimum (a change of 0.03 in $$\alpha$$ moves it
  by about $$10^{-5}$$) that R's rounding steers `optimize`. sc-droplets evaluates it as R does:
  R's own `lgamma` (ported from R's `nmath`), the C library's `log`, and R's `sum`, which
  accumulates in x87 long double on x86-64 (emulated exactly where the hardware has no such
  type). R built for Apple silicon or aarch64 sums in plain double, so its $$\alpha$$, and the
  p-values that follow, can differ in the last digits from R on x86-64; sc-droplets gives the
  x86-64 result everywhere.
- `known.empty` (1.20 on) is not offered.
- 1.0 – 1.2 with a Dirichlet-multinomial (1.2's default) would need R's `rgamma` drawn one value
  at a time from R's single stream; it is refused. Pass `--alpha inf` for 1.2's multinomial.

## Validation

Compared with `DropletUtils::emptyDrops` from bioconda (1.0.3, 1.2.1, 1.4.2, 1.6.1, 1.8.0,
1.10.3, 1.12.1, 1.14.2, 1.18.0, 1.20.0, 1.26.0, 1.30.0, and 1.28 by its own knee in 1.30) with
`set.seed(100)` on 10x's pbmc_1k_v3 raw matrix (6.8 million barcodes) and on a simulated matrix:
every p-value, the estimated $$\alpha$$, every FDR and every call identical
(`tests/test_validation.py`, outputs in `tests/oracles/dropletutils/`, made by
`oracles/dropletutils/run.sh`).

The pieces are checked one by one too: R's seeding against `set.seed` and dqrng, PCG32 against
its reference outputs, the Monte Carlo against DropletUtils' compiled `montecarlo_pval` with
every simulated value pinned between two observed ones ($$\alpha$$ infinite, 0.5, 50 and 2000),
edgeR's Good-Turing, R's `lgamma`, `sum`, `optimize` and `p.adjust` against R.

## References

- Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based
  single-cell RNA sequencing data, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- Amezquita et al. 2020, Orchestrating single-cell analysis with Bioconductor, *Nat Methods*
  17:137, [doi:10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x).
- Gale & Sampson 1995, Good-Turing frequency estimation without tears, *J Quant Linguist* 2:217,
  [doi:10.1080/09296179508590051](https://doi.org/10.1080/09296179508590051).
- O'Neill 2014, PCG: a family of simple fast space-efficient statistically good algorithms for
  random number generation, HMC-CS-2014-0905.
- Source: [`R/emptyDrops.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/emptyDrops.R),
  [`R/ambientProfileEmpty.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/ambientProfileEmpty.R),
  [`src/montecarlo_pval.cpp`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/src/montecarlo_pval.cpp)
  (1.30); edgeR's `R/goodTuring.R` and `src/good_turing.c`.
