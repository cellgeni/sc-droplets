---
title: DIEM
parent: Planned
grand_parent: Methods
nav_order: 2
method_ids: [diem]
status: planned
---

# DIEM: a mixture of debris and cell profiles
{: .no_toc }

DIEM's semi-supervised filter, written for single-nucleus data. The barcodes below a count floor
fix the profile of debris; those above it are clustered, and a multinomial mixture is fitted by
EM with the debris component held to the low-count barcodes. Clusters whose droplets detect few
genes join the debris, the genes enriched there give every droplet a debris score, and droplets
with a low score and enough genes are cells. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `diem` (planned) |
|---|---|
| Tool | diem 1.0 – 2.4.1 (GitHub `marcalva/diem` only), the README's workflow ending in `call_targets` |
| Category | mixture |
| Counts | UMIs per gene: each droplet's whole profile, over every row given |
| Needs expected cells | no |
| Deterministic | yes: it calls `set.seed(1)` before the PCA and before k-means. With more than one thread (the default is one) an unguarded shared update makes runs differ |
| Original licence | GPL-3 |
| Extra inputs | none; the matrix, raw |

## How it works

Release 2.4.1, as the README runs it:

1. **Droplets.** Barcodes with no count are dropped. The test set is the barcodes with total
   $$\ge$$ `min_counts` (100) and at least `min_genes` (0) genes; `top_n`, when given, raises the
   floor to the total of the `top_n`-th barcode. Every other barcode with a count is the debris
   set. Genes are those with a count anywhere (`cpm_thresh = 0`).
2. **Variable genes.** On the test droplets with at least 200 genes (the cluster set): per gene,
   $$\log_{10}(\text{mean} + 1)$$ and $$\log_{10}(\text{variance} + 1)$$ of the raw counts
   (variance with divisor $$n - 1$$); R's `loess` of the second on the first, `span = 0.3`, R's
   defaults otherwise; the 2,000 genes with the largest residuals.
3. **PCs.** Each droplet's counts divided by its total, times the median total, then `log1p`.
   30 PCs of the variable genes, centred and scaled, by `irlba::prcomp_irlba` after `set.seed(1)`
   (base `prcomp` when 30 is at least half the smaller dimension).
4. **Start.** Each PC scaled to unit variance; `set.seed(1)`; R's `kmeans` with 20 centres,
   Hartigan–Wong, at most 10 iterations, one start (the documented `nstart_init = 30` is not
   passed on). The debris set is component 1 and the k-means clusters are 2 to 21; test droplets
   with fewer than 200 genes start unassigned. From these hard labels $$Z$$,
   $$\alpha_{gk} \propto 2\cdot10^{-10} + \sum_i Z_{ik}x_{gi}$$ and $$\pi_k \propto \sum_i Z_{ik}$$:
   the pseudocount is passed by position into the prior's argument, so each adds $$10^{-10}$$.
5. **EM.** A multinomial mixture over every droplet with a count. Debris-set droplets stay fixed
   in component 1; test droplets are free. One E-step from the start, then repeat: M-step,
   $$\pi$$ from every row of $$Z$$ (the fixed ones too) and
   $$\alpha_{gk} \propto 10^{-10} + \sum_i Z_{ik}x_{gi}$$; E-step,
   $$Z_{ik} = \operatorname{softmax}_k\big(\sum_g x_{gi}\log\alpha_{gk} + \log\pi_k\big)$$. It stops
   when the mean over free droplets of $$\sum_k \lvert Z^{\text{new}}_{ik} - Z^{\text{old}}_{ik}\rvert$$
   is below $$10^{-4}$$, or after 1,000 iterations.
6. **Clusters.** Each test droplet goes to its most probable component, the first on ties.
7. **Debris clusters.** Component 1, plus the clusters whose test droplets detect fewer than 200
   genes on average. The second part has a bug: `which(tapply(genes, Cluster, mean) < 200)` gives
   positions in the alphabetical order of the labels ("1", "10", "11", …, "19", "2", "20", …),
   and those positions are then taken as labels. With 10 or more labels (21 by default) the wrong
   clusters can be flagged: a low-gene cluster "10", at position 2, flags cluster 2.
8. **Debris score.** Expression $$\log(1 + x_{gi}/n_i)$$ over the expressed genes. A Welch t-test
   per gene, test droplets in debris clusters against the rest (genes with mean above 0 in both),
   Benjamini–Hochberg; debris genes have adjusted $$p < 0.05$$ and $$\log_2$$ fold change above 0,
   at most the first 5,000 by mean difference. A droplet's sum over the debris genes, minus the
   smallest of the clusters' mean sums, divided by the mean of that over the droplets in debris
   clusters, is its score: the lowest cluster averages 0 and the debris droplets 1.
9. **Call.** `call_targets(thresh_score = 0.5)`: a test droplet is `Clean`, a cell, when its score
   is below 0.5 and it detects at least 200 genes; otherwise `Debris`. For single cells the README
   calls `call_targets(clusters = "debris", thresh_score = NULL)` instead: clean when outside
   the debris clusters (the same flagged set) with at least 200 genes. The debris set gets no call.

DIEM seeds itself, so a run repeats on the same R, irlba and BLAS. The PCs depend on irlba's
random start vector and on the BLAS, their signs included; k-means with one start depends on the
PCs; and the cluster labels feed the bug of step 7.

## Versions

| diem | Change to the call |
|---|---|
| 1.0 (`v1.0`, Nov 2019) | the preprint's method: a multinomial mixture; test set the top 10,000 barcodes with at least 100 counts and 100 genes; EM to $$10^{-8}$$; a cell has posterior $$\ge$$ 0.95 of a clean component and at least 200 genes |
| 2.0.0 – 2.2.0 (Feb 2020) | Dirichlet-multinomial mixtures; 30 k-means clusters with 30 starts; call on posterior $$\ge$$ 0.95 (`pp_thresh`) |
| 2.3.0 (untagged, Apr – Jun 2020) | the debris score and `call_targets(thresh_score = 0.5)`; multinomial by default; k-means with 20 centres and 30 starts, clusters under 10 droplets merged; at most 300 EM iterations; the bug of step 7. The published method (which commit is not known) |
| 2.4.0 (Feb 2022) | EM rewritten in C++; k-means with one start, no merging; at most 1,000 iterations |
| 2.4.1 (Dec 2022) | fixes an out-of-bounds read in the split of free and fixed droplets, and lost dimnames that make 2.4.0's debris score fail with recent Matrix |

2.4.1 is the repository's last commit, so a pipeline that installs from GitHub today runs it.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `min_counts` (`set_debris_test_set`) | 100 | test droplets at or above it; the rest make the debris set |
| `min_genes` (`set_debris_test_set`) | 0 | genes a test droplet needs |
| `top_n` (`set_debris_test_set`) | none | raise the floor to the `top_n`-th barcode's total |
| `cpm_thresh` (`filter_genes`) | 0 | genes kept: CPM above it |
| `min_genes` (`get_pcs`) | 200 | droplets clustered at the start |
| `n_var_genes`, `lss`, `n_pcs` (`get_pcs`) | 2,000, 0.3, 30 | variable genes, the loess span, PCs |
| `k_init` (`init`) | 20 | k-means clusters |
| `seedn` (`get_pcs`, `init`) | 1 | the seed before the PCA and k-means |
| `model` (`init`, `run_em`) | `"mltn"` | multinomial; `"DM"` fits a Dirichlet-multinomial |
| `max_iter`, `eps`, `psc` (`run_em`) | 1,000, $$10^{-4}$$, $$10^{-10}$$ | EM's stop and pseudocount |
| `threads` (`run_em`) | 1 | worker threads; more race |
| `thresh_genes` (`estimate_dbr_score`) | 200 | clusters below it on average join the debris |
| `thresh_p`, `thresh_logfc`, `max_genes` (`estimate_dbr_score`) | 0.05, 0, 5,000 | debris genes |
| `thresh_score` (`call_targets`) | 0.5 | cells score below it |
| `min_genes` (`call_targets`) | 200 | cells detect at least this many genes |

## Status

Planned for 0.2. Implementing it needs:

- R's `loess` at its defaults, for the variable genes;
- irlba's `prcomp_irlba` with its random start vector, under R's Mersenne-Twister (the package
  has the generator);
- R's `kmeans` (Hartigan–Wong, in Fortran) with its start from `sample.int`;
- the multinomial EM, Welch t-tests and Benjamini–Hochberg (the package has R's `p.adjust`);
- the bug of step 7, reproduced as it is.

The oracle is diem built from GitHub at `v2.4.1` (a C++ compiler with Eigen; conda-forge
`r-rcppeigen`, `r-rcppprogress`, `r-irlba`), with `threads = 1`. It is deterministic, but the PCs
and the one-start k-means make an exact match from the matrix alone fragile. Two routes: compare
the calls within a tolerance, or stage the oracle, taking the original's k-means labels as a
fixture and checking the EM, the score and the call exactly.

## References

- Alvarez, Rahmani, Jew et al. 2020, Enhancing droplet-based single-nucleus RNA-seq resolution
  using the semi-supervised machine learning classifier DIEM, *Scientific Reports* 10:11019,
  [doi:10.1038/s41598-020-67513-5](https://doi.org/10.1038/s41598-020-67513-5).
- Source: [marcalva/diem at v2.4.1](https://github.com/marcalva/diem/tree/v2.4.1) (`R/filter.R`,
  `R/var.R`, `R/normalize.R`, `R/init_clust.R`, `R/em.R`, `R/post.R`, `src/MultEM.cpp`).
- Research notes: `docs/_research/planned-rna.md` §2.
