---
title: valiDrops
parent: Planned
grand_parent: Methods
nav_order: 5
method_ids: [validrops]
status: planned
---

# valiDrops: a segmented knee, then quality and cluster filters
{: .no_toc }

valiDrops's `valiDrops`, for single cells and nuclei, in three stages. A segmented regression
of the log-log rank curve gives a count threshold; adaptive filters on the mitochondrial
fraction, genes against UMIs and the protein-coding fraction remove poor barcodes; then the
survivors are clustered, and only clusters with clear marker genes keep their barcodes. Planned
for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `validrops` (planned) |
|---|---|
| Tool | valiDrops 1.0.0 (GitHub `madsen-lab/valiDrops` only, tag `v1.0.0`; DESCRIPTION 0.1.0) |
| Category | knee |
| Counts | UMIs per barcode in stage 1; genes detected and UMIs per gene after |
| Needs expected cells | no |
| Deterministic | no: `normalmixEM`'s random starts, the subsample of the mitochondrial fallback and irlba's start vector, and the package sets no seed. Stage 1 alone is deterministic with segmented 2.0-3 or later |
| Original licence | GPL-3 |
| Extra inputs | gene names its bundled tables recognise (human, mouse, rat, zebrafish, worm, fly), to find mitochondrial, ribosomal and protein-coding genes |

## How it works

Species and annotation (symbol, Entrez, Ensembl, HGNC, MGI; version suffixes stripped) are
detected against bundled tables: mitochondrial genes are those on chromosome MT, ribosomal genes
match `^RPL` or `^RPS`.

**Stage 1, the threshold** (`rank_barcodes`):

1. Barcodes with a count; `rank(-total)` (ties averaged); one point per distinct total, at its
   average rank. $$x = \ln$$ rank and $$y = \ln$$ total, each smoothed by a centred moving average
   of width $$\lceil 2N^{1/3}\rceil$$, $$N$$ the number of distinct totals (`zoo::rollmean`).
2. For $$k = 2, \dots, 5$$ breakpoints: `segmented` on `lm(y ~ x)`, starting from $$k$$ values
   evenly spaced between the $$\alpha$$ and $$1 - \alpha$$ quantiles of $$x$$, with segmented's
   own `alpha` $$= 0.999\,\alpha$$ and 10 bootstrap restarts; $$\alpha = 0.001$$, raised by 0.001
   up to 0.05 while segmented fails with "psi values too close".
3. The model: the fewest breakpoints whose RMSE is at most 1.5 times the smallest.
4. Between consecutive segments with slopes $$s_i, s_{i+1}$$, the angle
   $$\theta_i = \arctan\frac{s_i - s_{i+1}}{1 + s_i s_{i+1}}$$ in degrees. The breakpoint is the
   one with the smallest angle, from the second on (the first is never chosen). The threshold is
   the distinct total whose $$\ln$$ rank lies nearest that breakpoint, returned as
   `exp(log(total))`; barcodes with total at or above it go on. As `exp(log(t))` need not equal
   $$t$$ in floating point, barcodes at exactly that total can fall either side: a port must
   compute the same `exp(log)`.

**Stage 2, quality** (`quality_metrics`, `quality_filter`). Per barcode: $$\ln$$ UMIs, $$\ln$$
genes, and the mitochondrial, ribosomal and protein-coding UMI fractions. Three filters, each on
the survivors of the one before:

1. **Mitochondrial.** 10 times: `mixtools::normalmixEM` of $$\ln$$ genes, 2 components, random
   start (`rexp`, `rnorm`); the barcodes of the higher-mean component; over $$t$$ from their
   median mitochondrial fraction to 1 in steps of 0.001, $$c(t)$$ the number at or below $$t$$;
   the threshold `inflection::uik(t, c)`. The median of the 10. If it exceeds 0.3: 10 times on a
   random subsample of $$\min(5000, 0.8n)$$ barcodes, `segmented(lm(ln genes ~ mito))` with 1 to
   5 breakpoints until the smallest is at most 0.3; the median. Barcodes at or below it stay.
2. **Distance.** `segmented(lm(ln genes ~ ln UMIs))` with 3 breakpoints (fewer if it fails);
   residuals within the median $$\pm\,5$$ times `robustbase::Sn` stay.
3. **Coding.** A protein-coding fraction within the median $$\pm\,3\,S_n$$ stays.

**Stage 3, clusters** (`expression_metrics`, `expression_filter`), on the protein-coding genes of
stage 2's barcodes: normalised to 10,000 and `log1p`; 5,000 genes by scry's binomial deviance;
z-scores; 10 PCs by irlba; Seurat's shared-nearest-neighbour graph and Louvain clusters at
resolution 0.1 (shallow) and at a deep resolution where the smallest cluster has about 5
barcodes (a grid 1 to 20, then $$\pm 0.9$$ in steps of 0.1). Each deep cluster is compared with
the barcodes outside its majority shallow cluster: genes detected in 10 % of either group with
$$\log_2$$ fold change at least 0.25, Wilcoxon (presto if installed, else `wilcox.test`),
Bonferroni over all expressed genes. Its top markers are at most 10 with FDR $$\le$$ 0.05, by
$$p$$ (the first gene alone when none is significant). With $$p_1, p_2$$ the shares of
barcodes detecting a marker inside and outside, a cluster is valid when:

- no marker has $$(p_1 - p_2)/p_1 < -0.01$$, and their mean $$(p_1 - p_2)/p_1$$ is at least 0.2;
- mean $$p_1$$ is at least 0.3 and mean $$p_2$$ at most 0.7;
- at least one marker is significant, and $$-\log_{10}$$ of the smallest FDR reaches a level from
  the data: the smaller of its median $$+\,3\,S_n$$ over the clusters with mean
  $$(p_1 - p_2)/p_1 \le 0.4$$, and the breakpoint of a one-breakpoint `segmented` fit of that
  mean on it;
- the share of genes differentially expressed exceeds 0.01;
- its median mitochondrial and ribosomal fractions are at most median $$+\,3\,S_n$$ across clusters.

The cells are the barcodes of valid clusters; `stageThree = FALSE` stops after stage 2.
`label_dead` (off by default) adds a dead, live or uncertain label from a ridge `cv.glmnet`, which
does not change the call. Seurat's `FindClusters` seeds itself with 0, but nothing else is seeded:
only a `set.seed` by the caller makes a run repeat.

## Versions

| valiDrops | Change to the call |
|---|---|
| before the tag (Nov 2022 – Sep 2023) | the three stages settle in Nov 2022; 0834284 (Sep 2023) makes stage 3 optional; 905cd3a and the next two commits accept a `DelayedMatrix` |
| 1.0.0 (`v1.0.0`, Sep 2023) | the only tag; `main` adds only a LICENSE file (Mar 2025) |

The call also depends on the installed segmented. Up to 2.0-2 its bootstrap restarts draw from
the caller's random numbers; from 2.0-3 (Feb 2024), with `seed = NULL`, it seeds itself from the
digits of the mean of $$y$$ and restores the caller's state, which makes stage 1 deterministic.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `rank_barcodes` | `TRUE` | run stage 1 |
| `psi.min`, `psi.max` | 2, 5 | the breakpoints tried |
| `alpha`, `alpha.max` | 0.001, 0.05 | where breakpoints may lie, and how far that may widen |
| `boot` | 10 | segmented's bootstrap restarts |
| `factor` | 1.5 | a model's RMSE may exceed the best by this factor |
| `mito.max`, `mito.nreps` | 0.3, 10 | the most a mitochondrial threshold may be before the fallback; repeats |
| `npsi` | 3 | breakpoints of the distance fit |
| `dist.threshold`, `coding.threshold` | 5, 3 | $$S_n$$ multiples of the distance and coding filters |
| `stageThree` | `TRUE` | run stage 3 |
| `nfeats`, `npcs`, `res.shallow`, `k.min`, `top.n` | 5,000, 10, 0.1, 5, 10 | stage 3's genes, PCs, shallow resolution, smallest deep cluster, markers |
| `min.target.pct`, `max.background.pct`, `min.diff.pct`, `min.de.frac` | 0.3, 0.7, 0.2, 0.01 | the marker rules |
| `mitochondrial_clusters`, `ribosomal_clusters` | 3, 3 | $$S_n$$ multiples for the clusters' mitochondrial and ribosomal medians |
| `species`, `annotation` | detected | the gene tables used |

## Status

Planned for 0.2. Implementing it needs:

- stage 1: segmented's fit, ported exactly, with its bootstrap restarts seeded from the decimal
  digits of the mean of $$y$$; zoo's moving average; R's averaged ranks;
- stage 2: mixtools' `normalmixEM` with R's random start, inflection's `uik` (already ported for
  [zUMIs](../zumis/cell-barcodes.md)), robustbase's $$S_n$$, and the bundled gene tables;
- stage 3: scry's deviance, irlba, Seurat's graph and Louvain clustering, and presto's Wilcoxon.

The oracle is valiDrops from GitHub at `v1.0.0` in R 4.0.2 or later, with Seurat, scry
(Bioconductor), presto (GitHub `immunogenomics/presto`), segmented 2.0-3 or later, inflection,
mixtools and robustbase, and `set.seed` before the call. Stage 1 is deterministic and a good exact
target. Stages 2 and 3 can only be compared within a tolerance, or stage by stage given the
original's intermediate outputs, such as its clusters.

## References

- Kavaliauskaite & Madsen 2023, Automatic quality control of single-cell and single-nucleus
  RNA-seq using valiDrops, *NAR Genomics and Bioinformatics* 5(4):lqad101,
  [doi:10.1093/nargab/lqad101](https://doi.org/10.1093/nargab/lqad101).
- Source: [madsen-lab/valiDrops at v1.0.0](https://github.com/madsen-lab/valiDrops/tree/v1.0.0)
  (`R/valiDrops.R`, `R/rank_barcodes.R`, `R/quality_metrics.R`, `R/quality_filter.R`,
  `R/expression_metrics.R`, `R/expression_filter.R`).
- Research notes: `docs/_research/planned-rna.md` §5.
