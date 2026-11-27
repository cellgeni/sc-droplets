---
title: CB2
parent: Planned
grand_parent: Methods
nav_order: 1
method_ids: [cb2]
status: planned
---

# CB2: clusters of candidates, then emptyDrops
{: .no_toc }

scCB2's `CB2FindCell`. Candidates between a lower bound and a knee are grouped into tight
clusters; a cluster whose counts correlate less with the ambient profile than ambient draws do is
a cell, all its members at once. The barcodes left over go to DropletUtils' emptyDrops. Testing
clusters gives power to small cells that emptyDrops tests one by one. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `cb2` (planned) |
|---|---|
| Tool | scCB2 1.0.0 – 1.24.0 (Bioconductor), `CB2FindCell`; the GitHub releases 0.99.x before |
| Category | statistical-test |
| Counts | UMIs per barcode |
| Needs expected cells | no |
| Deterministic | no: the cluster test's simulations run on worker processes that are never seeded |
| Original licence | GPL-3 |
| Extra inputs | none; the matrix, raw |

## How it works

1. **Features and barcodes.** Rows named like antibody or multiplexing features (`TotalSeqB`,
   `Cell Multiplexing Oligo`) are dropped; so are barcodes and genes with no count.
2. **Upper bound.** A modified [barcodeRanks](../dropletutils/barcoderanks.md): one point per run of
   equal totals at its mid-rank, `smooth.spline(df = 20)` of $$\log_{10}$$ total on $$\log_{10}$$
   rank fitted to *every* run (DropletUtils fits only those above `lower`). The knee is the least
   curvature between the steepest slope and the largest slope before it, and the inflection the
   steepest slope, both as fitted totals rounded half to even. With
   $$s = (\text{knee} - \text{inflection}) / 10$$ from the first fit, the fit is repeated with
   `lower` raised to the latest inflection plus $$s$$, until the knee stops moving or at most 3
   barcodes remain above the new bound. The last knee is `upper`.
3. **Classes.** Barcodes with total $$\ge$$ `upper` are retained as cells. Those $$\le$$ `lower`
   (100) are background; their summed counts give the ambient profile $$p$$ (edgeR's Good-Turing
   proportions). Those in between are candidates.
4. **Barnyard.** When every gene name holds a `_` (a species prefix), a candidate is tested only
   if one species holds at least 75 % of its UMIs.
5. **Genes.** The 30,000 genes with the largest $$p$$ (ties ranked together).
6. **Statistic.** Each candidate's Pearson correlation with $$p$$.
7. **Clusters.** Candidates sorted by total, in chunks of 1,000 (a last chunk under 500 joins the
   one before). In a chunk of $$m > 20$$ barcodes, complete-linkage hierarchical clustering on
   $$1 - $$ correlation, cut into $$\lfloor m/20 \rfloor$$ clusters. A cluster's statistic is its
   members' median correlation with $$p$$; its size, the total of the member nearest that median;
   its tightness, the mean correlation between its members.
8. **Tight clusters only.** A profile $$c$$: the retained cells' Good-Turing profile when there
   are at least 10 of them and it has at least the entropy of $$p$$, otherwise $$p$$. For sizes
   $$100b$$, $$b = 1, \dots, 10$$, the mean correlation between 100 multinomial draws of $$100b$$
   UMIs from $$c$$ sets a bar. A cluster tighter than the bar for its size is tested (sizes below
   300 use the bar for 200, sizes of 1,100 and more the bar for 900); the others return as single
   barcodes.
9. **Cluster test.** For each size, 1,000 multinomial draws from $$p$$ and their correlations with
   $$p$$; $$P = (\#\{\text{draws} \le \text{observed}\} + 1)/1001$$, since a low correlation is
   cell-like. Benjamini-Hochberg over the clusters; clusters with adjusted $$P \le$$
   `FDR_threshold` (0.01) are cells.
10. **Single barcodes.** [emptyDrops](../dropletutils/emptydrops.md) with its defaults on the
    remaining candidates, the background and the retained cells (`lower` and `retain = upper`
    passed on); FDR $$\le$$ `FDR_threshold` makes a cell.
11. The cells are emptyDrops' cells, the retained barcodes among them, and the cluster cells.
    `GetCellMat` then drops cells with more than `MTfilter` of their UMIs in mitochondrial genes.

Step 8's draws and step 10 use R's random numbers, reproducible with `set.seed`. Step 9 runs on a
cluster of worker processes started inside the function, whose random numbers are not seeded:
its p-values change from run to run, whatever the caller seeds.

## Versions

| scCB2 | Change to the call |
|---|---|
| 0.99.x (GitHub, 2019 – 2020) | the paper's code; 0.99.29 changed the sparse correlation, 0.99.30 added `TopNGene`, 0.99.31 the barnyard filter |
| 1.0.0 (Bioconductor 3.12) | first release; `CB2FindCell` is the same from here to 1.24.0 |
| 1.6.0 (Bioconductor 3.15) | `GetCellMat`'s `MTfilter` default 1 → 0.25: cells with more than 25 % mitochondrial UMIs dropped; `CheckBackgroundCutoff` added (advice only) |

The call also depends on the installed DropletUtils, whose emptyDrops decides every barcode
outside the clusters: a CB2 call is a pair (scCB2 release, DropletUtils release).

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `FDR_threshold` | 0.01 | the FDR cut, for clusters and single barcodes |
| `lower` | 100 | background at or below it; candidates above |
| `upper` | the knee (step 2) | retained as cells at or above it |
| `GeneExpressionOnly` | `TRUE` | drop antibody and multiplexing rows by name |
| `TopNGene` | 30,000 | genes kept for the correlations |
| `Ncores` | 2 | worker processes for the cluster test |
| `MTfilter` (`GetCellMat`) | 0.25 (1 before 1.6.0) | the largest mitochondrial fraction a cell may have |

## Status

Planned for 0.2. Implementing it needs:

- the knee of step 2, a variant of barcodeRanks' spline curvature on every run; the package
  already ports R's `smooth.spline`;
- R's `hclust` with complete linkage, including its tie handling, and `cutree`;
- R's `rmultinom` (binomial draws one gene at a time), because step 8's draws advance the stream
  emptyDrops then takes its seeds from;
- [emptyDrops](../dropletutils/emptydrops.md), already implemented, at the DropletUtils release
  the pipeline ran.

The oracle is scCB2 1.24.0 from Bioconductor with a pinned DropletUtils, `set.seed` as for the
DropletUtils oracle. As shipped, the cluster test cannot be reproduced exactly, so the comparison
is within the Monte Carlo noise band, unless the oracle seeds the workers itself (for instance
`Ncores = 1` with a worker profile that calls `set.seed`).

## References

- Ni, Chen, Brown & Kendziorski 2020, CB2 improves power of cell detection in droplet-based
  single-cell RNA sequencing data, *Genome Biology* 21:137,
  [doi:10.1186/s13059-020-02054-8](https://doi.org/10.1186/s13059-020-02054-8).
- Source: [scCB2 on Bioconductor](https://bioconductor.org/packages/scCB2/)
  ([git](https://git.bioconductor.org/packages/scCB2), `R/CB2FindCell.R`, `R/GetCellMat.R`;
  mirror [zijianni/scCB2](https://github.com/zijianni/scCB2)).
- Research notes: `docs/_research/planned-rna.md` §1.
