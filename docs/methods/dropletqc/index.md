---
title: DropletQC
parent: Methods
nav_order: 15
has_children: true
---

# DropletQC
{: .no_toc }

DropletQC (Muskovic & Powell 2021) flags empty droplets and damaged cells among the barcodes
another caller kept, from each barcode's nuclear fraction: the share of its reads in introns. An
empty droplet holds mostly spliced, cytoplasmic RNA and has a low fraction; a nucleus without its
cytoplasm has a high one.

1. TOC
{:toc}

## Which release does what

| DropletQC | What changed |
|---|---|
| before 1.0.0 (as `dropletQC`, Feb 2021) | `nuclear_fraction_tags` returned the fractions in alphabetical order (fixed in 57f6f44); renamed DropletQC in May 2021 |
| 1.0.0 (`v1.0.0`, Nov 2021) | the only tag, and the last commit |

The call also depends on the installed ks, whose binned estimator changed in 1.14.3: see the
method's *Versions*.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [empty droplets](empty-drops.md) | `dropletqc` | the valley after the first peak of the nuclear fraction's density; a rescue for barcodes with many UMIs |

## Not implemented

- The nuclear fraction itself, which DropletQC reads from a BAM (Cell Ranger's `RE` tags, or a
  GTF): a count matrix does not hold it. The method takes it from an `obs` column or a table, or
  uses the unspliced fraction of the package's splicing readers as an approximation.
- `identify_damaged_cells`, which needs mclust's EEI mixture with its hierarchical start and a
  cell type per barcode: a follow-up.
