---
title: DropletQC
parent: Planned
grand_parent: Methods
nav_order: 3
method_ids: [dropletqc]
status: planned
---

# DropletQC: the first valley of the nuclear fraction
{: .no_toc }

DropletQC's `identify_empty_drops`. Each barcode's nuclear fraction is the share of its reads
that fall in introns, read from the BAM. An empty droplet holds mostly spliced, cytoplasmic RNA,
so its fraction is low: DropletQC cuts the fraction at the first valley of a kernel density
estimate, then rescues barcodes with many UMIs. An optional second step flags damaged cells within
each cell type. DropletQC flags barcodes another caller has kept; on a whole raw matrix the rule
is untested. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `dropletqc` (planned) |
|---|---|
| Tool | DropletQC 1.0.0 (GitHub `powellgenomicslab/DropletQC` only, tag `v1.0.0`) |
| Category | density |
| Counts | the nuclear fraction (intronic over intronic and exonic reads) and UMIs per barcode |
| Needs expected cells | no |
| Deterministic | yes: no random numbers |
| Original licence | MIT |
| Extra inputs | a BAM with Cell Ranger's `CB` and `RE` tags (or a BAM and a GTF/GFF); the barcodes to flag; a cell type per barcode for the damaged cells |

## How it works

### Nuclear fraction

1. **Reading.** `nuclear_fraction_tags` cuts the genome (the sequences with mapped reads, from
   `idxstatsBam`) into `tiles` (100) ranges. For each, `scanBam` reads the `CB` and `RE` tags of
   every alignment overlapping it, with no filter: PCR duplicates, secondary alignments and
   multimappers all count, and an alignment across a tile boundary counts in both tiles. Records
   without `CB` or `RE`, or whose `CB` is not among the barcodes given, are dropped.
2. **Fraction.** Per barcode, $$\text{NF} = N/(N + E)$$, with $$N$$ the records tagged `RE = N`
   (intronic) and $$E$$ those tagged `RE = E` (exonic); intergenic records (`RE = I`) are
   ignored, and 0/0 gives 0. These are alignment records, not UMIs.
3. **Without `RE`.** `nuclear_fraction_annotation` takes exons (merged, strand ignored) and
   introns minus exons from a GTF/GFF, counts overlaps of at least 5 bases, and counts a read that
   overlaps both in both; it reads the BAM in 1,000 tiles.

The `outs` shortcut reads the barcodes from Cell Ranger's `filtered_feature_bc_matrix`; the
paper applies DropletQC after EmptyDrops (`lower` 500, FDR 1 %), Cell Ranger, CellBender or EmptyNN.

### Empty droplets

`identify_empty_drops` takes the NF and the UMIs of each barcode:

1. **Derivative.** `ks::kdde(nf, deriv.order = 1)`: a binned kernel estimate of the first
   derivative of the density of NF. Gaussian kernel; bandwidth $$h$$ from `ks::hpi(nf, nstage = 2,
   binned = TRUE, deriv.order = 1)`, the two-stage plug-in selector for the derivative; 401 grid
   points from $$\min - 4.7h$$ to $$\max + 4.7h$$ (ks 1.15.3's binning). A density estimate with
   its own bandwidth is computed for the plot only.
2. **Cut-off.** Whether the derivative is positive, along the grid, run-length encoded. The
   cut-off is the grid point at index (length of run 1 + length of run 2), 1-based: the last
   point of the second run. When the derivative starts positive, that is the first minimum after
   the first peak; when it starts at zero or below, it is the first peak, and nothing checks for
   it. Fewer than 4 runs gives a warning, nothing more.
3. **Call.** `empty_droplet` when NF is below the cut-off, `cell` otherwise. Then `cell` again
   when NF > `nf_rescue` (0.05) and UMIs > `umi_rescue` (1,000), both strict.

### Damaged cells

`identify_damaged_cells` (optional) takes, besides, a cell type per barcode. For the barcodes
called `cell`, per cell type: `mclust::Mclust(cbind(NF, log10 UMIs), G = 1:2, modelNames = "EEI")`,
the number of components chosen by BIC. A type has damaged cells only when there are two
components, the one with the higher NF mean has the lower UMI mean, the NF means differ by more
than `nf_sep` (0.15), and $$10^{\mu_{\text{UMI,damaged}}} < 10^{\mu_{\text{UMI,cell}}}(1 - 0.5)$$,
with 0.5 from `umi_sep_perc` (50). The members of the high-NF component (mclust's MAP
classification) become `damaged_cell`.

Nothing draws random numbers: ks's plug-in bandwidth and binned estimate, and mclust's
hierarchical start, are deterministic. The cut-off depends on the ks version to the grid point.

## Versions

| DropletQC | Change to the call |
|---|---|
| before 1.0.0 (as `dropletQC`, Feb 2021) | 57f6f44 returns the fractions in the input's order (they came back in alphabetical order); 2a9ebb6 adds the plot option; renamed DropletQC in May 2021 |
| 1.0.0 (`v1.0.0`, Nov 2021) | the only tag, and the last commit; DESCRIPTION still says 0.0.0.9000 |

One variant.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `tiles` (`nuclear_fraction_tags`) | 100 | ranges the BAM is read in; an alignment across a boundary counts twice |
| `cell_barcode_tag` | `"CB"` | the barcode tag |
| `region_type_tag`, `exon_tag`, `intron_tag` | `"RE"`, `"E"`, `"N"` | the region tag and its exonic and intronic values |
| `nf_rescue` (`identify_empty_drops`) | 0.05 | a rescued cell has a larger NF |
| `umi_rescue` (`identify_empty_drops`) | 1,000 | and more UMIs |
| `nf_sep` (`identify_damaged_cells`) | 0.15 | the least gap between the two NF means |
| `umi_sep_perc` (`identify_damaged_cells`) | 50 | how much lower, in percent, the damaged cells' mean UMIs must be |

## Status

Planned for 0.2. Implementing it needs:

- the nuclear fraction from reads. The package's splicing readers give an unspliced fraction from
  UMIs ([Inputs](../../reference/inputs.md#splicing)), counted by velocyto, STARsolo, kb or
  alevin-fry rules: a different quantity, usable as an approximation, not as the original. The
  original needs a BAM with `RE` tags, or per-barcode intronic and exonic read counts;
- ks's `hpi` (two-stage plug-in for the derivative, binned) and binned `kdde`, to the grid point,
  at a pinned ks (1.15.3 is current; older releases are still to be checked for changes to
  `hpi`). The package's [kernel-density valley](../generic/kde-valley.md) chooses its valley
  differently;
- for damaged cells, mclust's EEI mixture with BIC, and a cell type for every barcode. The
  package's [two-dimensional mixture](../generic/gmm-2d.md) is scikit-learn's, with full
  covariances.

The oracle is DropletQC from GitHub at `v1.0.0` in R with Bioconductor's Rsamtools,
GenomicAlignments, GenomicFeatures and rtracklayer and CRAN's ks and mclust. It splits in two.
The rule needs only (NF, UMIs) per barcode, so `identify_empty_drops` can be fed a table we
compute and checked exactly. The fraction from a BAM needs a BAM: none of the repository's public
datasets has one, so a small public 10x BAM is needed. What to apply it to is also open: on a raw
matrix the density is dominated by empty droplets.

## References

- Muskovic & Powell 2021, DropletQC: improved identification of empty droplets and damaged cells
  in single-cell RNA-seq data, *Genome Biology* 22:329,
  [doi:10.1186/s13059-021-02547-0](https://doi.org/10.1186/s13059-021-02547-0).
- Source: [powellgenomicslab/DropletQC at v1.0.0](https://github.com/powellgenomicslab/DropletQC/tree/v1.0.0)
  (`R/nuclear_fraction_tags.R`, `R/nuclear_fraction_annotation.R`, `R/identify_empty_drops.R`,
  `R/identify_damaged_cells.R`).
- Research notes: `docs/_research/planned-rna.md` §3.
