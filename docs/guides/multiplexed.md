---
title: Multiplexed and Flex runs
parent: Guides
nav_order: 6
---

# Multiplexed and Flex runs
{: .no_toc }

Cell Ranger calls cells of a multiplexed run sample by sample, then removes the barcodes of
crowded GEMs across all samples. To reproduce that, call each sample's raw matrix on its own and
run the high-occupancy filter on the union.

1. TOC
{:toc}

## Which runs this is about

| Run | What Cell Ranger calls on |
|---|---|
| Flex (Fixed RNA Profiling), several samples | each sample's barcodes: the GEM barcode joined with that sample's probe barcode |
| On-chip multiplexing (OCM, GEM-X) | each sample's barcodes, by overhang tag |
| CellPlex or hashtags | the whole run; the tags assign samples to the cells afterwards |
| several GEM groups in one `count` run | each GEM group (the `-1`, `-2` suffix) |
| a barnyard (two genomes) | each genome on its own features; the cells are the union |

For hashing and genetic demultiplexing, call cells on the whole run, as for any sample. The rest
of this guide is about Flex and OCM.

## Call each sample

`cellranger multi` writes a raw matrix per sample:
`outs/per_sample_outs/<sample>/count/sample_raw_feature_bc_matrix.h5`. Call each one with the
release that ran:

```bash
sc-droplets call cellranger per_sample_outs/LymphNode1/count/sample_raw_feature_bc_matrix.h5 --tool-version 7.2.0 --barcodes LymphNode1.tsv.gz
```

The chemistry is read from the file (`Fixed RNA Profiling`), and it sets the ambient range and the
FDR for 7.2 and later. Give `--chemistry` when the input has no such attribute, such as an
`.h5ad` or STARsolo output. When a sample uses more than one probe barcode or OCM tag, give their
number with `--n-probe-barcodes`: the ambient range grows with it.

## Remove crowded GEMs

Cell Ranger then joins the samples' cells. A GEM, the first 16 bases of the barcode, that holds
cells of more probe barcodes than its occupancy allows probably held several cells, and all its
cells are removed. The threshold comes from a simulation of GEM occupancy: Poisson up to 9.x,
negative binomial from 10.0. The command line does not run this step yet. In Python:

```python
import sc_droplets as scd
from sc_droplets.methods.cellranger import high_occupancy_gems

samples = ["LymphNode1", "LymphNode2", "Spleen1", "Spleen2"]
calls = {}
for name in samples:
    adata = scd.read(f"per_sample_outs/{name}/count/sample_raw_feature_bc_matrix.h5")
    result = scd.call(adata, "cellranger@7.2.0", inplace=False)
    calls[name] = set(result.barcodes[result.is_cell])
removed, info = high_occupancy_gems(set().union(*calls.values()), "7.2.0")
cells = {name: called - set(removed) for name, called in calls.items()}
```

On 10x's 4-plex lymph node and spleen run (Cell Ranger 7.2), this gives 10x's cells for all four
samples: the filter's threshold is 3 probe barcodes per GEM, and it removes 28 barcodes. Called on
its own, the first lymph node sample has 7 more cells than 10x's, exactly the barcodes of crowded
GEMs.

## Other methods on Flex samples

Any method runs on a sample's raw matrix. Methods written for single-sample 3′ data know nothing
of the probe barcodes or the crowded GEMs, so compare their calls with Cell Ranger's before the
filter. In the [comparison](../overview/comparison.md), STARsolo's EmptyDrops_CR, with Cell
Ranger 3's ambient ranks, calls 8,369 cells on the lymph node sample, where Cell Ranger 7.2 calls
7,586. Splitting a whole multiplexed raw matrix into samples by probe barcode is not offered yet.
