---
title: Read-based methods
parent: Guides
nav_order: 5
---

# Read-based methods
{: .no_toc }

Several tools rank barcodes by reads, not UMIs: UMI-tools' whitelist, salmon alevin, alevin-fry,
bustools, zUMIs and Drop-seq. A count matrix holds UMIs, so the package counts UMIs unless it is
given reads. To reproduce these tools, give it the reads.

1. TOC
{:toc}

## Why it matters

The number of reads per UMI varies from barcode to barcode with sequencing depth, duplication and
the share of reads that map. Reads therefore rank barcodes differently from UMIs, and draw a
different curve. A knee found on reads is not the knee on UMIs. Run on UMIs, a read-based method
is a reasonable rule of its own, but not the original's call.

## Give the reads

`--reads` takes either an `obs` column of an `.h5ad`, or a table with two columns, barcode and
read count. The table is tab-separated if its name contains `.tsv` or `.txt` and comma-separated
otherwise, and a header line is skipped. Barcodes must be written as in the matrix (with `-1` for
10x). The methods that can count reads take `--unit reads`:

```bash
sc-droplets call umitools-distance raw_feature_bc_matrix.h5 --reads reads.tsv --unit reads --barcodes cells.tsv.gz
sc-droplets call alevinfry-permit-list raw_feature_bc_matrix.h5 --reads reads.tsv --unit reads --expect-cells 3000 --json
sc-droplets compare raw_feature_bc_matrix.h5 -m umitools-distance:unit=reads -m umitools-distance --reads reads.tsv
```

In Python, `reads` is an `obs` column or an array:

```python
import pandas as pd
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
reads = pd.read_csv("reads.tsv", sep="\t", index_col=0).iloc[:, 0]
adata.obs["n_reads"] = reads.reindex(adata.obs_names).fillna(0).to_numpy()
scd.call(adata, "zumis-cell-barcodes", reads="n_reads", unit="reads")
```

## Where the reads come from

| Pipeline | Reads per barcode |
|---|---|
| Cell Ranger 3+ | `molecule_info.h5`: the sum of each barcode's molecules' read counts (only the reads Cell Ranger counted) |
| STARsolo | `CellReads.stats` (written with `--soloCellReadStats Standard`), the reads matched to each barcode |
| UMI-tools, zUMIs, Drop-seq | the tool's own per-barcode read counts, as it wrote them |

From `molecule_info.h5`, for a run with one GEM group:

```python
import h5py
import numpy as np
import pandas as pd

with h5py.File("molecule_info.h5") as f:
    barcodes = f["barcodes"][:].astype(str)
    reads = np.bincount(f["barcode_idx"][:], weights=f["count"][:], minlength=len(barcodes))
table = pd.DataFrame({"barcode": [b + "-1" for b in barcodes], "reads": reads.astype(int)})
table[table["reads"] > 0].to_csv("reads.tsv", sep="\t", index=False)
```

These are not exactly the reads the read-based tool counted. UMI-tools counts barcodes in the
first 100 million reads, before any mapping, and alevin-fry counts the reads it mapped. Read
counts from another pipeline give a close curve, not the same one; each method's *Differences
from the original* says what it counted.

## Which methods take reads

| Method | Counts by default in the original | `--unit reads` |
|---|---|---|
| [UMI-tools distance](../methods/umitools/distance.md) and [density](../methods/umitools/density.md) | reads (or UMIs with `--method umis`) | yes |
| [alevin initial whitelist](../methods/alevin/initial-whitelist.md) | reads | yes |
| [alevin-fry permit list](../methods/alevinfry/permit-list.md) | reads | yes |
| [bustools allowlist](../methods/bustools/allowlist.md) | reads | yes |
| [zUMIs](../methods/zumis/cell-barcodes.md) | reads | yes |
| [Drop-seq (dropbead)](../methods/dropseq/dropbead.md) | reads | yes |
| [Kneedle](../methods/kneedle/kneedle.md) | as given | yes |
| the [generic rules](../methods/generic/index.md) | | `umis`, `reads` or `genes` |
