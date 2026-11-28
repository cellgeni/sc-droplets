# sc-droplets

Every method for calling cells vs empty droplets in droplet single-cell data, in one place,
version by version.

Each pipeline decides which barcodes are cells from the barcode-rank ("knee") curve, and each
does it differently: Cell Ranger, STARsolo, DropletUtils' emptyDrops, CellBender, UMI-tools,
alevin-fry, kb, Drop-seq tools and more. Most changed their rule between releases. sc-droplets
reimplements them so that any of them, as any release of its tool ran it, can be run on the
same matrix, and compared.

- **Documentation and the review of methods**: <https://cellgeni.github.io/sc-droplets/>
- **Python**: works on AnnData.
- **Command line**: works on `.h5ad`, 10x `.h5`, MTX folders, and the output folders of Cell
  Ranger, STARsolo, kb and alevin-fry.

> Status: 0.1.0, the first release. It runs 29 methods: those of 13 tools and the generic rules,
> in 76 release-specific behaviours, each checked against its original tool. 22 more methods are
> described on the site, planned for 0.2 and 0.3 or documented only. See [STATUS.md](STATUS.md).

## Install

```bash
pip install sc-droplets
```

Python 3.12 or later. The command is `sc-droplets`; the import is `sc_droplets`.

## Quick start

```bash
sc-droplets methods                                  # what is available
sc-droplets inspect raw_feature_bc_matrix.h5         # a summary and a knee in the terminal
sc-droplets call top-n raw_feature_bc_matrix.h5 --n-cells 3000 --barcodes cells.tsv.gz
sc-droplets compare raw_feature_bc_matrix.h5 -m cellranger -m otsu -m gmm --plot knee.png
sc-droplets report raw_feature_bc_matrix.h5 -m cellranger -m otsu -m gmm -o report.html
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")       # raw, unfiltered counts
scd.call(adata, "threshold", min_counts=500)       # stored in adata.obs["threshold"]
scd.call(adata, "cellranger")
scd.pl.compare(adata)                               # every stored call, side by side
```

## Licence

GPL-3.0-or-later. Several methods are ports of GPL code (DropletUtils, alevin, zUMIs, dropEst).
Methods whose original code is under a licence that does not allow copying (Cell Ranger 2.x and
7.2 onwards) are implemented from a written specification on their documentation page, not from
their code.
