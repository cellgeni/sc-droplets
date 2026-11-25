---
title: Home
nav_order: 1
permalink: /
---

# sc-droplets

Every droplet single-cell pipeline has to decide which barcodes are cells and which are empty
droplets holding only ambient RNA. Nearly all of them decide from the same picture, the
barcode-rank curve or "knee plot", but each in its own way, and most have changed their rule
between releases. Cell Ranger alone has used at least six.

This site and the `sc-droplets` package collect those methods in one place:

- a [review](overview/index.md) that sorts them into families and lists what has been compared
  before;
- [one page per algorithm](methods/index.md), with its exact rule, its version history and how
  far our implementation reproduces the original;
- a Python package that runs any of them on AnnData, and a command line that runs them on the
  files pipelines write.

```
barcodes ranked by counts
   │
   ├── cells           the plateau at the top
   ├── the knee        where cells end
   ├── empty droplets  a shoulder of ambient RNA
   └── the tail        barcodes with a handful of counts
```

{: .important }
sc-droplets is in early development. The first release, 0.1.0, will implement every knee and
threshold heuristic and every statistical test, each checked against its original tool. Until
then the package has the generic [top N](methods/generic/top-n.md) and
[fixed threshold](methods/generic/threshold.md) rules, and the
[Overview](overview/index.md) lists what is coming.

## Install

```bash
pip install sc-droplets
```

Python 3.12 or later. The command is `sc-droplets`; the Python package is `sc_droplets`.

## Quick start

```bash
sc-droplets methods
sc-droplets inspect raw_feature_bc_matrix.h5
sc-droplets call threshold raw_feature_bc_matrix.h5 --min-counts 500 --barcodes cells.tsv.gz
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m threshold:min_counts=500 --plot knee.png
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
scd.call(adata, "threshold", min_counts=500)
scd.pl.knee(adata)
```

[Getting started](getting-started.md) walks through these.

## Versions of a method

A method is named after the tool and algorithm it comes from, and takes the release of that tool
whose behaviour you want: `dropletutils-emptydrops@1.26` is emptyDrops as DropletUtils 1.26 ran
it. Without a version, the latest known behaviour is used. Releases whose code is not public, and
so cannot be reproduced, are refused with a pointer to the method's page.

## Licence

GPL-3.0-or-later. Several methods are ports of GPL code. Methods whose original code may not be
copied (Cell Ranger 2.x and 7.2 onwards) are implemented from the specification on their page,
not from their code.
