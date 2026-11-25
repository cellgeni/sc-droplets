---
title: Methods
nav_order: 4
has_children: true
---

# Methods

One page per algorithm, grouped by the tool it comes from. Each page has the same sections:
what the method is at a glance, how it works step by step, which releases of its tool behaved
how, its parameters beside the original tool's names for them, how to run it, where our
implementation differs from the original, and how it was checked against the original.

The table below lists the methods the package runs. It is generated from the package
(`docs/_scripts/methods_table.py`); the [Overview](../overview/index.md) lists those still to
come.

<!-- methods-table:start -->
| Method | Tool | Category | Versions | Page |
|---|---|---|---|---|
| `threshold` | Generic | rank-threshold | — | [Fixed threshold](generic/threshold.md) |
| `top-n` | Generic | rank-threshold | — | [Top N barcodes](generic/top-n.md) |
<!-- methods-table:end -->

## Reading a method name

Method ids are `tool-algorithm`. A spec adds the release of the tool to reproduce and any
parameters:

```
cellranger-emptydrops@9.0
dropletutils-emptydrops@1.26:lower=200
top-n:n_cells=3000
```

On the command line, `sc-droplets call <id>` takes the release as `--tool-version` and the
parameters as options; `compare` and `plot` take specs with `-m`.
