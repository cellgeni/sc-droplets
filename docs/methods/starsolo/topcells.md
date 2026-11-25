---
title: TopCells
parent: STARsolo
grand_parent: Methods
nav_order: 2
method_ids: [starsolo-topcells]
status: implemented
---

# STARsolo TopCells
{: .no_toc }

`--soloCellFilter TopCells N` is STAR's fixed-number filter. As written it keeps N + 1 barcodes,
plus any that tie with the last, and before 2.7.6a it kept every barcode.

1. TOC
{:toc}

## At a glance

| Id | `starsolo-topcells` |
|---|---|
| Tool | STAR 2.7.3a and later (`--soloCellFilter TopCells N`) |
| Category | rank-threshold |
| Counts | UMIs of the detected barcodes |
| Needs expected cells | yes: N |
| Deterministic | yes |
| Original licence | MIT |
| Validated against | `STAR --runMode soloCellFiltering` 2.7.10a and 2.7.11b (exact) |

## How it works

With the detected totals sorted in decreasing order, $$u_{(0)} \ge u_{(1)} \ge \dots$$, and nCB
detected barcodes:

1. The cut is $$c = u_{(k)}$$ with $$k = \min(nCB - 1, N)$$: the (N+1)-th largest, as the index is
   0-based.
2. From 2.7.5a, $$c = \max(c, 1)$$.
3. Every detected barcode with total $$\ge c$$ is a cell: N + 1 of them, plus ties.

## Versions

| STAR | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 2.7.3a – 2.7.4a | `v1` | `2.7.3a` | $$k = \max(nCB - 1, N)$$: while N < nCB − 1 the cut is the smallest total, so every detected barcode is kept; with N ≥ nCB, STAR reads past its array (refused here) |
| 2.7.5a – 2.7.5c | `v2` | `2.7.5a` | as v1, cut floored at 1 |
| 2.7.6a – 2.7.11b | `v3` | `2.7.6a` | $$k = \min(nCB - 1, N)$$, floored at 1 |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--n-cells` | required | `TopCells N` | N |

## Run it

```bash
sc-droplets call starsolo-topcells Solo.out --n-cells 3000 --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("Solo.out")
scd.call(adata, "starsolo-topcells", n_cells=3000)
```

For exactly N cells, use the generic [top N](../generic/top-n.md).

## Differences from the original

None known. Where STAR 2.7.3a – 2.7.5c would read past the end of its array (N ≥ nCB), this
raises an error instead.

## Validation

Compared with `STAR --runMode soloCellFiltering --soloCellFilter TopCells 3000` (2.7.10a, 2.7.11b)
on two public 10x matrices; identical. The older variants are read from their source.

## References

- Source: [`SoloFeature_cellFiltering.cpp`](https://github.com/alexdobin/STAR/blob/2.7.11b/source/SoloFeature_cellFiltering.cpp)
  (2.7.11b), [2.7.5c](https://github.com/alexdobin/STAR/blob/2.7.5c/source/SoloFeature_cellFiltering.cpp),
  [2.7.6a](https://github.com/alexdobin/STAR/blob/2.7.6a/source/SoloFeature_cellFiltering.cpp).
