---
title: distance
parent: UMI-tools
grand_parent: Methods
nav_order: 1
method_ids: [umitools-distance]
status: implemented
---

# UMI-tools whitelist, distance knee
{: .no_toc }

`umi_tools whitelist --knee-method=distance`, the default from 1.0.0. On the cumulative count
curve, the knee is the point farthest from the straight line joining the curve's ends. The
search is repeated on the first three times as many barcodes until the knee stops moving, since
the long tail of barcodes with a read or two pushes the first knee too far out.

Why the rule lands where it does, drawn on real data:
[UMI-tools' distance knee](../../explained/knee/umitools.md), in *Explained*.

1. TOC
{:toc}

## At a glance

| Id | `umitools-distance` |
|---|---|
| Tool | UMI-tools 1.0.0 and later (`whitelist --knee-method=distance`) |
| Category | knee |
| Counts | reads per barcode in UMI-tools; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | MIT |
| Validated against | `umi_tools.whitelist_methods` 1.0.0 and 1.1.6 (exact) |

## How it works

Barcodes with a count of 0 are left out. The others are sorted by decreasing count, ties in
input order (UMI-tools sorts with `Counter.most_common`, which keeps the order barcodes were
first seen; here, the matrix's order). With counts $$c_0 \ge c_1 \ge \dots \ge c_{n-1}$$
(0-based), the cumulative counts are $$v_i = c_0 + \dots + c_i$$.

**The knee of the first $$m$$ points.** The points are $$P_i = (i, v_i)$$, $$i < m$$. With
$$u = (P_{m-1} - P_0) / \lVert P_{m-1} - P_0 \rVert$$ the direction of the chord, each point's
distance from the chord is the length of the part of $$P_i - P_0$$ perpendicular to it:

$$d_i = \left\lVert (P_i - P_0) - \big((P_i - P_0) \cdot u\big)\, u \right\rVert$$

and the knee is the first $$i$$ with the largest $$d_i$$.

1. $$k$$ = the knee of all $$n$$ points. If $$k = 0$$, UMI-tools stops with an error, and so
   does this.
2. Repeat, at most 100 times: $$k \leftarrow$$ the knee of the first $$\min(3k, n)$$ points;
   stop when $$k$$ no longer changes.
3. The cells are the first $$k + 1$$ barcodes of the sorted order. A barcode tied with the last
   cell may be left out: ties fall in input order.

With `--n-cells N` (`--set-cell-number`), the knee is still found (and step 1 can still fail),
but the cells are the barcodes with a count strictly above $$c_N$$, the $$(N+1)$$-th largest:
N cells, or fewer when $$c_{N-1} = c_N$$.

### Arithmetic

Counts and cumulative counts are 64-bit integers, as in UMI-tools, so
$$\lVert P_{m-1} - P_0 \rVert^2$$ is exact before its square root; the distances are computed in
the same order of floating-point operations as UMI-tools' numpy code.

## Versions

| UMI-tools | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| before 1.0.0 | — | — | no distance knee (see [density](density.md)) |
| 1.0.0 – 1.1.6 | `v1.0` | `1.0.0` | as above |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--n-cells` | none | `--set-cell-number` | call the barcodes above the (N+1)-th count instead of the knee |
| `--unit` | `umis` | `--method` | count `umis` (the matrix's totals) or `reads` (needs `--reads`) |

UMI-tools refuses `--expect-cells` with this rule; so does this method (it has no such option).

## Run it

```bash
sc-droplets call umitools-distance raw_feature_bc_matrix.h5 --barcodes whitelist.tsv
sc-droplets call umitools-distance raw_feature_bc_matrix.h5 --n-cells 1000 --json
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "umitools-distance")
result.n_cells, result.thresholds["knee_rank"]
```

The thresholds hold `knee_rank` (k + 1, the number of barcodes up to the knee) and the count of
the last cell in the unit used. `sc-droplets plot diagnostics` draws the rank curve with the
cells, and the cumulative curve with each pass's chord and knee.

## Differences from the original

- UMI-tools counts reads (or UMIs) itself from the first 100 million reads; this takes the
  counts from the matrix. Given the same counts in the same order, the calls are the same.
- UMI-tools' error correction of barcodes (`--error-correct-threshold`,
  `--ed-above-threshold`) is not part of this method; see [UMI-tools](index.md#counting).

## Validation

Compared with `umi_tools.whitelist_methods.getCellWhitelist(knee_method="distance")` from
UMI-tools 1.0.0 and 1.1.6 (each in its own environment from bioconda), given the per-barcode
UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices in matrix order, with and without
`--set-cell-number`: identical barcodes.

## References

- Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular
  Identifiers to improve quantification accuracy, *Genome Res* 27:491–499,
  [doi:10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116).
- Source: [`umi_tools/whitelist_methods.py`](https://github.com/CGATOxford/UMI-tools/blob/1.1.6/umi_tools/whitelist_methods.py)
  (`getKneeEstimateDistance`, 1.1.6) and
  [`whitelist.py`](https://github.com/CGATOxford/UMI-tools/blob/1.1.6/umi_tools/whitelist.py).
