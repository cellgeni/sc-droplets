---
title: initial whitelist
parent: alevin
grand_parent: Methods
nav_order: 1
method_ids: [alevin-initial-whitelist]
status: implemented
---

# alevin initial whitelist
{: .no_toc }

alevin's `sampleTrueBarcodes`: the barcodes it takes as cells before quantifying, and the band of
low-confidence barcodes below them. Each barcode gets one of three classes: `cell` (above the
knee, or within the expected or forced cells), `low_confidence` (the band), and `background`.
Only `cell` is called. alevin's naive Bayes classifier then decides the final cells from these
two groups; that step is not implemented here.

Why the rule lands where it does, drawn on real data:
[alevin's tangent](../../explained/knee/alevin.md), in *Explained*.

1. TOC
{:toc}

## At a glance

| Id | `alevin-initial-whitelist` |
|---|---|
| Tool | salmon 0.11.0 – 1.10.3 (`salmon alevin`, with `--expectCells`, `--forceCells` or `--keepCBFraction`, or none for the knee) |
| Category | knee |
| Counts | reads per barcode, before barcode correction, in alevin; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no (optional) |
| Deterministic | yes |
| Original licence | GPL-3.0 |
| Validated against | programs built from salmon 0.11.3, 0.13.1, 0.14.0, 1.1.0 and 1.10.3's own `Alevin.cpp` (exact, up to barcodes tied at the knee) |

## How it works

The barcodes with a count, sorted by decreasing count: $$f_0 \ge f_1 \ge \dots \ge f_{n-1}$$
(0-based). $$M = \min(n, \text{maxNumBarcodes})$$ (100,000 by default). Each rule gives $$k$$,
the number of cells.

**Knee** (with none of the options below; `getLeftBoundary`):

1. The log cumulative curve of the top $$M$$, summed from the smallest up:
   $$c_i = \log(f_{M-1} + f_{M-2} + \dots + f_{M-1-i})$$, $$i = 0, \dots, M - 1$$ (natural log).
   In 0.11.x the sum is shifted by one barcode, $$f_M + \dots + f_{M-i}$$ (an off-by-one fixed in
   0.12.0).
2. The tangent from the origin: the first $$j$$ with $$c_j \ne c_0$$ such that no later point lies
   below the line through the origin and $$(j, c_j)$$: $$i \cdot c_j / j \le c_i$$ for every
   $$i > j$$.
3. With $$m = M - j$$, $$k$$ is the number of barcodes with a count above $$f_m$$, less one: alevin
   steps back from rank $$m$$ past the barcodes tied with it, and stops at the first one above.
   If no $$j$$ qualifies, alevin stops with an error.

**Density correction** (`gaussianKDE`, after the knee):

1. $$\ell_i = \log_{10} f_i$$ for the counts above $$0.001\,f_0$$ ($$m'$$ of them). The kernel
   variance is $$\sigma^2 = 0.01 \cdot \sum_i (\ell_i - \bar\ell)^2 / (m' - 1)$$; alevin stops
   with an error if it is 0.
2. A grid of 10,000 points from $$\ell_0$$ down, in steps of
   $$d = (\ell_0 - \ell_{m'-1}) / 10000$$ (each point the previous one minus $$d$$), and the
   density $$D_t = \sum_i \exp(-(\ell_i - g_t)^2 / 2\sigma^2)$$ at each.
3. The local minima ($$D_{t-1} > D_t < D_{t+1}$$), from the highest count down. For each, the
   count $$\tau = 10^{\ell_0 - t d}$$ and $$b$$, the number of barcodes with
   $$f \ge \tau$$:
   - if $$b > k$$, the knee stands and the search stops;
   - if $$b < 0.1\,k$$, the next minimum is tried;
   - otherwise $$k = b$$.

**Expected cells** (`--expectCells` $$E$$, from 0.12.0; Cell Ranger 2's rule):
$$r = \max(1, \lfloor 0.01 E \rfloor)$$, cut $$= \max(1, 0.1 f_r)$$, and $$i$$ is the first rank
in $$[r, \min(n, 2E))$$ with $$f_i <$$ cut. Then $$k = i$$ (to 1.1.0) or $$k = i + 1$$ (from
1.2.0). If no rank qualifies (or, from 1.2.0, $$i + 1 = \min(n, 2E)$$), alevin stops with an error.

**Forced cells** (`--forceCells` $$F$$, from 0.12.0):

- 0.12.0 – 0.13.1: $$k = F$$.
- 0.14.0: $$k = F$$, then $$k$$ is lowered while $$f_k$$ is below `--freqThreshold` (10).
- 0.14.1 – 1.1.0: the same, starting from $$k = F - 1$$, so $$F$$ forced cells give $$F - 1$$.
- From 1.2.0: as 0.14.1, then $$k + 1$$.

`--keepCBFraction` $$p$$ (from 0.13.0) forces $$F = \min(\lfloor p n \rfloor, 737{,}000)$$, and from
0.14.1 $$\min(\lfloor p n \rfloor, M_\text{max})$$ with $$M_\text{max}$$ = `--maxNumBarcodes`.

**The low-confidence band.** Its width is $$\lfloor 0.5 k \rfloor$$, but at least
`--lowRegionMinNumBarcodes` (200) and at most 1,000. From 0.14.1 it stops at the last barcode.
With $$t = k$$ plus the width, the whitelist is every barcode with a count above $$f_t$$: alevin
steps back past the barcodes tied with $$f_t$$ "to imitate stable sort". The `cell` class is the
top $$k$$, and `low_confidence` the rest of the whitelist.

## Versions

| salmon | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.11.0 – 0.11.3 | `v0.11` | `0.11.0` | the knee only, with the off-by-one |
| 0.12.0 | `v0.12` | `0.12.0` | knee fixed; expected cells ($$k = i$$) and forced cells ($$k = F$$) |
| 0.13.0 – 0.13.1 | `v0.13` | `0.13.0` | `--keepCBFraction` (at most 737,000) |
| 0.14.0 | `v0.14` | `0.14.0` | forced cells below the read floor dropped |
| 0.14.1 – 1.1.0 | `v0.14.1` | `0.14.1` | forced cells start at $$F - 1$$; the band stops at the last barcode; `--keepCBFraction` at most `--maxNumBarcodes` |
| 1.2.0 – 1.10.3 | `v1.2` | `1.2.0` | forced cells $$F$$ again; expected cells $$k = i + 1$$ |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expect-cells` | none | `--expectCells` | the rule from expected cells |
| `--force-cells` | none | `--forceCells` | the rule from forced cells |
| `--keep-cb-fraction` | none | `--keepCBFraction` | force this fraction of the barcodes seen |
| `--freq-threshold` | 10 | `--freqThreshold` | forced cells need at least this count (0.14.0 on) |
| `--max-num-barcodes` | 100,000 | `--maxNumBarcodes` | the knee looks at this many barcodes at most |
| `--low-region-min-num-barcodes` | 200 | `--lowRegionMinNumBarcodes` | the band's least width |
| `--unit` | `umis` | reads per barcode | count UMIs (the matrix) or reads (`--reads`) |

## Run it

```bash
sc-droplets call alevin-initial-whitelist raw_feature_bc_matrix.h5 --json
sc-droplets call alevin-initial-whitelist raw_feature_bc_matrix.h5 --tool-version 1.1.0 --force-cells 3000
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "alevin-initial-whitelist", expect_cells=5000)
adata.obs["alevin_initial_whitelist_class"].value_counts()
```

The thresholds hold $$k$$ (`cells`), the knee before the density correction (`knee_rank`), the
correction's count where it found one (`gauss_rank`), the whitelist's size (`whitelist`), and the
band's width as alevin reports it (`low_confidence`). `sc-droplets plot diagnostics` draws the
rank curve with the cells and the cut, and the kernel density with its minima.

## Differences from the original

- alevin counts the reads of each barcode as sequenced, before barcode correction. A matrix
  holds corrected barcodes, and sc-droplets counts their UMIs by default; reads need `--reads`.
- When equal counts straddle the knee, which of them alevin calls cells depends on the order its
  hash map lists the barcodes in, and on the C++ library's `std::sort`. Neither can be had from a
  matrix: here the tied barcodes are taken in the matrix's order. The number of cells, the
  whitelist and every other barcode's class are as alevin's.
- In 0.11.x the knee reads one barcode past the top $$M$$; with $$M = n$$ that is past the end of
  alevin's array, so sc-droplets refuses 0.11.x when fewer than `--maxNumBarcodes` barcodes have a
  count.
- The density sums numpy's `exp`, barcode by barcode in alevin's order; alevin calls the C
  library's. The two can differ in the last bit, which could move a minimum only where two
  neighbouring grid points have densities equal to the last bits.
- `--whitelist` (a list given) and the classifier are not offered.

## Validation

`oracles/alevin/run.sh` builds a small C++ program per release from salmon's own `Alevin.cpp`
(fetched at the tag): `sort_indexes`, `getLeftBoundary`, `gaussianKDE` and `sampleTrueBarcodes` as
they are, with `AlevinOpts` and the logger replaced by stand-ins. It runs them on the per-barcode
UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices for 0.11.3, 0.13.1, 0.14.0, 1.1.0 and
1.10.3: the knee, expected cells, forced cells and `--keepCBFraction`. The number of cells, the
knee, the density correction and the whitelist are identical; the cells are identical up to the
barcodes tied at the knee.

## References

- Srivastava, Malik, Smith, Sudbery & Patro 2019, Alevin efficiently estimates accurate gene
  abundances from dscRNA-seq data, *Genome Biology* 20:65,
  [doi:10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y).
- Source: [`Alevin.cpp` at 1.10.3](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/Alevin.cpp)
  (`sort_indexes`, `gaussianKDE`, `getLeftBoundary`, `sampleTrueBarcodes`), and at
  [0.14.1](https://github.com/COMBINE-lab/salmon/blob/v0.14.1/src/Alevin.cpp),
  [0.12.0](https://github.com/COMBINE-lab/salmon/blob/v0.12.0/src/Alevin.cpp) and
  [0.11.3](https://github.com/COMBINE-lab/salmon/blob/v0.11.3/src/Alevin.cpp).
