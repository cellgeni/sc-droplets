---
title: cell barcodes
parent: zUMIs
grand_parent: Methods
nav_order: 1
method_ids: [zumis-cell-barcodes]
status: implemented
---

# zUMIs cell barcodes
{: .no_toc }

zUMIs' automatic cell detection (`cellBC` → `.cellBarcode_unknown` → `.FindBCcut`, 2.5.6b and
later): the unit invariant knee (`inflection::uik`) of the cumulative reads of the barcodes,
ranked by reads; every barcode with at least the knee's count is a cell.

1. TOC
{:toc}

## At a glance

| Id | `zumis-cell-barcodes` |
|---|---|
| Tool | zUMIs 2.5.6b and later (automatic barcode detection, with no `barcode_file` or `barcode_num`), with inflection 1.3.5 or later |
| Category | knee |
| Counts | reads per barcode in zUMIs; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | GPL-3.0 (zUMIs), GPL-2.0-or-later (inflection) |
| Validated against | zUMIs 2.5.6b and 2.9.7's `cellBC`, with inflection 1.3.5 and 1.3.7, in R (exact) |

## How it works

1. The barcodes with at least `nReadsperCell` reads (100), sorted by decreasing count
   $$n_1 \ge n_2 \ge \dots \ge n_N$$, and their cumulative reads $$s_i = \sum_{k \le i} n_k$$ (R's
   `cumsum`, in long double).
2. `uik(1:N, s / 1000)`:
   - `check_curve` decides whether the curve is concave: at ranks 1, and
     $$\lfloor q_{0.25} \rfloor$$, $$\lfloor q_{0.5} \rfloor$$, $$\lfloor q_{0.75} \rfloor$$ of
     $$1, \dots, N$$ (R's quantile type 7), and $$N$$, the trapezoid integrals of the curve minus
     its chord to the left and to the right; their signs give the curve's type.
   - For a concave curve (cumulative reads are), $$y$$ is negated. With the chord
     $$\ell(x) = y_1 + (y_N - y_1)(x - x_1)/(x_N - x_1)$$, the knee is $$x$$ at the first minimum
     of $$y - \ell(x)$$: the point farthest above the chord.
3. The cut is $$n$$ at the knee's rank; the cells are the barcodes with $$n \ge$$ the cut.
4. With fewer than 10 cells, zUMIs takes the top 100 barcodes instead (all of them if there are
   fewer).

## Versions

| zUMIs | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.0.2 – 2.5.4 | gap | | a Gaussian mixture initialised on an unseeded random subset: refused |
| 2.5.6b – 2.9.7 | `v2.5.6b` | `2.5.6b` | the UIK |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--n-reads-per-cell` | 100 | `barcodes: nReadsperCell` | barcodes with fewer reads are left out |
| `--unit` | `umis` | reads per barcode | count UMIs (the matrix) or reads (`--reads`) |

## Run it

```bash
sc-droplets call zumis-cell-barcodes raw_feature_bc_matrix.h5 --json
sc-droplets call zumis-cell-barcodes raw_feature_bc_matrix.h5 --n-reads-per-cell 50 --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "zumis-cell-barcodes")
result.thresholds["umis"], result.thresholds["knee_rank"]
```

The thresholds hold the cut (`umis` or `reads`) and the knee's rank (`knee_rank`).
`sc-droplets plot diagnostics` draws the rank curve with the cut.

## Differences from the original

- zUMIs counts reads per barcode; sc-droplets counts UMIs by default, reads with `--reads`.
- When fewer than 10 cells pass, zUMIs keeps the top 100 in its table's order (a stable sort of
  the barcodes as `fread` reads them); here equal counts keep the matrix's order.
- `barcode_file`, `barcode_num` and their combination with the automatic detection are not
  offered: they are lists and numbers given, not rules.

## Validation

`oracles/zumis/run.sh` sources zUMIs' `barcodeIDFUN.R` from tags 2.5.6b and 2.9.7, with
inflection 1.3.5 and 1.3.7 from CRAN, in R from conda-forge. It calls `cellBC` on a file of the
per-barcode UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices, with `nReadsperCell` 100 and
10. Identical cells.

## References

- Parekh, Ziegenhain, Vieth, Enard & Hellmann 2018, zUMIs - A fast and flexible pipeline to
  process RNA sequencing data with UMIs, *GigaScience* 7(6),
  [doi:10.1093/gigascience/giy059](https://doi.org/10.1093/gigascience/giy059).
- Christopoulos 2016, Introducing Unit Invariant Knee (UIK) As an Objective Choice for Elbow Point
  in Multivariate Data Analysis Techniques, SSRN,
  [doi:10.2139/ssrn.3043076](https://doi.org/10.2139/ssrn.3043076).
- Source: [`barcodeIDFUN.R`](https://github.com/sdparekh/zUMIs/blob/2.9.7/barcodeIDFUN.R)
  (`cellBC`, `.cellBarcode_unknown`, `.FindBCcut`); inflection's
  [`uik.R`](https://github.com/cran/inflection/blob/1.3.7/R/uik.R),
  [`check_curve.R`](https://github.com/cran/inflection/blob/1.3.7/R/check_curve.R),
  [`ede.R`](https://github.com/cran/inflection/blob/1.3.7/R/ede.R),
  [`findipl.R`](https://github.com/cran/inflection/blob/1.3.7/R/findipl.R).
