---
title: allowlist
parent: bustools
grand_parent: Methods
nav_order: 1
method_ids: [bustools-allowlist]
status: implemented
---

# bustools allowlist
{: .no_toc }

`bustools allowlist` (`whitelist` before 0.43.1), which `kb count --filter bustools` runs. It reads
the first barcodes of the BUS file, takes the mean count of the ten largest among them, and keeps
every barcode with at least a fraction of that mean. The fraction is the chance that a barcode has
a sequencing error, at 1 % per base: about 0.149 for 16 bases.

1. TOC
{:toc}

## At a glance

| Id | `bustools-allowlist` |
|---|---|
| Tool | bustools 0.39.3 and later (`bustools allowlist`, `whitelist`, `onlist`; `kb count --filter bustools`) |
| Category | rank-threshold |
| Counts | reads per barcode (the counts of its BUS records); here UMIs, or reads with `--unit reads` |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | BSD-2-Clause |
| Validated against | bustools 0.39.3, 0.41.0, 0.43.2 and 0.45.1 (exact) |

## How it works

A BUS file sorted by barcode holds one record per barcode, UMI and equivalence class, each with
the number of reads it collapses. kb runs the allowlist on the sorted and corrected file. Here
each UMI of a barcode is one record, and the barcode's count is its UMIs (or its reads, with
`--unit reads`).

1. The barcodes with a count, in the order of the file: their sequences in A < C < G < T order,
   which is the order of bustools' 2-bit code. $$L$$ is the barcode length, without a suffix such
   as `-1`.
2. The sample. bustools reads the records in blocks of 100,000. At the end of each block, every
   barcode whose records ended in it is in the sample with its count. So is the barcode still
   open at the block's end, with the count of its records so far. Reading stops after the first
   block that brings the sample to 200 entries, or at the end of the file. A barcode open across
   blocks can enter the sample twice. With 10x data the sample is every barcode in the first
   100,000 records, thousands of them, not 200.
3. The sample sorted by decreasing count (ties by barcode). With $$a$$ = the sum of its first ten
   counts, divided by 10 and truncated, the threshold is
   $$t = \lfloor a\,(1 - 0.99^{L}) \rfloor$$.
4. The cells are the barcodes with count $$\ge t$$. `--threshold` sets $$t$$ instead.

bustools prints a barcode twice when it is open at the end of the sample's last block and passes
the threshold; the list of cells is the same.

## Versions

| bustools | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.39.3 – 0.45.1 | `default` | | one rule |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--threshold` | none | `--threshold` | the threshold itself |
| `--unit` | `umis` | reads per barcode | count UMIs (the matrix) or reads (`--reads`) |

## Run it

```bash
sc-droplets call bustools-allowlist raw_feature_bc_matrix.h5 --json
sc-droplets call bustools-allowlist counts_unfiltered/ --barcodes allowlist.txt
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "bustools-allowlist")
result.thresholds["umis"], result.diagnostics["sample_size"]
```

The thresholds hold $$t$$ (in the unit counted) and the mean of the top ten (`top10_mean`).
`sc-droplets plot diagnostics` draws the rank curve with the cut.

## Differences from the original

- A matrix does not hold the BUS file's records. Each UMI stands for one record here, whereas
  bustools sees one record per UMI and equivalence class, so a UMI matching several classes
  counts several times. Where the first 100,000 records end, and so the sample, can differ
  slightly from the file kb used.
- The matrix's barcodes are already corrected: bustools, run by kb on the corrected file, sees
  the same barcodes.
- Before 0.42.0 the current barcode started uninitialised. From 0.42.0 it starts as 0, the code
  of a barcode of all A's: such a barcode, if first in the file, is counted one read short (as
  here). Before, it depended on the memory's contents.

## Validation

`oracles/bustools/run.sh` writes a BUS file from the matrix, one record per UMI in barcode order.
It runs bioconda's bustools 0.39.3, 0.41.0, 0.43.2 and 0.45.1 (`bustools whitelist`, the name every
release accepts) on 10x's pbmc_1k_v3 and pbmc4k raw matrices, with and without `--threshold`.
Identical lists of barcodes and thresholds.

## References

- Melsted, Ntranos & Pachter 2019, The barcode, UMI, set format and BUStools, *Bioinformatics*
  35:4472–4473, [doi:10.1093/bioinformatics/btz279](https://doi.org/10.1093/bioinformatics/btz279).
- Source: [`bustools_whitelist.cpp`](https://github.com/BUStools/bustools/blob/v0.45.1/src/bustools_whitelist.cpp)
  (0.45.1; the same rule since 0.39.3).
