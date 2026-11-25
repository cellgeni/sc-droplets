---
title: soloBasicCellFilter.awk
parent: STARsolo
grand_parent: Methods
nav_order: 4
method_ids: [starsolo-basic-awk]
status: implemented
---

# STARsolo soloBasicCellFilter.awk
{: .no_toc }

Before STAR filtered cells itself (2.7.3a), it shipped an awk script that filters a raw
`matrix.mtx`. Its cut is CellRanger2.2's, without rounding.

1. TOC
{:toc}

## At a glance

| Id | `starsolo-basic-awk` |
|---|---|
| Tool | STAR 2.7.2a and later, `extras/scripts/soloBasicCellFilter.awk` |
| Category | rank-threshold |
| Counts | the sum of each barcode's matrix entries |
| Needs expected cells | uses 3,000 unless told otherwise |
| Deterministic | yes |
| Original licence | MIT |
| Validated against | the script's source (unchanged since 2.7.2a); unit tests |

## How it works

With $$L$$ barcodes that have matrix entries and their totals sorted in increasing order
$$t_1 \le \dots \le t_L$$ (1-based, as awk's `asort`):

1. The robust maximum is $$m = t_{L - \lfloor (1 - q)\,\text{maxCells} \rfloor}$$, the 31st largest
   for the defaults.
2. The cut is $$c = m / \text{maxMinRatio}$$, unrounded.
3. With `exactCells` > 0 instead: $$c = t_1$$ if $$L < \text{exactCells}$$, else
   $$c = t_{L - \text{exactCells}}$$, keeping exactCells + 1 barcodes plus ties.
4. Every barcode with total $$\ge c$$ is a cell. awk reads a barcode with no entries as 0, so a
   cut of 0 keeps every barcode in the list.

## Versions

| STAR | Variant | Behaviour |
|---|---|---|
| 2.7.2a – 2.7.11b | `v1` | as above; the script has not changed |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--max-cells` | 3000 | `maxCells` | expected number of cells |
| `--max-percentile` | 0.99 | `maxPercentile` | robust maximum percentile |
| `--max-min-ratio` | 10 | `maxMinRatio` | robust maximum ÷ this is the cut |
| `--exact-cells` | 0 | `exactCells` | if above 0, keep this many cells (plus one, plus ties) instead |

## Run it

```bash
sc-droplets call starsolo-basic-awk Solo.out --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("Solo.out")
scd.call(adata, "starsolo-basic-awk")
```

## Differences from the original

None known.

## Validation

Unit tests check the cut against values worked out from the script by hand, including
`exactCells`, fewer barcodes than `exactCells`, and a cut of 0.

## References

- Source: [`soloBasicCellFilter.awk`](https://github.com/alexdobin/STAR/blob/2.7.11b/extras/scripts/soloBasicCellFilter.awk).
