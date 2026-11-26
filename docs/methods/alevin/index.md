---
title: alevin
parent: Methods
nav_order: 11
has_children: true
---

# alevin
{: .no_toc }

alevin (Srivastava et al. 2019) was salmon's droplet pipeline from 0.11.0 to 1.10.3. Before
quantifying, it counts the reads of each barcode and takes an initial whitelist from the rank
curve: the barcodes above a knee, which it treats as cells, and a band of low-confidence barcodes
below them. A naive Bayes classifier, trained on those two groups, then decides the final cells.
salmon 1.11.0 dropped alevin in favour of alevin-fry.

1. TOC
{:toc}

## Which release does what

| salmon | What changed |
|---|---|
| 0.11.0 (Jul 2018) | the knee (a tangent from the origin to the log cumulative curve), corrected by a valley of a kernel density |
| 0.12.0 (Dec 2018) | the knee reads one barcode further up the curve (an off-by-one fixed); `--expectCells` (Cell Ranger 2's rule) and `--forceCells` added |
| 0.14.0 (May 2019) | forced cells below `--freqThreshold` reads are dropped |
| 0.14.1 (Jun 2019) | `--forceCells F` keeps F − 1 barcodes; the low-confidence band no longer runs past the last barcode |
| 1.2.0 (Apr 2020) | `--forceCells F` keeps F again, and `--expectCells` one barcode more |
| 1.10.0 – 1.10.3 | wider integers; the rules are unchanged |

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [initial whitelist](initial-whitelist.md) | `alevin-initial-whitelist` | the knee with its density correction, expected cells or forced cells; then the low-confidence band |

## Not implemented

The naive Bayes classifier that turns the initial whitelist into alevin's final cells is planned
for 0.2.
