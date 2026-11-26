---
title: alevin-fry
parent: Methods
nav_order: 9
has_children: true
---

# alevin-fry
{: .no_toc }

alevin-fry (He et al. 2022) quantifies droplet data from reads mapped by salmon or piscem.
`alevin-fry generate-permit-list` decides which barcodes are cells before anything is counted:
it counts the reads of each barcode and keeps those above a threshold, found by one of four
rules. simpleaf runs it with whichever rule it is given.

1. TOC
{:toc}

## Which release does what

| alevin-fry | What changed |
|---|---|
| 0.1.0 (Dec 2020) | the knee (a distance knee like UMI-tools', with a slack of 5), `--expect-cells`, `--force-cells`; the rules have not changed for single-barcode data since |
| 0.4.2 (Oct 2021) | `--unfiltered-pl`: every barcode of a given list with at least `--min-reads` reads |
| 0.12.0 (Mar 2026) | multi-barcode input (10x Flex, sample and cell barcodes): its permit lists use rules of their own, one rank higher for the knee and a different `--expect-cells` |
| 0.18.0 (Aug 2026) | one rule set for both: multi-barcode input gets the single-barcode rules again |

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [permit list](permit-list.md) | `alevinfry-permit-list` | the knee, expected cells, forced cells or a read floor; above the threshold, a cell |
