---
title: bustools
parent: Methods
nav_order: 10
has_children: true
---

# bustools
{: .no_toc }

bustools (Melsted et al. 2019) works on BUS files, the barcode, UMI and equivalence-class records
that kallisto writes. `bustools allowlist` (named `whitelist` until 0.43.1, and still accepted
under that name) makes the list of barcodes taken as cells, and `kb count --filter bustools`
runs it to make a filtered matrix. It derives one count threshold from a sample at the start of
the file.

1. TOC
{:toc}

## Which release does what

| bustools | What changed |
|---|---|
| 0.39.3 (Jul 2019) | `whitelist` added; the rule has not changed since |
| 0.42.0 (Dec 2022) | a variable left uninitialised is set to 0: only a first barcode of all A's is affected |
| 0.43.1 (Nov 2023) | also named `allowlist` (and `onlist` from 0.44.0) |

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [allowlist](allowlist.md) | `bustools-allowlist` | a fraction of the mean of the top 10 counts, in a sample from the start of the file |
