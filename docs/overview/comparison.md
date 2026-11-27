---
title: Comparison on public data
parent: Overview
nav_order: 3
---

# Comparison on public data
{: .no_toc }

Every method in this package, run with its defaults on the raw matrices of five public 10x
samples, and set against the call Cell Ranger made when 10x processed each one. Made by
`docs/_scripts/comparison.py`; the tables behind the figures are in `docs/overview/comparison/`.

1. TOC
{:toc}

## How to read it

- **10x's filtered matrix** is the reference: the cells Cell Ranger called when 10x ran it. It
  is one method among the others, not the truth.
- **Cells (Jaccard)**: how many barcodes a method calls, and the Jaccard index of its cell set
  with 10x's (1: the same barcodes).
- **Separation scores**: the median mitochondrial fraction (and, for PBMC 1k, the unspliced
  fraction from kb's nac counts) in a method's cells and in the barcodes it left out with at
  least 10 gene-expression counts, and the AUROC between them. In whole-cell samples, empty droplets hold mostly
  cytoplasmic, spliced RNA: a call that keeps them has a lower unspliced fraction among its cells
  and a lower AUROC. These are indicators, not ground truth.
- Defaults matter. Cell Ranger 2 and 3 run here with the expected cells 10x gave them; every other
  method with its own defaults, which is how most pipelines run them. Methods that count reads by
  default count UMIs here: a raw matrix has no reads per barcode.

## What it shows

**Cell Ranger here calls the cells 10x's own runs called.** At the release 10x ran, it calls the
same barcodes on PBMC 4k, 1k and 5k. On the spinal cord it calls 15,266 cells where 10x called
15,257: Cell Ranger 10.0 picks its ambient droplets by numpy's default sort, and which of the
barcodes tied at that rank fall inside depends on the machine 10x ran on, so a few barcodes at
the FDR cut go the other way. On the Flex sample it calls 7 more than 10x. Cell Ranger runs its
high-occupancy filter on the four samples of the plex together and removes those 7 there; this
page has one sample. The copies of its EmptyDrops step agree with it on whole cells and nuclei
(Jaccard 0.99 or more): STARsolo's EmptyDrops_CR and DropletUtils' `emptyDropsCellRanger`. On
Flex they part: 8,369 cells for STARsolo's, 7,605 for DropletUtils'.

**On whole cells, most rules land close together.** On the three PBMC samples, 17 of the 26
methods are within Jaccard 0.9 of 10x's call on all three. Four more stay within 0.85 on all
three: emptyDrops, UMI-tools' distance knee, cellsweep and multi-Otsu. They keep the same large
cells and differ in how far down the curve they go. Where they disagree, the
[Visualising calls](../guides/visualising.md) guide's quality plot shows the barcodes in
question.

**On nuclei, the methods split in three.** The spinal cord's curve has two cliffs. The first,
at about 15,000 barcodes, drops from the nuclei to a plateau near 250 counts. That plateau is the
empty droplets, heavy with ambient RNA, as nuclei preparations are. The second, at about 150,000
barcodes, drops from the plateau to barcodes with a few counts. Above the first cliff the curve
bends gradually, with no sharp knee.

- The knee rules stop in that gradual bend: 1,140 to 6,028 cells. These are OrdMag (Cell
  Ranger 2 and STARsolo's copy), `barcodeRanks`, `defaultDrops`, UMI-tools, dropbead, Kneedle,
  CellBender's prior and alevin.
- Cell Ranger's EmptyDrops and its copies, cellsweep, bustools, zUMIs, Otsu, the density valley
  and the mixtures find the first cliff: 13,500 to 20,000.
- Three rules find the second cliff and take in the empty droplets: alevin-fry's permit list, MAD
  and multi-Otsu, 154,000 to 160,000 barcodes. The sample has 158,080 barcodes with more than 100
  counts.

**emptyDrops at its defaults reaches its Monte Carlo limit on nuclei.** It tests every barcode
above 100 counts (158,080 here) with 10,000 simulations. Its smallest p-value is therefore
1/10,001. After the Benjamini-Hochberg correction over that many tests, the 12,988 barcodes
that beat every simulation still sit above FDR 0.001, and none is called. Its cells are the
2,290 above its knee, which it keeps without a test. Its `limited` column flags those barcodes.
DropletUtils' advice is to raise `--niters`. Cell Ranger's own EmptyDrops tests fewer barcodes:
at most 20,000 in 3.x, and from 8.0 only those above its ambient droplets' counts. From 9.0 it
also runs 100,000 simulations. It does not run into this limit here.

**The separation scores show what the disputed barcodes are, not which call is right.**

- **Whole cells, mitochondrial fraction.** Every call's cells carry more mitochondrial RNA than
  the barcodes it leaves out (AUROC 0.59 to 0.73). The calls that go furthest down the curve score
  highest. The calls that stop highest score lowest, since they leave cells out: Kneedle,
  CellBender's prior, zUMIs and `barcodeRanks`.
- **Whole cells, unspliced fraction (PBMC 1k).** Every call that leaves barcodes out separates
  them from its cells: its cells are 0.31 to 0.32 unspliced, the barcodes it leaves out 0.14 to
  0.15 (AUROC 0.84 to 0.92). Here the deepest calls
  score lowest: Cell Ranger 0.89, multi-Otsu 0.85, cellsweep 0.84, against 0.92 for the calls
  that stop between 1,200 and 2,900 counts.
- **Why the two measures disagree.** The barcodes the deeper calls add are mostly mitochondrial
  RNA, which has no introns. Cell Ranger 3's 114 rescued barcodes (median 972 counts) are 47 %
  mitochondrial and 25 % unspliced, against 10 % and 32 % in its 1,108 cells above the knee. The
  123 barcodes cellsweep takes beyond Cell Ranger (median 294 counts) are 90 % and 3 %. Damaged
  cells or debris, these are the barcodes the calls disagree on.
- **Nuclei.** The mitochondrial score flips: Cell Ranger's cells score 0.11. Nuclei carry less
  mitochondrial RNA than the ambient pool, which comes from lysed cytoplasm, so a lower score is
  the cleaner call. The rules that take in the empty droplets score about 0.5 (alevin-fry 0.52,
  MAD 0.55): no separation. Multi-Otsu's 0.93 means little, since it leaves out almost nothing
  above the floor.

**Four of the methods are steps of a larger method, not cell calls.** Their numbers say that,
not that the tool calls badly.

- Kneedle, with the sensitivity search yascp runs, finds the first sharp bend at the top of the
  curve: 321 to 1,515 cells on whole cells.
- CellBender's expected cells are a prior its model then refines. From 0.3.0 the prior is
  deliberately conservative: about half of 10x's cells.
- zUMIs calls about half of 10x's cells on PBMC, but 0.93 of them on nuclei.
- dropEst's bounds come from derivatives over 5 % of the barcodes it is given. dropEst passes only
  the barcodes it kept. Given every barcode of a raw 10x matrix, that window spans tens of
  thousands of ranks, and the estimate runs to 118,000 to 240,000.

## The calls

<!-- comparison:start -->

Cells called, and in brackets the Jaccard index of the cell set with 10x's own call:

| Method | PBMC 4k, v2 chemistry, Cell Ranger 2.1 | PBMC 1k, v3 chemistry, Cell Ranger 3.0 | PBMC 5k, GEM-X, Cell Ranger 9.0 | mouse spinal cord nuclei, Cell Ranger 10.0 | mouse lymph node, Flex (1 of a 4-plex), Cell Ranger 7.2 |
|---|---|---|---|---|---|
| 10x's filtered matrix | 4,340 | 1,222 | 5,710 | 15,257 | 7,579 |
| Cell Ranger, the release 10x ran | 4,340 (1.00) | 1,222 (1.00) | 5,710 (1.00) | 15,266 (1.00) | 7,586 (1.00) |
| `cellranger` | 4,606 (0.94) | 1,219 (0.99) | 5,710 (1.00) | 15,266 (1.00) | 7,578 (1.00) |
| `cellranger-ordmag` | 4,341 (1.00) | 1,116 (0.91) | 5,605 (0.98) | 5,924 (0.39) | 6,432 (0.85) |
| `starsolo-cellranger2.2` | 4,297 (0.99) | 1,135 (0.92) | 5,589 (0.98) | 5,228 (0.34) | 7,151 (0.89) |
| `starsolo-emptydrops-cr` | 4,591 (0.94) | 1,226 (0.99) | 5,710 (1.00) | 15,375 (0.99) | 8,369 (0.90) |
| `starsolo-basic-awk` | 4,297 (0.99) | 1,135 (0.92) | 5,589 (0.98) | 5,228 (0.34) | 7,151 (0.89) |
| `dropletutils-barcoderanks` | 3,519 (0.81) | 873 (0.71) | 4,775 (0.84) | 2,290 (0.15) | 5,191 (0.68) |
| `dropletutils-emptydrops` | 4,402 (0.85) | 1,166 (0.91) | 5,695 (0.98) | 2,290 (0.15) | 7,269 (0.94) |
| `dropletutils-emptydropscellranger` | 4,592 (0.94) | 1,226 (0.99) | 5,710 (1.00) | 15,378 (0.99) | 7,605 (0.99) |
| `dropletutils-defaultdrops` | 4,296 (0.99) | 1,135 (0.92) | 5,589 (0.98) | 5,228 (0.34) | 6,100 (0.80) |
| `umitools-distance` | 4,287 (0.99) | 1,089 (0.89) | 5,566 (0.97) | 3,554 (0.23) | 6,551 (0.86) |
| `umitools-density` | 4,770 (0.91) | 1,169 (0.94) | 5,640 (0.99) | 2,925 (0.19) | 7,267 (0.93) |
| `kneedle` | 936 (0.22) | 321 (0.26) | 1,299 (0.23) | 4,389 (0.29) | 1,515 (0.20) |
| `cellsweep-empty` | 5,033 (0.86) | 1,345 (0.91) | 5,836 (0.98) | 14,181 (0.91) | 8,912 (0.85) |
| `cellbender-prior` | 2,383 (0.55) | 641 (0.52) | 2,868 (0.50) | 6,028 (0.40) | 4,379 (0.57) |
| `alevinfry-permit-list` | 4,519 (0.96) | 1,144 (0.93) | 5,640 (0.99) | 154,250 (0.10) | 7,239 (0.93) |
| `bustools-allowlist` | 4,493 (0.97) | 1,161 (0.93) | 5,672 (0.99) | 13,571 (0.88) | 7,265 (0.93) |
| `alevin-initial-whitelist` | 4,763 (0.91) | 1,169 (0.94) | 5,640 (0.99) | 1,140 (0.07) | 7,267 (0.93) |
| `dropseq-dropbead` | 4,353 (1.00) | 1,217 (0.97) | 5,535 (0.97) | 4,506 (0.30) | 6,173 (0.81) |
| `dropest-bounds` | 117,894 (0.04) | 125,227 (0.01) | 239,902 (0.02) | 181,638 (0.08) | 18,426 (0.41) |
| `zumis-cell-barcodes` | 2,474 (0.57) | 594 (0.49) | 2,208 (0.39) | 15,092 (0.93) | 4,085 (0.54) |
| `mad` | 4,541 (0.96) | 1,160 (0.93) | 5,589 (0.98) | 155,698 (0.10) | 7,650 (0.92) |
| `otsu` | 4,489 (0.97) | 1,143 (0.92) | 5,654 (0.99) | 13,781 (0.89) | 6,968 (0.91) |
| `multiotsu` | 4,845 (0.90) | 1,317 (0.93) | 5,753 (0.99) | 159,891 (0.10) | 8,663 (0.87) |
| `kde-valley` | 4,768 (0.91) | 1,164 (0.94) | 5,690 (0.99) | 14,277 (0.91) | 7,347 (0.93) |
| `gmm` | 4,413 (0.98) | 1,121 (0.91) | 5,632 (0.99) | 19,885 (0.77) | 6,890 (0.90) |
| `gmm-2d` | 4,475 (0.97) | 1,151 (0.92) | 5,610 (0.98) | 18,284 (0.83) | 6,544 (0.86) |

Where a method fails on a dataset, its row on that dataset's table says why.

### PBMC 4k, v2 chemistry, Cell Ranger 2.1

![Calls on pbmc4k_2_1]({{ '/assets/comparison/pbmc4k_2_1.png' | relative_url }})

| Method | Cells | Jaccard with 10x | Smallest cell total | mito: cells | mito: others | mito: AUROC | Seconds |
|---|---|---|---|---|---|---|---|
| `cellranger@2.1.0:expect_cells=5000` | 4,340 | 1.000 | 1386 | 0.028 | 0.000 | 0.65 | 2 |
| `cellranger` | 4,606 | 0.942 | 500 | 0.029 | 0.000 | 0.67 | 19.6 |
| `cellranger-ordmag` | 4,341 | 1.000 | 1383 | 0.028 | 0.000 | 0.65 | 1.2 |
| `starsolo-cellranger2.2` | 4,297 | 0.990 | 1505 | 0.028 | 0.000 | 0.65 | 0.1 |
| `starsolo-emptydrops-cr` | 4,591 | 0.940 | 500 | 0.029 | 0.000 | 0.67 | 3 |
| `starsolo-basic-awk` | 4,297 | 0.990 | 1505 | 0.028 | 0.000 | 0.65 | 0.1 |
| `dropletutils-barcoderanks` | 3,519 | 0.811 | 2858 | 0.028 | 0.000 | 0.64 | 0.1 |
| `dropletutils-emptydrops` | 4,402 | 0.854 | 102 | 0.030 | 0.000 | 0.67 | 12.6 |
| `dropletutils-emptydropscellranger` | 4,592 | 0.940 | 501 | 0.029 | 0.000 | 0.67 | 1 |
| `dropletutils-defaultdrops` | 4,296 | 0.990 | 1506 | 0.028 | 0.000 | 0.65 | 0.1 |
| `umitools-distance` | 4,287 | 0.988 | 1533 | 0.028 | 0.000 | 0.65 | 0.1 |
| `umitools-density` | 4,770 | 0.910 | 373 | 0.029 | 0.000 | 0.67 | 0.6 |
| `kneedle` | 936 | 0.216 | 5137 | 0.030 | 0.000 | 0.65 | 0.2 |
| `cellsweep-empty` | 5,033 | 0.862 | 141 | 0.029 | 0.000 | 0.67 | 0.1 |
| `cellbender-prior` | 2,383 | 0.549 | 3700 | 0.028 | 0.000 | 0.64 | 0.2 |
| `alevinfry-permit-list` | 4,519 | 0.960 | 892 | 0.029 | 0.000 | 0.66 | 0.1 |
| `bustools-allowlist` | 4,493 | 0.966 | 973 | 0.029 | 0.000 | 0.66 | 0.2 |
| `alevin-initial-whitelist` | 4,763 | 0.911 | 387 | 0.029 | 0.000 | 0.67 | 0.5 |
| `dropseq-dropbead` | 4,353 | 0.997 | 1327 | 0.028 | 0.000 | 0.65 | 0.1 |
| `dropest-bounds` | 117,894 | 0.037 | 4 | 0.000 | — | — | 0.2 |
| `zumis-cell-barcodes` | 2,474 | 0.570 | 3638 | 0.028 | 0.000 | 0.64 | 0.1 |
| `mad` | 4,541 | 0.956 | 841 | 0.029 | 0.000 | 0.66 | 0.1 |
| `otsu` | 4,489 | 0.967 | 983 | 0.029 | 0.000 | 0.66 | 0.1 |
| `multiotsu` | 4,845 | 0.896 | 257 | 0.029 | 0.000 | 0.67 | 0.1 |
| `kde-valley` | 4,768 | 0.910 | 375 | 0.029 | 0.000 | 0.67 | 0.3 |
| `gmm` | 4,413 | 0.983 | 1176 | 0.029 | 0.000 | 0.66 | 1.2 |
| `gmm-2d` | 4,475 | 0.969 | 947 | 0.029 | 0.000 | 0.66 | 0.1 |

### PBMC 1k, v3 chemistry, Cell Ranger 3.0

![Calls on pbmc_1k_v3_3_0]({{ '/assets/comparison/pbmc_1k_v3_3_0.png' | relative_url }})

| Method | Cells | Jaccard with 10x | Smallest cell total | unspliced: cells | unspliced: others | unspliced: AUROC | mito: cells | mito: others | mito: AUROC | Seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| `cellranger@3.0.0:expect_cells=1000` | 1,222 | 1.000 | 503 | 0.316 | 0.143 | 0.89 | 0.103 | 0.077 | 0.69 | 11.8 |
| `cellranger` | 1,219 | 0.993 | 503 | 0.316 | 0.143 | 0.89 | 0.103 | 0.077 | 0.69 | 31.6 |
| `cellranger-ordmag` | 1,116 | 0.909 | 2213 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.67 | 2.3 |
| `starsolo-cellranger2.2` | 1,135 | 0.919 | 1868 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.9 |
| `starsolo-emptydrops-cr` | 1,226 | 0.993 | 503 | 0.316 | 0.143 | 0.89 | 0.103 | 0.077 | 0.69 | 4.5 |
| `starsolo-basic-awk` | 1,135 | 0.919 | 1868 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.9 |
| `dropletutils-barcoderanks` | 873 | 0.714 | 5204 | 0.311 | 0.150 | 0.91 | 0.098 | 0.077 | 0.65 | 0.8 |
| `dropletutils-emptydrops` | 1,166 | 0.915 | 101 | 0.316 | 0.143 | 0.89 | 0.101 | 0.077 | 0.68 | 14.4 |
| `dropletutils-emptydropscellranger` | 1,226 | 0.993 | 503 | 0.316 | 0.143 | 0.89 | 0.103 | 0.077 | 0.69 | 2.3 |
| `dropletutils-defaultdrops` | 1,135 | 0.919 | 1868 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.9 |
| `umitools-distance` | 1,089 | 0.891 | 2868 | 0.316 | 0.143 | 0.92 | 0.101 | 0.077 | 0.67 | 0.8 |
| `umitools-density` | 1,169 | 0.939 | 1182 | 0.318 | 0.143 | 0.92 | 0.102 | 0.077 | 0.68 | 1 |
| `kneedle` | 321 | 0.263 | 8717 | 0.311 | 0.154 | 0.89 | 0.092 | 0.077 | 0.62 | 0.9 |
| `cellsweep-empty` | 1,345 | 0.909 | 122 | 0.312 | 0.143 | 0.84 | 0.105 | 0.077 | 0.71 | 0.8 |
| `cellbender-prior` | 641 | 0.525 | 6487 | 0.307 | 0.154 | 0.90 | 0.096 | 0.077 | 0.64 | 1 |
| `alevinfry-permit-list` | 1,144 | 0.925 | 1659 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.8 |
| `bustools-allowlist` | 1,161 | 0.933 | 1285 | 0.318 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.9 |
| `alevin-initial-whitelist` | 1,169 | 0.939 | 1182 | 0.318 | 0.143 | 0.92 | 0.102 | 0.077 | 0.68 | 0.8 |
| `dropseq-dropbead` | 1,217 | 0.970 | 645 | 0.317 | 0.143 | 0.90 | 0.102 | 0.077 | 0.69 | 0.9 |
| `dropest-bounds` | 125,227 | 0.010 | 1 | 0.200 | — | — | 0.000 | — | — | 0.9 |
| `zumis-cell-barcodes` | 594 | 0.486 | 6772 | 0.307 | 0.154 | 0.89 | 0.096 | 0.077 | 0.64 | 0.8 |
| `mad` | 1,160 | 0.933 | 1320 | 0.318 | 0.143 | 0.92 | 0.102 | 0.077 | 0.68 | 0.8 |
| `otsu` | 1,143 | 0.924 | 1687 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.8 |
| `multiotsu` | 1,317 | 0.928 | 181 | 0.313 | 0.143 | 0.85 | 0.104 | 0.077 | 0.71 | 0.8 |
| `kde-valley` | 1,164 | 0.935 | 1255 | 0.318 | 0.143 | 0.92 | 0.102 | 0.077 | 0.68 | 0.8 |
| `gmm` | 1,121 | 0.913 | 2138 | 0.317 | 0.143 | 0.92 | 0.101 | 0.077 | 0.67 | 1.8 |
| `gmm-2d` | 1,151 | 0.925 | 1204 | 0.318 | 0.143 | 0.92 | 0.101 | 0.077 | 0.68 | 0.8 |

### PBMC 5k, GEM-X, Cell Ranger 9.0

![Calls on pbmc5k_gemx_9_0]({{ '/assets/comparison/pbmc5k_gemx_9_0.png' | relative_url }})

| Method | Cells | Jaccard with 10x | Smallest cell total | mito: cells | mito: others | mito: AUROC | Seconds |
|---|---|---|---|---|---|---|---|
| `cellranger@9.0.0` | 5,710 | 1.000 | 512 | 0.029 | 0.000 | 0.73 | 49 |
| `cellranger` | 5,710 | 1.000 | 512 | 0.029 | 0.000 | 0.73 | 46.9 |
| `cellranger-ordmag` | 5,605 | 0.982 | 3530 | 0.029 | 0.000 | 0.72 | 2.3 |
| `starsolo-cellranger2.2` | 5,589 | 0.979 | 3933 | 0.029 | 0.000 | 0.72 | 0.2 |
| `starsolo-emptydrops-cr` | 5,710 | 0.999 | 512 | 0.029 | 0.000 | 0.73 | 6.8 |
| `starsolo-basic-awk` | 5,589 | 0.979 | 3933 | 0.029 | 0.000 | 0.72 | 0.4 |
| `dropletutils-barcoderanks` | 4,775 | 0.836 | 8027 | 0.028 | 0.000 | 0.72 | 0.2 |
| `dropletutils-emptydrops` | 5,695 | 0.982 | 103 | 0.029 | 0.000 | 0.73 | 52.7 |
| `dropletutils-emptydropscellranger` | 5,710 | 0.999 | 512 | 0.029 | 0.000 | 0.73 | 2.8 |
| `dropletutils-defaultdrops` | 5,589 | 0.979 | 3933 | 0.029 | 0.000 | 0.72 | 0.9 |
| `umitools-distance` | 5,566 | 0.975 | 4327 | 0.029 | 0.000 | 0.72 | 0.8 |
| `umitools-density` | 5,640 | 0.987 | 2616 | 0.029 | 0.000 | 0.72 | 1 |
| `kneedle` | 1,299 | 0.227 | 14759 | 0.031 | 0.000 | 0.71 | 0.8 |
| `cellsweep-empty` | 5,836 | 0.978 | 84 | 0.029 | 0.000 | 0.73 | 0.7 |
| `cellbender-prior` | 2,868 | 0.502 | 11148 | 0.028 | 0.000 | 0.71 | 0.7 |
| `alevinfry-permit-list` | 5,640 | 0.987 | 2616 | 0.029 | 0.000 | 0.72 | 0.3 |
| `bustools-allowlist` | 5,672 | 0.992 | 1433 | 0.029 | 0.000 | 0.72 | 0.9 |
| `alevin-initial-whitelist` | 5,640 | 0.987 | 2616 | 0.029 | 0.000 | 0.72 | 1 |
| `dropseq-dropbead` | 5,535 | 0.969 | 4795 | 0.028 | 0.000 | 0.72 | 0.8 |
| `dropest-bounds` | 239,902 | 0.024 | 1 | 0.000 | — | — | 1 |
| `zumis-cell-barcodes` | 2,208 | 0.387 | 12288 | 0.029 | 0.000 | 0.71 | 0.7 |
| `mad` | 5,589 | 0.979 | 3933 | 0.029 | 0.000 | 0.72 | 0.4 |
| `otsu` | 5,654 | 0.989 | 2103 | 0.029 | 0.000 | 0.72 | 0.3 |
| `multiotsu` | 5,753 | 0.993 | 296 | 0.029 | 0.000 | 0.73 | 0.9 |
| `kde-valley` | 5,690 | 0.995 | 849 | 0.029 | 0.000 | 0.72 | 1 |
| `gmm` | 5,632 | 0.986 | 2948 | 0.029 | 0.000 | 0.72 | 0.9 |
| `gmm-2d` | 5,610 | 0.982 | 2541 | 0.029 | 0.000 | 0.72 | 0.9 |

### mouse spinal cord nuclei, Cell Ranger 10.0

![Calls on spinal_cord_10_0]({{ '/assets/comparison/spinal_cord_10_0.png' | relative_url }})

| Method | Cells | Jaccard with 10x | Smallest cell total | mito: cells | mito: others | mito: AUROC | Seconds |
|---|---|---|---|---|---|---|---|
| `cellranger@10.0.0` | 15,266 | 0.999 | 500 | 0.003 | 0.017 | 0.11 | 144.2 |
| `cellranger` | 15,266 | 0.999 | 500 | 0.003 | 0.017 | 0.11 | 86.3 |
| `cellranger-ordmag` | 5,924 | 0.388 | 3096 | 0.001 | 0.016 | 0.07 | 12.5 |
| `starsolo-cellranger2.2` | 5,228 | 0.343 | 3623 | 0.001 | 0.016 | 0.07 | 9.6 |
| `starsolo-emptydrops-cr` | 15,375 | 0.992 | 500 | 0.003 | 0.017 | 0.12 | 14.8 |
| `starsolo-basic-awk` | 5,228 | 0.343 | 3623 | 0.001 | 0.016 | 0.07 | 4.5 |
| `dropletutils-barcoderanks` | 2,290 | 0.150 | 8374 | 0.001 | 0.015 | 0.06 | 1.6 |
| `dropletutils-emptydrops` | 2,290 | 0.150 | 8374 | 0.001 | 0.015 | 0.06 | 55.2 |
| `dropletutils-emptydropscellranger` | 15,378 | 0.991 | 501 | 0.003 | 0.017 | 0.12 | 5.9 |
| `dropletutils-defaultdrops` | 5,228 | 0.343 | 3623 | 0.001 | 0.016 | 0.07 | 0.5 |
| `umitools-distance` | 3,554 | 0.233 | 5683 | 0.001 | 0.016 | 0.06 | 0.4 |
| `umitools-density` | 2,925 | 0.192 | 6910 | 0.001 | 0.016 | 0.06 | 4.4 |
| `kneedle` | 4,389 | 0.288 | 4442 | 0.001 | 0.016 | 0.07 | 0.4 |
| `cellsweep-empty` | 14,181 | 0.906 | 744 | 0.002 | 0.017 | 0.11 | 0.3 |
| `cellbender-prior` | 6,028 | 0.395 | 3033 | 0.001 | 0.016 | 0.07 | 0.7 |
| `alevinfry-permit-list` | 154,250 | 0.099 | 128 | 0.015 | 0.016 | 0.52 | 0.4 |
| `bustools-allowlist` | 13,571 | 0.877 | 843 | 0.002 | 0.017 | 0.10 | 0.7 |
| `alevin-initial-whitelist` | 1,140 | 0.075 | 12113 | 0.001 | 0.015 | 0.06 | 5.4 |
| `dropseq-dropbead` | 4,506 | 0.295 | 4326 | 0.001 | 0.016 | 0.07 | 0.4 |
| `dropest-bounds` | 181,638 | 0.084 | 3 | 0.014 | — | — | 0.6 |
| `zumis-cell-barcodes` | 15,092 | 0.930 | 615 | 0.003 | 0.017 | 0.13 | 0.4 |
| `mad` | 155,698 | 0.098 | 121 | 0.015 | 0.014 | 0.55 | 0.3 |
| `otsu` | 13,781 | 0.888 | 808 | 0.002 | 0.017 | 0.11 | 0.3 |
| `multiotsu` | 159,891 | 0.095 | 25 | 0.015 | 0.000 | 0.93 | 0.4 |
| `kde-valley` | 14,277 | 0.910 | 728 | 0.002 | 0.017 | 0.11 | 4 |
| `gmm` | 19,885 | 0.767 | 100 | 0.004 | 0.017 | 0.23 | 1.7 |
| `gmm-2d` | 18,284 | 0.834 | 130 | 0.003 | 0.017 | 0.20 | 0.6 |

### mouse lymph node, Flex (1 of a 4-plex), Cell Ranger 7.2

![Calls on flex_lymphnode1_7_2]({{ '/assets/comparison/flex_lymphnode1_7_2.png' | relative_url }})

| Method | Cells | Jaccard with 10x | Smallest cell total | mito: cells | mito: others | mito: AUROC | Seconds |
|---|---|---|---|---|---|---|---|
| `cellranger@7.2.0` | 7,586 | 0.999 | 696 | 0.013 | 0.000 | 0.61 | 11.4 |
| `cellranger` | 7,578 | 0.998 | 696 | 0.013 | 0.000 | 0.61 | 42.9 |
| `cellranger-ordmag` | 6,432 | 0.847 | 4043 | 0.013 | 0.000 | 0.61 | 1.1 |
| `starsolo-cellranger2.2` | 7,151 | 0.893 | 6851 | 0.014 | 0.000 | 0.62 | 0.2 |
| `starsolo-emptydrops-cr` | 8,369 | 0.896 | 503 | 0.014 | 0.000 | 0.62 | 9.5 |
| `starsolo-basic-awk` | 7,151 | 0.893 | 6851 | 0.014 | 0.000 | 0.62 | 0.3 |
| `dropletutils-barcoderanks` | 5,191 | 0.684 | 5409 | 0.013 | 0.000 | 0.60 | 0.2 |
| `dropletutils-emptydrops` | 7,269 | 0.943 | 101 | 0.013 | 0.000 | 0.61 | 33.4 |
| `dropletutils-emptydropscellranger` | 7,605 | 0.995 | 503 | 0.013 | 0.000 | 0.61 | 2.6 |
| `dropletutils-defaultdrops` | 6,100 | 0.804 | 3999 | 0.013 | 0.000 | 0.61 | 0.4 |
| `umitools-distance` | 6,551 | 0.863 | 2928 | 0.013 | 0.000 | 0.61 | 0.3 |
| `umitools-density` | 7,267 | 0.925 | 1559 | 0.013 | 0.000 | 0.61 | 0.5 |
| `kneedle` | 1,515 | 0.200 | 8980 | 0.013 | 0.000 | 0.59 | 0.2 |
| `cellsweep-empty` | 8,912 | 0.850 | 155 | 0.014 | 0.000 | 0.62 | 0.2 |
| `cellbender-prior` | 4,379 | 0.571 | 16876 | 0.013 | 0.000 | 0.61 | 0.4 |
| `alevinfry-permit-list` | 7,239 | 0.925 | 1602 | 0.013 | 0.000 | 0.61 | 0.2 |
| `bustools-allowlist` | 7,265 | 0.926 | 1561 | 0.013 | 0.000 | 0.61 | 0.4 |
| `alevin-initial-whitelist` | 7,267 | 0.925 | 1559 | 0.013 | 0.000 | 0.61 | 0.5 |
| `dropseq-dropbead` | 6,173 | 0.814 | 3832 | 0.013 | 0.000 | 0.61 | 0.2 |
| `dropest-bounds` | 18,426 | 0.411 | 33 | 0.014 | 0.000 | 0.59 | 0.3 |
| `zumis-cell-barcodes` | 4,085 | 0.539 | 6314 | 0.013 | 0.000 | 0.60 | 0.2 |
| `mad` | 7,650 | 0.923 | 1040 | 0.013 | 0.000 | 0.61 | 0.2 |
| `otsu` | 6,968 | 0.910 | 2096 | 0.013 | 0.000 | 0.61 | 0.2 |
| `multiotsu` | 8,663 | 0.875 | 286 | 0.013 | 0.000 | 0.62 | 0.2 |
| `kde-valley` | 7,347 | 0.928 | 1430 | 0.013 | 0.000 | 0.61 | 0.5 |
| `gmm` | 6,890 | 0.902 | 2260 | 0.013 | 0.000 | 0.61 | 0.2 |
| `gmm-2d` | 6,544 | 0.862 | 2549 | 0.013 | 0.000 | 0.61 | 0.2 |

<!-- comparison:end -->
