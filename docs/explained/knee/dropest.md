---
title: dropEst's lagged slopes
parent: Knee geometry
grand_parent: Explained
nav_order: 8
explains: [dropest-bounds]
---

# dropEst's lagged slopes
{: .no_toc }

dropEst turns droplet reads into corrected counts, and its R package, dropestr, estimates how many
cells a library holds. `EstimateCellsNumber` gives an expected number and two bounds. The
barcodes within the lower bound are surely cells, those beyond the upper bound surely not, and
dropestr's quality classifier learns from the two groups. In the family's terms:

- it draws log count against log rank;
- instead of a derivative it takes a slope over a lag, a fixed number of ranks;
- the expected number of cells is the steepest of those slopes, the middle of the cliff.

The lag is 5 % of the barcodes it is given. That one choice makes the answer depend on what the
function is handed, rather than on the cells.

1. TOC
{:toc}

## The rule

With the barcodes sorted in decreasing order, $$u_1 \ge \dots \ge u_n$$, and the lag
$$\lambda = \operatorname{round}(0.05\,n)$$:

1. The lagged slope $$y_i = (\log u_{i+\lambda} - \log u_i) / (\log(i + \lambda) - \log i)$$: the
   slope of the chord across $$\lambda$$ ranks, placed at rank $$i + \lambda$$.
2. The lagged slope of that slope, $$y'$$, over another $$\lambda$$.
3. The **upper bound** is where the longest run of $$y' > 0$$ starts: the longest convex stretch.
4. The **expected number of cells** is the rank of the steepest $$y$$ before the upper bound, less
   $$\lambda/2$$, which moves the chord's position back to its middle.
5. The **lower bound** is 0.75 × the expected number.

## Why it lands where it does

A lagged slope is a derivative smoothed over $$\lambda$$ ranks: each $$y_i$$ is the chord across a
window, rather than a tangent. The steepest chord is the middle of the cliff, so the expected
number is an inflection, like [barcodeRanks'](barcoderanks.md). The difference is the window.
barcodeRanks smooths with a spline whose stiffness does not depend on the matrix's size. dropEst
smooths over 5 % of whatever barcodes it is given.

![Left: the lagged slope on PBMC 1k, from every barcode and from the barcodes with at least 20 UMIs. Right: dropEst's estimate against the barcodes it is given, on both samples]({{ '/assets/explained/dropest-lag.png' | relative_url }})

**From every barcode with a count**, PBMC 1k has 234,600 barcodes, and the lag is 11,730 ranks,
nine times the number of cells. A chord that starts among the cells spans the cliff and ten
thousand barcodes of shoulder, and is not steep. The steepest chord of 11,730 ranks is in the
error tail. It starts at rank 119,362, the last barcode with 2 UMIs, and falls to 1. Its middle,
the expected number, is rank 125,227.

**From the barcodes with at least 20 UMIs** there are 2,761 barcodes, and the lag is 138. Now the
steepest chord is the cliff, and the estimate is 1,318.

The right panel repeats this for every floor:

| barcodes given | PBMC 1k | spinal cord nuclei |
|---|---|---|
| every barcode with a count | 125,227 | 181,638 |
| at least 2 UMIs | 3,994 | 17,692 |
| at least 10 | 1,612 | 14,307 |
| at least 20 | 1,318 | 14,292 |
| at least 50 | 1,290 | 14,217 |
| at least 100 | 278 | 14,214 |
| at least 500 | 274 | 14,593 |
| 10x | 1,222 | 15,257 |

The rule works when the barcodes it is given hold the cliff and a stretch of what follows, but not
the long tail. Then the lag is a small fraction of the cliff's width, and the steepest chord is
the cliff:

- **On the nuclei** that holds for any floor from 5 to 500: the estimate stays between 13,936 and
  14,593. The ambient plateau keeps enough barcodes below the nuclei.
- **On PBMC 1k** it holds for floors from about 15 to 50. From 100 UMIs up, the barcodes left are
  almost only cells (1,367 of them), and nothing below the cliff remains. The longest convex
  stretch then starts at rank 360, inside the cells, so the search for the steepest slope stops
  there. The estimate is 278.

## Where it goes wrong

- **On a whole raw matrix**, where the lag covers tens of thousands of ranks: 125,227 on PBMC 1k
  and 181,638 on the nuclei. Both are features of the error tail. The
  [comparison on public data](../../overview/comparison.md) finds the same on every sample.
- **On a tightly filtered one**, where nothing below the cliff remains: 278 on PBMC 1k from 100
  UMIs up.
- **The bounds inherit the estimate.** The classes `high`, `unknown` and `low` that dropestr's
  [classifier](../../methods/dropest/classifier.md) trains on are cut from the same numbers.

## Details that matter

Each is on the [method page](../../methods/dropest/bounds.md).

- **The lag** is R's `round(0.05 * n)`, which rounds halves to even.
- **Natural logarithms**, although the base does not matter: it scales both axes alike, and every
  slope with them.
- **The longest run** is the first, if several are equally long; with no convex stretch at all,
  R's `which.max` falls back to the first run.
- **The upper bound inside the lag**, where `max − lag` is negative, makes R fail to index;
  the package stops with an error.

## Further reading

- Petukhov, Guo, Baryawno, Severe, Scadden, Samsonova & Kharchenko 2018, dropEst: pipeline for
  accurate estimation of molecular counts in droplet-based single-cell RNA-seq experiments,
  *Genome Biology* 19:78, [doi:10.1186/s13059-018-1449-6](https://doi.org/10.1186/s13059-018-1449-6).
- Source: [`dropestr/R/cell_number.R`](https://github.com/kharchenkolab/dropEst/blob/v0.8.6/dropestr/R/cell_number.R)
  (`EstimateCellsNumber`).
- The reference page: [dropEst cell-number bounds](../../methods/dropest/bounds.md).
