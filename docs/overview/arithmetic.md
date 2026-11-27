---
title: Reads, UMIs, ties and rounding
parent: Overview
nav_order: 2
---

# Reads, UMIs, ties and rounding
{: .no_toc }

A method is more than its formula. Its call also depends on:

- which number it counts per barcode;
- which barcodes it looks at;
- how it orders equal counts, and how it rounds;
- which random numbers it draws.

These details move a call by a barcode or by thousands, and a reimplementation that skips them
only resembles the original. This page collects them across the tools. Each method's page lists
its own under *Differences from the original*.

1. TOC
{:toc}

## What is counted

**UMIs or reads.** Most tools rank barcodes by UMIs: Cell Ranger, STARsolo, DropletUtils,
CellBender and cellsweep. UMI-tools, alevin, alevin-fry, bustools, zUMIs and Drop-seq count reads
by default. The two rank barcodes differently, because the number of reads per UMI varies from
barcode to barcode. A count matrix holds UMIs, so sc-droplets counts UMIs. `--reads` gives each
barcode's reads (an `obs` column, or a table of barcode and count), and `--unit reads` then
counts them, as the original would.

**Which features.** sc-droplets counts the Gene Expression features by default, as most tools do;
`--feature-type` chooses others. Cell Ranger's own runs count Gene Expression for the call but
look at antibodies for the aggregate filter; on a matrix with antibody or CRISPR features, a rule
that counts every feature gets other totals.

**Which genome.** In a run against two genomes, Cell Ranger calls each genome's cells on that
genome's features alone and takes the union. A rule run on both genomes' counts together finds a
different curve.

**Which matrix.** Exonic and intronic counts make different curves, most of all in nuclei.
STARsolo writes `Gene` (exonic) and `GeneFull` (exons and introns) side by side. Cell Ranger
counts intronic reads by default from 7.0. Velocyto's spliced, unspliced and ambiguous counts add
up to yet another total.

## Which barcodes enter

A rule that looks at the shape of the curve sees a different curve when its tail changes.

| Tool | The barcodes the rule sees |
|---|---|
| Cell Ranger, emptyDropsCellRanger | every barcode of the raw matrix; Cell Ranger 3.0's lists all 6,794,880 of the v3 list for PBMC 1k, most with no count |
| STARsolo | the barcodes with a count (its "detected" barcodes) |
| DropletUtils emptyDrops | every barcode; those up to `lower` (100) are the ambient ones, those above it are tested |
| DropletUtils barcodeRanks | runs of equal totals above `lower` (100) |
| UMI-tools | barcodes with a count, counted from the first 100 million reads |
| Kneedle in yascp | barcodes with at least 1,000 |
| zUMIs | barcodes with at least 100 reads |
| dropbead | the top 15,000 |
| dropEst | the barcodes dropEst kept; given every barcode of a raw 10x matrix, its estimate runs to 118,000 to 240,000 ([comparison](comparison.md#what-it-shows)) |
| sc-droplets' distribution rules | barcodes above `--min-counts` (100) |

Counting every barcode matters only where a rule reads ranks the barcodes with a count do not
reach. Cell Ranger's ambient droplets are ranks 45,000 to 90,000; with fewer barcodes holding a
count, zero-count barcodes become ambient.

## Ties

Equal counts are everywhere at the low end of the curve. On the spinal cord nuclei, 1,509
barcodes share the count at a boundary of Cell Ranger 10.0's ambient range. Two questions
decide where tied barcodes fall.

**How are equal counts ordered?**

| Tool | Order of equal counts |
|---|---|
| Cell Ranger, CellBender | numpy's default sort: before numpy 1.25 everywhere, and on processors without the vector instructions after it, its classic quicksort, which is not stable; with AVX-512 (or AVX2 in newer releases) a vectorised sort that leaves them elsewhere |
| UMI-tools | the order barcodes were first seen (`Counter.most_common`) |
| STARsolo | position in the barcode list |
| DropletUtils, dropbead, zUMIs | input order (R's `order` is stable) |
| alevin | the order of its hash map, then the C++ library's `std::sort` |

sc-droplets ports numpy's classic quicksort for Cell Ranger and CellBender
(`primitives.legacysort`). A Cell Ranger run on a processor with the vectorised sort can still
differ: on the spinal cord, Cell Ranger 10.0 here calls 15,266 cells and 10x's run 15,257, every
difference at the FDR cut.

**Is a tie at the cut in?** `≥` keeps every barcode tied with the cut: alevin-fry, STARsolo,
Cell Ranger's bootstrap cut. `>` drops them all: DropletUtils `defaultDrops`, UMI-tools' density
valley. A rank keeps some and not others, by the order above: Cell Ranger's OrdMag takes the top
$$n$$ in its sort order; UMI-tools' distance knee the first $$k + 1$$. barcodeRanks avoids the
question by collapsing each run of equal totals into one point at its mid-rank.

## Rounding

| Rounding | Where |
|---|---|
| halves away from zero | Python 2's `round` (Cell Ranger 2.x, 3.x), Rust's `round` (alevin-fry) |
| halves to even | Python 3 and numpy (Cell Ranger 7.2+), R's `round` (DropletUtils) |
| $$\lfloor x + 0.5 \rfloor$$ | STAR |
| truncation | STARsolo's minimum for candidates, $$\lfloor 0.01 \cdot \text{median} \rfloor$$ |

The rank of OrdMag's robust maximum is $$\operatorname{round}(0.01 N)$$. With $$N = 2{,}250$$ it is
22.5: 23 under Python 2, 22 under numpy. Cell Ranger 7.2 computes $$1 - 0.99$$ in floating point,
which is slightly above 0.01, so its product is 22.50000000000002 and rounds to 23 either way.

## Floating point

| Detail | Where it shows |
|---|---|
| the C library's `log` and `log10` against numpy's vectorised ones, which can differ in the last bit | DropletUtils' dispersion $$\alpha$$; barcodeRanks 1.30's knee is 5199.8422235753687 on Linux and 5199.8422235753533 on macOS, in R and here |
| R's `sum` and `cumsum` accumulate in long double on x86-64 | emptyDrops' $$\alpha$$ (its likelihood is flat enough that the rounding steers `optimize`), zUMIs' and dropbead's cumulative counts |
| $$10^y$$ correctly rounded, as R's `powl` gives it | R's knee for a run of total 1,852 is 1851.9999999999993, so those barcodes are above it |
| 32-bit floats | scikit-image's Otsu and multi-Otsu from 0.20; cellsweep's `log1p`; CellBender reading a float32 `.h5ad` |
| fused multiply-add | scikit-image on arm64 picks another multi-Otsu split in about 1 histogram in 300; R with OpenBLAS changes `smooth.spline`'s last bits |

sc-droplets calls the C library where the original does, and squares with `x * x` rather than
`x ** 2`.

## Random numbers

| Tool | Generator and seed |
|---|---|
| Cell Ranger 2.x – 8.x | numpy's legacy stream seeded 0: OrdMag's bootstraps, then (3.x) the simulations continue it; from 4.0 the simulations start a fresh stream for each test |
| Cell Ranger 9.0+ | PCG64 seeded 42, once per run, 100,000 simulations |
| STARsolo EmptyDrops_CR | `std::mt19937` seeded $$19760110\,(s + 1)$$ for simulation $$s$$, and libstdc++'s `discrete_distribution`; a STAR built with libc++ (clang, macOS) draws differently |
| DropletUtils | R's `set.seed`; from 1.4 a PCG32 stream per iteration, seeded from R's generator |
| scikit-learn's mixtures (here `gmm`) | `random_state`, 0 here |

With the generator, the seed and the order of draws reproduced, a Monte Carlo test gives the same
p-values: sc-droplets matches STAR and DropletUtils exactly. The same tool built against another
C++ library, or a release that changed its stream, moves only the borderline barcodes.

## Off by one

- Cell Ranger 2.x indexes its baseline from the end; for fewer than 50 expected cells that index
  is 0, which takes the smallest value. 7.1 fixed this, and 7.2 takes the $$(b + 1)$$-th largest.
- STARsolo's robust maximum is the 31st largest total (0-based index 30); DropletUtils'
  `emptyDropsCellRanger` takes the 30th.
- DropletUtils barcodeRanks 1.0 – 1.8 places its knee (left edge − 1) runs too early.
- alevin 0.11.x sums its cumulative curve one barcode off; 0.14.1 – 1.1.0 turns $$F$$ forced cells
  into $$F - 1$$; from 1.2.0 the expected-cells rule keeps one more.
- alevin-fry 0.12 – 0.17 sets Flex's knee one rank above the single-barcode knee.

Each is a variant of its method, chosen by the tool's release: see the
[methods catalogue](../methods/index.md).
