---
title: alevin's tangent
parent: Knee geometry
grand_parent: Explained
nav_order: 9
explains: [alevin-initial-whitelist]
---

# alevin's tangent
{: .no_toc }

alevin, salmon's single-cell mode from 0.11 to 1.10, first chose an *initial whitelist* of
barcodes from their counts alone, then refined it with a classifier. Its knee is not a chord at
all. It is a tangent drawn from the origin to an unusual curve: the logarithm of the cumulative
counts, summed from the smallest barcode up. The geometry looks unlike the rest of the family, but
it comes down to the same kind of rule, a count above a scaled mean:

> a barcode is a cell when its count exceeds the mean count of the barcodes below it, times the
> logarithm of their sum.

1. TOC
{:toc}

## The rule

1. Take the top $$M$$ barcodes (`--maxNumBarcodes`, 100,000), and sum them from the smallest up:
   $$C_i$$ is the sum of the $$i + 1$$ smallest. The curve is $$c_i = \ln C_i$$.
2. Draw lines from the origin. The knee is the first point $$j$$ such that no later point lies
   below the line through the origin and $$(j, c_j)$$: the line from the origin that touches the
   curve with all of the curve's remainder on or above it.
3. The cells are the barcodes above that point, counted from the top: about $$M - j$$.
4. A correction follows. alevin finds the valleys in a kernel density of $$\log_{10}$$ counts, and
   may lower the number to one of them, but never raise it.

## Why it lands where it does

![The tangent from the origin to alevin's curve on PBMC 1k, and the counts against the threshold the tangent amounts to, on PBMC 1k and the nuclei]({{ '/assets/explained/alevin-tangent.png' | relative_url }})

Read the curve from left to right: each step adds one more barcode, smallest first, so the cells
come last and make the curve shoot up at the right-hand end. The line from the origin through
$$(j, c_j)$$ has slope $$c_j / j$$. The first $$j$$ with no later point below its line is where
$$c_i / i$$ is least. For a smooth curve, that is where the curve's slope equals the line's:

$$c'(j) = \frac{c_j}{j}$$

This is the parallel-tangent condition of the [family page](index.md#the-chord-and-the-parallel-tangent),
with the origin as one end of the chord. Now translate it back into counts. Adding the barcode with
count $$f$$ raises $$\ln C$$ by about $$f / C$$, so $$c'(j) \approx f / C_j$$, and the condition
becomes

$$f \approx \frac{C_j}{j} \cdot \ln C_j$$

**the mean count of the barcodes below, times the log of their sum.** The tangent is a threshold
in disguise. The condition holds only approximately in the discrete case, but closely, as the script
behind these pages finds:

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| mean count of the barcodes below the knee | 7.6 | 672.1 |
| ln of their sum | 13.52 | 18.01 |
| the product | 102.3 | 12,106 |
| the count at alevin's knee | 105 | 12,113 |
| cells | 1,362 | 1,140 |
| 10x's cells | 1,222 | 15,257 |

In the right-hand panels the counts cross that threshold, the dashed line, exactly at alevin's
knee. The factor $$\ln C$$ is 13 to 18, so a barcode must hold well over ten times the average of
everything below it:

- **On PBMC 1k** the barcodes below are mostly empty droplets, averaging 7.6 UMIs, and the bar
  is about 100 UMIs: the foot of the cliff, 1,362 cells.
- **On the nuclei** the barcodes below the knee include most of the nuclei, averaging 672 UMIs.
  The bar is about 12,100 UMIs, and only the 1,140 largest nuclei clear it.

The logarithm has a side effect the chord rules do not share. A tangent from the origin survives
any rescaling of the axes, since a linear map fixes the origin. It does not survive a shift.
Counting reads instead of UMIs multiplies every $$C$$ by roughly a constant, the reads per UMI.
That adds a constant to $$\ln C$$, shifts the curve, and moves the tangent.

## Where it goes wrong

- **The nuclei**: 1,140 cells, the lowest of any method in the
  [comparison on public data](../../overview/comparison.md).
- **The correction can only lower the knee.** On PBMC 1k the density's valley gives 1,169, below
  the knee's 1,362, and replaces it: the final call is 1,169. On the nuclei the valley gives 2,924,
  above the knee, so the knee stands. That correction is a density rule, for that family's pages
  ([Same idea, different tools](../../overview/families.md#a-valley-in-the-density-of-log-counts)).
- **The window is the top 100,000.** `--maxNumBarcodes` sets $$M$$, and with it which barcodes
  count as "below".

## Details that matter

Each is on the [method page](../../methods/alevin/initial-whitelist.md).

- **0.11's off-by-one.** The sum is shifted by one barcode, which moves the knee by one: 1,361
  instead of 1,362 on PBMC 1k, 1,139 instead of 1,140 on the nuclei
  ([versions](../../methods/alevin/initial-whitelist.md#versions)).
- **Ties.** From the knee's rank alevin steps back past the barcodes tied with it, stopping at
  the first barcode above.
- **Reads, not UMIs.** alevin counts reads per barcode before barcode correction. This package
  counts UMIs by default, or reads with `--unit reads`.
- **Other rules.** `--expectCells` is Cell Ranger 2's robust maximum divided by ten, and
  `--forceCells` and `--keepCBFraction` fix the number. The whitelist also adds a low-confidence
  band of half as many barcodes again (between 200 and 1,000) for the classifier.
- **Historical.** salmon 1.11 removed alevin. Its rules are kept for the data processed with it.

## Further reading

- Srivastava, Malik, Smith, Sudbery & Patro 2019, Alevin efficiently estimates accurate gene
  abundances from dscRNA-seq data, *Genome Biology* 20:65,
  [doi:10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y).
- Source: `getLeftBoundary` in
  [`Alevin.cpp`](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/Alevin.cpp) (1.10.3).
- The reference page: [alevin initial whitelist](../../methods/alevin/initial-whitelist.md).
- The classifier that follows: [alevin's naive Bayes whitelist](../../methods/alevin/naive-bayes.md),
  planned.
