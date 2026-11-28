---
title: dropbead's diagonal
parent: Knee geometry
grand_parent: Explained
nav_order: 1
explains: [dropseq-dropbead]
---

# dropbead's diagonal
{: .no_toc }

Drop-seq's first papers chose the number of cells by eye. They plotted the cumulative fraction of
reads against the barcodes ranked by reads, and read off where the curve flattens. dropbead, the
Rajewsky lab's R package for Drop-seq data, turned that reading into a rule:
`estimateCellNumber`. It is the simplest member of the family, and the cleanest case of its
central identity. In the family's terms:

- it draws the cumulative fraction of the top 15,000 barcodes against rank / 15,000;
- it takes the point farthest from the diagonal of that square;
- it does nothing else: no refit, no smoothing.

By the [family's identity](index.md#on-the-cumulative-curve-a-chord-is-a-mean), the cells are
**the barcodes above the mean count of the top 15,000**.

1. TOC
{:toc}

## The rule

With the counts sorted in decreasing order and $$N$$ = `max.cells` = 15,000:

1. $$x_i = i / N$$, and the cumulative fraction $$s_i = (c_1 + \dots + c_i) / (c_1 + \dots + c_N)$$,
   for $$i = 1, \dots, N$$.
2. The distance of each point from the diagonal $$s = x$$: dropbead writes it as
   $$\sin(\arctan m)\,\lvert s_i - x_i \rvert$$, with a slope $$m$$ that mixes the two scales. The
   factor $$\sin(\arctan m)$$ is the same for every point, so it moves nothing. What counts is the
   height $$s_i - x_i$$.
3. The number of cells is the first $$i$$ with the greatest height.

![dropbead's cumulative curve of the top 15,000 barcodes of PBMC 1k; its height above the diagonal; and the counts against the window's mean]({{ '/assets/explained/knee-mean.png' | relative_url }})

## Why it lands where it does

Moving from barcode $$i - 1$$ to barcode $$i$$, the height changes by

$$(s_i - x_i) - (s_{i-1} - x_{i-1}) = \frac{c_i}{c_1 + \dots + c_N} - \frac{1}{N}$$

which is positive exactly when $$c_i$$ is above $$\bar c_N = (c_1 + \dots + c_N)/N$$, the mean of
the window. The counts are sorted, so the height climbs while the counts are above the mean and
falls after. The greatest height is at the last barcode above the mean. dropbead's rule is a
threshold:

$$\text{cells} = \#\{\, i \le N : c_i > \bar c_N \,\}$$

The script behind these pages checks this against dropbead's own function (ported, and validated
against dropbead in R) on both samples, at windows of 2,000, 10,000, 15,000 and 50,000. The two
agree every time.

On PBMC 1k the top 15,000 hold the 1,222 cells and about 13,800 empty droplets. The mean, 633
UMIs, lies inside the cliff, and dropbead calls 1,217. The diagonal runs from the origin, so the
mean includes the first barcode. (UMI-tools' chord starts at the first barcode and leaves it out;
the difference is a barcode or so.)

So the only real parameter is the window. Here is the answer as the window grows:

![dropbead's number of cells against the size of its window, on PBMC 1k and the spinal cord nuclei]({{ '/assets/explained/dropbead-window.png' | relative_url }})

The two samples behave in opposite ways:

- **PBMC 1k** settles. Once the window holds the cells and a few thousand empty droplets, the
  mean falls below the cliff and stays there: 1,134 cells at a window of 5,000, 1,217 at 15,000,
  1,365 at 100,000. Any window from 5,000 barcodes up gives between 1,134 and 1,365 cells.
- **The nuclei** never settle. Their bend is gradual, so every larger window lowers the mean a
  little and takes in more nuclei: 1,868 cells at 5,000, 4,506 at 15,000, 13,416 at 100,000.
  Past about 300,000 the window takes in the error tail, the mean drops below the ambient
  plateau, and the answer jumps to the second cliff.

Suppose a window holds $$K$$ cells averaging $$C$$ UMIs and $$N - K$$ empty droplets averaging
$$a$$. Its mean is $$(KC + (N - K)\,a)/N$$. When $$K$$ is a small share of $$N$$, that mean sits
far below $$C$$, at the cliff. When $$K$$ approaches $$N$$, it approaches $$C$$ itself, and the rule
keeps only the cells above the average cell. A window of 15,000 leaves a wide margin for a run of
a few thousand cells, and none for a run of fifteen thousand.

## Where it goes wrong

- **More cells than a small fraction of the window.** The nuclei have 15,257 cells by 10x's count,
  as many as the window. The window is almost all nuclei, its mean of 4,325 UMIs is an average
  nucleus, and dropbead calls the 4,506 above it.
- **Fewer barcodes than the window.** With fewer than 15,000 barcodes (a filtered matrix, or a
  small run) dropbead's sums are missing values and it returns no number. The package stops
  with an error instead.
- **The plot and the rule disagree.** dropbead's plot, `plotCumulativeFractionOfReads`, draws the
  top 10,000, so its knee is a different one: 1,179 cells on PBMC 1k, and 3,398 on the nuclei.

## Details that matter

Each is on the [method page](../../methods/dropseq/dropbead.md).

- **Reads, not UMIs.** Drop-seq counts reads per barcode (`BAMTagHistogram`). This package counts
  UMIs by default, or reads with `--unit reads`.
- **The first maximum wins** (R's `which.max`). When a count equals the mean exactly, the earlier
  barcode is the knee.
- **Ties** fall in the order of dropbead's input file; here, the matrix's order
  ([differences](../../methods/dropseq/dropbead.md#differences-from-the-original)).
- **Arithmetic.** R's `sum` and `cumsum` accumulate in long double, and the port does too.

## Further reading

- Alles, Karaiskos, Praktiknjo et al. 2017, Cell fixation and preservation for droplet-based
  single-cell transcriptomics, *BMC Biology* 15:44,
  [doi:10.1186/s12915-017-0383-5](https://doi.org/10.1186/s12915-017-0383-5).
- Macosko et al. 2015, *Cell* 161:1202–1214,
  [doi:10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002): the cumulative
  plot read by eye.
- Source: `estimateCellNumber` in
  [`R/Tools.R`](https://github.com/rajewsky-lab/dropbead/blob/73fd121/R/Tools.R).
- The reference page: [dropbead](../../methods/dropseq/dropbead.md), with its options and
  validation.
- [Same idea, different tools](../../overview/families.md#the-farthest-point-from-a-chord).
