---
title: alevin naive Bayes whitelist
parent: alevin
grand_parent: Methods
nav_order: 2
method_ids: [alevin-naive-bayes]
status: planned
---

# alevin naive Bayes whitelist: the top of the list against the band below
{: .no_toc }

salmon alevin's final whitelist (`performWhitelisting`, in `WhiteList.cpp`). After quantifying,
alevin trains a Gaussian naive Bayes classifier on its [initial whitelist](initial-whitelist.md):
the top half as cells, the low-confidence band below it as not cells, on features of each barcode
such as its mapping rate, deduplication rate and genes detected. The barcodes in between stay
when the classifier favours the cells. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `alevin-naive-bayes` (planned) |
|---|---|
| Tool | salmon 0.11.0 – 1.10.3 (`salmon alevin`, removed in 1.11.0) |
| Category | classifier |
| Counts | reads per barcode, raw and mapped, and alevin's per-gene estimates |
| Needs expected cells | no (optional, through the initial whitelist's `--expectCells`) |
| Deterministic | the classifier, yes. With more than one worker thread the feature rows can follow the threads' order of completion, so the calls near the class boundaries could change between runs (read in the source, not seen in a run) |
| Original licence | GPL-3 |
| Extra inputs | raw reads (unmapped included) and mapped reads per barcode; alevin's own gene estimates; optionally `--mrna` and `--rrna` gene lists |

## How it works

Releases 1.1.0 – 1.10.3. The input is the initial whitelist (the knee, or expected or forced
cells: the top $$k$$ by corrected read count) followed by the low-confidence band below it, of
width $$\lfloor 0.5k \rfloor$$ clamped to [`--lowRegionMinNumBarcodes` (200), 1,000], in
decreasing order of count: $$N$$ barcodes, $$L$$ of them in the band. Barcodes in the list with no
mapped read are dropped first, and $$L$$ shrinks for those at an index above $$N - L$$; one at
exactly $$N - L$$ is not counted (an off-by-one).

The classifier is skipped, and every barcode in the list kept, with `--whitelist` or
`--noWhitelist`, in the hash and naive equivalence-class modes, or when $$L < 200$$ or
$$N - L < 90$$.

1. **Features** (`CollapsedCellOptimizer`). With $$r$$ the barcode's corrected read count, $$m$$
   its mapped reads, alevin's per-gene estimates $$\alpha_g$$, $$u = \sum_g \alpha_g$$ its
   deduplicated UMIs and $$G = \#\{g : \alpha_g > 0\}$$ its expressed genes:
   - mapping rate $$m/r$$;
   - deduplication rate $$1 - u/m$$;
   - mean by max, $$(u/G)/\max_g \alpha_g$$;
   - $$G$$;
   - the number of genes with $$\alpha_g > u/G$$;
   - when the `--mrna` or `--rrna` list names more than one gene, the mitochondrial or ribosomal
     fraction $$\sum_{g \in M} \alpha_g / u$$.

   Each is printed to `featureDump.txt` with 6 significant digits and parsed back with `strtod`.
   With exactly one gene in a list, the classifier expects the column but the writer leaves it
   out, and alevin exits with code 84.
2. **Classes**, by position: true, the first $$T = \lfloor (N - L)/2 \rfloor$$ barcodes (the top
   half of the initial whitelist); ambiguous, the next $$N - L - T$$; false, the last $$L$$ (the
   band).
3. **Fit**, with two off-by-ones. The true class uses rows $$0, \dots, T - 2$$ and the false class
   rows $$N - L, \dots, N - 2$$, each dropping its last barcode, yet both divide by the full class
   size $$n_c$$: per feature, $$\mu = \sum x / n_c$$ and $$\sigma^2 = \sum (x - \mu)^2 / n_c$$
   over the rows used. The priors are integer divisions and come out 0, and a zero prior adds
   nothing to the log score: the classes are compared without priors.
4. **Score** per class:
   $$-\tfrac12\Big(0.5 + \sum_{j:\,\sigma_j^2 \ne 0}\big[\log(2\pi\sigma_j^2) + (x_j - \mu_j)^2/\sigma_j^2\big]\Big)$$.
   A feature with zero variance in a class is skipped in that class only.
5. **Whitelist.** The $$T$$ true barcodes and every ambiguous barcode whose true score is at least
   its false score (ties go to true). The band is never kept. `whitelist.txt` lists them in list
   order; `--dumpFeatures` also writes `predictions.txt`, the two scores of the kept ambiguous
   barcodes.

The classifier draws no random numbers. Upstream, the order in which alevin's hash map lists tied
barcodes sets the list's order and so the class boundaries, and alevin's EM sets the features.
One more hazard is read in the source but not seen in a run: the worker threads write
`featureDump.txt` in the order they finish (under a mutex), while `performWhitelisting` ignores
its barcode column and takes row $$i$$ as the list's barcode $$i$$. With one worker (`-p` at most
2) the orders agree; with more, neighbouring rows can swap, and a barcode can be classified on
its neighbour's features.

## Versions

| salmon | Change to the call |
|---|---|
| 0.11.0 – 0.13.1 | features computed in memory, not rounded: mapping rate, deduplication rate, mean by max and genes over the mean with the mean over **all** genes of the reference (4), plus mitochondrial and ribosomal fractions, plus each barcode's largest correlation with the true set (`--useCorrelation`, off by default). 0.11.3 only rewrites the prediction loop |
| 0.14.0 – 1.0.x | features read back from `featureDump.txt` (6 significant digits), the mean over **expressed** genes; expressed genes added, and the mean number of molecules per UMI-graph arborescence (6), plus mitochondrial and ribosomal fractions |
| 1.1.0 – 1.10.3 | the arborescence feature dropped (5), plus mitochondrial and ribosomal fractions |
| 1.11.0 (Mar 2026) | `salmon alevin` removed; 1.10.3 (Mar 2024) is the last release with it |

The fit and the scores are the same code from 0.11.0 to 1.10.3. The list the classifier sees also
moves with the [initial whitelist's versions](initial-whitelist.md#versions).

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--mrna FILE`, `--rrna FILE` | none | gene lists that add the mitochondrial and ribosomal fractions |
| `--lowRegionMinNumBarcodes` | 200 | the band's least width; a band under 200 skips the classifier |
| `--noWhitelist` | off | skip the classifier |
| `--whitelist FILE` | none | a list given; no classifier |
| `--useCorrelation` | off | add the correlation feature (0.11 – 0.13) |
| `--dumpFeatures` | off | also write `predictions.txt` |
| `-p` | salmon's default | threads; above 2, more than one worker, which can reorder the feature rows |
| `--expectCells`, `--forceCells`, `--keepCBFraction`, `--maxNumBarcodes`, `--freqThreshold` | see the [initial whitelist](initial-whitelist.md#parameters) | the list the classifier sees |

## Status

Planned for 0.2. Implementing it needs:

- the initial whitelist and its band: [already implemented](initial-whitelist.md);
- the classifier: small, but the 6-significant-digit round trip and both off-by-ones must be
  reproduced;
- the features. From a count matrix alone only mean by max, expressed genes, genes over the mean
  and the mitochondrial and ribosomal fractions can be computed, and only from alevin's own
  estimates are they alevin's. The mapping and deduplication rates need raw (unmapped included)
  and mapped reads per barcode: STARsolo's `CellReads.stats` or Cell Ranger's
  `molecule_info.h5` approximate them, not exactly.

The oracle has two levels. The classifier alone: a program built from `WhiteList.cpp`'s
functions copied out of each tag, as `oracles/alevin/harness.py` does for `Alevin.cpp`, fed a
feature table, the whitelist compared exactly. End to end: salmon alevin 1.10.3 from bioconda on
FASTQs, with `-p` at most 2 so that the feature rows keep the list's order.

## References

- Srivastava, Malik, Smith, Sudbery & Patro 2019, Alevin efficiently estimates accurate gene
  abundances from dscRNA-seq data, *Genome Biology* 20:65,
  [doi:10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y).
- Source: [`WhiteList.cpp` at 1.10.3](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/WhiteList.cpp)
  (`performWhitelisting`, `naive_bayes_fit`, `get_log_likelihood`, `classifyBarcodes`) and
  [`CollapsedCellOptimizer.cpp`](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/CollapsedCellOptimizer.cpp),
  with `WhiteList.cpp` at [1.0.0](https://github.com/COMBINE-lab/salmon/blob/v1.0.0/src/WhiteList.cpp)
  and [0.11.0](https://github.com/COMBINE-lab/salmon/blob/v0.11.0/src/WhiteList.cpp).
- Research notes: `docs/_research/planned-rna.md` §9.
