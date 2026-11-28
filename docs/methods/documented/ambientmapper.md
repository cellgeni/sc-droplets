---
title: AmbientMapper
parent: Documented only
grand_parent: Methods
nav_order: 7
method_ids: [ambientmapper]
status: docs-only
---

# AmbientMapper: empty, singlet or doublet from reads mapped to several genomes
{: .no_toc }

AmbientMapper is for libraries that pool several species or genotypes, mapped competitively to
each of their reference genomes; it was designed for scifi-ATAC. Each read gets a likelihood for
each genome, and an ambient profile over genomes is learnt from the barcodes with the fewest
reads. Per barcode, an empty model (all reads ambient), a singlet and a doublet are fitted, and
BIC chooses between them. It needs reads mapped to every reference, which a count matrix does not
carry. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `ambientmapper` (documented only) |
|---|---|
| Tool | AmbientMapper v0.1.0 (the only tag, 2026-09-20; Python, GitHub `gomezcan/ambientmapper`), `ambientmapper run` |
| Category | mixture |
| Counts | reads per barcode and genome, weighted by each read's likelihood for each genome |
| Needs expected cells | no |
| Deterministic | yes: only the barcode subsamples are random, seeded with 42 |
| Original licence | MIT |
| Extra inputs | one coordinate-sorted BAM per reference genome, of the same reads, with a barcode tag (`CB` or `BC`), `AS`, `NM`, MAPQ and optionally `XA` |

## How it works

At the tag v0.1.0, `ambientmapper run` runs extract, filter, chunks, assign and genotyping:

1. **Reads.** Per genome, a table of reads; barcodes with fewer than `min_barcode_freq` (10)
   reads in that genome's table are dropped.
2. **Winners.** Per read, the best and second-best genome by fewer mismatches (`NM`), then a
   higher alignment score (`AS`), then a higher MAPQ. Empirical distributions of $$\Delta$$AS and
   $$\Delta$$MAPQ (best less second), by decile of the winner's score, are learnt on 50,000
   barcodes sampled with `random.Random(42)`. They give each read p-values; at `alpha` 0.05, with
   a rule on NM and a dominance filter across genomes, the read is a "winner" for its best genome
   or "ambiguous".
3. **Read likelihoods.** Per read and genome, a score
   $$w_\text{AS}\,\mathrm{AS} + w_\text{MAPQ}\,\mathrm{MAPQ} - w_\text{NM}\,\mathrm{NM}$$ with
   weights 1 (the code uses the raw values; the README says robust z-scores). A softmax over the
   genomes with $$\beta = 1$$, with an ambient constant $$10^{-3}$$ added to its denominator,
   gives $$L(r, g)$$. With `winner_only`, $$L(r, g) = 1$$ for the read's best genome alone; the
   configuration's default turns it on, but `run` leaves it off unless given
   `--genotyping-winner-only`. A barcode's expected counts are $$C(b, g) = \sum_r L(r, g)$$.
4. **Ambient profile.** $$\eta(g)$$ starts from the 2 % of barcodes with the smallest effective
   count (`eta_seed_quantile`): their summed $$L$$ by genome, normalised. Then, twice
   (`eta_iters`), the barcodes whose top genome holds at most 0.6 of their mass and whose top two
   genomes are within a ratio of 1.5 are kept if their Jensen–Shannon divergence from $$\eta$$ is
   at or below the 95th percentile (`empty_tau_quantile`), and $$\eta$$ is recomputed from them.
   The README's "barcodes with fewer than 200 reads" does not match this code.
5. **Models.** Per barcode, with $$E_r = \sum_g \eta(g) L(r, g)$$, $$n$$ its reads and the top 3
   genomes by $$C(b, g)$$ as candidates:
   - empty: $$\ell_0 = \sum_r \log E_r$$, no parameter;
   - singlet $$g_1$$: $$\ell_1 = \max_\alpha \sum_r \log[(1 - \alpha) L(r, g_1) + \alpha E_r]$$,
     with $$\alpha$$ on 0, 0.02, …, 0.5; one parameter;
   - doublet $$g_1, g_2$$:
     $$\ell_2 = \max_{\alpha, \rho} \sum_r \log[(1 - \alpha)(\rho L(r, g_1) + (1 - \rho) L(r, g_2)) + \alpha E_r]$$,
     with $$\rho$$ on 0.1, 0.15, …, 0.9; two parameters;
   - each scored by $$\mathrm{BIC} = -2\ell + k \log n$$, $$k$$ its parameters.
6. **Call.** With $$p_1 \ge p_2$$ the two largest shares of $$C(b, \cdot)$$:
   - fewer than `min_reads` (5) reads: `low_reads`;
   - **empty** if $$\mathrm{BIC}_0 + 10 < \min(\mathrm{BIC}_1, \mathrm{BIC}_2)$$
     (`empty_bic_margin`) **and** $$p_1 \le 0.6$$ **and** $$p_1 / p_2 \le 1.5$$; optional gates
     on the reads and on the divergence from $$\eta$$ are off by default;
   - otherwise a singlet if $$\mathrm{BIC}_1 + 6 < \mathrm{BIC}_2$$ (`bic_margin`), a doublet if
     $$\mathrm{BIC}_2 + 6 < \mathrm{BIC}_1$$, else `ambiguous`;
   - a singlet is `single_clean` if $$p_1 \ge 0.6$$ and $$p_1 / p_2 \ge 2$$, else
     `dirty_singlet`; a doublet is `doublet` if $$p_2 \ge 0.2$$, else `weak_doublet`.

A cross-mapping correction (`xmap`) replaces $$L(r, g_1)$$ by a learnt mixture when a genome's
cross-mapping share exceeds 0.05. The configuration enables it by default (`--no-xmap` turns it
off); whether `run` supplies the estimates it needs was not checked. The README also lists a call
`indistinguishable`, which `genotyping.py` at the tag does not produce.

A reading of the rules, not measured: beyond the floor of 5 reads, the call looks at which genomes
a barcode's reads come from, not at how many there are. An empty barcode is one whose reads are
spread across genomes as the ambient profile is; a barcode with few reads, all from one genome,
is a singlet.

The only randomness is the barcode subsamples, seeded with 42; the grids and thresholds are fixed.

## Versions

| AmbientMapper | Change to the call |
|---|---|
| v0.1.0 (commit 7e77faf, 2026-09-20) | the only tag and the version the preprint used; the method above |
| 0.2 series (unreleased) | intermediate files move from TSV to Parquet; the CHANGELOG claims "no behavior change" in the calls |

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `min_barcode_freq` | 10 | the reads a barcode needs in a genome's table |
| `alpha` (assign) | 0.05 | the p-value cut that makes a read a winner |
| `beta` | 1 | the softmax's inverse temperature |
| `w_as`, `w_mapq`, `w_nm` | 1, 1, 1 | weights of AS, MAPQ and NM in a read's score |
| `ambient_const` | $$10^{-3}$$ | the ambient term in the softmax's denominator |
| `winner_only` | on in the configuration, off in `run` | each read counts for its best genome alone |
| `min_reads` | 5 | fewer reads: `low_reads` |
| `empty_bic_margin` | 10 | the BIC lead the empty model needs |
| `empty_top1_max`, `empty_ratio12_max` | 0.6, 1.5 | how evenly an empty barcode's reads must spread |
| `empty_reads_max`, `empty_jsd_max` | off | optional gates on an empty call |
| `bic_margin` | 6 | the BIC lead between singlet and doublet |
| `single_mass_min`, `ratio_top1_top2_min` | 0.6, 2.0 | a clean singlet's top share and ratio |
| `doublet_minor_min` | 0.2 | a doublet's second share |
| `eta_seed_quantile`, `eta_iters`, `empty_tau_quantile` | 0.02, 2, 0.95 | the ambient profile's seed and refinement |
| `alpha_grid`, `max_alpha`, `rho_grid` | 0.02, 0.5, 0.05 | the grids of $$\alpha$$ and $$\rho$$ |
| `topk_genomes` | 3 | candidate genomes per barcode |
| `xmap_*` | enabled; threshold 0.05 | the cross-mapping correction |

## Status

Documented only. AmbientMapper works on reads mapped to several references, not on a count
matrix, and its call needs every read's likelihood for each genome. It is deterministic and
MIT-licensed, so it could be reproduced, but only from the BAMs. scifi-demux, the scifi-ATAC
preprocessing it is published with, calls no cells: it extracts UMIs, trims, demultiplexes, maps,
marks duplicates and tags reads, with no knee, threshold or model on counts per barcode. It has no
page.

AmbientMapper's call can be set beside the package's. `final/<sample>_cells_calls.tsv.gz` gives,
per barcode, the `call`, `genome_1`, `genome_2`, `alpha`, `rho`, the BICs, `n_reads`,
`n_effective`, `p_top1` – `p_top3`, `entropy` and `status_flag`;
`final/<sample>_expected_counts_by_genome.csv` holds $$C(b, g)$$. Barcodes removed by
`min_barcode_freq` do not appear. Join the calls by barcode to a raw barcode × peak (or gene)
matrix of the same library read with `scd.read`, mark every call but `empty` and `low_reads` in a
boolean `obs` column (absent barcodes as not cells), and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

## References

- Gomez-Cano, Jiang, Welch & Marand 2026, Empirical estimation of ambient contamination in
  combinatorial single-cell methods using multi-reference mapping, bioRxiv (not peer-reviewed),
  [doi:10.64898/2026.09.11.750809](https://doi.org/10.64898/2026.09.11.750809).
- Source: [gomezcan/ambientmapper](https://github.com/gomezcan/ambientmapper/tree/v0.1.0) at
  v0.1.0, commit 7e77faf (`src/ambientmapper/genotyping.py`, `assign_streaming.py`, `pipeline.py`,
  `cli.py`, `README.md`, `CHANGELOG.md`, `LICENSE`);
  [gomezcan/scifi-demux](https://github.com/gomezcan/scifi-demux) at v0.1.3 and v0.1.4.
- Research notes: docs/_research/docs-only.md §8.
