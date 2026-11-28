---
title: Ambimux
parent: Documented only
grand_parent: Methods
nav_order: 6
method_ids: [ambimux]
status: docs-only
---

# Ambimux: empty, singlet or doublet from donor genotypes
{: .no_toc }

Ambimux models each droplet of a channel that pools several genotyped donors as empty, a singlet
or a doublet, with a share of its molecules from the ambient pool. It uses only the alleles that
the droplet's RNA and ATAC reads carry at known variants. An EM fit gives each droplet posterior
probabilities of the three; a droplet not called empty holds a nucleus. It needs read-level BAMs
and the donors' genotypes, which a count matrix does not carry. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `ambimux` (documented only) |
|---|---|
| Tool | Ambimux v0.5.0 (the only tag, 2025-09-19; C, built from source) |
| Category | mixture |
| Counts | informative molecules: RNA UMIs and ATAC fragments over a variant, with their alleles; totals only for the fixed-empty floor |
| Needs expected cells | no |
| Deterministic | in principle: no random numbers, but with several threads sums are added in thread order, so the last bits may differ |
| Original licence | MIT |
| Extra inputs | coordinate-sorted RNA and/or ATAC BAMs with `CB` and `UB` tags; a VCF of the pooled donors' genotypes; a GTF and a peak BED |

## How it works

At the tag v0.5.0:

1. **Counting.** RNA molecules are deduplicated by barcode and UMI, ATAC fragments deduplicated.
   Per barcode, the molecules that overlap a bi-allelic SNP of the VCF are informative, each base
   with its Phred quality (at least `--min-phredq`, 30) and a mapping quality of at least 30. Reads
   outside genes and peaks leave the model by default (`--mdl-reads intra`).
2. **Fixed empties.** Barcodes with fewer than `--out-min` (100) RNA UMIs **and** fewer than 100
   ATAC fragments are fixed as empty, with an ambient fraction of 1, and are not tested. They inform
   the ambient pool.
3. **Model.** For droplet $$d$$, $$H_d \in \{0, 1, 2\}$$ (empty, singlet, doublet) with prior
   $$\lambda$$, and donors $$S_d$$ with proportions $$\pi$$; a doublet's pair has probability
   $$\pi_{s_1}\pi_{s_2} / \sum_{i \ne j} \pi_i \pi_j$$. Each informative molecule comes from the
   ambient pool with probability $$\alpha_d$$ (per modality and assignment) or from the nucleus,
   split equally between a doublet's two; an empty droplet's molecules are all ambient. A base's
   probability given its source is $$(1 - \varepsilon)\,p(\text{allele}) + \varepsilon/2$$, with
   $$\varepsilon$$ from its base quality (0.01 if missing). $$p(\text{allele})$$ is the donor's
   genotype (0, 0.5 or 1, from `GT`) or, for the ambient source, the donors' genotypes averaged
   with weights $$\pi_\text{amb}$$. Molecules without a variant add nothing: the likelihood is
   allelic only.
4. **Fit.** First $$\pi_\text{amb}$$, from the fixed empties, iterated to convergence (at most
   10,000 iterations). Then EM. Per droplet and assignment, $$\alpha$$ maximises the marginal
   likelihood (a Newton-type step) under a Beta-like prior centred on the mean singlet $$\alpha$$,
   clipped to $$[10^{-4}, 1 - 10^{-4}]$$, with weight `--amb-prior-w` ($$10^{-8}$$). An inner loop
   updates $$\lambda$$, $$\pi$$ (pseudocount $$10^{-8}$$) and the ambient genotype until the
   relative change is below `--eps` ($$10^{-6}$$), at most `--max-iter` (100) times. An outer loop
   runs until the mean singlet $$\alpha$$ moves by less than `--eps`, at most 100 times. The
   start: $$\alpha = 0.5$$; $$\lambda$$ the fixed-empty share for empty, and 0.9 and 0.1 of the
   tested share for singlet and doublet.
5. **Output.** Per tested droplet, $$\log P(H, S, X)$$ for every assignment; the posteriors `PP0`,
   `PP1` and `PP2` of empty, singlet and doublet; and `best_type`, the single assignment
   $$(H, S)$$ with the highest joint probability.
6. **Call.** The program sets no threshold. The README recommends calling singlets with `PP1`
   above 0.9, 0.95 or 0.99. A cell-versus-empty call is `best_type` not empty, or $$1 - \text{PP0}$$
   above a threshold of the user's choice.

A reading of the model, not measured: since $$\alpha = 1$$ reproduces the empty likelihood, a
droplet whose informative molecules match the ambient genotype mixture is decided by $$\lambda$$.
A droplet is called non-empty when its alleles lean towards one donor, or two. A droplet with few
informative molecules carries little evidence either way, whatever its total count.

There are no random numbers. Droplets are processed in a fixed hash-based permutation (`permute`
in `math_util.c`, unseeded). With several threads, sums are accumulated under a lock in the order
the threads arrive, so the last bits may differ between runs (not checked).

## Versions

| Ambimux | Change to the call |
|---|---|
| 0.2.0 (2023-11), 0.4.0 (2024-01) | version strings in `main.c`; no tags; not read |
| 0.5.0 (commit 832de39, 2024-09-17) | the version string becomes 0.5.0 |
| commits still named v0.5.0 (2024-10-28, 2025-08-20) | the estimate of the ambient donor proportions changed: "refactor ambient", then "Add robust pi_as estimate" and "scale vareps by num mol"; 376 changed lines in `mod.c` |
| tag v0.5.0 (11f55af, 2025-09-19) | the only tag; the model above |

So "ambimux v0.5.0" names more than one model: a build from before 20 August 2025 estimates the
ambient proportions differently from the tag. Only the tag can be reproduced.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--out-min` | 100 | RNA UMIs and ATAC fragments below which a barcode is fixed empty |
| `--mdl-reads` | `intra` | reads in the model: inside genes and peaks (`intra`), outside (`inter`), or `all` |
| `--max-nh` | 0 (no limit) | reads with more alignments than this, by the `NH` tag (`--rna-nh-tag`, `--atac-nh-tag`), are skipped; 1 keeps uniquely mapped reads only |
| `--min-phredq` | 30 | the lowest base quality at a variant |
| `--rna-mapq`, `--atac-mapq` | 30 | the lowest mapping quality |
| `--eps` | $$10^{-6}$$ | the EM's convergence: the inner loop, the outer loop and each droplet's $$\alpha$$ |
| `--amb-eps` | $$10^{-6}$$ | documented as the convergence threshold of the droplets' $$\alpha$$, but at the tag the model never reads it (`--eps` is used), so it does not change the call |
| `--max-iter` | 100 | EM iterations |
| `--amb-prior-w` | $$10^{-8}$$ | the weight of the prior on $$\alpha$$ |
| `--tx-basic` | off | keep only transcripts and exons tagged `basic` in the GTF |
| `--exclude` | none | a BED of regions whose ATAC reads are discarded |
| `--samples` | every sample in the VCF | a file of the VCF's sample IDs to demultiplex; without it the VCF's samples are not subset |
| threshold on `PP0` or `PP1` | none (README: 0.9, 0.95, 0.99 for singlets) | the user's call |

These options were read at the tag (commit 11f55af): their help in `main.c` and their use in
`clopts.c`, `bam_proc.c`, `gtf_anno.c`, `variants.c` and `mod.c`.

## Status

Documented only. Ambimux needs read-level BAMs and the pooled donors' genotypes, which a count
matrix does not have, and it answers another question: which donor, and how much ambient, in each
droplet of a multi-donor multiome channel. Its EM is deterministic in principle and MIT-licensed,
so it could be ported, but it would sit outside the package's model of a matrix in and a call out.

Where a pool has genotypes, its call can be set beside the package's. `summary.txt` gives, per
tested barcode, the counts, `best_type`, `best_sample`, the ambient fractions `best_rna_ambient`
and `best_atac_ambient`, `PP0` – `PP2` and `LLK0` – `LLK2`. Fixed empties are absent from it and
count as empty. Ambimux also writes gene, peak and allele count matrices over all barcodes (MEX).
Mark the barcodes whose `best_type` is not empty, or whose $$1 - \text{PP0}$$ passes the chosen
threshold, in a boolean `obs` column of the raw gene matrix read with `scd.read`, and
`scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

## References

- Alvarez et al. 2025, Integrated ambient modeling and genetic demultiplexing of single-cell
  RNA+ATAC multiome experiments with Ambimux, bioRxiv (not peer-reviewed),
  [doi:10.1101/2025.08.21.671671](https://doi.org/10.1101/2025.08.21.671671).
- Source: [marcalva/ambimux](https://github.com/marcalva/ambimux/tree/v0.5.0) at v0.5.0, commit
  11f55af (`mod.c`, `mod.h`, `math_util.c`, `main.c`, `clopts.c`, `bam_proc.c`, `gtf_anno.c`,
  `variants.c`, `README.md`, `LICENSE.txt`); the simulator
  [marcalva/ambisim](https://github.com/marcalva/ambisim).
- Research notes: docs/_research/docs-only.md §7.
