---
title: SADE
parent: Planned
grand_parent: Methods
nav_order: 15
method_ids: [sade]
status: planned
---

# SADE: the entropy of a barcode's Tn5 insertions
{: .no_toc }

SADE (Shannon entropy for scATAC-seq denoising), from the Streets lab. Each barcode's Tn5
insertions are counted in fixed windows of the genome; a model of open windows as Poisson and
closed ones as empty gives the barcode an entropy, highest for cells. A smoothing spline on the
ranked log entropies picks a steep drop as the cutoff, and barcodes whose signal covers too much
of the genome are removed as debris. The repository has no licence and no release, so this page
describes it in its own words and names commits. Planned for 0.3.

1. TOC
{:toc}

## At a glance

| Id | `sade` (planned) |
|---|---|
| Tool | SADE, no releases: commits `eef42ee` (2026-09-30, the preprint's date) and `f83ca20` (2026-11-24) |
| Category | knee |
| Counts | Tn5 insertions (two per fragment) per genome window, on chr1 by default |
| Needs expected cells | no (`-n` targeted cells, optional, only narrows the choice) |
| Deterministic | yes: no random numbers; threads change only the row order of a table |
| Original licence | none: no licence file, so all rights are reserved. Read, described here, never copied |
| Extra inputs | a fragments file, the genome's contig sizes and the ENCODE blacklist |

## How it works

As of `f83ca20`. The input is a 10x-style fragments file (rows with other than 5 fields are
dropped), one or more chromosomes (chr1 by default) and a window size $$w$$ (3,000 bp).

1. **Insertions per window.** Each fragment gives two insertions, at its start and at its end
   minus 1; its read count is ignored. A chromosome of length $$L$$ has
   $$M = \lfloor L / w \rfloor + 1$$ windows, and an insertion at $$x$$ falls in window
   $$\lfloor x / w \rfloor$$. Insertions in a window that overlaps the blacklist are set aside,
   but such windows still count in $$M$$. Several chromosomes are joined into one vector.
2. **Entropy per barcode.** Let $$n_j$$ be the counts of the $$M_+$$ non-empty windows,
   $$a_j = \lfloor (n_j + 2) / 2 \rfloor$$ (described as a paired-end adjustment), $$f_v$$ the
   number of windows with $$a_j = v$$ and $$\bar a$$ the mean of the $$a_j$$.
   - If $$\bar a \le 1$$ the barcode gets no entropy and is dropped.
   - The zero-truncated Poisson rate: $$\lambda = W_0(-\bar a e^{-\bar a}) + \bar a$$, with
     $$W_0$$ the principal branch of Lambert's W; $$P_0 = e^{-\lambda}$$.
   - Empty windows expected inside the open region: $$M_0 = M_+ P_0 / (1 - P_0)$$.
   - Open-state probability, adjusted for depth: $$p = (M_0 + M_+) / (M \bar a)$$; if $$p > 1$$
     the barcode is dropped. The reported open fraction is $$\min(p\,\bar a, 1)$$.
   - $$H_{states} = -p \log_2 p - (1 - p) \log_2 (1 - p)$$.
   - With $$\pi_v = \mathrm{Pois}(v; \lambda)(1 - P_0)$$ (so the code has it, at `f83ca20`,
     `calculate_entropy.py`; a zero-truncated probability would divide by $$1 - P_0$$ instead):
     $$H_{open} = -\sum_v f_v \pi_v \log_2 \pi_v - M_0 P_0 \log_2 P_0$$.
   - $$H = H_{states} + p\,H_{open}\,(1 - P_0) / (M_+ \bar a)$$.
3. **Ranked curve.** Entropies in decreasing order, $$y_r = \log_{10} H_{(r)}$$ for
   $$r = 0, 1, \dots$$, the first 20,000 kept. A cubic smoothing spline of $$y$$ on $$r$$
   (SciPy's `make_splrep`, smoothing $$s = 1$$, SciPy 1.15 or later) gives the first derivative
   $$d_1$$ and the second, scaled by its largest absolute value, $$d_2$$, at every rank.
4. **Candidates.** Scanning ranks upwards until $$k$$ are found (3 by default): a region opens at
   the first rank with $$d_2 < 0$$ and closes at the next with $$d_2 > 0$$; its candidate is the
   rank of least $$d_1$$ (the steepest fall) inside it. With $$k = 1$$ the candidate is the
   shoulder: the rank where the spline, both axes scaled to [0, 1], lies farthest above the chord
   from the first point to the last. With $$k > 1$$ the shoulder is used only when a region
   opened, never closed, and no region was found. With no candidate, the cutoff is the number of
   ranks with $$y$$ above its median minus 3 MAD, at most the last rank.
5. **Cutoff.** The candidate with the least $$d_1$$ gives the cutoff rank $$R$$. With `-n N`, a
   rank outside $$[\lfloor N/5 \rfloor, 20{,}000]$$ gives way to the next-steepest candidate,
   until one fits or none is left.
6. **Call.** The $$R$$ barcodes with the highest entropy; the barcode at rank $$R$$, whose entropy
   is reported as the cutoff, is not among them. Ties there are broken by a whole-line text sort.
   Of these, barcodes with an open fraction above `-s` (0.5; values outside [0.1, 1] revert to
   0.5) are flagged as debris and removed.

The result depends on SciPy's `make_splrep` (new in 1.15), the blacklist files shipped in the
repository and the chromosomes chosen.

## Versions

There are no releases; the default branch moves weekly (92 commits in November 2026), and the
entropy's normalisation alone changed about ten times between 2026-11-18 and 2026-11-21. Two
states matter:

| SADE | Change to the call |
|---|---|
| `eef42ee` (2026-09-30), the last commit before the preprint's submission | insertions at start and end; no blacklist; $$a_j = \lfloor (n_j + 1)/2 \rfloor$$; $$p = (M_0 + M_+)/M$$; $$\pi_v = (f_v / M_+)(1 - P_0)$$, the observed frequencies; $$H_{open} = -\sum_v \pi_v \log_2 \pi_v - P_0 \log_2 P_0$$; $$H = H_{states} + p\,H_{open}(1 - P_0)/M_+$$; $$k = 2$$; no shoulder or median fallback (no candidate: an error); debris when $$p > 0.5$$ |
| `f83ca20` (2026-11-24) | as on this page: blacklist on by default; end − 1; depth-adjusted entropy; $$k = 3$$ in the driver script; fallbacks; `-n`, `-p` and `autosomes` |

The preprint (posted 2026-10-23) describes $$k = 2$$ and an open-region entropy normalised by the
total number of insertions. That normalisation matches `f83ca20`'s $$M_+ \bar a$$ better than
`eef42ee`'s $$M_+$$, so the preprint may describe a later state, or a branch (inferred). The
README and the tutorial also say $$k = 2$$.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `-c` | `chr1` | chromosomes: a list, or `autosomes` |
| `-w` | 3,000 | window size, bp |
| `-k` | 3 (2 in the README and the preprint) | candidates sought |
| `-s` | 0.5 | open fraction above which a barcode is debris |
| `-b` | `Y` | leave out insertions in blacklisted windows |
| `-n` | 0 (off) | targeted cells, which bound the cutoff rank |
| `-p` | 4 | threads; change only the row order |
| `-g` | — | genome: hg38, mm10 or hg19, shipped with the repository |
| fixed | 20,000 ranks; cubic spline, $$s = 1$$; 3 MAD | steps 3–4 |

## Status

Planned for 0.3. Until the authors tag a release and add a licence, an implementation can only
be written from this page. Implementing it needs:

- a fragments-file reader, the same one as [Cell Ranger ATAC](cellranger-atac.md) needs, and the
  genome's contig sizes and blacklist;
- SciPy's `lambertw`, Poisson pmf and `make_splrep`, exactly;
- a pinned commit, expected to go stale.

The oracle runs a clone at the pinned commit in a Python environment with the SciPy of the day
and one thread. The repository's tutorial, on 10x's atac_pbmc_500_nextgem fragments (Cell Ranger
ATAC 2.0.0) with hg38 and $$k = 3$$, scores 54,046 barcodes, finds candidates at ranks 157, 534
and 995, and cuts at rank 534, entropy 0.01763 (the notebook's saved output, from `f83ca20` or a
commit shortly before). Any 10x ATAC or ARC fragments file serves; the method is deterministic,
so the comparison is exact. Blockers: no licence, no release, and a preprint that matches no
commit exactly.

## References

- Yang, An, Lee, Gehring, Booeshaghi, Huang & Streets 2026, Denoising single-cell ATAC-seq data
  with Shannon entropy, *bioRxiv*,
  [doi:10.64898/2026.10.20.760589](https://doi.org/10.64898/2026.10.20.760589).
- Source (no licence; read, never copied): [streetslab/SADE](https://github.com/streetslab/SADE/tree/f83ca20)
  at commits `f83ca20` and `eef42ee` (`sade.sh`, `find_threshold.py`, `autothreshold.py`); the
  preprint's analyses in [streetslab/SADE_manuscript](https://github.com/streetslab/SADE_manuscript).
- Research notes: docs/_research/planned-atac.md §5.
