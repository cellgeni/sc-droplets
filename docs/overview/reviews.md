---
title: Reviews and benchmarks
parent: Overview
nav_order: 1
---

# Reviews and benchmarks
{: .no_toc }

1. TOC
{:toc}

We found no independent review or benchmark devoted to empty-droplet detection (searched
November 2026). What exists are comparisons inside papers that propose a method, run by its
authors, and broader papers on preprocessing or quality control that touch cell calling in
passing. This page lists them so their results can be read with that in mind.

## Comparisons inside method papers

| Paper | Compared with | Data and ground truth | Finding |
|---|---|---|---|
| Kavaliauskaite & Madsen 2023, valiDrops, *NAR Genom Bioinform* 5:lqad101, [doi:10.1093/nargab/lqad101](https://doi.org/10.1093/nargab/lqad101) | EmptyDrops, scCB2, DIEM, EmptyNN, CellBender, dropkick | 47 samples from 5 datasets, single cells and nuclei; 9 filtering and 12 integration metrics | valiDrops ranked best in 35 of 47 samples. The broadest comparison there is. |
| Yan, Zhao & Simon 2021, EmptyNN, *Patterns* 2:100311, [doi:10.1016/j.patter.2021.100311](https://doi.org/10.1016/j.patter.2021.100311) | Cell Ranger 2, Cell Ranger 3, EmptyDrops, CellBender | cell hashing and genetic demultiplexing as ground truth | AUROC: EmptyNN 94.7 %, Cell Ranger 3 90.2 %, CellBender 88.8 %, Cell Ranger 2 86.9 %, EmptyDrops 79.9 % |
| Ni et al. 2020, CB2, *Genome Biol* 21:137, [doi:10.1186/s13059-020-02054-8](https://doi.org/10.1186/s13059-020-02054-8) | EmptyDrops | simulations built from ten 10x datasets | CB2 found 24 % more cells on average (4–81 %) at the same FDR |
| Alvarez et al. 2020, DIEM, *Sci Rep* 10:11019, [doi:10.1038/s41598-020-67513-5](https://doi.org/10.1038/s41598-020-67513-5) | a Cell Ranger-style quantile cut-off, EmptyDrops | single-nucleus adipose, preadipocyte and brain samples, PBMC 68k; spliced-read fraction as a proxy for contamination | DIEM better on nuclei, similar on PBMC |
| Heiser et al. 2021, dropkick, *Genome Res* 31:1742, [doi:10.1101/gr.271908.120](https://doi.org/10.1101/gr.271908.120) | thresholds, EmptyDrops, CellBender (runtime) | simulated and real data | better recovery of rare types and removal of empty and noisy barcodes |
| Xi et al. 2023, SiftCell, *Cell Syst* 14:620, [doi:10.1016/j.cels.2023.06.002](https://doi.org/10.1016/j.cels.2023.06.002) | Cell Ranger cut-off, EmptyDrops, DIEM, CellBender | PBMC v3, E18 brain nuclei, a Drop-seq mixture; a shuffled null as silver standard | |
| Muskovic & Powell 2021, DropletQC, *Genome Biol* 22:329, [doi:10.1186/s13059-021-02547-0](https://doi.org/10.1186/s13059-021-02547-0) | Cell Ranger, EmptyDrops calls | nuclear fraction | their calls still hold empty droplets with much ambient RNA, and damaged cells |
| Megas, Lorenzi & Marioni 2024, EmptyDropsMultiome, *Genome Biol* 25:121, [doi:10.1186/s13059-024-03259-x](https://doi.org/10.1186/s13059-024-03259-x) | Cell Ranger ARC, EmptyDrops | multiome | Cell Ranger ARC misses more than half the nuclei in some samples |
| Petukhov et al. 2018, dropEst, *Genome Biol* 19:78, [doi:10.1186/s13059-018-1449-6](https://doi.org/10.1186/s13059-018-1449-6) | count thresholds | | |
| Srivastava et al. 2019, alevin, *Genome Biol* 20:65, [doi:10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y) | Cell Ranger | | |
| He et al. 2022, alevin-fry, *Nat Methods* 19:316, [doi:10.1038/s41592-022-01408-3](https://doi.org/10.1038/s41592-022-01408-3) | permit-list strategies | | |

## Other studies that touch cell calling

- Brüning et al. 2022, *GigaScience* 11:giac001,
  [doi:10.1093/gigascience/giac001](https://doi.org/10.1093/gigascience/giac001): Cell Ranger 6,
  STARsolo, alevin, alevin-fry and kallisto run with their default filtering; kallisto returned
  the most cells, many with few genes.
- Montserrat-Ayuso & Esteve-Codina 2024, *BMC Genomics* 25:1124,
  [doi:10.1186/s12864-024-11015-5](https://doi.org/10.1186/s12864-024-11015-5): nuclear fraction
  across atlas datasets.
- Fleming et al. 2023, CellBender, *Nat Methods*,
  [doi:10.1038/s41592-023-01943-7](https://doi.org/10.1038/s41592-023-01943-7): Fig. 1d overlays
  the exonic read fraction on a snRNA-seq rank curve, showing how the composition of droplets
  changes at the knee.
- Caglayan et al. 2022, *Neuron*,
  [doi:10.1016/j.neuron.2022.09.010](https://doi.org/10.1016/j.neuron.2022.09.010): the intronic
  read ratio separates non-nuclear ambient RNA from nuclei.
- Janssen et al. 2023, *Genome Biol* 24:140,
  [doi:10.1186/s13059-023-02978-x](https://doi.org/10.1186/s13059-023-02978-x): benchmarks ambient
  RNA removal, taking cell calls as given.
- Quality control after calling, with mixtures or MADs: miQC (Hippen et al. 2021,
  [doi:10.1371/journal.pcbi.1009290](https://doi.org/10.1371/journal.pcbi.1009290)), SampleQC
  (Macnair & Robinson 2023, [doi:10.1186/s13059-023-02859-3](https://doi.org/10.1186/s13059-023-02859-3)),
  ddqc (Subramanian et al. 2022, [doi:10.1186/s13059-022-02820-w](https://doi.org/10.1186/s13059-022-02820-w)),
  pipeComp (Germain et al. 2020, [doi:10.1186/s13059-020-02136-7](https://doi.org/10.1186/s13059-020-02136-7)).
- Reviews and tutorials with a section on cell calling: Luecken & Theis 2019
  ([doi:10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746)), OSCA (Amezquita et al.
  2020, [doi:10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x)), Heumos et
  al. 2023 ([doi:10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w)), the
  kallisto–bustools protocol (Sullivan et al. 2025,
  [doi:10.1038/s41596-024-01057-0](https://doi.org/10.1038/s41596-024-01057-0)).
- Blog posts: Tom Smith's two posts on estimating the number of true cell barcodes (CGAT blog,
  May 2017), which explain UMI-tools' density method; Lior Pachter, "Rotating the knee (plot) and
  related yoga" (June 2019), on transposed knee plots.

## What is missing

No study has run the methods side by side on the same samples, with the versions stated, and
with an independent ground truth. sc-droplets makes the first half easy: every method on one
matrix, each at a named version. The *Comparison on public data* page (coming with 0.1.0) will
use it on public samples, with unsupervised indicators such as the unspliced fraction of called
and uncalled barcodes.
