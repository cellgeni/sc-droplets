# Other methods, reviews and mixture models

Everything outside Cell Ranger, STARsolo, CellBender, cellsweep and DropletUtils. Heuristics were
read from source; DOIs checked against PubMed where given.

Two version changes affect the package directly:

- salmon 2.0 (Rust rewrite, June 2026) removed `salmon alevin`; single-cell quantification moved
  to alevin-fry. The last C++ release is 1.12.0 (`cpp` branch). Alevin's knee and naive-Bayes
  whitelisting are historical (0.14–1.x).
- kneed changed twice in 2026: 0.8.6 (March) pauses detection after a local minimum until the
  next local maximum; 0.9.0 (October) mirrors the curve by x value instead of by index.

## A. Reviews, benchmarks and overviews (2017–2026)

### Benchmark sections inside method papers (run by the method's authors)

| Citation | What it compares |
|---|---|
| Ni Z, Chen S, Brown J, Kendziorski C. CB2. *Genome Biol* 2020;21:137. doi:10.1186/s13059-020-02054-8 | Only EmptyDrops; simulations from ten 10x datasets; CB2 found 24 % more cells on average (4–81 %) at controlled FDR. |
| Alvarez M et al. DIEM. *Sci Rep* 2020;10:11019. doi:10.1038/s41598-020-67513-5 | snRNA-seq (adipose ×6, preadipocytes, mouse brain) + PBMC 68k vs a Cell-Ranger-style cut-off and EmptyDrops (FDR 0.05); spliced-read fraction as a contamination proxy. |
| Yan F, Zhao Z, Simon LM. EmptyNN. *Patterns* 2021;2:100311. doi:10.1016/j.patter.2021.100311 | Ground truth from hashing and genetic demultiplexing; AUROC EmptyNN 94.73 %, CR v2 86.88 %, CR v3 90.16 %, EmptyDrops 79.87 %, CellBender 88.80 %. |
| Heiser CN et al. dropkick. *Genome Res* 2021;31:1742–1752. doi:10.1101/gr.271908.120 | Simulated and real data vs thresholds and EmptyDrops (plus CellBender runtime). |
| Muskovic W, Powell JE. DropletQC. *Genome Biol* 2021;22:329. doi:10.1186/s13059-021-02547-0 | Cell Ranger and EmptyDrops calls still contain high-ambient empty droplets and damaged cells; nuclear fraction as an add-on. |
| Xi J, Park SR, Lee JH, Kang HM. SiftCell. *Cell Syst* 2023;14:620–628.e3. doi:10.1016/j.cels.2023.06.002 | vs Cell Ranger cut-off, EmptyDrops, DIEM and CellBender on PBMC v3, E18 brain nuclei and a Drop-seq mix. |
| Kavaliauskaite G, Madsen JGS. valiDrops. *NAR Genom Bioinform* 2023;5:lqad101. doi:10.1093/nargab/lqad101 | Broadest in-paper benchmark: EmptyDrops, scCB2, DIEM, EmptyNN, CellBender, dropkick on 47 samples from 5 datasets; 9 filtering + 12 integration metrics. |
| Megas S, Lorenzi V, Marioni JC. EmptyDropsMultiome. *Genome Biol* 2024;25:121. doi:10.1186/s13059-024-03259-x | Multiome vs CellRanger-arc and EmptyDrops. |
| Petukhov V et al. dropEst. *Genome Biol* 2018;19:78. doi:10.1186/s13059-018-1449-6 | KDE quality classifier vs UMI thresholds. |
| Srivastava A et al. alevin. *Genome Biol* 2019;20:65. doi:10.1186/s13059-019-1670-y | Knee + naive-Bayes whitelisting vs Cell Ranger calls. |
| He D et al. alevin-fry. *Nat Methods* 2022;19:316–322. doi:10.1038/s41592-022-01408-3 | Permit-list strategies. |

### Independent studies, QC papers and reviews

- Brüning RS et al. 2022, *GigaScience* 11:giac001, doi:10.1093/gigascience/giac001 — Cell Ranger
  6, STARsolo, Alevin, Alevin-fry, Kallisto with default filtering; Kallisto returned the most
  cells.
- Montserrat-Ayuso T, Esteve-Codina A. 2024, *BMC Genomics* 25:1124, doi:10.1186/s12864-024-11015-5
  — DropletQC nuclear fraction across HCA, Tabula Muris, Tabula Sapiens.
- Post-calling QC with mixtures or MADs: miQC (Hippen 2021, doi:10.1371/journal.pcbi.1009290),
  SampleQC (Macnair & Robinson 2023, doi:10.1186/s13059-023-02859-3), ddqc (Subramanian 2022,
  doi:10.1186/s13059-022-02820-w), pipeComp (Germain 2020, doi:10.1186/s13059-020-02136-7).
- Ambient-removal benchmark: Janssen P et al. 2023, *Genome Biol* 24:140,
  doi:10.1186/s13059-023-02978-x.
- Reviews with cell-calling sections: Luecken & Theis 2019 (doi:10.15252/msb.20188746), OSCA
  (Amezquita 2020, doi:10.1038/s41592-019-0654-x), Heumos 2023 (doi:10.1038/s41576-023-00586-w),
  kb protocol (Sullivan 2025, doi:10.1038/s41596-024-01057-0).
- Blogs: Tom Smith, CGAT blog "Estimating the number of true cell barcodes" parts 1–2
  (May 2017); Lior Pachter, "Rotating the knee (plot) and related yoga" (2019-06-24).
- 2025–2026 preprints: Ambimux, AmbientMapper, SADE (section C12).
- CellBender (Fleming 2023, doi:10.1038/s41592-023-01943-7) Fig. 1d overlays the exonic read
  fraction on the rat6k snRNA-seq rank curve; Caglayan et al. 2022 (*Neuron*,
  doi:10.1016/j.neuron.2022.09.010) uses the intronic read ratio to separate non-nuclear ambient
  RNA from nuclei.

**No dedicated independent review or benchmark of empty-droplet methods was found.** valiDrops is
the closest multi-method benchmark.

## B. Knee and barcode-rank heuristics

Several "different" methods are one computation:

- maximum distance to a chord on the cumulative curve: UMI-tools `distance`, alevin-fry `-k`,
  dropbead, zUMIs `uik` (vertical distance, same argmax); kneed refines the idea;
- KDE valley on log counts: UMI-tools `density`, alevin's Gauss correction, BAP, DropSift, DropletQC;
- "robust max / 10": alevin `--expectCells`, alevin-fry `-e`, ScaleBio TopCell, CeleScope `auto`,
  bustools in spirit.

### B1. `umitools_distance` — UMI-tools `whitelist --knee-method=distance` (default from 1.0.0)

Smith, Heger, Sudbery, *Genome Res* 2017, doi:10.1101/gr.209601.116; CGATOxford/UMI-tools (MIT),
`umi_tools/whitelist_methods.py`. `whitelist` added in 0.5.0 (Aug 2017, density only); 1.0.0 (Feb
2019) made `distance` default, added `--allow-threshold-error`, `--ed-above-threshold`. Current
1.1.6. Defaults: `--method=reads`, `--subset-reads=100,000,000`, `--error-correct-threshold=1`;
`--expect-cells` rejected with `distance`; `--set-cell-number` overrides.

```
c = counts sorted desc ; v = cumsum(c)
knee(v): line from (0,v0) to (n-1,v[n-1]); d_i = perpendicular distance of (i,v_i); return argmax d (first max)
idx = knee(v); if idx==0: error
repeat up to 100 times: prev=idx; idx = knee(v[:idx*3]); stop when idx==prev
cells = top idx+1 barcodes in Counter.most_common order (ties arbitrary)
if set_cell_number N: thr = c[N]; cells = {count > thr}   # strict, may give < N
```

### B2. `umitools_density` — `--knee-method=density` (pre-1.0 default)

Prefilter counts > 0.001 × max; `gaussian_kde(log10 counts, bw_method=0.1)` on a 10000-point grid;
local minima by `argrelextrema(np.less)`.

```
for each local min m, from highest count to lowest:
    n_pass = #barcodes with count > 10^x[m]
    if expect_cells E: accept the first m with 0.1*E < n_pass <= E
    else: accept the first m with m >= 0.2*10000 and (max(log) - x[m] > 0.5 or x[m] < max(log)/2)
cells = {count > 10^x[m]}; none accepted -> error unless --allow-threshold-error
```

### B3. salmon alevin (0.14–1.x) — `alevin_knee`, `alevin_expect_cells`, `alevin_force_cells`

`src/Alevin.cpp` `sampleTrueBarcodes` at v1.10.0; GPL-3 (salmon 1.x). Defaults:
`maxNumBarcodes=100000`, `lowRegionMinNumBarcodes=200`, low-region max 1000, `freqThreshold=10`;
counts are corrected-read frequencies.

Knee (`getLeftBoundary`), a tangent from the origin to the log cumulative curve:

```
M = min(1e5, n); f = top-M counts desc
y[i] = ln(sum of the i+1 smallest of the top-M)
for j: skip while y[j]==y[0]; slope = y[j]/j; accept the first j where slope*i <= y[i] for all i > j
knee = M-j, moved left past ties
```

Gauss correction: KDE on log10 counts (kernel variance 0.01 × sample variance, 10000-point grid
from the top); accept the first local minimum whose barcode count is in [0.1·knee, knee].
`--expectCells E`: `b = max(1, int(0.01*E)); cut = max(1, 0.1*f[b])`; first i in
[b, min(n, 2E)) with f[i] < cut → i+1 barcodes. `--forceCells F`: top F, minus barcodes with < 10
reads. Low-confidence region `clamp(0.5·knee, 200, 1000)` below the knee trains the classifier
(C11).

### B4. alevin-fry — `alevinfry_knee`, `_expect_cells`, `_force_cells`, `_unfiltered`

COMBINE-lab/alevin-fry (BSD-3) `src/knee_finding.rs`, `src/cellfilter.rs`; latest 0.19.0; counts
are reads per corrected barcode.

- `-k`: UMI-tools distance with slack multiplier **5** (not 3), slice `cfreq[0..min(len-1, idx*5))`,
  `>=` ties (last maximum wins), max 100 iterations; threshold `freq[knee]`, keep count ≥ threshold.
- `-e E`: `idx = min(n-1, round(0.99·E)); thr = max(1, round(freq[idx]/10))` — **0.99·E, not the
  0.01·E** of Cell Ranger v2 and alevin; much more permissive. Unchanged since v0.8.
- `-f F`: `thr = freq[F-1]`; keep ≥ thr (ties may give > F).
- `-u`: barcodes on an external list with `--min-reads` ≥ 10. `-b`: explicit list.
- simpleaf: exactly one of `--knee`, `--unfiltered-pl`, `--forced-cells`, `--expect-cells`,
  `--explicit-pl` is required; `--min-reads 10`.

### B5. Drop-seq

- `dropseq_cumfrac_knee`: Macosko 2015 (doi:10.1016/j.cell.2015.05.002) cumulative fraction of reads
  vs ranked barcodes, inflection chosen by eye, then `NUM_CORE_BARCODES`.
- `dropbead_knee`: Alles 2017 (doi:10.1186/s12915-017-0383-5), `estimateCellNumber(max.cells=15000)`:
  top 15000, `cs` = cumulative fraction (within the top 15000), `x = (1..N)/N`, return
  `argmax |cs - x|`.
- `dropseq_min_transcripts`: `SelectCellsByNumTranscripts MIN_TRANSCRIPTS_PER_CELL` (fixed UMI
  threshold). broadinstitute/Drop-seq (MIT).

### B6. `bustools_allowlist` and kb `--filter bustools`

Melsted 2019 (doi:10.1093/bioinformatics/btz279); BUStools/bustools (BSD-2)
`src/bustools_whitelist.cpp` (unchanged since June 2019; renamed allowlist Oct 2023):

```
read BUS records (sorted by barcode) until >= 200 distinct barcodes
avg = mean read count of the top 10 of those (integer division)
thr = int(avg * (1 - 0.99^bclen))          # 16 bp -> ~0.149*avg
output every barcode with reads >= thr     # or --threshold
```

The ~200 barcodes are the lexicographically first, so effectively a random sample. kb count
`--filter bustools` runs it then correct/sort/count into a filtered matrix.

### B7. `dropest_cell_number` — dropestr `EstimateCellsNumber` (GPL-3)

```
u = UMI counts desc; on log(rank) vs log(u)
lag = round(0.05*n); y = lagged first derivative; y2 = lagged derivative of y
max = exp(x2 at the start of the longest run where y2 > 0)
expected = round(exp(x[argmin y[1:max-lag]]) - lag/2)
min = round(0.75*expected)
```

Labels High (rank ≤ min), Low (rank ≥ max), Unknown — the classifier's training sets (C9).

### B8. `zumis_uik` — zUMIs (GPL-3), `barcodeIDFUN.R`

Parekh 2018 (doi:10.1093/gigascience/giy059).

```
keep barcodes with reads >= nReadsperCell (100); sort desc; cs = cumsum(reads)
ntop = inflection::uik(index, cs/1000); cut = reads[ntop]; cells = reads >= cut
fewer than 10 cells -> fall back to the top 100
```

`uik` (R `inflection` 1.3.7) = argmax vertical distance from the chord (one pass of UMI-tools
distance). Before commit 6d528b2 (Nov 2019, ~v2.5.x) zUMIs used `mclust::Mclust(log10(reads),
modelNames=c("E","V"))`, cut = 10^(qnorm(0.01, μ_top, σ)).

### B9. Kneedle

Satopää 2011 (doi:10.1109/ICDCSW.2011.20); `kneed` (BSD-3, 0.9.0). Defaults `S=1.0`,
`curve="concave"`, `direction="increasing"`, `interp_method="interp1d"`, `online=False`. Normalise
to [0,1]; flip to concave increasing; difference curve `y_d = y_n − x_n`; local maxima; threshold
`T = y_lmx − S·mean(Δx_n)`; knee = local max after which `y_d` drops below `T` before the next
maximum. Used by Sanger HGI `yascp` and `nf-cellbender` (`get_estimates_from_umi_counts.py`: linear
rank vs UMI, `convex`, `decreasing`, S from 5000 down by 100 until a knee), and the Crick `carmack`
scATAC caller (`concave`, `decreasing`, `online=True`).

### B10. `bap_density_knee` (scATAC beads)

Lareau 2019 (doi:10.1038/s41587-019-0147-6); caleblareau/bap (MIT) `bin/R/00_knee_CL.R`. R port of
UMI-tools density with prefilter 0.001 × mode and absolute bandwidth 0.1; fallbacks 500 / 0.

### B11. Vendor and pipeline heuristics

- BD Rhapsody (documentation only, closed): second derivative of the log10 cumulative-reads curve;
  minima 2.0×, 1.75×, 1.5× deeper than the mean; refined ATAC mode uses Gaussian mixtures.
- ScaleBio ScaleRna: default EmptyDrops-like (FDR 0.001); TopCell `minUTC=100`, `thr =
  percentile_99 / minCellRatio(10)`; fixed `UTC`.
- Singleron CeleScope: default `EmptyDrops_CR`; `auto` = Cell Ranger v2 ordmag (expected 8000).
- inDrops: fixed `absolute_min_reads=250`, `min_reads=750`.

### B12. Generic heuristics in practice

Top-N; fixed thresholds (UMIs, reads, genes); MAD outliers (`scuttle::isOutlier` nmads=3 on log;
zUMIs `calcMADS`; ddqc per-cluster MAD 2.0); Otsu / multi-Otsu (Otsu 1979,
doi:10.1109/TSMC.1979.4310076; dropkick uses multi-Otsu for training labels).

## C. Statistical and ML methods

- **C1 CB2 / scCB2** (Ni 2020; GPL-3; R). `FDR_threshold=0.01`, `lower=100`, `upper` automatic,
  `TopNGene=30000`. `upper` = barcodeRanks knee iterated with `lower <- inflection +
  (knee-inflection)/10` until stable; ≥ upper are cells; Good-Turing ambient from ≤ lower;
  Pearson correlation with ambient; hierarchical clusters (chunks of 1000) in (lower, upper);
  cluster test with 1000 multinomial sims; BH ≤ 0.01; the rest through emptyDrops.
- **C2 EmptyNN** (Yan 2021; no licence; R/keras) — docs only. PU learning: N = ≤ 100 counts
  (subsample 10k), U = > 100; top 2000 genes; 5 iterations × 5 folds of an MLP; mean > 0.5 → cell.
- **C3 DIEM** (Alvarez 2020; GPL-3; R/C++). Debris set < 100 counts; k-means on PCs of the test set
  (k_init 20); semi-supervised multinomial-mixture EM; clusters with mean genes < 200 → debris;
  debris score from debris genes; clean if score < 0.5 and genes ≥ 200.
- **C4 valiDrops** (Kavaliauskaite 2023; GPL-3; R). Segmented regression on the smoothed log-log
  rank curve; adaptive QC (mito via 2-component normalmixEM + UIK, segmented logFeatures ~ logUMI
  residuals, coding fraction); clustering-based expression filter.
- **C5 dropkick** (Heiser 2021; MIT; Python) — docs only. Multi-Otsu on arcsinh(n_genes) for
  labels; HVG features; glmnet elastic-net logistic regression (alpha 0.1).
- **C6 DropletQC** (Muskovic 2021; MIT; R). Nuclear fraction NF = intronic / (intronic + exonic)
  reads (BAM `RE` tag or velocyto). Empty if NF < first KDE density minimum, rescued if NF > 0.05
  and UMI > 1000; damaged cells by per-cell-type `mclust(G=1:2, "EEI")` on (NF, log10 UMI).
- **C7 SiftCell** (Xi 2023; GPL ≥ 2; R + xgboost). Shuffle UMI-to-barcode assignments as a null;
  XGBoost on PCs; multinomial mixture for ambient fraction.
- **C8 DropSift** (Broad/McCarroll; MIT; R; no paper). KDE on log10 UMI for exemplar bounds
  (`PitAfterHighestPeak` / `PitBetweenHighestPeaks`, bandwidth grid 0.25–3); HDR regions in
  pct_intronic; RBF SVM on log UMI, pct_intronic, pct_mt, empty gene-module score; uses CellBender
  `frac_contamination` by default.
- **C9 dropEst classifier** (GPL-3). KDE Bayes classifier on reads/UMI, UMIs/gene, one-UMI-gene
  fraction (+ optional intergenic, not-aligned, mito) after robust sparse PCA; needs BAM features.
- **C10 EmptyDropsMultiome** (Megas 2024; GPL ≥ 3; R). Per-modality DM EmptyDrops, combined
  p-values, k-means line zones; `FDR_multi < 0.001`.
- **C11 alevin naive-Bayes whitelist** (`src/WhiteList.cpp`). Features: mapping rate, dedup rate,
  gene UMI stats, genes expressed (+ mito, rRNA); true = top half of the knee set, ambiguous =
  bottom half, false = low-confidence region; Gaussian NB with effectively equal priors.
- **C12 2025–2026 preprints**: SADE (scATAC entropy + spline knee; doi:10.64898/2026.10.20.760589),
  AmbientMapper / scifi-demux (doi:10.64898/2026.09.11.750809), Ambimux (genotype-based, multiome;
  doi:10.1101/2025.08.21.671671).
- **C13 out of scope as cell callers**: SoupX, CellBouncer (doi:10.1101/2025.03.23.644821),
  post-calling QC (scAutoQC, miQC, SampleQC, ddqc, scPipe `detect_outlier`).

## D. Mixture-model approaches

In tools or the literature:

1. zUMIs legacy (≤ Nov 2019): `Mclust(log10(reads), modelNames=c("E","V"))`, cut at the 1st
   percentile of the top component.
2. cellsweep `mx_filter`: 2-component GMM on log1p totals, cut at the max-entropy point.
3. scPipe `detect_outlier` (Tian 2018, doi:10.1371/journal.pcbi.1006361): robust MCD or mclust on
   QC features.
4. DropletQC damaged cells: bivariate EEI GMM on (NF, log10 UMI).
5. valiDrops: 2-component normalmixEM on log features.
6. miQC: flexmix mixture of regressions, mito% ~ genes.
7. SampleQC: robust multivariate GMM per sample group.
8. sctk scAutoQC (Teichlab, GPL-3): per-metric GMM (1–10 components by BIC) on log1p metrics.
9. BD Rhapsody refined ATAC: GMMs on log transposase sites and FRiP.
10. Count-space mixtures: DIEM, SiftCell-Mix, dropEst KDE class densities.
11. K-means as a hard GMM: cellranger-arc (K=2 on log ATAC × log GEX), EmptyDropsMultiome.
12. Cell Ranger ATAC: 2-component NB mixture on fragments in peaks (see the Cell Ranger notes).

Reasonable additions of our own (generic, not from one tool):

| id | Feature space | Parameters |
|---|---|---|
| `gmm` | log10 total UMIs (or reads), barcodes ≥ `min_counts` | k ∈ {2,3} by BIC; cells = top-mean component; rule `posterior > 0.5` or zUMIs-style quantile 0.01 of the top component |
| `gmm-2d` | (log10 UMI, log10 genes), (log10 UMI, NF), (log10 UMI, logit mito) | full covariance; k by BIC ≤ 4 |
| `nb-mixture` | raw totals | 2-component NB EM (Cell Ranger ATAC style); threshold where posteriors cross or odds rule |
| `otsu`, `multiotsu` | log10 UMIs or arcsinh n_genes | 256 bins; 2 or 3 classes |
| `kde-valley` | log10 UMI (or NF) | bandwidth UMI-tools 0.1 / Scott / Silverman / DropSift grid; rule umitools, pit-after-highest-peak, pit-between-two-highest, DropletQC first derivative sign |

## Units and ties

UMI-tools, alevin, alevin-fry and bustools default to **reads**; most R tools use UMIs. Ties:
UMI-tools `>` with arbitrary order; alevin-fry `>=`; alevin skips ties.

## Summary table

| Method | Category | Language | Licence | Python port | Priority |
|---|---|---|---|---|---|
| UMI-tools distance / density | knee / KDE valley | Python | MIT | trivial | high |
| alevin knee, expect, force | tangent knee + KDE / ordmag / top-N | C++ (removed in salmon 2.0) | GPL-3 (1.x) | yes | medium |
| alevin NB whitelist | Gaussian NB | C++ | GPL-3 | needs read features | low-medium |
| alevin-fry knee, expect, force, unfiltered | distance knee (slack 5) / ordmag 0.99·E / top-N / min reads | Rust | BSD-3 | trivial | high |
| Drop-seq cumfrac, dropbead | distance to diagonal | Java/R | MIT / GPL | trivial | medium |
| bustools / kb `--filter` | sample-mean heuristic | C++/Python | BSD-2 | trivial | medium |
| dropEst cell number / classifier | log-log derivative / KDE Bayes | R | GPL-3 | yes / needs BAM features | medium / low-medium |
| zUMIs uik, legacy mclust | UIK knee / GMM | R | GPL-3 | yes | medium |
| Kneedle | knee | Python | BSD-3 | yes | high |
| BAP | KDE valley (ATAC) | R | MIT | yes | low |
| BD Rhapsody | second derivative | closed | proprietary | approximate only | low |
| CB2 | clustered multinomial test + ED | R | GPL-3 | yes | high |
| EmptyNN | PU-learning NN | R/keras | none | docs only | — |
| DIEM | multinomial mixture EM | R/C++ | GPL-3 | yes | high |
| valiDrops | segmented knee + QC + clustering | R | GPL-3 | moderate | medium |
| dropkick | multi-Otsu + elastic net | Python | MIT | docs only | — |
| DropletQC | NF KDE valley + bivariate GMM | R | MIT | needs spliced/unspliced | high |
| SiftCell | shuffle null + XGBoost | R + Python | GPL ≥ 2 | yes | medium |
| DropSift | KDE exemplars + RBF SVM | R | MIT | yes | medium-high (nuclei) |
| EmptyDropsMultiome | DM tests, 2 modalities | R | GPL ≥ 3 | needs ED | low (multiome) |
| SADE | entropy + spline knee (ATAC) | Python notebooks | none | yes | low |
