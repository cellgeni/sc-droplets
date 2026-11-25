# Status

Working notes on what is in progress and what comes next. Updated as work moves; the plan of
record for the current item is in the section below.

## Released

Nothing yet.

## In progress: 0.1.0

Every knee and threshold heuristic and every statistical test, each checked against its
original tool, with a page per algorithm and the review pages. Decisions (November 2026):

- One name throughout: PyPI `sc-droplets`, import `sc_droplets`, command `sc-droplets`.
- GPL-3.0-or-later. Cell Ranger 2.x (AGPL) and 7.2+ (10x licence) code is never copied: the
  method page is the specification and the implementation is written from it.
- A method spec is `id[@tool_version][:name=value,...]`; each `Variant` covers a half-open range
  of tool releases; `Gap` marks releases whose code is not public (Cell Ranger 4.0–7.1).
- Heavy models are documented, not implemented: CellBender's VAE, EmptyNN, dropkick.
- Validation: exact against the original for deterministic methods; for Monte Carlo methods,
  ≥ 99.5 % agreement of called barcodes with every disagreement inside the noise band.
- Visualisation is part of the package: knee, density, per-method diagnostics, compare, QC,
  splicing (unspliced fraction), and an HTML report.

Milestones, each its own pull request after the scaffold:

- [x] 0. Scaffold (straight to `main`): registry, counting, results, readers, `top-n`,
  `threshold`, `pl.knee` and the terminal knee, CLI (`methods`, `inspect`, `call`, `compare`,
  `plot knee`), simulator, docs skeleton with the Overview and the review of earlier
  comparisons, CI.
- [ ] 1. Shared primitives (Good-Turing variants, multinomial and Dirichlet-multinomial
  likelihoods, Monte Carlo engines per tool, BH variants, R compatibility helpers) + Cell Ranger
  (OrdMag v1/v2/v7, EmptyDrops v3/v7.2/v8/v9/v10.0/v10.1, targeted gradient, antibody and
  high-occupancy filters) + STARsolo (CellRanger2.2 v1/v2, TopCells, EmptyDrops_CR, the awk
  filter). Oracles: 10x public datasets per Cell Ranger version; `STAR --runMode
  soloCellFiltering` 2.7.10a+.
- [ ] 2. DropletUtils: barcodeRanks v1–v4 (port of R `smooth.spline`), emptyDrops at 1.0, 1.2,
  1.8, 1.10, 1.12, 1.20, 1.28, 1.30, emptyDropsCellRanger, defaultDrops. Oracle: R in pinned
  Bioconductor containers.
- [ ] 3. UMI-tools, alevin, alevin-fry, Drop-seq (dropbead), bustools, dropEst bounds, zUMIs,
  Kneedle (kneed 0.8.5 / 0.8.6 / 0.9.0, vendored), cellsweep empty step, CellBender priors.
- [ ] 4. Generic MAD / Otsu / multi-Otsu / KDE valley / GMM (1-D, 2-D) / NB mixture; splicing
  readers, `pl.splicing`, knee tracks, separation scores; `pl.density`, `pl.compare`, `pl.qc`,
  `report`; the comparison-on-public-data page.
- [ ] 5. Overview pages in full, *planned* and *docs only* pages, guides; release 0.1.0.

Needs the maintainer of the GitHub organisation (not doable from here): a PyPI pending trusted
publisher for `sc-droplets` (workflow `publish.yml`, environment `pypi`), the
`quay.io/cellgeni/sc-droplets` repository with `QUAY_USERNAME` / `QUAY_TOKEN` secrets, and
GitHub Pages from `main`, `/docs`.

## Follow-ups noticed

- anndata 0.13 lists `X` as the layer `None` in `adata.layers.keys()`; filter it wherever layer
  names are listed.

## Next candidates (after 0.1.0)

- 0.2: CB2, DIEM, DropletQC (with a BAM reader for the intron tag), DropSift, dropEst classifier,
  alevin naive-Bayes whitelist, valiDrops, SiftCell, miQC-style mixture regression.
- 0.3: ATAC and multiome: Cell Ranger ATAC (NB mixture, source 1.2), Cell Ranger ARC,
  EmptyDropsMultiome, BAP, SADE.
