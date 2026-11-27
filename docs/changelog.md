---
title: Changelog
nav_order: 8
---

# Changelog

## Unreleased

The scaffold: the method registry with tool versions and gaps, readers for 10x HDF5 (v2 and v3),
Matrix Market folders (Cell Ranger, STARsolo, kb, alevin-fry) and pipeline output folders, the
generic `top-n` and `threshold` methods, `scd.call`, `scd.compare`, `scd.agreement`, the knee
plot (`scd.pl.knee`, and in the terminal), the `sc-droplets` command (`methods`, `inspect`,
`call`, `compare`, `plot knee`), and this site with the Overview and the review of earlier
comparisons.

Then Cell Ranger and STARsolo, each checked against the original on public data:

- `cellranger` (2.x, 3.x, 7.2, 8.0, 8.0.1, 9.0, 10.0, 10.1) and `cellranger-ordmag`, with the
  antibody-aggregate filter, the 8.0+ global-minimum and mitochondrial filters, and the
  high-occupancy GEM filter for multiplexed Flex (`high_occupancy_gems`). 10x's own filtered
  matrices are reproduced exactly for 2.1, 3.0, 7.2 (4 Flex samples), 8.0 and 9.0, and up to 15
  borderline barcodes of 15,257 for 10.0.
- `starsolo-cellranger2.2`, `starsolo-topcells`, `starsolo-emptydrops-cr` and
  `starsolo-basic-awk`, every release where they changed; identical to `STAR --runMode
  soloCellFiltering` 2.7.10a and 2.7.11b.
- Shared building blocks reproducing each tool's arithmetic: both Simple Good-Turing variants,
  numpy's classic quicksort (for tie order), C++ `mt19937` and libstdc++'s
  `discrete_distribution`, and the Cell Ranger and STAR Monte Carlo engines.
- `scd.pl.diagnostics` and `sc-droplets plot diagnostics`: each call's own view.

Then DropletUtils, checked against every release from 1.0 to 1.30 run in R:

- `dropletutils-barcoderanks` (the knee by spline curvature, 1.0 with its indexing bug, 1.10 and
  1.12; by chord distance, 1.28; by curve tracing, 1.30), `dropletutils-emptydrops` (ten
  variants: multinomial or Dirichlet-multinomial, R's or PCG32 random numbers, the rounding,
  `by.rank` and `test.ambient` rules of each release), `dropletutils-emptydropscellranger` (with
  1.14's median bug) and `dropletutils-defaultdrops`. On 10x's pbmc_1k_v3 raw matrix every knee,
  p-value, estimated alpha, FDR and call is R's own, with the same `set.seed`.
- New building blocks: R's random numbers and dqrng's seeds; PCG32 and Boost.Random's
  alias-table, exponential and gamma sampling; R's `smooth.spline` and `predict`; R's `lgamma`,
  long-double `sum` and `cumsum`, `optimize`, `quantile` and `p.adjust`; edgeR's Good-Turing.
- The 10x HDF5 reader keeps the file's attributes (chemistry, Cell Ranger version).

Then UMI-tools:

- `umitools-distance` (`whitelist --knee-method=distance`, 1.0.0+) and `umitools-density`
  (`--knee-method=density`, 0.5.0+), with `--set-cell-number`, `--expect-cells` and
  `--allow-threshold-error`. Given the same per-barcode counts, the barcodes of UMI-tools 0.5.5,
  1.0.0 and 1.1.6's own functions on pbmc_1k_v3 and pbmc4k.
- `sc-droplets plot diagnostics` draws their cumulative curve and chord, or their density and
  valleys.

Then Kneedle:

- `kneedle`: kneed's `KneeLocator` on the barcode-rank curve, in each behaviour its releases
  have had (0.5.1, 0.7.0, 0.8.6, 0.9.0), by default as Sanger's yascp runs it (convex, decreasing,
  at least 1,000 UMIs, sensitivity lowered from 5,000 until a knee appears). Identical knees to
  kneed 0.6.0, 0.8.5, 0.8.6 and 0.9.0 themselves on pbmc_1k_v3 and pbmc4k.

Then cellsweep:

- `cellsweep-empty`: cellsweep's empty-droplet step, with the rule of each release (the second
  derivative of the smoothed log curve to 1.0.1, kneed's online knee from 1.0.2, with the
  behaviour of whichever kneed is installed), `mx_filter` and the manual cutoffs, in the dtype of
  cellsweep's matrix. Identical cutoffs to cellsweep 1.0.1 and 1.0.2 (with kneed 0.8.5, 0.8.6,
  0.9.0) themselves on pbmc_1k_v3 and pbmc4k.

Then CellBender:

- `cellbender-prior`: what `cellbender remove-background` decides before training, in each
  release's form (0.1.0; 0.2.x; 0.3.0 to 0.4.0): the priors on counts in cells and in empty
  droplets, and the probable cells, the other analysed droplets and the surely-empty droplets.
  Identical priors and classes to CellBender 0.1.0, 0.2.2, 0.3.0, 0.3.2 and 0.4.0 themselves on
  pbmc_1k_v3 and pbmc4k.

Then alevin-fry and bustools:

- `alevinfry-permit-list`: alevin-fry's `generate-permit-list` rules (the knee, expected cells,
  forced cells, the unfiltered floor), with the rules 0.12–0.17 apply to multi-barcode (Flex)
  input. Identical cells to programs built from alevin-fry 0.1.0, 0.11.2, 0.17.1 and 0.19.0's own
  source on pbmc_1k_v3 and pbmc4k.
- `bustools-allowlist`: bustools' allowlist (`kb --filter bustools`), with its sample of the first
  blocks of the BUS file. Identical lists to bustools 0.39.3, 0.41.0, 0.43.2 and 0.45.1 on BUS
  files written from pbmc_1k_v3 and pbmc4k.

Then alevin:

- `alevin-initial-whitelist`: salmon alevin's initial whitelist (0.11.0 – 1.10.3): the tangent
  knee with its kernel-density correction, Cell Ranger 2's expected cells, forced cells and
  `--keepCBFraction`, then the low-confidence band, in each of the six forms these took.
  Identical counts and whitelists to programs built from salmon 0.11.3, 0.13.1, 0.14.0, 1.1.0 and
  1.10.3's own `Alevin.cpp` on pbmc_1k_v3 and pbmc4k; the cells match up to barcodes tied at the
  knee.

Then Drop-seq, dropEst and zUMIs:

- `dropseq-dropbead`: dropbead's `estimateCellNumber` (0.3.1).
- `dropest-bounds`: dropestr's `EstimateCellsNumber` (0.6.0 – 0.7.9), with the High / Unknown / Low
  bands of its classifier.
- `zumis-cell-barcodes`: zUMIs' automatic cell detection with inflection's `uik` (2.5.6b on); the
  earlier Gaussian mixture is a documented gap (mclust initialises on an unseeded random subset).
- Identical results to the tools' own R functions, sourced at their commits and tags, on
  pbmc_1k_v3 and pbmc4k.
- A validation suite on public data (`-m integration`), and an oracle workflow that runs STAR.
