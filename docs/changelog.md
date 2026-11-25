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
- The 10x HDF5 reader keeps the file's attributes (chemistry, Cell Ranger version).
- A validation suite on public data (`-m integration`), and an oracle workflow that runs STAR.
