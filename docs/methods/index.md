---
title: Methods
nav_order: 4
has_children: true
---

# Methods

One page per algorithm, grouped by the tool it comes from. Each page has the same sections:
what the method is at a glance, how it works step by step, which releases of its tool behaved
how, its parameters beside the original tool's names for them, how to run it, where our
implementation differs from the original, and how it was checked against the original.

The table below lists the methods the package runs. It is generated from the package
(`docs/_scripts/methods_table.py`); the [Overview](../overview/index.md) lists those still to
come.

<!-- methods-table:start -->
| Method | Tool | Category | Versions | Page |
|---|---|---|---|---|
| `alevin-initial-whitelist` | alevin | knee | v0.11: ≥0.11.0, <0.12.0; v0.12: ≥0.12.0, <0.13.0; v0.13: ≥0.13.0, <0.14.0; v0.14: ≥0.14.0, <0.14.1; v0.14.1: ≥0.14.1, <1.2.0; v1.2: ≥1.2.0 | [alevin initial whitelist](alevin/initial-whitelist.md) |
| `alevinfry-permit-list` | alevin-fry | knee | v0.1: ≥0.1.0, <0.4.2; v0.4.2: ≥0.4.2, <0.12.0; v0.12: ≥0.12.0, <0.18.0; v0.18: ≥0.18.0 | [alevin-fry generate-permit-list](alevinfry/permit-list.md) |
| `bustools-allowlist` | bustools | rank-threshold | default: ≥0.39.3 | [bustools allowlist](bustools/allowlist.md) |
| `cellranger` | Cell Ranger | statistical-test | v2: ≥2.0, <3.0; v3: ≥3.0, <4.0; v7.2: ≥7.2, <8.0; v8.0: ≥8.0, <8.0.1; v8.0.1: ≥8.0.1, <9.0; v9.0: ≥9.0, <10.0; v10.0: ≥10.0, <10.1; v10.1: ≥10.1 | [Cell Ranger cell calling](cellranger/cellranger.md) |
| `cellranger-ordmag` | Cell Ranger | rank-threshold | v2: ≥2.0, <4.0; v7.2: ≥7.2, <8.0; v8.0: ≥8.0, <8.0.1; v8.0.1: ≥8.0.1, <9.0; v9.0: ≥9.0, <10.0; v10.0: ≥10.0, <10.1; v10.1: ≥10.1 | [Cell Ranger OrdMag](cellranger/ordmag.md) |
| `cellbender-prior` | CellBender | density | v0.1: ≥0.1.0, <0.2.0; v0.2: ≥0.2.0, <0.3.0; v0.3: ≥0.3.0 | [CellBender priors](cellbender/prior.md) |
| `cellsweep-empty` | cellsweep | knee | v0.1: ≥0.1.0, <1.0.2; v1.0.2: ≥1.0.2 | [cellsweep empty droplets](cellsweep/empty.md) |
| `dropletutils-barcoderanks` | DropletUtils | knee | v1.0: ≥1.0, <1.10; v1.10: ≥1.10, <1.12; v1.12: ≥1.12, <1.28; v1.28: ≥1.28, <1.30; v1.30: ≥1.30 | [DropletUtils barcodeRanks](dropletutils/barcoderanks.md) |
| `dropletutils-defaultdrops` | DropletUtils | rank-threshold | v1.0: ≥1.0 | [DropletUtils defaultDrops](dropletutils/defaultdrops.md) |
| `dropletutils-emptydrops` | DropletUtils | statistical-test | v1.0: ≥1.0, <1.2; v1.2: ≥1.2, <1.4; v1.4: ≥1.4, <1.6; v1.6: ≥1.6, <1.8; v1.8: ≥1.8, <1.10; v1.10: ≥1.10, <1.12; v1.12: ≥1.12, <1.20; v1.20: ≥1.20, <1.28; v1.28: ≥1.28, <1.30; v1.30: ≥1.30 | [DropletUtils emptyDrops](dropletutils/emptydrops.md) |
| `dropletutils-emptydropscellranger` | DropletUtils | statistical-test | v1.14: ≥1.14, <1.18.1; v1.18.1: ≥1.18.1 | [DropletUtils emptyDropsCellRanger](dropletutils/emptydropscellranger.md) |
| `threshold` | Generic | rank-threshold | — | [Fixed threshold](generic/threshold.md) |
| `top-n` | Generic | rank-threshold | — | [Top N barcodes](generic/top-n.md) |
| `kneedle` | Kneedle (kneed) | knee | v0.5.1: ≥0.5.1, <0.7.0; v0.7: ≥0.7.0, <0.8.6; v0.8.6: ≥0.8.6, <0.9.0; v0.9: ≥0.9.0 | [Kneedle on the barcode-rank curve](kneedle/kneedle.md) |
| `starsolo-basic-awk` | STARsolo | rank-threshold | v1: ≥2.7.2a | [STARsolo soloBasicCellFilter.awk](starsolo/basic-awk.md) |
| `starsolo-cellranger2.2` | STARsolo | rank-threshold | v1: ≥2.7.3a, <2.7.5a; v2: ≥2.7.5a, <2.7.8a; v3: ≥2.7.8a | [STARsolo CellRanger2.2](starsolo/cellranger2-2.md) |
| `starsolo-emptydrops-cr` | STARsolo | statistical-test | v1: ≥2.7.8a | [STARsolo EmptyDrops_CR](starsolo/emptydrops-cr.md) |
| `starsolo-topcells` | STARsolo | rank-threshold | v1: ≥2.7.3a, <2.7.5a; v2: ≥2.7.5a, <2.7.6a; v3: ≥2.7.6a | [STARsolo TopCells](starsolo/topcells.md) |
| `umitools-density` | UMI-tools | density | v0.5: ≥0.5.0 | [UMI-tools whitelist, density valley](umitools/density.md) |
| `umitools-distance` | UMI-tools | knee | v1.0: ≥1.0.0 | [UMI-tools whitelist, distance knee](umitools/distance.md) |
<!-- methods-table:end -->

## Reading a method name

Method ids are `tool-algorithm`. A spec adds the release of the tool to reproduce and any
parameters:

```
cellranger-emptydrops@9.0
dropletutils-emptydrops@1.26:lower=200
top-n:n_cells=3000
```

On the command line, `sc-droplets call <id>` takes the release as `--tool-version` and the
parameters as options; `compare` and `plot` take specs with `-m`.
