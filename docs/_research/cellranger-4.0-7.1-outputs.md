# Cell Ranger 4.0–7.1 from 10x's own outputs

November 2026. Cell Ranger 4.0 to 7.1 published no source (only binaries under 10x's EULA, which
this project does not use). 10x does publish, for its demonstration datasets, the raw and filtered
matrices those releases made, and most dataset pages say which `--expect-cells` the run used. This
note records what those outputs show about each release's cell calling, without reading any of
10x's code. Tag: **[OUT]**, inferred from the tool's published outputs.

## Where the datasets come from

10x's dataset catalogue answers `https://www.10xgenomics.com/api/search?document=dataset` (791
entries, `pipeline` gives the release); each dataset page's data
(`/_next/data/<build>/en/datasets/<slug>.json`) lists its `cf.10xgenomics.com` files and, often,
"Run with --expect-cells=N". The datasets below are in `tests/datasets.py`, SHA-256-pinned.

| Name | Release | Chemistry (HDF5) | `--expect-cells` | 10x cells |
|---|---|---|---|---|
| `pbmc10k_nextgem_4_0` | 4.0.0 | 3′ v3 | 10,000 | 10,985 |
| `neuron10k_nextgem_4_0` | 4.0.0 | 3′ v3 (mouse) | 10,000 | 11,316 |
| `pbmc_parent_4_0` | 4.0.0 | 3′ v3 | not given | 10,194 |
| `pbmc10k_5p_v2_5_0` | 5.0.0 | 5′ R2-only, + 17 antibodies | 10,000 | 10,548 |
| `brain_tumor_6_0` | 6.0.0 | 3′ v3 | 2,000 | 1,615 |
| `breast_cancer_6_0` | 6.0.0 | 3′ v3 | 7,500 | 5,680 |
| `hgmm1k_6_0` | 6.0.0 | 3′ v3, human + mouse | 1,000 | 1,112 |
| `hgmm1k_lt_6_0` | 6.0.0 | 3′ v3 LT, human + mouse | 1,000 | 808 |
| `pbmc10k_controller_6_1` | 6.1.0 | 3′ v3 | 10,000 | 11,485 |
| `pbmc10k_chromium_x_6_1` | 6.1.0 | 3′ v3 | 10,000 | 11,996 |
| `pbmc10k_chromium_x_6_1_2` | 6.1.2 | the same raw file | 10,000 | 11,996 (the same cells) |
| `pbmc500_lt_6_1` | 6.1.0 | 3′ v3 LT | 500 | 705 |
| `hgmm10k_controller_6_1` | 6.1.0 | 3′ v3, human + mouse | 10,000 | 10,578 |
| `mouse_brain_nuclei_7_0` | 7.0.0 | 3′ v3 nuclei | not given | 7,377 |
| `jejunum_nuclei_7_0` | 7.0.0 | 3′ v3 nuclei | not given | 4,392 |
| `pbmc5k_7_0_1` | 7.0.1 | 3′ v3 | not given | 5,140 |
| `flex_spleen_ffpe_7_1` | 7.1.0 | Flex, singleplex (`multi`) | not given | 8,268 |

The 4.0 raw matrices hold the whole 3′ v3 whitelist (6,794,880 barcodes); from 6.0 they hold only
barcodes with reads, as 6.0's release notes say.

## Release notes (archived, 1.0–7.1), on calling

- 5.0: "upgraded protein aggregation detection and filtering … directly using the protein counts".
- 6.0: raw matrices only hold barcodes with at least one read, "which results in minor changes to
  the results of the cell calling algorithm"; species annotation of EmptyDrops-rescued cells fixed.
- 6.1: "changed certain parameters in the cell calling algorithm"; EmptyDrops for multi-species
  experiments uses each species' UMI counts separately.
- 7.0: `--expect-cells` no longer needed: estimated.
- 7.1: the estimate's upper range restricted to the low end of the EmptyDrops range (45,000;
  262,000 before); OrdMag no longer calls every barcode when there are very few cells.

## What each part of the rule does to the result

`cellranger-4.0-7.1-survey.py` (beside this note) recombines the parts of the 3.x and 7.2 rules:
OrdMag (3.x or 7.2), the simulation's random stream (continuing the stage's stream after OrdMag's
bootstraps, as 3.x does, or a fresh `RandomState(0)` per test, as 7.2 does), the block refill (3.x,
or 7.2's early refill for gaps of 20–999), the candidate limit (20,000 or none), the ambient ranks
(45,000–90,000 or LT's 4,500–9,000), EmptyDrops across genomes or per genome, and the antibody
aggregate count rule.

"own stream" below is 3.x with the simulation on a fresh `RandomState(0)` (3.x OrdMag, 3.x refill);
it is the `v4` rule, `v5` with the aggregate filter.

| Dataset | 3.x | own stream | 7.2 (with the page's expect) | 7.2 (estimated) |
|---|---|---|---|---|
| pbmc10k 4.0 | identical | identical | identical | 2 extra (OrdMag) |
| neuron10k 4.0 | 8 missing, 23 extra (p 0.007–0.012) | **identical** | 1 missing, 35 extra | 5 missing, 55 extra |
| pbmc parent 4.0 (default 3,000) | 1 missing, 6 extra (p 0.008–0.011) | **identical** | — | 2 missing, 124 extra |
| 5′ PBMC 5.0 | 5 extra (ranks 1–13: aggregates) | identical with aggregates removed | identical | identical |
| brain tumour 6.0 | 1 missing, 1 extra | 1 missing, 2 extra (p 0.008–0.013) | 0 missing, 1 extra | 1 missing, 2 extra |
| breast cancer 6.0 | 1 missing, 1 extra | **identical** | 1 missing | 1 extra |
| hgmm 1k 6.0 | 1 extra | **identical** | — | 125 missing (per genome) |
| hgmm 1k LT 6.0 | 147 missing (no ambient barcodes at 45,000–90,000: OrdMag only) | **identical** with LT ranks (7.2 refill: 2 missing) | — | 97 missing, 2 extra |
| PBMC Controller 6.1 | 1 extra (p 0.0077) | 1 extra (p 0.0093) | 1 extra | 7 extra |
| PBMC Chromium X 6.1 | 1 extra (p 0.0093) | **identical** (7.2 refill: 1 extra) | 1 extra | 2 extra |
| PBMC 500 LT 6.1 | 31 missing (OrdMag only) | **identical** with LT ranks | — | identical |
| hgmm 10k 6.1 | 95 missing, 523 extra | across genomes 70 missing, 546 extra; per genome 11 missing, 1 extra (any per-genome combination: 4–11 missing, 0–2 extra, some missing with p = 1) | — | 17 missing, 1 extra |
| mouse brain nuclei 7.0 | 77 missing, 1 extra (3,000) | — | — | 2 missing, 5 extra (p 0.009–0.014); cap 262,144 the same |
| jejunum nuclei 7.0 | 1 missing, 57 extra (3,000) | — | — | 1 missing, 3 extra (p 0.007–0.010) |
| PBMC 5k 7.0.1 | identical with 5,000 | identical with 5,000 | identical with 5,000 | identical |
| Flex spleen 7.1 | 5 missing, 1,171 extra | — | — | 4 missing, 1 extra (p 0.009–0.011) |

Conclusions:

- **4.0**: 3.x's algorithm with the simulation on a fresh stream seeded 0. The neuron and parent
  datasets rule out 3.x's shared stream and 7.2's OrdMag; the refill and the candidate limit make
  no difference on them (fewer than 20,000 candidates). Implemented as `cellranger` `v4`, and
  `cellranger-ordmag` `v2` extended to 5.x.
- **5.0**: as 4.0, after the aggregate count rule (as 7.2's). The dataset cannot tell the stream
  apart (every combination agrees once the aggregates are gone); `v5` keeps 4.0's.
- **6.0–6.1**: the same rule reproduces five of seven single-genome datasets, the LT ones with
  LT's ambient ranks, so the LT range predates 7.2. 6.1 tests each genome on its own, but not in
  either form tried (some barcodes 10x called look ambient in their genome's test). Two
  single-genome datasets stay 1–3 barcodes off at the FDR cut. Not explained by:
  other tie orders for the ambient droplets (stable sorts either way), sorting over the whole
  whitelist or over barcodes with counts only, other `--expect-cells`, 7.2's refill, 7.2's OrdMag.
- **7.0–7.1**: closest to 7.2's rule (estimated cells, own stream), off by 4–7 barcodes at the FDR
  cut on three of four datasets; 7.0's larger cap does not change the estimates here.

6.0–7.1 remain a `Gap`.
