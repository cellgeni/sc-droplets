---
title: empty droplets
parent: DropletQC
grand_parent: Methods
nav_order: 1
method_ids: [dropletqc]
status: implemented
---

# DropletQC empty droplets
{: .no_toc }

DropletQC's `identify_empty_drops`. An empty droplet holds mostly spliced, cytoplasmic RNA, so its
nuclear fraction, the share of its reads in introns, is low. DropletQC estimates the derivative of
the density of the fractions, cuts at the valley after the first peak, and rescues the barcodes
below the cut whose fraction and UMI count both clear a floor. It flags the barcodes another caller
kept; it is not meant for a whole raw matrix.

1. TOC
{:toc}

## At a glance

| Id | `dropletqc` |
|---|---|
| Tool | DropletQC 1.0.0 (GitHub `powellgenomicslab/DropletQC`, tag `v1.0.0`), with the ks it runs on |
| Category | density |
| Counts | the nuclear fraction, and UMIs per barcode for the rescue |
| Needs expected cells | no |
| Deterministic | yes: no random numbers |
| Original licence | MIT (ks: GPL-2 or GPL-3) |
| Validated against | DropletQC's own `identify_empty_drops`, with ks 1.14.2 and 1.15.3, in R |

## How it works

DropletQC is given barcodes, each with a nuclear fraction $$f$$ and a UMI count $$u$$.

1. **The derivative.** `ks::kdde(f, deriv.order = 1)` with ks's defaults: a binned kernel estimate
   of the first derivative of the density of $$f$$, with a Gaussian kernel.
   - The bandwidth $$h$$ is `ks::hpi(f, nstage = 2, binned = TRUE, deriv.order = 1)`, Wand and
     Jones' two-stage plug-in for the first derivative. With $$\psi_r$$ the density functionals
     and $$\phi^{(r)}$$ the derivatives of the standard normal density at 0:
     $$\hat\psi_{10}$$ is the normal-scale value at the sample's standard deviation;
     $$g_8 = \big(2\phi^{(8)}(0) / (-\hat\psi_{10}\, n)\big)^{1/11}$$;
     $$\hat\psi_8$$ is the binned estimate with bandwidth $$g_8$$;
     $$g_6 = \big(2\phi^{(6)}(0) / (-\hat\psi_8\, n)\big)^{1/9}$$;
     $$\hat\psi_6$$ the binned estimate with $$g_6$$; and
     $$h = \big(3\,\psi_2^{NS} / (\hat\psi_6\, n)\big)^{1/7}$$, with $$\psi_2^{NS}$$ the
     standard normal's.
   - A binned estimate: 401 grid points from $$\min f - (4 + r)g$$ to $$\max f + (4 + r)g$$ for
     a functional of order $$r$$ (from $$\min f - 4.7h$$ to $$\max f + 4.7h$$ for the
     derivative), linear binning of the values onto the grid, and the counts convolved with the
     kernel's derivative sampled at the grid spacing out to $$(4 + r)$$ bandwidths, by FFT.
2. **The cut.** The derivative's sign along the grid, run-length encoded. The cut is the grid
   point at index (length of the first run + length of the second), 1-based: the last point of
   the second run. When the derivative starts positive, that is the valley after the first peak;
   when it starts at zero or below, it is the first peak. With one run there is no cut and no
   empty droplet; with fewer than 4, DropletQC warns that there may be no empty droplets.
3. **The call.** `empty_droplet` when $$f$$ is below the cut, `cell` otherwise; then `cell` again
   when $$f >$$ `nf_rescue` (0.05) and $$u >$$ `umi_rescue` (1,000), both strict.

### The nuclear fraction

DropletQC computes it from a BAM: per barcode, $$f = N/(N + E)$$ over the alignment records with
Cell Ranger's tag `RE = N` (intronic) and `RE = E` (exonic), reads not UMIs, duplicates and
multimappers included; 0/0 is 0 (`nuclear_fraction_tags`; `nuclear_fraction_annotation` uses a
GTF instead). A count matrix does not hold it. `--fraction` takes it as an `obs` column or a table
of barcode and fraction, for instance DropletQC's own output; by default (`auto`) the method uses
an `obs` column `nuclear_fraction`, or else the unspliced fraction of the package's splicing
readers (velocyto, STARsolo Velocyto, kb nac, alevin-fry USA; [Inputs](../../reference/inputs.md#splicing)),
which counts UMIs, not reads: an approximation, and the call says so.

### Which barcodes

DropletQC flags the barcodes it is given: in its paper those an earlier caller kept (EmptyDrops
with `lower` 500, Cell Ranger, CellBender, EmptyNN), in its tutorial Cell Ranger's filtered
barcodes. `--among` names a stored call (or any boolean `obs` column) whose cells are given;
without it, every barcode with at least `--min-counts` (500) UMIs is. The barcodes not given are
not cells (class `not given`).

## Versions

| DropletQC | Variant | `--tool-version` | Change to the call |
|---|---|---|---|
| 1.0.0 | `v1.0` | `1.0.0` | the only release |

The estimate depends on the installed ks, whose 1.14.3 changed two lines of the binned estimator
(`--ks-version`):

| ks | Change |
|---|---|
| 1.13.2 – 1.14.2 | the kernel reaches $$\lceil (4+r)h(M-1)/(b-a) \rceil$$ grid points; the counts sit at offset $$(N+1)/2$$ of the convolution |
| 1.14.3 – 1.15.3 | $$\lceil (4+r)h/\delta \rceil$$, the same number rounded on another path; the counts at $$\operatorname{round}(N/2) + 1$$, which with R's halves-to-even rounding is one cell later whenever $$(N-1)/2$$ is odd |

On 4,500 test values the bandwidth moves from 0.019629 (ks 1.14.2) to 0.019654 (1.15.3) and the
cut by one grid point.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--fraction` | `auto` | the first column of `nf_umi` | each barcode's nuclear fraction: an `obs` column, a table of barcode and fraction, or `auto` (`nuclear_fraction`, else `unspliced_fraction`) |
| `--among` | none | the barcodes in `nf_umi` | the barcodes given: a stored call's cells, or a boolean `obs` column |
| `--min-counts` | 500 | | without `--among`, the barcodes given have at least this many UMIs |
| `--nf-rescue` | 0.05 | `nf_rescue` | a rescued barcode has a larger fraction ... |
| `--umi-rescue` | 1000 | `umi_rescue` | ... and more UMIs |
| `--ks-version` | 1.15.3 | | the ks release whose estimator to follow; 1.13.2 or later |

## Run it

```bash
sc-droplets call dropletqc raw_feature_bc_matrix.h5 --splicing velocyto.loom --json
sc-droplets call cellranger raw_feature_bc_matrix.h5 -o calls.h5ad
sc-droplets call dropletqc calls.h5ad --among cellranger --fraction nuclear_fraction.csv --barcodes cells.tsv.gz
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
scd.add_metrics(adata, splicing="velocyto.loom")
scd.call(adata, "cellranger")
result = scd.call(adata, "dropletqc", among="cellranger")
result.thresholds["nuclear_fraction"], result.thresholds["bandwidth"]
adata.obs["dropletqc_class"].value_counts()
```

The classes are `cell`, `rescued` (below the cut, kept by the rescue), `empty_droplet` and `not
given`; `obs["dropletqc_nuclear_fraction"]` holds the fraction used. The thresholds hold the cut
(`nuclear_fraction`), the bandwidth, the rescue floors and the number of barcodes given.
`sc-droplets plot diagnostics` draws the derivative with the cut and the fraction against the
UMIs.

## Differences from the original

- The nuclear fraction is not computed from a BAM here; given as DropletQC computed it, the call
  is DropletQC's. The unspliced fraction is an approximation (UMIs, not reads; the reader's
  rules, not Cell Ranger's `RE` tags).
- The estimate is an FFT convolution in both, R's there and numpy's here: the bandwidth and the
  derivative agree to the last bits, not in them. A fraction would have to fall within those bits
  of the cut to be called differently.
- A barcode given without a fraction gets 0, as DropletQC's 0/0.
- `identify_damaged_cells` is not offered yet.

## Validation

`oracles/dropletqc/run.sh` sources DropletQC's `R/identify_empty_drops.R` at `v1.0.0` and runs it
in R from conda-forge, once with ks 1.14.2 and once with ks 1.15.3, on seven tables of nuclear
fractions and UMIs (`make_inputs.py`):

- 10x's pbmc_1k_v3 (Cell Ranger 3.0), with kb's nascent fraction: the 1,222 barcodes of 10x's
  filtered matrix (as DropletQC's paper uses it, after a caller), and the 1,241 with at least 500
  UMIs;
- simulated fractions: two peaks (two mixtures), one peak (DropletQC's warning, a cut at the last
  grid point), 25 barcodes, and fractions of small counts with many ties.

Every barcode gets DropletQC's status under both ks releases; the bandwidth and the cut agree to
$$4 \times 10^{-15}$$ relative. ks's change shows in the outputs: on the first simulated mixture
the bandwidth is the same under both and the cut moves one grid point (0.29186 and 0.29437); on
PBMC 1k's filtered barcodes the cut is 0.1239 with ks 1.14.2 and 0.1262 with 1.15.3.

## References

- Muskovic & Powell 2021, DropletQC: improved identification of empty droplets and damaged cells
  in single-cell RNA-seq data, *Genome Biology* 22:329,
  [doi:10.1186/s13059-021-02547-0](https://doi.org/10.1186/s13059-021-02547-0).
- Wand & Jones 1995, *Kernel Smoothing*, Chapman & Hall: the plug-in selector.
- Source: [`R/identify_empty_drops.R`](https://github.com/powellgenomicslab/DropletQC/blob/v1.0.0/R/identify_empty_drops.R)
  at `v1.0.0`; ks's `R/selector.R`, `R/kdde.R`, `R/kfe.R`, `R/normal.R`, `R/binning.R` and
  `src/ks.c` ([CRAN](https://cran.r-project.org/package=ks)).
