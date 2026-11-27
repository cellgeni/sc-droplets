---
title: Commands
parent: Reference
nav_order: 1
---

# Commands
{: .no_toc }

Results and data go to stdout; progress, notes and errors go to stderr. A failure prints one
line starting `Error:` and exits with code 1.

1. TOC
{:toc}

## `sc-droplets`

```bash
sc-droplets --version
sc-droplets --help
```

| Option | Meaning |
|---|---|
| `--version`, `-V` | Print the installed version |

## `methods`

```bash
sc-droplets methods
sc-droplets methods top-n
```

Without an argument, a table of every method. With a method id, its variants (the releases of
its tool and what each did), its parameters with defaults and the original tool's names for
them, and its references.

## `inspect`

```bash
sc-droplets inspect raw_feature_bc_matrix.h5
sc-droplets inspect Solo.out --solo-feature GeneFull --width 80
sc-droplets inspect raw_feature_bc_matrix.h5 --splicing velocyto.loom
```

A summary of a raw matrix (barcodes, features by type, totals at ranks 1 to 100,000, the median
unspliced and mitochondrial fractions when they can be had) and its knee drawn in the terminal.
Warns when the matrix looks filtered.

| Option | Meaning |
|---|---|
| `--feature-type` | Feature type to count (repeatable); `all` for every feature. Default: Gene Expression when present |
| `--layer` | Count this layer instead of `X` |
| `--genome` | Count only this genome's features |
| `--solo-feature` | STARsolo feature folder (`Gene`, `GeneFull`, `GeneFull_Ex50pAS`, `GeneFull_ExonOverIntron`, `Velocyto`) |
| `--splicing` | Where to read spliced and unspliced counts, matched by barcode (see [Inputs](inputs.md#splicing)); default: the input's own layers |
| `--width` | Width of the terminal plot in characters (default 60) |

## `call`

```bash
sc-droplets call top-n raw_feature_bc_matrix.h5 --n-cells 3000 -o calls.h5ad
sc-droplets call threshold raw_feature_bc_matrix.h5 --min-counts 500 --table calls.tsv --json
```

One subcommand per method; its own options are listed on the method's page and by
`sc-droplets call <id> --help`. Every subcommand also takes:

| Option | Meaning |
|---|---|
| `--output`, `-o` | Write the input with the call added, as `.h5ad` |
| `--cells-only` | With `-o`, write only the called cells |
| `--barcodes` | Write the cell barcodes, one per line (gzipped if the name ends in `.gz`) |
| `--table` | Write one row per barcode: total, call, class and scores (`.csv`, or `.tsv`) |
| `--json` | Print the call's summary as JSON on stdout |
| `--tool-version`, `-t` | Release of the original tool to reproduce (methods that follow a tool) |
| `--feature-type` | As for `inspect`; default: what the method counts |
| `--layer` | As for `inspect` |
| `--genome` | As for `inspect` |
| `--reads` | Read counts per barcode, for methods that rank reads: an obs column, or a CSV/TSV of barcode and count |
| `--seed` | Random seed; default: the seed the original tool uses |
| `--key` | Key to store the call under in `obs` and `uns` (default: the method id) |
| `--solo-feature` | As for `inspect` |
| `--splicing` | As for `inspect`: the unspliced fraction (and the mitochondrial fraction) go into `obs` of `-o`'s output |

Without `-o`, `--barcodes`, `--table` or `--json` nothing is kept, and the command says so.

## `compare`

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m threshold:min_counts=500
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m top-n:n_cells=2000 --json
sc-droplets compare raw_feature_bc_matrix.h5 -m otsu -m gmm --splicing velocyto.loom --table calls.csv
```

Runs several methods on one sample and prints how many cells each called and the smallest total
among them. Each `-m` is a method spec, `id[@tool_version][:name=value,...]`; each call is stored
under its own key (the id, the version and the spec's parameters).

With the unspliced fraction (the input's layers, or `--splicing`), each call also gets its
median in the call's cells, in the barcodes it left out with at least 10 gene-expression counts, and the AUROC
between the two; likewise the mitochondrial fraction when genes are named `MT-`. These separation
scores are unsupervised indicators of a call's quality, not ground truth.

| Option | Meaning |
|---|---|
| `--method`, `-m` | A method spec (repeatable, at least one) |
| `--output`, `-o` | Write the input with every call added, as `.h5ad` |
| `--table` | Per-barcode table: the total, one column per call, and the fractions measured |
| `--plot` | Knee plot of every call (`.png`, `.pdf`, `.svg`) |
| `--json` | Print the summary (with the separation scores) and the pairwise Jaccard index of the cell sets as JSON |
| `--feature-type`, `--layer`, `--genome`, `--reads`, `--solo-feature`, `--splicing` | As for `call` |

## `plot knee`

```bash
sc-droplets plot knee raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -o knee.png
sc-droplets plot knee calls.h5ad --key top_n -o knee.pdf --rotated
sc-droplets plot knee raw_feature_bc_matrix.h5 -m otsu --splicing velocyto.loom --color unspliced_fraction --track unspliced_fraction -o knee.png
```

The barcode-rank curve with each call's cut-off. Given an `.h5ad` with stored calls and neither
`-m` nor `--key`, it draws every stored call. `--color` colours each barcode by a measure;
`--track` adds a panel under the curve with the measure's median in log-spaced bins of rank, so a
change in composition at the knee shows.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--rotated` | Totals on x, ranks on y |
| `--cumulative` | Cumulative share of counts against rank |
| `--color` | Colour the barcodes by this obs column (`unspliced_fraction`, `mito_fraction`, ...) |
| `--track` | A panel with this obs column's median along rank (repeatable) |
| `--title` | Figure title |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature`, `--splicing` | As for `inspect` |

## `plot density`

```bash
sc-droplets plot density raw_feature_bc_matrix.h5 -m otsu -m kde-valley -o density.png
sc-droplets plot density calls.h5ad -o density.pdf --weighted
```

The histogram of $$\log_{10}$$ totals with each call's cut-off. Given an `.h5ad` with stored
calls and neither `-m` nor `--key`, it draws every stored call.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--min-counts` | Leave out barcodes below this total (1) |
| `--bins` | Histogram bins (100) |
| `--weighted` | Weight each barcode by its total |
| `--title` | Figure title |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature` | As for `inspect` |

## `plot compare`

```bash
sc-droplets plot compare raw_feature_bc_matrix.h5 -m cellranger -m otsu -m gmm -o compare.png
sc-droplets plot compare calls.h5ad -o compare.pdf
```

Several calls side by side. Left: the barcode-rank curve with each call's cut, and under it one
strip per call, shaded by the share of each log-spaced bin of rank the call took, so where the
calls part is plain. Right: the number of cells per call, and the pairwise Jaccard index of
their cell sets.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature` | As for `inspect` |

## `plot qc`

```bash
sc-droplets plot qc raw_feature_bc_matrix.h5 -m cellranger -m gmm -o qc.png
```

Genes detected against count, and the mitochondrial fraction against count when genes are named
`MT-`, each barcode coloured by how many calls took it. With calls that disagree, a third panel
compares the genes per count of the barcodes some calls took with those every call took and
those none took in the band of totals where the calls disagree.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--min-counts` | Leave out barcodes below this total (10) |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature` | As for `inspect` |

## `plot splicing`

```bash
sc-droplets plot splicing raw_feature_bc_matrix.h5 -m cellranger -m otsu --splicing velocyto.loom -o splicing.png
sc-droplets plot splicing counts_unfiltered -m gmm -m kde-valley -o splicing.pdf
```

The unspliced fraction against the calls: against log10 count with each barcode coloured by how
many calls took it; in each call's cells and in the barcodes it left out (with the AUROC between
them); and, with several calls, in the barcodes every call took, some took and none took. In
whole-cell samples empty droplets sit low, cells higher, damaged cells high.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--metric` | The obs column to draw (default `unspliced_fraction`; `mito_fraction`) |
| `--min-counts` | Leave out barcodes below this total (10): too few molecules for a fraction |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature`, `--splicing` | As for `inspect` |

## `plot diagnostics`

```bash
sc-droplets plot diagnostics raw_feature_bc_matrix.h5 -m starsolo-emptydrops-cr -o diagnostics.png
sc-droplets plot diagnostics calls.h5ad --key cellranger_3_0_0 -o diagnostics.pdf
```

One call's own view. For a threshold method: the rank curve, its cells, its cut and the robust
maximum the cut came from (with the spread of OrdMag's bootstraps). For a test against the ambient
profile: three panels, the classes on the rank curve with the ambient ranks shaded, the
candidates' log-likelihood against their total, and the histogram of their p-values.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure; required |
| `--method`, `-m` | Run this method spec and draw it |
| `--key` | Draw this stored call |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature` | As for `inspect` |

Exactly one of `--method` and `--key`.

## `report`

```bash
sc-droplets report raw_feature_bc_matrix.h5 -m cellranger -m otsu -m gmm -o report.html
sc-droplets report calls.h5ad --splicing velocyto.loom -o report.html --title "PBMC 4k"
```

One self-contained HTML page for a sample: a table of the calls (with the separation scores
when the unspliced or mitochondrial fraction can be had), the knee, density, agreement, quality
and splicing plots, and for each call its parameters, thresholds, notes and diagnostic view.
The figures are inline SVG (dense scatters as images inside them); nothing is loaded from
elsewhere, so the page can be mailed or archived. The page names the input file, not its path.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The page to write; required |
| `--method`, `-m` | Run this method spec and report it (repeatable) |
| `--key` | Report this stored call (repeatable); default: every stored call |
| `--title` | The page's title (default: the input's file name) |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature`, `--splicing` | As for `inspect` |
