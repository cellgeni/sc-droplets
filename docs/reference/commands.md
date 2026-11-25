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
```

A summary of a raw matrix (barcodes, features by type, totals at ranks 1 to 100,000, whether
spliced/unspliced counts are present) and its knee drawn in the terminal. Warns when the matrix
looks filtered.

| Option | Meaning |
|---|---|
| `--feature-type` | Feature type to count (repeatable); `all` for every feature. Default: Gene Expression when present |
| `--layer` | Count this layer instead of `X` |
| `--genome` | Count only this genome's features |
| `--solo-feature` | STARsolo feature folder (`Gene`, `GeneFull`, `GeneFull_Ex50pAS`, `GeneFull_ExonOverIntron`, `Velocyto`) |
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

Without `-o`, `--barcodes`, `--table` or `--json` nothing is kept, and the command says so.

## `compare`

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m threshold:min_counts=500
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m top-n:n_cells=2000 --json
```

Runs several methods on one sample and prints how many cells each called and the smallest total
among them. Each `-m` is a method spec, `id[@tool_version][:name=value,...]`; each call is stored
under its own key (the id, the version and the spec's parameters).

| Option | Meaning |
|---|---|
| `--method`, `-m` | A method spec (repeatable, at least one) |
| `--output`, `-o` | Write the input with every call added, as `.h5ad` |
| `--table` | Per-barcode table: the total and one column per call |
| `--plot` | Knee plot of every call (`.png`, `.pdf`, `.svg`) |
| `--json` | Print the summary and the pairwise Jaccard index of the cell sets as JSON |
| `--feature-type`, `--layer`, `--genome`, `--reads`, `--solo-feature` | As for `call` |

## `plot knee`

```bash
sc-droplets plot knee raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -o knee.png
sc-droplets plot knee calls.h5ad --key top_n -o knee.pdf --rotated
```

The barcode-rank curve with each call's cut-off. Given an `.h5ad` with stored calls and neither
`-m` nor `--key`, it draws every stored call.

| Option | Meaning |
|---|---|
| `--output`, `-o` | The figure (`.png`, `.pdf`, `.svg`); required |
| `--method`, `-m` | Run this method spec and draw it (repeatable) |
| `--key` | Draw this stored call (repeatable) |
| `--rotated` | Totals on x, ranks on y |
| `--cumulative` | Cumulative share of counts against rank |
| `--title` | Figure title |
| `--feature-type`, `--layer`, `--genome`, `--solo-feature` | As for `inspect` |
