---
title: Developers
nav_order: 6
---

# Developers
{: .no_toc }

How the package is built and how a method is added. `AGENTS.md` in the repository says the same
for coding agents; the two are kept in step.

1. TOC
{:toc}

## Set up

```bash
git clone https://github.com/cellgeni/sc-droplets
cd sc-droplets
uv sync --extra dev
uv run sc-droplets --help
uv run pytest -m "not integration"
```

## Layout

| Path | Holds |
|---|---|
| `src/sc_droplets/registry.py` | `MethodInfo`, `Variant`, `Gap`, `Method`, `param`, `@register`, specs and tool versions |
| `src/sc_droplets/counts.py` | `Droplets`: what a method sees (counts of the chosen features, totals, reads) |
| `src/sc_droplets/result.py` | `CallResult` and how a call is stored on AnnData |
| `src/sc_droplets/api.py` | `call`, `compare`, `agreement`, `list_methods` |
| `src/sc_droplets/io/` | readers (10x HDF5, MTX folders, pipeline folders) and writers |
| `src/sc_droplets/methods/` | one module per tool page |
| `src/sc_droplets/primitives/` | each tool's arithmetic, shared: Good-Turing variants, numpy's classic quicksort, C++, R and Boost random numbers, R's `smooth.spline`, `lgamma`, long-double sums and statistics, the Monte Carlo engines |
| `oracles/` | scripts that run the original tools on public data (`.github/workflows/oracles.yml`); their output goes to `tests/oracles/`. `oracles/dropletutils/` runs every DropletUtils release from bioconda with pixi, `oracles/umitools/` UMI-tools' own whitelist functions, `oracles/kneedle/` kneed's `KneeLocator`, `oracles/cellsweep/` cellsweep's `infer_empty_droplets`, `oracles/cellbender/` CellBender's `SingleCellRNACountsDataset` from each release's tagged source, `oracles/bustools/` bustools' `whitelist` on a BUS file written from the matrix, `oracles/alevinfry/` programs built from alevin-fry's source at each tag, `oracles/alevin/` programs built from salmon's `Alevin.cpp` at each tag |
| `src/sc_droplets/pl/` | plots |
| `src/sc_droplets/cli.py` | the command; `call <id>` subcommands are built from the registry |
| `docs/` | this site; `docs/_research/` holds the source-reading notes the pages start from |
| `tests/` | the unit suite; `integration` tests download public data |

## Adding a method

Write the page first. The page is the specification, and for code under a licence that does not
allow copying (Cell Ranger 2.x, 7.2 and later), it is the only thing the implementation may be
written from.

1. **Page.** Add `docs/methods/<tool>/<algorithm>.md` with front matter `method_ids` and
   `status`, and the sections every method page has, in order: At a glance, How it works,
   Versions, Parameters, Run it, Differences from the original, Validation, References. Read the
   original's source at each release where it changed, and write down every constant, rounding
   rule and tie break.
2. **Code.** In `src/sc_droplets/methods/<tool>.py`, a `Method` subclass with a `MethodInfo`
   (one `Variant` per behaviour, each with the range of tool releases it covers, and a `Gap` for
   releases that cannot be reproduced), a `Params` dataclass built with `param`, and `run`.
   Decorate it with `@register` and add the module to `methods/__init__.py`.
3. **Oracle.** Check it against the original tool on the same matrix: exact for deterministic
   methods, within Monte Carlo noise for the rest. Correctness is never shown by the
   implementation agreeing with itself. Add the cases to `tests/test_validation.py`: public data
   goes in `tests/datasets.py`; a tool that has to be run goes in `oracles/`.
4. **Table.** Regenerate the table on the Methods page:

   ```bash
   uv run python docs/_scripts/methods_table.py
   ```

`tests/test_docs_are_accurate.py` fails until the page exists, lists the method's options, and
the table is current.

## Tests

```bash
uv run pytest -m "not integration"
uv run pytest tests/test_registry.py::test_specs_split_into_id_version_and_parameters
uv run pytest -m "not integration and not perf" --cov=sc_droplets --cov-fail-under=90
uv run python tests/datasets.py                       # download the public datasets (about 1 GB)
uv run pytest -m integration -n auto --dist worksteal tests/test_validation.py # every method against its original tool, on every core
```

CI runs the validation suite as two jobs (`--shard 1/2`, `--shard 2/2`), each on every core of its
runner. The `Oracles` workflow runs a tool's job only when that tool's scripts, or a helper they
use, changed (`.github/scripts/changed-oracles.sh`); started by hand, it runs them all.

House style: test names read as sentences; a regression test's docstring says what used to go
wrong; comments explain why.

## Previewing the site

```bash
cd docs
bundle install
bundle exec jekyll serve
```

## Branches

`dev` collects work in progress: each milestone or fix is a branch from `dev` and a pull request
into `dev`, merged once its checks pass. `main` holds releases only, and the site is built from
it, so documentation reaches the site with a release.

## Releases

1. A branch `release/X.Y.Z` from `dev`; the version in `pyproject.toml` and `uv.lock`.
2. A `## X.Y.Z` section in [Changelog](changelog.md).
3. A pull request into `main`. Merge, then push a tag `X.Y.Z` (no `v`), and merge `main` back
   into `dev`. The tag runs the suite, checks the tag against
   `pyproject.toml`, publishes to PyPI and pushes `quay.io/cellgeni/sc-droplets:X.Y.Z`.
