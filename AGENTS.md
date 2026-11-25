# AGENTS.md

Guidance for coding agents working in this repository. `CLAUDE.md` imports this file. Edit this
one, not that.

`sc-droplets` collects every method for calling cells vs empty droplets in droplet single-cell
data, version by version: a review site (one page per algorithm), a Python package on AnnData,
and a command line. The point is fidelity: `cellranger-emptydrops@9.0` must call the cells Cell
Ranger 9.0 called on the same matrix. A method that only resembles its original is a bug.

Naming:
- PyPI distribution, repository and command: `sc-droplets`.
- Import package: `sc_droplets` (`import sc_droplets as scd`).

## Commands

```bash
uv sync --extra dev                                  # set up
uv run sc-droplets --help                            # run the CLI from source
uv run pytest -m "not integration"                   # unit suite; what CI gates on
uv run pytest tests/test_api.py::test_top_n_breaks_ties_by_input_order_or_keeps_them_all
uv run pytest -m "not integration and not perf" --cov=sc_droplets --cov-fail-under=90
uv run python docs/_scripts/methods_table.py         # regenerate the table on docs/methods/index.md
```

Markers: `integration` (downloads public data or builds oracle environments; needs the network),
`perf` (cost guards, their own CI job), `slow`.

## Architecture

- `registry.py`: `ToolVersion` (orders `2.7.10a`, `cellranger-7.2.0`, `v0.3.2`), `Variant` (a
  behaviour over tool releases [since, until)), `Gap` (releases that cannot be reproduced),
  `MethodInfo`, `Method` (base class: `info`, a `Params` dataclass built with `param`, `run`),
  `@register`, `parse_spec` (`id[@version][:name=value,...]`), `default_key`.
  - Parameter values resolve as: dataclass default < `variant.defaults` < what the caller gave.
    `None` from the caller means "not set".
  - A method parameter may not take a name in `RESERVED_PARAMS` (the CLI's own options).
- `counts.py`: `Droplets`, the only thing a method sees: CSR counts of the selected features
  (Gene Expression by default, as most tools count), float64 totals, optional reads, obs.
- `result.py`: `CallResult` and `store`: `obs[k]`, `obs[k + "_class"]`, `obs[k + "_" + col]`,
  `uns["sc_droplets"][k]`. Everything in `uns` must survive `write_h5ad` (`h5ad_safe`).
- `api.py`: `call`, `compare`, `agreement`, `list_methods`, `stored`.
- `io/`: `read` finds the raw matrix in any file or pipeline folder (`find_matrix`); `tenx.py`
  (10x HDF5 v2/v3), `mtx.py` (Cell Ranger, STARsolo, Velocyto, kb incl. nac, alevin-fry incl.
  USA); `write.py` stages every output beside its name and renames it into place.
- `methods/`: one module per tool page, listed in `methods/__init__.py::_MODULES`.
  `cellranger.py` holds the chemistry tables (ambient range and FDR per release and kit family),
  OrdMag v2/v7, the EmptyDrops steps, and the aggregate and high-occupancy filters.
- `primitives/`: the arithmetic of each tool, reused across methods: `goodturing` (Cell Ranger's
  and STAR's), `legacysort` (numpy's classic quicksort, with a `window` so a slice is cheap: tie
  order decides Cell Ranger's ambient droplets), `stdrandom` (C++ `mt19937` and libstdc++'s
  `discrete_distribution`), `simulate` (Monte Carlo engines: Cell Ranger 3.x–8.x with its block
  refills, 9.0+ PCG64, STAR), `multinomial`, `pvalues`, `rounding`.
- `pl/`: plots take an AnnData, draw stored calls or run specs, return the Axes, accept `ax=`.
  `terminal.py` draws the knee in braille for `inspect`.
- `cli.py`: typer. `call <id>` subcommands are generated from the registry
  (`_add_call_command`): options from `Params`, plus the shared ones. Results and data to
  stdout; messages to stderr (`console`); a failure is one `Error:` line and exit 1 (`_fail`).
  Help text is Rich markup: escape `[` (`rich.markup.escape`, or `\\[` in literals).

## Rules the code depends on

- **The page is the spec.** Write `docs/methods/<tool>/<algorithm>.md` before the code. For
  Cell Ranger 2.x (AGPLv2) and 7.2+ (10x licence), the implementation is written from the page
  only; never copy, paste or transliterate their code. GPL, MIT and BSD sources (DropletUtils,
  STAR, UMI-tools, CellBender, kneed, Cell Ranger 3.0.2) may be ported closely, keeping their
  notices. `docs/_research/` holds the source-reading notes the pages start from.
- **Every behaviour change of the original is a `Variant`.** Ranges are half-open
  [since, until). Releases without public source are a `Gap`, which `resolve` refuses.
- **Reproduce the original's arithmetic.** Its RNG and seed (`MethodInfo.default_seed`), sort
  stability and tie order, rounding (Python 2 `round` vs `np.round`, R's half-to-even), and
  off-by-ones. Document each on the page under "Differences from the original" when it cannot
  be reproduced.
- **Correctness is checked against an oracle**: the original tool's output on the same matrix
  (exact for deterministic methods; for Monte Carlo methods ≥ 99.5 % of called barcodes agree and
  every disagreement is within the noise band). Self-consistency alone is not accepted.
- **Inputs are raw.** Methods need the empty droplets; `looks_filtered` warns otherwise.
- anndata 0.13 lists `X` as the layer `None` in `layers.keys()`: filter it out.

## Tests

- `tests/conftest.py`: `sim` (session) / `adata` (a copy) are a small simulated raw matrix;
  writers for every input layout (`write_10x_h5_v3`, `write_10x_h5_v2`, `write_10x_mtx`,
  `write_starsolo`, `write_kb`, `write_alevin_fry`).
- Readers are checked against scanpy (`scanpy.read_10x_h5`, `read_10x_mtx`).
- `tests/test_validation.py` (marker `integration`) is the oracle suite: Cell Ranger against 10x's
  own filtered matrices (`tests/datasets.py` downloads and caches them in `$SC_DROPLETS_DATA`),
  STARsolo against STAR's calls in `tests/oracles/starsolo/` (made by `oracles/starsolo/run.sh`
  in the `Oracles` workflow, which runs on branch pushes that touch `oracles/`). A new method
  adds its cases here.
- CLI tests use a module-level `CliRunner`; strip ANSI codes and collapse whitespace before
  matching (`_out`).
- `tests/test_docs_are_accurate.py`: every `sc-droplets …` line in README and `docs/` parses
  against the command tree; every `scd.` name exists; the API reference's signatures match;
  every registered method has its page (front matter `method_ids`, `status`) with the sections
  in order and every option listed; the methods table is current; pages have `title`,
  `nav_order` and an existing `parent`; relative links and anchors exist.
- House style: test names read as sentences; a regression test's docstring says what used to go
  wrong; comments explain why.

## Documentation

Just the Docs (`just-the-docs/just-the-docs@v0.10.1`), built by GitHub Pages from `docs/` on
`main`; <https://cellgeni.github.io/sc-droplets/>. `docs/developers.md` covers the same ground as
this file, for people; keep the two in step.

- Every page needs `title` and `nav_order`, and `parent` (plus `grand_parent` under a tool) to
  appear in the sidebar.
- Method pages: front matter `method_ids: [id, ...]` and `status: implemented | planned |
  docs-only`; sections in order: At a glance, How it works, Versions, Parameters, Run it,
  Differences from the original, Validation, References.
- Formulas: `$$...$$`, inline and display alike (kramdown turns them into `\(..\)` / `\[..\]`, which
  KaTeX renders). A literal `\(` in a page shows as "("; the docs test refuses it.
- Link between pages as relative `page.md`.
- Preview: `cd docs && bundle install && bundle exec jekyll serve`.

## Status

`STATUS.md` records what is released, what is in progress (with its decisions and steps),
follow-ups noticed along the way, and the next candidates. Update it as work moves.

## Git and releases

- `dev` collects work in progress; `main` holds releases only. Each milestone or fix is a branch
  from `dev` and a pull request into `dev`, merged (by the agent too) once its checks pass and
  its review comments are dealt with.
- GitHub Pages builds the site from `main`, so documentation changes reach the site with a
  release.
- Release: branch `release/X.Y.Z` from `dev`; version in `pyproject.toml` and `uv.lock`; a
  `## X.Y.Z` section in `docs/changelog.md`; a pull request into `main`; merge; push tag `X.Y.Z`
  (no `v`); merge `main` back into `dev`. `publish.yml` checks the tag
  against `pyproject.toml`, runs `tests.yml`, publishes to PyPI (trusted publishing), then
  `quay-on-tag.yml` builds, smoke-tests and pushes `quay.io/cellgeni/sc-droplets`.
- Workflows pin third-party actions by commit SHA with the version in a comment, and pass tags
  and inputs to scripts through `env:`.
