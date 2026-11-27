"""The `sc-droplets` command.

Conventions (as adata-cli): results and data go to stdout, everything else
(progress, notes, errors) to stderr, so stdout pipes cleanly. A failure is one
`Error:` line on stderr and exit code 1.

`sc-droplets call <id>` has one subcommand per registered method, built from
the method's `Params`: adding a method to the registry adds its command.
"""

from __future__ import annotations

import enum
import inspect
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Type

import typer
from rich.console import Console
from rich.markup import escape
from rich.table import Table

from sc_droplets import registry

app = typer.Typer(
    help="Call cells vs empty droplets with any published method, as any version of its tool ran it.",
    rich_markup_mode="rich",
    no_args_is_help=True,
)
call_app = typer.Typer(
    help="Run one method. Each method is a subcommand with its own options; "
    "[bold]sc-droplets methods[/] lists them.",
    rich_markup_mode="rich",
    no_args_is_help=True,
)
plot_app = typer.Typer(help="Draw calls.", rich_markup_mode="rich", no_args_is_help=True)
app.add_typer(call_app, name="call")
app.add_typer(plot_app, name="plot")

console = Console(stderr=True, soft_wrap=True)
out_console = Console(soft_wrap=True)


def _fail(error: BaseException) -> typer.Exit:
    message = error.args[0] if isinstance(error, KeyError) and error.args else error
    console.print(f"[bold red]Error:[/] {escape(str(message))}")
    return typer.Exit(code=1)


def _version_callback(value: bool) -> None:
    if value:
        from sc_droplets import __version__

        out_console.print(__version__, highlight=False)
        raise typer.Exit()


@app.callback()
def _root(
    version: bool = typer.Option(
        False, "--version", "-V", help="Show the installed version and exit.", callback=_version_callback, is_eager=True
    ),
) -> None:
    """Call cells vs empty droplets with any published method."""


# ---------------------------------------------------------------------------
# Shared options and helpers
# ---------------------------------------------------------------------------


def _input_argument():
    return typer.Argument(
        ...,
        help="A raw (unfiltered) matrix: .h5ad, .zarr, 10x .h5, an MTX folder, or a "
        "Cell Ranger outs/, STARsolo Solo.out/, kb or alevin-fry output folder.",
        exists=True,
        dir_okay=True,
        file_okay=True,
        readable=True,
    )


def _feature_type_option():
    return typer.Option(
        None,
        "--feature-type",
        help="Feature type to count (repeatable), as in var['feature_types']; 'all' counts every "
        "feature. Default: what the method counts, usually Gene Expression.",
    )


def _layer_option():
    return typer.Option(None, "--layer", help="Count this layer instead of X (.h5ad input).")


def _genome_option():
    return typer.Option(None, "--genome", help="Count only this genome's features (barnyard).")


def _solo_feature_option():
    return typer.Option("Gene", "--solo-feature", help="STARsolo feature folder to read from Solo.out/.")


def _reads_option():
    return typer.Option(
        None,
        "--reads",
        help="Read counts per barcode, for methods that rank reads: an obs column, or a CSV/TSV "
        "of barcode and count.",
    )


def _splicing_option():
    return typer.Option(
        None,
        "--splicing",
        help="Where to read spliced and unspliced counts, matched by barcode: a velocyto .loom, "
        "STARsolo Velocyto/raw or Solo.out (GeneFull − Gene without Velocyto), kb --workflow nac, "
        "alevin-fry USA, an .h5ad with spliced and unspliced layers. Default: the input's own layers.",
    )


def _add_metrics(adata, splicing: Optional[Path]) -> List[str]:
    """The unspliced and mitochondrial fractions in obs, from --splicing or the input's layers."""
    from sc_droplets.metrics import add_metrics

    with console.status("Reading the spliced and unspliced counts" if splicing else "Measuring the barcodes"):
        added = add_metrics(adata, splicing)
    if added:
        console.print(f"Added {', '.join(added)} to obs")
    return added


def _load(path: Path, solo_feature: str = "Gene", genome: Optional[str] = None):
    from sc_droplets.io import read

    with console.status(f"Reading {escape(str(path))}"):
        adata = read(path, solo_feature=solo_feature, genome=genome)
    fmt = adata.uns.get("sc_droplets_input", {}).get("format", Path(path).suffix.lstrip(".") or "folder")
    console.print(f"Read {adata.n_obs:,} barcodes × {adata.n_vars:,} features ({escape(str(fmt))})")
    return adata


def _attach_reads(adata, reads: Optional[str]) -> Optional[str]:
    """The obs column holding read counts: as named, or loaded from a table."""
    if reads is None:
        return None
    if reads in adata.obs.columns:
        return reads
    path = Path(reads)
    if not path.exists():
        raise ValueError(f"--reads {reads}: neither an obs column nor a file")
    import pandas as pd

    sep = "\t" if ".tsv" in path.name or ".txt" in path.name else ","
    table = pd.read_csv(path, sep=sep, header=None, comment="#")
    if table.shape[1] < 2:
        raise ValueError(f"{path}: expected two columns, barcode and read count")
    if not str(table.iloc[0, 1]).replace(".", "", 1).isdigit():
        table = table.iloc[1:]
    counts = pd.Series(table.iloc[:, 1].astype(float).to_numpy(), index=table.iloc[:, 0].astype(str))
    adata.obs["n_reads"] = counts.reindex(adata.obs_names).fillna(0.0).to_numpy()
    return "n_reads"


def _feature_types(values: Optional[List[str]]):
    if not values:
        return None
    return "all" if values == ["all"] else list(values)


def _report_call(result) -> None:
    from sc_droplets.counts import looks_filtered

    version = f" as {result.tool_version}" if result.tool_version else ""
    console.print(
        f"[bold]{escape(result.method)}[/]{escape(version)} (variant {escape(result.variant)}): "
        f"{result.n_cells:,} cells of {len(result.barcodes):,} barcodes"
    )
    for name, value in result.thresholds.items():
        console.print(f"  threshold ({escape(name)}): {value:g}")
    reason = looks_filtered(result.totals)
    if reason:
        console.print(f"[yellow]Warning:[/] {escape(reason)}")
    for note in result.notes:
        console.print(f"[yellow]Note:[/] {escape(note)}")


def _write(adata, result, output, cells_only, barcodes, table, inputs) -> None:
    from sc_droplets.io.write import check_outputs, write_barcodes, write_h5ad, write_table

    check_outputs([output, barcodes, table], inputs)
    if barcodes is not None:
        write_barcodes(result.barcodes[result.is_cell], barcodes)
        console.print(f"Wrote {result.n_cells:,} cell barcodes to {escape(str(barcodes))}")
    if table is not None:
        write_table(result.table, table)
        console.print(f"Wrote the per-barcode table to {escape(str(table))}")
    if output is not None:
        target = adata[adata.obs[result.key].to_numpy()].copy() if cells_only else adata
        write_h5ad(target, output)
        console.print(f"Wrote {target.n_obs:,} barcodes to {escape(str(output))}")


# ---------------------------------------------------------------------------
# methods
# ---------------------------------------------------------------------------


@app.command("methods")
def methods_command(
    method: Optional[str] = typer.Argument(None, help="A method id, to show its versions and parameters."),
) -> None:
    """List every method, or describe one.

    Examples:
        sc-droplets methods
        sc-droplets methods top-n
    """
    try:
        if method is None:
            table = Table(title="Methods", show_lines=False)
            for column in ("id", "tool", "title", "category", "versions"):
                table.add_column(column)
            for cls in registry.all_methods():
                info = cls.info
                versions = ", ".join(v.name for v in info.variants) if info.versioned else "—"
                table.add_row(info.id, info.tool, info.title, info.category, versions)
            out_console.print(table)
            return
        cls = registry.get(registry.parse_spec(method).id)
        _describe(cls)
    except Exception as e:
        raise _fail(e)


def _describe(cls: Type[registry.Method]) -> None:
    info = cls.info
    out_console.print(f"[bold]{escape(info.id)}[/]: {escape(info.title)} ({escape(info.tool)})")
    out_console.print(escape(info.summary))
    out_console.print(f"Category: {info.category}; counts {info.counts}; "
                      f"{'deterministic' if info.deterministic else 'Monte Carlo'}")
    out_console.print(f"Docs: {info.docs_url}")
    if info.versioned:
        table = Table(title="Variants")
        for column in ("variant", "tool versions", "what it does"):
            table.add_column(column)
        for v in info.variants:
            table.add_row(v.name, v.versions, v.summary)
        for gap in info.gaps:
            table.add_row("(gap)", f"≥{gap.since}, <{gap.until}", gap.reason)
        out_console.print(table)
    specs = cls.param_specs()
    if specs:
        table = Table(title="Parameters")
        for column in ("option", "default", "original", "meaning"):
            table.add_column(column)
        for spec in specs:
            table.add_row(spec.flag, _default_text(cls, spec), spec.original or "", spec.help)
        out_console.print(table)
    for ref in info.references:
        out_console.print(f"Reference: {escape(ref.text)}" + (f" {ref.link}" if ref.link else ""))


def _default_text(cls: Type[registry.Method], spec: registry.ParamSpec) -> str:
    """The default as help shows it, per variant where variants differ."""
    per_variant = {v.name: v.defaults.get(spec.name, spec.default) for v in cls.info.variants}
    values = set(map(repr, per_variant.values()))
    if len(values) == 1:
        value = next(iter(per_variant.values()))
        return "required" if value is None and spec.required else str(value)
    return "; ".join(f"{name}: {value}" for name, value in per_variant.items())


# ---------------------------------------------------------------------------
# inspect
# ---------------------------------------------------------------------------


@app.command("inspect")
def inspect_command(
    input: Path = _input_argument(),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    solo_feature: str = _solo_feature_option(),
    splicing: Optional[Path] = _splicing_option(),
    width: int = typer.Option(60, "--width", min=20, help="Width of the terminal knee plot, in characters."),
) -> None:
    """Summarise a raw matrix and draw its knee in the terminal.

    Examples:
        sc-droplets inspect raw_feature_bc_matrix.h5
        sc-droplets inspect Solo.out --solo-feature GeneFull
        sc-droplets inspect raw_feature_bc_matrix.h5 --splicing velocyto.loom
    """
    try:
        import numpy as np

        from sc_droplets.counts import Droplets, feature_type_column, looks_filtered
        from sc_droplets.pl.terminal import knee_text

        adata = _load(input, solo_feature, genome)
        types = _feature_types(feature_type)
        if types is None:
            column = feature_type_column(adata.var)
            types = "Gene Expression" if column is None or "Gene Expression" in set(adata.var[column]) else "all"
        droplets = Droplets.from_anndata(adata, layer=layer, feature_types=types, genome=genome)
        totals = droplets.totals
        nonzero = totals[totals > 0]
        ordered = np.sort(nonzero)[::-1]
        column = feature_type_column(adata.var)
        table = Table(title=f"{escape(str(input))}", show_header=False)
        table.add_column("what")
        table.add_column("value")
        table.add_row("barcodes", f"{adata.n_obs:,} ({nonzero.size:,} with any count)")
        if column is not None:
            kinds = adata.var[column].astype(str).value_counts()
            table.add_row("features", ", ".join(f"{escape(k)} {n:,}" for k, n in kinds.items()))
        else:
            table.add_row("features", f"{adata.n_vars:,}")
        counted = "every feature" if droplets.feature_types is None else ", ".join(droplets.feature_types)
        table.add_row("counted", escape(counted) + (f" (layer {escape(layer)})" if layer else ""))
        table.add_row("total counts", f"{totals.sum():,.0f}")
        for rank in (1, 10, 100, 1000, 10000, 100000):
            if rank <= ordered.size:
                table.add_row(f"total at rank {rank:,}", f"{ordered[rank - 1]:,.0f}")
        from sc_droplets.metrics import MIN_COUNTS, MITO, UNSPLICED

        added = _add_metrics(adata, splicing)
        measured = totals >= MIN_COUNTS
        for column, label in ((UNSPLICED, "unspliced fraction"), (MITO, "mitochondrial fraction")):
            if column in added:
                values = adata.obs[column].to_numpy(dtype=float)[measured]
                values = values[~np.isnan(values)]
                where = f"median {np.median(values):.3f} of {values.size:,} barcodes with ≥ {MIN_COUNTS} counts" if values.size else "no barcode measured"
                table.add_row(label, where)
            elif column == UNSPLICED:
                table.add_row(label, "none: no spliced/unspliced layers (see --splicing)")
        out_console.print(table)
        for line in knee_text(totals, width=width):
            out_console.print(line, highlight=False)
        reason = looks_filtered(totals)
        if reason:
            console.print(f"[yellow]Warning:[/] {escape(reason)}")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


# ---------------------------------------------------------------------------
# call <method>
# ---------------------------------------------------------------------------

_COMMON = (
    "input",
    "output",
    "cells_only",
    "barcodes",
    "table",
    "json_out",
    "feature_type",
    "layer",
    "genome",
    "reads",
    "seed",
    "key",
    "solo_feature",
    "splicing",
    "tool_version",
)


def _common_parameters(versioned: bool) -> List[inspect.Parameter]:
    P = inspect.Parameter
    kind = P.POSITIONAL_OR_KEYWORD
    params = [
        P("input", kind, default=_input_argument(), annotation=Path),
        P(
            "output",
            kind,
            default=typer.Option(None, "--output", "-o", help="Write the input, with the call added, as .h5ad."),
            annotation=Optional[Path],
        ),
        P(
            "cells_only",
            kind,
            default=typer.Option(False, "--cells-only", help="With -o, write only the called cells."),
            annotation=bool,
        ),
        P(
            "barcodes",
            kind,
            default=typer.Option(None, "--barcodes", help="Write the cell barcodes, one per line (.gz to compress)."),
            annotation=Optional[Path],
        ),
        P(
            "table",
            kind,
            default=typer.Option(None, "--table", help="Write a per-barcode table: total, call, scores (.csv or .tsv)."),
            annotation=Optional[Path],
        ),
        P(
            "json_out",
            kind,
            default=typer.Option(False, "--json", help="Print the call's summary as JSON on stdout."),
            annotation=bool,
        ),
        P("feature_type", kind, default=_feature_type_option(), annotation=Optional[List[str]]),
        P("layer", kind, default=_layer_option(), annotation=Optional[str]),
        P("genome", kind, default=_genome_option(), annotation=Optional[str]),
        P("reads", kind, default=_reads_option(), annotation=Optional[str]),
        P(
            "seed",
            kind,
            default=typer.Option(None, "--seed", help="Random seed; default: the one the original tool uses."),
            annotation=Optional[int],
        ),
        P(
            "key",
            kind,
            default=typer.Option(None, "--key", help="obs/uns key to store the call under (default: the method id)."),
            annotation=Optional[str],
        ),
        P("solo_feature", kind, default=_solo_feature_option(), annotation=str),
        P("splicing", kind, default=_splicing_option(), annotation=Optional[Path]),
    ]
    if versioned:
        params.append(
            P(
                "tool_version",
                kind,
                default=typer.Option(
                    None,
                    "--tool-version",
                    "-t",
                    help="Release of the original tool to reproduce (default: the latest known).",
                ),
                annotation=Optional[str],
            )
        )
    return params


def _method_parameters(cls: Type[registry.Method]) -> List[inspect.Parameter]:
    params = []
    for spec in cls.param_specs():
        help_text = spec.help
        default = _default_text(cls, spec)
        help_text += f" [default: {default}]" if default != "required" else " [required]"
        if spec.original:
            help_text += f" (original: {spec.original})"
        if spec.choices:
            choice = enum.Enum(f"{cls.__name__}_{spec.name}", [(c, c) for c in spec.choices], type=str)
            annotation: Any = Optional[choice]
            option = typer.Option(None, spec.flag, help=escape(help_text), show_default=False)
        elif spec.type is bool:
            annotation = Optional[bool]
            name = spec.name.replace("_", "-")
            option = typer.Option(None, f"--{name}/--no-{name}", help=escape(help_text), show_default=False)
        else:
            annotation = Optional[spec.type]
            option = typer.Option(None, spec.flag, help=escape(help_text), show_default=False)
        params.append(inspect.Parameter(spec.name, inspect.Parameter.POSITIONAL_OR_KEYWORD, default=option, annotation=annotation))
    return params


def _add_call_command(cls: Type[registry.Method]) -> None:
    info = cls.info

    def command(**kwargs: Any) -> None:
        common = {name: kwargs.pop(name, None) for name in _COMMON}
        _run_call(cls, common, kwargs)

    command.__signature__ = inspect.Signature(_common_parameters(info.versioned) + _method_parameters(cls))
    command.__name__ = "call_" + info.id.replace("-", "_").replace(".", "_")
    lines = [f"{info.title} ({info.tool}). {info.summary}"]
    if info.versioned:
        lines.append("Versions: " + "; ".join(f"{v.name} = {v.versions}" for v in info.variants) + ".")
    lines.append(f"Docs: {info.docs_url}")
    call_app.command(info.id, help="\n\n".join(lines), short_help=info.title)(command)


def _run_call(cls: Type[registry.Method], common: Dict[str, Any], given: Dict[str, Any]) -> None:
    from sc_droplets.api import call

    try:
        adata = _load(common["input"], common["solo_feature"], common["genome"])
        if common["splicing"] is not None:
            _add_metrics(adata, common["splicing"])
        reads = _attach_reads(adata, common["reads"])
        result = call(
            adata,
            cls.info.id,
            tool_version=common.get("tool_version"),
            key_added=common["key"],
            layer=common["layer"],
            feature_types=_feature_types(common["feature_type"]),
            genome=common["genome"],
            reads=reads,
            random_state=common["seed"],
            **{k: v for k, v in given.items() if v is not None},
        )
        _report_call(result)
        _write(adata, result, common["output"], common["cells_only"], common["barcodes"], common["table"], [common["input"]])
        if common["json_out"]:
            out_console.print(result.to_json(), highlight=False)
        elif not any(common[k] for k in ("output", "barcodes", "table")):
            console.print("Nothing written: pass -o, --barcodes, --table or --json to keep the call.")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


for _cls in registry.all_methods():
    _add_call_command(_cls)


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------


@app.command("compare")
def compare_command(
    input: Path = _input_argument(),
    method: List[str] = typer.Option(
        ...,
        "--method",
        "-m",
        help="A method spec, id\\[@tool_version]\\[:name=value,...], e.g. top-n:n_cells=3000 (repeatable).",
    ),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write the input with every call added, as .h5ad."),
    table: Optional[Path] = typer.Option(None, "--table", help="Per-barcode table: total and one column per call."),
    plot: Optional[Path] = typer.Option(None, "--plot", help="Knee plot of every call (.png, .pdf or .svg)."),
    json_out: bool = typer.Option(False, "--json", help="Print the summary as JSON on stdout."),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    reads: Optional[str] = _reads_option(),
    solo_feature: str = _solo_feature_option(),
    splicing: Optional[Path] = _splicing_option(),
) -> None:
    """Run several methods on one sample and compare their calls.

    Parameters go after a colon in each spec; every call is stored under its
    own key. With splicing (the input's layers or --splicing), each call also
    gets the median unspliced fraction in its cells and in the barcodes it
    left out, and the AUROC between them; likewise the mitochondrial fraction.

    Examples:
        sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m threshold:min_counts=500
        sc-droplets compare raw_feature_bc_matrix.h5 -m otsu -m gmm --splicing velocyto.loom --table calls.csv
    """
    from sc_droplets.api import agreement, compare
    from sc_droplets.io.write import check_outputs, write_h5ad, write_table

    try:
        check_outputs([output, table, plot], [input])
        adata = _load(input, solo_feature, genome)
        metrics = _add_metrics(adata, splicing)
        reads_column = _attach_reads(adata, reads)
        summary = compare(
            adata,
            method,
            layer=layer,
            feature_types=_feature_types(feature_type),
            genome=genome,
            reads=reads_column,
        )
        jaccard = agreement(adata, list(summary.index))
        if json_out:
            out_console.print(
                json.dumps(
                    {"calls": json.loads(summary.to_json(orient="index")), "jaccard": json.loads(jaccard.to_json())},
                    indent=2,
                ),
                highlight=False,
            )
        else:
            grid = Table(title="Calls")
            scored = [m for m in metrics if f"{m}_auroc" in summary.columns]
            for column in ("key", "variant", "cells", "smallest cell total"):
                grid.add_column(column)
            for metric in scored:
                short = "unspliced" if metric.startswith("unspliced") else "mito"
                grid.add_column(f"{short}: cells / others (AUROC)")
            for key, row in summary.iterrows():
                low = "—" if row["min_cell_total"] is None else f"{row['min_cell_total']:g}"
                scores = [f"{row[f'{m}_cells']:.3f} / {row[f'{m}_others']:.3f} ({row[f'{m}_auroc']:.2f})" for m in scored]
                grid.add_row(str(key), str(row["variant"]), f"{row['n_cells']:,}", low, *scores)
            out_console.print(grid)
        if table is not None:
            from sc_droplets.pl._data import curve_for

            frame = adata.obs[list(summary.index) + metrics].copy()
            # The totals the calls were made on (their layer and feature types).
            curve = curve_for(adata, list(summary.index), layer=layer, feature_types=_feature_types(feature_type))
            frame.insert(0, "total", curve.totals)
            frame.index.name = "barcode"
            write_table(frame, table)
            console.print(f"Wrote the per-barcode table to {escape(str(table))}")
        if plot is not None:
            _save_knee(adata, list(summary.index), plot)
        if output is not None:
            write_h5ad(adata, output)
            console.print(f"Wrote {escape(str(output))}")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


def _save_knee(adata, keys: Sequence[str], path: Path, **kwargs: Any) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from sc_droplets.pl import knee

    ax = knee(adata, list(keys), **kwargs)
    if kwargs.get("color") or kwargs.get("track"):
        # A colour bar shared by the curve and its tracks defeats tight_layout; trim on saving instead.
        ax.figure.savefig(path, dpi=150, bbox_inches="tight")
    else:
        ax.figure.tight_layout()
        ax.figure.savefig(path, dpi=150)
    plt.close(ax.figure)
    console.print(f"Wrote {escape(str(path))}")


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------


def _plot_keys(adata, input: Path, key, method, layer, types, genome) -> List[str]:
    """The calls a plot draws: stored keys asked for, then method specs run; every stored call when
    neither is given."""
    keys: List[str] = list(key or [])
    stored = adata.uns.get("sc_droplets", {})
    missing = [k for k in keys if k not in stored]
    if missing:
        raise KeyError(f"no stored call {', '.join(missing)} in {input}; stored: {', '.join(stored) or 'none'}")
    if method:
        from sc_droplets.api import call

        for spec in method:
            keys.append(call(adata, spec, layer=layer, feature_types=types, genome=genome).key)
    return keys or list(stored)


@plot_app.command("knee")
def plot_knee(
    input: Path = _input_argument(),
    output: Path = typer.Option(..., "--output", "-o", help="Figure to write (.png, .pdf or .svg)."),
    method: Optional[List[str]] = typer.Option(
        None, "--method", "-m", help="Run this method spec and draw it (repeatable); see compare."
    ),
    key: Optional[List[str]] = typer.Option(
        None, "--key", help="Draw this call stored in an .h5ad input (repeatable). Default: every stored call."
    ),
    rotated: bool = typer.Option(False, "--rotated", help="Totals on x, ranks on y."),
    cumulative: bool = typer.Option(False, "--cumulative", help="Cumulative share of counts against rank."),
    color: Optional[str] = typer.Option(
        None, "--color", help="Colour the barcodes by this obs column (unspliced_fraction, mito_fraction...)."
    ),
    track: Optional[List[str]] = typer.Option(
        None, "--track", help="A panel under the curve with this obs column's median along rank (repeatable)."
    ),
    title: Optional[str] = typer.Option(None, "--title", help="Figure title."),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    solo_feature: str = _solo_feature_option(),
    splicing: Optional[Path] = _splicing_option(),
) -> None:
    """Barcode-rank (knee) plot with each call's cut-off.

    Examples:
        sc-droplets plot knee raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -o knee.png
        sc-droplets plot knee calls.h5ad -o knee.pdf --rotated
        sc-droplets plot knee raw_feature_bc_matrix.h5 -m otsu --splicing velocyto.loom --color unspliced_fraction --track unspliced_fraction -o knee.png
    """
    from sc_droplets.io.write import check_outputs

    try:
        check_outputs([output], [input])
        adata = _load(input, solo_feature, genome)
        if color or track or splicing:
            _add_metrics(adata, splicing)
        types = _feature_types(feature_type)
        keys = _plot_keys(adata, input, key, method, layer, types, genome)
        _save_knee(
            adata, keys, output, rotated=rotated, cumulative=cumulative, title=title, layer=layer, feature_types=types,
            color=color, track=list(track or []) or None,
        )
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


@plot_app.command("density")
def plot_density(
    input: Path = _input_argument(),
    output: Path = typer.Option(..., "--output", "-o", help="Figure to write (.png, .pdf or .svg)."),
    method: Optional[List[str]] = typer.Option(
        None, "--method", "-m", help="Run this method spec and draw it (repeatable); see compare."
    ),
    key: Optional[List[str]] = typer.Option(
        None, "--key", help="Draw this call stored in an .h5ad input (repeatable). Default: every stored call."
    ),
    min_counts: float = typer.Option(1, "--min-counts", help="Leave out barcodes below this total.", min=0),
    bins: int = typer.Option(100, "--bins", help="Histogram bins.", min=2),
    weighted: bool = typer.Option(False, "--weighted", help="Weight each barcode by its total."),
    title: Optional[str] = typer.Option(None, "--title", help="Figure title."),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    solo_feature: str = _solo_feature_option(),
) -> None:
    """Histogram of log10 totals with each call's cut-off.

    Examples:
        sc-droplets plot density raw_feature_bc_matrix.h5 -m otsu -m kde-valley -o density.png
        sc-droplets plot density calls.h5ad -o density.pdf --weighted
    """
    from sc_droplets.io.write import check_outputs

    try:
        check_outputs([output], [input])
        adata = _load(input, solo_feature, genome)
        types = _feature_types(feature_type)
        keys = _plot_keys(adata, input, key, method, layer, types, genome)
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from sc_droplets.pl import density

        ax = density(adata, keys, layer=layer, feature_types=types, min_counts=min_counts, bins=bins,
                     weighted=weighted, title=title)
        ax.figure.tight_layout()
        ax.figure.savefig(output, dpi=150)
        plt.close(ax.figure)
        console.print(f"Wrote {escape(str(output))}")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


@plot_app.command("splicing")
def plot_splicing(
    input: Path = _input_argument(),
    output: Path = typer.Option(..., "--output", "-o", help="Figure to write (.png, .pdf or .svg)."),
    method: Optional[List[str]] = typer.Option(
        None, "--method", "-m", help="Run this method spec and draw it (repeatable); see compare."
    ),
    key: Optional[List[str]] = typer.Option(
        None, "--key", help="Draw this call stored in an .h5ad input (repeatable). Default: every stored call."
    ),
    metric: str = typer.Option("unspliced_fraction", "--metric", help="The obs column to draw (mito_fraction...)."),
    min_counts: float = typer.Option(10, "--min-counts", help="Leave out barcodes below this total.", min=0),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    solo_feature: str = _solo_feature_option(),
    splicing: Optional[Path] = _splicing_option(),
) -> None:
    """The unspliced fraction against the calls: against count, in each call's cells and the
    barcodes left out, and where the calls disagree.

    Examples:
        sc-droplets plot splicing raw_feature_bc_matrix.h5 -m cellranger -m otsu --splicing velocyto.loom -o splicing.png
        sc-droplets plot splicing counts_unfiltered -m gmm -m kde-valley -o splicing.pdf
    """
    from sc_droplets.io.write import check_outputs

    try:
        check_outputs([output], [input])
        adata = _load(input, solo_feature, genome)
        _add_metrics(adata, splicing)
        types = _feature_types(feature_type)
        keys = _plot_keys(adata, input, key, method, layer, types, genome)
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from sc_droplets.pl import splicing as draw

        figure = draw(adata, keys, metric=metric, min_counts=min_counts, layer=layer, feature_types=types)
        figure.savefig(output, dpi=150)
        plt.close(figure)
        console.print(f"Wrote {escape(str(output))}")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


@plot_app.command("diagnostics")
def plot_diagnostics(
    input: Path = _input_argument(),
    output: Path = typer.Option(..., "--output", "-o", help="Figure to write (.png, .pdf or .svg)."),
    method: Optional[str] = typer.Option(None, "--method", "-m", help="Run this method spec and draw it."),
    key: Optional[str] = typer.Option(None, "--key", help="Draw this call stored in an .h5ad input."),
    feature_type: Optional[List[str]] = _feature_type_option(),
    layer: Optional[str] = _layer_option(),
    genome: Optional[str] = _genome_option(),
    solo_feature: str = _solo_feature_option(),
) -> None:
    """One call's own view: its cut and robust maximum, or its classes, likelihoods and p-values.

    Examples:
        sc-droplets plot diagnostics raw_feature_bc_matrix.h5 -m starsolo-emptydrops-cr -o diagnostics.png
        sc-droplets plot diagnostics calls.h5ad --key cellranger -o diagnostics.pdf
    """
    from sc_droplets.io.write import check_outputs

    try:
        if (method is None) == (key is None):
            raise ValueError("give exactly one of --method and --key")
        check_outputs([output], [input])
        adata = _load(input, solo_feature, genome)
        if key is not None:
            if key not in adata.uns.get("sc_droplets", {}):
                raise KeyError(f"no stored call {key} in {input}")
            target = key
        else:
            from sc_droplets.api import call

            target = call(adata, method, layer=layer, feature_types=_feature_types(feature_type), genome=genome).key
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from sc_droplets.pl import diagnostics

        figure = diagnostics(adata, target)
        figure.savefig(output, dpi=150)
        plt.close(figure)
        console.print(f"Wrote {escape(str(output))}")
    except typer.Exit:
        raise
    except Exception as e:
        raise _fail(e)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Run the CLI on `argv`, or on the command line when it is None."""
    app(args=list(argv) if argv is not None else None)
