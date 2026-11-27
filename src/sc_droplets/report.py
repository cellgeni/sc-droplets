"""One self-contained HTML page for a sample and its calls: `scd.report`.

The figures are the package's own plots, inlined as SVG; scatter layers with many points are
rasterised inside the SVG so the page stays small. Nothing is loaded from elsewhere: the page can
be mailed, archived or opened offline.
"""

from __future__ import annotations

import html
import io
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np

from sc_droplets.result import UNS_KEY

#: Scatter layers with more points than this are drawn as an image inside the SVG.
RASTER_POINTS = 2000

_CSS = """
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; margin: 2em auto;
       max-width: 1200px; padding: 0 1em; color: #222; }
h1 { font-size: 1.6em; margin-bottom: 0.2em; }
h2 { font-size: 1.25em; margin-top: 2em; border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }
h3 { font-size: 1.05em; margin-top: 1.5em; }
.meta { color: #666; font-size: 0.9em; }
table { border-collapse: collapse; margin: 0.8em 0; font-size: 0.88em; }
th, td { border: 1px solid #ddd; padding: 0.3em 0.6em; text-align: left; vertical-align: top; }
th { background: #f5f5f5; }
td.number { text-align: right; font-variant-numeric: tabular-nums; }
figure { margin: 1em 0; }
figure svg { max-width: 100%; height: auto; }
figcaption { color: #555; font-size: 0.88em; margin-top: 0.3em; }
nav a { margin-right: 1em; }
.note { color: #8a5a00; }
"""


def _svg(figure) -> str:
    """A matplotlib figure as inline SVG, heavy scatter layers rasterised, the figure closed."""
    import matplotlib.pyplot as plt

    for ax in figure.axes:
        for collection in ax.collections:
            offsets = collection.get_offsets()
            if offsets is not None and len(offsets) > RASTER_POINTS:
                collection.set_rasterized(True)
    buffer = io.StringIO()
    figure.savefig(buffer, format="svg", bbox_inches="tight", dpi=110)
    plt.close(figure)
    text = buffer.getvalue()
    return text[text.index("<svg"):]


def _cell(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return '<td class="number">—</td>'
    if isinstance(value, (bool, np.bool_)):
        return f"<td>{'yes' if value else 'no'}</td>"
    if isinstance(value, (int, np.integer)):
        return f'<td class="number">{int(value):,}</td>'
    if isinstance(value, (float, np.floating)):
        return f'<td class="number">{float(value):.4g}</td>'
    if isinstance(value, (list, tuple, np.ndarray)):
        return f"<td>{html.escape(', '.join(str(v) for v in np.ravel(value)))}</td>"
    return f"<td>{html.escape(str(value))}</td>"


def _table(header: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in header)
    body = "".join("<tr>" + "".join(_cell(v) for v in row) + "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _figure(anchor: str, svg: str, caption: str) -> str:
    return f'<figure id="{anchor}">{svg}<figcaption>{html.escape(caption)}</figcaption></figure>'


def _summary_rows(adata: Any, keys: Sequence[str]) -> List[List[Any]]:
    rows = []
    for key in keys:
        record = adata.uns[UNS_KEY][key]
        rows.append([key, record.get("method", ""), record.get("tool_version", "") or "", record.get("variant", ""),
                     int(record.get("n_cells", 0)), record.get("min_cell_total")])
    return rows


def report(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    path: Union[None, str, Path] = None,
    *,
    title: Optional[str] = None,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
) -> str:
    """Build a self-contained HTML report of a sample's calls; write it to `path` when given, and
    return it.

    `calls` takes stored keys or method specs (a spec not yet stored is run first); None reports
    every stored call. The page holds: a table of the calls (with the separation scores when
    `obs` has the unspliced or mitochondrial fraction, `scd.add_metrics`); the knee, density,
    comparison and quality plots, and the splicing plot when the unspliced fraction is there;
    and for each call its parameters, thresholds, notes and diagnostic view.
    """
    from sc_droplets import __version__, pl
    from sc_droplets.metrics import MITO, UNSPLICED, separation
    from sc_droplets.pl._data import resolve_calls
    from sc_droplets.registry import get as get_method

    keys = resolve_calls(adata, calls)
    if not keys:
        raise ValueError("nothing to report: give method specs or stored keys")
    source: Dict[str, Any] = dict(adata.uns.get("sc_droplets_input", {}))
    name = title or (Path(str(source["path"])).name if "path" in source else "sc-droplets report")
    drawing = {"layer": layer, "feature_types": feature_types}

    sections: List[str] = []
    nav: List[str] = []

    def section(anchor: str, heading: str, body: str) -> None:
        nav.append(f'<a href="#{anchor}">{html.escape(heading)}</a>')
        sections.append(f'<section id="{anchor}"><h2>{html.escape(heading)}</h2>{body}</section>')

    header = ["key", "method", "tool version", "variant", "cells", "smallest cell total"]
    rows = _summary_rows(adata, keys)
    for metric in (UNSPLICED, MITO):
        if metric in adata.obs.columns:
            scores = separation(adata, keys, metric=metric)
            short = "unspliced" if metric == UNSPLICED else "mito"
            header += [f"{short}: cells", f"{short}: others", f"{short}: AUROC"]
            for row, key in zip(rows, keys):
                row += [scores.loc[key, f"{metric}_cells"], scores.loc[key, f"{metric}_others"],
                        scores.loc[key, f"{metric}_auroc"]]
    note = ("<p>Separation scores: the median in each call's cells and in the barcodes it left out with "
            "at least 10 counts, and the AUROC between them (0.5: no separation). Indicators, not ground "
            "truth.</p>") if len(header) > 6 else ""
    section("summary", "Calls", _table(header, rows) + note)

    plots = [("knee", "Knee", pl.knee(adata, keys, **drawing).figure,
              "The barcode-rank curve with each call's smallest cell total."),
             ("density", "Density", pl.density(adata, keys, **drawing).figure,
              "The histogram of log10 totals with each call's cut.")]
    if len(keys) > 1:
        plots.append(("compare", "Agreement", pl.compare(adata, keys, **drawing),
                      "Where the calls part along the rank curve, how many cells each called, and the Jaccard "
                      "index of their cell sets."))
    plots.append(("qc", "Quality", pl.qc(adata, keys, **drawing),
                  "Genes detected (and the mitochondrial fraction) against count, coloured by how many calls "
                  "took each barcode."))
    if UNSPLICED in adata.obs.columns:
        plots.append(("splicing", "Splicing", pl.splicing(adata, keys, **drawing),
                      "The unspliced fraction against count, in each call's cells and the barcodes it left out, "
                      "and where the calls disagree."))
    for anchor, heading, figure, caption in plots:
        section(anchor, heading, _figure(f"{anchor}-figure", _svg(figure), caption))

    parts = []
    for key in keys:
        record = adata.uns[UNS_KEY][key]
        method = str(record.get("method", ""))
        try:
            docs = get_method(method).info.docs_url
        except Exception:
            docs = None
        link = f' (<a href="{html.escape(docs)}">method page</a>)' if docs else ""
        parts.append(f'<h3 id="call-{html.escape(key)}">{html.escape(key)}: {html.escape(str(record.get("title", method)))}{link}</h3>')
        params = dict(record.get("params", {}))
        thresholds = dict(record.get("thresholds", {}))
        parts.append(_table(["parameter", "value"], sorted(params.items())))
        if thresholds:
            parts.append(_table(["threshold", "value"], sorted(thresholds.items())))
        notes = [str(n) for n in np.atleast_1d(record.get("notes", [])) if str(n)]
        if notes:
            parts.append("<ul>" + "".join(f'<li class="note">{html.escape(n)}</li>' for n in notes) + "</ul>")
        parts.append(_figure(f"call-{key}-figure", _svg(pl.diagnostics(adata, key)), f"{key}: the method's own view."))
    section("calls", "Each call", "".join(parts))

    shape = f"{adata.n_obs:,} barcodes × {adata.n_vars:,} features"
    # The file's name, not its path: a report is for sharing, and a path tells of the machine.
    origin = f" from {html.escape(Path(str(source['path'])).name)} ({html.escape(str(source.get('format', '')))})" if "path" in source else ""
    when = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    page = (
        "<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        f"<title>{html.escape(name)}</title><style>{_CSS}</style></head><body>"
        f"<h1>{html.escape(name)}</h1>"
        f'<p class="meta">{shape}{origin}. sc-droplets {html.escape(__version__)}, {when}.</p>'
        f"<nav>{''.join(nav)}</nav>{''.join(sections)}</body></html>\n"
    )
    if path is not None:
        from sc_droplets.io.write import write_text

        write_text(page, Path(path))
    return page
