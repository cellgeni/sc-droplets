"""zUMIs' automatic cell detection (2.5.6b and later): `inflection::uik` on the cumulative reads.

Page: docs/methods/zumis/cell-barcodes.md. Ported from zUMIs' `barcodeIDFUN.R` (`cellBC`,
`.cellBarcode_unknown`, `.FindBCcut`; GPL-3) and from the R package inflection 1.3.5 - 1.3.7
(`uik`, `check_curve`, `findipl`, `ede`, `lin2`; GPL >= 2), as sc-droplets is GPL-3. R's `sum` and
`cumsum` accumulate in long double (`rmath.ld_sum`, `ld_cumsum`); `quantile` is type 7.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.primitives.rmath import ld_cumsum, ld_sum
from sc_droplets.primitives.rstats import quantile7
from sc_droplets.registry import Gap, Method, MethodInfo, Reference, Variant, param, register

N_READS_PER_CELL = 100
FEWEST = 10
FALLBACK = 100


class ZumisError(ValueError):
    """Where zUMIs (or inflection) stops with an error."""


def lin2(x1, y1, x2, y2, x):
    """inflection's `lin2`: the line through (x1, y1) and (x2, y2) at x."""
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def findipl(x: np.ndarray, y: np.ndarray, j: int):
    """inflection's `findipl` at the 1-based point j: the trapezoid integrals of the curve minus its
    chord, left (sl) and right (sr) of j. At the ends R's indexing leaves an empty product: 0."""
    n = x.size
    if j >= 2:
        fl = y[:j] - lin2(x[0], y[0], x[j - 1], y[j - 1], x[:j])
        sl = ld_sum(np.diff(x[:j]) * (0.5 * (fl[:-1] + fl[1:])))
    else:
        sl = 0.0
    if j <= n - 1:
        fr = y[j - 1:] - lin2(x[j - 1], y[j - 1], x[-1], y[-1], x[j - 1:])
        sr = ld_sum(np.diff(x[j - 1:]) * (0.5 * (fr[:-1] + fr[1:])))
    else:
        sr = 0.0
    return sl, sr


def check_curve(x: np.ndarray, y: np.ndarray) -> int:
    """inflection's `check_curve`: 1 if the curve is concave (or concave then convex), else 0."""
    n = x.size
    points = [1] + [int(quantile7(np.arange(1, n + 1, dtype=np.float64), p)) for p in (0.25, 0.5, 0.75)] + [n]
    integrals = [findipl(x, y, j) for j in points]
    left = np.sign([integrals[k][0] for k in (1, 2, 3, 4)])
    right = np.sign([integrals[k][1] for k in (0, 1, 2, 3)])
    if np.unique(left).size == 1:
        cleft = "concave" if left[0] > 0 else "convex"
    else:
        cleft = "concave" if left[0] > 0 else "convex"
    if np.unique(right).size == 1:
        cright = "concave" if right[0] > 0 else "convex"
    else:
        cright = "concave" if right[3] > 0 else "convex"
    return 1 if cleft == "concave" else 0  # concave_convex and concave: 1; convex_concave and convex: 0


def uik(x: np.ndarray, y: np.ndarray) -> float:
    """inflection's `uik`: x at the first minimum of y minus its chord (y negated if concave)."""
    if x.size <= 3:
        raise ZumisError("Method is not applicable for such a small vector (inflection::uik)")
    if check_curve(x, y) == 1:
        y = -y
    lf = y - lin2(x[0], y[0], x[-1], y[-1], x)
    return float(x[int(np.argmin(lf))])


@register
class CellBarcodes(Method):
    info = MethodInfo(
        id="zumis-cell-barcodes",
        title="zUMIs cell barcodes",
        tool="zUMIs",
        category="knee",
        summary="zUMIs' automatic cell detection: the unit invariant knee (inflection::uik) of the cumulative "
                "reads; every barcode with at least the knee's count.",
        docs="methods/zumis/cell-barcodes.md",
        counts="either",
        variants=(Variant("v2.5.6b", since="2.5.6b", summary="inflection::uik on the cumulative reads"),),
        gaps=(
            Gap("0.0.1", "2.5.6b",
                "zUMIs 0.0.1 had no automatic detection, and 0.0.2 - 2.5.4 fitted mclust::Mclust, which from "
                "mclust 5.4 initialises on a random subset of 2,000 barcodes with no seed: not reproducible"),
        ),
        latest_known="2.9.7",
        references=(
            Reference(
                "Parekh, Ziegenhain, Vieth, Enard & Hellmann 2018, zUMIs - A fast and flexible pipeline to "
                "process RNA sequencing data with UMIs, GigaScience 7(6)",
                doi="10.1093/gigascience/giy059",
            ),
            Reference(
                "Christopoulos 2016, Introducing Unit Invariant Knee (UIK) As an Objective Choice for Elbow Point "
                "in Multivariate Data Analysis Techniques, SSRN",
                doi="10.2139/ssrn.3043076",
            ),
        ),
        source="https://github.com/sdparekh/zUMIs/blob/2.9.7/barcodeIDFUN.R",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        n_reads_per_cell: int = param(N_READS_PER_CELL, "Barcodes with fewer reads are left out",
                                      original="barcodes: nReadsperCell", minimum=0)
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="reads per barcode", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        kept = np.flatnonzero(values >= params.n_reads_per_cell)
        order = kept[np.argsort(-values[kept], kind="stable")]  # data.table's order(-n) is stable
        n = values[order]
        if n.size == 0:
            raise ZumisError(f"no barcode has --n-reads-per-cell {params.n_reads_per_cell} reads")
        cs = ld_cumsum(n.astype(np.float64))
        index = np.arange(1, n.size + 1, dtype=np.float64)
        knee = int(uik(index, cs / 1000))
        cut = n[knee - 1]
        is_cell = np.zeros(values.size, dtype=bool)
        notes = []
        if np.sum(n >= cut) < FEWEST:
            is_cell[order[:FALLBACK]] = True
            notes.append(f"fewer than {FEWEST} cells: zUMIs takes the top {min(FALLBACK, n.size)} barcodes instead")
        else:
            is_cell[order[n >= cut]] = True
        return result(type(self), variant, params, droplets, is_cell,
                      thresholds={params.unit: float(cut), "knee_rank": float(knee)}, notes=notes)
