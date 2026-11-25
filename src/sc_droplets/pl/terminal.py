"""A knee plot in Unicode braille, for a terminal with no display.

Each character cell holds 2 × 4 dots, so a 60 × 14 character plot has a
120 × 56 dot canvas. Both axes are log10.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

# Dot bits of a braille cell, by (column, row) within the cell.
_BITS = ((0x01, 0x02, 0x04, 0x40), (0x08, 0x10, 0x20, 0x80))


def _format(value: float) -> str:
    if value >= 1e6:
        return f"{value / 1e6:.0f}M" if value >= 1e7 else f"{value / 1e6:.1f}M"
    if value >= 1e3:
        return f"{value / 1e3:.0f}k" if value >= 1e4 else f"{value / 1e3:.1f}k"
    return f"{value:.0f}"


def knee_text(
    totals: np.ndarray,
    *,
    width: int = 60,
    height: int = 14,
    marks: Optional[Dict[str, float]] = None,
) -> List[str]:
    """Lines of a braille knee plot of `totals` (any order, zeros ignored).

    `marks` maps a label to a total; each is drawn as a dotted horizontal line
    and named under the plot.
    """
    values = np.sort(np.asarray(totals, dtype=np.float64))[::-1]
    values = values[values > 0]
    if values.size == 0:
        return ["(no barcode has any counts)"]
    w, h = width * 2, height * 4
    top = np.log10(values[0])
    bottom = 0.0 if values[-1] >= 1 else np.log10(values[-1])
    span_y = max(top - bottom, 1e-9)
    span_x = max(np.log10(values.size), 1e-9)

    canvas = np.zeros((h, w), dtype=bool)

    def row_of(total: float) -> int:
        return int(round((top - np.log10(total)) / span_y * (h - 1)))

    previous = None
    for col in range(w):
        rank = 10 ** (col / (w - 1) * span_x)
        index = min(values.size - 1, int(round(rank)) - 1)
        row = row_of(values[index])
        rows = [row] if previous is None else range(min(previous, row), max(previous, row) + 1)
        for r in rows:
            canvas[r, col] = True
        previous = row
    legend = []
    for label, total in (marks or {}).items():
        if total is None or not np.isfinite(total) or total <= 0:
            continue
        r = min(h - 1, max(0, row_of(total)))
        canvas[r, ::3] = True
        legend.append(f"  ···· {label}: {_format(total)}")

    lines = []
    label_width = max(len(_format(values[0])), len(_format(values[-1]))) + 1
    for cell_row in range(height):
        chars = []
        for cell_col in range(width):
            code = 0x2800
            for dx in range(2):
                for dy in range(4):
                    if canvas[cell_row * 4 + dy, cell_col * 2 + dx]:
                        code |= _BITS[dx][dy]
            chars.append(chr(code))
        if cell_row == 0:
            label = _format(values[0])
        elif cell_row == height - 1:
            label = _format(10**bottom)
        else:
            label = ""
        lines.append(f"{label:>{label_width}} ┤" + "".join(chars))
    axis = " " * label_width + " └" + "─" * width
    ticks = " " * (label_width + 2) + "1" + " " * (width - 1 - len(_format(values.size))) + _format(values.size)
    lines.extend([axis, ticks, " " * (label_width + 2) + "barcode rank (log) vs counts (log)"])
    lines.extend(legend)
    return lines
