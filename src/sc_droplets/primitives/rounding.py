"""Rounding as each original tool does it.

Python 2's `round`, C's `round` and R's `round` on .5 differ from Python 3's
and numpy's, which round half to even. Cell Ranger up to 3.x ran on Python 2;
7.2 onwards calls `np.round`. STAR is C++.
"""

from __future__ import annotations

import math

import numpy as np


def round_half_away(x: float) -> float:
    """Python 2 `round(x)` and C `round(x)`: halves go away from zero.

    Exact, as Python 2.7's correctly rounded `round`: `0.49999999999999994`
    rounds to 0, where `floor(x + 0.5)` would give 1.
    """
    if not math.isfinite(x):
        return x
    magnitude = abs(x)
    whole = math.floor(magnitude)
    # magnitude - whole is exact for any double below 2**52.
    y = whole + 1.0 if magnitude - whole >= 0.5 else whole
    return math.copysign(y, x)


def round_half_even(x):
    """numpy's `np.round` (and Python 3's `round`): halves go to the even neighbour."""
    return np.round(x)
