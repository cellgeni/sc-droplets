"""Plots of calls: `scd.pl.knee` (more arrive with each milestone).

Every plot takes an AnnData, draws the calls stored on it (or runs the method
specs it is given), returns the matplotlib Axes, and accepts `ax=` to draw
into an existing figure.
"""

from sc_droplets.pl.knee import knee
from sc_droplets.pl.terminal import knee_text

__all__ = ["knee", "knee_text"]
