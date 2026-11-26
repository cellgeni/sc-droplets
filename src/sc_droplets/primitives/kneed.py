"""Kneedle as `kneed.KneeLocator` computes it, in each behaviour its releases have had.

Page: docs/methods/kneedle/kneedle.md. Behaviours (kneed releases):

- "0.5.1" (0.5.1 – 0.6.0): extrema with >=/<=; polynomial by scikit-learn; a convex curve is
  mirrored in x as well as y; the walk stops at x == 1.0.
- "0.7" (0.7.0 – 0.8.5): polynomial by numpy.polyfit with a degree; only y is transformed.
  (0.8.3 stops the walk at the last index instead of at x == 1.0; with distinct x only the last
  point is at 1.0, so the two agree.)
- "0.8.6": after a local minimum, no knee is detected until the next local maximum.
- "0.9": no interpolation step; a concave decreasing or convex increasing curve is mirrored in x.

Ported from kneed (https://github.com/arvkevi/kneed), under the BSD 3-Clause licence:

    Copyright (c) 2017, Kevin Arvai
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are
    permitted provided that the following conditions are met:
    1. Redistributions of source code must retain the above copyright notice, this list of
       conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice, this list of
       conditions and the following disclaimer in the documentation and/or other materials
       provided with the distribution.
    3. Neither the name of the copyright holder nor the names of its contributors may be used
       to endorse or promote products derived from this software without specific prior
       written permission.
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
    TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

BEHAVIOURS = ("0.5.1", "0.7", "0.8.6", "0.9")
CURVES = ("convex", "concave")
DIRECTIONS = ("increasing", "decreasing")


@dataclass
class Kneedle:
    """Kneedle's curves for one (x, y); `find` walks them for a sensitivity."""

    x: np.ndarray
    x_normalized: np.ndarray
    y_difference: np.ndarray
    maxima: np.ndarray
    minima: np.ndarray
    curve: str
    direction: str
    behaviour: str
    knees: List[int] = field(default_factory=list)

    def thresholds(self, S: float) -> np.ndarray:
        """T_k = d at the k-th maximum - S |mean step of x_n|, as kneed computes it."""
        return self.y_difference[self.maxima] - (S * np.abs(np.diff(self.x_normalized).mean()))

    def find(self, S: float, online: bool = False) -> Optional[int]:
        """kneed's `find_knee`: the knee as an index into `x`, or None.

        `knees` is left holding every knee the walk found, in order.
        """
        self.knees = []
        if not self.maxima.size:
            return None
        Tmx = self.thresholds(S)
        n = self.x_normalized.shape[0]
        is_max = np.zeros(n, dtype=bool)
        is_max[self.maxima] = True
        is_min = np.zeros(n, dtype=bool)
        is_min[self.minima] = True
        gate = self.behaviour in ("0.8.6", "0.9")
        maxima_seen = 0
        threshold = threshold_index = None
        detection_active = True
        knee = None
        for i in range(int(self.maxima[0]), n):
            if self.behaviour == "0.5.1":
                if self.x_normalized[i] == 1.0:
                    break
            elif i == n - 1:
                break
            j = i + 1
            if is_max[i]:
                threshold = Tmx[maxima_seen]
                threshold_index = i
                maxima_seen += 1
                detection_active = True
            if is_min[i]:
                threshold = 0.0
                if gate:
                    detection_active = False
            if (detection_active or not gate) and self.y_difference[j] < threshold:
                knee = self._knee_index(threshold_index)
                if knee not in self.knees:
                    self.knees.append(knee)
                if not online:
                    return knee
        return knee if self.knees else None

    def _knee_index(self, threshold_index: int) -> int:
        """kneed's knee for the maximum at `threshold_index`: that point, or for a convex
        increasing or concave decreasing curve, the one at the mirrored index."""
        n = self.x.shape[0]
        mirrored = (self.curve == "convex") != (self.direction == "decreasing")
        return n - 1 - threshold_index if mirrored else threshold_index

    def maximum_of(self, knee_index: int) -> int:
        """The maximum of the difference curve a knee came from (the mapping is its own inverse)."""
        return self._knee_index(knee_index)


def _normalize(a: np.ndarray, as_float: bool) -> np.ndarray:
    if as_float:
        a = np.asarray(a, dtype=float)
    return (a - min(a)) / (max(a) - min(a))


def prepare(
    x,
    y,
    *,
    curve: str = "concave",
    direction: str = "increasing",
    interp_method: str = "interp1d",
    polynomial_degree: int = 7,
    behaviour: str = "0.9",
) -> Kneedle:
    """Steps 1–5 of KneeLocator (smoothing, scaling, transform, difference curve, extrema)."""
    from scipy.signal import argrelextrema

    if behaviour not in BEHAVIOURS:
        raise ValueError(f"unknown kneed behaviour {behaviour!r}; one of {BEHAVIOURS}")
    if curve not in CURVES or direction not in DIRECTIONS:
        raise ValueError("Please check that the curve and direction arguments are valid.")
    x = np.array(x)
    y = np.array(y)
    if x.shape != y.shape:
        raise ValueError("x and y must be the same length.")
    if interp_method == "interp1d":
        if behaviour == "0.9":
            smooth = y.astype(float)
        else:
            from scipy import interpolate

            smooth = interpolate.interp1d(x, y)(x)
    elif interp_method == "polynomial":
        if behaviour == "0.5.1":
            if polynomial_degree != 7:
                raise ValueError("kneed before 0.7.0 fits a polynomial of degree 7 only")
            from sklearn.linear_model import LinearRegression
            from sklearn.preprocessing import PolynomialFeatures

            features = PolynomialFeatures(7)
            model = LinearRegression()
            model.fit(features.fit_transform(x.reshape(-1, 1)), y)
            smooth = model.predict(features.fit_transform(x.reshape(-1, 1)))
        else:
            smooth = np.poly1d(np.polyfit(x, y, polynomial_degree))(x)
    else:
        raise ValueError(f"{interp_method} is an invalid interp_method parameter, use either 'interp1d' or 'polynomial'")

    as_float = behaviour == "0.9"
    x_n = _normalize(x, as_float)
    y_n = _normalize(smooth, as_float)
    if behaviour == "0.5.1":
        if curve == "convex":
            x_n = x_n.max() - x_n
            y_n = y_n.max() - y_n
        if direction == "decreasing":
            y_n = np.flip(y_n, axis=0)
        if curve == "convex":
            x_n = np.flip(x_n, axis=0)
            y_n = np.flip(y_n, axis=0)
    else:
        if direction == "decreasing":
            if curve == "concave":
                y_n = np.flip(y_n)
            else:
                y_n = y_n.max() - y_n
        elif curve == "convex":
            y_n = np.flip(y_n.max() - y_n)
        if behaviour == "0.9" and (direction, curve) in (("decreasing", "concave"), ("increasing", "convex")):
            x_n = np.flip(1 - x_n)
    difference = y_n - x_n
    return Kneedle(
        x=x,
        x_normalized=x_n,
        y_difference=difference,
        maxima=argrelextrema(difference, np.greater_equal)[0],
        minima=argrelextrema(difference, np.less_equal)[0],
        curve=curve,
        direction=direction,
        behaviour=behaviour,
    )


def knee_locator(x, y, S: float = 1.0, curve: str = "concave", direction: str = "increasing",
                 interp_method: str = "interp1d", online: bool = False, polynomial_degree: int = 7,
                 behaviour: str = "0.9") -> Optional[float]:
    """`KneeLocator(...).knee`: the x value of the knee, or None."""
    found = prepare(x, y, curve=curve, direction=direction, interp_method=interp_method,
                    polynomial_degree=polynomial_degree, behaviour=behaviour)
    index = found.find(S, online)
    return None if index is None else found.x[index]
