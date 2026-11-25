"""numpy's classic `argsort(kind="quicksort")`, reproduced so ties fall where they fell.

Cell Ranger picks its ambient droplets, and the last few OrdMag cells, by
position in `np.argsort(counts)`. Among equal counts (and ranks 45,000 to
90,000 are nearly all ties) the order is whatever numpy's quicksort left.
Recent numpy releases sort with SIMD code on some processors (AVX2/AVX-512
on x86-64), which leaves ties in a different order; numpy before 1.25 (and
numpy on processors without those units) used the median-of-three quicksort
below, with insertion sort for partitions of 16 or fewer.

This reimplements that algorithm (numpy `npysort/quicksort`, BSD licence),
vectorising each partition step: the left scan stops at values >= pivot,
the right scan at values <= pivot, and the k-th stop of each swaps until the
scans cross. Insertion sort keeps equal values in order, so the small
partitions are finished by one stable sort. numpy 1.17 added a depth limit
that falls back to heapsort; Cell Ranger 3.x predates it, and ordinary count
vectors never reach it.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

SMALL_QUICKSORT = 15


def legacy_argsort(values: np.ndarray, window: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """`np.argsort(values, kind="quicksort")` as numpy < 1.25 computed it.

    With `window=(lo, hi)`, only positions lo..hi-1 of the result are
    guaranteed (the rest are left partly sorted): quicksort sorts each
    partition independently, so partitions outside the window are skipped.
    Callers that use a slice of the order (the top N, ranks 45,000–90,000)
    pass it and get the same slice, much faster.
    """
    v = np.asarray(values)
    n = v.shape[0]
    tosort = np.arange(n, dtype=np.int64)
    if n < 2:
        return tosort
    lo, hi = (0, n) if window is None else (max(0, window[0]), min(n, window[1]))
    if lo >= hi:
        return tosort

    def wanted(left: int, right: int) -> bool:
        return left <= right and right >= lo and left < hi

    small_segments = []
    stack = [(0, n - 1)]
    while stack:
        pl, pr = stack.pop()
        if not wanted(pl, pr):
            continue
        while pr - pl > SMALL_QUICKSORT:
            pm = pl + ((pr - pl) >> 1)
            if v[tosort[pm]] < v[tosort[pl]]:
                tosort[pm], tosort[pl] = tosort[pl], tosort[pm]
            if v[tosort[pr]] < v[tosort[pm]]:
                tosort[pr], tosort[pm] = tosort[pm], tosort[pr]
            if v[tosort[pm]] < v[tosort[pl]]:
                tosort[pm], tosort[pl] = tosort[pl], tosort[pm]
            vp = v[tosort[pm]]
            pj_start = pr - 1
            tosort[pm], tosort[pj_start] = tosort[pj_start], tosort[pm]
            segment = v[tosort[pl : pr + 1]]
            # Left scan from pl+1 stops at values >= vp (the pivot at pr-1 stops it at worst);
            # right scan from pr-2 down stops at values <= vp (v[pl] <= vp stops it).
            left = np.flatnonzero(segment[1 : pj_start - pl + 1] >= vp) + pl + 1
            right = (pj_start - 1) - np.flatnonzero(segment[pj_start - 1 - pl :: -1] <= vp)
            m = min(left.size, right.size)
            crossed = np.flatnonzero(left[:m] >= right[:m])
            k = int(crossed[0]) if crossed.size else m
            if k:
                a = left[:k]
                b = right[:k]
                tosort[a], tosort[b] = tosort[b].copy(), tosort[a].copy()
            # After k swaps the left scan stops at the next candidate, or earlier at
            # the last swapped right position (which now holds a value >= pivot).
            pi = int(left[k]) if k == 0 else int(min(left[k], right[k - 1]))
            tosort[pi], tosort[pj_start] = tosort[pj_start], tosort[pi]
            left_wanted = wanted(pl, pi - 1)
            right_wanted = wanted(pi + 1, pr)
            if left_wanted and right_wanted:
                if pi - pl < pr - pi:
                    stack.append((pi + 1, pr))
                    pr = pi - 1
                else:
                    stack.append((pl, pi - 1))
                    pl = pi + 1
            elif left_wanted:
                pr = pi - 1
            elif right_wanted:
                pl = pi + 1
            else:
                pl, pr = 0, -1
                break
        if pr > pl and wanted(pl, pr):
            small_segments.append((pl, pr))
    if small_segments:
        starts = np.array([s for s, _ in small_segments])
        ends = np.array([e for _, e in small_segments]) + 1
        lengths = ends - starts
        positions = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)])
        segment_id = np.repeat(np.arange(len(small_segments)), lengths)
        members = tosort[positions]
        # Insertion sort is stable: order by segment, then value, keeping current order.
        order = np.lexsort((np.arange(positions.size), v[members], segment_id))
        tosort[positions] = members[order]
    return tosort
