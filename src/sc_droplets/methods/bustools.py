"""bustools `allowlist` (`whitelist` before 0.43.1): a fraction of the mean of the top ten counts
in a sample from the start of the BUS file.

Page: docs/methods/bustools/allowlist.md. Ported from bustools' `src/bustools_whitelist.cpp`
(0.45.1; the same rule since 0.39.3), under the BSD 2-Clause licence:

    Copyright (c) 2018, BUStools. All rights reserved.
    Redistribution and use in source and binary forms, with or without modification, are
    permitted provided that the following conditions are met: 1. Redistributions of source code
    must retain the above copyright notice, this list of conditions and the following
    disclaimer. 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation and/or other
    materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
    AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.

bustools streams the sorted BUS file; the threshold comes from the barcodes of its first blocks
of 100,000 records. A matrix has no records, so each UMI of a barcode stands for one (its reads,
with --unit reads, spread evenly over them), in bustools' barcode order.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "bustools"
LATEST = "0.45.1"
BLOCK = 100_000  # records read at a time
SAMPLE = 200  # barcodes wanted in the sample
TOP = 10
ERROR_RATE = 0.01


class BustoolsError(ValueError):
    """Where bustools cannot run (or would read past its sample)."""


def sequences(barcodes) -> Tuple[np.ndarray, int]:
    """The barcodes as bustools codes them: the sequence before any `-` suffix, A/C/G/T only."""
    seqs = np.array([str(b).split("-", 1)[0] for b in barcodes], dtype=object)
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise BustoolsError(f"bustools needs barcodes of one length; found lengths {sorted(lengths)}")
    if not set("".join(seqs)) <= set("ACGT"):
        raise BustoolsError("bustools codes barcodes in two bits: only A, C, G and T")
    if len(set(seqs)) != len(seqs):
        raise BustoolsError("barcodes repeat once their suffix is dropped (several libraries?): bustools would merge them")
    return seqs, lengths.pop()


@dataclass
class Sample:
    """What bustools read before deciding: the sample (barcode position, count) and its blocks."""

    entries: List[Tuple[int, int]]
    blocks: int


def sample(records: np.ndarray, counts: np.ndarray, first_is_zero: bool) -> Sample:
    """The first phase of `bustools_whitelist`: whole blocks of records until the sample has 200
    entries. `records` and `counts` are per barcode, in file order; a barcode's count is spread
    evenly over its records. `first_is_zero`: the first barcode is all A's (code 0), which the
    initial `curr_bc = 0` mistakes for a barcode already open (bustools 0.42.0+)."""
    ends = np.cumsum(records)
    total = int(ends[-1]) if ends.size else 0
    starts = ends - records

    def part(j: int, a: int, b: int) -> int:
        """The count of barcode j's records a..b-1 (local positions)."""
        q, r = divmod(int(counts[j]), int(records[j]))
        return q * (b - a) + min(b, r) - min(a, r)

    entries: List[Tuple[int, int]] = []
    current, count = (0, -1) if first_is_zero else (None, -1)
    blocks = 0
    while True:
        lo, hi = blocks * BLOCK, min((blocks + 1) * BLOCK, total)
        if lo >= total:
            break
        blocks += 1
        first = bisect.bisect_right(ends, lo)
        last = bisect.bisect_right(ends, hi - 1)
        for j in range(first, last + 1):
            a, b = max(int(starts[j]), lo) - int(starts[j]), min(int(ends[j]), hi) - int(starts[j])
            if j == current:
                count += part(j, a, b)
            else:
                if count != -1:
                    entries.append((current, count))
                current, count = j, part(j, a, b)
        if count != -1:
            entries.append((current, count))  # the barcode open at the block's end, so far
        if len(entries) >= SAMPLE:
            break
    return Sample(entries, blocks)


def threshold_of(entries: List[Tuple[int, int]], barcode_length: int) -> Tuple[int, int]:
    """(threshold, mean of the top ten): the sample by count, decreasing, ties by barcode."""
    if len(entries) < TOP:
        raise BustoolsError(f"the sample has {len(entries)} barcodes; bustools reads past its end below {TOP}")
    top = sorted(entries, key=lambda e: (-e[1], e[0]))[:TOP]
    total = sum(c for _, c in top)
    if total > 2**31 - 1:
        raise BustoolsError("the top ten counts overflow bustools' 32-bit sum")
    mean = total // TOP
    return int(mean * (1 - math.pow(1 - ERROR_RATE, barcode_length))), mean


@register
class Allowlist(Method):
    info = MethodInfo(
        id="bustools-allowlist",
        title="bustools allowlist",
        tool=TOOL,
        category="rank-threshold",
        summary=(
            "bustools' allowlist (kb --filter bustools): every barcode with at least 1 - 0.99^L of the "
            "mean of the top ten counts in a sample from the start of the BUS file."
        ),
        docs="methods/bustools/allowlist.md",
        counts="either",
        variants=(Variant("default", since="0.39.3", summary="one rule since whitelist was added"),),
        latest_known=LATEST,
        references=(
            Reference(
                "Melsted, Ntranos & Pachter 2019, The barcode, UMI, set format and BUStools, Bioinformatics 35:4472-4473",
                doi="10.1093/bioinformatics/btz279",
            ),
        ),
        source="https://github.com/BUStools/bustools/blob/v0.45.1/src/bustools_whitelist.cpp",
        original_license="BSD-2-Clause",
    )

    @dataclass
    class Params:
        threshold: Optional[int] = param(None, "The threshold itself (bustools' sample is then not read)",
                                         original="--threshold", minimum=1)
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="reads per barcode",
                          choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        umis = droplets.totals
        present = np.flatnonzero(values > 0)
        if not np.all(values[present] == np.floor(values[present])):
            raise BustoolsError("bustools counts whole reads; these counts have fractions")
        seqs, length = sequences(np.asarray(droplets.barcodes)[present])
        order = np.argsort(seqs.astype(str), kind="stable")  # A < C < G < T: the 2-bit code's order
        file_index = present[order]
        counts = values[file_index].astype(np.int64)
        records = np.maximum(umis[file_index], 1).astype(np.int64)
        diagnostics = {"unit": params.unit, "barcode_length": length}
        thresholds = {}
        totals = counts.copy()
        if params.threshold is not None:
            cut = int(params.threshold)
        else:
            first_is_zero = bool(seqs.size) and seqs[order[0]] == "A" * length
            found = sample(records, counts, first_is_zero)
            cut, mean = threshold_of(found.entries, length)
            if first_is_zero:
                totals[0] -= 1  # its first record went to the "already open" branch
            thresholds["top10_mean"] = float(mean)
            diagnostics.update({"sample_size": len(found.entries), "sample_blocks": found.blocks})
        thresholds[params.unit] = float(cut)
        is_cell = np.zeros(values.shape[0], dtype=bool)
        is_cell[file_index[totals >= cut]] = True
        return result(type(self), variant, params, droplets, is_cell, thresholds=thresholds, diagnostics=diagnostics)
