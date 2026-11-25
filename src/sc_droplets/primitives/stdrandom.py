"""C++ standard-library random numbers, as STAR's Linux builds draw them.

STARsolo's EmptyDrops_CR draws from `std::discrete_distribution` with a
`std::mt19937` engine. Reproducing its simulations needs both bit for bit:

- `std::mt19937(seed)` is Matsumoto and Nishimura's MT19937 seeded with
  `init_genrand(seed mod 2**32)`, which is what numpy's legacy
  `RandomState(seed)` does; its raw 32-bit outputs are
  `randint(0, 2**32, dtype=uint32)` (the full range bypasses rejection).
- libstdc++'s `discrete_distribution` normalises the weights by their
  sequential sum, takes their sequential partial sums, sets the last to 1,
  and returns the first index whose cumulative probability is >= u, with
  u from `generate_canonical<double, 53>`: two engine outputs a then b,
  u = (a + b·2**32) / 2**64, and 1 - 2**-53 if that rounds to 1.

libc++ (clang, macOS) implements `discrete_distribution` differently, so this
matches STAR binaries built with GCC, which the released Linux binaries are.
"""

from __future__ import annotations

import numpy as np

_TWO_32 = 4294967296.0
_TWO_64 = 18446744073709551616.0
_JUST_BELOW_ONE = np.nextafter(1.0, 0.0)


def mt19937_raw(seed: int, n: int) -> np.ndarray:
    """The first `n` 32-bit outputs of `std::mt19937(seed)`."""
    state = np.random.RandomState(int(seed) % (1 << 32))
    return state.randint(0, 1 << 32, size=n, dtype=np.uint32)


def generate_canonical(raw: np.ndarray) -> np.ndarray:
    """libstdc++'s `generate_canonical<double, 53>` from consecutive engine outputs.

    `raw` holds pairs (a, b), the first of each pair the low word.
    """
    raw = np.asarray(raw, dtype=np.float64)
    low = raw[0::2]
    high = raw[1::2]
    # sum = 0.0 + a * 1.0, then sum += b * 2**32: the same two IEEE operations.
    u = (low + high * _TWO_32) / _TWO_64
    u[u >= 1.0] = _JUST_BELOW_ONE
    return u


class LibstdcxxDiscrete:
    """`std::discrete_distribution<>` as libstdc++ builds and samples it."""

    def __init__(self, weights: np.ndarray):
        weights = np.asarray(weights, dtype=np.float64)
        if weights.size < 2:
            # libstdc++ clears a one-weight distribution, which then always returns 0.
            self.cumulative = np.ones(1)
            return
        total = np.cumsum(weights)[-1]  # std::accumulate: left to right
        probabilities = weights / total
        cumulative = np.cumsum(probabilities)  # std::partial_sum
        cumulative[-1] = 1.0
        self.cumulative = cumulative

    def sample_from_raw(self, raw: np.ndarray) -> np.ndarray:
        """Indices drawn from engine outputs, two per draw (lower_bound)."""
        return np.searchsorted(self.cumulative, generate_canonical(raw), side="left")
