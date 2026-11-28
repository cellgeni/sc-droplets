"""R's default random numbers, and dqrng's seeds drawn from them.

DropletUtils seeds its Monte Carlo from R's random number generator, so with
the value given to `set.seed` its simulations are reproducible:

- `set.seed(s)` with the default kind (Mersenne-Twister, R >= 1.7) scrambles
  the seed with 50 rounds of the LCG s = 69069·s + 1 (mod 2**32), then fills
  the 624-word state with the next 624 values of the same LCG; the position
  is set to 624, so the first draw regenerates the state. From there it is
  Matsumoto and Nishimura's MT19937 with its usual tempering, which numpy's
  `MT19937` bit generator reproduces when handed that state.
- `unif_rand()` is the 32-bit output times 2.3283064365386963e-10 (2**-32),
  nudged into (0, 1) if it falls on 0.
- `R_unif_index(2**32)` with the default "Rejection" sample kind (R >= 3.6)
  builds a 32-bit integer from three `unif_rand()` draws, 16 bits each
  (`rbits(32)`), and keeps the low 32 bits: the second and third.
- dqrng's `generateSeedVectors(n, nwords=2)` makes each seed from two such
  integers, and `convert_seed<uint64_t>` joins them, the first as the high
  word.

Written from R's `src/main/RNG.c` (GPL-2 or later, Copyright (C) 1997-2024
The R Core Team) and dqrng's `inst/include/R_randgen.h` and `convert_seed.h`
(AGPL-3 or GPL-2 by their headers' choice; Copyright (C) 2018-2019 Ralf
Stubner and Aaron Lun).
"""

from __future__ import annotations

import numpy as np

_LCG = 69069
_MASK32 = 0xFFFFFFFF
_TO_UNIT = 2.3283064365386963e-10
_I2_32M1 = 2.328306437080797e-10  # 1 / (2**32 - 1), R's i2_32m1


class RMersenneTwister:
    """R's Mersenne-Twister after `set.seed(seed)`."""

    def __init__(self, seed: int):
        s = int(seed) & _MASK32
        for _ in range(50):
            s = (_LCG * s + 1) & _MASK32
        words = []
        for _ in range(625):
            s = (_LCG * s + 1) & _MASK32
            words.append(s)
        # words[0] is .Random.seed's position slot, which FixupSeeds resets to 624.
        self._bits = np.random.MT19937()
        self._bits.state = {
            "bit_generator": "MT19937",
            "state": {"key": np.array(words[1:], dtype=np.uint32), "pos": 624},
        }

    def raw(self, n: int) -> np.ndarray:
        """The next `n` tempered 32-bit outputs."""
        return self._bits.random_raw(int(n)).astype(np.uint64)

    def unif_rand(self, n: int) -> np.ndarray:
        """The next `n` values of R's `unif_rand()`."""
        u = self.raw(n).astype(np.float64) * _TO_UNIT
        u[u <= 0.0] = 0.5 * _I2_32M1
        u[1.0 - u <= 0.0] = 1.0 - 0.5 * _I2_32M1
        return u

    def unif_index_u32(self, n: int) -> np.ndarray:
        """`n` draws of `R_unif_index(2**32)` (rejection sampling), as uint64."""
        u = self.unif_rand(3 * int(n)).reshape(-1, 3)
        # rbits(32): v = (v1·65536 + v2)·65536 + v3, masked to 32 bits.
        high = np.floor(u[:, 1] * 65536.0).astype(np.uint64)
        low = np.floor(u[:, 2] * 65536.0).astype(np.uint64)
        return (high << np.uint64(16)) | low


def generate_seed_vectors(rng: RMersenneTwister, n: int) -> np.ndarray:
    """dqrng's `generateSeedVectors(n, nwords=2)`, joined into 64-bit seeds."""
    words = rng.unif_index_u32(2 * int(n)).reshape(-1, 2)
    return (words[:, 0] << np.uint64(32)) | words[:, 1]
