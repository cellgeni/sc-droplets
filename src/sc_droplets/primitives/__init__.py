"""Building blocks the methods share, each written to reproduce one tool's arithmetic.

- `rounding`: Python 2 / C `round` (half away from zero) beside numpy's
  half-to-even.
- `goodturing`: Simple Good-Turing as Cell Ranger and as STAR compute it.
- `multinomial`: multinomial log-likelihoods of sparse count rows.
- `pvalues`: Monte Carlo p-values and the Benjamini-Hochberg variants.
- `stdrandom`: C++ `std::mt19937` and libstdc++'s `discrete_distribution`,
  for STARsolo's simulations.
- `simulate`: the Monte Carlo engines, one per tool and release.
"""
