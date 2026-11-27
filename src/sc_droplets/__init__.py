"""sc-droplets: every cell vs empty-droplet calling method, version by version.

    import sc_droplets as scd
    adata = scd.read("raw_feature_bc_matrix.h5")
    scd.call(adata, "top-n", n_cells=3000)
    scd.pl.knee(adata)

Documentation: https://cellgeni.github.io/sc-droplets/
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("sc-droplets")
except PackageNotFoundError:  # pragma: no cover - running from a source tree without install
    __version__ = "0.0.0+unknown"

from sc_droplets import pl
from sc_droplets.api import agreement, call, compare, list_methods, stored
from sc_droplets.io import read
from sc_droplets.metrics import add_metrics, separation, unspliced_fraction
from sc_droplets.registry import UnknownMethod, UnsupportedVersion
from sc_droplets.result import CallResult
from sc_droplets.simulate import simulate_droplets

__all__ = [
    "__version__",
    "CallResult",
    "UnknownMethod",
    "UnsupportedVersion",
    "add_metrics",
    "agreement",
    "call",
    "compare",
    "list_methods",
    "pl",
    "read",
    "separation",
    "simulate_droplets",
    "stored",
    "unspliced_fraction",
]
