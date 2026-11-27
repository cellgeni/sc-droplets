"""Fit scikit-learn's GaussianMixture to a dataset's barcodes, as the `gmm` and `gmm-2d` rules do.

    python gmm.py FEATURES.tsv OUT_DIR

FEATURES.tsv is `barcode<TAB>umis<TAB>genes` (make_features.py). The candidates have at least 100
UMIs; their log10 (the C library's, Python's math module) is the one dimension, with log10 genes
the second. For 2 and 3 components, at scikit-learn's default tolerance and converged (tol 1e-8,
max_iter 10000), seed 0. Writes OUT_DIR/gmm.tsv (`dimensions<TAB>components<TAB>tol<TAB>max_iter
<TAB>means<TAB>covariances<TAB>weights<TAB>iterations<TAB>cells<TAB>log-likelihood`, numbers as
Python's repr, comma-separated; cells are those above 0.5 for the component with the highest
mean log10 UMIs) and OUT_DIR/versions.txt.
"""

import math
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import sklearn
from sklearn.mixture import GaussianMixture

rows = [line.split("\t") for line in Path(sys.argv[1]).read_text().splitlines()]
umis = np.array([float(r[1]) for r in rows])
genes = np.array([float(r[2]) for r in rows])
keep = umis >= 100
x = np.array([math.log10(v) for v in umis[keep].tolist()])
y = np.array([math.log10(v) for v in genes[keep].tolist()])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)


def numbers(values):
    return ",".join(repr(float(v)) for v in np.ravel(values))


lines = []
for dimensions, data in ((1, x.reshape(-1, 1)), (2, np.column_stack([x, y]))):
    for components in (2, 3):
        for tol, max_iter in ((1e-3, 100), (1e-8, 10000)):
            model = GaussianMixture(n_components=components, covariance_type="full", tol=tol, max_iter=max_iter,
                                    random_state=0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model.fit(data)
            top = int(np.argmax(model.means_[:, 0]))
            cells = int((model.predict_proba(data)[:, top] > 0.5).sum())
            loglik = float(model.score(data) * data.shape[0])
            lines.append("\t".join([str(dimensions), str(components), repr(tol), str(max_iter), numbers(model.means_),
                                    numbers(model.covariances_), numbers(model.weights_), str(model.n_iter_),
                                    str(cells), repr(loglik)]))
            print(dimensions, components, tol, model.n_iter_, cells, round(loglik, 4))
(out / "gmm.tsv").write_text("\n".join(lines) + "\n")
(out / "versions.txt").write_text(
    f"scikit-learn {sklearn.__version__}\nnumpy {np.__version__}\npython {platform.python_version()}\n"
    f"platform {platform.system()}-{platform.machine()}\n"
)
