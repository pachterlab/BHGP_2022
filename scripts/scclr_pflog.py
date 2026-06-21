#!/usr/bin/env python3
"""Helpers for the current scclr PFlog transform.

PFlog is the shifted CLR

    center(log1p(4 * alpha * x))

computed by scclr.  The earlier depth-targeted helper
``center(log1p(4 * alpha * x * sbar / s_i))`` is intentionally not used here.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp


def import_scclr(extra_path: str | os.PathLike[str] | None = None):
    path = extra_path or os.environ.get("SC_CLR_PATH")
    if path:
        sys.path.insert(0, str(Path(path)))
    import scclr

    _assert_correct_pflog(scclr)
    return scclr


def _assert_correct_pflog(scclr) -> None:
    X = sp.csr_matrix([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]])
    res = scclr.normalize(X, alpha=0.5)
    got = res.sparse.toarray()
    expected = np.array(
        [
            [np.log1p(2.0), 0.0, np.log1p(6.0)],
            [0.0, np.log1p(4.0), 0.0],
        ]
    )
    if not np.allclose(got, expected):
        src = getattr(scclr, "__file__", "<unknown>")
        raise RuntimeError(
            "Imported scclr does not implement current PFlog "
            f"center(log1p(4*alpha*x)); got incompatible module at {src}. "
            "Set SC_CLR_PATH=/path/to/scclr/python or update the environment."
        )


def normalize_pflog(X, alpha: float | None = None, scclr_path=None):
    scclr = import_scclr(scclr_path)
    if alpha is None:
        return scclr.normalize(X, target="auto")
    return scclr.normalize(X, alpha=float(alpha))


def to_dense(sclr, rows=None, dtype=np.float32) -> np.ndarray:
    L = sclr.sparse if rows is None else sclr.sparse[rows]
    center = np.asarray(sclr.row_center, dtype=float).ravel()
    if rows is not None:
        center = center[rows]
    return (L.toarray() - center[:, None]).astype(dtype, copy=False)


def gene_var(sclr) -> np.ndarray:
    """Per-gene variance of Z = L - row_center without densifying Z."""
    L = sclr.sparse
    m = np.asarray(sclr.row_center, dtype=float).ravel()
    n = L.shape[0]
    mean_l = np.asarray(L.sum(0)).ravel() / n
    var_l = np.asarray(L.multiply(L).sum(0)).ravel() / n - mean_l**2
    cov_lm = np.asarray(L.T.dot(m)).ravel() / n - mean_l * m.mean()
    return var_l + m.var() - 2.0 * cov_lm


def cov_gene(sclr) -> float:
    gvar = gene_var(sclr)
    return float(np.sqrt(np.var(gvar)) / np.mean(gvar))
