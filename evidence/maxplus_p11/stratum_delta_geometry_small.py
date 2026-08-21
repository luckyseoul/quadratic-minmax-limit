#!/usr/bin/env python3
"""Compute profile-stratum four-point residual norms at p=5 or p=7.

This is the small-prime cross-check for ``stratum_delta_geometry.py``.  It
uses the independent exhaustive Max+ caches, classifies every eps=+1 vector
by its number of active square directions, and reports the residual geometry
of the resulting probability measures.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/storage/e1work/scripts")
from kgen import square_coords  # type: ignore


def m4_flat_norm_sq(p: int) -> Fraction:
    return Fraction(
        (p - 1) * (p + 1) * (p * p + 1) * (3 * p * p + 17),
        24 * (p * p - 5),
    )


def active_counts(y: np.ndarray, p: int) -> np.ndarray:
    _dirs, _forms, coords = square_coords(p)
    eps = y[:, 0]
    yf = y[:, 1:]
    out = np.zeros(len(y), dtype=np.int8)
    for tmap in np.asarray(coords, dtype=np.int64):
        line_sums = np.stack(
            [yf[:, tmap == value].sum(axis=1) for value in range(p)], axis=1
        )
        out += np.any(line_sums != eps[:, None], axis=1)
    return out


def pair_gram(y: np.ndarray, iu: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    features = y[:, iu[0]] * y[:, iu[1]]
    # BLAS float64 is exact here: every output is an integer of magnitude at
    # most 5726, far below the 53-bit mantissa limit.
    ff = features.astype(np.float64)
    return np.rint(ff.T @ ff).astype(np.int64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=int, choices=(5, 7))
    args = parser.parse_args()
    p = args.p
    n = p * p + 1
    source = Path(f"/tmp/maxplus_p{p}.npy")
    yall = np.load(source)
    yall = np.rint(yall[yall[:, 0] == 1]).astype(np.int8)
    kvals = active_counts(yall, p)
    labels = tuple(f"k{k}" for k in sorted(set(map(int, kvals))))
    groups = {label: yall[kvals == int(label[1:])] for label in labels}

    iu = np.triu_indices(n, 1)
    pair_id = -np.ones((n, n), dtype=np.int32)
    pair_id[iu] = np.arange(len(iu[0]), dtype=np.int32)
    pair_id[(iu[1], iu[0])] = pair_id[iu]
    quads = np.asarray(list(itertools.combinations(range(n), 4)), dtype=np.int16)
    a, b, c, d = (quads[:, j].astype(np.int32) for j in range(4))
    pairings = (
        (pair_id[a, b], pair_id[c, d]),
        (pair_id[a, c], pair_id[b, d]),
        (pair_id[a, d], pair_id[b, c]),
    )

    moments: dict[str, np.ndarray] = {}
    consistency: dict[str, int] = {}
    for label in labels:
        gram = pair_gram(groups[label], iu)
        values = [gram[idx] for idx in pairings]
        consistency[label] = int(
            max(np.max(np.abs(values[0] - values[1])), np.max(np.abs(values[0] - values[2])))
        )
        if consistency[label] != 0:
            raise RuntimeError(f"pairing mismatch for {label}")
        moments[label] = values[0].astype(np.float64) / len(groups[label])

    flat = float(m4_flat_norm_sq(p))
    residual = np.empty((len(labels), len(labels)), dtype=np.float64)
    for i, left in enumerate(labels):
        for j, right in enumerate(labels):
            residual[i, j] = float(np.dot(moments[left], moments[right])) - flat

    counts = np.asarray([len(groups[k]) for k in labels], dtype=np.int64)
    weights = counts / counts.sum()
    global_delta_sq = float(weights @ residual @ weights)
    norms = np.sqrt(np.maximum(np.diag(residual), 0.0))
    cosines = residual / np.outer(norms, norms)
    report = {
        "p": p,
        "source": str(source),
        "counts": dict(zip(labels, counts.tolist())),
        "weights": dict(zip(labels, weights.tolist())),
        "pairing_consistency_max_integer_error": consistency,
        "m4_flat_norm_sq": str(m4_flat_norm_sq(p)),
        "delta_norm_sq_by_stratum": dict(zip(labels, np.diag(residual).tolist())),
        "delta_inner_product_matrix_order": labels,
        "delta_inner_product_matrix": residual.tolist(),
        "delta_cosine_matrix": cosines.tolist(),
        "weighted_global_delta_norm_sq": global_delta_sq,
    }
    output = Path(f"/mnt/storage/e1work/maxplus_p11/stratum_delta_geometry_p{p}.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
