#!/usr/bin/env python3
"""Exceptional-constituent projection shells at p=5,7 and a p=11 witness.

This is a diagnostic, not a general theorem.  It records two failed pointwise
routes to QVAR:

* p=7 has a zero exceptional-projection shell;
* even after excluding zero, p=11 has a nonzero shell below the mean target
  ``3n``.  Census row 11_453_817 has squared projection norm ``4304/15``,
  while ``3n=366``.

The p=11 run uses the untracked complete census and Phi matrix under
``E1WORK_P11`` (default ``/mnt/storage/e1work/maxplus_p11``).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15588 import maxplus, phi_matrix, z_basis  # noqa: E402


def _clusters(values: np.ndarray, tol: float = 1e-6) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    lo = 0
    for hi in range(1, len(values) + 1):
        if hi == len(values) or values[hi] - values[hi - 1] > tol:
            out.append((lo, hi))
            lo = hi
    return out


def _exceptional_forms(p: int, Phi: np.ndarray) -> tuple[np.ndarray, float]:
    values, vectors = np.linalg.eigh(Phi)
    d = (p * p + 1) // 2
    lo, hi = next(block for block in _clusters(values) if block[1] - block[0] == d)
    forms = np.einsum(
        "tr,tij->rij", vectors[:, lo:hi], z_basis(p), optimize=True
    )
    return forms, float(values[lo:hi].mean())


def _norms(rows: np.ndarray, forms: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.float64)
    coordinates = np.einsum("ai,rij,aj->ar", rows, forms, rows, optimize=True)
    return np.sum(coordinates * coordinates, axis=1)


def _rational(x: float) -> Fraction:
    return Fraction(float(x)).limit_denominator(100_000)


def _full_small_prime(p: int) -> dict:
    rows = maxplus(p)
    forms, eigenvalue = _exceptional_forms(p, phi_matrix(p))
    norms = _norms(rows, forms)
    rounded = np.round(norms, 8)
    values, counts = np.unique(rounded, return_counts=True)
    shells = [
        {"norm_sq": str(_rational(value)), "count": int(count)}
        for value, count in zip(values, counts)
    ]
    nonzero = norms[norms > 1e-7]
    quantum = Fraction(128, p * p - 1)
    return {
        "p": p,
        "n": p * p + 1,
        "row_count": len(rows),
        "lambda_exceptional": eigenvalue,
        "mean_projection_norm_sq": float(norms.mean()),
        "target_3n": 3 * (p * p + 1),
        "shells": shells,
        "min_nonzero_norm_sq": str(_rational(nonzero.min())),
        "observed_quantum": str(quantum),
        "all_shells_are_quantum_multiples": all(
            (_rational(value) / quantum).denominator == 1 for value in values
        ),
    }


def _p11_witness() -> dict:
    work = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    phi_path = work / "phiZ_p11.npy"
    rows_path = work / "maxplus_p11_eps1.npy"
    if not phi_path.is_file() or not rows_path.is_file():
        return {
            "p": 11,
            "skipped": True,
            "reason": f"missing {phi_path} or {rows_path}",
        }

    index = 11_453_817
    rows = np.load(rows_path, mmap_mode="r")
    row = np.asarray(rows[index : index + 1])
    forms, eigenvalue = _exceptional_forms(11, np.load(phi_path))
    norm = _rational(_norms(row, forms)[0])
    quantum = Fraction(128, 11 * 11 - 1)
    return {
        "p": 11,
        "n": 122,
        "census_index": index,
        "row_sha256": hashlib.sha256(row.tobytes()).hexdigest(),
        "lambda_exceptional": eigenvalue,
        "norm_sq": str(norm),
        "target_3n": 366,
        "nonzero": norm > 0,
        "below_target": norm < 366,
        "observed_quantum": str(quantum),
        "is_quantum_multiple": (norm / quantum).denominator == 1,
    }


def main() -> dict:
    small = {str(p): _full_small_prime(p) for p in (5, 7)}
    witness = _p11_witness()
    result = {
        "title": "Exceptional projection shell diagnostic",
        "small_prime_full_shells": small,
        "p11_nonzero_below_3n_witness": witness,
        "pointwise_exceptional_floor_false": any(
            shell["norm_sq"] == "0"
            for shell in small["7"]["shells"]
        ),
        "nonzero_shell_ge_3n_false": bool(
            not witness.get("skipped", False) and witness["below_target"]
        ),
        "interpretation": (
            "QVAR is an ensemble mean bound.  Neither a pointwise floor nor a "
            "lower bound of 3n on every nonzero exceptional shell can prove it."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(f"wrote {path}")
    return result


if __name__ == "__main__":
    main()
