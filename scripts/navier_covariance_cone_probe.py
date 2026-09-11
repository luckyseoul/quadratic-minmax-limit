#!/usr/bin/env python3
"""Finite diagnostic for the Navier-style PSD covariance cone.

For a signing A and M=A/sqrt(n), the all-law source proof uses the two
correlations built from |M|+M and |M|-M.  This script asks whether choosing a
different common PSD majorizer P with P+M,P-M positive semidefinite can
improve the diagonal-dependent arcsine baseline.  It is deliberately a
finite floating-point scout: any improvement is a candidate for an exact
all-orders lemma, while no improvement is not a theorem.
"""

from __future__ import annotations

import argparse
import json
import math

import cvxpy as cp
import numpy as np


def signing(order: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    upper = rng.choice((-1.0, 1.0), size=(order, order))
    matrix = np.triu(upper, 1)
    return matrix + matrix.T


def score(diagonal: np.ndarray) -> float:
    return float(sum(
        1.0 / math.sqrt(diagonal[i] * diagonal[j])
        for i in range(len(diagonal)) for j in range(i + 1, len(diagonal))
    ))


def abs_majorizer(matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.abs(values)) @ vectors.T


def weighted_majorizer(matrix: np.ndarray, weights: np.ndarray) -> np.ndarray:
    order = len(matrix)
    variable = cp.Variable((order, order), symmetric=True)
    problem = cp.Problem(
        cp.Minimize(weights @ cp.diag(variable)),
        [variable + matrix >> 0, variable - matrix >> 0],
    )
    problem.solve(solver="CLARABEL", tol_gap_abs=1e-9, tol_feas=1e-9,
                  tol_gap_rel=1e-9, max_iter=500)
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(problem.status)
    value = np.asarray(variable.value, dtype=float)
    min_eigenvalue = min(np.linalg.eigvalsh(value + matrix).min(),
                         np.linalg.eigvalsh(value - matrix).min())
    if min_eigenvalue < -2e-6:
        raise RuntimeError(("PSD residual", min_eigenvalue))
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=10)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--rays", type=int, default=16)
    arguments = parser.parse_args()
    if arguments.order < 3 or arguments.rays < 1:
        raise ValueError("order >= 3 and rays >= 1 required")

    matrix = signing(arguments.order, arguments.seed) / math.sqrt(arguments.order)
    canonical = abs_majorizer(matrix)
    canonical_score = score(np.diag(canonical))
    rng = np.random.default_rng(arguments.seed ^ 0x4E4156494552)
    candidates = []
    for ray in range(arguments.rays):
        weights = np.exp(rng.uniform(-2.0, 2.0, arguments.order))
        candidate = weighted_majorizer(matrix, weights)
        candidates.append(score(np.diag(candidate)))
    best = max(candidates)
    print(json.dumps({
        "kind": "finite_floating_cone_scout",
        "order": arguments.order,
        "seed": arguments.seed,
        "rays": arguments.rays,
        "abs_majorizer_score": canonical_score,
        "best_weighted_score": best,
        "relative_improvement": best / canonical_score - 1.0,
        "candidate_scores": candidates,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
