#!/usr/bin/env python3
"""Compare the p=11 square-circle operator with the exact Max+ Phi.

This is a finite representation-theory reconnaissance tool.  It rebuilds
the same deterministic orthonormal basis of

    Z = {W: PWP=W, diag(W)=0}

used by the full p=11 spectrum computation, constructs the square-circle
evaluation matrix from the complete second dual shell, and tests how its
two nonzero eigenspaces decompose under ``Phi``.  No finite p=11 spectrum is
promoted to a general-p theorem by this script.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def canonical(vector: np.ndarray) -> tuple[int, ...]:
    nonzero = np.flatnonzero(vector)
    if len(nonzero) == 0:
        raise ValueError("cannot orient the zero vector")
    oriented = vector if vector[nonzero[0]] > 0 else -vector
    return tuple(int(value) for value in oriented)


def clusters(values: np.ndarray, tolerance: float = 1e-7) -> list[dict]:
    groups: list[list[float]] = []
    for value in np.sort(values):
        if not groups or abs(value - groups[-1][-1]) > tolerance * max(
            1.0, abs(value)
        ):
            groups.append([float(value)])
        else:
            groups[-1].append(float(value))
    return [
        {
            "eigenvalue": float(sum(group) / len(group)),
            "multiplicity": len(group),
        }
        for group in groups
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shell", type=Path, required=True)
    parser.add_argument("--phi", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    p = 11
    n = p * p + 1
    started = time.perf_counter()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    archive = np.load(args.shell)
    second = archive["ambient_numerator"][
        archive["scaled_norm"] == 2 * (p - 1)
    ]
    complete = {canonical(vector): vector for vector in second}
    point_pairs = set()
    for i in range(n):
        for j in range(i + 1, n):
            z = np.zeros(n, dtype=np.int64)
            z[i] = 1
            z[j] = -int(C[i, j])
            point_pairs.add(canonical(p * z + C @ z))
    circle_keys = complete.keys() - point_pairs
    circles = np.asarray([complete[key] // 2 for key in circle_keys], dtype=np.float64)
    if circles.shape != (p * n // 2, n):
        raise ArithmeticError(f"unexpected circle matrix {circles.shape}")
    print(f"circle vectors: {circles.shape}", flush=True)

    eigenvalues, eigenvectors = np.linalg.eigh(C.astype(np.float64))
    Vplus = eigenvectors[:, np.isclose(eigenvalues, p, atol=1e-8)]
    d = Vplus.shape[1]
    upper = np.triu_indices(d)
    off_diagonal = upper[0] != upper[1]

    diagonal_constraints = np.empty((n, len(upper[0])), dtype=np.float64)
    for index, row in enumerate(Vplus):
        outer = np.outer(row, row)[upper]
        outer[off_diagonal] *= np.sqrt(2.0)
        diagonal_constraints[index] = outer
    _left, singular, vh = np.linalg.svd(diagonal_constraints, full_matrices=True)
    rank = int(np.count_nonzero(singular > 1e-9))
    Zbasis = vh[rank:].T
    if Zbasis.shape != (d * (d + 1) // 2, n * (n - 6) // 8):
        raise ArithmeticError(f"unexpected Z basis {Zbasis.shape}")

    coordinates = circles @ Vplus
    outer = coordinates[:, :, None] * coordinates[:, None, :]
    features = outer[:, upper[0], upper[1]]
    features[:, off_diagonal] *= np.sqrt(2.0)
    evaluations = features @ Zbasis
    gram = evaluations @ evaluations.T
    gram = (gram + gram.T) / 2.0
    gram_values, gram_vectors = np.linalg.eigh(gram)
    expected_low = p**3 * (p - 1)
    expected_high = p**3 * (p + 1)
    low_indices = np.flatnonzero(np.isclose(gram_values, expected_low, atol=1e-5))
    high_indices = np.flatnonzero(np.isclose(gram_values, expected_high, atol=1e-5))
    zero_indices = np.flatnonzero(np.abs(gram_values) < 1e-5)
    if (len(zero_indices), len(low_indices), len(high_indices)) != (122, 305, 244):
        raise ArithmeticError(
            "unexpected circle Gram multiplicities: "
            f"{len(zero_indices), len(low_indices), len(high_indices)}"
        )

    low_basis = evaluations.T @ gram_vectors[:, low_indices] / np.sqrt(expected_low)
    high_basis = evaluations.T @ gram_vectors[:, high_indices] / np.sqrt(expected_high)
    phi = np.load(args.phi)
    if phi.shape != (Zbasis.shape[1], Zbasis.shape[1]):
        raise ArithmeticError(f"unexpected Phi shape {phi.shape}")

    def compression_report(basis: np.ndarray) -> dict:
        compressed = basis.T @ phi @ basis
        residual = phi @ basis - basis @ compressed
        values = np.linalg.eigvalsh((compressed + compressed.T) / 2.0)
        return {
            "dimension": basis.shape[1],
            "orthonormality_error": float(
                np.linalg.norm(basis.T @ basis - np.eye(basis.shape[1]), ord=2)
            ),
            "phi_preservation_error": float(np.linalg.norm(residual, ord=2)),
            "phi_preservation_relative": float(
                np.linalg.norm(residual, ord=2) / np.linalg.norm(phi, ord=2)
            ),
            "phi_spectrum": clusters(values, tolerance=2e-7),
            "phi_min": float(values[0]),
            "phi_max": float(values[-1]),
        }

    low = compression_report(low_basis)
    high = compression_report(high_basis)
    result = {
        "experiment": "r1_p11_circle_phi_intertwiner",
        "status": "finite_p11_reconnaissance_not_general_proof",
        "p": p,
        "circle_count": len(circles),
        "Z_dimension": Zbasis.shape[1],
        "circle_gram_spectrum": clusters(gram_values, tolerance=1e-7),
        "predicted_circle_gram_spectrum": [
            {"eigenvalue": 0, "multiplicity": n},
            {"eigenvalue": expected_low, "multiplicity": d * (p - 1) // 2},
            {"eigenvalue": expected_high, "multiplicity": d * (p - 3) // 2},
        ],
        "low_circle_eigenspace": low,
        "high_circle_eigenspace": high,
        "elapsed_seconds": time.perf_counter() - started,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
