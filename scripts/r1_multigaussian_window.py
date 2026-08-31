#!/usr/bin/env python3
"""Finite multi-Gaussian anti-cancellation probe for the R1 theta route.

This is a laboratory tool, not a proposition and not an R1 proof.  It uses
PARI/GP to enumerate exact low shells of

    L* = P Z^n,   P=(I+C/p)/2,

forms the degree-four harmonic shell operator from Proposition 15.631, and
tests signed combinations of Gaussian parameters.  The combination is
normalized to have coefficient one on the first primal odd-coset shell.
An SDP then asks whether the *enumerated* dual shells have a positive
operator margin while early higher primal shells are suppressed.

Two controls prevent a false positive:

* the known first dual shell is separated from the higher-shell operator;
* the SDP is rerun with that first-shell coefficient forced to zero.

Any positive margin is finite-shell numerical evidence only.  A proof would
still need certified bounds for both omitted dual shells and the remaining
primal odd-coset tail.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def rank_L(p: int) -> int:
    return (p * p + 1) // 2


def lattice_det(p: int) -> int:
    return 2 * p ** (((p + 1) // 2) ** 2)


def paley_conference(p: int) -> np.ndarray:
    return np.rint(paley_conference_prime_power(p)).astype(np.int64)


def gp_matrix(A: np.ndarray) -> str:
    rows = [",".join(str(int(x)) for x in row) for row in A]
    return "[" + ";".join(rows) + "]"


def enumerate_dual_shells(
    p: int,
    scaled_bound: int,
    max_half_vectors: int = 1_000_000,
) -> dict[int, np.ndarray]:
    """Enumerate one vector from each +/- pair with 2p||u||^2 <= bound.

    PARI works in dual coordinates.  If ``B`` is a saturated basis of L,
    ``G=B^T B``, and ``Q=2p G^-1``, then a dual coordinate vector ``k`` maps
    to the ambient numerator ``B Q k`` of ``u=(B Q k)/(2p)``.
    """
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    program = f"""
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
Q={2*p}*G^-1;
if(denominator(Q)!=1,error("nonintegral scaled dual Gram"));
DM=qfminim(Q,{int(scaled_bound)},{int(max_half_vectors)});
V=DM[3];
print("COUNT=",DM[1]);
print("STORED=",matsize(V)[2]);
print("MAXNORM=",DM[2]);
for(j=1,matsize(V)[2],k=V[,j];s=k~*Q*k;unum=B*Q*k;print("VECTOR=",s,"|",Vec(unum)));
quit;
"""
    proc = subprocess.run(
        ["gp", "-q", "-s", "1G"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in proc.stderr:
        raise RuntimeError(f"PARI/GP failed: {proc.stderr.strip()}")

    count = stored = maxnorm = None
    rows: dict[int, list[list[int]]] = {}
    for line in proc.stdout.splitlines():
        if line.startswith("COUNT="):
            count = int(line.split("=", 1)[1])
        elif line.startswith("STORED="):
            stored = int(line.split("=", 1)[1])
        elif line.startswith("MAXNORM="):
            maxnorm = int(line.split("=", 1)[1])
        elif line.startswith("VECTOR="):
            label, raw = line.split("|", 1)
            scaled_norm = int(label.split("=", 1)[1])
            vector = ast.literal_eval(raw)
            rows.setdefault(scaled_norm, []).append(vector)
    if count is None or stored is None or maxnorm is None:
        raise RuntimeError(f"could not parse PARI output: {proc.stdout[-2000:]}")
    if count != 2 * stored:
        raise RuntimeError(
            "PARI enumeration was truncated or did not return +/- pairs: "
            f"count={count}, stored={stored}"
        )
    if stored >= max_half_vectors:
        raise RuntimeError("PARI reached the vector-storage limit")

    out = {
        s: np.asarray(vectors, dtype=np.int64) / float(2 * p)
        for s, vectors in rows.items()
    }
    for s, vectors in out.items():
        norms = np.einsum("vi,vi->v", vectors, vectors)
        if not np.allclose(norms, s / (2 * p), atol=1e-10, rtol=1e-10):
            raise ArithmeticError(f"ambient norm check failed on shell {s}")
    return out


def admissible_W_basis(
    C: np.ndarray, p: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ambient W basis, +p coordinates, and intrinsic W basis."""
    eigenvalues, eigenvectors = np.linalg.eigh(C.astype(np.float64))
    U = eigenvectors[:, np.isclose(eigenvalues, p, atol=1e-8)]
    d = U.shape[1]
    if d != rank_L(p):
        raise ArithmeticError("wrong +p eigenspace rank")

    sym: list[np.ndarray] = []
    for a in range(d):
        E = np.zeros((d, d), dtype=np.float64)
        E[a, a] = 1.0
        sym.append(E)
        for b in range(a + 1, d):
            E = np.zeros((d, d), dtype=np.float64)
            E[a, b] = E[b, a] = 1.0 / math.sqrt(2.0)
            sym.append(E)
    sym_array = np.stack(sym)
    diagonal_map = np.einsum("ia,kab,ib->ik", U, sym_array, U, optimize=True)
    _u, singular, vh = np.linalg.svd(diagonal_map, full_matrices=True)
    tolerance = max(diagonal_map.shape) * singular[0] * np.finfo(float).eps
    rank = int(np.count_nonzero(singular > tolerance))
    coefficient_nullspace = vh[rank:]
    intrinsic = np.einsum(
        "rk,kab->rab", coefficient_nullspace, sym_array, optimize=True
    )
    ambient = np.einsum("ia,rab,jb->rij", U, intrinsic, U, optimize=True)

    gram = np.einsum("rij,sij->rs", ambient, ambient, optimize=True)
    if np.max(np.abs(gram - np.eye(len(ambient)))) > 2e-10:
        raise ArithmeticError("W basis is not Frobenius orthonormal")
    if np.max(np.abs(np.diagonal(ambient, axis1=1, axis2=2))) > 2e-10:
        raise ArithmeticError("W basis violates the zero-diagonal constraint")
    return ambient, U, intrinsic


def harmonic_operator(
    vectors: np.ndarray,
    W_basis: np.ndarray,
    eigenspace_basis: np.ndarray,
    intrinsic_W_basis: np.ndarray,
    chunk_size: int = 4096,
) -> np.ndarray:
    """Matrix T with alpha^T T alpha = sum_x H_{sum alpha_i W_i}(x)."""
    if len(vectors) == 0:
        return np.zeros((len(W_basis), len(W_basis)), dtype=np.float64)
    X = np.asarray(vectors, dtype=np.float64)
    n = X.shape[1]
    d = n // 2
    coordinates = X @ eigenspace_basis
    flat_W = intrinsic_W_basis.reshape(len(W_basis), d * d)

    first = np.zeros((len(W_basis), len(W_basis)), dtype=np.float64)
    for lo in range(0, len(X), chunk_size):
        z = coordinates[lo : lo + chunk_size]
        quadratic_features = np.einsum("vi,vj->vij", z, z, optimize=True)
        quadratic_values = quadratic_features.reshape(len(z), d * d) @ flat_W.T
        first += quadratic_values.T @ quadratic_values

    norm_sq = np.einsum("vi,vi->v", X, X)
    weighted_covariance = coordinates.T @ (norm_sq[:, None] * coordinates)
    W_covariance = np.einsum(
        "kij,jm->kim", intrinsic_W_basis, weighted_covariance, optimize=True
    )
    second = W_covariance.reshape(len(W_basis), d * d) @ flat_W.T
    second = (second + second.T) / 2.0
    identity = np.eye(len(W_basis))
    third = float(np.dot(norm_sq, norm_sq)) * identity
    out = (
        first
        - 4.0 * second / (d + 4)
        + 2.0 * third / ((d + 2) * (d + 4))
    )
    return (out + out.T) / 2.0


def phased_dual_shell_operators(
    p: int,
    shells: dict[int, np.ndarray],
    W_basis: np.ndarray,
    eigenspace_basis: np.ndarray,
    intrinsic_W_basis: np.ndarray,
) -> dict[int, np.ndarray]:
    """Operators for sum phase(u) H_W(u/2) on complete +/- shells."""
    out = {}
    for scaled_norm, half_shell in shells.items():
        phase = -1.0 if scaled_norm & 1 else 1.0
        # qfminim stores one vector per +/- pair; H is even and degree four.
        out[scaled_norm] = phase * harmonic_operator(
            half_shell,
            W_basis,
            eigenspace_basis,
            intrinsic_W_basis,
        ) / 8.0
    return out


def maxplus_shell(p: int) -> np.ndarray:
    path = Path(f"/tmp/maxplus_p{p}.npy")
    if not path.is_file():
        from enum_maxplus import enum_maxplus

        enum_maxplus(p, workers=min(8, max(1, (os.cpu_count() or 2) - 1)))
    return np.load(path).astype(np.float64)


def primal_lattice_minimum_shell(
    p: int,
    max_half_vectors: int = 100_000,
) -> np.ndarray:
    """Enumerate the complete signed norm-(p+1) shell of L exactly."""
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    program = f"""
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
DM=qfminim(G,{p + 1},{int(max_half_vectors)});
V=DM[3];
print("COUNT=",DM[1]);
print("STORED=",matsize(V)[2]);
print("MAXNORM=",DM[2]);
for(j=1,matsize(V)[2],k=V[,j];y=B*k;print("VECTOR=",k~*G*k,"|",Vec(y)));
quit;
"""
    proc = subprocess.run(
        ["gp", "-q", "-s", "1G"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in proc.stderr:
        raise RuntimeError(f"PARI/GP failed: {proc.stderr.strip()}")
    count = stored = maxnorm = None
    vectors: list[list[int]] = []
    norms: list[int] = []
    for line in proc.stdout.splitlines():
        if line.startswith("COUNT="):
            count = int(line.split("=", 1)[1])
        elif line.startswith("STORED="):
            stored = int(line.split("=", 1)[1])
        elif line.startswith("MAXNORM="):
            maxnorm = int(line.split("=", 1)[1])
        elif line.startswith("VECTOR="):
            label, raw = line.split("|", 1)
            norms.append(int(label.split("=", 1)[1]))
            vectors.append(ast.literal_eval(raw))
    if count is None or stored is None or maxnorm is None:
        raise RuntimeError(f"could not parse primal qfminim output: {proc.stdout[-2000:]}")
    if count != 2 * stored or stored >= max_half_vectors:
        raise RuntimeError(
            "primal minimum-shell enumeration truncated: "
            f"count={count}, stored={stored}"
        )
    if maxnorm != p + 1 or any(norm != p + 1 for norm in norms):
        raise ArithmeticError(
            f"L has a vector below the expected norm {p + 1}: max={maxnorm}"
        )
    positive = np.asarray(vectors, dtype=np.float64)
    if not np.all(np.einsum("vi,vi->v", positive, positive) == p + 1):
        raise ArithmeticError("primal minimum-shell ambient norm check failed")
    return np.concatenate([positive, -positive], axis=0)


def two_L_volume(p: int) -> float:
    return (2.0 ** rank_L(p)) * math.sqrt(float(lattice_det(p)))


def transformed_low_operator(
    p: int,
    t: float,
    shell_operators: dict[int, np.ndarray],
) -> np.ndarray:
    """Operator for e^(pi*t*n) times the enumerated dual transform."""
    d = rank_L(p)
    n = p * p + 1
    factor = math.exp(math.pi * t * n) * t ** (-(d / 2.0 + 4.0))
    factor /= two_L_volume(p)
    out = np.zeros_like(next(iter(shell_operators.values())))
    for scaled_norm, operator in shell_operators.items():
        weight = math.exp(-math.pi * scaled_norm / (8.0 * p * t))
        out += weight * operator
    return factor * out


def transformed_first_scalar(p: int, t: float) -> float:
    d = rank_L(p)
    n = p * p + 1
    factor = math.exp(math.pi * t * n) * t ** (-(d / 2.0 + 4.0))
    factor /= two_L_volume(p)
    return (
        factor
        * math.exp(-math.pi / (8.0 * t))
        / (8.0 * (d + 2))
    )


def spectrum_clusters(matrix: np.ndarray, tolerance: float = 1e-7) -> list[dict]:
    """Cluster a symmetric numerical spectrum by absolute tolerance."""
    values = np.linalg.eigvalsh(matrix)
    clusters: list[list[float]] = []
    for value in values:
        if not clusters or abs(value - clusters[-1][-1]) > tolerance:
            clusters.append([float(value)])
        else:
            clusters[-1].append(float(value))
    return [
        {"eigenvalue": float(sum(cluster) / len(cluster)), "multiplicity": len(cluster)}
        for cluster in clusters
    ]


def exceptional_channel(
    primal_operator: np.ndarray,
    shell_operators: dict[int, np.ndarray],
    expected_dimension: int,
    tolerance: float = 1e-7,
) -> dict:
    """Extract the PSL exceptional block from the first primal operator.

    In the audited p=5 and p=7 cases the exceptional Weil constituent is the
    unique eigenspace of the first-primal-shell operator having dimension
    ``rank_L(p)``.  This is a numerical block identifier, not a general-p
    theorem.  We explicitly check that every enumerated dual-shell operator
    preserves the block and is scalar on it.
    """
    values, vectors = np.linalg.eigh(primal_operator)
    clusters: list[list[int]] = []
    for index, value in enumerate(values):
        if not clusters or abs(value - values[clusters[-1][-1]]) > tolerance:
            clusters.append([index])
        else:
            clusters[-1].append(index)
    candidates = [indices for indices in clusters if len(indices) == expected_dimension]
    if len(candidates) != 1:
        raise ArithmeticError(
            "could not identify a unique rank-dimensional exceptional block: "
            f"multiplicities={[len(indices) for indices in clusters]}"
        )
    indices = candidates[0]
    basis = vectors[:, indices]
    projector = basis @ basis.T
    identity = np.eye(expected_dimension)
    shell_scalars: dict[int, float] = {}
    scalar_errors: dict[int, float] = {}
    preservation_errors: dict[int, float] = {}
    for scaled_norm, operator in shell_operators.items():
        compression = basis.T @ operator @ basis
        scalar = float(np.trace(compression) / expected_dimension)
        shell_scalars[scaled_norm] = scalar
        scalar_errors[scaled_norm] = float(
            np.linalg.norm(compression - scalar * identity, ord=2)
        )
        preservation_errors[scaled_norm] = float(
            np.linalg.norm((np.eye(len(operator)) - projector) @ operator @ basis, ord=2)
        )
    max_error = max([*scalar_errors.values(), *preservation_errors.values()])
    if max_error > 2e-7:
        raise ArithmeticError(f"exceptional block scalarity failed: {max_error}")
    return {
        "dimension": expected_dimension,
        "first_primal_eigenvalue": float(np.mean(values[indices])),
        "first_primal_multiplicity": len(indices),
        "dual_shell_scalars": shell_scalars,
        "max_scalarity_error": max(scalar_errors.values()),
        "max_preservation_error": max(preservation_errors.values()),
    }


def transformed_scalar_channel(
    p: int,
    t: float,
    shell_scalars: dict[int, float],
) -> float:
    """Transformed enumerated dual sum on one invariant scalar channel."""
    d = rank_L(p)
    n = p * p + 1
    factor = math.exp(math.pi * t * n) * t ** (-(d / 2.0 + 4.0))
    factor /= two_L_volume(p)
    return factor * sum(
        math.exp(-math.pi * scaled_norm / (8.0 * p * t)) * scalar
        for scaled_norm, scalar in shell_scalars.items()
    )


def solve_scalar_coercive_window(
    t_values: np.ndarray,
    low_values: np.ndarray,
    first_scalars: np.ndarray,
    primal_cancellations: int,
    leakage: float,
    l1_cap: float,
) -> dict:
    """Maximize first-shell scalar minus the absolute known scalar tail."""
    from scipy.optimize import linprog

    count = len(t_values)
    # Variables are (a_1,...,a_J, u_1,...,u_J, rho), rho>=|tail.a|.
    size = 2 * count + 1
    rho_index = 2 * count
    tail_values = low_values - first_scalars
    scale = max(float(np.max(np.abs(first_scalars))), float(np.max(np.abs(tail_values))))
    if not math.isfinite(scale) or scale <= 0.0:
        raise ArithmeticError("invalid exceptional-channel scale")
    first_normalized = first_scalars / scale
    tail_normalized = tail_values / scale
    objective = np.zeros(size)
    objective[:count] = -first_normalized
    objective[rho_index] = 1.0
    upper_rows: list[np.ndarray] = []
    upper_bounds: list[float] = []
    for j in range(count):
        row = np.zeros(size)
        row[j], row[count + j] = 1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
        row = np.zeros(size)
        row[j], row[count + j] = -1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
    row = np.zeros(size)
    row[count : 2 * count] = 1.0
    upper_rows.append(row)
    upper_bounds.append(l1_cap)
    row = np.zeros(size)
    row[:count] = -first_normalized
    upper_rows.append(row)
    upper_bounds.append(0.0)
    for sign in (-1.0, 1.0):
        row = np.zeros(size)
        row[:count] = sign * tail_normalized
        row[rho_index] = -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)

    primal_rows = []
    for k in range(1, primal_cancellations + 1):
        coefficients = np.exp(-8.0 * math.pi * t_values * k)
        primal_rows.append(coefficients)
        row = np.zeros(size)
        row[:count] = coefficients
        upper_rows.append(row)
        upper_bounds.append(leakage)
        upper_rows.append(-row)
        upper_bounds.append(leakage)

    equality = np.zeros((1, size))
    equality[0, :count] = 1.0
    result = linprog(
        objective,
        A_ub=np.asarray(upper_rows),
        b_ub=np.asarray(upper_bounds),
        A_eq=equality,
        b_eq=np.asarray([1.0]),
        bounds=(
            [(None, None)] * count
            + [(0.0, None)] * count
            + [(0.0, None)]
        ),
        method="highs",
    )
    if not result.success or result.x is None:
        return {
            "status": str(result.message),
            "success": False,
            "solver": "scipy-highs",
        }
    weights = np.asarray(result.x[:count], dtype=np.float64)
    first_value = float(first_scalars @ weights)
    tail_value = float(tail_values @ weights)
    return {
        "status": "optimal",
        "success": True,
        "solver": "scipy-highs",
        "coercive_margin": first_value - abs(tail_value),
        "weights": weights.tolist(),
        "weight_l1": float(np.abs(weights).sum()),
        "primal_tail_coefficients": [
            float(row @ weights) for row in primal_rows
        ],
        "first_dual_scalar": first_value,
        "known_higher_dual_scalar": tail_value,
        "known_total_scalar": first_value + tail_value,
        "first_shell_coercive_against_known_scalar_tail": bool(
            first_value > abs(tail_value)
        ),
    }


def simultaneous_channels(
    operators: list[np.ndarray],
) -> tuple[float, np.ndarray, float]:
    """Certify commuting shell operators and return normalized eigenchannels."""
    scale = max(float(np.linalg.norm(M, ord=2)) for M in operators)
    if not math.isfinite(scale) or scale <= 0.0:
        raise ArithmeticError("invalid transformed-operator scale")
    normalized = [M / scale for M in operators]
    reference = sum((math.sqrt(2.0) + j) * M for j, M in enumerate(normalized))
    _values, basis = np.linalg.eigh(reference)
    columns = []
    error = 0.0
    for operator in normalized:
        transformed = basis.T @ operator @ basis
        columns.append(np.diag(transformed))
        off_diagonal = transformed - np.diag(np.diag(transformed))
        error = max(error, float(np.max(np.abs(off_diagonal))))
    if error > 2e-7:
        raise ArithmeticError(
            f"dual shell operators failed simultaneous diagonalization: {error}"
        )
    return scale, np.column_stack(columns), error


def solve_window_sdp(
    p: int,
    t_values: np.ndarray,
    low_operators: list[np.ndarray],
    first_scalars: np.ndarray,
    primal_cancellations: int,
    leakage: float,
    l1_cap: float,
    cancel_first_dual: bool = False,
) -> dict:
    from scipy.optimize import linprog

    count = len(t_values)
    operator_scale, channel_matrix, simultaneous_error = simultaneous_channels(
        low_operators
    )
    # Variables are (a_1,...,a_J, u_1,...,u_J, gamma), with u_j>=|a_j|.
    size = 2 * count + 1
    gamma_index = 2 * count
    objective = np.zeros(size)
    objective[gamma_index] = -1.0
    upper_rows: list[np.ndarray] = []
    upper_bounds: list[float] = []
    for j in range(count):
        row = np.zeros(size)
        row[j], row[count + j] = 1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
        row = np.zeros(size)
        row[j], row[count + j] = -1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
    row = np.zeros(size)
    row[count : 2 * count] = 1.0
    upper_rows.append(row)
    upper_bounds.append(l1_cap)

    primal_rows = []
    for k in range(1, primal_cancellations + 1):
        row = np.exp(-8.0 * math.pi * t_values * k)
        primal_rows.append(row)
        full = np.zeros(size)
        full[:count] = row
        upper_rows.append(full)
        upper_bounds.append(leakage)
        upper_rows.append(-full)
        upper_bounds.append(leakage)

    for channel in channel_matrix:
        row = np.zeros(size)
        row[:count] = -channel
        row[gamma_index] = 1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)

    equality_rows = []
    equality_bounds = []
    row = np.zeros(size)
    row[:count] = 1.0
    equality_rows.append(row)
    equality_bounds.append(1.0)
    if cancel_first_dual:
        row = np.zeros(size)
        row[:count] = first_scalars / np.max(np.abs(first_scalars))
        equality_rows.append(row)
        equality_bounds.append(0.0)

    result = linprog(
        objective,
        A_ub=np.asarray(upper_rows),
        b_ub=np.asarray(upper_bounds),
        A_eq=np.asarray(equality_rows),
        b_eq=np.asarray(equality_bounds),
        bounds=(
            [(None, None)] * count
            + [(0.0, None)] * count
            + [(None, None)]
        ),
        method="highs",
    )
    if not result.success or result.x is None:
        return {
            "status": str(result.message),
            "success": False,
            "cancel_first_dual": cancel_first_dual,
            "solver": "scipy-highs",
            "simultaneous_diagonalization_error": simultaneous_error,
        }
    weights = np.asarray(result.x[:count], dtype=np.float64)
    matrix = sum(weights[j] * low_operators[j] for j in range(count))
    first = float(first_scalars @ weights)
    tail = matrix - first * np.eye(matrix.shape[0])
    eig = np.linalg.eigvalsh(matrix)
    tail_eig = np.linalg.eigvalsh(tail)
    primal_values = [float(row @ weights) for row in primal_rows]
    return {
        "status": "optimal",
        "success": True,
        "solver": "scipy-highs",
        "simultaneous_diagonalization_error": simultaneous_error,
        "cancel_first_dual": cancel_first_dual,
        "gamma": float(result.x[gamma_index]) * operator_scale,
        "weights": weights.tolist(),
        "weight_l1": float(np.abs(weights).sum()),
        "primal_tail_coefficients": primal_values,
        "first_dual_scalar": first,
        "known_higher_dual_min_eigenvalue": float(tail_eig[0]),
        "known_higher_dual_max_abs_eigenvalue": float(
            np.max(np.abs(tail_eig))
        ),
        "known_total_min_eigenvalue": float(eig[0]),
        "known_total_max_eigenvalue": float(eig[-1]),
        "first_shell_coercive_against_known_tail": bool(
            first > np.max(np.abs(tail_eig))
        ),
    }


def solve_coercive_window(
    p: int,
    t_values: np.ndarray,
    low_operators: list[np.ndarray],
    first_scalars: np.ndarray,
    primal_cancellations: int,
    leakage: float,
    l1_cap: float,
) -> dict:
    """Maximize first-shell coefficient minus known-tail operator norm."""
    from scipy.optimize import linprog

    count = len(t_values)
    operator_scale, channels, simultaneous_error = simultaneous_channels(
        low_operators
    )
    first_channels = first_scalars / operator_scale
    tail_channels = channels - first_channels[None, :]
    # Variables are (a_1,...,a_J, u_1,...,u_J, rho), rho>=0.
    size = 2 * count + 1
    rho_index = 2 * count
    objective = np.zeros(size)
    objective[:count] = -first_channels
    objective[rho_index] = 1.0
    upper_rows: list[np.ndarray] = []
    upper_bounds: list[float] = []
    for j in range(count):
        row = np.zeros(size)
        row[j], row[count + j] = 1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
        row = np.zeros(size)
        row[j], row[count + j] = -1.0, -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
    row = np.zeros(size)
    row[count : 2 * count] = 1.0
    upper_rows.append(row)
    upper_bounds.append(l1_cap)
    row = np.zeros(size)
    row[:count] = -first_channels
    upper_rows.append(row)
    upper_bounds.append(0.0)

    primal_rows = []
    for k in range(1, primal_cancellations + 1):
        row = np.exp(-8.0 * math.pi * t_values * k)
        primal_rows.append(row)
        full = np.zeros(size)
        full[:count] = row
        upper_rows.append(full)
        upper_bounds.append(leakage)
        upper_rows.append(-full)
        upper_bounds.append(leakage)
    for channel in tail_channels:
        row = np.zeros(size)
        row[:count] = channel
        row[rho_index] = -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)
        row = np.zeros(size)
        row[:count] = -channel
        row[rho_index] = -1.0
        upper_rows.append(row)
        upper_bounds.append(0.0)

    equality = np.zeros((1, size))
    equality[0, :count] = 1.0
    result = linprog(
        objective,
        A_ub=np.asarray(upper_rows),
        b_ub=np.asarray(upper_bounds),
        A_eq=equality,
        b_eq=np.asarray([1.0]),
        bounds=(
            [(None, None)] * count
            + [(0.0, None)] * count
            + [(0.0, None)]
        ),
        method="highs",
    )
    if not result.success or result.x is None:
        return {
            "status": str(result.message),
            "success": False,
            "solver": "scipy-highs",
        }
    weights = np.asarray(result.x[:count], dtype=np.float64)
    matrix = sum(weights[j] * low_operators[j] for j in range(count))
    first_value = float(first_scalars @ weights)
    tail_matrix = matrix - first_value * np.eye(matrix.shape[0])
    tail_eigenvalues = np.linalg.eigvalsh(tail_matrix)
    total_eigenvalues = np.linalg.eigvalsh(matrix)
    tail_norm = float(np.max(np.abs(tail_eigenvalues)))
    margin = first_value - tail_norm
    return {
        "status": "optimal",
        "success": True,
        "solver": "scipy-highs",
        "simultaneous_diagonalization_error": simultaneous_error,
        "coercive_margin": margin,
        "weights": weights.tolist(),
        "weight_l1": float(np.abs(weights).sum()),
        "primal_tail_coefficients": [
            float(row @ weights) for row in primal_rows
        ],
        "first_dual_scalar": first_value,
        "known_higher_dual_operator_norm": tail_norm,
        "known_higher_dual_min_eigenvalue": float(tail_eigenvalues[0]),
        "known_higher_dual_max_eigenvalue": float(tail_eigenvalues[-1]),
        "known_total_min_eigenvalue": float(total_eigenvalues[0]),
        "known_total_max_eigenvalue": float(total_eigenvalues[-1]),
        "first_shell_coercive_against_known_tail": bool(margin > 1e-7),
    }


def audit(
    p: int,
    scaled_bound: int,
    t_values: np.ndarray,
    primal_cancellations: int,
    leakage: float,
    l1_cap: float,
) -> dict:
    C = paley_conference(p)
    W_basis, eigenspace_basis, intrinsic_W_basis = admissible_W_basis(C, p)
    shells = enumerate_dual_shells(p, scaled_bound)
    shell_ops = phased_dual_shell_operators(
        p,
        shells,
        W_basis,
        eigenspace_basis,
        intrinsic_W_basis,
    )
    first = shell_ops[p]
    expected_first = np.eye(len(W_basis)) / (8.0 * (rank_L(p) + 2))
    first_error = float(np.max(np.abs(first - expected_first)))
    if first_error > 2e-9:
        raise ArithmeticError(f"first dual shell identity failed: {first_error}")

    maxplus = maxplus_shell(p)
    primal_first = harmonic_operator(
        maxplus,
        W_basis,
        eigenspace_basis,
        intrinsic_W_basis,
    )
    primal_eigenvalues = np.linalg.eigvalsh(primal_first)

    lattice_minimum = primal_lattice_minimum_shell(p)
    lattice_minimum_operator = harmonic_operator(
        lattice_minimum,
        W_basis,
        eigenspace_basis,
        intrinsic_W_basis,
    )

    exceptional_operators = {-1: lattice_minimum_operator, **shell_ops}
    exceptional = exceptional_channel(
        primal_first,
        exceptional_operators,
        expected_dimension=rank_L(p),
    )
    lattice_minimum_exceptional_scalar = exceptional["dual_shell_scalars"].pop(-1)

    low = [transformed_low_operator(p, float(t), shell_ops) for t in t_values]
    first_scalars = np.asarray(
        [transformed_first_scalar(p, float(t)) for t in t_values]
    )
    exceptional_low = np.asarray(
        [
            transformed_scalar_channel(
                p,
                float(t),
                exceptional["dual_shell_scalars"],
            )
            for t in t_values
        ]
    )
    ordinary = solve_window_sdp(
        p,
        t_values,
        low,
        first_scalars,
        primal_cancellations,
        leakage,
        l1_cap,
        cancel_first_dual=False,
    )
    control = solve_window_sdp(
        p,
        t_values,
        low,
        first_scalars,
        primal_cancellations,
        leakage,
        l1_cap,
        cancel_first_dual=True,
    )
    coercive = solve_coercive_window(
        p,
        t_values,
        low,
        first_scalars,
        primal_cancellations,
        leakage,
        l1_cap,
    )
    exceptional_coercive = solve_scalar_coercive_window(
        t_values,
        exceptional_low,
        first_scalars,
        primal_cancellations,
        leakage,
        l1_cap,
    )
    shell_summary = {
        str(s): {
            "norm": s / (2 * p),
            "signed_vector_count": 2 * len(shells[s]),
            "phase": -1 if s & 1 else 1,
            "operator_min_eigenvalue": float(np.linalg.eigvalsh(shell_ops[s])[0]),
            "operator_max_eigenvalue": float(np.linalg.eigvalsh(shell_ops[s])[-1]),
            "operator_spectrum": spectrum_clusters(shell_ops[s]),
        }
        for s in sorted(shells)
    }
    return {
        "p": p,
        "n": p * p + 1,
        "rank": rank_L(p),
        "W_dimension": len(W_basis),
        "scaled_dual_bound": scaled_bound,
        "dual_shells": shell_summary,
        "first_dual_shell_identity_max_error": first_error,
        "maxplus_count": len(maxplus),
        "actual_first_primal_harmonic_operator": {
            "min_eigenvalue": float(primal_eigenvalues[0]),
            "max_eigenvalue": float(primal_eigenvalues[-1]),
            "positive_semidefinite": bool(primal_eigenvalues[0] >= -1e-8),
            "spectrum": spectrum_clusters(primal_first),
        },
        "exceptional_channel": exceptional,
        "primal_lattice_minimum_shell": {
            "norm": p + 1,
            "signed_vector_count": len(lattice_minimum),
            "exceptional_scalar": lattice_minimum_exceptional_scalar,
            "operator_spectrum": spectrum_clusters(lattice_minimum_operator),
        },
        "t_values": t_values.tolist(),
        "primal_cancellations": primal_cancellations,
        "leakage_cap": leakage,
        "weight_l1_cap": l1_cap,
        "window": ordinary,
        "first_shell_coercive_window": coercive,
        "exceptional_channel_coercive_window": exceptional_coercive,
        "forced_first_dual_cancellation_control": control,
        "proved_R1": False,
        "finite_shell_diagnostic_only": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=5)
    parser.add_argument("--scaled-bound", type=int, default=20)
    parser.add_argument(
        "--t-values",
        default="0.04,0.05,0.06,0.07,0.085,0.10,0.12,0.15",
    )
    parser.add_argument("--primal-cancellations", type=int, default=3)
    parser.add_argument("--leakage", type=float, default=1e-3)
    parser.add_argument("--l1-cap", type=float, default=50.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    t_values = np.asarray([float(x) for x in args.t_values.split(",")])
    result = audit(
        args.p,
        args.scaled_bound,
        t_values,
        args.primal_cancellations,
        args.leakage,
        args.l1_cap,
    )
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
