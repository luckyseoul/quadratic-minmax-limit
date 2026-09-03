#!/usr/bin/env python3
"""Exact small-order calculator for the original MO two-half geometry.

For a symmetric signing A with M = Phi(A), and a skew signing R, define

    U(x,y) = |Q_A(x) + Q_A(y)|,
    W(x,y) = |x^T R y|,
    D(x,y) = 2M - U(x,y).

The equal-endpoint Hadamard doubling target is exactly

    max_{x,y} (W-D) <= 2(sqrt(2)-1)M + error.

This script verifies stored optimal A and R instances, optionally re-solves
the finite orientation minimax problem, writes a JSON record, and plots the
upper envelope W_max(D).  It deliberately stops at n=8: the calculation is a
geometric probe and a finite certificate, not evidence of asymptotics.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = ROOT / "evidence" / "original_mo_two_half_geometry.json"
DEFAULT_PLOT = ROOT / "evidence" / "original_mo_two_half_geometry.png"


def _negative_cycle_five() -> np.ndarray:
    matrix = np.ones((5, 5), dtype=np.int64) - np.eye(5, dtype=np.int64)
    for i in range(5):
        j = (i + 1) % 5
        matrix[i, j] = matrix[j, i] = -1
    return matrix


INSTANCES: dict[int, dict[str, Any]] = {
    5: {
        "m": 4,
        "A": _negative_cycle_five().tolist(),
        "r_bits": "0000000000",
        "expected_B": 16,
        "expected_excess": 8,
        "expected_optimal_count": 780,
        "method": "complete enumeration of all 2^10 skew signings",
    },
    6: {
        "m": 5,
        "A": [
            [0, 1, 1, 1, 1, 1],
            [1, 0, 1, -1, -1, 1],
            [1, 1, 0, 1, -1, -1],
            [1, -1, 1, 0, 1, -1],
            [1, -1, -1, 1, 0, 1],
            [1, 1, -1, -1, 1, 0],
        ],
        "r_bits": "100111001110000",
        "expected_B": 18,
        "expected_excess": 8,
        "expected_optimal_count": 1760,
        "method": "complete enumeration of all 2^15 skew signings",
    },
    7: {
        "m": 9,
        "A": [
            [0, 1, 1, 1, 1, 1, 1],
            [1, 0, -1, 1, 1, 1, 1],
            [1, -1, 0, 1, 1, -1, -1],
            [1, 1, 1, 0, -1, -1, 1],
            [1, 1, 1, -1, 0, 1, -1],
            [1, 1, -1, -1, 1, 0, -1],
            [1, 1, -1, 1, -1, -1, 0],
        ],
        "r_bits": "001000001011100100100",
        "expected_B": 22,
        "expected_excess": 4,
        "method": "complete CP-SAT minimax model; OPTIMAL required",
    },
    8: {
        "m": 10,
        "A": [
            [0, 1, 1, 1, 1, 1, 1, 1],
            [1, 0, -1, -1, -1, 1, 1, -1],
            [1, -1, 0, -1, -1, -1, 1, 1],
            [1, -1, -1, 0, 1, 1, -1, 1],
            [1, -1, -1, 1, 0, -1, 1, -1],
            [1, 1, -1, 1, -1, 0, -1, 1],
            [1, 1, 1, -1, 1, -1, 0, -1],
            [1, -1, 1, 1, -1, 1, -1, 0],
        ],
        "r_bits": "0101001000000110001101011101",
        "expected_B": 28,
        "expected_excess": 8,
        "method": "complete CP-SAT minimax model; OPTIMAL required",
    },
}


def boolean_states(n: int, fix_first: bool = False) -> np.ndarray:
    """Return Boolean rows safely, avoiding unsigned subtraction."""
    masks = np.arange(1 << n, dtype=np.uint64)[:, None]
    shifts = np.arange(n, dtype=np.uint64)[None, :]
    bits = ((masks >> shifts) & np.uint64(1)).astype(np.int64)
    states = 2 * bits - 1
    return states[states[:, 0] == 1] if fix_first else states


def edge_pairs(n: int) -> list[tuple[int, int]]:
    return list(itertools.combinations(range(n), 2))


def rademacher_abs_mean(k: int) -> float:
    """Return E|eps_1+...+eps_k| from its exact binomial formula."""
    if k < 1:
        raise ValueError("k must be positive")
    return k * math.comb(k - 1, (k - 1) // 2) / (1 << (k - 1))


def sharp_influence_constant(n: int, m: int) -> float:
    """K_n = max_A Inf_1(Q_A/Phi(A)) = n*mu_(n-1)/m_n."""
    if n < 2 or m < 1:
        raise ValueError("expected n>=2 and m>=1")
    return n * rademacher_abs_mean(n - 1) / m


def skew_from_bits(n: int, bits: str) -> np.ndarray:
    edges = edge_pairs(n)
    if len(bits) != len(edges) or set(bits) - {"0", "1"}:
        raise ValueError(f"expected {len(edges)} binary digits for n={n}")
    matrix = np.zeros((n, n), dtype=np.int64)
    for (i, j), bit in zip(edges, bits, strict=True):
        matrix[i, j] = 1 if bit == "1" else -1
        matrix[j, i] = -matrix[i, j]
    return matrix


def bits_from_skew(matrix: np.ndarray) -> str:
    return "".join("1" if matrix[i, j] == 1 else "0" for i, j in edge_pairs(len(matrix)))


def quadratic_values(matrix: np.ndarray, states: np.ndarray) -> np.ndarray:
    return np.einsum("bi,ij,bj->b", states, matrix, states, dtype=np.int64) // 2


def directed_halfcut_norms(A: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Return Phi of every directed half-cut neighbor, indexed by cut mask.

    With S=A*R viewed as a tournament and U the mask, the neighbor reverses
    exactly the cut arcs directed from U to its complement.  Its matrix is
    (A + D A D + D R - R D)/2, where D is -1 on U and +1 off U.
    """
    n = len(A)
    states = boolean_states(n, fix_first=True)
    norms = np.empty(1 << n, dtype=np.int64)
    for mask in range(1 << n):
        signs = np.asarray(
            [-1 if mask & (1 << i) else 1 for i in range(n)],
            dtype=np.int64,
        )
        switched = A * signs[:, None] * signs[None, :]
        commutator = signs[:, None] * R - R * signs[None, :]
        neighbor = (A + switched + commutator) // 2
        norms[mask] = int(np.max(np.abs(quadratic_values(neighbor, states))))
    return norms


def directed_cut_halves(
    A: np.ndarray, R: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[int, int]:
    """Return the signed outward/inward A-energies F,G for U={x!=y}."""
    n = len(A)
    in_cut_side = x != y
    tournament = A * R
    outward = 0
    inward = 0
    for i, j in edge_pairs(n):
        if in_cut_side[i] == in_cut_side[j]:
            continue
        if in_cut_side[i]:
            u, v = i, j
        else:
            u, v = j, i
        contribution = int(A[u, v] * y[u] * y[v])
        if tournament[u, v] == 1:
            outward += contribution
        else:
            inward += contribution
    return outward, inward


def validate_signings(A: np.ndarray, R: np.ndarray) -> None:
    n = len(A)
    if A.shape != (n, n) or R.shape != (n, n):
        raise AssertionError("non-square input")
    if not np.array_equal(A, A.T) or np.any(np.diag(A)):
        raise AssertionError("A is not symmetric with zero diagonal")
    if not np.array_equal(R, -R.T) or np.any(np.diag(R)):
        raise AssertionError("R is not skew with zero diagonal")
    mask = ~np.eye(n, dtype=bool)
    if not np.all(np.abs(A[mask]) == 1) or not np.all(np.abs(R[mask]) == 1):
        raise AssertionError("off-diagonal entries are not signs")


def analyze_orientation(A: np.ndarray, M: int, R: np.ndarray) -> dict[str, Any]:
    validate_signings(A, R)
    states = boolean_states(len(A), fix_first=True)
    q = quadratic_values(A, states)
    phi = int(np.max(np.abs(q)))
    U = np.abs(q[:, None] + q[None, :])
    W = np.abs(states @ R @ states.T)
    D = 2 * M - U
    score = U + W
    B = int(np.max(score))
    excess = int(np.max(W - D))
    if B - 2 * M != excess:
        raise AssertionError("B-2M and max(W-D) disagree")
    halfcut_norms = directed_halfcut_norms(A, R)
    halfcut_max = int(np.max(halfcut_norms))
    if B != 2 * halfcut_max:
        raise AssertionError("fourth-phase and directed half-cut norms disagree")

    envelope: dict[int, int] = {}
    point_counts: dict[int, int] = {}
    for d in sorted(int(value) for value in np.unique(D)):
        mask = D == d
        envelope[d] = int(np.max(W[mask]))
        point_counts[d] = int(np.count_nonzero(mask))

    wi, wj = np.argwhere(score == B)[0]
    x = states[int(wi)]
    y = states[int(wj)]
    u = (x + y) // 2
    v = (x - y) // 2
    I = int((u @ A @ u + v @ A @ v) // 2)
    X = int(u @ A @ v)
    C = int(-(u @ R @ v))
    qx, qy = int(q[int(wi)]), int(q[int(wj)])
    cross = int(x @ R @ y)
    if (qx, qy) != (I + X, I - X):
        raise AssertionError("two-half internal identity failed")
    if cross != 2 * C:
        raise AssertionError("two-half skew identity failed")
    if abs(I) + abs(X) > M:
        raise AssertionError("endpoint diamond failed")
    epsilon = M - max(abs(qx), abs(qy))
    separation = abs(qx - qy)
    witness_D = int(2 * M - abs(qx + qy))
    if witness_D != separation + 2 * epsilon:
        raise AssertionError("endpoint slack decomposition failed")
    outward, inward = directed_cut_halves(A, R, x, y)
    if cross != 2 * (inward - outward):
        raise AssertionError("directed cut-half identity failed")
    cancellation = 4 * min(abs(outward), abs(inward))
    if outward * inward >= 0:
        cancellation *= -1
    if abs(cross) - separation != cancellation:
        raise AssertionError("directed cut cancellation identity failed")

    target = 2 * (math.sqrt(2) - 1) * M
    return {
        "n": len(A),
        "M": M,
        "phi_A": phi,
        "B": B,
        "max_W_minus_D": excess,
        "target_excess": target,
        "target_margin": target - excess,
        "zero_error_doubling_passes": excess <= target,
        "r_bits": bits_from_skew(R),
        "directed_halfcut": {
            "max_phi": halfcut_max,
            "identity_twice_max_phi_equals_B": True,
            "maximizing_cut_masks": [
                format(mask, f"0{len(A)}b")
                for mask in np.flatnonzero(halfcut_norms == halfcut_max)
            ],
            "phi_histogram": {
                str(value): int(np.count_nonzero(halfcut_norms == value))
                for value in sorted(int(v) for v in np.unique(halfcut_norms))
            },
        },
        "envelope_Wmax_by_D": {str(d): envelope[d] for d in envelope},
        "pair_count_by_D": {str(d): point_counts[d] for d in envelope},
        "maximizing_witness": {
            "x": x.tolist(),
            "y": y.tolist(),
            "Qx": qx,
            "Qy": qy,
            "xRy": cross,
            "D": witness_D,
            "endpoint_separation": separation,
            "epsilon": epsilon,
            "D_equals_separation_plus_twice_epsilon": True,
            "outward_half_F": outward,
            "inward_half_G": inward,
            "opposite_sign_halves": outward * inward < 0,
            "cross_minus_separation": abs(cross) - separation,
            "u": u.tolist(),
            "v": v.tolist(),
            "I": I,
            "X": X,
            "C": C,
            "H": int(M - abs(I) - abs(X)),
        },
    }


def _orientation_coefficients(states: np.ndarray) -> np.ndarray:
    edges = edge_pairs(states.shape[1])
    rows = []
    for x in states:
        for y in states:
            rows.append([int(x[i] * y[j] - x[j] * y[i]) for i, j in edges])
    return np.asarray(rows, dtype=np.int64)


def solve_exhaustive(A: np.ndarray, M: int, chunk_size: int = 4096) -> dict[str, Any]:
    """Enumerate every R.  Intended only for n<=6."""
    n = len(A)
    edges = edge_pairs(n)
    total = 1 << len(edges)
    states = boolean_states(n, fix_first=True)
    q = quadratic_values(A, states)
    U = np.abs(q[:, None] + q[None, :]).reshape(-1)
    coefficients = _orientation_coefficients(states)
    best_B: int | None = None
    best_bits = ""
    optimal_count = 0
    for start in range(0, total, chunk_size):
        stop = min(total, start + chunk_size)
        masks = np.arange(start, stop, dtype=np.uint64)[:, None]
        shifts = np.arange(len(edges), dtype=np.uint64)[None, :]
        signs = 2 * (((masks >> shifts) & np.uint64(1)).astype(np.int64)) - 1
        cross = coefficients @ signs.T
        scores = np.max(U[:, None] + np.abs(cross), axis=0)
        chunk_best = int(np.min(scores))
        if best_B is None or chunk_best < best_B:
            best_B = chunk_best
            optimal_count = int(np.count_nonzero(scores == chunk_best))
            index = int(np.flatnonzero(scores == chunk_best)[0])
            best_bits = "".join("1" if value == 1 else "0" for value in signs[index])
        elif chunk_best == best_B:
            optimal_count += int(np.count_nonzero(scores == chunk_best))
    assert best_B is not None
    return {
        "status": "EXHAUSTIVE",
        "objective_B": best_B,
        "best_r_bits": best_bits,
        "optimal_orientation_count": optimal_count,
        "orientations_checked": total,
        "max_W_minus_D": best_B - 2 * M,
    }


def solve_cpsat(A: np.ndarray, M: int, hint_bits: str, workers: int, max_time: float) -> dict[str, Any]:
    """Solve the complete finite minimax model; OPTIMAL is mandatory."""
    from ortools.sat.python import cp_model

    n = len(A)
    edges = edge_pairs(n)
    states = boolean_states(n, fix_first=True)
    q = quadratic_values(A, states)
    model = cp_model.CpModel()
    bits = [model.new_bool_var(f"r_{i}_{j}") for i, j in edges]
    # R and -R have the same objective.
    model.add(bits[0] == 0)
    B = model.new_int_var(0, 2 * M + n * n, "B")
    for xi, x in enumerate(states):
        for yi, y in enumerate(states):
            U = abs(int(q[xi] + q[yi]))
            coeffs = [int(x[i] * y[j] - x[j] * y[i]) for i, j in edges]
            expr = sum(c * (2 * bit - 1) for c, bit in zip(coeffs, bits, strict=True))
            model.add(B >= U + expr)
            model.add(B >= U - expr)
    model.minimize(B)
    for bit, value in zip(bits, hint_bits, strict=True):
        model.add_hint(bit, int(value))

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if max_time > 0:
        solver.parameters.max_time_in_seconds = max_time
    status = solver.solve(model)
    status_name = solver.status_name(status)
    result = {
        "status": status_name,
        "objective_B": int(round(solver.objective_value)) if status in (cp_model.FEASIBLE, cp_model.OPTIMAL) else None,
        "best_bound_B": int(round(solver.best_objective_bound)),
    }
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"n={n} finite minimax solve was {status_name}, not OPTIMAL: {result}")
    # Keep the persisted certificate independent of parallel search order.
    # Timing, branch/conflict counts, and the particular optimal witness can
    # all vary across otherwise identical CP-SAT replays.  The pinned stored
    # witness is checked separately by analyze_orientation; here only the
    # proved optimum and bound belong in the reproducible record.
    result["max_W_minus_D"] = result["objective_B"] - 2 * M
    return result


def solve_instance(n: int, workers: int, max_time: float) -> dict[str, Any]:
    instance = INSTANCES[n]
    A = np.asarray(instance["A"], dtype=np.int64)
    M = int(instance["m"])
    if n <= 6:
        return solve_exhaustive(A, M)
    return solve_cpsat(A, M, str(instance["r_bits"]), workers, max_time)


def build_record(run_solver: bool, workers: int, max_time: float) -> dict[str, Any]:
    records: dict[str, Any] = {}
    for n, instance in INSTANCES.items():
        A = np.asarray(instance["A"], dtype=np.int64)
        M = int(instance["m"])
        stored = analyze_orientation(A, M, skew_from_bits(n, str(instance["r_bits"])))
        if stored["phi_A"] != M:
            raise AssertionError(f"stored A at n={n} has Phi={stored['phi_A']}, expected {M}")
        if stored["B"] != instance["expected_B"] or stored["max_W_minus_D"] != instance["expected_excess"]:
            raise AssertionError(f"stored n={n} result disagrees with pinned values")
        row: dict[str, Any] = {
            "classification": "exhaustive finite certificate",
            "orientation_optimization_method": instance["method"],
            "A": instance["A"],
            "stored_orientation": stored,
        }
        if run_solver:
            solved = solve_instance(n, workers, max_time)
            if solved["objective_B"] != stored["B"]:
                raise AssertionError(f"n={n} stored R is not optimal")
            expected_count = instance.get("expected_optimal_count")
            if expected_count is not None and solved.get("optimal_orientation_count") != expected_count:
                raise AssertionError(f"n={n} optimal orientation count changed")
            row["independent_optimization"] = solved
        records[str(n)] = row

    exact_table_path = ROOT / "evidence" / "exact_m_table.json"
    exact_table_payload = exact_table_path.read_bytes()
    with exact_table_path.open(encoding="utf-8") as handle:
        exact = json.load(handle)
    influence_sequence = {
        key: {
            "m": int(value["m"]),
            "mu_n_minus_1": rademacher_abs_mean(int(key) - 1),
            "K_n": sharp_influence_constant(int(key), int(value["m"])),
        }
        for key, value in exact.items()
        if key.isdigit()
    }

    return {
        "schema": "original-mo-two-half-geometry-v1",
        "status": "finite geometric evidence; original MO convergence remains open",
        "exact_m_table_sha256": hashlib.sha256(exact_table_payload).hexdigest(),
        "definitions": {
            "Q_A(x)": "sum_{i<j} A_ij x_i x_j",
            "M": "Phi(A)=max_x |Q_A(x)|",
            "U": "|Q_A(x)+Q_A(y)|",
            "W": "|x^T R y|",
            "D": "2M-U",
            "epsilon(x,y)": "M-max(|Q_A(x)|,|Q_A(y)|)",
            "B(A,R)": "max_{x,y}(U+W)",
            "identity": "B(A,R)-2M=max_{x,y}(W-D)",
            "endpoint_slack_identity": "D=|Q_A(x)-Q_A(y)|+2epsilon(x,y)",
            "directed_cancellation_identity": "W-|Q_A(x)-Q_A(y)| is 4min(|F|,|G|) for FG<0 and -4min(|F|,|G|) for FG>=0",
            "directed_halfcut_identity": "B(A,R)/2=max_U Phi(A with outward S=A*R half-cut flipped)",
            "zero_error_target": "max(W-D)<=2(sqrt(2)-1)M",
            "sharp_influence_identity": "m_n=n*mu_(n-1)/K_n, where mu_k=E|sum_(1..k) eps_i|",
        },
        "exact_small_order_influence_sequence": influence_sequence,
        "solver_replayed": run_solver,
        "instances": records,
    }


def plot_record(record: dict[str, Any], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with (ROOT / "evidence" / "exact_m_table.json").open(encoding="utf-8") as handle:
        exact = json.load(handle)

    fig, axes = plt.subplots(3, 2, figsize=(12, 13), constrained_layout=True)
    ns = sorted(int(key) for key in exact if key.isdigit())
    alphas = [float(exact[str(n)]["alpha"]) for n in ns]
    ax = axes[0, 0]
    alpha_line = ax.plot(ns, alphas, "o-", color="#244a7c", linewidth=1.8, label=r"$\alpha_n$")
    ax.axhline(1 / math.pi, color="#6b7280", linestyle="--", label=r"proved lower bound $1/\pi$")
    ax.axhline(0.5, color="#8b5cf6", linestyle=":", label=r"upper benchmark $1/2$")
    for n, alpha in zip(ns, alphas, strict=True):
        ax.annotate(str(exact[str(n)]["m"]), (n, alpha), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=8)
    ax.set(title=r"Exact values and equivalent influence constant", xlabel="n", ylabel=r"$\alpha_n$")
    ax.grid(alpha=0.25)
    influence_ax = ax.twinx()
    influence = [record["exact_small_order_influence_sequence"][str(n)]["K_n"] for n in ns]
    influence_line = influence_ax.plot(
        ns,
        influence,
        "s--",
        color="#b45309",
        linewidth=1.2,
        markersize=4,
        label=r"$K_n=n\mu_{n-1}/m_n$",
    )
    influence_ax.set_ylabel(r"sharp influence constant $K_n$", color="#b45309")
    influence_ax.tick_params(axis="y", labelcolor="#b45309")
    lines = alpha_line + influence_line + ax.lines[1:]
    ax.legend(lines, [line.get_label() for line in lines], fontsize=7, loc="best")

    ax = axes[0, 1]
    calc_ns = sorted(int(key) for key in record["instances"])
    ratios = [record["instances"][str(n)]["stored_orientation"]["max_W_minus_D"] / INSTANCES[n]["m"] for n in calc_ns]
    colors = ["#2c7a7b" if value <= 2 * (math.sqrt(2) - 1) else "#c2413b" for value in ratios]
    ax.bar([str(n) for n in calc_ns], ratios, color=colors)
    ax.axhline(2 * (math.sqrt(2) - 1), color="black", linestyle="--", label=r"target $2(\sqrt{2}-1)$")
    ax.set(title="Optimized vertical excess", xlabel="n", ylabel=r"$\min_R\max(W-D)/M$")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)

    for ax, n in zip(axes[1:].flat, calc_ns, strict=True):
        result = record["instances"][str(n)]["stored_orientation"]
        envelope = {int(d): int(w) for d, w in result["envelope_Wmax_by_D"].items()}
        ds = sorted(envelope)
        ws = [envelope[d] for d in ds]
        M = int(result["M"])
        upper_x = max(ds) if ds else 2 * M
        line_x = np.linspace(0, upper_x, 200)
        line_y = line_x + 2 * (math.sqrt(2) - 1) * M
        ax.plot(ds, ws, "o-", color="#244a7c", linewidth=1.8, label=r"exact $W_{\max}(D)$")
        ax.plot(line_x, line_y, "--", color="#c2413b", label="required diagonal")
        ax.fill_between(line_x, line_y, max(max(ws) + 1, max(line_y) + 1), color="#fecaca", alpha=0.25)
        verdict = "passes" if result["zero_error_doubling_passes"] else "fails"
        ax.set(
            title=f"n={n}, M={M}: g={result['max_W_minus_D']} ({verdict})",
            xlabel=r"internal slack $D=2M-|Q(x)+Q(y)|$",
            ylabel=r"skew interaction $W=|x^TRy|$",
        )
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

    fig.suptitle("Original MO problem: exact two-half diagonal envelopes", fontsize=16)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true", help="independently optimize every stored finite instance")
    parser.add_argument("--workers", type=int, default=8, help="CP-SAT workers for n=7,8")
    parser.add_argument("--max-time", type=float, default=180.0, help="per-instance CP-SAT time limit; 0 means unlimited")
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args()

    record = build_record(args.solve, args.workers, args.max_time)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(record, indent=2, sort_keys=True) + "\n"
    args.json.write_text(payload, encoding="utf-8")
    plot_record(record, args.plot)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    print(json.dumps({
        "json": str(args.json),
        "plot": str(args.plot),
        "json_sha256": digest,
        "solver_replayed": args.solve,
        "B": {n: record["instances"][str(n)]["stored_orientation"]["B"] for n in INSTANCES},
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
