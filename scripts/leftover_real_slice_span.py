#!/usr/bin/env python3
"""Real edge-feature span on a fixed Max-minus pair slice.

This is a diagnostic for the minus-slice branch of residual (ii).  Put

    f_ab(y) = C_ab y_a y_b,   U = {y in Max-minus : f_01(y)=-1}.

If a graph G has S_G=-2 on all of U, its edge-incidence vector lies in
the real annihilator of {f(y)-f(y0): y in U}.  The script computes that
annihilator numerically, compares it with the always-constant star space
plus the distinguished edge, and checks the p=5 leftover witness.

The singular-value gap is only a discovery diagnostic; any theorem based
on it must subsequently receive an exact algebraic rank proof.
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

from e1_gmin_m4_prop15406 import WITNESS  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


DEFAULT_MINUS = {
    5: Path("/tmp/maxminus_p5.npy"),
    7: Path("/tmp/maxminus_p7.npy"),
}


def edge_data(n: int) -> tuple[list[tuple[int, int]], np.ndarray, np.ndarray]:
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    left = np.fromiter((a for a, _b in edges), dtype=np.int32)
    right = np.fromiter((b for _a, b in edges), dtype=np.int32)
    return edges, left, right


def known_constant_space(
    n: int, edges: list[tuple[int, int]], distinguished: tuple[int, int]
) -> np.ndarray:
    """Columns: n vertex stars and the distinguished edge."""
    idx = {edge: j for j, edge in enumerate(edges)}
    K = np.zeros((len(edges), n + 1), dtype=np.float64)
    for j, (a, b) in enumerate(edges):
        K[j, a] = 1.0
        K[j, b] = 1.0
    K[idx[distinguished], n] = 1.0
    return K


def wedge_slice_generators(
    C: np.ndarray, edges: list[tuple[int, int]]
) -> np.ndarray:
    """One constant star, f_01, and f_0k+Delta_k f_1k for k>=2."""
    n = C.shape[0]
    idx = {edge: j for j, edge in enumerate(edges)}
    G = np.zeros((len(edges), n), dtype=np.float64)
    # The vertex-0 star scores identically -p on Max-minus.
    for k in range(1, n):
        G[idx[(0, k)], 0] = 1.0
    G[idx[(0, 1)], 1] = 1.0
    for col, k in enumerate(range(2, n), start=2):
        delta = int(C[0, 1]) * int(C[0, k]) * int(C[1, k])
        G[idx[(0, k)], col] = 1.0
        G[idx[(1, k)], col] = float(delta)
    return G


def run_prime(p: int, minus_path: Path) -> dict:
    t0 = time.time()
    Y = np.sign(np.load(minus_path, mmap_mode="r")).astype(np.int8)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    n = C.shape[0]
    edges, left, right = edge_data(n)
    signs = C[left, right].astype(np.int8)
    F = (Y[:, left] * Y[:, right] * signs).astype(np.int8)
    e_index = edges.index((0, 1))
    U = F[:, e_index] == -1
    FU = F[U]
    D = ((FU[1:].astype(np.int16) - FU[0].astype(np.int16)) // 2).astype(
        np.float64
    )

    # The right singular vectors furnish the annihilator directly.  Full
    # matrices are unnecessary and prohibitively large for p=7.
    singular = np.linalg.svd(D, full_matrices=False, compute_uv=False)
    scale = singular[0] if singular.size else 1.0
    tol = max(D.shape) * np.finfo(np.float64).eps * scale
    rank = int(np.count_nonzero(singular > tol))
    nullity = int(D.shape[1] - rank)

    full_singular = np.linalg.svd(F.astype(np.float64), full_matrices=False, compute_uv=False)
    full_scale = full_singular[0] if full_singular.size else 1.0
    full_tol = max(F.shape) * np.finfo(np.float64).eps * full_scale
    full_rank = int(np.count_nonzero(full_singular > full_tol))
    slice_function_dim = full_rank - rank

    K = known_constant_space(n, edges, (0, 1))
    known_rank = int(np.linalg.matrix_rank(K))
    known_residual = float(np.max(np.abs(D @ K)))
    G = wedge_slice_generators(C, edges)
    FG = F.astype(np.float64) @ G
    generator_function_rank = int(np.linalg.matrix_rank(FG))
    generator_u_spread = float(np.max(np.abs(D @ G)))

    row = {
        "p": p,
        "n": n,
        "ensemble_rows": int(Y.shape[0]),
        "U_rows": int(U.sum()),
        "edges": len(edges),
        "variation_rank": rank,
        "annihilator_dim": nullity,
        "full_score_function_rank": full_rank,
        "constant_on_U_function_dim": slice_function_dim,
        "expected_span_1_fe_dim": 2,
        "constant_on_U_is_only_span_1_fe": slice_function_dim == 2,
        "wedge_slice_generator_function_rank": generator_function_rank,
        "wedge_slice_generators_complete": generator_function_rank
        == slice_function_dim,
        "wedge_slice_generator_U_variation": generator_u_spread,
        "known_star_plus_edge_dim": known_rank,
        "extra_annihilator_dim": nullity - known_rank,
        "known_space_annihilated_max_error": known_residual,
        "sv_largest": float(singular[0]) if singular.size else 0.0,
        "sv_smallest_positive": float(singular[rank - 1]) if rank else 0.0,
        "sv_first_zero": float(singular[rank]) if rank < singular.size else None,
        "numeric_tolerance": float(tol),
        "full_numeric_tolerance": float(full_tol),
        "seconds": round(time.time() - t0, 3),
    }

    if p == 5:
        idx = {edge: j for j, edge in enumerate(edges)}
        x = np.zeros(len(edges), dtype=np.float64)
        for edge in WITNESS:
            x[idx[edge]] = 1.0
        scores_u = FU.astype(np.int16) @ x.astype(np.int16)
        # Distance from x to the star+edge space detects whether this known
        # leftover is already explained by the obvious constants.
        coeff, *_ = np.linalg.lstsq(K, x, rcond=None)
        residual = x - K @ coeff
        row["p5_witness"] = {
            "size": int(x.sum()),
            "U_score_values": sorted({int(v) for v in scores_u.tolist()}),
            "in_annihilator": bool(np.max(np.abs(D @ x)) == 0),
            "distance_to_star_plus_edge": float(np.linalg.norm(residual)),
        }
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", type=int, nargs="+", default=[5, 7])
    ap.add_argument("--minus", type=Path, help="override cache for one prime")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    rows = {}
    for p in args.primes:
        path = args.minus if args.minus is not None else DEFAULT_MINUS[p]
        row = run_prime(p, path)
        rows[str(p)] = row
        print(json.dumps(row, indent=2), flush=True)
    out = {
        "title": "real edge-feature annihilator on a Max-minus pair slice",
        "diagnostic_only": True,
        "rows": rows,
    }
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
