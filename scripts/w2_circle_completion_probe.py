#!/usr/bin/env python3
"""Probe the local sign-completion lemma behind the GQR-circle W2 route.

Let C be the Paley conference matrix of order p^2+1 and let S be one
nonsquare F_p-subline through infinity.  There is a signed incidence vector
v_S, supported on S, with C v_S = -p v_S.  A Max- sign vector y can be
flipped on S and remain Max- exactly when y|_S = +/- v_S|_S.

For all such completions in a stored Max+ ensemble (transported to Max-),
this script counts, for every pair outside S, how often

    C_ij y_i y_j = -1.

If every count is positive, PSL transitivity shows that every nonsquare
circle disjoint from a fixed edge is realised by a difference inside that
edge's U-slice.  Pairs inside S are checked separately and should always
have the U sign.

The scan is chunked and reads only p+1 columns on its first pass.  It is
intended for NUKA's local p=11 ensemble as well as the small local caches.
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
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import square_line_matrix  # noqa: E402
from e1_gmin_m4_prop15603 import direction_line_matrix  # noqa: E402


DEFAULT_CACHES = {
    5: Path("/mnt/storage/e1work/qvar_nuka_caches/maxplus_p5.npy"),
    7: Path("/mnt/storage/e1work/qvar_nuka_caches/maxplus_p7.npy"),
    11: Path("/home/nick/e1work/maxplus_p11/maxplus_p11_eps1.npy"),
}


def _order(x: int, q: int, mul) -> int:
    y, order = x, 1
    while y != 1:
        y = mul(y, x)
        order += 1
        if order > q:
            return 0
    return order


def anti_transport(p: int, C: np.ndarray):
    """Return y=d*A[pi] transporting current-label Max+ rows to Max-."""
    q, mul, _add, _chi, _frob, _norm, _ia, _ib = field_ctx(p)
    gen = next(x for x in range(2, q) if _order(x, q, mul) == q - 1)
    pi = np.zeros(q + 1, dtype=np.int64)
    pi[0] = 0
    for x in range(q):
        pi[1 + x] = 1 + mul(x, gen)
    d = np.zeros(q + 1, dtype=np.int8)
    d[0] = 1
    d[1:] = -C[pi[0], pi[1:]] * C[0, 1:]
    return pi, d


def signed_sparse_eigenvector(p: int, C: np.ndarray):
    """One nonsquare circle and its exact {0,+/-1} -p eigenvector."""
    circle = direction_line_matrix(p, square=False)[0].astype(np.uint8)
    supp = np.flatnonzero(circle)
    M = C[:, supp].astype(np.float64).copy()
    M[supp, np.arange(len(supp))] += p
    _u, singular, vh = np.linalg.svd(M, full_matrices=False)
    raw = vh[-1]
    signs = np.sign(raw).astype(np.int8)
    if np.any(signs == 0):
        raise RuntimeError("zero coordinate in sparse eigenvector")
    v = np.zeros(C.shape[0], dtype=np.int8)
    v[supp] = signs
    if not np.array_equal(C @ v.astype(np.int64), -p * v.astype(np.int64)):
        raise RuntimeError(
            f"failed exact sparse eigenvector check; smallest singular={singular[-1]}"
        )
    return circle, v


def completion_indices(
    A: np.ndarray,
    pi: np.ndarray,
    d: np.ndarray,
    v: np.ndarray,
    chunk: int,
):
    supp = np.flatnonzero(v)
    source = pi[supp]
    target = v[supp]
    hits = []
    for lo in range(0, len(A), chunk):
        hi = min(lo + chunk, len(A))
        y = A[lo:hi, source].astype(np.int8, copy=False) * d[supp]
        keep = np.all(y == target, axis=1) | np.all(y == -target, axis=1)
        if np.any(keep):
            hits.append(np.flatnonzero(keep).astype(np.int64) + lo)
        if lo == 0 or (lo // chunk) % 8 == 0:
            found = sum(len(x) for x in hits)
            print(f"  scan {hi}/{len(A)} completions={found}", flush=True)
    if not hits:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(hits)


def pair_counts(X: np.ndarray, C: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Upper-triangle counts of C_ij X_i X_j = -1 on coords."""
    Z = X[:, coords].astype(np.int32)
    gram = Z.T @ Z
    signs = C[np.ix_(coords, coords)].astype(np.int64)
    counts = (len(X) - signs * gram.astype(np.int64)) // 2
    iu = np.triu_indices(len(coords), 1)
    return counts[iu]


def run(p: int, cache: Path, chunk: int) -> dict:
    t0 = time.time()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    A = np.load(cache, mmap_mode="r")
    if A.shape[1] != p * p + 1:
        raise ValueError(f"cache shape {A.shape} does not match p={p}")
    pi, d = anti_transport(p, C)
    circle, v = signed_sparse_eigenvector(p, C)
    supp = np.flatnonzero(circle)
    outside = np.flatnonzero(circle == 0)
    print(
        f"p={p} cache={cache} rows={len(A)} circle={len(supp)} outside={len(outside)}",
        flush=True,
    )
    idx = completion_indices(A, pi, d, v, chunk)
    if len(idx) == 0:
        raise RuntimeError("no sign completions found")
    X = A[idx][:, pi].astype(np.int8) * d
    if not np.all((X[:, supp] == v[supp]).all(axis=1) | (X[:, supp] == -v[supp]).all(axis=1)):
        raise RuntimeError("completion reload/transport mismatch")
    outside_counts = pair_counts(X, C, outside)
    inside_counts = pair_counts(X, C, supp)

    # Independent code-side rank check for p small enough to enumerate all
    # sublines with the existing exact helper.
    circle_rank = None
    n_eligible = None
    if p <= 11:
        from e1_gmin_m4_prop15406 import gf2_rref
        from walsh_subline_dual import all_sublines

        all_circles = all_sublines(p)
        square = square_line_matrix(p)
        parity = (
            all_circles.astype(np.int32) @ square.astype(np.int32).T
        ) & 1
        nonsquare = all_circles[parity.max(axis=1) == 0]
        eligible = nonsquare[(nonsquare[:, 0] ^ nonsquare[:, 1]) == 0]
        circle_rank = int(gf2_rref(eligible.copy())[2])
        n_eligible = int(len(eligible))

    def distribution(values: np.ndarray) -> dict[str, int]:
        keys, counts = np.unique(values, return_counts=True)
        return {str(int(k)): int(c) for k, c in zip(keys, counts)}

    rec = {
        "p": p,
        "cache": str(cache),
        "cache_rows": int(len(A)),
        "circle_support": supp.tolist(),
        "circle_weight": int(len(supp)),
        "sparse_eigenvector_exact": True,
        "n_completions_in_cache": int(len(X)),
        "inside_pair_count_distribution": distribution(inside_counts),
        "outside_pair_count_distribution": distribution(outside_counts),
        "min_outside_U_completions": int(outside_counts.min()),
        "all_inside_pairs_always_U": bool(np.all(inside_counts == len(X))),
        "all_outside_pairs_have_U_completion": bool(np.all(outside_counts > 0)),
        "eligible_circle_rank": circle_rank,
        "eligible_circle_count": n_eligible,
        "target_slice_rank": (p * p - 1) // 2,
        "seconds": round(time.time() - t0, 3),
    }
    print(json.dumps(rec, indent=2), flush=True)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int, choices=sorted(DEFAULT_CACHES))
    ap.add_argument("--cache", type=Path)
    ap.add_argument("--chunk", type=int, default=1_000_000)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    cache = args.cache or DEFAULT_CACHES[args.p]
    rec = run(args.p, cache, args.chunk)
    if args.output:
        args.output.write_text(json.dumps(rec, indent=2) + "\n")
        print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
