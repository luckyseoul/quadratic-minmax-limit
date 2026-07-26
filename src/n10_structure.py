#!/usr/bin/env python3
"""
Structural gap campaign at n=10 (Paley p=3, Φ=15, m_10=13).

Goals (ranked attack #2 in HANDOFF §5.1):
  1. Maximizer balance of Paley C_10 on every edge (Prop 15.21).
  2. Exact min Φ at Hamming distance k from Paley (undirected edge flips).
  3. Min Hamming distance (over Seidel switching class) from Paley to any
     matrix with Φ ≤ 13 (i.e. any exact optimum).
  4. Spectral fingerprint of found optima: op, δ=tr(A^4)-n(n-1)^2, ρ, r.

All Φ evaluations are exact (n=10 half-cube). Parallel via ProcessPool.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    form_Q,
    halfspace_boolean_vector,
    is_conference_matrix,
    op_norm_sym,
    paley_conference_prime_power,
    phi,
    rho_from_boolean_evec,
)

# BLAS single-thread per process (set before heavy numpy in workers)
for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")

N = 10
M_TARGET = 13  # exact m_10
PHI_PALEY = 15  # (1/2)*n*p with p=3


def _edges(n: int = N) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def _flip_edge(A: np.ndarray, i: int, j: int) -> None:
    A[i, j] *= -1.0
    A[j, i] *= -1.0


def _hamming_edges(A: np.ndarray, B: np.ndarray) -> int:
    """Undirected edge Hamming distance (number of disagreeing upper-triangle signs)."""
    n = A.shape[0]
    d = 0
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] * B[i, j] < 0:
                d += 1
    return d


def switch(A: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Seidel switching: A' = D A D with D=diag(eps), eps in {±1}^n."""
    return (eps[:, None] * A) * eps[None, :]


def all_switchings_of(C: np.ndarray, fix_first: bool = True) -> List[np.ndarray]:
    """
    Full switching class of C (up to global sign of eps if fix_first).
    For n=10: 2^9 = 512 matrices.
    """
    n = C.shape[0]
    out = []
    # eps[0] fixed +1 if fix_first
    nfree = n - 1 if fix_first else n
    for mask in range(1 << nfree):
        eps = np.ones(n, dtype=np.float64)
        for b in range(nfree):
            if mask & (1 << b):
                idx = b + 1 if fix_first else b
                eps[idx] = -1.0
        out.append(switch(C, eps))
    return out


def exact_phi_and_maximizers(A: np.ndarray) -> Tuple[float, List[np.ndarray]]:
    """Exact Φ and all maximizers with x[0]=+1 (half-cube)."""
    n = A.shape[0]
    npat = 1 << (n - 1)
    best = 0.0
    maxes: List[np.ndarray] = []
    for xb in range(npat):
        x = np.ones(n, dtype=np.float64)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = form_Q(A, x)
        aq = abs(q)
        if aq > best + 1e-12:
            best = aq
            maxes = [x.copy()]
        elif abs(aq - best) < 1e-12:
            maxes.append(x.copy())
    return best, maxes


def maximizer_balance_report(C: np.ndarray) -> Dict:
    """Prop 15.21 balance check on every undirected edge of conference C."""
    M, maxes = exact_phi_and_maximizers(C)
    edges = _edges(C.shape[0])
    unbalanced = []
    balanced_count = 0
    for p, r in edges:
        # edge balanced if some maximizer has ε(x)=C_pr x_p x_r = -sign(Q_C(x))
        ok = False
        for x in maxes:
            q = form_Q(C, x)
            s = 1.0 if q >= 0 else -1.0
            eps = C[p, r] * x[p] * x[r]
            if eps == -s or (abs(q) < 1e-12 and eps == -1.0):
                ok = True
                break
        if ok:
            balanced_count += 1
        else:
            unbalanced.append((int(p), int(r)))
    return {
        "phi": M,
        "n_maximizers_halfcube": len(maxes),
        "n_edges": len(edges),
        "n_balanced": balanced_count,
        "n_unbalanced": len(unbalanced),
        "unbalanced_edges": unbalanced,
        "single_flip_local_opt": len(unbalanced) == 0,
    }


def phi_after_flips(C: np.ndarray, flip_idxs: Tuple[int, ...], edges: List[Tuple[int, int]]) -> float:
    A = C.copy()
    for t in flip_idxs:
        i, j = edges[t]
        _flip_edge(A, i, j)
    return phi(A)


def _kflip_shard(args) -> Tuple[int, float, Optional[Tuple[int, ...]]]:
    """Worker: min Φ over a list of k-subsets of edge indices."""
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    C_list, edges, combos = args
    C = np.array(C_list, dtype=np.float64)
    edges_t = [(int(a), int(b)) for a, b in edges]
    best = float("inf")
    best_combo = None
    for combo in combos:
        A = C.copy()
        for t in combo:
            i, j = edges_t[t]
            _flip_edge(A, i, j)
        v = phi(A)
        if v < best:
            best = v
            best_combo = tuple(int(x) for x in combo)
    return len(combos), best, best_combo


def min_phi_at_distance_k(
    C: np.ndarray, k: int, workers: int
) -> Dict:
    """Exact min Φ among all matrices at undirected Hamming distance k from C."""
    edges = _edges(C.shape[0])
    ne = len(edges)
    if k == 0:
        return {"k": 0, "n_combos": 1, "min_phi": phi(C), "best_edges": []}
    combos = list(combinations(range(ne), k))
    n_combos = len(combos)
    # shard
    shard_size = max(1, math.ceil(n_combos / max(1, workers * 4)))
    shards = []
    C_list = C.tolist()
    edges_list = edges
    for s in range(0, n_combos, shard_size):
        shards.append((C_list, edges_list, combos[s : s + shard_size]))
    best = float("inf")
    best_combo = None
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_kflip_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            _, b, combo = fut.result()
            if b < best:
                best = b
                best_combo = combo
    best_edges = [edges[t] for t in (best_combo or ())]
    return {
        "k": k,
        "n_combos": n_combos,
        "min_phi": best,
        "best_edges": best_edges,
        "time_s": time.time() - t0,
        "undercuts_paley": best < PHI_PALEY - 1e-9,
        "reaches_m10": best <= M_TARGET + 1e-9,
    }


def _sa_find_optima_worker(args) -> Dict:
    """SA minimising exact Φ; record best matrix and its Φ."""
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    seed, iters, n = args
    rng = np.random.default_rng(seed)
    edges = _edges(n)
    ne = len(edges)
    A = np.zeros((n, n), dtype=np.float64)
    for i, j in edges:
        s = float(rng.choice([-1.0, 1.0]))
        A[i, j] = A[j, i] = s
    cur = phi(A)
    best = cur
    best_A = A.copy()
    T0 = max(2.0, cur * 0.15)
    for t in range(iters):
        T = T0 * (1.0 - t / iters) + 1e-9
        e = int(rng.integers(0, ne))
        i, j = edges[e]
        _flip_edge(A, i, j)
        new = phi(A)
        d = new - cur
        if d <= 0 or rng.random() < math.exp(-d / T):
            cur = new
            if cur < best:
                best = cur
                best_A = A.copy()
        else:
            _flip_edge(A, i, j)  # revert
    return {"phi": float(best), "A": best_A.tolist(), "seed": seed}


def spectral_fingerprint(A: np.ndarray) -> Dict:
    n = A.shape[0]
    op = op_norm_sym(A)
    trA4 = float(np.trace(A @ A @ A @ A))
    delta = trA4 - n * (n - 1) ** 2
    ph = phi(A)
    # ρ = max|xAx| / (n op) = 2Φ/(n op)
    rho = (2.0 * ph) / (n * op) if op > 0 else float("nan")
    r = rho * op / math.sqrt(n - 1)
    return {
        "phi": ph,
        "op": op,
        "op_over_sqrt_n1": op / math.sqrt(n - 1),
        "delta": delta,
        "rho": rho,
        "r": r,
        "is_conference": bool(is_conference_matrix(A)),
    }


def min_switch_hamming(A: np.ndarray, switch_class: List[np.ndarray]) -> int:
    return min(_hamming_edges(A, S) for S in switch_class)


def run_campaign(
    out_dir: Path,
    workers: int = 86,
    k_max: int = 4,
    sa_workers: int = 86,
    sa_iters: int = 8000,
) -> Dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    report: Dict = {
        "n": N,
        "m_10": M_TARGET,
        "Phi_Paley": PHI_PALEY,
        "workers": workers,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    t_all = time.time()

    C = paley_conference_prime_power(3)
    assert is_conference_matrix(C)
    x_hs = halfspace_boolean_vector(3)
    assert abs(rho_from_boolean_evec(C, x_hs) - 1.0) < 1e-9
    report["paley_is_conference"] = True
    report["paley_rho"] = 1.0
    report["paley_phi_exact"] = phi(C)

    # --- 1. maximizer balance ---
    print("[1/4] maximizer balance on Paley C_10 ...", flush=True)
    bal = maximizer_balance_report(C)
    report["maximizer_balance"] = bal
    print(
        f"  Φ={bal['phi']} maximizers={bal['n_maximizers_halfcube']} "
        f"balanced={bal['n_balanced']}/{bal['n_edges']} "
        f"local_opt_1flip={bal['single_flip_local_opt']}",
        flush=True,
    )

    # --- 2. k-flip neighborhoods ---
    print(f"[2/4] exact min Φ at Hamming distance k=0..{k_max} ...", flush=True)
    k_table = []
    for k in range(0, k_max + 1):
        # combinatorial cost: C(45,k); k=5 is 1.2M, fine; k=6 is 8M still ok
        row = min_phi_at_distance_k(C, k, workers=workers if k >= 2 else 1)
        k_table.append(row)
        print(
            f"  k={k}: min_Φ={row['min_phi']}  combos={row['n_combos']}  "
            f"reaches_m={row.get('reaches_m10')}  t={row.get('time_s', 0):.2f}s",
            flush=True,
        )
    report["k_flip_table"] = k_table

    # extend k=5 if k_max>=4 and time allows
    if k_max >= 4:
        print("[2b] k=5 neighborhood ...", flush=True)
        row5 = min_phi_at_distance_k(C, 5, workers=workers)
        k_table.append(row5)
        report["k_flip_table"] = k_table
        print(
            f"  k=5: min_Φ={row5['min_phi']} combos={row5['n_combos']} "
            f"reaches_m={row5.get('reaches_m10')} t={row5.get('time_s', 0):.2f}s",
            flush=True,
        )

    # --- 3. SA hunt for Φ=13 optima + distance to switching class ---
    print(f"[3/4] SA search for optima (workers={sa_workers}, iters={sa_iters}) ...", flush=True)
    t0 = time.time()
    switch_class = all_switchings_of(C, fix_first=True)
    report["switch_class_size"] = len(switch_class)

    jobs = [(10_000 + w, sa_iters, N) for w in range(sa_workers)]
    results = []
    with ProcessPoolExecutor(max_workers=sa_workers) as ex:
        futs = [ex.submit(_sa_find_optima_worker, j) for j in jobs]
        for fut in as_completed(futs):
            results.append(fut.result())
    phis = [r["phi"] for r in results]
    best_phi = min(phis)
    optima = [r for r in results if r["phi"] <= M_TARGET + 1e-9]
    near = [r for r in results if r["phi"] <= PHI_PALEY - 1e-9]

    fingerprints = []
    distances = []
    for r in sorted(results, key=lambda z: z["phi"])[:20]:
        A = np.array(r["A"], dtype=np.float64)
        fp = spectral_fingerprint(A)
        d = min_switch_hamming(A, switch_class)
        fp["hamming_to_paley_switch_class"] = d
        fp["seed"] = r["seed"]
        fingerprints.append(fp)
        if r["phi"] <= M_TARGET + 1e-9:
            distances.append(d)

    report["sa"] = {
        "workers": sa_workers,
        "iters": sa_iters,
        "time_s": time.time() - t0,
        "best_phi_found": best_phi,
        "n_hits_m13": len(optima),
        "n_undercut_paley": len(near),
        "phi_histogram": {
            str(v): sum(1 for p in phis if abs(p - v) < 1e-9)
            for v in sorted(set(phis))
        },
        "min_hamming_among_m13_hits": (min(distances) if distances else None),
        "max_hamming_among_m13_hits": (max(distances) if distances else None),
        "mean_hamming_among_m13_hits": (
            float(np.mean(distances)) if distances else None
        ),
        "top_fingerprints": fingerprints,
    }
    print(
        f"  best_Φ={best_phi}  hits_m13={len(optima)}/{sa_workers}  "
        f"hamming_min={report['sa']['min_hamming_among_m13_hits']}  "
        f"t={report['sa']['time_s']:.1f}s",
        flush=True,
    )

    # --- 4. greedy multi-start from Paley: can we reach 13? ---
    print("[4/4] multi-start greedy descent from switched Paley + random ...", flush=True)
    descent = _greedy_descent_wave(C, switch_class, workers=workers, restarts=workers * 2)
    report["greedy_descent"] = descent
    print(
        f"  greedy best_Φ={descent['best_phi']}  "
        f"min_hamming_at_best={descent.get('min_hamming_at_best')}",
        flush=True,
    )

    report["total_time_s"] = time.time() - t_all

    # summary verdict
    k_reach = None
    for row in k_table:
        if row.get("reaches_m10"):
            k_reach = row["k"]
            break
    report["verdict"] = {
        "single_edge_local_opt_paley": bal["single_flip_local_opt"],
        "min_k_exact_flip_reaches_m13": k_reach,
        "sa_min_hamming_to_switch_class_for_m13": report["sa"][
            "min_hamming_among_m13_hits"
        ],
        "note": (
            "If min_k is large and SA Hamming is large, optima sit far from the "
            "Paley switching class — E(1) cannot be proved by local edge-flip alone."
        ),
    }

    out_path = out_dir / "n10_structure.json"
    out_path.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {out_path}", flush=True)
    return report


def _greedy_one(args) -> Dict:
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    A0_list, seed, max_steps = args
    A = np.array(A0_list, dtype=np.float64)
    rng = np.random.default_rng(seed)
    edges = _edges(A.shape[0])
    ne = len(edges)
    cur = phi(A)
    for _ in range(max_steps):
        # try random edge; accept if improves
        order = rng.permutation(ne)
        improved = False
        for e in order:
            i, j = edges[int(e)]
            _flip_edge(A, i, j)
            new = phi(A)
            if new < cur - 1e-12:
                cur = new
                improved = True
                break
            _flip_edge(A, i, j)
        if not improved:
            break
    return {"phi": float(cur), "A": A.tolist()}


def _greedy_descent_wave(
    C: np.ndarray, switch_class: List[np.ndarray], workers: int, restarts: int
) -> Dict:
    # starts: half from switched Paley (with 0-3 random flips), half random
    jobs = []
    rng = np.random.default_rng(12345)
    edges = _edges(N)
    for r in range(restarts):
        if r < restarts // 2:
            S = switch_class[r % len(switch_class)].copy()
            # scramble k random edges
            k = int(rng.integers(0, 6))
            for _ in range(k):
                i, j = edges[int(rng.integers(0, len(edges)))]
                _flip_edge(S, i, j)
            jobs.append((S.tolist(), 2000 + r, 200))
        else:
            A = np.zeros((N, N))
            for i, j in edges:
                s = float(rng.choice([-1.0, 1.0]))
                A[i, j] = A[j, i] = s
            jobs.append((A.tolist(), 3000 + r, 200))

    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_greedy_one, j) for j in jobs]
        for fut in as_completed(futs):
            results.append(fut.result())

    best_phi = min(r["phi"] for r in results)
    bests = [r for r in results if abs(r["phi"] - best_phi) < 1e-9]
    hammings = [
        min_switch_hamming(np.array(r["A"], dtype=np.float64), switch_class)
        for r in bests
    ]
    # also fingerprint one best
    A_best = np.array(bests[0]["A"], dtype=np.float64)
    return {
        "restarts": restarts,
        "best_phi": best_phi,
        "n_at_best": len(bests),
        "min_hamming_at_best": min(hammings) if hammings else None,
        "max_hamming_at_best": max(hammings) if hammings else None,
        "fingerprint_best": spectral_fingerprint(A_best),
    }


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    out = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    k_max = int(os.environ.get("K_MAX", "4"))
    sa_iters = int(os.environ.get("SA_ITERS", "8000"))
    report = run_campaign(
        out_dir=out,
        workers=workers,
        k_max=k_max,
        sa_workers=workers,
        sa_iters=sa_iters,
    )
    print(json.dumps(report["verdict"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
