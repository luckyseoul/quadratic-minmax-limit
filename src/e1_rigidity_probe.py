#!/usr/bin/env python3
"""
E(1) rigidity probe on ρ=1 family.

Theory target (Prop 15.20): if every near-minimiser of r is within
k=o(n) edge flips of the conference switching class, then
m_n = Φ(C) - o(n^{3/2}).

This script:
1. Starts from Paley C (ρ=1), runs SA minimising exact-rescorable proxies
2. Measures Hamming distance to the switching class of C (min over switches)
3. Records r-proxy via phi_local + exact Φ when n is small enough

Does NOT claim a proof. Writes JSON evidence for structural constraints.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    is_conference_matrix,
    op_norm_sym,
    paley_conference_prime_power,
    phi,
    phi_local,
)


def switching_distance(A: np.ndarray, C: np.ndarray, n_switch: int = 200, seed: int = 0) -> int:
    """
    Approximate min over diagonal switchings ε of Hamming edge-distance
    between A^ε and C. Exact min over 2^{n-1} is hard; we use:
    - greedy coordinate ascent on agreement
    - random multi-start
    """
    n = A.shape[0]
    rng = np.random.default_rng(seed)
    # Fix C; switch A by ε: A_ij -> ε_i A_ij ε_j
    # Agree count = sum_{i<j} 1[ε_i A_ij ε_j == C_ij]
    # Equivalent: max_ε sum_{i<j} ε_i ε_j (A_ij C_ij)
    B = A * C  # Hadamard; diagonal 0

    def agree_from_eps(eps: np.ndarray) -> int:
        # number of edges where A^eps == C
        M = eps[:, None] * B * eps[None, :]
        # M_ij = ε_i A_ij C_ij ε_j; equals +1 iff A^ε and C agree
        iu = np.triu_indices(n, 1)
        return int(np.sum(M[iu] > 0))

    e = n * (n - 1) // 2
    best_agree = -1
    for trial in range(n_switch):
        if trial == 0:
            eps = np.ones(n)
        else:
            eps = rng.choice([-1.0, 1.0], size=n)
            eps[0] = 1.0
        # greedy flip coordinates to improve agreement (max cut style on B)
        improved = True
        while improved:
            improved = False
            order = rng.permutation(n)
            for i in order:
                if i == 0:
                    continue
                # delta if flip eps[i]
                # contribution of edges (i,j): currently eps_i B_ij eps_j
                # after flip: -eps_i B_ij eps_j, change = -2 eps_i B_ij eps_j
                row = B[i] * eps
                cur_sum = float(eps[i] * np.sum(row) - eps[i] * B[i, i])  # off-diag
                # flipping i changes total sum_{edges} by -2 * sum_{j!=i} eps_i B_ij eps_j
                s = float(np.dot(B[i], eps) - B[i, i] * eps[i])
                delta_agree_proxy = -2.0 * eps[i] * s
                if delta_agree_proxy > 1e-9:
                    eps[i] *= -1
                    improved = True
        ag = agree_from_eps(eps)
        if ag > best_agree:
            best_agree = ag
    return e - best_agree  # Hamming edge distance


def run_one(args: Tuple[int, int, int]) -> Dict:
    p, seed, steps = args
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Phi = 0.5 * n * p
    rng = np.random.default_rng(seed)
    # init: Paley with random switch + random flips
    s = rng.choice([-1.0, 1.0], size=n)
    A = s[:, None] * C * s[None, :]
    k0 = int(rng.integers(0, max(1, n // 2)))
    for _ in range(k0):
        i, j = rng.integers(0, n, size=2)
        if i != j:
            A[i, j] *= -1
            A[j, i] *= -1

    def energy(M: np.ndarray) -> float:
        # minimise local phi as SA energy (not a certified UB)
        return float(phi_local(M, restarts=8, rng=rng))

    cur_e = energy(A)
    best_e = cur_e
    best_A = A.copy()
    T0, T1 = float(Phi), 0.5
    for t in range(steps):
        T = T0 * (T1 / T0) ** (t / max(steps - 1, 1))
        i, j = rng.integers(0, n, size=2)
        if i == j:
            continue
        A[i, j] *= -1
        A[j, i] *= -1
        new_e = energy(A)
        if new_e <= cur_e or rng.random() < np.exp((cur_e - new_e) / max(T, 1e-9)):
            cur_e = new_e
            if new_e < best_e:
                best_e = new_e
                best_A = A.copy()
        else:
            A[i, j] *= -1
            A[j, i] *= -1

    # exact Φ for finalist (n=10 cheap; n=26 heavy but doable)
    exact = float(phi(best_A))
    op = float(op_norm_sym(best_A))
    r = 2.0 * exact / (n * np.sqrt(n - 1))
    ham = switching_distance(best_A, C, n_switch=80, seed=seed)
    return {
        "p": p,
        "n": n,
        "seed": seed,
        "phi_local_best": best_e,
        "phi_exact": exact,
        "Phi_paley": Phi,
        "gap": Phi - exact,
        "op": op,
        "r": r,
        "hamming_to_C_switch_approx": ham,
        "is_conference": bool(is_conference_matrix(best_A)),
    }


def main() -> None:
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    W = max(1, (os.cpu_count() or 2) - 2)
    jobs = [(p, s, steps) for s in range(reps)]
    t0 = time.time()
    rows: List[Dict] = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        for r in ex.map(run_one, jobs, chunksize=1):
            rows.append(r)
            print(
                f"seed={r['seed']} exactΦ={r['phi_exact']} gap={r['gap']} "
                f"r={r['r']:.4f} ham≈{r['hamming_to_C_switch_approx']}",
                flush=True,
            )
    gaps = [r["gap"] for r in rows]
    hams = [r["hamming_to_C_switch_approx"] for r in rows]
    report = {
        "p": p,
        "n": p * p + 1,
        "reps": reps,
        "steps": steps,
        "workers": W,
        "min_exact_phi": min(r["phi_exact"] for r in rows),
        "max_gap_under_paley": max(gaps),
        "min_gap": min(gaps),
        "mean_hamming_approx": float(np.mean(hams)),
        "min_hamming_approx": int(min(hams)),
        "rows": rows,
        "elapsed_sec": time.time() - t0,
        "note": (
            "If min_exact_phi < Phi_paley, certified undercut (E(1) gap). "
            "Hamming is approximate min over switches. Not a proof of E(1)."
        ),
    }
    out = Path(__file__).resolve().parents[1] / "evidence" / f"e1_rigidity_p{p}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: report[k] for k in report if k != "rows"}, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
