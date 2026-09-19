#!/usr/bin/env python3
"""Finite census of Proposition 6.6 residue (6.20) and Paley-skew diamond slack.

Implication this would establish: for each tested signing A of order n<=14,
whether the unshielded set (6.20) is empty, and if not, whether the explicit
Paley-skew R of Prop 6.6 still satisfies the mixed-state diamond
    |Q_A(x)+Q_A(y)| + |x^T R y|  <=  2 sqrt(2) Phi(A)
on that set (up to the n^{3/2} Dini slack, which at these n is order-1).
A nonempty residue on which Paley-R fails is a finite witness that the
global estimate (6.24) is the live gap. Emptiness is a finite observation
only; it does not close the ray.

Backend: one NumPy process per signing. Pair matrices at n=14 are 8192^2
GEMMs (~0.5 GiB). GPU unused (driver down). Serial inner GEMM is the
natural unit; fan-out is over independent signings via ProcessPool.
"""
from __future__ import annotations

import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from paley_structural_lift_family import paley_conference  # noqa: E402

A_CONST = (math.sqrt(2.0) - 1.0) / math.pi
RHO = A_CONST * A_CONST


def next_prime_3mod4(n: int) -> int:
    def is_prime(p: int) -> bool:
        if p < 2:
            return False
        if p % 2 == 0:
            return p == 2
        d = 3
        while d * d <= p:
            if p % d == 0:
                return False
            d += 2
        return True

    q = n if n % 2 else n + 1
    if q % 4 != 3:
        q += 2 if (q % 4 == 1) else 0
        while q % 4 != 3:
            q += 2
    while not is_prime(q):
        q += 4
    return q


def paley_tournament(q: int) -> np.ndarray:
    """Skew Paley tournament of prime order q = 3 mod 4."""
    if q % 4 != 3:
        raise ValueError("need q = 3 mod 4")
    chi = np.zeros(q, dtype=np.int64)
    for k in range(1, q):
        chi[k] = 1 if pow(k, (q - 1) // 2, q) == 1 else -1
    T = np.zeros((q, q), dtype=np.int64)
    for i in range(q):
        for j in range(q):
            if i != j:
                T[i, j] = chi[(j - i) % q]
    return T


def balance_tournament(T: np.ndarray) -> np.ndarray:
    """Degree-balance by outdegree transfers, Prop 6.6. Mutates a copy."""
    R = T.copy()
    n = R.shape[0]
    # Potential: sum |row-sum| (odd n) or sum max(|row-sum|-1, 0) (even n).
    def row_sums(M):
        return M.sum(axis=1)

    guard = n * n
    while guard > 0:
        guard -= 1
        s = row_sums(R)
        if n % 2:
            pot = int(np.abs(s).sum())
            if pot == 0:
                break
        else:
            pot = int(np.maximum(np.abs(s) - 1, 0).sum())
            if pot == 0:
                break
        i = int(np.argmax(s))
        j = int(np.argmin(s))
        if i == j:
            break
        if R[i, j] == 1:
            R[i, j] = -1
            R[j, i] = 1
        else:
            # reverse a directed two-path i -> k -> j
            mids = np.where((R[i] == 1) & (R[:, j] == 1))[0]
            mids = mids[(mids != i) & (mids != j)]
            if mids.size == 0:
                break
            k = int(mids[0])
            R[i, k] = -1
            R[k, i] = 1
            R[k, j] = -1
            R[j, k] = 1
    return R


def paley_skew_R(n: int, z: np.ndarray) -> np.ndarray:
    q = next_prime_3mod4(n)
    T = paley_tournament(q)
    Rn = T[:n, :n].copy()
    Rn = balance_tournament(Rn)
    D = np.diag(z.astype(np.int64))
    return D @ Rn @ D


def boolean_cube(n: int) -> np.ndarray:
    """(2^{n-1}, n) with x_0 = +1."""
    m = 1 << (n - 1)
    idx = np.arange(m, dtype=np.int64)
    X = np.empty((m, n), dtype=np.int8)
    X[:, 0] = 1
    for c in range(1, n):
        X[:, c] = (1 - 2 * ((idx >> (c - 1)) & 1)).astype(np.int8)
    return X


def cycle5() -> np.ndarray:
    n = 5
    A = np.ones((n, n), dtype=np.int64)
    np.fill_diagonal(A, 0)
    for i in range(n):
        A[i, (i + 1) % n] = -1
        A[(i + 1) % n, i] = -1
    return A


def random_signing(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    U = rng.integers(0, 2, size=(n, n)) * 2 - 1
    A = np.triu(U, 1)
    A = A + A.T
    np.fill_diagonal(A, 0)
    return A.astype(np.int64)


def census_one(spec: dict) -> dict:
    t0 = time.time()
    kind = spec["kind"]
    n = spec["n"]
    if kind == "conference":
        A = paley_conference(spec["q"])
        assert A.shape[0] == n
        label = f"C{n}"
    elif kind == "cycle5":
        A = cycle5()
        label = "C5-cycle"
    elif kind == "npy":
        A = np.load(spec["path"])
        n = int(A.shape[0])
        spec = dict(spec)
        spec["n"] = n
        label = spec.get("label", Path(spec["path"]).stem)
    elif kind == "random":
        A = random_signing(n, spec["seed"])
        label = f"rand n={n} seed={spec['seed']}"
    else:
        raise ValueError(kind)

    X = boolean_cube(n).astype(np.float64)
    Q = 0.5 * np.sum((X @ A.astype(np.float64)) * X, axis=1)
    absQ = np.abs(Q)
    M = float(absQ.max())
    z_idx = int(np.argmax(absQ))
    z = X[z_idx]
    # Hamming to {z, -z}: n even would include antipode in the x_0=+1 cube
    # only if z_0 = -1, which we forbade. Distance to -z is n - d(z).
    dots_z = X @ z
    h = np.minimum((n - dots_z) / 2.0, (n + dots_z) / 2.0)

    inner = X @ X.T  # (m, m)
    d = (n - inner) / 2.0
    cut = d * (n - d)
    hx = h[:, None]
    hy = h[None, :]
    Qsum = np.abs(Q[:, None] + Q[None, :])

    thr_h = RHO * n
    thr_cut = RHO * n * n
    thr_prod = (RHO / 4.0) * n * n
    n32 = n ** 1.5
    two_sqrt2_M = 2.0 * math.sqrt(2.0) * M
    energy_thr = two_sqrt2_M - n32
    energy_live = energy_thr > 0.0

    residue = (
        (hx > thr_h)
        & (hy > thr_h)
        & (cut > thr_cut)
        & ((hx * hy) > thr_prod)
        & (Qsum > energy_thr)
    )
    # skip i=j? diagonal has d=0 so cut=0, already excluded by cut>thr
    n_res = int(residue.sum())
    m_states = X.shape[0]
    n_pairs = m_states * m_states

    out = {
        "label": label,
        "kind": kind,
        "n": n,
        "Phi": M,
        "alpha": M / n32,
        "energy_threshold": energy_thr,
        "energy_condition_live": energy_live,
        "rho_n": thr_h,
        "n_states": m_states,
        "n_pairs": n_pairs,
        "n_residue": n_res,
        "residue_frac": n_res / n_pairs,
        "max_|Qx+Qy|_over_n32": float(Qsum.max() / n32),
        "2sqrt2_alpha": two_sqrt2_M / n32,
    }

    if n_res > 0:
        R = paley_skew_R(n, z)
        C = X @ R.astype(np.float64) @ X.T
        diamond = Qsum + np.abs(C)
        slack = diamond - two_sqrt2_M
        res_slack = slack[residue]
        out.update(
            {
                "paley_q": int(next_prime_3mod4(n)),
                "max_diamond_res_over_n32": float(diamond[residue].max() / n32),
                "max_slack_res": float(res_slack.max()),
                "max_slack_res_over_n32": float(res_slack.max() / n32),
                "frac_res_slack_positive": float((res_slack > 1e-9).mean()),
                "max_|xRy|_res_over_n32": float(np.abs(C[residue]).max() / n32),
                "max_slack_all_over_n32": float(slack.max() / n32),
            }
        )
    out["seconds"] = round(time.time() - t0, 3)
    return out


SPECS = [
    {"kind": "cycle5", "n": 5},
    {"kind": "conference", "n": 6, "q": 5},
    {"kind": "conference", "n": 10, "q": 9},
    {"kind": "conference", "n": 14, "q": 13},
    {
        "kind": "npy",
        "n": 8,
        "label": "exact_m8",
        "path": str(
            ROOT / "evidence" / "residue_620_census_20260919" / "A8_exact.npy"
        ),
    },
    {"kind": "random", "n": 8, "seed": 1},
    {"kind": "random", "n": 8, "seed": 2},
    {"kind": "random", "n": 11, "seed": 1},
    {"kind": "random", "n": 12, "seed": 1},
    {"kind": "random", "n": 13, "seed": 1},
    {"kind": "random", "n": 14, "seed": 1},
]


def main() -> None:
    outdir = ROOT / "evidence" / "residue_620_census_20260919"
    outdir.mkdir(parents=True, exist_ok=True)
    workers = min(len(SPECS), 10)
    print(f"backend=ProcessPool workers={workers} (one signing each); GPU unused", flush=True)
    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(census_one, s): s for s in SPECS}
        for fut in as_completed(futs):
            rec = fut.result()
            results.append(rec)
            print(json.dumps(rec, sort_keys=True), flush=True)
    results.sort(key=lambda r: (r["n"], r["label"]))
    payload = {
        "rho": RHO,
        "a": A_CONST,
        "workers": workers,
        "seconds": round(time.time() - t0, 3),
        "results": results,
    }
    path = outdir / "census.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print("WROTE", path, flush=True)


if __name__ == "__main__":
    main()
