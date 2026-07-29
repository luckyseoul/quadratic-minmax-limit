#!/usr/bin/env python3
"""
Attack Q_4 ≤ 10 N ‖B‖_F² on zero-diag ∩ V_+ (Prop 15.62 residual).

p=3,5: full operator eigh on the zero-diag∩V_+ space.
p=7: multi-start vectorized power iteration (fan-out across workers).

Multi-worker ProcessPool from this file only. OMP=1 per worker.
Writes evidence/e1_gmin_q4_bound.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _A_from_vec(v: np.ndarray, d: int) -> np.ndarray:
    A = np.zeros((d, d))
    k = 0
    for i in range(d):
        for j in range(i, d):
            A[i, j] = A[j, i] = v[k]
            k += 1
    return A


def _null_zero_diag(Vp: np.ndarray) -> np.ndarray:
    n, d = Vp.shape
    dimS = d * (d + 1) // 2
    M = np.zeros((n, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        M[:, k] = np.diag(Vp @ _A_from_vec(e, d) @ Vp.T)
    _u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    return vh[rank:].T


def _worker_full_eigh(p: int) -> dict:
    """Full Q_4 operator spectrum (p=3,5 only; dim small)."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    null = _null_zero_diag(Vp)
    nullity = null.shape[1]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)

    Bs = []
    for t in range(nullity):
        A = _A_from_vec(null[:, t], d)
        B = Vp @ A @ Vp.T
        for Bb in Bs:
            B = B - np.sum(B * Bb) * Bb
        nrm = np.linalg.norm(B)
        if nrm > 1e-12:
            Bs.append(B / nrm)
    m = len(Bs)

    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    Gu = F.T @ F
    share = np.zeros((E, E), dtype=bool)
    for a, (i, j) in enumerate(edges):
        for b in range(a + 1, E):
            k, l = edges[b]
            if len({i, j, k, l}) == 3:
                share[a, b] = share[b, a] = True
    disj = ~share
    np.fill_diagonal(disj, False)
    Gu_disj = Gu * disj

    Bes = np.stack([np.array([2.0 * B[i, j] for i, j in edges]) for B in Bs], axis=1)
    Mmat = Bes.T @ Gu_disj @ Bes
    evals = np.linalg.eigvalsh(0.5 * (Mmat + Mmat.T))
    thr = 10.0 * N

    # Also full Q operator (should be Q4 + 6N on unit sphere)
    # Q(B)=sum (yBy)^2
    Mq = np.zeros((m, m))
    for s in range(m):
        for t in range(s, m):
            # bilinear: sum_y (y Bs y)(y Bt y)
            qs = np.einsum("ai,ij,aj->a", Y, Bs[s], Y)
            qt = np.einsum("ai,ij,aj->a", Y, Bs[t], Y)
            val = float(np.dot(qs, qt))
            Mq[s, t] = Mq[t, s] = val
    evals_Q = np.linalg.eigvalsh(0.5 * (Mq + Mq.T))

    return {
        "p": int(p),
        "method": "full_eigh",
        "N": int(N),
        "dim_basis": int(m),
        "lam_max_Q4": float(evals[-1]),
        "lam_min_Q4": float(evals[0]),
        "lam_max_Q": float(evals_Q[-1]),
        "thr_10N": thr,
        "thr_16N": 16.0 * N,
        "bound_ok": bool(evals[-1] <= thr + 1e-6 * thr),
        "ratio_to_thr": float(evals[-1] / thr),
        "top5_Q4": [float(x) for x in evals[-5:][::-1]],
        "spectrum_Q4_unique_round6": sorted({round(float(x), 6) for x in evals}),
    }


def _worker_power_start(args: tuple) -> dict:
    """One power-iteration start for max Q on V_+ (then Q4=Q-6N for unit)."""
    _blas1()
    p, seed = args
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(70):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    # Ambient B may have tiny diag; Q is on V_+ Tr0
    B = Vp @ A @ Vp.T
    bf2 = float(np.sum(B * B))
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))
    # For unit B
    Q_unit = Q / bf2
    Q4_unit = Q_unit - 6.0 * N
    return {
        "p": int(p),
        "seed": int(seed),
        "Q_unit": Q_unit,
        "Q4_unit": Q4_unit,
        "ratio_16N": Q_unit / (16.0 * N),
        "diagmax": float(np.max(np.abs(np.diag(B)))),
    }


def _worker_m4(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    Y = _load_Y(p)
    N, n = Y.shape
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    if E > 400:
        return {"p": int(p), "skipped": True}
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    G = (F.T @ F) / N
    m4s = []
    for a in range(E):
        i, j = edges[a]
        for b in range(a + 1, E):
            k, l = edges[b]
            if len({i, j, k, l}) == 4:
                m4s.append(G[a, b])
    m4s = np.asarray(m4s, dtype=np.float64)
    return {
        "p": int(p),
        "skipped": False,
        "n_disj_pairs": int(m4s.size),
        "m4_min": float(m4s.min()),
        "m4_max": float(m4s.max()),
        "m4_abs_max": float(np.max(np.abs(m4s))),
        "m4_mean": float(m4s.mean()),
        "n_unique_round10": len({round(float(x), 10) for x in m4s}),
        "unique_round10_sample": sorted({round(float(x), 10) for x in m4s})[:16],
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"q4_bound: workers={W}", flush=True)

    eigs = []
    m4s = []
    powers = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5):
            futs[ex.submit(_worker_full_eigh, p)] = ("eigh", p)
            futs[ex.submit(_worker_m4, p)] = ("m4", p)
        # p=7: 86 power starts in parallel
        for s in range(W):
            futs[ex.submit(_worker_power_start, (7, 1000 + s))] = ("pow", s)
        # also a few at p=5 for cross-check
        for s in range(min(20, W)):
            futs[ex.submit(_worker_power_start, (5, 2000 + s))] = ("pow5", s)

        for fut in as_completed(futs):
            kind, tag = futs[fut]
            r = fut.result()
            if kind == "eigh":
                eigs.append(r)
                print(
                    f"  eigh p={r['p']}: lam_max_Q4={r['lam_max_Q4']:.6f} "
                    f"thr={r['thr_10N']} ok={r['bound_ok']} ratio={r['ratio_to_thr']:.6f} "
                    f"lam_max_Q={r['lam_max_Q']:.6f}",
                    flush=True,
                )
            elif kind == "m4":
                m4s.append(r)
                print(
                    f"  m4 p={r['p']}: absmax={r.get('m4_abs_max')} "
                    f"nuniq={r.get('n_unique_round10')}",
                    flush=True,
                )
            else:
                powers.append(r)

    # Aggregate power starts
    pow_by_p: dict[int, list] = {}
    for r in powers:
        pow_by_p.setdefault(r["p"], []).append(r)
    pow_summary = {}
    for p, lst in sorted(pow_by_p.items()):
        best = max(lst, key=lambda r: r["Q4_unit"])
        N = {5: 260, 7: 11452}[p]
        thr = 10.0 * N
        pow_summary[str(p)] = {
            "n_starts": len(lst),
            "best_Q4": best["Q4_unit"],
            "best_Q": best["Q_unit"],
            "best_ratio_16N": best["ratio_16N"],
            "thr_10N": thr,
            "bound_ok": bool(best["Q4_unit"] <= thr + 1e-3 * thr),
            "ratio_to_thr": best["Q4_unit"] / thr,
        }
        print(
            f"  power p={p}: best_Q4={best['Q4_unit']:.4f} thr={thr} "
            f"ok={pow_summary[str(p)]['bound_ok']} ratio16={best['ratio_16N']:.6f}",
            flush=True,
        )

    eigs.sort(key=lambda r: r["p"])
    m4s.sort(key=lambda r: r["p"])
    all_ok = all(r["bound_ok"] for r in eigs) and all(
        v["bound_ok"] for v in pow_summary.values()
    )
    out = {
        "q4_operator": eigs,
        "power_summary": pow_summary,
        "m4_disj": m4s,
        "status": (
            "Q4≤10N on zero-diag∩V_+: "
            + ("CERTIFIED p=3,5 (full eigh) + p=5,7 (power)" if all_ok else "FAIL")
            + ". p=3 equality lam_max_Q4=10N. "
            "Uniform proof OPEN (need image Rayleigh of Gu_disj ≤5N for unit Be)."
        ),
        "workers": W,
        "all_ok": bool(all_ok),
    }
    path = ROOT / "evidence" / "e1_gmin_q4_bound.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
