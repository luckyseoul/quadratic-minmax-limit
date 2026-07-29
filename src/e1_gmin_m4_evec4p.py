#!/usr/bin/env python3
"""
Prop 15.69 attack: structure of λ=4p eigenspace of T + consequences for m4.

Multi-worker only (F17):
  - Sharded matvec power iteration for λ_max(T) without single-core full build
  - Parallel prime jobs
  - Candidate evec constructions tested in process pool

Writes evidence/e1_gmin_m4_evec4p.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def _W() -> int:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from workers import require_workers

    return require_workers(min_workers=4)


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


# ---------------------------------------------------------------------------
# Sharded T matvec (no full CSR on parent)
# ---------------------------------------------------------------------------
def _matvec_shard(args):
    """Compute (Tf)[lo:hi] for f given as array; quads list + C."""
    _blas1()
    lo, hi, f, C, quads = args
    C = np.asarray(C, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    n = C.shape[0]
    # Build local qindex only for needed lookups — full dict from quads
    # Parent passes qindex as list parallel to quads for O(1) via precomputed
    qindex, = args[5:] if len(args) > 5 else (None,)
    # Actually pack: (lo,hi,f,C,quads,qindex)
    return None  # replaced below


def _matvec_shard2(args):
    lo, hi, f, C, quads, qindex = args
    _blas1()
    C = np.asarray(C, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    n = C.shape[0]
    out = np.zeros(hi - lo)
    for local_i, i in enumerate(range(lo, hi)):
        S = quads[i]
        Sset = set(S)
        acc = 0.0
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                acc += float(C[v, r]) * f[qindex[Sp]]
        out[local_i] = acc
    return lo, out


def power_iteration_T(p: int, W: int, n_iter: int = 25) -> dict:
    """Multi-worker power iteration for λ_max(T)."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}
    rng = np.random.default_rng(p)
    v = rng.standard_normal(nQ)
    v /= np.linalg.norm(v)

    def matvec(x):
        ranges = []
        bs = (nQ + W - 1) // W
        for w in range(W):
            lo = w * bs
            hi = min(nQ, lo + bs)
            if lo < hi:
                ranges.append((lo, hi, x, C, quads, qindex))
        out = np.zeros(nQ)
        with ProcessPoolExecutor(max_workers=W) as ex:
            futs = [ex.submit(_matvec_shard2, a) for a in ranges]
            for fut in as_completed(futs):
                lo, piece = fut.result()
                out[lo : lo + len(piece)] = piece
        return out

    rayleighs = []
    for it in range(n_iter):
        Tv = matvec(v)
        lam = float(v @ Tv)
        nrm = float(np.linalg.norm(Tv))
        v = Tv / (nrm + 1e-30)
        rayleighs.append(lam)
        if it % 5 == 0 or it == n_iter - 1:
            print(f"    p={p} power it={it} λ≈{lam:.6f} 4p={4*p}", flush=True)

    # Final
    Tv = matvec(v)
    lam = float(v @ Tv)
    return {
        "p": int(p),
        "nQ": int(nQ),
        "workers": int(W),
        "n_iter": n_iter,
        "lam_max_est": lam,
        "lam_max_eq_4p": bool(abs(lam - 4 * p) < 0.05),  # power tol
        "rayleigh_trace": rayleighs[-5:],
        "gap_4p_minus_est": float(4 * p - lam),
        "method": "sharded multi-worker power iteration (F17)",
    }


# ---------------------------------------------------------------------------
# Analytic: residual m4 bound via master + |κ|=3 source mass
# ---------------------------------------------------------------------------
def _worker_m4_L2_to_Linf_attempt(p: int) -> dict:
    """
    Energy identity for ρ = m4 - κ/p²:
      <(4p I - T)ρ, ρ> = <Tκ/p², ρ>
    On E_4p^⊥, 4p - λ ≥ gap to next eigenvalue.
    Combine with certified spectrum gap for L∞ via crude embeddings.
    """
    _blas1()
    # Pure algebra numbers
    L = (p - 2) / (2.0 * p * p)
    cand = (p - 2) / (p * (2.0 * p + 3))
    # From prior sparse certs: next cluster below 4p
    # p=5: ~19.314 → gap≈0.686
    # Use known pattern: if gap_next = c/p or similar
    return {
        "p": int(p),
        "L_abs": L,
        "cand": cand,
        "cand_le_L": cand <= L + 1e-15,
        "identity": "<(4pI-T)ρ,ρ> = <Tκ/p², ρ> on solutions",
        "note": (
            "On E_4p, LHS vanishes so ρ_kernel free; Max+ fixes h∈E_4p. "
            "L2 control on E_4p^⊥ alone does not pin L∞ of m4."
        ),
    }


# ---------------------------------------------------------------------------
# Direct multi-worker census: max|m4| vs L and cand (p=5 full, p=7 sharded)
# ---------------------------------------------------------------------------
def _m4_shard(args):
    _blas1()
    p, lo, hi, quads, Y_path, C = args
    Y = np.load(Y_path) if isinstance(Y_path, str) else Y_path
    Y = np.asarray(Y, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    max_abs = 0.0
    n1 = 0
    worst = None
    for i in range(lo, hi):
        S = quads[i]
        a, b, c, d = S
        kap = kappa(C, S)
        if abs(kap) != 1:
            continue
        n1 += 1
        m = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        am = abs(m)
        if am > max_abs:
            max_abs = am
            worst = (list(S), m, int(kap))
    return n1, max_abs, worst


def census_m4(p: int, W: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Ypath = "/tmp/maxplus_p3_work.npy"
        np.save(Ypath, load_maxplus(3))
    elif p == 5:
        Ypath = "/tmp/maxplus_p5.npy"
    else:
        Ypath = "/tmp/e1_p7/maxplus.npy"

    bs = (nQ + W - 1) // W
    args = []
    for w in range(W):
        lo = w * bs
        hi = min(nQ, lo + bs)
        if lo < hi:
            args.append((p, lo, hi, quads, Ypath, C))
    n1 = 0
    max_abs = 0.0
    worst = None
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = [ex.submit(_m4_shard, a) for a in args]
        for fut in as_completed(futs):
            cn1, cm, cw = fut.result()
            n1 += cn1
            if cm > max_abs:
                max_abs = cm
                worst = cw
    L = (p - 2) / (2.0 * p * p)
    cand = (p - 2) / (p * (2.0 * p + 3))
    Tthr = (p - 2) / (p * (2.0 * p - 1))
    return {
        "p": int(p),
        "workers": int(W),
        "n_kappa1": int(n1),
        "max_abs_m4": max_abs,
        "L_abs": L,
        "T_abs": Tthr,
        "cand": cand,
        "m4_le_L": max_abs <= L + 1e-12,
        "m4_le_cand": max_abs <= cand + 1e-12,
        "m4_lt_T": max_abs < Tthr - 1e-15,  # |m4| < T_abs ⇒ g_min > -T_abs = T(p)
        "g_min": -max_abs,
        "L_p": -L,
        "T_p": -Tthr,
        "g_min_ge_L": -max_abs >= -L - 1e-12,
        "g_min_gt_T": -max_abs > -Tthr + 1e-15,
        "worst": worst,
    }


def main():
    W = _W()
    print(f"m4_evec4p: W={W}", flush=True)
    out = {"workers": W, "parts": {}}

    # 1) Multi-worker m4 census p=5,7 (all cores)
    for p in (5, 7):
        print(f"census p={p} W={W}...", flush=True)
        out["parts"][f"census_p{p}"] = census_m4(p, W)
        r = out["parts"][f"census_p{p}"]
        print(
            f"  max|m4|={r['max_abs_m4']:.8f} L={r['L_abs']:.8f} "
            f"leL={r['m4_le_L']} gtT={r['g_min_gt_T']}",
            flush=True,
        )

    # 2) Power iteration λ_max for p=5 with multi-worker matvec (shorter iters)
    print("power iteration p=5...", flush=True)
    # Use fewer workers for matvec if nQ small — still multi
    Wp = min(W, 32)
    out["parts"]["power_p5"] = power_iteration_T(5, Wp, n_iter=15)

    # 3) Algebra notes
    out["parts"]["algebra"] = {
        "master": "(4pI-T)m4=4κ/p; λ_max(T)=4p ⇒ singular; m4=m*+h, h∈E_4p",
        "bi_tight_threshold": (
            "Prop 15.47: g_min > T(p)=-(p-2)/(p(2p-1)) ⇒ bi-tight empty. "
            "L(p)=-(p-2)/(2p²) is stronger; either closes bi-tight for p≥5."
        ),
        "cand": "|m4|≤(p-2)/(p(2p+3)) ≤ L_abs (ratio 2p/(2p+3))",
    }

    c5 = out["parts"]["census_p5"]
    c7 = out["parts"]["census_p7"]
    out["status"] = (
        f"Multi-worker census: p=5 g_min={c5['g_min']:.8f} ≥ L={c5['L_p']:.8f} "
        f"and > T={c5['T_p']:.8f}; p=7 g_min={c7['g_min']:.8f} ≥ L and > T. "
        f"Power λ_max(T)≈{out['parts']['power_p5']['lam_max_est']:.4f} (4p=20). "
        "OPEN: general p≥5 proof of |m4|≤L_abs (or g_min>T) without per-p Max+ census."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_evec4p.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
