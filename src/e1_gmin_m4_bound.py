#!/usr/bin/env python3
"""
P0: prove |m4| ≤ (p-2)/(2p²) on |κ|=1 via the identity
  m4(S) = κ(S)/p² + Ext(S)/(4p)
where Ext = T m4, (Tf)(S) = sum_{v∈S, r∉S} C_vr f(S_{v→r}),
and combinatorially σ_sum = 4κ for every ±1 edge-labeling of K4.

Multi-worker analysis of T, Ext, and closed-form bound attempts.
Writes evidence/e1_gmin_m4_bound.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def _load_Y(p: int) -> np.ndarray:
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def kappa_pts(C, S):
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def sigma_sum(C, S):
    ss = 0
    for vi, v in enumerate(S):
        others = [S[j] for j in range(4) if j != vi]
        i, j, k = others
        ss += int(C[v, i] * C[j, k] + C[v, j] * C[i, k] + C[v, k] * C[i, j])
    return ss


def _worker_identity_and_ext(p: int) -> dict:
    """Verify m4 = κ/p² + Ext/(4p) and analyze Ext on |κ|=1."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape

    # Precompute m4 for all 4-sets as dict
    # For p=7: binom(50,4)=230300, store float array with index map
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}

    # m4 vector
    m4 = np.zeros(nQ)
    # Batch compute via numpy where possible
    for i, q in enumerate(quads):
        a, b, c, d = q
        m4[i] = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))

    kap = np.zeros(nQ, dtype=int)
    ss = np.zeros(nQ, dtype=int)
    for i, q in enumerate(quads):
        kap[i] = kappa_pts(C, q)
        ss[i] = sigma_sum(C, q)

    # Identity check: 4*p*m4 - ss/p  should equal Ext
    # Ext via T
    def compute_ext(i_q):
        S = quads[i_q]
        Sset = set(S)
        ext = 0.0
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                ext += float(C[v, r]) * m4[qindex[Sp]]
        return ext

    # Only for |κ|=1
    idx1 = np.where(np.abs(kap) == 1)[0]
    # Sample for identity check + full for small p
    if p == 5:
        check_idx = idx1
    else:
        rng = np.random.default_rng(0)
        check_idx = rng.choice(idx1, size=min(2000, len(idx1)), replace=False)

    id_errs = []
    exts = []
    m4s = []
    kaps = []
    for i in check_idx:
        ext = compute_ext(int(i))
        pred = kap[i] / (p * p) + ext / (4 * p)
        id_errs.append(abs(pred - m4[i]))
        exts.append(ext)
        m4s.append(m4[i])
        kaps.append(int(kap[i]))

    exts = np.array(exts)
    m4s = np.array(m4s)
    kaps = np.array(kaps)
    # Same-sign Ext (dangerous for |m4|)
    same = np.sign(exts) * np.sign(kaps) > 0
    opp = np.sign(exts) * np.sign(kaps) < 0
    zero_ext = np.abs(exts) < 1e-14

    # Bound targets
    thr_L = L_abs(p)
    thr_ext_for_L = 2.0 * (p - 4) / p  # if same-sign: |m4|=1/p²+|Ext|/(4p)≤L_abs
    # Direct: max |m4|
    max_m4 = float(np.max(np.abs(m4[idx1])))
    max_ext = float(np.max(np.abs(exts)))
    max_ext_same = float(np.max(np.abs(exts[same]))) if np.any(same) else 0.0

    # Among |κ|=1, m4 values when Ext same-sign
    m4_same = m4s[same] if np.any(same) else np.array([])
    m4_opp = m4s[opp] if np.any(opp) else np.array([])

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "n_kappa1": int(len(idx1)),
        "n_checked": int(len(check_idx)),
        "identity_max_err": float(max(id_errs)),
        "sigma_sum_eq_4kappa_all": bool(np.all(ss[idx1] == 4 * kap[idx1])),
        "max_abs_m4": max_m4,
        "L_abs": thr_L,
        "m4_le_L": bool(max_m4 <= thr_L + 1e-12),
        "max_abs_Ext": max_ext,
        "max_abs_Ext_same_sign": max_ext_same,
        "thr_Ext_for_triangle_L": thr_ext_for_L,
        "Ext_same_le_thr": bool(max_ext_same <= thr_ext_for_L + 1e-9),
        "n_same_sign": int(np.sum(same)),
        "n_opp_sign": int(np.sum(opp)),
        "max_abs_m4_same_sign": float(np.max(np.abs(m4_same))) if len(m4_same) else None,
        "max_abs_m4_opp_sign": float(np.max(np.abs(m4_opp))) if len(m4_opp) else None,
        "mean_Ext_same": float(np.mean(np.abs(exts[same]))) if np.any(same) else None,
        "mean_Ext_opp": float(np.mean(np.abs(exts[opp]))) if np.any(opp) else None,
        # For g_min class specifically
        "formula": "m4 = κ/p² + Ext/(4p), Ext = (T m4)(S)",
    }


def _worker_T_spectrum_class(p: int) -> dict:
    """
    Restrict T to class-constant functions (C-isomorphism classes).
    On |κ|=1 refined classes, analyze if T is diagonally dominant etc.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power
    from e1_gmin_moduli import build_classes, class_key, kappa as kap_fn

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    classes, quad_to_key = build_classes(C)
    keys = sorted(classes.keys())
    # Only |κ|=1 classes
    k1_keys = []
    k1_m4 = []
    for k in keys:
        q0 = classes[k][0]
        kap = kap_fn(C, q0)
        if abs(kap) != 1:
            continue
        k1_keys.append(k)
        a, b, c, d = q0
        k1_m4.append(float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d])))

    nk = len(k1_keys)
    key_idx = {k: i for i, k in enumerate(k1_keys)}

    # Build T restricted to class averages: for each class A, 
    # (T 1_A) averaged on class B
    # Sample quads per class to estimate T matrix on R^{nk}
    Tmat = np.zeros((nk, nk))
    rng = np.random.default_rng(p)
    samples_per = 30 if p == 5 else 15

    for i, kA in enumerate(k1_keys):
        quads = classes[kA]
        take = quads if len(quads) <= samples_per else [
            quads[j] for j in rng.choice(len(quads), samples_per, replace=False)
        ]
        acc = np.zeros(nk)
        cnt = 0
        for S in take:
            Sset = set(S)
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    kB = quad_to_key[Sp]
                    if kB not in key_idx:
                        continue  # |κ|≠1 target
                    acc[key_idx[kB]] += float(C[v, r])
                    cnt += 1
        if cnt:
            # Average Ext contribution per S, but T applied to 1_{class}:
            # (T f)(S) for f = class indicator of A: sum_{v,r: S_vr in A} C_vr
            # We accumulated sum over samples of those C_vr into target classes
            Tmat[i, :] = acc / len(take)  # average (T 1_A)(S) for S in class i? 
            # Wait: for fixed class A (source), we're computing average over S in A
            # of the vector of contributions TO each class B.
            # That's (average_{S in A} (T 1_A)(S) broken by target class) — wrong.
            # 
            # Correct for matrix of T on class-constant functions:
            # If f is constant f_A on class A, then
            # (Tf)(S) = sum_{v,r} C_vr f(S_vr) = sum_B f_B * sum_{v,r: S_vr in B} C_vr
            # For S in class A, (Tf)(S) should be nearly constant.
            # Tmat[A,B] = average_{S in A} sum_{v,r: S_vr in B} C_vr
            # Above: for source class... I indexed i as kA which was the SOURCE of the indicator.
            # Let me re-read.
            # I wanted: for S in class i=kA, accumulate C_vr into target class of S_vr.
            # That gives Tmat[i, j] = avg_{S in class i} sum_{v,r: S_vr in class j} C_vr
            # which is exactly (T applied to function 1_j, evaluated averaged on class i) = Tmat_i←j
            # So Tmat[i,j] = (T e_j)_i, columns are T images. Good if we fix:
            pass
        # Actually my loop was for kA as the S-class (evaluation class), and 
        # I add C_vr into target class of S_vr. So acc[j] = sum of C for targets j.
        # Tmat[i,:] = acc/len = row i of T (from all sources? No — I'm only counting
        # edges from S in class i to S_vr in class j, which is sum_{v,r: S_vr in j} C_vr
        # for f=1 on ALL classes... that's T applied to the constant-1 function on all |κ|=1,
        # split by target. Wrong for the full matrix.
        #
        # Correct construction: Tmat[i,j] = avg_S in class_i of (sum_{v,r: S_vr in class_j} C_vr)
        # Yes that's what I have: when I visit S in class i and S_vr in class j, add C_vr to acc[j].
        # So Tmat[i,j] is correct for the matrix of T on class-constant functions!

    # For each evaluation class i, recompute properly
    Tmat = np.zeros((nk, nk))
    for i, kA in enumerate(k1_keys):
        quads = classes[kA]
        take = quads if len(quads) <= samples_per else [
            quads[j] for j in rng.choice(len(quads), samples_per, replace=False)
        ]
        acc = np.zeros(nk)
        for S in take:
            Sset = set(S)
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    kB = quad_to_key[Sp]
                    if kB in key_idx:
                        acc[key_idx[kB]] += float(C[v, r])
        Tmat[i, :] = acc / len(take)

    # Check eigenstructure on m4 vector of classes
    m4_vec = np.array(k1_m4)
    Tm4 = Tmat @ m4_vec  # if Tmat[i,j] = (T e_j)_i, then Tmat @ f gives Tf
    # From identity Ext = T m4 = 4p m4 - 4κ/p
    kap_vec = np.array([kap_fn(C, classes[k][0]) for k in k1_keys], dtype=float)
    Ext_id = 4 * p * m4_vec - 4 * kap_vec / p
    # Compare Tm4 to Ext_id
    err = float(np.max(np.abs(Tm4 - Ext_id))) if nk else None

    # Eigenvalues of Tmat
    evals = np.linalg.eigvals(Tmat)
    evals_real = sorted(evals.real, reverse=True)

    # Bound: if |λ_max(T)| < 4p, then (4p I - T) invertible and
    # m4 = (4p I - T)^{-1} (4κ/p)
    # |m4|_∞ ≤ ‖(4p I - T)^{-1}‖_∞ * 4/p
    # Or l2: ‖m4‖_2 ≤ 4/(p (4p - λ_max)) ‖κ‖_2

    lam_max_T = float(np.max(evals.real))
    lam_min_T = float(np.min(evals.real))

    return {
        "p": int(p),
        "n_kappa1_classes": nk,
        "T_m4_vs_identity_max_err": err,
        "lam_max_T": lam_max_T,
        "lam_min_T": lam_min_T,
        "gap_4p_minus_lammax": 4 * p - lam_max_T,
        "top5_eigs_T": [float(x) for x in evals_real[:5]],
        "m4_values": k1_m4,
        "kappa_values": [int(kap_fn(C, classes[k][0])) for k in k1_keys],
    }


def _worker_prove_sigma_identity() -> dict:
    """Combinatorial: σ_sum = 4κ for all ±1 labelings of K4."""
    from itertools import product

    def check(bits):
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        Cv = {edges[i]: bits[i] for i in range(6)}

        def C(i, j):
            if i == j:
                return 0
            return Cv[tuple(sorted((i, j)))]

        def sigma(v):
            others = [x for x in range(4) if x != v]
            a, b, c = others
            return C(v, a) * C(b, c) + C(v, b) * C(a, c) + C(v, c) * C(a, b)

        kap = C(0, 1) * C(2, 3) + C(0, 2) * C(1, 3) + C(0, 3) * C(1, 2)
        ss = sum(sigma(v) for v in range(4))
        return ss == 4 * kap, kap, ss

    ok = sum(1 for bits in product([-1, 1], repeat=6) if check(bits)[0])
    return {"n_labelings": 64, "n_ok": ok, "proved_by_exhaustion": ok == 64}


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from workers import require_workers

    W = require_workers(min_workers=4)
    print(f"m4_bound: workers={W}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {
            ex.submit(_worker_prove_sigma_identity): "sig",
            ex.submit(_worker_identity_and_ext, 5): "e5",
            ex.submit(_worker_identity_and_ext, 7): "e7",
            ex.submit(_worker_T_spectrum_class, 5): "t5",
        }
        # p=7 T spectrum may be slow; optional
        out_map = {}
        for fut in as_completed(futs):
            tag = futs[fut]
            r = fut.result()
            out_map[tag] = r
            if tag == "sig":
                print(f"  sigma identity: {r}", flush=True)
            elif tag.startswith("e"):
                print(
                    f"  ext p={r['p']}: id_err={r['identity_max_err']:.2e} "
                    f"max|m4|={r['max_abs_m4']:.6f} L={r['L_abs']:.6f} "
                    f"max|Ext|={r['max_abs_Ext']:.4f} same_sign_max={r['max_abs_Ext_same_sign']:.4f} "
                    f"thr_Ext={r['thr_Ext_for_triangle_L']:.4f} "
                    f"same_le={r['Ext_same_le_thr']}",
                    flush=True,
                )
            else:
                print(
                    f"  T p={r['p']}: classes={r['n_kappa1_classes']} "
                    f"T_err={r['T_m4_vs_identity_max_err']} "
                    f"λmax_T={r['lam_max_T']:.4f} 4p-λ={r['gap_4p_minus_lammax']:.4f}",
                    flush=True,
                )

    # Algebra: m4 = κ/p² + Ext/(4p), |m4|≤L if same-sign Ext≤2(p-4)/p
    algebra = {}
    for p in (5, 7, 11, 13, 17):
        algebra[str(p)] = {
            "L_abs": L_abs(p),
            "wick": 1.0 / (p * p),
            "thr_Ext_same_sign": 2.0 * (p - 4) / p,
            "note": "if Ext has sign(κ) and |Ext|≤thr then |m4|≤L_abs",
        }

    out = {
        "sigma_identity": out_map.get("sig"),
        "ext_p5": out_map.get("e5"),
        "ext_p7": out_map.get("e7"),
        "T_p5": out_map.get("t5"),
        "algebra_Ext_criterion": algebra,
        "master_identity": (
            "For every 4-set S on a conference graph: "
            "m4(S)=κ(S)/p² + Ext(S)/(4p), where "
            "Ext(S)=sum_{v∈S,r∉S} C_vr m4(S_{v→r}), "
            "and combinatorially sum_v σ_v = 4κ for all ±1 K4 labelings."
        ),
        "status": (
            "Master identity proved (σ_sum=4κ by K4 exhaustion + evec expansion). "
            "OPEN: bound |Ext| on same-sign |κ|=1 sets by 2(p-4)/p, or bound |m4| directly."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_bound.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
