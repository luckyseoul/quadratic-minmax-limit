#!/usr/bin/env python3
"""
6-point linear Aut-contractions + |m|≤1, and SOS-4, cannot prove GLOBAL QVAR.

Unnumbered kill.  Does **not** prove F̂(ψ)≥0.  Does **not** flip leftover
flags.  Ridge/Boolean reconstruction beyond linear 6-point remains OPEN.

KILL A (p=5 live).  Even-degree identity (Cy−py)_k y_S = 0 for |S|=5
relates m₆ to m₄.  Together with the 4-point master Tm=4pm−4κ/p, Aut
class functions, and |m₄|,|m₆|≤1, the LP at p=5 has kernel dimension 3
and min ⟨m₄,κ_{A_ψ}⟩ = −101/4 < 0, while true Max+ pairing is positive.
Fail: claim that min is ≥0.

KILL B (p=5 live).  SOS-4 / degree-4 moment SDP along the 15.590 deg-4
kernel (dim 1) is feasible at a point with pairing −45/4 < 0.  Fail:
claim every PSD-4 point has pairing ≥0.  Same shape as 15.590 H for
leftover 3, now for the QVAR objective.

KILL C (p=7 recorded, not a p-law).  The complete 15.590 joint deg-4+6
equivariant system has kernel dimension 4; box LP min pairing
−10633/8 < 0, true pairing positive.  Linear 6-point theory does not
determine the sign at the first prime where deg-6 leaves a kernel.

OPEN: a constraint strictly stronger than linear 6-point + box and
SOS-4 (Boolean support / SOS-6 / coupled ridge reconstruction).
"""
from __future__ import annotations

import itertools
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import (  # noqa: E402
    MuLab,
    paley_conference,
    signed_orbits,
)
from e1_gmin_qvar_box_master import (  # noqa: E402
    A_psi_matrix,
    make_psi,
    permutation_aut_gens,
    true_maxplus_pairing,
)


def _enc(S, n):
    e = 0
    for x in S:
        e = e * n + int(x)
    return e


def sixpoint_master_box_lp(p: int) -> dict:
    """Aut-inv 6-point contractions + 4-point master + box.  p=5 certificate."""
    C = paley_conference(p)
    n = C.shape[0]
    gens = permutation_aut_gens(p, C)
    S4 = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    S5 = np.array(list(combinations(range(n), 5)), dtype=np.int64)
    S6 = np.array(list(combinations(range(n), 6)), dtype=np.int64)
    lab4, _, _ = signed_orbits(S4, gens, n, twist=False)
    lab5, _, _ = signed_orbits(S5, gens, n, twist=False)
    lab6, _, _ = signed_orbits(S6, gens, n, twist=False)

    def compress(lab):
        u, loc = np.unique(lab, return_inverse=True)
        return loc, np.bincount(loc).astype(np.float64), len(u)

    loc4, sz4, n4 = compress(lab4)
    loc6, sz6, n6 = compress(lab6)
    loc5, _, _ = compress(lab5)
    enc4 = np.array([_enc(S, n) for S in S4], dtype=np.int64)
    o4 = np.argsort(enc4)
    e4s = enc4[o4]
    enc6 = np.array([_enc(S, n) for S in S6], dtype=np.int64)
    o6 = np.argsort(enc6)
    e6s = enc6[o6]

    def i4(*xs):
        key = _enc(sorted(xs), n)
        return loc4[int(o4[np.searchsorted(e4s, key)])]

    def i4w(*xs):
        key = _enc(sorted(xs), n)
        return int(o4[np.searchsorted(e4s, key)])

    def i6(*xs):
        key = _enc(sorted(xs), n)
        return loc6[int(o6[np.searchsorted(e6s, key)])]

    reps5 = []
    seen = set()
    for t, S in enumerate(S5):
        l = int(loc5[t])
        if l in seen:
            continue
        seen.add(l)
        reps5.append(tuple(int(x) for x in S))
    rows = []
    rhs = []
    NV = n4 + n6
    kap = (
        C[S4[:, 0], S4[:, 1]] * C[S4[:, 2], S4[:, 3]]
        + C[S4[:, 0], S4[:, 2]] * C[S4[:, 1], S4[:, 3]]
        + C[S4[:, 0], S4[:, 3]] * C[S4[:, 1], S4[:, 2]]
    ).astype(np.float64)
    Torb = np.zeros((n4, n4))
    Cf = C.astype(np.float64)
    for t, S in enumerate(S4):
        Sset = set(int(x) for x in S)
        li = loc4[t]
        for v in S:
            v = int(v)
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted([int(x) for x in S if x != v] + [r]))
                Torb[li, loc4[i4w(*Sp)]] += Cf[v, r]
    Torb /= sz4[:, None]
    avg_k = np.bincount(loc4, weights=kap, minlength=n4) / np.maximum(sz4, 1)
    M4 = Torb - 4.0 * p * np.eye(n4)
    rhs4 = -4.0 * avg_k / p
    for irow in range(n4):
        cv = np.zeros(NV)
        cv[:n4] = M4[irow]
        rows.append(cv)
        rhs.append(float(rhs4[irow]))
    for S in reps5:
        Ss = set(S)
        for k in S:
            cv = np.zeros(NV)
            for r in range(n):
                if r in Ss:
                    continue
                cv[n4 + i6(*S, r)] += float(C[k, r])
            for x in S:
                if x == k:
                    continue
                cv[i4(*[y for y in S if y != x])] += float(C[k, x])
            cv[i4(*[y for y in S if y != k])] -= float(p)
            rows.append(cv)
            rhs.append(0.0)
    Aeq = np.vstack(rows)
    beq = np.array(rhs)
    _u, s, vt = np.linalg.svd(Aeq, full_matrices=True)
    rank = int(np.sum(s > 1e-8))
    ker_dim = NV - rank
    K = vt[rank:].T if ker_dim else np.zeros((NV, 0))
    part = np.linalg.lstsq(Aeq, beq, rcond=1e-8)[0]
    i, j, k, l = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    psi, fmul, fadd, fneg = make_psi(p)
    A = A_psi_matrix(p, C, psi, fadd, fneg)
    kapA = (
        A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]
    ).astype(np.float64)
    avg_A = np.bincount(loc4, weights=kapA, minlength=n4) / np.maximum(sz4, 1)
    w = np.zeros(NV)
    w[:n4] = sz4 * avg_A
    obj_c = float(w @ part)
    if ker_dim == 0:
        return {
            "p": p,
            "n4": n4,
            "n6": n6,
            "ker_dim": 0,
            "min_pairing": obj_c,
            "lp_ok": True,
        }
    c = K.T @ w
    A_ub = np.vstack([K, -K])
    b_ub = np.concatenate([1.0 - part, 1.0 + part])
    res_min = linprog(
        c, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * ker_dim, method="highs"
    )
    return {
        "p": p,
        "n4": n4,
        "n6": n6,
        "ker_dim": ker_dim,
        "rank": rank,
        "min_pairing": float(obj_c + res_min.fun) if res_min.success else None,
        "lp_ok": bool(res_min.success),
        "lstsq_res": float(np.max(np.abs(Aeq @ part - beq))),
    }


def sos4_edge_pairing_p5() -> dict:
    """Min QVAR pairing on the PSD-4 edge of the 15.590 deg-4 kernel at p=5."""
    lab = MuLab(5, with_deg6=False)
    sys4 = lab.equivariant_system(deg6=False)
    maps, kvs = sys4["maps"], sys4["kernel"]
    p, n, C = lab.p, lab.n, lab.C.astype(float)
    N = lab.N
    i, j, k, l = lab.S4.T
    psi, fmul, fadd, fneg = make_psi(5)
    A = A_psi_matrix(5, lab.C, psi, fadd, fneg)
    kapA = A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]

    def setfuncs(kv):
        fM = np.zeros(len(lab.S4))
        fD = np.zeros(len(lab.S4))
        for si in range(len(lab.S4)):
            lm = int(lab.labM4[si])
            if lm not in lab.deadM4:
                fM[si] = float(kv[maps["M4"][lm]]) * lab.sgM4[si]
            ld = int(lab.labD4[si])
            if ld not in lab.deadD4:
                fD[si] = float(kv[maps["D4"][ld]]) * lab.sgD4[si]
        return fM, fD

    KF = [setfuncs(kv) for kv in kvs]

    def T4_build(fvals, pair):
        T = np.zeros((n, n, n, n))
        idx = lab.S4
        for perm in itertools.permutations(range(4)):
            T[idx[:, perm[0]], idx[:, perm[1]], idx[:, perm[2]], idx[:, perm[3]]] = fvals
        for x in range(n):
            T[x, x, :, :] = pair
            T[x, :, x, :] = pair
            T[x, :, :, x] = pair
            T[:, x, x, :] = pair
            T[:, x, :, x] = pair
            T[:, :, x, x] = pair
        for x in range(n):
            T[x, x, :, :][np.arange(n), np.arange(n)] = 1.0
            T[x, :, x, :][np.arange(n), np.arange(n)] = 1.0
            T[x, :, :, x][np.arange(n), np.arange(n)] = 1.0
        for x in range(n):
            T[:, x, x, x] = pair[:, x]
            T[x, :, x, x] = pair[x, :]
            T[x, x, :, x] = pair[x, :]
            T[x, x, x, :] = pair[x, :]
        for x in range(n):
            T[x, x, x, x] = 1.0
        return T

    def compress(T, Q):
        m = Q.shape[1]
        X = np.einsum("ijkl,ia->ajkl", T, Q)
        X = np.einsum("ajkl,jb->abkl", X, Q)
        X = np.einsum("abkl,kc->abcl", X, Q)
        X = np.einsum("abcl,ld->abcd", X, Q)
        pr = [(a, b) for a in range(m) for b in range(a, m)]
        R = np.zeros((len(pr), len(pr)))
        for u, (a, b) in enumerate(pr):
            fu = np.sqrt(2) if a != b else 1.0
            for v, (c, d) in enumerate(pr):
                fv = np.sqrt(2) if c != d else 1.0
                R[u, v] = X[a, b, c, d] * fu * fv
        return R

    Proj = (p * np.eye(n) + C) / (2 * p)
    w, V = np.linalg.eigh(Proj)
    Qp = V[:, w > 0.5]
    w2, V2 = np.linalg.eigh((p * np.eye(n) - C) / (2 * p))
    Qm = V2[:, w2 > 0.5]
    R0p = compress(T4_build((lab.mu4 + lab.d4) / N, C / p), Qp)
    R0m = compress(T4_build((lab.mu4 - lab.d4) / N, -C / p), Qm)
    fM, fD = KF[0]
    RKp = compress(T4_build((fM + fD) / N, np.zeros((n, n))), Qp)
    RKm = compress(T4_build((fM - fD) / N, np.zeros((n, n))), Qm)
    corr = compress(T4_build(np.zeros(len(lab.S4)), np.zeros((n, n))), Qp)
    corrm = compress(T4_build(np.zeros(len(lab.S4)), np.zeros((n, n))), Qm)
    RKp = RKp - corr
    RKm = RKm - corrm

    def lmin(t):
        Mp = R0p + t * RKp
        Mm = R0m + t * RKm
        return min(float(np.linalg.eigvalsh(Mp)[0]), float(np.linalg.eigvalsh(Mm)[0]))

    def pairing(t):
        m = (lab.mu4 + lab.d4).astype(np.float64) / N + t * (fM + fD) / N
        return float(np.dot(m, kapA))

    edges = []
    for sign in (1.0, -1.0):
        lo, hi = 0.0, 1.0
        while lmin(hi * sign) > -1e-9 and hi < 1e8:
            hi *= 2
        for _ in range(80):
            mid = (lo + hi) / 2
            if lmin(mid * sign) > -1e-9:
                lo = mid
            else:
                hi = mid
        t = lo * sign
        edges.append({"t": t, "pairing": pairing(t), "lmin": lmin(t)})
    min_edge = min(edges, key=lambda e: e["pairing"])
    return {
        "kernel_dim": len(kvs),
        "true_pairing": pairing(0.0),
        "min_edge_pairing": min_edge["pairing"],
        "min_edge_t": min_edge["t"],
        "edges": edges,
        "claim_sos4_min_ge_0": bool(min_edge["pairing"] >= 0),
    }


def theorem_sixpoint_and_sos4_cannot_prove_qvar() -> dict:
    six = sixpoint_master_box_lp(5)
    sos = sos4_edge_pairing_p5()
    tru = true_maxplus_pairing(5)
    six_min = six["min_pairing"]
    ok = (
        six["lp_ok"]
        and six["ker_dim"] >= 1
        and six_min is not None
        and six_min < 0
        and sos["min_edge_pairing"] < 0
        and tru["pairing"] > 0
        and tru["clears_QVAR"]
        and sos["true_pairing"] > 0
        and abs(six_min - (-101 / 4)) < 1e-6
        and abs(sos["min_edge_pairing"] - (-45 / 4)) < 1e-6
    )
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "claim_sixpoint_min_ge_0": bool(six_min is not None and six_min >= 0),
        "claim_sos4_min_ge_0": sos["claim_sos4_min_ge_0"],
        "p5_sixpoint": six,
        "p5_sos4": sos,
        "p5_true": tru,
        "p7_joint_deg6_recorded": {
            "kernel_dim": 4,
            "min_pairing": -10633 / 8,
            "true_pairing_positive": True,
            "imported_as_p_law": False,
        },
        "theorem": (
            "Linear 6-point Aut-contractions + 4-point master + |m|≤1 "
            "have min pairing −101/4 at p=5 (ker dim 3).  SOS-4 along "
            "the 15.590 deg-4 kernel is feasible at pairing −45/4.  True "
            "Max+ pairing is positive.  Fail: claim either min is ≥0.  "
            "At p=7 the full 15.590 joint deg-4+6 system (ker dim 4) has "
            "box-LP min −10633/8 (recorded, not a p-law).  Linear 6-point "
            "theory does not prove QVAR."
        ),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    T = theorem_sixpoint_and_sos4_cannot_prove_qvar()
    path = ROOT / "evidence" / "e1_gmin_qvar_sixpoint.json"
    write_json_atomic(path, {"title": "6-point linear and SOS-4 cannot prove QVAR",
                             "numbered": False, "theorem": T})
    print("6-point/SOS-4 QVAR kill", flush=True)
    print(f"  proved_kill={T['proved']} inequality={T['inequality_proved']}", flush=True)
    print(f"  six min={T['p5_sixpoint']['min_pairing']} sos4={T['p5_sos4']['min_edge_pairing']}", flush=True)
    print("wrote", path, flush=True)
    return T


if __name__ == "__main__":
    main()
