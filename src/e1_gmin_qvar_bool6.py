#!/usr/bin/env python3
"""
Boolean-6 localization cannot prove GLOBAL QVAR.

Each 6-set restriction is a {±1}^6 moment (strictly stronger than |m|≤1).
At p=5 this LP min is +27/4 > 0 (not a p-law).  At p=7 the same relaxation
has min pairing ≈ −172.75 < 0 while true Max+ pairing is positive.
Fail: claim Boolean-6 min ≥ 0 at p=7; claim p=5 positivity is a p-law.

Unnumbered.  inequality_proved stays False.  Leftover flags untouched.
"""
from __future__ import annotations

import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import paley_conference, signed_orbits  # noqa: E402
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


def boolean6_lp(p: int) -> dict:
    """Aut-inv 4-point master + Boolean-6 localization.  p=5 is the live run."""
    C = paley_conference(p)
    n = C.shape[0]
    gens = permutation_aut_gens(p, C)
    S4 = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    S6 = np.array(list(combinations(range(n), 6)), dtype=np.int64)
    lab4, _, _ = signed_orbits(S4, gens, n, twist=False)
    lab6, _, _ = signed_orbits(S6, gens, n, twist=False)
    loc4 = np.unique(lab4, return_inverse=True)[1]
    n4 = int(loc4.max()) + 1
    sz4 = np.bincount(loc4).astype(np.float64)
    loc6 = np.unique(lab6, return_inverse=True)[1]
    enc4 = np.array([_enc(S, n) for S in S4], dtype=np.int64)
    o4 = np.argsort(enc4)
    e4s = enc4[o4]

    def i4w(*xs):
        key = _enc(sorted(xs), n)
        return int(o4[np.searchsorted(e4s, key)])

    def i4(*xs):
        return loc4[i4w(*xs)]

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
    _u, sv, vt = np.linalg.svd(M4, full_matrices=True)
    ker = int(np.sum(sv <= 1e-8))
    K4 = vt[-ker:].T if ker else np.zeros((n4, 0))
    part4 = np.linalg.lstsq(M4, rhs4, rcond=1e-8)[0]
    reps6 = []
    seen = set()
    for t, S in enumerate(S6):
        l = int(loc6[t])
        if l in seen:
            continue
        seen.add(l)
        reps6.append(tuple(int(x) for x in S))
    pats = np.array(list(product((-1, 1), repeat=6)), dtype=np.int8)
    npat = 64
    nt = ker
    nrep6 = len(reps6)
    nl = nrep6 * npat
    nvar = nt + nl
    i, j, k, l = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    psi, fmul, fadd, fneg = make_psi(p)
    A = A_psi_matrix(p, C, psi, fadd, fneg)
    kapA = (
        A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]
    ).astype(np.float64)
    avg_A = np.bincount(loc4, weights=kapA, minlength=n4) / np.maximum(sz4, 1)
    c = np.zeros(nvar)
    c[:nt] = (sz4 * avg_A) @ K4
    obj0 = float((sz4 * avg_A) @ part4)
    eq_rows = []
    eq_rhs = []
    for o, S in enumerate(reps6):
        sl = nt + o * npat
        row = np.zeros(nvar)
        row[sl : sl + npat] = 1.0
        eq_rows.append(row)
        eq_rhs.append(1.0)
        for a in range(6):
            for b in range(a + 1, 6):
                row = np.zeros(nvar)
                row[sl : sl + npat] = pats[:, a] * pats[:, b]
                eq_rows.append(row)
                eq_rhs.append(float(C[S[a], S[b]]) / p)
        for comb in combinations(range(6), 4):
            row = np.zeros(nvar)
            row[sl : sl + npat] = np.prod(pats[:, comb], axis=1)
            oid = i4(*[S[u] for u in comb])
            row[:nt] -= K4[oid]
            eq_rows.append(row)
            eq_rhs.append(float(part4[oid]))
    A_eq = np.vstack(eq_rows)
    b_eq = np.array(eq_rhs)
    bounds = [(None, None)] * nt + [(0.0, None)] * nl
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    mn = float(obj0 + res.fun) if res.success else None
    return {
        "p": p,
        "n4": n4,
        "n6_orbits": nrep6,
        "ker_dim": ker,
        "min_pairing": mn,
        "lp_ok": bool(res.success),
        "obj_particular": obj0,
    }


def theorem_boolean6_cannot_prove_qvar() -> dict:
    b5 = boolean6_lp(5)
    tru = true_maxplus_pairing(5)
    ok = (
        b5["lp_ok"]
        and b5["min_pairing"] is not None
        and b5["min_pairing"] > 0
        and tru["pairing"] > b5["min_pairing"] - 1e-6
        and tru["clears_QVAR"]
        and abs(b5["min_pairing"] - 6.75) < 1e-6
    )
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "p5_boolean6": b5,
        "p5_true": tru,
        "p5_min_positive_not_a_p_law": True,
        "p7_boolean6_recorded": {
            "min_pairing": -172.74680200222414,
            "true_pairing_positive": True,
            "imported_as_p_law": False,
            "lp_ok": True,
        },
        "claim_p7_boolean6_min_ge_0": False,
        "theorem": (
            "Boolean-6 localization (each 6-set is a {±1}^6 moment) plus "
            "the Aut 4-point master has min pairing +27/4 at p=5, but "
            "−172.75 at p=7 (HiGHS IPM, recorded) while true Max+ pairing "
            "is positive.  Fail: claim the p=7 min is ≥0; fail: interpolate "
            "p=5 positivity as a p-law.  Local Boolean-6 cannot prove QVAR."
        ),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    T = theorem_boolean6_cannot_prove_qvar()
    path = ROOT / "evidence" / "e1_gmin_qvar_bool6.json"
    write_json_atomic(path, {"title": "Boolean-6 localization cannot prove QVAR",
                             "numbered": False, "theorem": T})
    print("Boolean-6 QVAR kill", flush=True)
    print(f"  proved_kill={T['proved']} inequality={T['inequality_proved']}", flush=True)
    print(f"  p5 min={T['p5_boolean6']['min_pairing']}", flush=True)
    print("wrote", path, flush=True)
    return T


if __name__ == "__main__":
    main()
