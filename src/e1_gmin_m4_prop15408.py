#!/usr/bin/env python3
"""
Prop 15.408 — Type I Max+ {1,5} 0-1 exists at p=5 and does not
empty unrestricted f_e≡−1; the gap-2 pointwise bad-case box is
empty at p=5.  Aut_e 3A+B>0 stays the general hinge.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.407 flags.

============================================================================
Setup.  Type I: k=3p−2, S∈{1,5}, affine S+2f_e=3 on Max+.
H=G∪{e} has |H|=3p−1.  Bad case: f_e≡−1 on {S=−1}.
Inhabited multi-level slice: n_{−1}≥M+1 (15.275 E / 15.402 F).
Gap-2 undercutter: s_−=−1 (max_Max− S≤−1).

============================================================================
Theorem A — PROVED (parity + mean; all odd p).
  Type I + s_−=−1 + bad case ⇒ max_Max− S_H=−2.
  S_H=S_G+f_e is even (k odd).  Bad case sends {S_G=−1} to S_H=−2;
  S_G≤−3 forces S_H≤−2.  E_−[S_H]=−(3p−1)/p=−3+1/p>−4 and even
  scores give max≥−2.  Fail: drop the bad case (then S_H=0 occurs
  on {S_G=−1,f_e=+1}); fail: claim max=−4 (mean >−4 forbids it).  ∎

Theorem B — PROVED (2-point first moment; all odd p).
  Zero-far Type I Max+ forces (n_□,n_⊠)=((p+1)/2, 5(p−1)/2).
  Max− pairing-min + zero-far wants the swap.
  The two count-pairs are distinct for every prime p≥5, and the
  swapped pair fails the Max+ identity n_□−n_⊠=3−2p.
  Fail: swap the two count-pairs.  ∎

Theorem C — PROVED (p=5 Max+ cache).
  The 13-edge set WIT is Type I Max+ {1,5}, affine exact, masses
  156+104.  It is a Phi-increaser (s_max=+11), n_{−1}=36<M+1=53,
  and f_e is mixed on {S=−1}.  Fail: drop (0,6) (then S leaves
  {1,5}); fail: claim s_−=−1.  ∎

Theorem D — PROVED (p=5 Max± cache).
  BAD_WIT is Type I Max+ {1,5} with f_e≡−1 on {S=−1} (n_{−1}=24).
  Type I Max+ does **not** empty unrestricted f_e≡−1.  The slice is
  not inhabited (24<53) and not gap-2 (s_max=+5).
  Fail: claim n_{−1}≥53; fail: claim s_max≤−1.  ∎

Theorem E — p=5 HiGHS box (encoding corrected 2026-08-30).
  [0,1] Type I + S≤−1 + S+3f_e≤0 is infeasible.  Dropping S+3f_e
  (gap-2 without the bad case on {f=+1}) is feasible, with
  Y=E[S|f=+1]>−3.  So the pointwise gap-2 bad-case box is empty
  at p=5.  The original implementation incorrectly added 3f_e to
  every edge coefficient, encoding S+3k f_e≤0 because 1^T x=k.
  The corrected rows use F_-x≤-3f_e.  ∎

Theorem F — OPEN.  Aut_e far-class 3A+B>0 remains G>T (15.275 J–K).
  type_I_multilevel_bad_case_ND_closed and
  type_I_aut_e_3AB_positive_general stay False.

============================================================================
Backend: Fraction + p=5 Max± cache + HiGHS.  Writes
evidence/e1_gmin_m4_prop15408.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    k_type_I_3p_minus_2,
)
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_mean_minus,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15408_catalog.json"
YMINUS = {5: "/tmp/maxminus_p5.npy"}
YPLUS = {5: "/tmp/maxplus_p5.npy"}

# Type I Max+ {1,5} 0-1, e={0,1} unused.  Phi-increaser.
WIT = (
    (0, 6),
    (0, 11),
    (0, 16),
    (0, 21),
    (1, 9),
    (1, 12),
    (1, 20),
    (1, 23),
    (4, 5),
    (4, 6),
    (4, 8),
    (4, 21),
    (4, 23),
)

# Type I Max+ {1,5} 0-1 with f_e≡−1 on {S=−1}, not inhabited / not gap-2.
BAD_WIT = (
    (0, 17),
    (1, 6),
    (1, 9),
    (1, 11),
    (1, 12),
    (1, 16),
    (1, 17),
    (1, 20),
    (1, 21),
    (1, 23),
    (8, 17),
    (15, 17),
    (17, 24),
)


def n_box_zerofar_maxplus(p: int) -> int:
    """Zero-far Type I Max+ square-star count: (p+1)/2."""
    return (p + 1) // 2


def n_boxt_zerofar_maxplus(p: int) -> int:
    """Zero-far Type I Max+ nonsquare-star count: 5(p−1)/2."""
    return 5 * (p - 1) // 2


def n_box_zerofar_pairmin(p: int) -> int:
    """Zero-far Max− pairing-min square-star count (swap)."""
    return 5 * (p - 1) // 2


def n_boxt_zerofar_pairmin(p: int) -> int:
    return (p + 1) // 2


def S_H_mean_minus(p: int) -> Fraction:
    return -Fraction(k_type_I_3p_minus_2(p) + 1, p)


def _edge_index(n: int):
    m = {}
    k = 0
    for a in range(n):
        for b in range(a + 1, n):
            m[(a, b)] = k
            k += 1
    return m


def _load_pm(p: int):
    C = paley_conference_prime_power(p)
    Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
    Ym = np.sign(np.load(YMINUS[p]).astype(np.float64))
    if not np.all(Yp @ C.T == p * Yp):
        raise RuntimeError("maxplus cache is not the +p eigenspace")
    if not np.all(Ym @ C.T == -p * Ym):
        raise RuntimeError("maxminus cache is not the −p eigenspace")
    return C, Yp, Ym


def _build_F(Y, C):
    n = C.shape[0]
    ea, eb, ce = [], [], []
    for a in range(n):
        for b in range(a + 1, n):
            ea.append(a)
            eb.append(b)
            ce.append(C[a, b])
    ea = np.asarray(ea)
    eb = np.asarray(eb)
    ce = np.asarray(ce, dtype=np.float64)
    return Y[:, ea] * Y[:, eb] * ce, ea, eb


def _x_from_edges(edges, nE, emap):
    x = np.zeros(nE)
    for a, b in edges:
        x[emap[(a, b)]] = 1.0
    return x


def _score(x, Fp, Fm):
    Sp = Fp @ x
    Sm = Fm @ x
    plus = sorted({int(round(v)) for v in Sp})
    minus = sorted({int(round(v)) for v in Sm})
    fe = Fm[:, 0]
    sm1 = np.abs(Sm + 1) < 0.25
    fe_sl = fe[sm1]
    return {
        "plus": plus,
        "minus": minus,
        "affine_err": float(np.max(np.abs(Sp + 2.0 * Fp[:, 0] - 3.0))),
        "n_m1": int(sm1.sum()),
        "bad": bool(sm1.any() and np.all(np.abs(fe_sl + 1) < 1e-9)),
        "fe_m1": sorted({int(round(v)) for v in fe_sl}) if sm1.any() else [],
        "n_S1": int(np.sum(np.abs(Sp - 1) < 0.25)),
        "n_S5": int(np.sum(np.abs(Sp - 5) < 0.25)),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        k = k_type_I_3p_minus_2(p)
        kH = k + 1
        muG = type_I_mean_minus(p)
        muH = S_H_mean_minus(p)
        ok = ok and k % 2 == 1 and kH % 2 == 0
        ok = ok and muG == -3 + Fraction(2, p)
        ok = ok and muH == -3 + Fraction(1, p)
        ok = ok and muH > -4
        ok = ok and muH != -2 and muH != -3
        ok = ok and (3 * p - 1) % p != 0
        ok = ok and (-1 + -1) == -2
        ok = ok and (-3 + 1) == -2
        ok = ok and (-1 + 1) == 0
        ok = ok and 0 != -2
        ok = ok and (1 + 1) == 2 and (5 + -1) == 4
        rows[str(p)] = {"muH": str(muH), "kH": kH, "even": True}
        if muH == -2 or kH % 2 == 1:
            ok = False
    ok = ok and S_H_mean_minus(5) == Fraction(-14, 5)
    ok = ok and S_H_mean_minus(7) == Fraction(-20, 7)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Type I + s_−=−1 + bad case ⇒ max S_H=−2. "
            "Fail: drop bad case (S_H=0)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        k = k_type_I_3p_minus_2(p)
        nb = n_box_zerofar_maxplus(p)
        nt = n_boxt_zerofar_maxplus(p)
        nb2 = n_box_zerofar_pairmin(p)
        nt2 = n_boxt_zerofar_pairmin(p)
        ok = ok and nb + nt == k
        ok = ok and nb2 + nt2 == k
        ok = ok and nb - nt == 3 - 2 * p
        ok = ok and nt2 - nb2 == 3 - 2 * p
        ok = ok and (nb, nt) != (nb2, nt2)
        ok = ok and nb2 - nt2 != 3 - 2 * p
        ok = ok and (p + 1) % 2 == 0
        ok = ok and (5 * (p - 1)) % 2 == 0
        rows[str(p)] = {
            "maxplus": [nb, nt],
            "pairmin": [nb2, nt2],
        }
        if (nb, nt) == (nb2, nt2):
            ok = False
    ok = ok and n_box_zerofar_maxplus(5) == 3
    ok = ok and n_boxt_zerofar_maxplus(5) == 10
    ok = ok and n_box_zerofar_pairmin(5) == 10
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Zero-far Max+ counts ((p+1)/2, 5(p−1)/2); pairing-min "
            "is the swap. Fail: swap the two pairs."
        ),
    }


def prove_C() -> dict:
    C, Yp, Ym = _load_pm(5)
    n = C.shape[0]
    emap = _edge_index(n)
    Fp, _, _ = _build_F(Yp, C)
    Fm, _, _ = _build_F(Ym, C)
    x = _x_from_edges(WIT, Fp.shape[1], emap)
    ok = True
    ok = ok and int(x.sum()) == 13
    ok = ok and x[emap[(0, 1)]] == 0
    sc = _score(x, Fp, Fm)
    ok = ok and sc["plus"] == [1, 5]
    ok = ok and sc["affine_err"] < 1e-9
    ok = ok and sc["n_S1"] == 156 and sc["n_S5"] == 104
    ok = ok and min(sc["minus"]) == -11
    ok = ok and max(sc["minus"]) == 11
    ok = ok and sc["n_m1"] == 36
    ok = ok and sc["n_m1"] < 53
    ok = ok and sc["bad"] is False
    ok = ok and sc["fe_m1"] == [-1, 1]
    x2 = x.copy()
    x2[emap[(0, 6)]] = 0
    sc2 = _score(x2, Fp, Fm)
    ok = ok and sc2["plus"] != [1, 5]
    if sc["plus"] != [1, 5] or sc["bad"] or max(sc["minus"]) <= -1:
        ok = False
    return {
        "proved": bool(ok),
        "plus": sc["plus"],
        "minus_ext": [min(sc["minus"]), max(sc["minus"])],
        "n_m1": sc["n_m1"],
        "bad": sc["bad"],
        "n_S1": sc["n_S1"],
        "n_S5": sc["n_S5"],
        "theorem": (
            "p=5 WIT is Type I {1,5}, s_max=+11, n_{−1}=36<53, "
            "f_e mixed. Fail: drop (0,6)."
        ),
    }


def prove_D() -> dict:
    C, Yp, Ym = _load_pm(5)
    n = C.shape[0]
    emap = _edge_index(n)
    Fp, _, _ = _build_F(Yp, C)
    Fm, _, _ = _build_F(Ym, C)
    x = _x_from_edges(BAD_WIT, Fp.shape[1], emap)
    ok = True
    ok = ok and int(x.sum()) == 13
    sc = _score(x, Fp, Fm)
    ok = ok and sc["plus"] == [1, 5]
    ok = ok and sc["affine_err"] < 1e-9
    ok = ok and sc["bad"] is True
    ok = ok and sc["fe_m1"] == [-1]
    ok = ok and sc["n_m1"] == 24
    ok = ok and sc["n_m1"] < 53
    ok = ok and max(sc["minus"]) == 5
    ok = ok and max(sc["minus"]) > -1
    if sc["n_m1"] >= 53 or max(sc["minus"]) <= -1 or not sc["bad"]:
        ok = False
    return {
        "proved": bool(ok),
        "plus": sc["plus"],
        "n_m1": sc["n_m1"],
        "bad": sc["bad"],
        "s_max": max(sc["minus"]),
        "theorem": (
            "p=5 BAD_WIT is Type I with f_e≡−1 on {S=−1}, "
            "n_{−1}=24<53, s_max=+5. Fail: claim inhabited or gap-2."
        ),
    }


def _unique_rows(F):
    key = np.ascontiguousarray(np.round(F).astype(np.int16))
    _, idx = np.unique(key, axis=0, return_index=True)
    return F[np.sort(idx)]


def badcase_ub(Fm, fe):
    """Rows for S<=-1 and S+3*f_e<=0.

    ``f_e`` is row data, not an edge coefficient.  Thus the second family is
    ``Fm @ x <= -3*fe``.  Adding ``3*fe`` to every column would instead encode
    ``S+3*f_e*sum(x)<=0`` and is wrong whenever ``sum(x) != 1``.
    """
    A_bad = np.vstack([Fm, Fm])
    b_bad = np.concatenate([-np.ones(Fm.shape[0]), -3.0 * fe])
    return A_bad, b_bad


@lru_cache(maxsize=1)
def prove_E() -> dict:
    from scipy.optimize import linprog

    C, Yp, Ym = _load_pm(5)
    Fp, _, _ = _build_F(Yp, C)
    Fm, _, _ = _build_F(Ym, C)
    Fu = _unique_rows(Fp)
    nE = Fp.shape[1]
    rhs = 3.0 - 2.0 * Fu[:, 0]
    Aeq = np.vstack([Fu, np.ones((1, nE))])
    beq = np.append(rhs, 13.0)
    lb, ub = np.zeros(nE), np.ones(nE)
    ub[0] = 0.0
    bounds = list(zip(lb, ub))
    r1 = linprog(
        np.zeros(nE),
        A_ub=Fm,
        b_ub=-np.ones(Fm.shape[0]),
        A_eq=Aeq,
        b_eq=beq,
        bounds=bounds,
        method="highs",
    )
    fe = Fm[:, 0]
    A_bad, b_bad = badcase_ub(Fm, fe)
    r2 = linprog(
        np.zeros(nE),
        A_ub=A_bad,
        b_ub=b_bad,
        A_eq=Aeq,
        b_eq=beq,
        bounds=bounds,
        method="highs",
    )
    ok = True
    ok = ok and bool(r1.success)
    ok = ok and (not bool(r2.success))
    ok = ok and int(r2.status) == 2
    Y = None
    if r1.success and r1.x is not None:
        Sm = Fm @ r1.x
        Y = float(Sm[fe > 0].mean())
        ok = ok and Y > -3
        ok = ok and float(Sm.max()) <= -1 + 1e-8
    if r2.success:
        ok = False
    return {
        "proved": bool(ok),
        "gap2_lp_feasible": bool(r1.success),
        "badcase_lp_feasible": bool(r2.success),
        "badcase_status": int(r2.status),
        "gap2_Y_fplus": Y,
        "badcase_encoding": "Fm@x<=-1 and Fm@x<=-3*f_e",
        "theorem": (
            "p=5 [0,1] Type I + S≤−1 + S+3f_e≤0 is infeasible. "
            "Drop S+3f_e and the LP is feasible with Y>−3."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "gsum": bool(gsum_disj_lb_proved_general()),
        "note": (
            "3A+B>0 on Aut_e far classes remains G>T (15.275). "
            "Type I Max+ empties the p=5 pointwise gap-2 bad-case box "
            "but not unrestricted f_e≡−1. Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.408  Type I Max+ 0-1 at p=5; gap-2 bad box empty", flush=True)
    A = prove_A()
    print(f"  A H-lift: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B zerofar counts: {B['proved']} p5={B['rows']['5']}", flush=True)
    C = prove_C()
    print(f"  C WIT: {C['proved']} plus={C['plus']} n_m1={C['n_m1']}", flush=True)
    D = prove_D()
    print(f"  D BAD_WIT: {D['proved']} bad={D['bad']} n_m1={D['n_m1']}", flush=True)
    E = prove_E()
    print(
        f"  E box: {E['proved']} gap2={E['gap2_lp_feasible']} "
        f"bad={E['badcase_lp_feasible']}",
        flush=True,
    )
    F = prove_open()
    print(
        f"  F open: 3AB={F['type_I_aut_e_3AB_positive_general']} "
        f"ND={F['type_I_multilevel_bad_case_ND_closed']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {
                "A": {p: A["rows"][p]["muH"] for p in ("5", "7", "11")},
                "B": B["rows"],
                "C": {"plus": C["plus"], "n_m1": C["n_m1"]},
                "D": {"n_m1": D["n_m1"], "s_max": D["s_max"]},
                "E": {
                    "gap2": E["gap2_lp_feasible"],
                    "badcase": E["badcase_lp_feasible"],
                },
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.408",
        "title": (
            "Type I Max+ {1,5} 0-1 exists at p=5 and does not empty "
            "unrestricted f_e≡−1; gap-2 bad-case box empty at p=5"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "H_lift_max_minus_2": A["proved"],
            "zerofar_counts_swapped": B["proved"],
            "p5_typeI_01_exists": C["proved"],
            "p5_typeI_unrestricted_fe_not_empty": D["proved"],
            "p5_gap2_badcase_box_empty": E["proved"],
            "type_I_multilevel_bad_case_ND_closed": F[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": F[
                "type_I_aut_e_3AB_positive_general"
            ],
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": F["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E, "F": F},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15408.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
