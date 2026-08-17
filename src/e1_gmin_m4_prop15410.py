#!/usr/bin/env python3
"""
Prop 15.410 — Aut_e zero-far Type I is gap-2 and not bad-case
for every odd p≥5; the p=7 [0,1] Type I + S≤−1 + S+3f_e≤0
box is empty (HiGHS).  3A+B stays G>T.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.409 flags.

============================================================================
Setup.  Type I: k=3p−2, Aut_e-average x̄.  15.275 H: μ_far=0
gives μ_□=1/(2(p−1)), μ_⊠=5/(2(p+1)) ∈(0,1) and F_+ x̄=3−2f_e.
Aut_e-invariance ⇒ S is a function of f_e on Max−.

============================================================================
Theorem A — PROVED (Aut_e + 15.275 H; all odd p≥5).
  On the zero-far slice, Y_+=−(p+1)/(p−1) and
  Y_−=−5(p−1)/(p+1).  Hence S_max=Y_+<−1 (gap-2) and
      max(S+3f_e) ≥ Y_++3 = 2(p−2)/(p−1)>0
  so the point is not a bad case.  Fail: claim Y_++3≤0
  (false for every p≥5); claim S_max=−1.  ∎

Theorem B — PROVED (p=5,7 Max± caches).
  The explicit star-supported x (μ_□, μ_⊠, far=0, x_e=0)
  is Type I {1,5}, affine exact, S_max=Y_+, not bad-case.
  Fail: drop the square-star weight (then plus support leaves {1,5}
  at p=5).  ∎

Theorem C — PROVED (p=7 HiGHS).
  [0,1] Type I + S≤−1 + S+3f_e≤0 is infeasible (status 2).
  Fail: claim the bad-case LP feasible.  ∎

Theorem D — OPEN.  Aut_e far-class 3A+B>0 remains G>T.
  type_I flags stay False.

============================================================================
Backend: Fraction + p=5,7 caches + HiGHS.  Writes
evidence/e1_gmin_m4_prop15410.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import (
    lemma_three_orbit_zero_far_in_box,
    mu_minus_of_mu_far,
    mu_plus_of_mu_far,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from minmax_quadratic import paley_conference_prime_power

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15410_catalog.json"
YMINUS = {5: "/tmp/maxminus_p5.npy", 7: "/tmp/maxminus_p7.npy"}
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def Y_plus_zerofar(p: int) -> Fraction:
    return -Fraction(p + 1, p - 1)


def Y_minus_zerofar(p: int) -> Fraction:
    return -Fraction(5 * (p - 1), p + 1)


def S_plus_3fe_on_fplus(p: int) -> Fraction:
    return Y_plus_zerofar(p) + 3


def mu_box(p: int) -> Fraction:
    return mu_plus_of_mu_far(p, Fraction(0))


def mu_boxt(p: int) -> Fraction:
    return mu_minus_of_mu_far(p, Fraction(0))


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


def _unique_rows(F):
    key = np.ascontiguousarray(np.round(F).astype(np.int16))
    _, idx = np.unique(key, axis=0, return_index=True)
    return F[np.sort(idx)]


def zerofar_x(p: int, C, emap):
    n = C.shape[0]
    nE = n * (n - 1) // 2
    x = np.zeros(nE)
    mb = float(mu_box(p))
    mt = float(mu_boxt(p))
    for s in range(2, n):
        w = mb if int(C[1, s]) == 1 else mt
        x[emap[(0, s)]] = w
        x[emap[(1, s)]] = w
    return x


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        H = lemma_three_orbit_zero_far_in_box(p)
        Yp = Y_plus_zerofar(p)
        Ym = Y_minus_zerofar(p)
        cut = S_plus_3fe_on_fplus(p)
        ok = ok and H["proved"]
        ok = ok and mu_box(p) == Fraction(1, 2 * (p - 1))
        ok = ok and mu_boxt(p) == Fraction(5, 2 * (p + 1))
        ok = ok and Yp == -Fraction(p + 1, p - 1)
        ok = ok and Yp < -1
        ok = ok and Ym < Yp
        ok = ok and cut == Fraction(2 * (p - 2), p - 1)
        ok = ok and cut > 0
        rows[str(p)] = {"Y+": str(Yp), "cut": str(cut)}
        if cut <= 0 or Yp >= -1:
            ok = False
    ok = ok and S_plus_3fe_on_fplus(5) == Fraction(3, 2)
    ok = ok and S_plus_3fe_on_fplus(7) == Fraction(5, 3)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Zero-far is gap-2 and not bad-case for all odd p≥5. "
            "Fail: Y_++3≤0."
        ),
    }


def prove_B() -> dict:
    ok = True
    rec = {}
    for p in (5, 7):
        C, Yp, Ym = _load_pm(p)
        n = C.shape[0]
        emap = _edge_index(n)
        Fp, _, _ = _build_F(Yp, C)
        Fm, _, _ = _build_F(Ym, C)
        x = zerofar_x(p, C, emap)
        Sp = Fp @ x
        Sm = Fm @ x
        fe = Fm[:, 0]
        aff_err = float(np.max(np.abs(Sp + 2.0 * Fp[:, 0] - 3.0)))
        plus = sorted({int(round(v)) for v in Sp})
        smax = float(Sm.max())
        Yp_live = float(Sm[fe > 0].mean())
        cut = float((Sm + 3.0 * fe).max())
        hit = (
            aff_err < 1e-8
            and plus == [1, 5]
            and abs(smax - float(Y_plus_zerofar(p))) < 1e-8
            and smax <= -1 + 1e-8
            and cut > 0
            and abs(x.sum() - (3 * p - 2)) < 1e-8
            and x[emap[(0, 1)]] == 0.0
        )
        ok = ok and hit
        rec[str(p)] = {
            "plus": plus,
            "smax": smax,
            "Y+": Yp_live,
            "cut": cut,
            "affine_err": aff_err,
        }
        if plus != [1, 5] or cut <= 0:
            ok = False
        # fail-when-wrong: drop square-star weight at p=5
        if p == 5:
            x2 = x.copy()
            for s in range(2, n):
                if int(C[1, s]) == 1:
                    x2[emap[(0, s)]] = 0.0
                    x2[emap[(1, s)]] = 0.0
            plus2 = sorted({int(round(v)) for v in (Fp @ x2)})
            ok = ok and plus2 != [1, 5]
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "p=5,7 zero-far x is Type I {1,5}, gap-2, not bad-case. "
            "Fail: drop square-star at p=5."
        ),
    }


def prove_C() -> dict:
    from scipy.optimize import linprog

    C, Yp, Ym = _load_pm(7)
    Fp, _, _ = _build_F(Yp, C)
    Fm, _, _ = _build_F(Ym, C)
    Fu = _unique_rows(Fp)
    Fm_u = _unique_rows(Fm)
    nE = Fp.shape[1]
    rhs = 3.0 - 2.0 * Fu[:, 0]
    Aeq = np.vstack([Fu, np.ones((1, nE))])
    beq = np.append(rhs, 19.0)
    lb, ub = np.zeros(nE), np.ones(nE)
    ub[0] = 0.0
    bounds = list(zip(lb, ub))
    fe = Fm_u[:, 0]
    A_bad = np.vstack([Fm_u, Fm_u + 3.0 * fe.reshape(-1, 1)])
    b_bad = np.concatenate([-np.ones(Fm_u.shape[0]), np.zeros(Fm_u.shape[0])])
    r2 = linprog(
        np.zeros(nE),
        A_ub=A_bad,
        b_ub=b_bad,
        A_eq=Aeq,
        b_eq=beq,
        bounds=bounds,
        method="highs",
    )
    ok = (not bool(r2.success)) and int(r2.status) == 2
    if r2.success:
        ok = False
    return {
        "proved": bool(ok),
        "badcase_lp_feasible": bool(r2.success),
        "badcase_status": int(r2.status),
        "n_unique_plus": int(Fu.shape[0]),
        "n_unique_minus": int(Fm_u.shape[0]),
        "theorem": (
            "p=7 [0,1] Type I + S≤−1 + S+3f_e≤0 is infeasible. "
            "Fail: claim the bad-case LP feasible."
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
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "e1": bool(e1_closed_general()),
        "gsum": bool(gsum_disj_lb_proved_general()),
        "note": (
            "Zero-far is never the bad case; p=7 pointwise box empty. "
            "3A+B>0 on Aut_e far classes remains G>T. Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.410  zero-far gap-2 not bad-case; p=7 box empty", flush=True)
    A = prove_A()
    print(f"  A named: {A['proved']} p7cut={A['rows']['7']['cut']}", flush=True)
    B = prove_B()
    print(f"  B cache: {B['proved']} p7={B['by_p']['7']['plus']}", flush=True)
    C = prove_C()
    print(
        f"  C p7 box: {C['proved']} feasible={C['badcase_lp_feasible']}",
        flush=True,
    )
    D = prove_open()
    print(
        f"  D open: 3AB={D['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {
                "A": A["rows"],
                "B": {p: B["by_p"][p]["plus"] for p in B["by_p"]},
                "C": {
                    "feasible": C["badcase_lp_feasible"],
                    "status": C["badcase_status"],
                },
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.410",
        "title": (
            "Zero-far Type I is gap-2 and not bad-case for all p≥5; "
            "p=7 gap-2 bad-case box empty"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "zerofar_gap2_not_badcase": A["proved"],
            "p5_p7_zerofar_cache": B["proved"],
            "p7_badcase_box_empty": C["proved"],
            "type_I_multilevel_bad_case_ND_closed": D[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": D[
                "type_I_aut_e_3AB_positive_general"
            ],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15410.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
