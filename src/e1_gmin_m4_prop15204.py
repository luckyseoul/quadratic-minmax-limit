#!/usr/bin/env python3
"""
Prop 15.204 — Ker(Gsum) PSD characterization; free-e dual over scheme⊕cross;
cost_D·p census; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182, 15.190–203)
  Dual-eq empty ⇔ free-e max κ_e < 2−α on ker(Gsum)∩box, α=6/(p n).
  scheme⊕cross ⊆ ker (15.202).  D(C) dual cost_D bounds free-e over sc
  (15.203).  Open: ker=sc general and/or cost_D<2−α general, or |μ|/K₄.

PROVED Max+-free
  A. PSD ker characterization.  Gsum = E₊[f fᵀ]+E₋[f fᵀ] ≽ 0 on edge
     space, with f_e(y)=C_e y_i y_j.  For any edge-weight κ,
        κᵀ Gsum κ = E₊[q_κ²] + E₋[q_κ²],   q_κ(y) := ∑_e κ_e f_e(y)
                  = ½ yᵀ (κ ⊙ C) y   (zero-diag convention).
     Hence κ ∈ ker(Gsum) ⇔ κᵀ Gsum κ = 0 ⇔ q_κ ≡ 0 on Max+ ∪ Max−.  ∎

  B. Consequence for free-e.  ker(Gsum) = {κ : yᵀ(κ⊙C)y = 0 ∀ y∈Max±}.
     scheme⊕cross ⊆ ker ⇒ free-e max on true ker ≤ free-e max on sc,
     with equality whenever ker = scheme⊕cross.  ∎

  C. Free-e dual over sc (lower box).  As in 15.202, when only the lower
     box is active at optimum,
        max_{κ∈sc∩box} κ_e = α · S*_sc,
     where
        S*_sc = min{ ∑_{g≠e} w_g : w_e=1, w≥0, w ⊥ scheme⊕cross }.
     Equivalently L_scᵀ w = 0 with W_e=1, W≥0 (edge weights).  ∎

  D. Sufficient dual sum bound.  For all primes p≥5:
        8/p < 2 − α = 2 − 6/(p n).
     Proof: 8/p + 6/(p n) < 2 ⇔ (8n+6)/(p n) < 2 ⇔ 8n+6 < 2 p n
        ⇔ 8(p²+1)+6 < 2p(p²+1) ⇔ 0 < 2p³ − 8p² − 2p + 14
        ⇔ 0 < p³ − 4p² − p + 7.  At p=5 value 27>0; derivative
        3p²−8p−1 > 0 on [5,∞).  ∎
     Therefore any feasible sc-dual with sum_ne ≤ (4/3) n gives
        free-e max_sc ≤ α·(4/3)n = 8/p < 2−α,
     hence blocks dual-eq whenever ker=sc.  ∎

  E. Recall (15.201): α(n−2)<2−α for p≥5, so S*_sc ≤ n−2 also blocks
     dual-eq when ker=sc.  ∎

CERTIFIED (not general)
  F. Dual LP over sc: S*_sc = 123/7 (p=5), 3912/113 (p=7); free-e max
     = α S* matches 369/455, 11736/19775 (15.202).  Both ≤ α(n−2) and
     ≤ (3/2)·scheme, so dual-eq empty at p=5,7 if ker=sc (which holds
     by dim match at those p).
  G. cost_D(p)·p ∈ [6.5, 7.7] for primes p∈{5,7,11,13,17,19}; all have
     cost_D < 2−α and sum_ne ≤ (4/3)n (p=7 closest: 64.03 ≤ 66.67).
  H. Quadratic rigidity on V+: for M∈Sym(V+), uᵀ M u = 0 on all Max+
     ⇔ M is the ++ block of a scheme matrix (dim ker Q₊ = n−1 =
     rank scheme→Sym₊₊), residual outside scheme <1e−14 at p=3,5.
     Combined with (A)–(B) this reconfirms ker=sc at p=3,5 independently
     of the dim(Gsum) count in 15.202.

OPEN (blocks residual i / predicate flip)
  I. Prove for all primes p≥5 one of:
     (i)   ker(Gsum)=scheme⊕cross (e.g. extend H: ker Q₊ = scheme_image
           for all p, plus matching − side / single f),
     (ii)  S*_sc ≤ n−2 (or ≤ 3(n−2)/4), Max+-free closed form,
     (iii) cost_D(p) ≤ 8/p (or sum_ne(D(C)) ≤ (4/3)n) for all p≥5,
           together with (i),
     (iv)  |μ₄|≤2/n on |κ|=1, or K₄≤Wick_hi.
     Then flip gsum_disj_lb_proved_general → residual_i → E1 → L only if
     residual (ii) full allows.  Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15204.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15203 import dual_construction_cost  # noqa: E402


def theorem_psd_ker_characterization() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Gsum ≽ 0 and κᵀ Gsum κ = E₊[q_κ²]+E₋[q_κ²] with "
            "q_κ(y)=½ yᵀ(κ⊙C)y.  Hence ker(Gsum)={κ : q_κ≡0 on Max±}."
        ),
    }


def theorem_free_e_over_sc_bounds_true() -> dict:
    return {
        "proved": True,
        "theorem": (
            "scheme⊕cross ⊆ ker(Gsum) ⇒ free-e max on ker ≤ free-e max on sc, "
            "with equality when ker=sc.  Dual-eq blocked if free-e max_sc < 2−α "
            "and ker=sc."
        ),
    }


def theorem_eight_over_p_beats_need() -> dict:
    """8/p < 2−α for all primes p≥5 (Fraction algebra)."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        al = alpha_particular(p)
        lhs = Fraction(8, p)
        need = 2 - al
        # cubic p³ − 4p² − p + 7
        poly = p**3 - 4 * p**2 - p + 7
        beats = lhs < need
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "8_over_p": str(lhs),
                "need": str(need),
                "beats": beats,
                "poly": poly,
            }
        if not beats or poly <= 0:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: 8/p < 2−α.  Hence any sc-dual with "
            "sum_ne ≤ (4/3)n yields free-e max_sc ≤ 8/p < 2−α."
        ),
        "by_p_sample": sample,
    }


def free_e_dual_Sstar_sc(p: int) -> dict:
    """
    Solve S*_sc = min{∑ w − w_e : w_e=1, w≥0, L_scᵀ w = 0} by linprog.
    Max+-free once ker=sc is known; here used as census for p=5,7.
    """
    import numpy as np
    from scipy.optimize import linprog
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    d = n // 2
    evals, evecs = np.linalg.eigh(C)
    order = np.argsort(evals)
    Vm = evecs[:, order[:d]]
    Vp = evecs[:, order[d:]]
    ii, jj = np.triu_indices(n, k=1)
    m = len(ii)
    n_scheme = n - 1
    n_cross = d * d
    nparam = n_scheme + n_cross
    L = np.zeros((m, nparam))
    for v in range(n - 1):
        f = np.zeros(n)
        f[v] = 1.0
        f[n - 1] = -1.0
        L[:, v] = f[ii] + f[jj]
    for a in range(d):
        for b in range(d):
            k = n_scheme + a * d + b
            B = np.zeros((d, d))
            B[a, b] = 1.0
            A = Vp @ B @ Vm.T + Vm @ B.T @ Vp.T
            L[:, k] = C[ii, jj] * A[ii, jj]
    # min 1ᵀw s.t. Lᵀ w = 0, w_0 = 1, w ≥ 0
    c = np.ones(m)
    A_eq = np.vstack([L.T, np.eye(1, m)])
    b_eq = np.zeros(nparam + 1)
    b_eq[-1] = 1.0
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    if not res.success:
        return {"p": p, "ok": False, "message": res.message}
    Sstar = float(res.fun - 1.0)
    al = alpha_particular(p)
    free_e_ub = al * Fraction(Sstar).limit_denominator(10**9)
    # match known exact values when close
    known = {5: Fraction(123, 7), 7: Fraction(3912, 113)}
    S_exact = known.get(p)
    if S_exact is not None and abs(Sstar - float(S_exact)) < 1e-8:
        Sstar_out = str(S_exact)
        free_e_ub = al * S_exact
    else:
        Sstar_out = str(Fraction(Sstar).limit_denominator(10**6))
    need = 2 - al
    return {
        "p": p,
        "ok": True,
        "Sstar": Sstar_out,
        "Sstar_float": Sstar,
        "free_e_ub": str(free_e_ub),
        "need": str(need),
        "blocks_if_ker_sc": float(free_e_ub) < float(need) - 1e-12,
        "n_minus_2": n - 2,
        "Sstar_le_n_minus_2": Sstar <= (n - 2) + 1e-9,
        "proved_general": False,
    }


def certify_cost_D_times_p(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19]
    rows = {}
    all_ok = True
    for p in primes:
        d = dual_construction_cost(p)
        n = n_of(p)
        sn = d["sum_ne"]
        cost = d["cost"]
        need = d["need"]
        row = {
            "cost": cost,
            "need": need,
            "cost_times_p": cost * p,
            "sum_ne": sn,
            "four_thirds_n": 4 * n / 3,
            "sum_ne_le_4_3_n": sn <= 4 * n / 3 + 1e-6,
            "cost_le_8_over_p": cost <= 8 / p + 1e-6,
            "blocks_if_ker_sc": d["blocks_if_ker_eq_sc"],
        }
        rows[str(p)] = row
        if not (row["blocks_if_ker_sc"] and row["sum_ne_le_4_3_n"]):
            all_ok = False
    return {
        "certified_range": primes,
        "all_block_if_ker_sc_and_sum_ok": all_ok,
        "by_p": rows,
        "proved_general": False,
        "note": (
            "cost_D·p ≈ 6.5–7.7 on this range; sum_ne≤(4/3)n holds. "
            "General cost_D≤8/p still open."
        ),
    }


def free_e_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_204() -> bool:
    return free_e_bound_proved_general()


def hinge_status_204() -> dict:
    return {
        "psd_ker_characterization": True,
        "free_e_sc_bounds_true_ker": True,
        "eight_over_p_beats_need": True,
        "ker_eq_sc_general": False,
        "Sstar_le_n_minus_2_general": False,
        "cost_D_le_8_over_p_general": False,
        "residual_i_closed": residual_i_closed_via_204(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove ker=sc for all p≥5 and (S*_sc≤n−2 or cost_D≤8/p or "
            "sum_ne≤(4/3)n), or |μ|≤2/n / K₄≤Wick_hi Max+-free."
        ),
    }


def main() -> dict:
    psd = theorem_psd_ker_characterization()
    scb = theorem_free_e_over_sc_bounds_true()
    eight = theorem_eight_over_p_beats_need()
    duals = {str(p): free_e_dual_Sstar_sc(p) for p in (5, 7)}
    costc = certify_cost_D_times_p()
    out = {
        "title": (
            "Prop 15.204 ker PSD characterization; free-e dual over sc; "
            "residual (i) OPEN"
        ),
        "proved": {
            "psd_ker_characterization": psd["proved"],
            "free_e_sc_bounds_true_ker": scb["proved"],
            "eight_over_p_beats_need": eight["proved"],
            "ker_eq_sc_general": False,
            "Sstar_le_n_minus_2_general": False,
            "cost_D_le_8_over_p_general": False,
            "free_e_bound_general": False,
            "residual_i_closed": residual_i_closed_via_204(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "psd_ker": psd,
            "free_e_sc": scb,
            "eight_over_p": eight,
        },
        "certified": {
            "Sstar_sc_dual": duals,
            "cost_D_times_p": costc,
        },
        "hinge": hinge_status_204(),
        "F3": (
            "PSD ker form and 8/p algebra are Max+-free, but ker=sc and "
            "S*/cost_D general bounds remain open.  Do not flip residual/E1/L."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15204.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.204 ker PSD characterization; free-e dual over sc")
    print(f"  psd ker: {psd['proved']}")
    print(f"  8/p < 2-α p>=5: {eight['proved']}")
    for pk, v in duals.items():
        if not v.get("ok"):
            print(f"  p={pk}: dual fail {v}")
            continue
        print(
            f"  p={pk}: S*={v['Sstar']} free-e ub={v['free_e_ub']} "
            f"need={v['need']} blocks_if_ker_sc={v['blocks_if_ker_sc']} "
            f"S*<=n-2={v['Sstar_le_n_minus_2']}"
        )
    print(
        f"  cost_D census all_ok={costc['all_block_if_ker_sc_and_sum_ok']} "
        f"range={costc['certified_range']}"
    )
    print(f"  residual_i closed: {residual_i_closed_via_204()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
