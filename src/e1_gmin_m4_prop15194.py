#!/usr/bin/env python3
"""
Prop 15.194 — Row negative-mass obstruction for dual-eq ker-box (residual i).

Does **not** flip residual_i / gsum_disj_lb / E1 / L. Soft-close forbidden.

SETUP (15.182 dual-eq normal form + 15.192 Gsum row algebra)
  Gsum_ee = 2, adjacent (wedge) entries = 0, row sum = n.
  Hence for disj partners f ⊥ e:
        ∑_{f ⊥ e} a_f = n − 2,   a_f := Gsum_ef.
  Dual-eq: κ ∈ ker(Gsum), κ_e = 2 − α, −α ≤ κ_g ≤ 1 − α (g ≠ e),
  α = 6/(p n), and 1ᵀκ = 0.

PROVED Max+-free Fraction (four-line box min, independent of full ker):
  A. Negative / positive mass.  N_e := ∑_{f⊥e, a_f<0} |a_f|,
     P_e := ∑_{f⊥e, a_f>0} a_f, so P_e − N_e = n − 2.
  B. e-row of Gsum κ = 0 (wedges 0):
        2 κ_e + ∑_{f⊥e} a_f κ_f = 0.
  C. Box min of ∑ a_f κ_f over −α ≤ κ_f ≤ 1−α independently:
        ∑ a_f κ_f ≥ −α P_e − (1−α) N_e.
     Hence
        κ_e ≤ (α P_e + (1−α) N_e)/2 = (α(n−2) + N_e)/2.
  D. Dual-eq needs κ_e = 2 − α.  Sufficient obstruction:
        (α(n−2) + N_e)/2 < 2 − α
     ⇔  N_e < 4 − α n = 4 − 6/p.
     So if N_e < 4 − 6/p for every edge e and every prime p≥5, the
     continuous dual-eq ker-box is empty and residual (i) closes.
     **No** pointwise Gsum≥−1/p and **no** |μ₄|≤2/n required.

CERTIFIED (not general; exact Fraction / float census):
  E. Paley N_e at a fixed edge (edge-transitive):
        p=3: N_e = 16/3 > thr = 2
        p=5: N_e = 384/65 > thr = 14/5
        p=7: N_e ≈ 13.077 > thr = 22/7
     Pure N_e < 4−6/p is **FALSE** for actual Paley Gsum at p=3,5,7.
     The independent-box bound is too weak: ub κ_e ≈ 3.5 > need at p=5.
  F. Refined row+mass LP (e-row + 1ᵀκ=0 + box; no full ker):
        p=3: max κ_e = 14/5 = 2.8 ≥ need 9/5   (does **not** block)
        p=5: max κ_e ≈ 1.189 < need 127/65≈1.954  (blocks)
        p=7: max κ_e ≈ 1.166 < need 347/175≈1.983 (blocks)
     Matches the auditor's ~1.19 / ~1.17.  So e-row + mass alone already
     kills dual-eq at p=5,7; full ker is even tighter (15.190/192).

OPEN (blocks residual i / predicate flip):
  G. Pure N_e < 4−6/p is **not** a viable general target (fails census).
  H. Prove Max+-free that the **row+mass** max
        max{κ_e : (Gsum κ)_e=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)}
     is < 2−α for all primes p≥5
     (or full ker-box empty, or |μ|≤2/n on |κ|=1).
  I. Optional refined L¹ target: incorporate mass via a Lagrange shift
        ∑ (a_f − λ) κ_f, optimise λ, get a weighted negative-mass bound
     that is Max+-free and beats need for p≥5.

Until H: residual_i False; gsum_disj_lb False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15194.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
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


# ---------------------------------------------------------------------------
# A–D: four-line obstruction (Max+-free Fraction algebra)
# ---------------------------------------------------------------------------


def N_e_threshold(p: int) -> Fraction:
    """4 − 6/p: pure independent-box obstruction threshold."""
    return 4 - Fraction(6, p)


def kappa_e_indep_upper_bound(p: int, N_e: Fraction) -> Fraction:
    """(α(n−2) + N_e)/2 from independent box min of ∑ a_f κ_f."""
    n = n_of(p)
    alpha = alpha_particular(p)
    return (alpha * (n - 2) + N_e) / 2


def dual_eq_need_kappa_e(p: int) -> Fraction:
    return 2 - alpha_particular(p)


def theorem_row_negative_mass_obstruction() -> dict:
    """
    Proved: if N_e < 4−6/p for every edge e, dual-eq ker-box is empty.
    Algebra only; does not claim N_e < thr for actual Gsum.
    """
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        thr = N_e_threshold(p)
        alpha = alpha_particular(p)
        n = n_of(p)
        # algebra identity: (α(n-2)+N)/2 < 2-α  ⇔  N < 4-αn = 4-6/p
        # check at boundary N = thr - eps conceptually via identity
        # α n = 6/p
        id_ok = alpha * n == Fraction(6, p)
        # at N = thr: ub = (α(n-2) + 4 - 6/p)/2 = (α(n-2) + 4 - α n)/2
        # = (4 - 2α)/2 = 2 - α = need  (equality, not strict)
        N_eq = thr
        ub_eq = kappa_e_indep_upper_bound(p, N_eq)
        need = dual_eq_need_kappa_e(p)
        # strict inequality N < thr ⇒ ub < need
        N_strict = thr - Fraction(1, 10**9)
        ub_strict = kappa_e_indep_upper_bound(p, N_strict)
        rows[str(p)] = {
            "thr": str(thr),
            "ub_at_thr": str(ub_eq),
            "need": str(need),
            "equality_at_thr": ub_eq == need,
            "strict_blocks": ub_strict < need,
            "alpha_n_eq_6_over_p": id_ok,
        }
        if not (id_ok and ub_eq == need and ub_strict < need):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Under Gsum_ee=2, wedges=0, row sum n: if for every edge e the "
            "disj negative mass N_e < 4−6/p, then independent box min of the "
            "e-row forces κ_e ≤ (α(n−2)+N_e)/2 < 2−α, so dual-eq κ_e=2−α "
            "is impossible. Continuous ker-box empty ⇒ residual (i) closes. "
            "No pointwise Gsum LB and no |μ₄| bound required."
        ),
        "proof_sketch": [
            "P_e − N_e = n−2",
            "2κ_e + ∑ a_f κ_f = 0",
            "∑ a_f κ_f ≥ −α P_e − (1−α) N_e  (box min)",
            "κ_e ≤ (α(n−2)+N_e)/2",
            "κ_e < 2−α  ⇔  N_e < 4−6/p",
        ],
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13") if k in rows},
        "n_primes": len(primes),
        "note": (
            "Sufficient condition only. Census (E) shows N_e > thr for "
            "Paley at p=3,5,7, so the pure N_e target does not hold there."
        ),
    }


def theorem_threshold_algebra(p: int) -> dict:
    """Single-p check of the N_e ↔ κ_e translation."""
    thr = N_e_threshold(p)
    need = dual_eq_need_kappa_e(p)
    ub = kappa_e_indep_upper_bound(p, thr)
    return {
        "p": p,
        "N_e_threshold": str(thr),
        "need_kappa_e": str(need),
        "ub_at_threshold": str(ub),
        "ub_equals_need_at_thr": ub == need,
        "proved": ub == need,
    }


# ---------------------------------------------------------------------------
# E–F: census N_e and row+mass LP
# ---------------------------------------------------------------------------


def _maxplus_minus(p: int):
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    if p == 7:
        cache = Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7")
        if not (cache / "maxplus.npy").is_file():
            alt = Path("/tmp/e1_p7")
            cache = alt if (alt / "maxplus.npy").is_file() else cache
        if (cache / "maxplus.npy").is_file():
            Mp = np.load(cache / "maxplus.npy").astype(float)
            Mm = np.load(cache / "maxminus.npy").astype(float)
            return C, Mp, Mm
        return C, None, None
    if p > 7:
        return C, None, None
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    return C, Mp, Mm


def census_N_e_exact(p: int, edge: tuple[int, int] = (0, 1)) -> dict:
    """
    Exact Fraction N_e for Paley at edge (default (0,1)).
    p=7 uses float then limit_denominator if Max+ cache present.
    """
    import numpy as np

    C, Mp, Mm = _maxplus_minus(p)
    if Mp is None:
        return {"p": p, "skipped": True, "reason": "no Max± for this p"}
    n = int(C.shape[0])
    edges = list(combinations(range(n), 2))
    e = tuple(sorted(edge))
    if e not in edges:
        e = (0, 1)
    disj = [f for f in edges if len(set(e) & set(f)) == 0]
    Np, Nm = len(Mp), len(Mm)

    # Prefer exact integer dots when Max± are ±1
    use_exact = p <= 5
    if use_exact:
        Mp_i = Mp.astype(np.int8)
        Mm_i = Mm.astype(np.int8)
        C_i = C.astype(np.int8)
        Ce = int(C_i[e[0], e[1]])
        fe_p = Ce * Mp_i[:, e[0]] * Mp_i[:, e[1]]
        fe_m = Ce * Mm_i[:, e[0]] * Mm_i[:, e[1]]
        N_e = Fraction(0)
        P_e = Fraction(0)
        value_mult: dict[str, int] = {}
        for f in disj:
            i, j = f
            Cf = int(C_i[i, j])
            ff_p = Cf * Mp_i[:, i] * Mp_i[:, j]
            ff_m = Cf * Mm_i[:, i] * Mm_i[:, j]
            sp = int(np.dot(fe_p.astype(np.int64), ff_p.astype(np.int64)))
            sm = int(np.dot(fe_m.astype(np.int64), ff_m.astype(np.int64)))
            a = Fraction(sp, Np) + Fraction(sm, Nm)
            value_mult[str(a)] = value_mult.get(str(a), 0) + 1
            if a < 0:
                N_e += -a
            elif a > 0:
                P_e += a
    else:
        Ce = float(C[e[0], e[1]])
        fe_p = Ce * Mp[:, e[0]] * Mp[:, e[1]]
        fe_m = Ce * Mm[:, e[0]] * Mm[:, e[1]]
        N_f = 0.0
        P_f = 0.0
        value_mult = {}
        for f in disj:
            i, j = f
            Cf = float(C[i, j])
            a = float(
                np.mean(fe_p * (Cf * Mp[:, i] * Mp[:, j]))
                + np.mean(fe_m * (Cf * Mm[:, i] * Mm[:, j]))
            )
            # bucket for report
            key = str(Fraction(a).limit_denominator(100000))
            value_mult[key] = value_mult.get(key, 0) + 1
            if a < -1e-12:
                N_f += -a
            elif a > 1e-12:
                P_f += a
        N_e = Fraction(N_f).limit_denominator(1000000)
        P_e = Fraction(P_f).limit_denominator(1000000)

    thr = N_e_threshold(p)
    alpha = alpha_particular(p)
    ub = kappa_e_indep_upper_bound(p, N_e)
    need = dual_eq_need_kappa_e(p)
    # P − N should equal n−2
    return {
        "p": p,
        "edge": list(e),
        "n_disj": len(disj),
        "N_e": str(N_e),
        "P_e": str(P_e),
        "P_minus_N": str(P_e - N_e),
        "n_minus_2": n - 2,
        "row_sum_disj_ok": (P_e - N_e) == (n - 2)
        or abs(float(P_e - N_e) - (n - 2)) < 1e-6,
        "threshold": str(thr),
        "N_e_lt_threshold": N_e < thr,
        "indep_ub_kappa_e": str(ub),
        "need_kappa_e": str(need),
        "indep_ub_blocks_dual_eq": ub < need,
        "value_multiplicities": {
            k: value_mult[k]
            for k in sorted(value_mult.keys(), key=lambda s: Fraction(s))
        },
        "exact": use_exact,
        "proved_general": False,
    }


def certify_row_mass_lp_max_kappa_e(p: int) -> dict:
    """
    max κ_e s.t. (Gsum κ)_e = 0, 1ᵀκ = 0, box on g≠e.
    Certified census only — not a general proof.
    """
    import numpy as np
    from scipy.optimize import linprog

    C, Mp, Mm = _maxplus_minus(p)
    if Mp is None:
        return {"p": p, "skipped": True, "reason": "no Max± for this p"}
    n = int(C.shape[0])
    edges = list(combinations(range(n), 2))
    E = len(edges)
    ei = edges.index((0, 1))
    Fp = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mp])
    Fm = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mm])
    Gsum = (Fp.T @ Fp) / len(Mp) + (Fm.T @ Fm) / len(Mm)
    alpha = float(alpha_particular(p))
    need = float(dual_eq_need_kappa_e(p))
    c = np.zeros(E)
    c[ei] = -1.0
    A_eq = np.vstack([Gsum[ei : ei + 1], np.ones((1, E))])
    b_eq = np.zeros(2)
    bounds = [
        (None, None) if i == ei else (-alpha, 1.0 - alpha) for i in range(E)
    ]
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    max_ke = float(-res.fun) if res.success else None
    # independent ub for comparison
    Ne_c = census_N_e_exact(p)
    indep_ub = float(Fraction(Ne_c["indep_ub_kappa_e"])) if not Ne_c.get("skipped") else None
    return {
        "p": p,
        "max_kappa_e": max_ke,
        "max_kappa_e_frac": (
            str(Fraction(max_ke).limit_denominator(100000))
            if max_ke is not None
            else None
        ),
        "need": need,
        "need_frac": str(dual_eq_need_kappa_e(p)),
        "blocks_dual_eq": max_ke is not None and max_ke < need - 1e-9,
        "indep_ub_kappa_e": indep_ub,
        "success": bool(res.success),
        "proved_general": False,
        "note": (
            "e-row + mass + box only (not full ker). "
            "Blocks dual-eq at p=5,7; not at p=3."
        ),
    }


# ---------------------------------------------------------------------------
# Predicates
# ---------------------------------------------------------------------------


def row_negative_mass_bound_proved_general() -> bool:
    """
    True only when Max+-free proof that N_e < 4−6/p for all p≥5 is shipped.
    Census shows the inequality fails at p=5,7 — stay False.
    """
    return False


def row_mass_lp_blocks_dual_eq_general() -> bool:
    """True only with Max+-free proof of row+mass max κ_e < 2−α for p≥5."""
    return False


def residual_i_closed_via_row_mass() -> bool:
    """Residual (i) via this prop — False until a general bound is proved."""
    return bool(
        row_negative_mass_bound_proved_general()
        or row_mass_lp_blocks_dual_eq_general()
    )


def hinge_status_194() -> dict:
    return {
        "row_negative_mass_sufficient_lemma": True,  # algebra proved
        "pure_N_e_lt_thr_for_paley": False,  # fails census p=3,5,7
        "row_mass_lp_blocks_p5_p7": True,  # certified
        "row_mass_lp_blocks_general": False,
        "residual_i_closed": residual_i_closed_via_row_mass(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove Max+-free that row+mass max κ_e < 2−α for all p≥5 "
            "(or full ker-box empty, or |μ|≤2/n). Pure N_e < 4−6/p fails "
            "Paley census at p=3,5,7 and is not a viable target as stated."
        ),
    }


def main() -> dict:
    lemma = theorem_row_negative_mass_obstruction()
    thr_checks = {str(p): theorem_threshold_algebra(p) for p in (5, 7, 11, 13)}
    census = {}
    for p in (3, 5, 7):
        census[str(p)] = census_N_e_exact(p)
    row_lp = {}
    for p in (3, 5, 7):
        row_lp[str(p)] = certify_row_mass_lp_max_kappa_e(p)

    # sanity: pure N_e does not block; row+mass blocks p≥5 when available
    pure_fails = all(
        not census[str(p)].get("N_e_lt_threshold", True)
        for p in (3, 5)
        if not census[str(p)].get("skipped")
    )
    row_blocks_ge5 = all(
        row_lp[str(p)].get("blocks_dual_eq") is True
        for p in (5, 7)
        if not row_lp[str(p)].get("skipped")
    )

    out = {
        "title": (
            "Prop 15.194 row negative-mass obstruction: sufficient lemma proved; "
            "pure N_e target fails Paley census; row+mass LP blocks p=5,7; "
            "residual (i) still OPEN"
        ),
        "proved": {
            "row_negative_mass_sufficient_lemma": lemma["proved"],
            "threshold_algebra": all(v["proved"] for v in thr_checks.values()),
            "pure_N_e_lt_4_minus_6_over_p_general": False,
            "row_mass_lp_blocks_general": False,
            "residual_i_closed": residual_i_closed_via_row_mass(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "sufficient_lemma": lemma,
            "threshold_checks": thr_checks,
        },
        "certified": {
            "N_e_census": census,
            "row_mass_lp": row_lp,
            "pure_N_e_fails_at_p3_p5": pure_fails,
            "row_mass_blocks_at_p5_p7": row_blocks_ge5,
        },
        "hinge": hinge_status_194(),
        "F3": (
            "Do not flip residual (i)/E1/L: pure N_e < 4−6/p fails census; "
            "row+mass LP block at p=5,7 is evidence, not a general proof."
        ),
        "audit_note": (
            "Sufficient lemma is the four-line box min. Pure L¹ target is "
            "too weak for Paley. Prefer proving row+mass max κ_e < 2−α "
            "Max+-free (census ~1.19/1.17 at p=5,7 vs need ~1.95/1.98)."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15194.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.194 row negative-mass obstruction")
    print(f"  sufficient lemma proved: {lemma['proved']}")
    for p in (3, 5, 7):
        c = census[str(p)]
        if c.get("skipped"):
            print(f"  N_e p={p}: skipped")
            continue
        print(
            f"  N_e p={p}: {c['N_e']} thr={c['threshold']} "
            f"lt_thr={c['N_e_lt_threshold']} indep_blocks={c['indep_ub_blocks_dual_eq']}"
        )
    for p in (3, 5, 7):
        r = row_lp[str(p)]
        if r.get("skipped"):
            print(f"  row+mass LP p={p}: skipped")
            continue
        print(
            f"  row+mass LP p={p}: max_ke={r.get('max_kappa_e_frac')} "
            f"need={r.get('need_frac')} blocks={r.get('blocks_dual_eq')}"
        )
    print(f"  residual_i via row mass: {residual_i_closed_via_row_mass()}")
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
