#!/usr/bin/env python3
"""
Prop 15.196 — Spectral row-energy bound for mass-corrected dual-eq (residual i).

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.192–195)
  a_f = Gsum_ef for f ⊥ e; wedges 0; ∑_{f⊥e} a_f = n−2.
  Q_e := ∑_{f⊥e} a_f²  (disj row energy).
  Dual-eq needs mass-min ∑ a κ ≤ −2M, M=2−α (15.195).
  Gsum ≽ 0 ⇒ every entry a_f ≥ −2 (2×2 principal: diag 2).

PROVED Max+-free
  A. Entry floor. Gsum_ee=2, Gsum ≽ 0 ⇒ a_f ≥ −2 for all f.  ∎

  B. Spectral energy identity. Write e_* = E^{-1}1 + e_⊥ with e_⊥⊥1,
     E=#edges.  Then e_*^T Gsum e_* = 2 = n/E + e_⊥^T Gsum e_⊥, so
        e_⊥^T Gsum e_⊥ = 2 − n/E = 2(n−2)/(n−1).
     And ||Gsum e_*||² = Q_e + 4 = n²/E + e_⊥^T Gsum² e_⊥.
     On 1^⊥, Gsum ≼ λ₂ I (λ₂ = second-largest eigenvalue of Gsum), hence
     Gsum² ≼ λ₂ Gsum on 1^⊥, and
        Q_e + 4 ≤ n²/E + λ₂ · 2(n−2)/(n−1)
                 = 2n/(n−1) + 2 λ₂ (n−2)/(n−1)
                 = 2( n + λ₂(n−2) ) / (n−1).
     Thus
        Q_e ≤ Q_ub(λ₂) := 2(n + λ₂(n−2))/(n−1) − 4.  ∎

  C. Mass-min under a≥−2 and ∑a=n−2, ∑a²≤Q.
     For two-level a with k entries at −2 and the rest equal (the
     sum-constrained extreme for negative mass), the mass-min of ∑aκ
     is an explicit Fraction function of (p,k,Q).  If this worst-case
     mass-min is > −2M for every admissible k, dual-eq is impossible.
     (Proved by the 15.195 criterion + finite check of admissible k.)

  D. Implication for residual (i).  If for every prime p≥5 and every edge e
        Q_e ≤ Q_*(p)
     where Q_*(p) is large enough that (C) blocks dual-eq (e.g. Q_*(p)=10
     works for all primes p∈[5,47] checked; critical Q_max≈10.86 at p=5),
     then dual-eq ker-box is empty and residual (i) closes.  ∎

CERTIFIED (not general)
  E. Paley Q_e (edge-transitive, one edge):
        p=3: Q=112/9≈12.444  (dual-eq feasible; critical boundary)
        p=5: Q=34512/4225≈8.169
        p=7: Q=55009/8220≈6.692
     Both p=5,7 have Q < 10 and mass-min under true Q blocks (15.195).
  F. Paley λ₂(Gsum):
        p=3: 8 (mult n)
        p=5: 88/13≈6.769 (mult n)
        p=7: 2160/409≈5.281 (mult n)
     Spectral Q_ub(λ₂) equals actual Q at p=3 (tight); strict at p=5,7.
     At p=5, Q_ub(λ₂)=144/13≈11.077 > critical Q_max≈10.86 — spectral
     ub alone does **not** quite block at p=5 (needs λ₂≲6.208 or actual Q).

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free for all primes p≥5 one of:
     (i)  Q_e ≤ 10 (or ≤ Q_*(p) with Q_* from (C)), or
     (ii) λ₂(Gsum) ≤ (6n−7)/(n−2) = 6 + 5/(n−2)  (⇒ Q_ub≤10), or
     (iii) |μ|≤2/n on |κ|=1, or full ker-box empty.
     Census suggests Q decreasing and λ₂ decreasing — promising but open.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15196.json
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
from e1_gmin_m4_prop15194 import dual_eq_need_kappa_e  # noqa: E402


def theorem_entry_floor_minus_2() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Gsum_ee=2 and Gsum≽0 ⇒ every off-diagonal entry a≥−2 "
            "(principal 2×2 minor)."
        ),
    }


def Q_ub_from_lambda2(p: int, lam2: Fraction) -> Fraction:
    """Q ≤ 2(n + λ₂(n−2))/(n−1) − 4."""
    n = n_of(p)
    return 2 * (n + lam2 * (n - 2)) / (n - 1) - 4


def theorem_spectral_Q_ub() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Q_e + 4 = ||Gsum e_*||² ≤ 2(n + λ₂(n−2))/(n−1), hence "
            "Q_e ≤ 2(n + λ₂(n−2))/(n−1) − 4, where λ₂ is the second-largest "
            "eigenvalue of Gsum. Uses only Gsum1=n1, Gsum_ee=2, and "
            "Gsum² ≼ λ₂ Gsum on 1^⊥."
        ),
        "formula": "Q_ub(λ₂)=2(n+λ₂(n-2))/(n-1)-4",
    }


def lambda2_for_Q_ub_le_10(p: int) -> Fraction:
    """λ₂ ≤ (6n−7)/(n−2) ⇒ Q_ub ≤ 10."""
    n = n_of(p)
    return Fraction(6 * n - 7, n - 2)


def theorem_lambda2_threshold_for_Q10() -> dict:
    primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        thr = lambda2_for_Q_ub_le_10(p)
        # check Q_ub(thr) == 10
        qub = Q_ub_from_lambda2(p, thr)
        rows[str(p)] = {
            "lambda2_thr": str(thr),
            "Q_ub_at_thr": str(qub),
            "equals_10": qub == 10,
        }
        if qub != 10:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "λ₂ ≤ (6n−7)/(n−2) = 6+5/(n−2) ⇒ Q_ub(λ₂) ≤ 10 for all p≥5. "
            "If also Q_e ≤ Q_ub and Q≤10 forces mass-min > −2M, residual (i) closes."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13") if k in rows},
    }


def mass_min_twolevel(
    k: int, r: Fraction, n_disj: int, n_wedge: int, M: Fraction, alpha: Fraction
) -> Fraction:
    """Mass-min ∑aκ for k copies of −2, n_disj−k of r, n_wedge of 0."""
    lo, hi = -alpha, 1 - alpha
    width = hi - lo
    B = -M
    n_tot = n_disj + n_wedge
    rem = B - n_tot * lo
    items = [(-2, k), (Fraction(0), n_wedge), (r, n_disj - k)]
    items.sort(key=lambda t: t[0])
    sum_a = lo * (k * (-2) + (n_disj - k) * r)
    for cost, cnt in items:
        if rem <= 0 or cnt == 0:
            continue
        take_total = min(rem, width * cnt)
        sum_a += cost * take_total
        rem -= take_total
    return sum_a


def worst_mass_min_under_Q(p: int, Qmax: Fraction) -> dict:
    """
    Worst two-level mass-min under a≥−2, ∑a=n−2, ∑a²≤Qmax.
    Returns whether this worst-case still exceeds −2M (blocks dual-eq).
    """
    n = n_of(p)
    alpha = alpha_particular(p)
    M = dual_eq_need_kappa_e(p)
    target = -2 * M
    n_wedge = 2 * (n - 2)
    n_disj = (n - 2) * (n - 3) // 2
    S = n - 2
    kmax = min(n_disj - 1, int(Qmax / 4) + 5)
    worst = None
    worst_k = None
    for k in range(0, kmax + 1):
        r = Fraction(S + 2 * k, n_disj - k)
        if r < -2 or r > 2:
            continue
        q = 4 * k + (n_disj - k) * r * r
        if q > Qmax:
            continue
        s = mass_min_twolevel(k, r, n_disj, n_wedge, M, alpha)
        if worst is None or s < worst:
            worst = s
            worst_k = k
    return {
        "p": p,
        "Qmax": str(Qmax),
        "worst_mass_min": str(worst) if worst is not None else None,
        "target_minus_2M": str(target),
        "blocks_dual_eq": worst is not None and worst > target,
        "worst_k": worst_k,
    }


def theorem_Q_le_10_blocks_p_ge_5() -> dict:
    """
    Checkable Fraction: for all primes p∈[5,47], Qmax=10 forces
    two-level worst mass-min > −2M.  (Finite check; pattern stable.)
    """
    primes = [p for p in range(5, 48) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        r = worst_mass_min_under_Q(p, Fraction(10))
        rows[str(p)] = r
        if not r["blocks_dual_eq"]:
            ok = False
    return {
        "proved_for_primes_5_to_47": ok,
        "proved_general_all_primes": False,  # finite check only
        "theorem": (
            "For every prime p with 5≤p≤47: under a≥−2, ∑a=n−2, ∑a²≤10, "
            "the two-level mass-min is > −2(2−α).  If this extends to all "
            "p≥5 and Q_e≤10 Max+-free, residual (i) closes."
        ),
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13", "47") if k in rows},
        "note": "Finite check is evidence for the implication shape; general p needs proof.",
    }


def certified_paley_Q_and_lambda2() -> dict:
    """Census values — not general."""
    return {
        "3": {
            "Q": str(Fraction(112, 9)),
            "lambda2": str(Fraction(8)),
            "Q_ub_lambda2": str(Q_ub_from_lambda2(3, Fraction(8))),
            "blocks_under_actual_Q": False,
        },
        "5": {
            "Q": str(Fraction(34512, 4225)),
            "lambda2": str(Fraction(88, 13)),
            "Q_ub_lambda2": str(Q_ub_from_lambda2(5, Fraction(88, 13))),
            "blocks_under_actual_Q": True,
            "blocks_under_Q10": True,
            "blocks_under_Q_ub": worst_mass_min_under_Q(
                5, Q_ub_from_lambda2(5, Fraction(88, 13))
            )["blocks_dual_eq"],
        },
        "7": {
            "Q": str(Fraction(55009, 8220)),
            "lambda2": str(Fraction(2160, 409)),
            "Q_ub_lambda2": str(Q_ub_from_lambda2(7, Fraction(2160, 409))),
            "blocks_under_actual_Q": True,
            "blocks_under_Q10": True,
            "blocks_under_Q_ub": worst_mass_min_under_Q(
                7, Q_ub_from_lambda2(7, Fraction(2160, 409))
            )["blocks_dual_eq"],
        },
        "proved_general": False,
    }


def Q_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_spectral_Q() -> bool:
    return Q_bound_proved_general()


def hinge_status_196() -> dict:
    return {
        "entry_floor_minus_2": True,
        "spectral_Q_ub_formula": True,
        "Q_le_10_blocks_checked_p_le_47": True,
        "Q_bound_general": False,
        "lambda2_bound_general": False,
        "residual_i_closed": residual_i_closed_via_spectral_Q(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove Max+-free Q_e≤10 (or λ₂≤6+5/(n−2), or |μ|≤2/n, or "
            "ker-box empty) for all primes p≥5."
        ),
    }


def main() -> dict:
    A = theorem_entry_floor_minus_2()
    B = theorem_spectral_Q_ub()
    C = theorem_lambda2_threshold_for_Q10()
    D = theorem_Q_le_10_blocks_p_ge_5()
    cert = certified_paley_Q_and_lambda2()
    out = {
        "title": (
            "Prop 15.196 spectral row-energy bound for mass-corrected dual-eq; "
            "Q≤10 would close residual (i); general Q/λ₂ bound OPEN"
        ),
        "proved": {
            "entry_floor_minus_2": A["proved"],
            "spectral_Q_ub_formula": B["proved"],
            "lambda2_threshold_for_Q10": C["proved"],
            "Q_le_10_blocks_primes_5_to_47": D["proved_for_primes_5_to_47"],
            "Q_bound_general": False,
            "residual_i_closed": residual_i_closed_via_spectral_Q(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "entry_floor": A,
            "spectral_Q_ub": B,
            "lambda2_thr": C,
            "Q10_blocks": D,
        },
        "certified": cert,
        "hinge": hinge_status_196(),
        "F3": (
            "Q≤10 (or λ₂≤6+5/(n−2)) Max+-free would close residual (i) via "
            "15.195 mass-min. Census Q≈8.17/6.69 at p=5,7. No soft-close."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15196.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.196 spectral row-energy / mass-corrected dual-eq")
    print(f"  entry floor a≥−2: {A['proved']}")
    print(f"  spectral Q_ub formula: {B['proved']}")
    print(f"  λ₂ thr for Q_ub≤10: {C['proved']}")
    print(f"  Q≤10 blocks p=5..47: {D['proved_for_primes_5_to_47']}")
    for pk in ("5", "7"):
        c = cert[pk]
        print(
            f"  Paley p={pk}: Q={c['Q']} λ2={c['lambda2']} "
            f"Q_ub={c['Q_ub_lambda2']} blocks_Q10={c['blocks_under_Q10']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_spectral_Q()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
