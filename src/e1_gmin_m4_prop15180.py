#!/usr/bin/env python3
"""
Prop 15.180 — Residual (i) attack log: dual-eq support constraints; hinge OPEN.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

CONTEXT
  Residual (ii) CLOSED (15.179 freeze). Residual (i) Type I k=3p−2 remains
  the unique denseness-path ND gap for E(1).

PROVED Max+-free (this prop + 15.176–178):
  A. Dual equality (e∉G) forces (Gsum x)_e = need = 6/p−4 (15.176).
  B. Wedge partners contribute 0 ⇒ need = sum of n_d disj Gsums only.
  C. PSD Gsum≥−2 ⇒ dual eq impossible for n_d∈{0,1}, all primes p≥5 (15.178).
  D. n_d=2 wedge partners (Gsum between them 0): box≥−2√2, kills dual eq
     for all primes p≥7 (15.178).
  E. Score parity under dual affine: for each f_e=±1, exist S_W with
     |S_W|≤n_w, |S_W − (3−2f_e)|≤n_d, S_W ≡ n_w (mod 2). Filters some
     (n_d,n_w) but leaves n_d≥2 open for all p≥5.
  F. Under dual equality, sum of Gsums over disj *pairs within G* equals
        Q_pairs = 2(13p−12)/p − 2k = 30 − 6p − 24/p < 0 (p≥5).
     Hence G must contain at least one disj edge-pair (else Q_pairs=0).
  G. Farkas: Gsum≥−1/p on all disj ⇒ need < k(−1/p) for all primes p≥5
     (15.176). Equivalent |μ₄|≤1/(2p) on |κ|=1 (15.177).

CERTIFIED (not general):
  H. p=5: min Gsum=−6/65 > −1/p > μ_*; dual eq obstructed by census.
  I. p=7: max|m₄|_{κ1}<1/(2p) (prior census).

OPEN (blocks residual i / E1 / L):
  J. Prove disj Gsum≥−1/p (or |μ₄|≤1/(2p) on |κ|=1, or any μ>μ_*, or
     dual-eq impossible for all remaining supports n_d≥2) Max+-free for
     all primes p≥5.
  Dead / too weak for general close: weak CS, foreign-star CS, K4 U-PSD,
  absolute T bootstrap, scheme-min −12/(pn), non-uniform box with only
  PSD max≤2, freeze (S_H not constant for Type I).

Until J: gsum_disj_lb_proved_general False; type_I closed False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15180.json
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)
from e1_gmin_m4_prop15176 import (  # noqa: E402
    mu_minus_1_over_p,
    mu_star_threshold,
    residual_i_farkas_if_mu,
    theorem_minus_1_over_p_suffices,
)
from e1_gmin_m4_prop15178 import (  # noqa: E402
    dual_eq_box_lower_psd,
    nd2_wedge_partners_lb,
    theorem_nd2_wedge_kill_p_ge_7,
    theorem_nd_le_1_impossible,
)


def Q_pairs_dual_eq(p: int) -> Fraction:
    """Sum of Gsum over unordered pairs within G under dual equality."""
    k = k_type_I_3p_minus_2(p)
    ES2 = Fraction(13 * p - 12, p)  # both Max± under dual
    # E+[S²]+E-[S²] = 2k + sum_{pairs} Gsum
    return 2 * ES2 - 2 * k


def theorem_Q_pairs_negative(primes: list[int] | None = None) -> dict:
    """Q_pairs = 30−6p−24/p < 0 for all primes p≥5."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        qp = Q_pairs_dual_eq(p)
        # closed form check
        closed = Fraction(30) - 6 * p - Fraction(24, p)
        match = qp == closed
        neg = qp < 0
        rows[str(p)] = {
            "Q_pairs": str(qp),
            "closed_form": str(closed),
            "match": match,
            "negative": neg,
        }
        if not (match and neg):
            ok = False
    return {
        "proved": ok,
        "formula": "30 - 6p - 24/p",
        "theorem": (
            "Under dual equality ES2+=(ES2−)= (13p−12)/p, "
            "sum of Gsum over pairs in G equals 2 ES2−2k=30−6p−24/p<0 for p≥5. "
            "Hence G contains at least one vertex-disjoint edge pair."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_primes": len(primes),
    }


def score_parity_viable_nd(p: int) -> dict:
    """n_d values compatible with dual affine scores S∈{1,5} and |S_D|≤n_d."""
    k = k_type_I_3p_minus_2(p)
    viable = []
    for nd in range(0, k + 1):
        nw = k - nd
        ok1 = okm = False
        for sw in range(-nw, nw + 1, 2):
            if abs(1 - sw) <= nd:
                ok1 = True
            if abs(5 - sw) <= nd:
                okm = True
        if ok1 and okm:
            viable.append(nd)
    need = dual_equality_correlation_need(p)
    open_nd = [nd for nd in viable if not (need < dual_eq_box_lower_psd(nd))]
    return {
        "p": p,
        "viable_nd": viable,
        "open_after_psd": open_nd,
        "n_open": len(open_nd),
        "min_open_nd": min(open_nd) if open_nd else None,
    }


def theorem_open_core_nd_ge_2(primes: list[int] | None = None) -> dict:
    """For all primes p≥5, dual-eq support after PSD still has some n_d≥2."""
    if primes is None:
        primes = [p for p in range(5, 40) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        r = score_parity_viable_nd(p)
        has_open = r["n_open"] > 0 and r["min_open_nd"] is not None and r["min_open_nd"] >= 2
        rows[str(p)] = {
            "n_open": r["n_open"],
            "min_open_nd": r["min_open_nd"],
            "open_core_nonempty": has_open,
        }
        if not has_open:
            ok = False
    return {
        "proved_open_core_nonempty": ok,
        "note": (
            "PSD kills n_d≤1; score parity leaves n_d≥2 viable. "
            "This records the open dual-eq support, not a close."
        ),
        "by_p": rows,
    }


def hinge_status_180() -> dict:
    return {
        "residual_ii_closed": True,  # 15.179
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "nd_le_1_impossible": True,
        "nd2_wedge_p_ge_7_impossible": True,
        "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(),
        "bound_proved_general": False,
        "open": (
            "Prove disj Gsum≥−1/p or |μ₄|≤1/(2p) on |κ|=1 (or dual-eq "
            "impossible on all n_d≥2 supports) for all primes p≥5 Max+-free."
        ),
        "mu_star_p5": str(mu_star_threshold(5)),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    qp = theorem_Q_pairs_negative(primes)
    nd1 = theorem_nd_le_1_impossible(primes)
    nd2 = theorem_nd2_wedge_kill_p_ge_7(primes)
    core = theorem_open_core_nd_ge_2([p for p in range(5, 40) if is_prime(p)])
    sample_farkas = {
        str(p): residual_i_farkas_if_mu(p, mu_minus_1_over_p(p))
        for p in (5, 7, 11, 13)
    }
    out = {
        "title": (
            "Prop 15.180 residual (i) dual-eq support constraints; "
            "Gsum LB / |μ₄| bound still OPEN"
        ),
        "proved": {
            "Q_pairs_formula_negative": qp["proved"],
            "nd_le_1_impossible": nd1["proved"],
            "nd2_wedge_p_ge_7": nd2["proved_for_primes_ge_7"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(primes),
            "open_core_nd_ge_2_nonempty": core["proved_open_core_nonempty"],
            "mu4_or_gsum_bound_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "Q_pairs": qp,
            "nd_le_1": nd1,
            "nd2_wedge": nd2,
            "open_core": core,
            "farkas_if_minus_1_over_p": sample_farkas,
            "nd2_wedge_lb": nd2_wedge_partners_lb(),
        },
        "hinge": hinge_status_180(),
        "F3": "Do not flip residual (i)/E1/L without general Gsum or μ₄ LB",
        "dead_approaches": [
            "weak CS Gsum≥4/p²−2",
            "foreign-star CS / K4 U-PSD",
            "absolute T row-sum bootstrap (4p<d1)",
            "scheme-min −12/(pn) (15.158; false at p=3)",
            "Type I freeze (S_H∈{2,4} not constant)",
            "non-uniform box with only PSD |Gsum|≤2",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15180.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.180 residual (i) support constraints")
    print(f"  Q_pairs<0 proved={qp['proved']}")
    print(f"  n_d≤1 kill={nd1['proved']}; n_d=2 wedge p≥7={nd2['proved_for_primes_ge_7']}")
    print(f"  open core n_d≥2 nonempty={core['proved_open_core_nonempty']}")
    print(f"  gsum/E1/L still open (bound_proved_general=False)")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
