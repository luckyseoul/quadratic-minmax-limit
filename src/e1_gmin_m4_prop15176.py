#!/usr/bin/env python3
"""
Prop 15.176 — Correct dual-equality Farkas threshold (e ∉ G).

Corrects a wrong simplification in Prop 15.172.

SETUP (15.169–170, Max+-free Fraction):
  Freeness-fail Type I: edge set G with |G|=k=3p−2, distinguished freeness
  edge e ∉ G, affine S_G + 2 f_e = 3 on Max+.
  Dual equality (bad case): S_G = 3−2f_e on Max+, S_G = −3−2f_e on Max−.
  Then
      (Gsum x)_e = ∑_{e'∈G} Gsum_{e e'} = 6/p − 4,   x = 1_G
  with **no** diagonal term Gsum_ee (e ∉ G).

WRONG (15.172): claimed obstruction for all off-diag Gsum ≥ μ with μ > −2/p
  by writing need_off = (6/p−4)−2 = 6/p−6 = (k−1)(−2/p). That subtracts a
  diagonal as if e∈G and uses k−1 partners. **Invalid for e∉G.**

CORRECT threshold:
  ∑_{e'∈G} Gsum_{e e'} ≥ k · μ   if every Gsum_{e e'} ≥ μ
  (wedges contribute 0 ≥ μ when μ≤0, so only disj partners in G can bind).
  Obstruction when need < k·μ, i.e.
      μ > μ_*(p) := (6/p − 4)/k = (6 − 4p)/(p(3p − 2)).
  Asymptotically μ_* ∼ −4/(3p).

SUFFICIENT proved-target LBs (any one closes residual-(i) Farkas algebra):
  (i)  Gsum_ab ≥ −1/p  for all disj ab  (since −1/p > μ_* for all primes p≥5)
  (ii) Gsum_ab ≥ −12/(p n)  (candidate; > μ_* for p≥5; still unproved general)
  (iii) any fixed μ(p) > μ_*(p)

Equivalent form of (i) via 4-set moments (certified algebra + census):
  On a 4-set S, μ₄=(E₊+E₋)[y_a y_b y_c y_d]/2 and κ=∑_{3 matchings} C_e C_{e'}.
  For |κ|=1 the three matching Gsums are {2μ₄,2μ₄,−2μ₄} (up to global sign),
  so min Gsum on S equals −2|μ₄|.
  For |κ|=3 all three Gsums share sign and are positive when |μ₄| is large (census).
  Hence global min disj Gsum = −2 max{|μ₄| : |κ|=1}.
  Therefore Gsum ≥ −1/p  ⇔  |μ₄| ≤ 1/(2p) on every |κ|=1 four-set.
  Census: p=3 max|μ₄|=1/3 > 1/6 (bound false, out of residual scope);
           p=5 max|μ₄|=3/65 < 1/10 (bound holds).

NOT sufficient:
  Gsum ≥ −2/p  (i.e. H ≥ −1/p): −2/p < μ_* for p≥5, so need ≰ k·(−2/p).

CERTIFIED: census min Gsum at p=5 is −6/65 > −1/p > μ_*; candidate equals min.
OPEN: prove any sufficient disj Gsum LB Max+-free for all primes p≥5
  (preferred: |μ₄|≤1/(2p) on |κ|=1, or Gsum≥−1/p).
  gsum_disj_lb_proved_general stays False until then.
  Residual (ii) uses the same e∉G box-sum pattern with its own need(k).
  Do NOT soft-close E1/L / package while this is open.

Does **not** soft-close E1/L.

Writes evidence/e1_gmin_m4_prop15176.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    gsum_pairwise_lower_bound,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)


def mu_star_threshold(p: int) -> Fraction:
    """
    Correct Farkas threshold: dual equality obstructed if every Gsum_{e,e'}
    for e'∈G is ≥ μ with μ > μ_*(p) = (6/p−4)/k.
    """
    k = k_type_I_3p_minus_2(p)
    need = dual_equality_correlation_need(p)
    return need / k


def mu_minus_1_over_p(p: int) -> Fraction:
    return Fraction(-1, p)


def mu_minus_2_over_p(p: int) -> Fraction:
    """Wrong 15.172 threshold (not sufficient for e∉G)."""
    return Fraction(-2, p)


def mu_beats_star(p: int, mu: Fraction) -> bool:
    """Strict: μ > μ_* ⇒ k·μ > need ⇒ box sum > need."""
    return mu > mu_star_threshold(p)


def residual_i_farkas_if_mu(p: int, mu: Fraction) -> dict:
    """e∉G Farkas: need < k·μ."""
    k = k_type_I_3p_minus_2(p)
    need = dual_equality_correlation_need(p)
    box = k * mu
    return {
        "p": p,
        "k": k,
        "need": str(need),
        "mu": str(mu),
        "box_k_mu": str(box),
        "need_lt_box": need < box,
        "mu_star": str(mu_star_threshold(p)),
        "mu_gt_star": mu > mu_star_threshold(p),
    }


def compare_thresholds(p: int) -> dict:
    star = mu_star_threshold(p)
    m1 = mu_minus_1_over_p(p)
    m2 = mu_minus_2_over_p(p)
    cand = gsum_pairwise_lower_bound(p)
    return {
        "p": p,
        "mu_star": str(star),
        "mu_minus_1_over_p": str(m1),
        "mu_minus_2_over_p": str(m2),
        "mu_candidate_minus_12_pn": str(cand),
        "minus_1_over_p_beats_star": mu_beats_star(p, m1),
        "minus_2_over_p_beats_star": mu_beats_star(p, m2),
        "candidate_beats_star": mu_beats_star(p, cand),
        "farkas_if_minus_1_over_p": residual_i_farkas_if_mu(p, m1)["need_lt_box"],
        "farkas_if_minus_2_over_p": residual_i_farkas_if_mu(p, m2)["need_lt_box"],
        "farkas_if_candidate": residual_i_farkas_if_mu(p, cand)["need_lt_box"],
        "wrong_172_claims_minus_2_over_p_suffices": False,
    }


def theorem_minus_1_over_p_suffices(primes: list[int] | None = None) -> bool:
    """Proved Fraction: −1/p > μ_*(p) for all primes p≥5."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    for p in primes:
        # -1/p > (6-4p)/(p(3p-2))  ⇔  -1 > (6-4p)/(3p-2)
        # ⇔  -(3p-2) > 6-4p  ⇔  -3p+2 > 6-4p  ⇔  p > 4
        if not (p > 4 and mu_beats_star(p, mu_minus_1_over_p(p))):
            return False
        if not residual_i_farkas_if_mu(p, mu_minus_1_over_p(p))["need_lt_box"]:
            return False
    return True


def theorem_minus_2_over_p_fails(primes: list[int] | None = None) -> bool:
    """Proved: −2/p ≤ μ_* (never strict beat) for primes p≥5."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    return all(not mu_beats_star(p, mu_minus_2_over_p(p)) for p in primes)


def hinge_status_176() -> dict:
    return {
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "correct_mu_star": "(6/p-4)/k = (6-4p)/(p(3p-2))",
        "sufficient_targets": [
            "disj Gsum ≥ -1/p (Fraction-sufficient for residual i Farkas)",
            "disj Gsum ≥ -12/(p n) (candidate; algebra OK; unproved general)",
            "any μ > μ_*(p)",
        ],
        "insufficient": "Gsum ≥ -2/p (H≥-1/p) — does not beat μ_*",
        "bound_proved_general": False,
        "open": (
            "Prove disj Gsum ≥ -1/p (or any μ>μ_*) Max+-free for all primes "
            "p≥5, then set gsum_disj_lb_proved_general True."
        ),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): compare_thresholds(p) for p in (5, 7, 11, 13, 17, 19, 23)}
    t1 = theorem_minus_1_over_p_suffices(primes)
    t2 = theorem_minus_2_over_p_fails(primes)
    out = {
        "title": (
            "Prop 15.176 correct e∉G Farkas threshold μ_*; "
            "−1/p suffices, −2/p does not; LB still OPEN"
        ),
        "proved": {
            "e_not_in_G_setup": True,
            "mu_star_formula": True,
            "minus_1_over_p_suffices_for_residual_i": t1,
            "minus_2_over_p_does_not_suffice": t2,
            "candidate_beats_star_sample": all(
                compare_thresholds(p)["candidate_beats_star"]
                for p in (5, 7, 11, 13, 17)
            ),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "by_p": rows,
            "n_primes_checked_sufficiency": len(primes),
        },
        "hinge": hinge_status_176(),
        "F3": "Do not flip residual/E1/L without proved disj Gsum LB beating μ_*",
        "retraction": (
            "Prop 15.172 claim 'μ>-2/p obstructs dual equality' is invalid "
            "for e∉G; superseded by μ_*(p) here."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15176.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.176 correct Farkas threshold (e∉G)")
    print(f"  μ_* formula OK; −1/p suffices={t1}; −2/p fails={t2}")
    print(f"  sample p=5: {rows['5']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
