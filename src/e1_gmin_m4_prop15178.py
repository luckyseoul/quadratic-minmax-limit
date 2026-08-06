#!/usr/bin/env python3
"""
Prop 15.178 — Star identity for μ₄κ; dual-equality n_d case kills; hinge OPEN.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

SETUP
  Dual equality (e∉G, |G|=k=3p−2) forces (Gsum x)_e = need = 6/p−4 (15.176).
  Write n_d = #{e'∈G : e' disjoint from e}, n_w = k−n_d (wedge partners).
  Wedge Gsum ≡ 0, so (Gsum x)_e = sum of n_d disj Gsums only.

PROVED Max+-free:
  A. Star identity. For every edge e and balanced μ₄,
        ∑_{S∋e} μ₄(S) κ(S) = 3(n−2)/2 = 3(p²−1)/2.
     Proof. ∑_S μ₄ κ = n(n−1)(n−2)/8 (average of Max± closed forms
     ∑ m₄κ = n(n−1)(n−2)/8 from Prop 15.92/99, same for both eigenspaces).
     Each 4-set is counted in binom(4,2)=6 edges, each edge has the same
     star sum by edge-transitivity of the count (or double-count):
        binom(n,2) · star = 6 · ∑_S μ₄κ ⇒ star = 3(n−2)/2. ∎
  B. PSD box (any conference Gsum ≽ 0, diag 2): every Gsum_ab ≥ −2.
     Hence (Gsum x)_e ≥ −2 n_d.
  C. Dual-equality n_d ≤ 1 impossible for all primes p≥5:
        need = 6/p−4 ≤ 6/5−4 = −14/5 < −2 ≤ −2 n_d  for n_d∈{0,1}
     (n_d=0: sum=0 > need). Strict need < −2 n_d ⇒ contradiction.
  D. n_d=2 wedge partners (the two disj partners share a vertex, so
     Gsum between them is 0). Principal 3×3 on {e,f,g}:
        [[2,g1,g2],[g1,2,0],[g2,0,2]] ≽ 0  ⇒  g1²+g2² ≤ 4
     ⇒ g1+g2 ≥ −2√2. Obstructs dual equality when
        −2√2 > 6/p−4  ⇔  p > 3/(2−√2) ≈ 5.121,
     i.e. all primes p≥7. (At p=5: −2√2 < −14/5, no kill.)
  E. Farkas algebra (15.176–177): Gsum≥−1/p on all disj (⇔ |μ₄|≤1/(2p)
     on |κ|=1) still obstructs for every n_d≤k, all primes p≥5.

CERTIFIED (not general):
  F. Star identity matches census p=3,5.
  G. p=5,7: max|μ₄|_{|κ|=1} ≤ 1/(2p) (p=5: 3/65; p=7: 109/2863 census).

OPEN:
  H. n_d≥2 matching-heavy configurations for all p≥5 (and n_d=2 at p=5)
     still require a pointwise disj Gsum LB beating μ_*=(6/p−4)/k, or an
     alternate residual. Preferred: |μ₄|≤1/(2p) on |κ|=1 (15.177).
  Until H: gsum_disj_lb_proved_general False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15178.json
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
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402


def sum_mu4_kappa_global(p: int) -> Fraction:
    """∑_S μ₄ κ = n(n−1)(n−2)/8 (balanced; Max± closed forms agree)."""
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2), 8)


def star_mu4_kappa(p: int) -> Fraction:
    """∑_{S∋e} μ₄ κ = 3(n−2)/2 = 3(p²−1)/2."""
    n = n_of(p)
    return Fraction(3 * (n - 2), 2)


def theorem_star_identity(primes: list[int] | None = None) -> dict:
    """Max+-free double count: star = 3(n−2)/2 from global ∑μ₄κ."""
    if primes is None:
        primes = [p for p in range(3, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        global_sum = sum_mu4_kappa_global(p)
        # binom(n,2) * star = 6 * global
        star_from_global = 6 * global_sum / Fraction(n * (n - 1), 2)
        star = star_mu4_kappa(p)
        match = star_from_global == star == Fraction(3 * (p * p - 1), 2)
        rows[str(p)] = {
            "global_sum_mu4_kappa": str(global_sum),
            "star": str(star),
            "star_from_double_count": str(star_from_global),
            "equals_3_halves_p2_minus_1": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "∑_S μ₄κ = n(n−1)(n−2)/8; each 4-set hits 6 edges ⇒ "
            "star_e = 3(n−2)/2 = 3(p²−1)/2."
        ),
        "by_p": rows,
    }


def psd_gsum_lb() -> Fraction:
    """Gsum ≽ 0, diag 2 ⇒ every entry ≥ −2."""
    return Fraction(-2)


def dual_eq_box_lower_psd(n_d: int) -> Fraction:
    """(Gsum x)_e ≥ −2 n_d (wedges contribute 0)."""
    return psd_gsum_lb() * n_d


def theorem_nd_le_1_impossible(primes: list[int] | None = None) -> dict:
    """Dual equality impossible when n_d ∈ {0,1} for all primes p≥5."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        need = dual_equality_correlation_need(p)
        kill0 = Fraction(0) > need  # n_d=0
        kill1 = dual_eq_box_lower_psd(1) > need  # −2 > need
        # need = 6/p−4; need < −2 ⇔ 6/p < 2 ⇔ p > 3, true for p≥5
        rows[str(p)] = {
            "need": str(need),
            "lb_nd0": "0",
            "lb_nd1": str(dual_eq_box_lower_psd(1)),
            "kill_nd0": kill0,
            "kill_nd1": kill1,
        }
        if not (kill0 and kill1 and need < -2):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For p≥5, need=6/p−4 < −2. PSD ⇒ Gsum≥−2 so box≥−2n_d. "
            "n_d=0 ⇒ box=0>need; n_d=1 ⇒ box≥−2>need. Dual equality impossible."
        ),
        "by_p": rows,
    }


def nd2_wedge_partners_lb() -> float:
    """Min g1+g2 for 3×3 PSD with middle 0: −2√2."""
    return -2.0 * math.sqrt(2.0)


def theorem_nd2_wedge_kill_p_ge_7(primes: list[int] | None = None) -> dict:
    """
    n_d=2 with wedge partners (Gsum_fg=0): box ≥ −2√2.
    Obstructs dual equality for all primes p≥7.
    """
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    lb = nd2_wedge_partners_lb()
    # Algebra: −2√2 > 6/p−4 ⇔ 4−2√2 > 6/p ⇔ p > 6/(4−2√2)=3/(2−√2)
    thr = 3.0 / (2.0 - math.sqrt(2.0))
    rows = {}
    ok_ge7 = True
    for p in primes:
        need = float(dual_equality_correlation_need(p))
        kills = lb > need
        rows[str(p)] = {
            "need": need,
            "lb_neg_2_sqrt2": lb,
            "kills_dual_eq": kills,
            "expect_kill": p >= 7,
        }
        if p >= 7 and not kills:
            ok_ge7 = False
        if p == 5 and kills:
            ok_ge7 = False  # must NOT kill at p=5
    return {
        "proved_for_primes_ge_7": ok_ge7,
        "threshold_p": thr,
        "lb": lb,
        # p=5: need=−14/5=−2.8 > −2√2≈−2.828 ⇒ wedge-n_d=2 does NOT kill
        "p5_need_gt_lb": float(dual_equality_correlation_need(5)) > lb,
        "theorem": (
            "n_d=2 wedge partners: [[2,g1,g2],[g1,2,0],[g2,0,2]]≽0 ⇒ "
            "g1²+g2²≤4 ⇒ g1+g2≥−2√2. Beats need=6/p−4 iff p>3/(2−√2)≈5.12 "
            "(all primes p≥7). p=5 open for this configuration."
        ),
        "by_p": rows,
    }


def certify_star_census(p: int) -> dict:
    """Census star sum of μ₄κ for one edge (p≤5)."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        raise ValueError("census p<=5")
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    e = (0, 1)
    star = 0.0
    for c, d in combinations([i for i in range(n) if i not in e], 2):
        a, b = e
        kap = (
            C[a, b] * C[c, d]
            + C[a, c] * C[b, d]
            + C[a, d] * C[b, c]
        )
        mp = float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
        mm = float(np.mean(Mm[:, a] * Mm[:, b] * Mm[:, c] * Mm[:, d]))
        mu = 0.5 * (mp + mm)
        star += mu * kap
    expected = float(star_mu4_kappa(p))
    return {
        "p": p,
        "star_empirical": star,
        "star_proved": expected,
        "match": bool(abs(star - expected) < 1e-9),
        "not_general_proof": True,
    }


def hinge_status_178() -> dict:
    return {
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "proved_nd_le_1_impossible": True,
        "proved_nd2_wedge_p_ge_7": True,
        "proved_star_identity": True,
        "bound_proved_general": False,
        "open": (
            "Prove |μ₄|≤1/(2p) on |κ|=1 (or any disj Gsum μ>μ_*) for all "
            "primes p≥5 Max+-free — covers remaining n_d≥2 (and p=5 n_d=2). "
            "Then flip gsum_disj_lb_proved_general."
        ),
        "mu_star_p5": str(mu_star_threshold(5)),
        "mu4_thr_p5": str(mu4_threshold(5)),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    star = theorem_star_identity([p for p in range(3, 80) if is_prime(p)])
    nd1 = theorem_nd_le_1_impossible(primes)
    nd2 = theorem_nd2_wedge_kill_p_ge_7(primes)
    cert3 = certify_star_census(3)
    cert5 = certify_star_census(5)
    # p=5: −2√2 does not beat need
    p5_nd2_open = not (
        nd2_wedge_partners_lb() > float(dual_equality_correlation_need(5))
    )
    out = {
        "title": (
            "Prop 15.178 star identity + dual-eq n_d case kills; "
            "pointwise Gsum LB still OPEN"
        ),
        "proved": {
            "star_identity": star["proved"],
            "nd_le_1_dual_eq_impossible": nd1["proved"],
            "nd2_wedge_dual_eq_impossible_p_ge_7": nd2["proved_for_primes_ge_7"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(primes),
            "mu4_bound_proved_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "star": star,
            "nd_le_1": nd1,
            "nd2_wedge": nd2,
            "p5_nd2_wedge_still_open": p5_nd2_open,
            "sample_farkas_if_minus_1_over_p": {
                str(p): residual_i_farkas_if_mu(p, mu_minus_1_over_p(p))
                for p in (5, 7, 11, 13)
            },
        },
        "certified_census": {"p3": cert3, "p5": cert5},
        "hinge": hinge_status_178(),
        "F3": (
            "Do not flip residual/E1/L without general |μ4|≤1/(2p) or "
            "Gsum≥−1/p (or alternate covering all n_d)."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15178.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.178 star identity + n_d case kills")
    print(f"  star proved={star['proved']} census p5 match={cert5['match']}")
    print(f"  n_d≤1 impossible={nd1['proved']}")
    print(f"  n_d=2 wedge p≥7={nd2['proved_for_primes_ge_7']} (p5 open={p5_nd2_open})")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
