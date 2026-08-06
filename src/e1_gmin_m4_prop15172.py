#!/usr/bin/env python3
"""
Prop 15.172 — Denseness-path remainder: Max+-free Gsum structure + Farkas threshold.

Does **not** soft-close L. Ships non-theater math progress on the Gsum hinge.

PROVED (Fraction / linear algebra, Max+-free):
  A. Gsum := E_+[f fᵀ] + E_−[f fᵀ] on edges, f_e = C_ij y_i y_j.
  B. Gsum_ee = 2; Gsum_ab = 0 on wedges (share a vertex); Gsum 1 = n 1.
  C. #disj neighbors of e is (n−2)(n−3)/2; sum of Gsum over them = n−2;
     hence average disj entry = 2/(n−3).
  D. Unique triangular-scheme matrix G0 with (B)+(C) has constant disj entry
     γ=2/(n−3) and eigenvalues {n (×1), 0 (×n−1), 2(n−2)/(n−3) (×N−n)};
     so G0 ≽ 0. Actual Gsum = G0 + mean-zero C-refinement (15.158: not a scheme).
  E. Matching PSD: for any matching M of size r≥2,
     average_{pairs in M} Gsum ≥ −2/(r−1).
  F. **Correct dual-equality Farkas threshold.** Dual equality forces
     (Gsum x)_e = 6/p−4 with x=1_G (includes diagonal 2), i.e.
     sum_{e'∈G\\{e}} Gsum_ee' = 6/p−6.
     Off-diag entries ≥ μ (wedges are 0 ≥ μ when μ≤0) ⇒ sum ≥ (k−1)μ
     with k=3p−2. Obstruction need_off < (k−1)μ rearranges to the
     p-independent form **μ > −2/p**. Equality μ=−2/p makes need_off=(k−1)μ
     (borderline; needs saturation analysis).
  G. Candidate scheme LB −12/(p n) is **false at p=3** as a universal-prime
     bound (census min = −2/3 = −2/p < −12/(3·10)=−2/5). So scheme-min cannot
     be the general proof (also blocked by 15.158).

CERTIFIED (finite Max+ census, not general proof):
  H. p=3: disj min Gsum = −2/p (tight for the −2/p threshold).
  I. p=5: disj min Gsum = −6/65 = −12/(p n) > −2/p; Farkas strict under either LB.

OPEN (hinge remains):
  J. Prove for all primes p≥5 that disj Gsum_ab ≥ −2/p (or any μ > −2/p),
     Max+-free — **or** alternate residual-(i)/(ii) avoiding dual-equality Farkas.
  Until then: gsum_disj_lb_proved_general() stays False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15172.json
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
    gsum_disj_lb_proved_general,
    gsum_pairwise_lower_bound,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)


def avg_disj_gsum_identity(p: int) -> dict:
    """Max+-free: avg disj Gsum = 2/(n−3) from row sum + wedge zeros."""
    n = n_of(p)
    n_disj = (n - 2) * (n - 3) // 2
    row_off = n - 2  # sum over all e'≠e equals n−2; wedges contribute 0
    avg = Fraction(row_off, n_disj)  # = 2/(n-3)
    assert avg == Fraction(2, n - 3)
    return {
        "p": p,
        "n": n,
        "n_disj_neighbors": n_disj,
        "sum_disj_Gsum": row_off,
        "avg_disj": str(avg),
        "proved": True,
        "theorem": (
            "Gsum_ee=2, Gsum1=n1, wedge Gsum=0 ⇒ sum over disj neighbors = n−2; "
            "avg = 2/(n−3)."
        ),
    }


def triangular_G0_spectrum(p: int) -> dict:
    """
    G0 = 2I + γ A_disj, γ=2/(n−3), is the unique Johnson J(n,2) scheme matrix
    with diag 2, wedge 0, row sum n. Eigenvalues from triangular graph T(n).
    """
    n = n_of(p)
    N = n * (n - 1) // 2
    gamma = Fraction(2, n - 3)
    # A_disj eigs: deg_disj (×1); −(n−3) (×n−1); +1 (×N−n)
    # (from A_wedge eigs deg_w, n−4, −2 and A_disj = J−I−A_wedge on 1^⊥)
    deg_d = (n - 2) * (n - 3) // 2
    eigs = [
        (2 + gamma * deg_d, 1),  # = n
        (2 + gamma * (-(n - 3)), n - 1),  # = 0
        (2 + gamma * 1, N - n),  # = 2(n−2)/(n−3)
    ]
    assert eigs[0][0] == n
    assert eigs[1][0] == 0
    assert eigs[2][0] == Fraction(2 * (n - 2), n - 3)
    psd = all(val >= 0 for val, _ in eigs)
    return {
        "p": p,
        "n": n,
        "N_edges": N,
        "gamma": str(gamma),
        "eigenvalues": [{"value": str(v), "mult": m} for v, m in eigs],
        "G0_PSD": psd,
        "proved": psd,
        "theorem": (
            "Triangular scheme matrix with diag 2, wedge 0, row sum n has "
            "constant disj entry 2/(n−3) and spectrum {n, 0, 2(n−2)/(n−3)}; G0≽0."
        ),
    }


def matching_average_lb(r: int) -> Fraction:
    """For matching of size r≥2: average Gsum on pairs ≥ −2/(r−1) (PSD)."""
    if r < 2:
        raise ValueError("r≥2")
    return Fraction(-2, r - 1)


def matching_psd_bound(p: int) -> dict:
    n = n_of(p)
    r_max = n // 2
    lb = matching_average_lb(r_max)
    thr = farkas_threshold_mu(p)
    return {
        "p": p,
        "max_matching_size": r_max,
        "avg_Gsum_on_max_matching_LB": str(lb),
        "farkas_threshold_mu": str(thr),
        "max_matching_avg_beats_threshold": lb > thr,
        "proved": True,
        "note": (
            "Matching bound controls averages inside matchings, not pointwise "
            "min over all disj pairs. Pointwise still open."
        ),
    }


def farkas_threshold_mu(p: int) -> Fraction:
    """
    Dual equality impossible if every off-diagonal Gsum ≥ μ with μ > −2/p.

    Derivation: need_off = 6/p − 6 = (k−1)μ_eq with μ_eq = −2/p, k=3p−2.
    """
    return Fraction(-2, p)


def dual_equality_farkas_if_mu(p: int, mu: Fraction) -> dict:
    """Obstruction under pairwise off-diag LB μ (includes wedges ≥ μ if μ≤0)."""
    k = k_type_I_3p_minus_2(p)
    need_total = dual_equality_correlation_need(p)
    need_off = need_total - 2  # remove diagonal Gsum_ee=2
    box = (k - 1) * mu
    obstruct = need_off < box
    thr = farkas_threshold_mu(p)
    return {
        "p": p,
        "k": k,
        "need_total": str(need_total),
        "need_off": str(need_off),
        "mu": str(mu),
        "box_sum_LB_off": str(box),
        "need_off_lt_box": obstruct,
        "threshold_mu": str(thr),
        "mu_gt_threshold": mu > thr,
        "borderline_eq": mu == thr,
    }


def candidate_lb_false_at_p3() -> dict:
    """Census: at p=3, min disj Gsum = −2/3 < candidate −12/(pn)=−2/5."""
    p = 3
    n = n_of(p)
    cand = gsum_pairwise_lower_bound(p)
    # Exact census value for Paley n=10 (full Max± enum is tiny / deterministic)
    min_census = Fraction(-2, 3)  # certified below via boolean_evecs when run
    return {
        "p": p,
        "n": n,
        "candidate_minus_12_over_pn": str(cand),
        "census_min_disj": str(min_census),
        "candidate_false_at_p3": min_census < cand,
        "equals_minus_2_over_p": min_census == Fraction(-2, p),
        "theorem": (
            "Candidate disj LB −12/(pn) fails at p=3 (min=−2/p). "
            "Cannot be a universal-all-primes scheme identity."
        ),
    }


def certify_gsum_min_small_p(p: int) -> dict:
    """Finite Max± census of min disj Gsum (p=3,5 only in tests; not general)."""
    import numpy as np
    from itertools import combinations

    from minmax_quadratic import paley_conference_prime_power
    from e1_bitight_infeas import boolean_evecs

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    if p <= 5:
        Mp = boolean_evecs(C, float(p), 0)
        Mm = boolean_evecs(C, -float(p), 0)
    else:
        raise ValueError("certify_gsum_min_small_p only for p<=5 in this module")
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)

    def G_of(M: np.ndarray) -> np.ndarray:
        f = np.array(
            [[C[a, b] * y[a] * y[b] for a, b in edges] for y in M],
            dtype=np.float64,
        )
        return (f.T @ f) / len(M)

    Gsum = G_of(Mp) + G_of(Mm)
    disj_vals = []
    for a, b in combinations(range(E), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 4:
            disj_vals.append(float(Gsum[a, b]))
    mn = min(disj_vals)
    thr = float(farkas_threshold_mu(p))
    cand = float(gsum_pairwise_lower_bound(p))
    # exact fractions for known p
    exact = {
        3: Fraction(-2, 3),
        5: Fraction(-6, 65),
    }.get(p)
    if exact is not None:
        assert abs(mn - float(exact)) < 1e-12
        mn_frac = exact
    else:
        mn_frac = Fraction(mn).limit_denominator(100000)
    return {
        "p": p,
        "n": n,
        "n_Max_plus": int(len(Mp)),
        "n_Max_minus": int(len(Mm)),
        "min_disj_Gsum": str(mn_frac),
        "minus_2_over_p": str(farkas_threshold_mu(p)),
        "candidate_minus_12_over_pn": str(gsum_pairwise_lower_bound(p)),
        "min_ge_minus_2_over_p": mn >= thr - 1e-12,
        "min_gt_minus_2_over_p": mn > thr + 1e-12,
        "min_eq_candidate": abs(mn - cand) < 1e-12,
        "certified_census_only": True,
        "not_general_proof": True,
    }


def hinge_status() -> dict:
    """Honest status of denseness-path residual hinge after 15.172."""
    thr_ok = all(
        dual_equality_farkas_if_mu(p, Fraction(-2, p) + Fraction(1, 10**9))[
            "need_off_lt_box"
        ]
        for p in (5, 7, 11, 13, 17, 19, 23)
    )
    # exact equality case: mu=-2/p does NOT give strict obstruction
    border = dual_equality_farkas_if_mu(5, Fraction(-2, 5))
    return {
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "farkas_threshold_is_minus_2_over_p": True,
        "strict_mu_above_threshold_obstructs_sample": thr_ok,
        "borderline_mu_eq_threshold_p5": border,
        "candidate_scheme_lb_valid_all_primes": False,  # fails p=3
        "open": (
            "Prove disj Gsum ≥ −2/p (strict for p≥5) Max+-free, or alternate "
            "residual (i)/(ii). E1/L remain OPEN."
        ),
        "progress": [
            "avg disj = 2/(n−3) and G0 spectrum proved Max+-free",
            "Farkas threshold simplified to μ > −2/p",
            "candidate −12/(pn) falsified as universal at p=3",
            "p=3,5 census mins certified",
        ],
    }


def main() -> dict:
    rows_G0 = {str(p): triangular_G0_spectrum(p) for p in (3, 5, 7, 11)}
    rows_avg = {str(p): avg_disj_gsum_identity(p) for p in (3, 5, 7, 11)}
    rows_match = {str(p): matching_psd_bound(p) for p in (3, 5, 7)}
    thr_algebra = {
        str(p): dual_equality_farkas_if_mu(p, Fraction(-2, p) + Fraction(1, p * 1000))
        for p in (5, 7, 11, 13, 17)
    }
    # census
    cert3 = certify_gsum_min_small_p(3)
    cert5 = certify_gsum_min_small_p(5)
    out = {
        "title": (
            "Prop 15.172 Gsum structure + Farkas threshold −2/p "
            "(denseness remainder; L still OPEN)"
        ),
        "proved": {
            "avg_disj_identity": all(r["proved"] for r in rows_avg.values()),
            "G0_spectrum_PSD": all(r["G0_PSD"] for r in rows_G0.values()),
            "matching_psd_average_lb": True,
            "farkas_threshold_mu_eq_minus_2_over_p": True,
            "candidate_lb_false_at_p3": candidate_lb_false_at_p3()[
                "candidate_false_at_p3"
            ],
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed_general": False,
            "L_closed": False,
        },
        "algebra": {
            "avg_disj": rows_avg,
            "G0": rows_G0,
            "matching": rows_match,
            "farkas_if_strict_above_minus_2_over_p": thr_algebra,
            "candidate_p3": candidate_lb_false_at_p3(),
        },
        "certified_census": {"p3": cert3, "p5": cert5},
        "hinge": hinge_status(),
        "F3": "Do not flip gsum_disj_lb_proved_general without Max+-free general LB",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15172.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.172 denseness remainder (Gsum structure)")
    print(f"  avg disj identity: {out['proved']['avg_disj_identity']}")
    print(f"  G0 PSD: {out['proved']['G0_spectrum_PSD']}")
    print(f"  Farkas threshold μ > −2/p: {out['proved']['farkas_threshold_mu_eq_minus_2_over_p']}")
    print(f"  candidate false at p=3: {out['proved']['candidate_lb_false_at_p3']}")
    print(f"  census p3 min={cert3['min_disj_Gsum']} p5 min={cert5['min_disj_Gsum']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()} (still False)")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
