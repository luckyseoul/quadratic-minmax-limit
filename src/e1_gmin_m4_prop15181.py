#!/usr/bin/env python3
"""
Prop 15.181 — Vertex-star constraints on Gsum; κ-counts through an edge;
residual-(i) open core refined. Hinge still OPEN.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

SETUP
  Conference matrix C of order n=p²+1: C=Cᵀ, C_ii=0, C_ij=±1 (i≠j),
  C²=(n−1)I, eigenvalues ±p. Max± = boolean ±p-eigenvectors.
  Edge features f_e = C_e y_i y_j for e={i,j}.
  Gsum = E₊[f fᵀ] + E₋[f fᵀ] (15.170).

PROVED Max+-free (any conference C, Cy=±p y with y_i=±1):
  A. Vertex star of features. For every vertex i and every Max+ vector,
        ∑_{e∋i} f_e = p.
     Proof. ∑_{j≠i} C_ij y_i y_j = y_i (Cy)_i = y_i (p y_i) = p.
     Similarly Max−: ∑_{e∋i} f_e = −p.
  B. Gsum star-sum. For every vertex i and every edge f,
        ∑_{e∋i} Gsum_{e f} = 2.
     Proof. E₊[f_f · p] + E₋[f_f · (−p)] = p E₊[f_f] − p E₋[f_f]
     = p(1/p) − p(−1/p) = 2  (using E±[f_f]=±1/p from E[y_a y_b]=±C_ab/p).
  C. Consequence for f not incident to i: the n−3 quantities
        h_x := Gsum({i,x}, f),  x∉{i}∪f,
     sum to 2. Principal block on those star edges is 2I (pairwise wedges),
     so [[2I,h],[hᵀ,2]]≽0 ⇒ ‖h‖²≤4 (Schur). Entrywise |h_x|≤2.
  D. κ-counts through a fixed edge e (C² algebra only).
     Let D=binom(n−2,2) = # of 4-sets containing e.
     Write n₁, n₃ for how many have |κ|∈{1,3}.
     Then ∑ κ² = n₁ + 9 n₃ = 3(n−2)(n−5)/2, and n₁+n₃=D, so
        n₃ = (n−2)(n−6)/8,   n₁ = 3(n−2)²/8.
     Proof of ∑κ²: expand κ²=3+2(cross); the three cross double-counts
     equal −2(n−2), −(n−2) respectively via (C²)_{uv}=0 (details in
     theorem_kappa_counts_through_edge). ∎
  E. Farkas algebra unchanged (15.176–178): Gsum≥−1/p on disj still
     kills dual equality for all n_d; n_d≤1 and n_d=2 wedge (p≥7) already
     killed by PSD alone.

CERTIFIED (not general):
  F. p=5: min_x (Gsum x)_e over 0≤x≤1, 1ᵀx=k, x_e=0 equals −1.2 > need=−2.8
     (actual Gsum row; dual-eq linear box empty). Matches pointwise min −6/65.
  G. Matching n_d=2 free 3×3 PSD allows g1+g2≥−4 (worse than wedge −2√2);
     matching-heavy is the hard PSD case.

OPEN (blocks residual i / E1 / L):
  H. Prove disj Gsum≥−1/p or |μ₄|≤1/(2p) on |κ|=1 (or dual-eq impossible
     on all remaining n_d≥2 matching-heavy supports, including p=5 n_d=2)
     Max+-free for all primes p≥5.
  Vertex-star + κ-counts refine the frame but do not yet force H.

Until H: gsum_disj_lb_proved_general False; type_I closed False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15181.json
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
    nd2_wedge_partners_lb,
    theorem_nd2_wedge_kill_p_ge_7,
    theorem_nd_le_1_impossible,
)


def vertex_star_feature_identity() -> dict:
    """
    Max+-free: ∑_{e∋i} f_e = p on Max+, = −p on Max−.
    From Cy=±p y and y_i²=1.
    """
    return {
        "proved": True,
        "max_plus": "sum_{e ni i} f_e = p",
        "max_minus": "sum_{e ni i} f_e = -p",
        "proof": (
            "sum_{j≠i} C_ij y_i y_j = y_i (C y)_i = y_i (λ y_i) = λ "
            "with λ=±p on Max±."
        ),
    }


def gsum_star_sum_identity() -> dict:
    """Max+-free: ∑_{e∋i} Gsum_{e f} = 2 for every vertex i and edge f."""
    return {
        "proved": True,
        "identity": "sum_{e ni i} Gsum_{e,f} = 2  for all vertices i, edges f",
        "proof": (
            "E+[f_f * p] + E-[f_f * (-p)] = p/p - p*(-1/p) = 2, "
            "using E±[f_f]=±1/p."
        ),
    }


def star_schur_on_pending(p: int) -> dict:
    """
    For f not at i: n−3 pending Gsums h_x sum to 2 and ‖h‖²≤4.
    Single-coordinate lower bound (equal remainder): root of
      L² + (2−L)²/(m−1) ≤ 4, m=n−3.
    Too weak for residual i (≫ need for small L target −1/p).
    """
    n = n_of(p)
    m = n - 3
    # binary search L minimal with L²+(2-L)²/(m-1) ≤ 4, L∈[-2,2]
    lo, hi = -2.0, 2.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if mid * mid + (2 - mid) ** 2 / (m - 1) <= 4.0 + 1e-15:
            hi = mid  # feasible, try smaller
        else:
            lo = mid
    L = hi
    need = float(dual_equality_correlation_need(p))
    return {
        "p": p,
        "m": m,
        "sum_h": 2,
        "norm_h_sq_le": 4,
        "single_coord_lb_approx": L,
        "need": need,
        "single_coord_beats_need": L > need,
        "note": "Single pending Gsum LB from one vertex star; not residual-i close.",
    }


def n3_through_edge(p: int) -> Fraction:
    """n₃ = (n−2)(n−6)/8 for 4-sets through a fixed edge with |κ|=3."""
    n = n_of(p)
    return Fraction((n - 2) * (n - 6), 8)


def n1_through_edge(p: int) -> Fraction:
    """n₁ = 3(n−2)²/8 for 4-sets through a fixed edge with |κ|=1."""
    n = n_of(p)
    return Fraction(3 * (n - 2) * (n - 2), 8)


def sum_kappa_sq_through_edge(p: int) -> Fraction:
    """∑ κ² over 4-sets through e = 3(n−2)(n−5)/2."""
    n = n_of(p)
    return Fraction(3 * (n - 2) * (n - 5), 2)


def theorem_kappa_counts_through_edge(primes: list[int] | None = None) -> dict:
    """
    Max+-free C² algebra: n₃=(n−2)(n−6)/8, n₁=3(n−2)²/8,
    ∑κ²=3(n−2)(n−5)/2 for 4-sets containing a fixed edge.
    """
    if primes is None:
        primes = [p for p in range(3, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        D = Fraction((n - 2) * (n - 3), 2)
        s2 = sum_kappa_sq_through_edge(p)
        # cross total = ∑κ² − 3D = −3(n−2)
        cross = s2 - 3 * D
        n3 = (s2 - D) / 8
        n1 = D - n3
        match = (
            n3 == n3_through_edge(p)
            and n1 == n1_through_edge(p)
            and cross == Fraction(-3 * (n - 2))
            and n1 + n3 == D
            and n1 + 9 * n3 == s2
        )
        rows[str(p)] = {
            "D": str(D),
            "sum_kappa_sq": str(s2),
            "cross_total": str(cross),
            "n1": str(n1),
            "n3": str(n3),
            "match_closed_forms": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Fix edge e={a,b}, W=V\\e, D=binom(|W|,2). For κ on 4-sets e∪{c,d}: "
            "κ²=3+2(cross12+cross13+cross23). Using (C²)_{uv}=0 for u≠v: "
            "sum_W C_ac C_bc=0 ⇒ s23=sum_{c<d} 2(C_ac C_bc)(C_ad C_bd)=−(n−2); "
            "sum_{c≠d} C_cd C_ac C_bd=−C_ab(n−2) ⇒ s12+s13=−2(n−2). "
            "Cross total −3(n−2). Hence ∑κ²=3D−3(n−2)=3(n−2)(n−5)/2. "
            "With n1+n3=D and n1+9n3=∑κ²: n3=(n−2)(n−6)/8, n1=3(n−2)²/8."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_primes": len(primes),
        "all_match": ok,
    }


def theorem_matching_nd2_psd_lb() -> dict:
    """
    Free 3×3 PSD on {e,f,g} with all three pairs disj (matching partners):
    min g1+g2 = −4 (achieved at g1=g2=−2, Gsum_fg=2).
    Strictly weaker than wedge (Gsum_fg=0 ⇒ −2√2). Matching is the hard case.
    """
    # Check eigenvalues at the claimed minimizer
    import numpy as np

    a = b = -2.0
    c = 2.0
    M = np.array([[2.0, a, b], [a, 2.0, c], [b, c, 2.0]])
    evals = np.linalg.eigvalsh(M)
    return {
        "proved": bool(evals.min() >= -1e-12),
        "min_g1_plus_g2": -4.0,
        "achieved_at": {"g1": -2.0, "g2": -2.0, "g_fg": 2.0},
        "wedge_lb": nd2_wedge_partners_lb(),
        "matching_weaker_than_wedge": -4.0 < nd2_wedge_partners_lb(),
        "theorem": (
            "Matching n_d=2: [[2,g1,g2],[g1,2,c],[g2,c,2]]≽0 allows "
            "g1+g2=−4 (c=2). Wedge c=0 forces g1+g2≥−2√2. "
            "PSD alone does not kill matching n_d=2 at p=5 (need=−2.8∈(−4,−2√2))."
        ),
    }


def theorem_nd2_schur_all_configs(primes: list[int] | None = None) -> dict:
    """
    Record: free Schur with B=2I (all partner pairs wedge — impossible for
    simple graphs at n_d=2 except the wedge-partner geometry already treated)
    vs matching free B.
    Does not claim a new kill beyond 15.178.
    """
    if primes is None:
        primes = [p for p in range(5, 40) if is_prime(p)]
    rows = {}
    for p in primes:
        need = float(dual_equality_correlation_need(p))
        rows[str(p)] = {
            "need": need,
            "wedge_lb_neg2sqrt2": nd2_wedge_partners_lb(),
            "matching_lb_neg4": -4.0,
            "wedge_kills": nd2_wedge_partners_lb() > need,
            "matching_psd_kills": -4.0 > need,  # never for p≥5
        }
    return {
        "proved_new_kill": False,
        "note": (
            "Matching n_d=2 PSD floor −4 never beats need=6/p−4>−4. "
            "Wedge kill for p≥7 remains the only n_d=2 PSD kill (15.178)."
        ),
        "by_p": rows,
    }


def certify_ksparse_box_p5() -> dict:
    """
    Census p=5: LP min (Gsum x)_e over 0≤x≤1, sum x=k, x_e=0 is −6/5 > need.
    Dual-equality linear constraint empty at p=5. Not a general proof.
    """
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    edges = list(combinations(range(n), 2))
    E = len(edges)
    e0 = (0, 1)
    i0 = edges.index(e0)

    def feats(M):
        F = np.zeros((M.shape[0], E))
        for j, (a, b) in enumerate(edges):
            F[:, j] = C[a, b] * M[:, a] * M[:, b]
        return F

    G = (feats(Mp).T @ feats(Mp)) / len(Mp) + (feats(Mm).T @ feats(Mm)) / len(
        Mm
    )
    row = G[i0]
    k = k_type_I_3p_minus_2(p)
    need = float(dual_equality_correlation_need(p))
    # Greedy: most negative k entries among non-e0 (wedges are 0)
    order = np.sort(row)
    # exclude diag contribution: row[i0]=2, but x_e0=0 so drop index i0
    mask = np.ones(E, dtype=bool)
    mask[i0] = False
    vals = np.sort(row[mask])
    lp_min = float(vals[:k].sum())  # put mass on k most negative
    return {
        "p": 5,
        "k": k,
        "need": need,
        "lp_min_approx_greedy": lp_min,
        "need_lt_lp_min": need < lp_min,
        "dual_eq_linear_empty_p5": need < lp_min,
        "pointwise_min_gsum": float(vals[0]),
        "certified_minus_6_over_65": abs(vals[0] + 6 / 65) < 1e-9,
        "not_general_proof": True,
    }


def hinge_status_181() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "vertex_star_proved": True,
        "kappa_counts_proved": True,
        "nd_le_1_impossible": True,
        "nd2_wedge_p_ge_7": True,
        "matching_nd2_psd_only_neg4": True,
        "bound_proved_general": False,
        "open": (
            "Prove disj Gsum≥−1/p or |μ₄|≤1/(2p) on |κ|=1, or dual-eq kill "
            "on all matching-heavy n_d≥2 (esp. p=5 n_d=2), Max+-free p≥5."
        ),
        "mu_star_p5": str(mu_star_threshold(5)),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    star_f = vertex_star_feature_identity()
    star_g = gsum_star_sum_identity()
    kap = theorem_kappa_counts_through_edge(
        [p for p in range(3, 80) if is_prime(p)]
    )
    match = theorem_matching_nd2_psd_lb()
    nd2s = theorem_nd2_schur_all_configs([p for p in range(5, 40) if is_prime(p)])
    nd1 = theorem_nd_le_1_impossible(primes)
    nd2w = theorem_nd2_wedge_kill_p_ge_7(primes)
    cert5 = certify_ksparse_box_p5()
    schur_sample = {str(p): star_schur_on_pending(p) for p in (5, 7, 11, 13)}
    out = {
        "title": (
            "Prop 15.181 vertex-star + κ-counts through edge; "
            "residual (i) matching-heavy core still OPEN"
        ),
        "proved": {
            "vertex_star_features": star_f["proved"],
            "gsum_star_sum": star_g["proved"],
            "kappa_counts_through_edge": kap["proved"],
            "matching_nd2_psd_floor_neg4": match["proved"],
            "nd_le_1_impossible": nd1["proved"],
            "nd2_wedge_p_ge_7": nd2w["proved_for_primes_ge_7"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(primes),
            "mu4_or_gsum_bound_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "vertex_star_features": star_f,
            "gsum_star_sum": star_g,
            "kappa_counts": kap,
            "matching_nd2_psd": match,
            "nd2_config_record": nd2s,
            "star_schur_pending_sample": schur_sample,
            "nd_le_1": nd1,
            "nd2_wedge": nd2w,
        },
        "certified_census": {"p5_ksparse_box": cert5},
        "hinge": hinge_status_181(),
        "F3": (
            "Do not flip residual/E1/L without general |μ4|≤1/(2p) or "
            "Gsum≥−1/p or alternate covering all n_d≥2 matching-heavy."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15181.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.181 vertex-star + κ-counts; residual (i) OPEN")
    print(f"  vertex star proved={star_f['proved']}; gsum star={star_g['proved']}")
    print(f"  kappa counts proved={kap['proved']}")
    print(f"  matching nd2 PSD floor -4 proved={match['proved']}")
    print(f"  p5 dual-eq linear empty={cert5['dual_eq_linear_empty_p5']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
