#!/usr/bin/env python3
"""
Prop 15.229 — Size-3 Cy-expansion vanishes on |κ|=1 (Max+-free);
(p⁴−1)μ+2φ = R̄₄; residual-(i) still OPEN on size-4 remainder.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.191, 15.177)
  For y∈Max±: ∏_{i∈S}(Cy)_i = (±p)^4 ∏_{i∈S} y_i = p⁴ χ_S.
  Expand LHS over maps f:S→[n]; group by |im f|.  Take E_Max±.
  On |κ|=1: size1+size2 = −2φ (15.191 E, Max+-free).  Size-4 =: Q₄^±.
  Target: control μ=½(m₄⁺+m₄⁻) so |μ|≤1/(2p) on |κ|=1.

PROVED Max+-free
  A. Unconstrained triple sum vanishes for domain doubletons.
     Fix domain doubleton {i,j}⊂S and singles k,l.  Then
        ∑_{u,v,w∈[n]} C_{iu}C_{ju}C_{kv}C_{lw}C_{vw}
          = (n−1) C_{kl} (C²)_{ij} = 0
     for i≠j, by C²=(n−1)I (off-diagonal zero).  ∎

  B. Collision calculus for size-3 (domain type 2+1+1).
     Size-3 weight (with factor C_{vw} from the Max± moment
     E[y_v y_w]=±C_{vw}/p) equals unconstrained minus collisions.
     C_{vv}=0 kills v=w and u=v=w.  The remaining collisions are
        u=v≠w:  (n−1) C_{il}C_{jl}C_{kl},
        u=w≠v:  (n−1) C_{ik}C_{jk}C_{lk}.
     Hence for each domain doubleton assignment,
        size3_{(i,j;k,l)} = −(n−1) C_{kl}(C_{ik}C_{jk}+C_{il}C_{jl}).  ∎
     (Sign of π=±C/p factors out globally as ±1/p and does not affect vanishing.)

  C. Full size-3 moment weight is −(n−1)·s, with
        s = ∑_{edges {k,l} of K4} C_{kl}(C_{ik}C_{jk}+C_{il}C_{jl})
     where {i,j} is the complementary edge.  ∎

  D. s=0 on every |κ|=1 labelling (64-exhaustion).
     s depends only on the six ±1 edge labels of K4.  Exhaustion of all 64
     zero-diag symmetric labelings: s=0 whenever |κ|=1 (24+24 labelings);
     s=±12 when |κ|=3 (8+8 labelings).  Max+-free pure K4 algebra.  ∎

  E. Size-3 Cy-expansion vanishes on |κ|=1.
     Combine (B)–(D): size-3 contribution to E[∏(Cy)_i] is 0 on every
     |κ|=1 four-set, for both Max+ and Max−.  ∎

  F. Functional equation.
     On |κ|=1, for both Max±:
        p⁴ m₄^± = −2φ + per(C[S,S]) m₄^± + R^±
     where per(C[S,S]) equals the derangement permanent (=1 on |κ|=1 by
     15.191 A: fixed points give C_{ii}=0).  R^± collects size-4 maps with
     image ≠ S.  Hence
        (p⁴ − 1) m₄^± = −2φ + R^±,
        (p⁴ − 1) μ + 2φ = R̄₄ := ½(R⁺+R⁻).  ∎

  G. Particular majorant consequence (conditional).
     |μ| ≤ (2|φ| + |R̄₄|)/(p⁴−1).  Under |φ|≤2(p−2), if
        |R̄₄| ≤ (p⁴−1)/(2p) − 4(p−2)
     then |μ|≤1/(2p).  At p=5 the RHS is 624/10−12=50.4; at p=7:
     2400/14−20≈151.4.  ∎

CERTIFIED
  H. s=0 on |κ|=1 and s=±12 on |κ|=3 for all 64 labelings; collision
     formula matches direct triple loops at p=3,5; size3 weight 0 on
     all |κ|=1 four-sets checked p=3 (full) and p=5 (sample).
  I. Census: p=5 has μ=f4 and satisfies the functional equation; p=7
     max|μ|=109/2863≤1/(2p).

OPEN (blocks residual i / predicate flip)
  J. Bound |R̄₄|≤(p⁴−1)/(2p)−4(p−2) on |κ|=1 for all primes p≥5
     (or otherwise prove |μ|≤1/(2p) / Path C / K₄ / free-e+ker=sc).
     Soft-close forbidden.

Until J: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15229.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import theorem_phi_bound_general  # noqa: E402
from e1_gmin_m4_prop15191 import (  # noqa: E402
    theorem_derangement_perm_on_kappa1,
    theorem_size1_plus_size2,
)


def kappa_of_C6(C6: tuple[int, ...]) -> int:
    c01, c02, c03, c12, c13, c23 = C6
    return c01 * c23 + c02 * c13 + c03 * c12


def inner_s_of_C6(C6: tuple[int, ...]) -> int:
    """s = ∑_{edges e={k,l}} C_e (C_{ik}C_{jk}+C_{il}C_{jl}), complement {i,j}."""
    c01, c02, c03, c12, c13, c23 = C6
    s = 0
    s += c01 * (c02 * c03 + c12 * c13)
    s += c02 * (c01 * c03 + c12 * c23)
    s += c03 * (c01 * c02 + c13 * c23)
    s += c12 * (c01 * c13 + c02 * c23)
    s += c13 * (c01 * c12 + c03 * c23)
    s += c23 * (c02 * c12 + c03 * c13)
    return s


def theorem_s_vanishes_on_kappa1() -> dict:
    """64-exhaustion: s=0 on |κ|=1; s=±12 on |κ|=3."""
    by_k: dict[int, set[int]] = {}
    ok_k1 = True
    for bits in product((-1, 1), repeat=6):
        k = kappa_of_C6(bits)
        s = inner_s_of_C6(bits)
        by_k.setdefault(k, set()).add(s)
        if abs(k) == 1 and s != 0:
            ok_k1 = False
    return {
        "proved": ok_k1 and by_k.get(1) == {0} and by_k.get(-1) == {0},
        "theorem": (
            "For every zero-diag symmetric ±1 labelling of K4: if |κ|=1 then "
            "s=0; if |κ|=3 then s∈{±12}.  Exhaustion of 64 labelings."
        ),
        "by_kappa": {str(k): sorted(v) for k, v in sorted(by_k.items())},
        "depends_on": [],
    }


def theorem_size3_collision_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For domain doubleton {i,j} and singles k,l: the size-3 "
            "C-weight with factor C_{vw} equals "
            "−(n−1)C_{kl}(C_{ik}C_{jk}+C_{il}C_{jl}), by C²=(n−1)I and "
            "C_{vv}=0 killing the v=w collision.  Summing over the six "
            "doubletons: full size-3 weight = −(n−1)·s."
        ),
        "depends_on": ["C2_eq_n_minus_1_I", "C_diag_zero"],
    }


def theorem_size3_vanishes_kappa1() -> dict:
    t_s = theorem_s_vanishes_on_kappa1()
    t_c = theorem_size3_collision_formula()
    return {
        "proved": t_s["proved"] and t_c["proved"],
        "theorem": (
            "On every |κ|=1 four-set of any symmetric conference: the size-3 "
            "contribution to E_Max±[∏(Cy)_i] vanishes.  Proof: collision "
            "formula ⇒ weight=−(n−1)s; s=0 by 64-exhaustion on |κ|=1.  "
            "Both Max+ (π=C/p) and Max− (π=−C/p)."
        ),
        "depends_on": [
            "theorem_s_vanishes_on_kappa1",
            "theorem_size3_collision_formula",
            "prop_15_189_pairwise_pi",
        ],
    }


def R_bar_budget_for_mu_thr(p: int) -> Fraction:
    """(p⁴−1)/(2p) − 4(p−2): |R̄|≤ this ⇒ |μ|≤1/(2p) under |φ|≤2(p−2)."""
    return Fraction(p**4 - 1, 2 * p) - 4 * (p - 2)


def theorem_functional_equation() -> dict:
    t3 = theorem_size3_vanishes_kappa1()
    t12 = theorem_size1_plus_size2()
    tper = theorem_derangement_perm_on_kappa1()
    return {
        "proved": t3["proved"] and t12["proved"] and tper["proved"],
        "theorem": (
            "On |κ|=1: p⁴ m₄^± = −2φ + m₄^± + R^± with R^± the size-4 "
            "non-self-image contribution, hence (p⁴−1)μ + 2φ = R̄₄.  "
            "Uses size1+2=−2φ, size3=0, derangement permanent=1."
        ),
        "equation": "(p^4 - 1) * mu + 2 * phi = R_bar_4",
        "depends_on": [
            "theorem_size3_vanishes_kappa1",
            "prop_15_191_size1_plus_size2",
            "prop_15_191_derangement_permanent",
        ],
    }


def theorem_conditional_R_bound_closes() -> dict:
    ok = True
    sample = {}
    phi_ok = theorem_phi_bound_general()["proved"]
    for p in range(5, 201):
        if not is_prime(p):
            continue
        bud = R_bar_budget_for_mu_thr(p)
        # need bud ≥ 0 so the implication is usable
        if bud < 0:
            ok = False
            break
        # check implication algebra: (2*2(p-2) + bud)/(p^4-1) ≤ 1/(2p)
        lhs_num = 4 * (p - 2) + bud
        # should equal (p^4-1)/(2p)
        if lhs_num != Fraction(p**4 - 1, 2 * p):
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31):
            sample[str(p)] = {
                "R_bar_budget": str(bud),
                "implies_mu_le_thr": True,
            }
    return {
        "proved_implication": ok and phi_ok,
        "bound_proved_general": False,
        "theorem": (
            "If |R̄₄| ≤ (p⁴−1)/(2p)−4(p−2) on every |κ|=1 four-set for all "
            "primes p≥5 (and |φ|≤2(p−2)), then |μ|≤1/(2p) ⇒ residual-(i) "
            "Farkas.  The R̄₄ bound is OPEN."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_functional_equation",
            "prop_15_186_phi_bound",
        ],
    }


def residual_i_closed_via_229() -> bool:
    return False


def hinge_status_229() -> dict:
    return {
        "s_vanishes_kappa1": theorem_s_vanishes_on_kappa1()["proved"],
        "size3_collision": theorem_size3_collision_formula()["proved"],
        "size3_vanishes": theorem_size3_vanishes_kappa1()["proved"],
        "functional_equation": theorem_functional_equation()["proved"],
        "conditional_R_closes": theorem_conditional_R_bound_closes()[
            "proved_implication"
        ],
        "R_bound_general": theorem_conditional_R_bound_closes()[
            "bound_proved_general"
        ],
        "residual_i_closed_via_229": residual_i_closed_via_229(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤(p⁴−1)/(2p)−4(p−2) on |κ|=1; "
            "or Aut-SOS/Path C/K4/free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_s_vanishes_on_kappa1()
    b = theorem_size3_collision_formula()
    c = theorem_size3_vanishes_kappa1()
    d = theorem_functional_equation()
    e = theorem_conditional_R_bound_closes()
    status = hinge_status_229()
    out = {
        "prop": "15.229",
        "title": (
            "Size-3 Cy vanishes on |κ|=1; (p⁴−1)μ+2φ=R̄₄; "
            "R̄₄ bound OPEN"
        ),
        "proved": {
            "s_vanishes_kappa1": a["proved"],
            "size3_collision_formula": b["proved"],
            "size3_vanishes_kappa1": c["proved"],
            "functional_equation": d["proved"],
            "conditional_R_closes": e["proved_implication"],
            "R_bound_general": e["bound_proved_general"],
            "residual_i_closed": residual_i_closed_via_229(),
        },
        "algebra": {
            "size3_weight": "-(n-1)*s",
            "s_on_kappa1": 0,
            "s_on_kappa3": "±12",
            "functional_eq": "(p^4-1)*mu + 2*phi = R_bar_4",
            "R_bar_budget_p5": str(R_bar_budget_for_mu_thr(5)),
            "R_bar_budget_p7": str(R_bar_budget_for_mu_thr(7)),
            "mu4_thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_229(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Size-3 vanishes Max+-free (C² + 64-exhaust).  "
            "Functional equation reduces residual-(i) to an R̄₄ bound.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15229.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.229  residual_i={residual_i_closed_via_229()}")
    print(
        f"  s_vanish={a['proved']} size3=0={c['proved']} "
        f"func_eq={d['proved']} R_bound={e['bound_proved_general']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
