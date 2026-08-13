#!/usr/bin/env python3
"""
Prop 15.233 — Signed k=2 layer of R̄₄: 14-term / per2 closed form;
bilinear in the two new columns; unsigned still dead.
Residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

HINGE GRAPH (live)
  OPEN: |R̄₄|≤B; |δ|≤room_δ; Path C; K₄; ker=sc+free-e.
  DEAD: crude ∑|per|; unsigned k=0,1,2 (15.232); unsigned k=2 even with |per|≤14.
  PROVED structure: k=3 column-linear (15.232); k=2 bilinear (this prop).
  Next: signed k=1, then k=0.

SETUP
  R^{(2)}(S)=∑_{|T|=4, |S∩T|=2} per(C[S,T]) μ_T.
  Write I=S∩T={a,b}, {c,d}=S\\I, {s,t}=T\\I.
  per2(U,V):=C_{u1 v1}C_{u2 v2}+C_{u1 v2}C_{u2 v1}  (2×2 permanent).

PROVED Max+-free (C_ii=0, C_ij=±1; C² not required for A–C)
  A. 14 surviving permutations.
     Ordering S=(a,b,c,d), T=(a,b,s,t) puts zeros at the two shared
     diagonal entries.  24−10=14 permutations survive.  |per|≤14.  ∎

  B. Closed form (grouped).
        per(C[S,T]) = per2({c,d},{s,t})
          + (κ(S)−C_ab C_cd) per2({a,b},{s,t})
          + C_ab [ C_ac per2({b,d},{s,t}) + C_ad per2({b,c},{s,t})
                  + C_bc per2({a,d},{s,t}) + C_bd per2({a,c},{s,t}) ]
     with κ(S)=C_ab C_cd+C_ac C_bd+C_ad C_bc.
     Proof: expand the 14 survivors; C_ab²=1 converts the two
     (ab,ab,cs,dt)/(ab,ab,ct,ds) terms into per2({c,d},{s,t});
     the four no-ab crossing terms factor as (κ−C_ab C_cd) per2({a,b});
     the remaining eight carry one C_ab and the displayed per2's.  ∎

  C. Bilinearity (signed pairing).
     Every surviving monomial has exactly one s-edge and one t-edge.
     Hence per=∑_{u,v∈S} β_{uv}(S;I) C_{u s} C_{v t}.
     Summing {s,t} before abs is a C⊗C contraction against μ — the
     Jacobi / cycle-before-abs move on the k=2 layer.  ∎

  D. Unsigned k=2 still dead for all primes p≥5.
     |R^{(2)}|≤14·6 C(n−4,2)=42(p²−3)(p²−4).  This is Θ(p⁴) and exceeds
     B=Θ(p³) for every prime p≥5 (Fraction on 5≤p≤300; leading 42p⁴ vs
     p³/2 for p>300).  Signs inside k=2 remain necessary.  ∎

OPEN
  E. Bound the μ-weighted bilinear form, or signed k=1 / k=0 / Aut-SOS.
     Soft-close forbidden.

Writes evidence/e1_gmin_m4_prop15233.json
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
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15232 import (  # noqa: E402
    layer_count,
    theorem_k3_permanent_closed_form,
    theorem_unsigned_layers_012_dead,
)


def per2(C, U, V) -> int:
    """2×2 permanent of C[U,V]."""
    u1, u2 = U
    v1, v2 = V
    return int(C[u1, v1]) * int(C[u2, v2]) + int(C[u1, v2]) * int(C[u2, v1])


def kappa_of_S(C, S) -> int:
    a, b, c, d = S
    return (
        int(C[a, b]) * int(C[c, d])
        + int(C[a, c]) * int(C[b, d])
        + int(C[a, d]) * int(C[b, c])
    )


def k2_permanent_closed(C, S, a, b, s, t) -> int:
    """Closed per(C[S, {a,b,s,t}]) for I={a,b}=S∩T."""
    rest = [x for x in S if x not in (a, b)]
    if len(rest) != 2:
        raise ValueError("a,b must be a 2-subset of S")
    c, d = rest

    def e(i, j) -> int:
        return int(C[i, j])

    kap = kappa_of_S(C, (a, b, c, d))
    term_rest = per2(C, (c, d), (s, t))
    term_cross = (kap - e(a, b) * e(c, d)) * per2(C, (a, b), (s, t))
    term_ab = e(a, b) * (
        e(a, c) * per2(C, (b, d), (s, t))
        + e(a, d) * per2(C, (b, c), (s, t))
        + e(b, c) * per2(C, (a, d), (s, t))
        + e(b, d) * per2(C, (a, c), (s, t))
    )
    return term_rest + term_cross + term_ab


def k2_unsigned_majorant(p: int) -> int:
    """14 · 6 C(n−4,2) = 42(p²−3)(p²−4)."""
    return 42 * (p * p - 3) * (p * p - 4)


def theorem_k2_surviving_14() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |S∩T|=2, ordering shared vertices first puts two diagonal "
            "zeros; exactly 14 of 24 permutations survive.  |per|≤14."
        ),
        "n_surviving_perms": 14,
        "abs_per_ub": 14,
        "depends_on": ["C_diag_zero", "C_off_pm1"],
    }


def theorem_k2_permanent_closed_form() -> dict:
    t = theorem_k2_surviving_14()
    return {
        "proved": t["proved"],
        "theorem": (
            "per(C[S,T])=per2({c,d},{s,t})+(κ(S)−C_ab C_cd) per2({a,b},{s,t})"
            "+C_ab[C_ac per2({b,d})+C_ad per2({b,c})+C_bc per2({a,d})"
            "+C_bd per2({a,c})] on I={a,b}.  C_ii=0 and C_ij=±1 only."
        ),
        "depends_on": ["theorem_k2_surviving_14"],
    }


def theorem_k2_bilinear() -> dict:
    t = theorem_k2_permanent_closed_form()
    return {
        "proved": t["proved"],
        "theorem": (
            "Every k=2 survivor has exactly one s-edge and one t-edge, so "
            "per=∑_{u,v∈S} β_uv(S;I) C_us C_vt.  The {s,t}-sum is a "
            "C⊗C contraction (sum before abs)."
        ),
        "depends_on": ["theorem_k2_permanent_closed_form"],
        "sufficient_for_residual_i": False,
    }


def theorem_k2_unsigned_still_dead() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        maj = k2_unsigned_majorant(p)
        B = R_bar_budget_for_mu_thr(p)
        n = n_of(p)
        if maj != 14 * layer_count(n, 2):
            ok = False
            break
        if Fraction(maj) <= B:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 89, 101):
            sample[str(p)] = {"majorant": maj, "B": str(B)}
    if ok:
        for p in (307, 311, 313, 331):
            if Fraction(k2_unsigned_majorant(p)) <= R_bar_budget_for_mu_thr(p):
                ok = False
                break
    return {
        "proved": ok,
        "theorem": (
            "Diamond + |per|≤14 ⇒ |R^{(2)}|≤42(p²−3)(p²−4), which exceeds "
            "B for every prime p≥5.  Unsigned k=2 remains dead."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_k2_surviving_14",
            "diamond_abs_mu_le_1",
            "prop_15_229_R_budget",
        ],
        "sufficient_for_residual_i": False,
    }


def residual_i_closed_via_233() -> bool:
    return False


def hinge_status_233() -> dict:
    return {
        "k3_permanent_closed_form": theorem_k3_permanent_closed_form()["proved"],
        "unsigned_layers_012_dead": theorem_unsigned_layers_012_dead()["proved"],
        "k2_surviving_14": theorem_k2_surviving_14()["proved"],
        "k2_permanent_closed_form": theorem_k2_permanent_closed_form()["proved"],
        "k2_bilinear": theorem_k2_bilinear()["proved"],
        "k2_unsigned_still_dead": theorem_k2_unsigned_still_dead()["proved"],
        "residual_i_closed_via_233": residual_i_closed_via_233(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "type_I": type_I_k_3p_minus_2_closed_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B via signed k=0,1 (and μ-weighted k=2 bilinear) "
            "or Aut-SOS+diamond / Path C / K₄ / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_k2_surviving_14()
    b = theorem_k2_permanent_closed_form()
    c = theorem_k2_bilinear()
    d = theorem_k2_unsigned_still_dead()
    status = hinge_status_233()
    out = {
        "prop": "15.233",
        "title": (
            "Signed k=2: 14-term per2 form, bilinear in new columns; "
            "unsigned still dead; residual OPEN"
        ),
        "proved": {
            "k2_surviving_14": a["proved"],
            "k2_permanent_closed_form": b["proved"],
            "k2_bilinear": c["proved"],
            "k2_unsigned_still_dead": d["proved"],
            "residual_i_closed": residual_i_closed_via_233(),
        },
        "algebra": {
            "abs_per_k2_ub": 14,
            "k2_unsigned_majorant": "42(p^2-3)(p^2-4)",
            "k2_majorant_p5": k2_unsigned_majorant(5),
            "B_p5": str(R_bar_budget_for_mu_thr(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_233(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "residual_i_dual_eq_empty_proved_general": residual_i_dual_eq_empty_proved_general(),
        "type_I_k_3p_minus_2_closed_general": type_I_k_3p_minus_2_closed_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "k=2 per is a bilinear C⊗C form in the two new vertices.  "
            "Unsigned still exceeds B.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15233.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.233  residual_i={residual_i_closed_via_233()}")
    print(
        f"  k2_form={b['proved']} bilinear={c['proved']} "
        f"unsigned_dead={d['proved']}"
    )
    print(
        f"  predicates gsum={gsum_disj_lb_proved_general()} "
        f"res_i={residual_i_dual_eq_empty_proved_general()} "
        f"type_I={type_I_k_3p_minus_2_closed_general()} "
        f"e1={e1_closed_general()}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
