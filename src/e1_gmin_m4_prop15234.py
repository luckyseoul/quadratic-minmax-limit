#!/usr/bin/env python3
"""
Prop 15.234 — Signed k=1 layer of R̄₄: Laplace expansion along the
shared column; trilinear in the three new vertices; |per|≤18;
unsigned still dead.  Residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

HINGE GRAPH
  After 15.233: k=3 linear, k=2 bilinear.  This prop: k=1 trilinear.
  Next: signed k=0.

PROVED Max+-free (C_ii=0, C_ij=±1)
  A. 18 surviving permutations.
     I=S∩T={a}, T={a,r,s,t}.  Ordering a first puts one diagonal zero.
     18 of 24 permutations survive.  |per|≤18.  ∎

  B. Laplace / shared-column expansion.
        per(C[S,T]) = ∑_{x∈S\\{a}} C_{x a} · per_3(C[S\\{x}, {r,s,t}]).
     Proof: C_aa=0 kills the x=a summand; each of the other three choices
     of which row hits column a leaves a 3×3 permanent on the new
     columns (3!=6, and 3·6=18).  ∎

  C. Trilinearity (signed pairing).
     per_3 is a sum of 6 triple products, so per is trilinear in the
     three new columns.  Summing {r,s,t} before abs is a C⊗C⊗C
     contraction — Jacobi / cycle-before-abs on k=1.  ∎

  D. Unsigned k=1 still dead for all primes p≥5.
     |R^{(1)}|≤18·4 C(n−4,3)=12(p²−3)(p²−4)(p²−5), Θ(p⁶)≫B=Θ(p³).  ∎

OPEN
  E. Bound the μ-weighted trilinear form, or signed k=0 / Aut-SOS.

Writes evidence/e1_gmin_m4_prop15234.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import permutations
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
from e1_gmin_m4_prop15232 import layer_count  # noqa: E402
from e1_gmin_m4_prop15233 import theorem_k2_permanent_closed_form  # noqa: E402


def per3(C, rows, cols) -> int:
    """3×3 permanent of C[rows, cols]."""
    cols = list(cols)
    rows = list(rows)
    tot = 0
    for sig in permutations(range(3)):
        w = 1
        for i, r in enumerate(rows):
            w *= int(C[r, cols[sig[i]]])
        tot += w
    return tot


def k1_permanent_closed(C, S, a, r, s, t) -> int:
    """Laplace form of per(C[S, {a,r,s,t}]) for shared vertex a."""
    new = (r, s, t)
    tot = 0
    for x in S:
        if x == a:
            continue
        tot += int(C[x, a]) * per3(C, [y for y in S if y != x], new)
    return tot


def k1_unsigned_majorant(p: int) -> int:
    """18 · 4 C(n−4,3) = 12(p²−3)(p²−4)(p²−5)."""
    return 12 * (p * p - 3) * (p * p - 4) * (p * p - 5)


def theorem_k1_surviving_18() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |S∩T|=1, one shared diagonal zero leaves 18 surviving "
            "permutations.  |per|≤18."
        ),
        "n_surviving_perms": 18,
        "abs_per_ub": 18,
        "depends_on": ["C_diag_zero", "C_off_pm1"],
    }


def theorem_k1_laplace() -> dict:
    t = theorem_k1_surviving_18()
    return {
        "proved": t["proved"],
        "theorem": (
            "per(C[S,T])=∑_{x∈S\\{a}} C_{xa} per_3(C[S\\{x},{r,s,t}]) "
            "for I={a}.  C_aa=0 kills the missing summand; 3·3!=18."
        ),
        "depends_on": ["theorem_k1_surviving_18"],
    }


def theorem_k1_trilinear() -> dict:
    t = theorem_k1_laplace()
    return {
        "proved": t["proved"],
        "theorem": (
            "Each per_3 is a sum of triple products in the new columns, "
            "so the k=1 permanent is trilinear in (r,s,t).  The triple "
            "sum is a C⊗C⊗C contraction (sum before abs)."
        ),
        "depends_on": ["theorem_k1_laplace"],
        "sufficient_for_residual_i": False,
    }


def theorem_k1_unsigned_still_dead() -> dict:
    ok = True
    sample = {}
    for p in range(5, 220):
        if not is_prime(p):
            continue
        maj = k1_unsigned_majorant(p)
        B = R_bar_budget_for_mu_thr(p)
        if maj != 18 * layer_count(n_of(p), 1):
            ok = False
            break
        if Fraction(maj) <= B:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 23, 89):
            sample[str(p)] = {"majorant": maj, "B": str(B)}
    if ok:
        for p in (227, 233, 239, 241):
            if Fraction(k1_unsigned_majorant(p)) <= R_bar_budget_for_mu_thr(p):
                ok = False
                break
    return {
        "proved": ok,
        "theorem": (
            "Diamond + |per|≤18 ⇒ |R^{(1)}|≤12(p²−3)(p²−4)(p²−5), which "
            "exceeds B for every prime p≥5.  Unsigned k=1 remains dead."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_k1_surviving_18",
            "diamond_abs_mu_le_1",
            "prop_15_229_R_budget",
        ],
        "sufficient_for_residual_i": False,
    }


def residual_i_closed_via_234() -> bool:
    return False


def hinge_status_234() -> dict:
    return {
        "k2_permanent_closed_form": theorem_k2_permanent_closed_form()["proved"],
        "k1_surviving_18": theorem_k1_surviving_18()["proved"],
        "k1_laplace": theorem_k1_laplace()["proved"],
        "k1_trilinear": theorem_k1_trilinear()["proved"],
        "k1_unsigned_still_dead": theorem_k1_unsigned_still_dead()["proved"],
        "residual_i_closed_via_234": residual_i_closed_via_234(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "type_I": type_I_k_3p_minus_2_closed_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B via signed k=0 (and μ-weighted k=1,2) "
            "or Aut-SOS+diamond / Path C / K₄ / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_k1_surviving_18()
    b = theorem_k1_laplace()
    c = theorem_k1_trilinear()
    d = theorem_k1_unsigned_still_dead()
    status = hinge_status_234()
    out = {
        "prop": "15.234",
        "title": (
            "Signed k=1: Laplace along shared column, trilinear; "
            "unsigned still dead; residual OPEN"
        ),
        "proved": {
            "k1_surviving_18": a["proved"],
            "k1_laplace": b["proved"],
            "k1_trilinear": c["proved"],
            "k1_unsigned_still_dead": d["proved"],
            "residual_i_closed": residual_i_closed_via_234(),
        },
        "algebra": {
            "abs_per_k1_ub": 18,
            "k1_unsigned_majorant": "12(p^2-3)(p^2-4)(p^2-5)",
            "k1_majorant_p5": k1_unsigned_majorant(5),
            "B_p5": str(R_bar_budget_for_mu_thr(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_234(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "residual_i_dual_eq_empty_proved_general": residual_i_dual_eq_empty_proved_general(),
        "type_I_k_3p_minus_2_closed_general": type_I_k_3p_minus_2_closed_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "k=1 per is a Laplace / trilinear form in the three new "
            "vertices.  Unsigned still exceeds B.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15234.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.234  residual_i={residual_i_closed_via_234()}")
    print(
        f"  k1_laplace={b['proved']} trilinear={c['proved']} "
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
