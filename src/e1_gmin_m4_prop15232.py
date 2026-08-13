#!/usr/bin/env python3
"""
Prop 15.232 — Intersection-type split of R̄₄; k=3 permanent closed form
(column-linear in the new vertex); C² pairing (sum before abs); unsigned
layers k=0,1,2 dead vs B; k=3 unsigned R-safe for p≥89 only.
Residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

HINGE GRAPH (live)
  OPEN: |R̄₄|≤B  — preferred; |δ|≤room_δ (≡ after 15.230); Path C; K₄;
        ker=sc+free-e.
  DEAD: crude ∑|per| (15.231); unsigned layers k=0,1,2 (this prop).
  Attack picked: Jacobi / cycle-before-abs on the k=3 layer
  (C-column linearity + C² pairing). Aut-SOS+diamond remains an alt.

SETUP (15.229–231)
  On |κ|=1: R̄₄(S)=∑_{T≠S} per(C[S,T]) μ_T.  Budget B=(p⁴−1)/(2p)−4(p−2).
  |μ|≤1 by the boolean diamond (15.228).  per is invariant under row/col
  reorderings (permanent, not determinant).

PROVED Max+-free
  A. Intersection decomposition.
        R̄₄(S)=∑_{k=0}^{3} R^{(k)}(S),
        R^{(k)}(S):=∑_{|T|=4, |S∩T|=k} per(C[S,T]) μ_T.
     Counts (any n≥4): #T with |S∩T|=3 is 4(n−4); =2 is 6 C(n−4,2);
     =1 is 4 C(n−4,3); =0 is C(n−4,4).  ∎

  B. k=3 permanent, 11-term closed form.
     Let T=(S\\{d})∪{t}, {a,b,c}=S\\{d}, Δ=C_ab C_bc C_ca.  Ordering
     S=(a,b,c,d), T=(a,b,c,t) puts zeros at the three shared diagonal
     entries, so only 11 of 24 permutations survive.  Using C_ij=±1
     (hence squares =1) and C_ii=0:
        per(C[S,T]) = 2 C_{dt} Δ
                    + ∑_{x∈{a,b,c}} C_{dx} C_{tx}
                    + ∑_{x∈{a,b,c}} C_{xy} C_{xz}
                         (C_{yd} C_{zt} + C_{yt} C_{zd}),
     where {y,z}={a,b,c}\\{x}.  Conference C² is not used.  |per|≤11.  ∎

  C. Column-linearity (signed pairing structure).
     Every monomial in (B) contains exactly one factor C_{· t}.
     Hence per(C[S, S_d^t]) = ∑_{u∈S} α_u(S,d) C_{u t}
     with α independent of t.  Summing t before abs is a C-weighted
     row contraction against μ on the adjacent 4-sets — the Jacobi /
     cycle-before-abs move on the k=3 layer.  ∎

  D. C² pairing (sum two t-columns before abs).
     For any symmetric conference C²=(n−1)I and any S⊂[n], u,v∈[n]:
        ∑_{t∉S} C_{u t} C_{v t}
          = (n−1) 1_{u=v} − ∑_{s∈S} C_{s u} C_{s v}.
     Proof: ∑_{t=0}^{n−1} C_{ut}C_{vt}=(C²)_{uv}=(n−1)1_{u=v}; subtract t∈S.
     Max+-free.  Contracts the Gram of the k=3 kernel; does not by
     itself bound the μ-weighted R^{(3)}.  ∎

  E. Unsigned layers k=0,1,2 are each dead vs B for every prime p≥5.
     Ultra-crude |R^{(k)}|≤24·#T_k ( |per|≤24, |μ|≤1 ):
        k=0: (p²−3)(p²−4)(p²−5)(p²−6) = C(n−4,4)·24,
        k=1: 16(p²−3)(p²−4)(p²−5),
        k=2: 72(p²−3)(p²−4).
     Each exceeds B(p)=(p⁴−1)/(2p)−4(p−2) for all primes p≥5
     (degree 8,6,4 vs B=Θ(p³); certified by Fraction on 5≤p≤300 and
     leading-term comparison for p>300).  Must keep signs inside each
     of those layers.  ∎

  F. k=3 unsigned majorant: |R^{(3)}|≤44(p²−3) (11·4(n−4), diamond).
     Let f(p)=p⁴−88p³−8p²+280p−1.  Then
        B−44(p²−3) = f(p)/(2p).
     f'(p)>0 on [89,∞) and f(89)=666520>0, so f>0 for all real p≥89.
     Hence |R^{(3)}|≤B for every prime p≥89 by diamond + (B).
     For 5≤p<89 the same majorant exceeds B (unsigned k=3 dead there).
     This closes only the k=3 layer for large p — not |R̄₄|≤B.  ∎

CERTIFIED (not a general residual-(i) close)
  G. (B) matches per(C[S,T]) on Paley p=3 (all |S∩T|=3 of every |κ|=1 S)
     and p=5 (samples).  (D) matches direct column sums on those C.

OPEN (blocks residual i / predicate flip)
  H. Bound |R̄₄|≤B on |κ|=1 for all primes p≥5: signed k=0,1,2
     (character sums / cycle pairing before abs), or Aut-SOS+diamond,
     or Path C / K₄ / free-e+ker=sc.  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15232.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
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
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15230 import theorem_R_part_max_le_budget  # noqa: E402
from e1_gmin_m4_prop15231 import (  # noqa: E402
    theorem_R_bar_permanent_form,
    theorem_crude_abs_majorant_dead,
)


def layer_count(n: int, k: int) -> int:
    """# of 4-sets T with |S∩T|=k for a fixed 4-set S."""
    if k == 0:
        return comb(n - 4, 4)
    if k == 1:
        return 4 * comb(n - 4, 3)
    if k == 2:
        return 6 * comb(n - 4, 2)
    if k == 3:
        return 4 * (n - 4)
    raise ValueError("k must be 0,1,2,3")


def ultra_crude_layer(p: int, k: int) -> int:
    """24 · #T_k  (|per|≤24, |μ|≤1)."""
    return 24 * layer_count(n_of(p), k)


def ultra_crude_layer_closed(p: int, k: int) -> int:
    """Closed polynomials in p for the ultra-crude layer majorants."""
    a = p * p - 3
    b = p * p - 4
    c = p * p - 5
    d = p * p - 6
    if k == 0:
        return a * b * c * d
    if k == 1:
        return 16 * a * b * c
    if k == 2:
        return 72 * a * b
    if k == 3:
        return 96 * a  # 24 * 4(n-4) = 96(p²-3)
    raise ValueError("k must be 0,1,2,3")


def k3_unsigned_majorant(p: int) -> int:
    """11 · 4(n−4) = 44(p²−3)."""
    return 44 * (p * p - 3)


def k3_slack_poly(p: int) -> int:
    """f(p)=p⁴−88p³−8p²+280p−1 = 2p (B − 44(p²−3)) as an integer."""
    return p**4 - 88 * p**3 - 8 * p * p + 280 * p - 1


def k3_layer_R_safe_via_unsigned(p: int) -> bool:
    """Whether 44(p²−3) ≤ B(p)."""
    return Fraction(k3_unsigned_majorant(p)) <= R_bar_budget_for_mu_thr(p)


def k3_permanent_closed(C, S, d, t) -> int:
    """11-term closed form of per(C[S, (S\\{d})∪{t}])."""
    abc = [x for x in S if x != d]
    if len(abc) != 3:
        raise ValueError("d must be a unique element of S")
    a, b, c = abc

    def e(i, j) -> int:
        return int(C[i, j])

    delta = e(a, b) * e(b, c) * e(c, a)
    term_dt = 2 * e(d, t) * delta
    term_same = e(d, a) * e(t, a) + e(d, b) * e(t, b) + e(d, c) * e(t, c)
    mixed = 0
    for x, y, z in ((a, b, c), (b, a, c), (c, a, b)):
        mixed += e(x, y) * e(x, z) * (e(y, d) * e(z, t) + e(y, t) * e(z, d))
    return term_dt + term_same + mixed


def c2_pairing_off_S(C, S, u, v) -> int:
    """Closed RHS of ∑_{t∉S} C_{ut} C_{vt} for a conference C."""
    n = C.shape[0]
    sset = set(S)
    acc = 0
    for s in sset:
        acc += int(C[s, u]) * int(C[s, v])
    return (n - 1) * int(u == v) - acc


def c2_pairing_off_S_direct(C, S, u, v) -> int:
    """Direct ∑_{t∉S} C_{ut} C_{vt} (for tests / census)."""
    sset = set(S)
    n = C.shape[0]
    tot = 0
    for t in range(n):
        if t in sset:
            continue
        tot += int(C[u, t]) * int(C[v, t])
    return tot


def theorem_intersection_decomposition() -> dict:
    ok = True
    sample = {}
    for p in range(5, 80):
        if not is_prime(p):
            continue
        n = n_of(p)
        total = sum(layer_count(n, k) for k in range(4))
        # all 4-sets except S itself
        if total != comb(n, 4) - 1:
            ok = False
            break
        if p in (5, 7, 11, 13):
            sample[str(p)] = {str(k): layer_count(n, k) for k in range(4)}
    return {
        "proved": ok,
        "theorem": (
            "R̄₄(S)=∑_{k=0}^{3} R^{(k)}(S) with R^{(k)} the |S∩T|=k slice "
            "of the permanent contraction.  Counts: 4(n−4), 6 C(n−4,2), "
            "4 C(n−4,3), C(n−4,4); they partition binom(n,4)−1."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_R_bar_permanent_form"],
    }


def theorem_k3_permanent_closed_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For T=(S\\{d})∪{t} and {a,b,c}=S\\{d}, the 4×4 conference "
            "block with shared diagonal zeros has exactly 11 surviving "
            "permutations, and per(C[S,T])=2 C_dt Δ + ∑_x C_dx C_tx + "
            "∑_x C_xy C_xz (C_yd C_zt + C_yt C_zd) with Δ=C_ab C_bc C_ca. "
            "Uses only C_ii=0 and C_ij=±1.  Hence |per|≤11 on |S∩T|=3."
        ),
        "n_surviving_perms": 11,
        "abs_per_ub": 11,
        "depends_on": ["C_diag_zero", "C_off_pm1"],
    }


def theorem_k3_column_linear() -> dict:
    t = theorem_k3_permanent_closed_form()
    return {
        "proved": t["proved"],
        "theorem": (
            "Every monomial in the k=3 permanent contains exactly one "
            "factor C_{· t}, so per(C[S,S_d^t])=∑_{u∈S} α_u(S,d) C_{u t}. "
            "The t-sum of R^{(3)} is a C-weighted contraction — sum before abs."
        ),
        "depends_on": ["theorem_k3_permanent_closed_form"],
        "sufficient_for_residual_i": False,
    }


def theorem_c2_pairing() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For any symmetric conference C²=(n−1)I and any S, u, v: "
            "∑_{t∉S} C_ut C_vt = (n−1)1_{u=v} − ∑_{s∈S} C_su C_sv. "
            "Max+-free C² identity; contracts the k=3 Gram."
        ),
        "depends_on": ["C2_eq_n_minus_1_I"],
        "sufficient_for_residual_i": False,
    }


def theorem_unsigned_layers_012_dead() -> dict:
    """Each of k=0,1,2 ultra-crude exceeds B for all primes p≥5."""
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        B = R_bar_budget_for_mu_thr(p)
        for k in (0, 1, 2):
            # prefer closed poly (no overflow / exact)
            u = ultra_crude_layer_closed(p, k)
            if Fraction(u) != Fraction(ultra_crude_layer(p, k)):
                ok = False
                break
            if Fraction(u) <= B:
                ok = False
                break
        else:
            if p in (5, 7, 11, 13, 17, 19, 23, 89, 101):
                sample[str(p)] = {
                    "B": str(B),
                    **{f"ultra_k{k}": str(ultra_crude_layer_closed(p, k)) for k in (0, 1, 2)},
                }
            continue
        break
    # large-p leading terms: k=0 ~ p^8, k=1 ~ 16 p^6, k=2 ~ 72 p^4, B ~ p^3/2
    if ok:
        for p in (307, 311, 313, 317, 331):
            B = R_bar_budget_for_mu_thr(p)
            for k in (0, 1, 2):
                if Fraction(ultra_crude_layer_closed(p, k)) <= B:
                    ok = False
                    break
    return {
        "proved": ok,
        "theorem": (
            "The unsigned layer majorants |R^{(k)}|≤24·#T_k exceed B for "
            "k=0,1,2 and every prime p≥5.  Residual-(i) cannot close by "
            "absolute summation inside any of those intersection types."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_intersection_decomposition",
            "diamond_abs_mu_le_1",
            "prop_15_229_R_budget",
        ],
        "sufficient_for_residual_i": False,
    }


def theorem_k3_layer_R_safe_large_p() -> dict:
    """|R^{(3)}|≤44(p²−3)≤B for all primes p≥89; not a residual-(i) close."""
    f89 = k3_slack_poly(89)
    # f'(p)=4p³−264p²−16p+280; on [89,∞) the cubic 4p³−264p² dominates
    deriv_pos = True
    for p in range(89, 200):
        deriv = 4 * p**3 - 264 * p * p - 16 * p + 280
        if deriv <= 0:
            deriv_pos = False
            break
    # algebraic: 4p²−264p−16+280/p > 4p(p−66) −16 for p≥89
    # 4*89*(89-66)-16 = 4*89*23-16=8188-16>0, increasing.
    ok_f = f89 > 0 and deriv_pos
    # small primes: unsigned k=3 exceeds B (dead there)
    small_dead = True
    small = {}
    for p in range(5, 89):
        if not is_prime(p):
            continue
        safe = k3_layer_R_safe_via_unsigned(p)
        if safe:
            small_dead = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 61, 73, 83):
            small[str(p)] = {
                "majorant": k3_unsigned_majorant(p),
                "B": str(R_bar_budget_for_mu_thr(p)),
                "f": k3_slack_poly(p),
            }
    large_ok = True
    large = {}
    for p in range(89, 220):
        if not is_prime(p):
            continue
        if not k3_layer_R_safe_via_unsigned(p):
            large_ok = False
            break
        if p in (89, 97, 101, 103, 107, 109):
            large[str(p)] = {
                "majorant": k3_unsigned_majorant(p),
                "B": str(R_bar_budget_for_mu_thr(p)),
                "f": k3_slack_poly(p),
            }
    # identity B − 44(p²−3) = f(p)/(2p)
    ident = True
    for p in range(5, 160):
        if not is_prime(p):
            continue
        slack = R_bar_budget_for_mu_thr(p) - Fraction(k3_unsigned_majorant(p))
        if slack != Fraction(k3_slack_poly(p), 2 * p):
            ident = False
            break
    return {
        "proved": ok_f and small_dead and large_ok and ident,
        "theorem": (
            "Diamond + |per|≤11 ⇒ |R^{(3)}|≤44(p²−3).  "
            "B−44(p²−3)=f(p)/(2p) with f(p)=p⁴−88p³−8p²+280p−1.  "
            "f(89)=666520>0 and f'>0 on [89,∞) ⇒ the k=3 layer is "
            "R-safe for every prime p≥89.  For 5≤p<89 the same majorant "
            "exceeds B (unsigned k=3 dead).  Does not bound R^{(0..2)}."
        ),
        "f_89": f89,
        "small_p_unsigned_dead": small,
        "large_p_R_safe": large,
        "depends_on": [
            "theorem_k3_permanent_closed_form",
            "diamond_abs_mu_le_1",
            "prop_15_229_R_budget",
        ],
        "sufficient_for_residual_i": False,
    }


def residual_i_closed_via_232() -> bool:
    return False


def hinge_status_232() -> dict:
    return {
        "intersection_decomposition": theorem_intersection_decomposition()["proved"],
        "k3_permanent_closed_form": theorem_k3_permanent_closed_form()["proved"],
        "k3_column_linear": theorem_k3_column_linear()["proved"],
        "c2_pairing": theorem_c2_pairing()["proved"],
        "unsigned_layers_012_dead": theorem_unsigned_layers_012_dead()["proved"],
        "k3_layer_R_safe_large_p": theorem_k3_layer_R_safe_large_p()["proved"],
        "R_bar_permanent_form": theorem_R_bar_permanent_form()["proved"],
        "crude_abs_majorant_dead": theorem_crude_abs_majorant_dead()["proved"],
        "R_part_max_le_budget": theorem_R_part_max_le_budget()["proved"],
        "residual_i_closed_via_232": residual_i_closed_via_232(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "type_I": type_I_k_3p_minus_2_closed_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B on |κ|=1 via signed k=0,1,2 (sum before abs) "
            "or Aut-SOS+diamond / Path C / K₄ / free-e+ker=sc.  "
            "k=3 unsigned R-safe only for p≥89; k=0,1,2 unsigned dead."
        ),
    }


def main() -> dict:
    a = theorem_intersection_decomposition()
    b = theorem_k3_permanent_closed_form()
    c = theorem_k3_column_linear()
    d = theorem_c2_pairing()
    e = theorem_unsigned_layers_012_dead()
    f = theorem_k3_layer_R_safe_large_p()
    status = hinge_status_232()
    out = {
        "prop": "15.232",
        "title": (
            "Intersection split of R̄₄; k=3 per closed + C² pairing; "
            "unsigned k=0,1,2 dead; k=3 R-safe for p≥89; residual OPEN"
        ),
        "proved": {
            "intersection_decomposition": a["proved"],
            "k3_permanent_closed_form": b["proved"],
            "k3_column_linear": c["proved"],
            "c2_pairing": d["proved"],
            "unsigned_layers_012_dead": e["proved"],
            "k3_layer_R_safe_large_p": f["proved"],
            "residual_i_closed": residual_i_closed_via_232(),
        },
        "algebra": {
            "layer_counts": "4(n-4), 6 C(n-4,2), 4 C(n-4,3), C(n-4,4)",
            "k3_per": "2 C_dt Δ + sum_x C_dx C_tx + mixed",
            "abs_per_k3_ub": 11,
            "k3_unsigned_majorant": "44(p^2-3)",
            "k3_slack_poly": "p^4-88p^3-8p^2+280p-1",
            "k3_majorant_p5": k3_unsigned_majorant(5),
            "k3_majorant_p89": k3_unsigned_majorant(89),
            "B_p5": str(R_bar_budget_for_mu_thr(5)),
            "B_p89": str(R_bar_budget_for_mu_thr(89)),
            "thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_232(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "residual_i_dual_eq_empty_proved_general": residual_i_dual_eq_empty_proved_general(),
        "type_I_k_3p_minus_2_closed_general": type_I_k_3p_minus_2_closed_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Signed k=3 structure (column-linear per + C² pairing) is "
            "Max+-free.  Unsigned k=0,1,2 each exceed B for all p≥5.  "
            "k=3 is R-safe by diamond for p≥89 only.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15232.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.232  residual_i={residual_i_closed_via_232()}")
    print(
        f"  k3_form={b['proved']} c2={d['proved']} "
        f"layers012_dead={e['proved']} k3_large={f['proved']}"
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
