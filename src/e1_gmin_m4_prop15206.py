#!/usr/bin/env python3
"""
Prop 15.206 — Local |κ|-extension degrees n₃=(n−6)/4 on |κ|=1 centres
(Max+-free conference algebra); residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP
  C symmetric conference order n (Cᵀ=C, zero diag, off-diag ±1, C²=(n−1)I).
  Paley: n=p²+1.  For 4-set S and centre a∈S, rest T=S\\{a}:
    S_a→r := (S\\{a}) ∪ {r} (r∉S),
    κ(S) = ∑_{three matchings} C_e C_{e'}.
  n₃(a,S) := #{r∉S : |κ(S_a→r)|=3},
  n₁(a,S) := #{r∉S : |κ(S_a→r)|=1}.
  Always n₁+n₃ = n−4.
  One-center residuals (15.77): S₃ = ∑_{r:|κ'|=3} C_{ar} ρ(S_a→r),
  and on |κ|=1, star_a=+1: S₁+S₃ = pρ − 2/p² (⇒ m₄ bound via joint).

PROVED Max+-free (any symmetric conference, C²=(n−1)I)
  A. Pair sum.  For b≠c: ∑_r C_{rb} C_{rc} = 0.  Restricting off a triangle
     △={b,c,d}: ∑_{r∉△} C_{rb} C_{rc} = −C_{bd} C_{cd}.  ∎

  B. Sign-pattern count.  Fix triangle △={b,c,d}, set
        e_b := C_{cd}, e_c := C_{bd}, e_d := C_{bc},
        e := (e_b, e_c, e_d) ∈ {±1}³,
        N(ε) := #{r∉△ : (C_{rb},C_{rc},C_{rd})=ε}.
     Fourier on ({±1})³ (or expand ∏(1+ε_i C_{r v_i})/8):
        N(e)+N(−e)
          = (1/4)[ (n−3) + ∑_{pairs ij} e_i e_j s_{ij} ]
     with s_{ij}=∑_{r∉△} C_{r v_i} C_{r v_j} = −e_i e_j by (A).
     Each pair term e_i e_j s_{ij} = −1, three pairs ⇒ −3.
     Hence N(e)+N(−e) = (1/4)(n−3−3) = (n−6)/4.  ∎
     (Requires 4∣(n−6); for Paley n=p²+1, p odd: p²≡1 (mod 8) ⇒ n≡2
     (mod 8) ⇒ n−6≡4 (mod 8) ⇒ 4∣(n−6).)  ∎

  C. Link to κ.  For S={a}∪△ and r∉△:
        κ(S_a→r) = C_{rb} e_b + C_{rc} e_c + C_{rd} e_d,
     so |κ(S_a→r)|=3 ⇔ (C_{rb},C_{rc},C_{rd}) ∈ {e,−e}.  ∎

  D. Subtract the centre.  Vertex a lies in the complement of △.  Moreover
        (C_{ab},C_{ac},C_{ad})=±e  ⇔  |κ(S)|=3.
     Therefore, counting r∉S = r∉△ and r≠a:
        n₃(a,S) = N(e)+N(−e) − 1_{|κ(S)|=3}
                 = (n−6)/4           if |κ(S)|=1,
                 = (n−6)/4 − 1       if |κ(S)|=3.  ∎

  E. Corollary on |κ|=1.  For every |κ|=1 four-set and every centre,
        n₃ = (n−6)/4,   n₁ = n−4−n₃ = (3n−10)/4.
     For Paley n=p²+1: n₃=(p²−5)/4, n₁=(3p²−7)/4.  ∎

  F. Consequence for residual-(i) one-center attack.  On |κ|=1,
        |S₃| ≤ ∑_{|κ'|=3} |ρ'| ≤ n₃ · max_{|κ|=3}|ρ|.
     Combined with star·S₁≤0 and S₁+S₃=pρ−2/p² (15.77–78), a
     Max+-free bound max_{|κ|=3}|ρ| ≤ B(p) with
        n₃ B(p) ≤ p ρ_cand − 2/p² = p M_cand − 1/p − 2/p²
     would force ρ≤ρ_cand ⇒ m₄≤M_cand ⇒ |μ₄|≤1/(2p) (15.205)
     ⇒ Gsum≥−1/p ⇒ residual-(i) Farkas.  At p=5 the joint budget is
     negative, so one also needs strongly negative S₁ on maximisers
     (certified sharp at p=5).  Crude B=O(1) is far too weak (cancellation
     in signed S₃ is essential).  ∎

CERTIFIED
  G. n₃∈{(n−6)/4} on all |κ|=1 centres and n₃=(n−6)/4−1 on |κ|=3
     centres, full enumeration of centres for p=3,5 and large samples p=7.

OPEN (blocks residual i / predicate flip)
  H. Prove |μ₄|≤M_cand (e.g. via star·S₁≤0 + tight signed S₃ bound using n₃),
     or free-e≤r(p)·scheme, or K₄≤Wick_hi / ker=sc+cost_D.
     The degree formula (E) is Max+-free but does not alone close residual (i).

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15206.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15205 import M_cand  # noqa: E402


def n3_formula_kappa1(p: int) -> int:
    """(n−6)/4 for Paley n=p²+1."""
    n = n_of(p)
    assert (n - 6) % 4 == 0
    return (n - 6) // 4


def n1_formula_kappa1(p: int) -> int:
    n = n_of(p)
    return n - 4 - n3_formula_kappa1(p)


def theorem_n3_degree() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Any symmetric conference C (C²=(n−1)I, zero diag, ±1 off): "
            "for 4-set S and centre a, n₃(a,S)=(n−6)/4 if |κ(S)|=1, "
            "and (n−6)/4−1 if |κ(S)|=3.  Proof: pair sums from C² give "
            "N(e)+N(−e)=(n−6)/4 for opposite-triangle sign patterns e; "
            "centre a matches ±e iff |κ(S)|=3.  Paley: 4∣(n−6)."
        ),
        "corollaries": {
            "kappa1_n3": "(n−6)/4",
            "kappa1_n1": "(3n−10)/4",
            "paley_n3": "(p²−5)/4",
            "paley_n1": "(3p²−7)/4",
        },
    }


def theorem_paley_divisibility() -> dict:
    """4∣(n−6) for n=p²+1, p odd."""
    ok = True
    sample = {}
    for p in range(3, 60, 2):
        if not is_prime(p):
            continue
        n = n_of(p)
        div = (n - 6) % 4 == 0
        # p odd ⇒ p²≡1 (mod 8) ⇒ n≡2 (mod 8) ⇒ n−6≡4 (mod 8) ⇒ 4|n−6
        if p in (3, 5, 7, 11):
            sample[str(p)] = {"n": n, "n_minus_6_mod4": (n - 6) % 4, "ok": div}
        if not div:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For odd p: p²≡1 (mod 8), n=p²+1≡2 (mod 8), n−6≡4 (mod 8), "
            "hence 4∣(n−6) so (n−6)/4 is an integer."
        ),
        "by_p_sample": sample,
    }


def certify_n3_counts(primes: list[int] | None = None, max_sets: int = 3000) -> dict:
    """Census n₃ on |κ|=1 and |κ|=3 centres."""
    from minmax_quadratic import paley_conference_prime_power

    if primes is None:
        primes = [3, 5, 7]
    rows = {}
    all_ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(int)
        n = C.shape[0]
        base = (n - 6) // 4
        ok1 = ok3 = 0
        bad = 0
        n_s = 0
        for S in combinations(range(n), 4):
            a, b, c, d = S
            kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
            ak = abs(kap)
            for centre in (a,):  # one centre enough per set for speed; full for small p
                if p <= 5:
                    centres = S
                else:
                    centres = (a,)
                for cen in centres:
                    rest = [x for x in S if x != cen]
                    bb, cc, dd = rest
                    eb, ec, ed = C[cc, dd], C[bb, dd], C[bb, cc]
                    n3 = 0
                    for r in range(n):
                        if r in S:
                            continue
                        kr = C[r, bb] * eb + C[r, cc] * ec + C[r, dd] * ed
                        if abs(kr) == 3:
                            n3 += 1
                    expect = base if ak == 1 else base - 1
                    if n3 == expect:
                        if ak == 1:
                            ok1 += 1
                        else:
                            ok3 += 1
                    else:
                        bad += 1
            n_s += 1
            if p >= 7 and n_s >= max_sets:
                break
        rows[str(p)] = {
            "base": base,
            "formula_n3_k1": n3_formula_kappa1(p),
            "ok_kappa1_centres": ok1,
            "ok_kappa3_centres": ok3,
            "bad": bad,
            "match_formula": bad == 0,
        }
        if bad:
            all_ok = False
    return {
        "certified": all_ok,
        "proved_general": True,  # algebra is the proof; census is check
        "by_p": rows,
        "note": "Algebra proves general; census is integrity check.",
    }


def joint_budget_cand(p: int) -> Fraction:
    """p·ρ_cand − 2/p² = p M_cand − 1/p − 2/p² (need S₁+S₃ ≤ this if using joint)."""
    mc = M_cand(p)
    return p * mc - Fraction(1, p) - Fraction(2, p * p)


def theorem_crude_S3_insufficient() -> dict:
    """Document that n₃·O(1) does not beat joint budget (need signed cancel)."""
    rows = {}
    for p in (5, 7, 11, 13):
        n3 = n3_formula_kappa1(p)
        jb = joint_budget_cand(p)
        # crude |m4|≤(p-2)/p on |κ|=3 ⇒ |ρ|≤(p-2)/p+3/p²
        rho3 = Fraction(p - 2, p) + Fraction(3, p * p)
        s3ub = n3 * rho3
        rows[str(p)] = {
            "n3": n3,
            "joint_budget": str(jb),
            "crude_S3_ub": str(s3ub),
            "crude_beats_budget": s3ub <= jb,
        }
    return {
        "proved": True,
        "theorem": (
            "n₃=(n−6)/4 alone with |ρ|≤(p−2)/p+3/p² does not imply "
            "S₃≤joint_budget (crude_S3_ub ≫ budget).  Signed cancellation "
            "in S₃=∑ C_ar ρ' is required for the M_cand path."
        ),
        "by_p": rows,
    }


def free_e_bound_proved_general() -> bool:
    return False


def mu_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_206() -> bool:
    return free_e_bound_proved_general() or mu_bound_proved_general()


def hinge_status_206() -> dict:
    return {
        "n3_degree_proved": True,
        "crude_S3_insufficient": True,
        "Mcand_via_S3_general": False,
        "residual_i_closed": residual_i_closed_via_206(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Signed S₃ bound or star·S₁≤0+joint for |μ|≤M_cand; or free-e/K₄/ker."
        ),
    }


def main() -> dict:
    th = theorem_n3_degree()
    div = theorem_paley_divisibility()
    cert = certify_n3_counts()
    crude = theorem_crude_S3_insufficient()
    out = {
        "title": (
            "Prop 15.206 n₃=(n−6)/4 on |κ|=1 centres Max+-free; "
            "residual (i) OPEN"
        ),
        "proved": {
            "n3_degree": th["proved"],
            "paley_divisibility": div["proved"],
            "crude_S3_insufficient": crude["proved"],
            "Mcand_mu_bound_general": False,
            "residual_i_closed": residual_i_closed_via_206(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "n3_degree": th,
            "divisibility": div,
            "crude_S3": crude,
        },
        "certified": {"n3_census": cert},
        "hinge": hinge_status_206(),
        "F3": (
            "n₃ formula is Max+-free and general, but residual (i) still needs "
            "|μ|≤M_cand / free-e / K₄ / ker.  Do not flip predicates."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15206.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.206 n3=(n-6)/4 on |kappa|=1 centres")
    print(f"  n3 degree proved: {th['proved']}")
    print(f"  4|(n-6) Paley: {div['proved']}")
    print(f"  census ok: {cert['certified']}")
    for pk, v in cert["by_p"].items():
        print(
            f"  p={pk}: base={v['base']} ok_k1={v['ok_kappa1_centres']} "
            f"ok_k3={v['ok_kappa3_centres']} bad={v['bad']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_206()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
