#!/usr/bin/env python3
"""
Prop 15.239 — Per φ = (2n+1)φ and Per κ = (n−1)²κ − 6φ (any conference);
unconditional Per μ_part Cy-FE on Paley n=p²+1; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  Symmetric conference C: Cᵀ=C, C_ii=0, C_ij=±1, C²=(n−1)I.
  (Per f)_S := ∑_{|T|=4} per(C[S,T]) f(T)
            = ∑_{w₀,…,w₃ all distinct} (∏_{i=0}^{3} C_{s_i,w_i}) f({w's})
  (domain S ordered as s₀..s₃).  φ,κ the usual 4-set functions.
  Paley: n=p²+1 so (n−1)²=p⁴.

PROVED Max+-free (any symmetric conference of order n≥4)
  A. Per φ = (2n+1) φ on every 4-set S.
     Expand (Per φ)_S via ∑_x ∑_{inj w: avoid x} ∏_i (C_{s_i,w_i} C_{x,w_i})
     and ∏_{i<j}(1−δ_{w_i w_j}).  Pair-expansion partitions of {0,1,2,3}:
       • (1,1,1,1) and (1,1,2): vanish (singleton W-factors force x=s_i=s_j).
       • (4): Möbius sign −6; ∑_x W₄(x)=(n−1)φ  ⇒ contrib −6(n−1)φ.
         (W₄(x)=∑_{u≠x}∏_i C_{s_i,u}; ∑_x W₄=(n−1)∑_u∏C= (n−1)φ.)
       • (2,2): three partitions, each sign +1; for {{0,1},{2,3}},
         W₂({i,j},x)=−C_{s_i,x}C_{s_j,x}, product sum_x = ∑_x∏_i C_{s_i,x}=φ.
         Contrib 3φ.
       • (1,3): four partitions, each sign +2; W₁({j},x)=(n−1)1_{x=s_j},
         so sum_x W₃ W₁=(n−1)·(∑_{u≠s_j}∏_{i∈S} C_{s_i,u})=(n−1)φ.
         Contrib 8(n−1)φ.
     Total: [−6(n−1)+3+8(n−1)]φ=(2n+1)φ.  ∎
     (W₁,W₂ from C²=(n−1)I and C_ii=0 only.)

  B. Per κ = (n−1)² κ − 6 φ on every 4-set S.
     Write κ({w})=∑_{three matchings m} C_{w_a w_b}C_{w_c w_d}.
     (Per κ)_S = ∑_m ∑_{inj w} (∏_i C_{s_i,w_i}) C_{w_a w_b}C_{w_c w_d}.
     For each matching m=(ij|kl):
       Unrestricted sum = (n−1)² C_{s_i s_j} C_{s_k s_l}
       (factor as [∑_{u,v} C_{s_i u}C_{s_j v}C_{uv}]²-style:
        ∑_{u,v} C_{au}C_{bv}C_{uv}=(n−1)C_{ab} by C²).
     Sum over three matchings: (n−1)² κ(S).
     Non-injective remainder by image size:
       Size 1: C_{w_a w_b}=C_{uu}=0 ⇒ 0.
       Size 2 (need w_a≠w_b, w_c≠w_d else matching factor 0):
         For m=01|23: sum_{u≠v} A(u,v)B(u,v) with
           A=C_{s0 u}C_{s1 v}+C_{s0 v}C_{s1 u} (and B for 2,3)
         double-counts each w once as (u,v) and (v,u), so size2_m
         = ½∑_{u≠v}AB = −½∑_u A(u,u)B(u,u).
         A(u,u)=2C_{s0 u}C_{s1 u}, hence ∑_u AB_diag=4∑_u∏C_{s_i u}=4φ,
         size2_m=−2φ.  Three matchings: −6φ.
       Size 3 (doubleton one from {0,1} and one from {2,3}):
         Four doubleton slots (0,2),(0,3),(1,2),(1,3).  For (0,2):
           w0=w2=u, w1=a, w3=b distinct: weight sums to φ by C²
           (the (n−1)δ term vanishes; remainder −∑∏C = −φ, sign flip ⇒ φ).
         Four slots ⇒ size3_m=4φ.  Three matchings: 12φ.
     Non-inj total = −6φ+12φ=6φ.
     Hence Per κ = (n−1)²κ − 6φ.  ∎

  C. Per μ_part Cy-FE (Paley n=p²+1, primes p≥5).
     μ_part=aκ+bφ, a=(p²−1)/den, b=−2/den, den=p²(p²−5).
     By (A)+(B) with (n−1)²=p⁴:
        Per μ_part = a(p⁴κ−6φ)+b(2n+1)φ = p⁴ μ_part + 2φ
     on every 4-set (φ-coeff identity −6a+b(2n+1)=p⁴b+2 as in 15.238 A).
     On |κ|=1, Cy-FE is Per μ=p⁴μ+2φ, so Per δ=p⁴δ for μ=μ_part+δ.  ∎

CERTIFIED
  D. (A)(B) match full permanent operator on Paley p=3 (all 4-sets) and
     p=5 (all |κ| types sampled / full k1 types).

OPEN (blocks residual i / predicate flip)
  E. Bound |μ|≤1/(2p) (or 2/n / m4₂ / K₄ thr / Path C / free-e+ker=sc)
     for all primes p≥5.  The Per eigenforms close the particular Cy channel
     but do not by themselves control δ.  Soft-close forbidden.

Until E: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15239.json
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
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15238 import (  # noqa: E402
    mu_part_coeffs,
    p4_mupart_plus_2phi,
    per_kappa_form,
    per_mupart_via_eigen,
    per_phi_form,
    theorem_per_mupart_identity,
)


def theorem_Per_phi() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Any symmetric conference of order n≥4: (Per φ)_S=(2n+1)φ_S "
            "for every 4-set S.  Proof: inj-expansion + pair IE; types "
            "(1^4) and (2,1^2) vanish; (4)⇒−6(n−1)φ; (2,2)⇒3φ; "
            "(3,1)⇒8(n−1)φ; sum=(2n+1)φ.  Uses only C²=(n−1)I and C_ii=0."
        ),
        "depends_on": ["C2_eq_n_minus_1_I", "C_diag_zero"],
    }


def theorem_Per_kappa() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Any symmetric conference of order n≥4: (Per κ)_S=(n−1)²κ_S−6φ_S "
            "for every 4-set S.  Proof: κ=sum of 3 matchings; each unrestricted "
            "injection-free sum is (n−1)² times that matching of C[S,S]; "
            "total (n−1)²κ.  Size-1 collisions 0; size-2 per matching −2φ "
            "(A,B double-count + C²); size-3 per matching 4φ (four cross "
            "doubletons each φ); non-inj total 6φ.  Hence Per κ=(n−1)²κ−6φ."
        ),
        "depends_on": ["C2_eq_n_minus_1_I", "C_diag_zero", "C_symmetric"],
    }


def theorem_Per_mupart_unconditional() -> dict:
    """Paley: Per μ_part = p⁴ μ_part + 2φ on all 4-sets."""
    t_phi = theorem_Per_phi()
    t_kap = theorem_Per_kappa()
    t_id = theorem_per_mupart_identity()
    ok = t_phi["proved"] and t_kap["proved"] and t_id["proved"]
    # recheck identity with (n-1)^2 = p^4
    for p in range(5, 200):
        if not is_prime(p):
            continue
        n = n_of(p)
        if (n - 1) ** 2 != p**4:
            ok = False
            break
        for kappa, phi in ((1, -2), (1, 2 * (p - 2)), (-1, 2)):
            # per forms with (n-1)^2
            a, b = mu_part_coeffs(p)
            Pmp = a * ((n - 1) ** 2 * kappa - 6 * phi) + b * ((2 * n + 1) * phi)
            if Pmp != p4_mupart_plus_2phi(p, kappa, phi):
                ok = False
                break
        if not ok:
            break
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5 (Paley n=p²+1): Per μ_part = p⁴ μ_part + 2φ "
            "on every 4-set.  Combines (A)(B) with (n−1)²=p⁴ and the "
            "coefficient identity of 15.238 A."
        ),
        "depends_on": [
            "theorem_Per_phi",
            "theorem_Per_kappa",
            "theorem_per_mupart_identity",
        ],
    }


def theorem_delta_Per_on_kappa1() -> dict:
    t = theorem_Per_mupart_unconditional()
    return {
        "proved_structure": t["proved"],
        "bound_proved_general": False,
        "theorem": (
            "On |κ|=1 for Paley: Cy-FE + (C) ⇒ (Per δ)_S = p⁴ δ_S for "
            "μ=μ_part+δ.  Does not bound |δ| or |μ|."
        ),
        "depends_on": [
            "theorem_Per_mupart_unconditional",
            "prop_15_229_functional_eq",
        ],
    }


def residual_i_closed_via_239() -> bool:
    return False


def hinge_status_239() -> dict:
    return {
        "Per_phi": theorem_Per_phi()["proved"],
        "Per_kappa": theorem_Per_kappa()["proved"],
        "Per_mupart_unconditional": theorem_Per_mupart_unconditional()[
            "proved"
        ],
        "delta_Per_structure": theorem_delta_Per_on_kappa1()[
            "proved_structure"
        ],
        "residual_i_closed_via_239": residual_i_closed_via_239(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|μ|≤1/(2p) or 2/n or m4₂/K4/Path C/free-e+ker=sc for all p≥5"
        ),
    }


def main() -> dict:
    a = theorem_Per_phi()
    b = theorem_Per_kappa()
    c = theorem_Per_mupart_unconditional()
    d = theorem_delta_Per_on_kappa1()
    status = hinge_status_239()
    out = {
        "prop": "15.239",
        "title": (
            "Per φ=(2n+1)φ and Per κ=(n−1)²κ−6φ Max+-free; "
            "Per μ_part Cy-FE unconditional; residual OPEN"
        ),
        "proved": {
            "Per_phi": a["proved"],
            "Per_kappa": b["proved"],
            "Per_mupart_unconditional": c["proved"],
            "delta_Per_structure": d["proved_structure"],
            "residual_i_closed": residual_i_closed_via_239(),
        },
        "algebra": {
            "Per_phi": "(2*n+1)*phi",
            "Per_kappa": "(n-1)^2*kappa - 6*phi",
            "Paley": "(n-1)^2 = p^4",
            "Per_mupart": "p^4*mu_part + 2*phi",
            "thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_239(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Permanent eigenforms on span{κ,φ} are Max+-free for every "
            "symmetric conference.  Particular Cy channel closed on Paley. "
            "Residual-(i) still needs a pointwise |μ| bound.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15239.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.239  residual_i={residual_i_closed_via_239()}")
    print(
        f"  Per_phi={a['proved']} Per_kappa={b['proved']} "
        f"mupart={c['proved']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
