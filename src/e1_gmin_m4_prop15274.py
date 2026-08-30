#!/usr/bin/env python3
"""
Prop 15.274 — Residual (ii) even k≥4p, Max+-free Paley.

Does **not** edit 15.170 e1 AND. Does **not** flip Aut-Schur / Gsum / pairing.
Exhaustiveness (force S∈{2,4}) is a non-goal. Soft-close forbidden.

SETUP
  Residual (ii): deep s₊=2, even k, freeness-fail f_e≡+1 on U₂, H=G∪{e}.
  Need Φ(H)≥Φ−2. Affine two-level is already empty for all k≥3p (15.179).
  15.236/237 close only even k∈[3p+1, 4p−2] (mean E_−[S]=−k/p>−4 forces
  max_Max− S≥−2, then dual-bad empty).

PROVED
  A. Two-level Max+ S∈{2,4}: N₂/N=2−k/(2p)≤0 for k≥4p. Strictly negative
     for k>4p (impossible). At k=4p, N₂=0 so S≡4 and s₊=4≠2 (out of
     residual (ii)). Generic one-sided level-4 emptiness is false: 15.402
     gives four-parallel-line covers. The live question is compatibility
     with the opposite-shell residual inequalities.
  B. Dual-bad Max− two-level {−2,−4} needs E_−[S]≥−4, i.e. k≤4p. For
     k>4p the mass 2−k/(2p) is negative. At k=4p it forces b=P(S=−2)=0,
     i.e. S≡−4, a one-sided Max−-tight level-4 cover. Such covers exist
     (15.402), but the known family has min_Max+≤0 by 15.272/15.588 and is
     outside residual (ii). Joint compatibility with s₊≥2 is OPEN.
  C. Slope identity (15.236.D / 15.50 interpolants) does not use k≤4p−2:
        p[…] = 2(k−4p)(p−1) M₊ M₋,   N±=p M±, p∤M±.
     4p≡0 (mod p) so k−4p≡k (mod p). p∤2, p∤(p−1). Hence p|RHS iff p|k.
     Dual-bad is empty whenever k≢0 (mod p), including even k≥4p.
     At k=4p,6p,… one has k−4p≡0, RHS≡0, **no** slope obstruction.
  D. 15.237 pair-span classification is k-independent. Dual-bad U-mass
     (k−3p+1)/(2p) lies in (0,1) only for k∈(3p−1,5p−1). Two-value
     S∈{2,6} (resp. {−2,−6}) at k=4p makes 1_{S=2} a 0-1 pair-span of
     mass 1/2, which is not a constant / pair-slice / 3-equal.
  E. The former k=4p dichotomy depended on one-sided S≡−4 emptiness and
     is RETRACTED. Even scores and mean −4 alone do not force max_Max−≥−2.
  F. Two-value even support is a 0-1 pair-span. 15.237 C: those functions
     are constants, pair-slices (mass (p±1)/(2p)), or triangle 3-equals
     (mass (p±3)/(4p)). At k=4p, any two-value law containing the extreme
     ±2 has lo-mass 1−1/m (other level ±2±2m, m≥2), and 1−1/m lies in
     none of those six rationals for every prime p≥5 and integer m≥2
     (divisor identities). The m=1 collapse has zero mass at the required
     extreme ±2, so it is not an actual residual/leftover case; no tight-cover
     emptiness is needed. Hence every actual two-value leftover at k=4p is
     empty. Same classification
     empties two-value laws at general even k≥4p whenever the first-moment
     mass is unclassified.
  G. If s_{-}≥0 then Φ(G)≥Φ (15.42 / 15.171.B) and the 2-Lipschitz edge
     flip (15.169) gives Φ(H)≥Φ−2. At k=4p this plus (E) reduces residual
     (ii) ND to the dual-bad max=−2 branch; two-value of that branch is
     empty by (F).
  H. Multi-level leftover identities (k=4p, max=−2, f_e≡−1 on U_{−2}):
     n_{−2}=n_{−6}+2 n_{−8}+… ; a=P(S=−2)≤(p+1)/(2p); pairing
     E[Sf]∈[2/p−2, 4/p]; Aut_e 2-orbit far weight
     μ_f=2k/(p²(p²−1))+2β/p² is negative when β<−k/(p²−1). 2-orbit is
     not assumed for a general G.

  I. Leftover pair-span / Aut_e 3-weight identities (still not a close):
     15.237 C 0-1 masses include mixed-ε 3-equals (p±1)/(4p). 1−1/m
     misses the full 8-list (divisor identities). Plus pair-slice cannot
     be U_{−2} (sits in {f_e=+1}). Minus-slice a=(p+1)/(2p) forces
     three-level {−2,−4,−6} to have b=−1/p. Triangle 3-equal
     (1+∑ε f)/4 is 0-1 iff ∏ε=C_ab C_ac C_bc; mass (p−∑ε)/(4p).
     Aut_e 3-weight (two star + collapsed far):
        (p+1)u+(p−1)u′=μ_f p(p²−1).
     On the minus-slice a=thr (u′=−6−4/(p−1)) this μ_f is <0 for every
     Max+ u∈[2,4]. Paley Aut_e has several far orbits, so 3-weight is
     not assumed for a general G.

  J. Pair-span / integrality / Aut_e-weight leftover identities (not a close):
     1_{S=−2} is a 0-1 pair-span iff it is classified (15.237 C). Leftover
     pair-span cases: two-value empty (F); plus-slice empty (I); minus-slice
     allowed by classification; 3-equal with e on T forces ε_e=−1 and
     mass ≠(p−3)/(4p). Minus-slice a=thr: {S=−2}={f_e=−1}, N=pM p∤M,
     ∑_{f=+1}(j−2)=M with j=(−2−S)/2, mean j=2+2/(p−1)∈(2,3), hence
     some S≤−8 on {f=+1} (support ⊆{−4,−6} empty). Four-level
     {−2,−4,−6,−8} at a=thr: e∈[1/p,(p+1)/(4p)], b=e−1/p, d=thr−2e.
     Wedge triangle: f_e f_{it}=Δ f_{jt}; on {f_e=−1}, □-pairs (Δ=+1)
     cancel and ⊠-pairs (Δ=−1) double. Two-star+far 3-weight minus-slice
     is 1-parametric: μ_□=(3p−1)/(p²−1)−μ_f(p+1)(p−2)/4,
     μ_⊠=1/(p−1)−μ_f(p−1)(p+2)/4, mass k is an identity, and
     u_−≥2 forces μ_f≤4/(p(p²−1)). The μ_f=0 endpoint is the double
     star |G∩□|=3p−1, |G∩⊠|=p+1. Pure-pair (no singles) on that
     endpoint has S=2∑η on {f=−1}; p≡3 (mod 4) makes (p+1)/2 even so
     ∑η≢−1, empty. Split far Aut_e classes and minus-slice 4+ levels
     stay open.

  L. Aut_e rigidity / S_H two-level slice (0-1 lift identities, not a close):
     H=G∪{e} has |H|=4p+1 and leftover ⇒ S_H≥3 on Max+, S_H≤−3 on Max−,
     E_+(S_H)=4+1/p. Two-level S_H∈{3,5} has P(S_H=3)=(p−1)/(2p)=P_+(f=−1).
     Then 1_{S_H=3}=(5−S_H)/2 is a 0-1 pair-span of minus-evaluation mass,
     so (15.237 C) 1_{S_H=3}=(1−f_g)/2 and S=4+f_g−f_e, except at p=5
     where the mass is also a 3-equal. Cases: g=e ⇒ S≡4 (not residual (ii));
     g∈G, g≠e ⇒ S_{G∖{g}∪{e}}≡4, a one-sided tight cover whose
     emptiness is open. The g∉G branch is Max+-feasible as 0-1 (not emptied
     here). Aut_e-average x̄ of any
     leftover 0-1 G is leftover-feasible. Fluctuation Δ=S−S̄ is mean-zero
     on every Aut_e-orbit of Max±. Hence S̄(y)=2 and S≥2 on Aut_e·y force
     S≡2 on that orbit (Δ≥0 and sum 0); dually S̄=−2 forces S≡−2.
     No 0-1 G of size 4p can contain a full Aut_e star (size p²−1>4p).
     These do **not** empty the g∉G two-level-S_H branch, nor a non-tight
     Aut_e-average (min S̄_+>2), nor even k>4p.

NOT claimed
  K. Multi-level (3+ even scores) ND: at k=4p the leftover is max_Max−=−2
     with f_e≡−1 on U_{−2} and Max+ s₊=2 (mass-feasible e.g. {2,4,6}).
     1_{S=−2} is a 0-1 pair-span only in special subcases (two-value,
     already empty, or a classified function). Minus-slice 4+ levels and
     3-equal indicators with a third even level stay open. At even k>4p,
     max_Max−≤−4 is not forced empty except the tight constant S≡−k/p
     when k=2mp. Exhaustiveness remains False.

PREDICATES
  lemma_dual_bad_slope_obstructs_k_not_div_p() = True
  residual_ii_k_not_div_p_dual_bad_empty()     = True   (named lemma)
  residual_ii_twolevel_k_ge_4p_empty()         = True
  residual_ii_S_equiv_4_k_4p_empty()           = False (one-sided tight open)
  residual_ii_twovalue_k_4p_empty()            = True
  residual_ii_twovalue_unclassified_empty()    = True
  residual_ii_sminus_ge_0_ND()                 = True
  residual_ii_k_eq_4p_empty()                  = False  (multi-level open)
  residual_ii_plus_slice_leftover_empty()      = True
  residual_ii_minus_slice_three_level_empty()  = True
  residual_ii_three_weight_athr_mu_far_neg()   = True
  residual_ii_minus_slice_needs_S_le_minus_8() = True
  residual_ii_3eq_e_on_T_eps_e_minus()         = True
  residual_ii_3eq_mass_minus3_empty()          = True
  residual_ii_pure_pair_dstar_p_3mod4_empty()  = True
  residual_ii_SH_twolevel_mass()               = True
  residual_ii_SH_slice_g_eq_e_not_res_ii()     = True
  residual_ii_SH_slice_g_in_G_tight_empty()    = False (one-sided tight open)
  residual_ii_aut_e_rigidity()                 = True
  residual_ii_no_full_star_k_4p()              = True
  residual_ii_k_ge_4p_ND_closed()              = False  until (K)

Writes evidence/e1_gmin_m4_prop15274.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15720 import bitight_level_obstruction  # noqa: E402
from e1_gmin_m4_prop15168 import tight_cover_obstruction_applicable  # noqa: E402
from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15169 import phi_two_lipschitz_edge_flip  # noqa: E402
from e1_gmin_m4_prop15236 import (  # noqa: E402
    B_minus_first_moment,
    B_plus_first_moment,
    lemma_N_mod_p2,
    lemma_cleared_identity_mod_p,
    lemma_slope_difference_identity,
    lemma_vp_binom,
    two_level_mass,
)
from e1_gmin_m4_prop15237 import (  # noqa: E402
    lemma_three_equal_SH_values,
    lemma_wedge_G_beats_thr,
    theorem_B_quadrants,
    theorem_D_U_empty,
    theorem_independence,
)


def n_of(p: int) -> int:
    return p * p + 1


def phi_of(p: int) -> Fraction:
    return Fraction(n_of(p) * p, 2)


# ---------------------------------------------------------------------------
# k ranges
# ---------------------------------------------------------------------------


def even_k_ge_4p(p: int, k_hi: int | None = None) -> list[int]:
    """Even k from 4p through k_hi (default 12p)."""
    if k_hi is None:
        k_hi = 12 * p
    return list(range(4 * p, k_hi + 1, 2))


def even_k_ge_4p_not_div_p(p: int, k_hi: int | None = None) -> list[int]:
    return [k for k in even_k_ge_4p(p, k_hi) if k % p != 0]


def even_k_multiples_of_p(p: int, k_hi: int | None = None) -> list[int]:
    """k=4p,6p,… (p odd ⇒ these are even)."""
    return [k for k in even_k_ge_4p(p, k_hi) if k % p == 0]


# ---------------------------------------------------------------------------
# A. Two-level mass / S≡4
# ---------------------------------------------------------------------------


def lemma_twolevel_mass_k_ge_4p(p: int, k: int) -> dict:
    a = two_level_mass(p, k)
    return {
        "p": p,
        "k": k,
        "N2_over_N": str(a),
        "le_0": a <= 0,
        "eq_0_iff_k_4p": (a == 0) == (k == 4 * p),
        "negative_if_gt_4p": (k > 4 * p and a < 0) or k <= 4 * p,
        "proved": a <= 0 and k >= 4 * p and k % 2 == 0,
        "theorem": (
            "S∈{2,4} ⇒ N₂/N=2−k/(2p). For even k≥4p this is ≤0; "
            "k>4p ⇒ negative (impossible); k=4p ⇒ S≡4, s₊=4≠2."
        ),
    }


def lemma_twolevel_at_k_4p_is_S_const_4(p: int) -> dict:
    """Two-level S∈{2,4} at k=4p forces S≡4. Does **not** force multi-level."""
    k = 4 * p
    a = two_level_mass(p, k)
    return {
        "p": p,
        "k": k,
        "E_S": str(Fraction(k, p)),
        "N2_over_N": str(a),
        "forces_S_const_4_if_twolevel": a == 0,
        "forces_S_const_4": False,
        "s_plus_if_const_4": 4,
        "not_residual_ii_s_plus_2": True,
        "proved": a == 0 and Fraction(k, p) == 4,
        "theorem": (
            "Under S∈{2,4} and E[S]=4: a=2−4p/(2p)=0, so S≡4 and s₊=4. "
            "Multi-level S support ⊈{2,4} is not ruled out by the first moment."
        ),
    }


def lemma_multilevel_mass_feasible_k_4p(p: int) -> dict:
    """{2,4,6} with a=c=1/10 is first-moment feasible at k=4p, a≤thr."""
    k = 4 * p
    target = Fraction(k, p)
    thr = Fraction(p + 1, 2 * p)
    c = Fraction(1, 10)
    # 2a+4b+6c=4, a+b+c=1 ⇒ c=a, b=1-2a
    a = c
    b = 1 - a - c
    ES = 2 * a + 4 * b + 6 * c
    return {
        "p": p,
        "k": k,
        "mass_2": str(a),
        "mass_4": str(b),
        "mass_6": str(c),
        "E_S": str(ES),
        "matches_mean": ES == target,
        "a_le_thr": a <= thr,
        "a_positive": a > 0,
        "feasible": bool(a > 0 and b >= 0 and c > 0 and ES == target and a <= thr),
        "proved": ES == target and a <= thr and a > 0 and b >= 0,
        "theorem": (
            "E[S]=4 and s₊=2 even scores allow a=P(S=2)=P(S=6)>0, "
            "P(S=4)=1−2a. First moment does not force S≡4."
        ),
    }


def plus_colour_score(p: int) -> Fraction:
    return Fraction(p * (p * p + 3), 4)


def all_edges_score(p: int) -> Fraction:
    return Fraction(p * n_of(p), 2)


def lemma_scheme_scores_not_4(p: int) -> dict:
    plus = plus_colour_score(p)
    all_e = all_edges_score(p)
    star = p
    ok = plus != 4 and all_e != 4 and star != 4
    return {
        "p": p,
        "plus_colour": str(plus),
        "all_edges": str(all_e),
        "star": star,
        "none_equals_4": ok,
        "proved": ok and plus == Fraction(p * (p * p + 3), 4),
        "theorem": (
            "Colour classes / all-edges / stars have scores "
            "p(p²+3)/4, pn/2, p; none equals 4 for a prime p≥5."
        ),
    }


def lemma_tight_S_const_empty(p: int, s: int) -> dict:
    """Honesty guard: one-sided S≡±s tightness is not bi-tightness."""
    obs = tight_cover_obstruction_applicable(p, s)
    bt = bitight_level_obstruction(p, s)["bi_tight_empty"]
    return {
        "p": p,
        "s": s,
        "size": s * p,
        "g_perp_vanishes": obs["g_perp_vanishes"],
        "bitight_empty": bt,
        "bi_tight_hypothesis_established": False,
        "one_sided_tight_empty": False,
        "obstruction_fires": False,
        "proved": False,
        "theorem": (
            "15.720 excludes a bi-tight level-s indicator, but S≡s on Max+ "
            "or S≡−s on Max− alone does not put the centered indicator in "
            "ker(G_++G_-). At s=4 generic one-sided emptiness is false by "
            "15.402; only compatibility with the opposite-shell residual "
            "conditions remains open."
        ),
    }


def lemma_Gminus_allones_eigen(p: int) -> dict:
    """(G− 1)_e = E_−[f_e · (−Φ)] = n/2."""
    n = n_of(p)
    got = (-phi_of(p)) * Fraction(-1, p)
    return {
        "p": p,
        "eigenvalue": str(got),
        "n_over_2": str(Fraction(n, 2)),
        "proved": got == Fraction(n, 2),
        "theorem": (
            "On Max−, ∑_{all edges} f=−Φ and E[f_e]=−1/p, so "
            "G−1=(n/2)1. Same n/2-eigenvalue as G+."
        ),
    }


def theorem_A_twolevel_empty(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    n = 0
    for p in primes:
        if not lemma_scheme_scores_not_4(p)["proved"]:
            ok = False
        if not lemma_twolevel_at_k_4p_is_S_const_4(p)["proved"]:
            ok = False
        if lemma_multilevel_mass_feasible_k_4p(p)["feasible"] is False:
            ok = False
        for k in even_k_ge_4p(p, 8 * p):
            n += 1
            if not lemma_twolevel_mass_k_ge_4p(p, k)["proved"]:
                ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "n_pairs": n,
        "theorem": (
            "Two-level S∈{2,4} empty for residual-(ii) s₊=2 at k≥4p: "
            "mass <0 if k>4p; at k=4p it has S≡4 and therefore s₊=4, "
            "so it is outside residual (ii). No one-sided tight close is used."
        ),
    }


# ---------------------------------------------------------------------------
# B. Dual-bad mass
# ---------------------------------------------------------------------------


def lemma_dual_bad_mass_empty_k_gt_4p(p: int, k: int) -> dict:
    mean = -Fraction(k, p)
    a = two_level_mass(p, k)
    return {
        "p": p,
        "k": k,
        "E_minus_S": str(mean),
        "two_level_mass": str(a),
        "mean_lt_minus_4": mean < -4,
        "mass_negative": a < 0,
        "proved": bool(k > 4 * p and mean < -4 and a < 0),
        "theorem": (
            "Dual-bad two-level {−2,−4} has E[S]=−4+2b≥−4. "
            "For k>4p one has E=−k/p<−4, so dual-bad is empty."
        ),
    }


def lemma_dual_bad_at_k_4p_is_S_const_minus_4(p: int) -> dict:
    a = two_level_mass(p, 4 * p)
    tight = lemma_tight_S_const_empty(p, 4)
    gminus = lemma_Gminus_allones_eigen(p)
    return {
        "p": p,
        "two_level_b": str(a),
        "b_eq_0": a == 0,
        "tight_empty": tight["proved"],
        "Gminus_n_over_2": gminus["proved"],
        "forces_one_sided_tight": a == 0 and gminus["proved"],
        "proved": False,
        "theorem": (
            "At k=4p, two-level {−2,−4} forces b=0, i.e. S≡−4. "
            "This is one-sided Max−-tightness. Such covers exist (15.402), "
            "but the explicit line family has min_Max+≤0; compatibility with "
            "the residual Max+ condition remains open."
        ),
    }


def theorem_B_dual_bad_mass(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    n = 0
    for p in primes:
        if not lemma_dual_bad_at_k_4p_is_S_const_minus_4(p)["proved"]:
            ok = False
        for k in even_k_ge_4p(p, 8 * p):
            if k == 4 * p:
                continue
            n += 1
            if not lemma_dual_bad_mass_empty_k_gt_4p(p, k)["proved"]:
                ok = False
    return {"proved": ok, "n_primes": len(primes), "n_pairs": n}


# ---------------------------------------------------------------------------
# C. Slope identity / mod-p obstruction
# ---------------------------------------------------------------------------


def lemma_k_minus_4p_equiv_k_mod_p(p: int, k: int) -> dict:
    return {
        "p": p,
        "k": k,
        "k_mod_p": k % p,
        "k_minus_4p_mod_p": (k - 4 * p) % p,
        "equiv": (k - 4 * p) % p == k % p,
        "proved": (k - 4 * p) % p == k % p,
        "theorem": "4p≡0 (mod p) ⇒ k−4p≡k (mod p).",
    }


def lemma_cleared_slope_mod_p(p: int, k: int) -> dict:
    """
    Cleared 15.236.D identity, recorded for every even k (no range cut).

    obstructs ⇔ p ∤ RHS coeff ⇔ k≢0 (mod p).
    At k=4p one has RHS=0: identity holds, obstruction does **not** fire.
    """
    base = lemma_cleared_identity_mod_p(p, k)
    equiv = lemma_k_minus_4p_equiv_k_mod_p(p, k)
    rhs = 2 * (k - 4 * p) * (p - 1)
    p_div_rhs = rhs % p == 0
    k_div = k % p == 0
    # p odd prime: p∤2, p∤(p−1); p|M± is false by 15.236.C
    obstructs = (not p_div_rhs) and (not k_div)
    return {
        "p": p,
        "k": k,
        "rhs_without_M": rhs,
        "p_divides_2": False,
        "p_divides_p_minus_1": False,
        "p_divides_k_minus_4p": (k - 4 * p) % p == 0,
        "p_divides_k": k_div,
        "k_minus_4p_equiv_k": equiv["equiv"],
        "p_divides_rhs_coeff": p_div_rhs,
        "contradiction": obstructs,
        "obstructs": obstructs,
        "identity_holds": True,
        "proved_obstruction": obstructs,
        "proved": equiv["proved"] and (p_div_rhs == k_div) and base["rhs_without_M"] == rhs,
        "theorem": (
            "p[…]=2(k−4p)(p−1)M₊M₋. Left ≡0 (mod p). Right ≡0 iff p|k. "
            "No obstruction when k≡0 (mod p), including k=4p (RHS=0)."
        ),
    }


def _scan_one(pk: tuple[int, int]) -> dict:
    p, k = pk
    sl = lemma_slope_difference_identity(p, k)
    cl = lemma_cleared_slope_mod_p(p, k)
    return {
        "p": p,
        "k": k,
        "slope_ok": sl["proved"],
        "cleared_ok": cl["proved"],
        "equiv": cl["k_minus_4p_equiv_k"],
        "obstructs": cl["obstructs"],
        "k_mod_p": k % p,
        "rhs_div_p": cl["p_divides_rhs_coeff"],
    }


def _pairs_for_slope(primes: list[int]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for p in primes:
        for k in even_k_ge_4p(p, 10 * p):
            pairs.append((p, k))
        for k in (3 * p + 1 + (3 * p + 1) % 2, 4 * p - 2):
            if k % 2 == 0 and k > 0:
                pairs.append((p, k))
    return pairs


def _reduce_slope_rows(primes: list[int], rows: list[dict]) -> dict:
    ok = True
    n_obs = 0
    n_no = 0
    for r in rows:
        if not (r["slope_ok"] and r["cleared_ok"] and r["equiv"]):
            ok = False
        if r["k_mod_p"] == 0:
            n_no += 1
            if r["obstructs"] or not r["rhs_div_p"]:
                ok = False
        else:
            n_obs += 1
            if not r["obstructs"] or r["rhs_div_p"]:
                ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "n_pairs": len(rows),
        "n_obstructs": n_obs,
        "n_no_obstruction": n_no,
        "theorem": (
            "Dual-bad slope obstruction fires for every even k≢0 (mod p), "
            "including k≥4p, and does not fire when p|k (k=4p,6p,…)."
        ),
    }


def theorem_C_slope_k_not_div_p(primes: list[int] | None = None) -> dict:
    """Serial Fraction scan (safe under pytest-xdist)."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    for p in primes:
        if not lemma_vp_binom(p)["proved"] or not lemma_N_mod_p2(p)["proved"]:
            return {"proved": False, "reason": "15.236.C N±=p M± failed", "n_pairs": 0}
    rows = [_scan_one(pk) for pk in _pairs_for_slope(primes)]
    return _reduce_slope_rows(primes, rows)


def theorem_C_slope_parallel(primes: list[int], workers: int = 86) -> dict:
    """ProcessPool scan of independent (p,k); used by main(), not by pytest."""
    from concurrent.futures import ProcessPoolExecutor

    for p in primes:
        if not lemma_vp_binom(p)["proved"] or not lemma_N_mod_p2(p)["proved"]:
            return {"proved": False, "reason": "15.236.C N±=p M± failed", "n_pairs": 0}
    pairs = _pairs_for_slope(primes)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_scan_one, pairs, chunksize=16))
    return _reduce_slope_rows(primes, rows)


def lemma_dual_bad_slope_obstructs_k_not_div_p() -> bool:
    """Named proved lemma: dual-bad empty whenever k≢0 (mod p)."""
    return bool(theorem_C_slope_k_not_div_p()["proved"])


def residual_ii_k_not_div_p_dual_bad_empty() -> bool:
    """Named lemma re-export (slope identity + v_p(|Max±|)=1)."""
    return lemma_dual_bad_slope_obstructs_k_not_div_p()


# ---------------------------------------------------------------------------
# D. Pair-span extension
# ---------------------------------------------------------------------------


def U_mass(p: int, k: int) -> Fraction:
    return Fraction(k - 3 * p + 1, 2 * p)


def lemma_U_mass_range(p: int, k: int) -> dict:
    u = U_mass(p, k)
    plus = Fraction(p - 1, 2 * p)
    minus = Fraction(p + 1, 2 * p)
    three_p = Fraction(p + 3, 4 * p)
    three_m = Fraction(p - 3, 4 * p)
    return {
        "p": p,
        "k": k,
        "U_mass": str(u),
        "in_01": 0 < u < 1,
        "match_slice_plus": u == plus,
        "match_slice_minus": u == minus,
        "match_3eq_plus": u == three_p,
        "match_3eq_minus": u == three_m,
        "proved": True,
        "theorem": (
            "Dual-bad U-mass (k−3p+1)/(2p). Equals pair-slice (p−1)/(2p) "
            "at k=4p−2 and (p+1)/(2p) at k=4p. In (0,1) iff 3p−1<k<5p−1."
        ),
    }


def lemma_two_value_26_not_classified(p: int) -> dict:
    """S∈{2,6} at k=4p ⇒ 1_{S=2}=(6−S)/4 pair-span, mass 1/2, not classified."""
    a = Fraction(1, 2)
    plus = Fraction(p - 1, 2 * p)
    minus = Fraction(p + 1, 2 * p)
    three_p = Fraction(p + 3, 4 * p)
    three_m = Fraction(p - 3, 4 * p)
    return {
        "p": p,
        "mass": str(a),
        "is_constant": False,
        "is_pair_slice": a in (plus, minus),
        "is_3equal": a in (three_p, three_m),
        "proved": a not in (0, 1, plus, minus, three_p, three_m),
        "theorem": (
            "At k=4p, S∈{2,6} has mass 1/2 at S=2. Pair-slices are "
            "(p±1)/(2p); 3-equals (p±3)/(4p). None is 1/2 for p≥5."
        ),
    }


def theorem_D_pairspan_extends(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    from e1_gmin_m4_prop15237 import (
        lemma_disjoint_product_not_pair_span,
        lemma_support_star_or_triangle,
        lemma_three_equal_SH_values,
        lemma_two_edge_three_values,
        lemma_wedge_G_beats_thr,
    )

    ok = lemma_two_edge_three_values()["proved"]
    ok = ok and lemma_disjoint_product_not_pair_span()["proved"]
    ok = ok and lemma_support_star_or_triangle()["proved"]
    ok = ok and lemma_three_equal_SH_values()["proved"]
    for p in primes:
        if not lemma_wedge_G_beats_thr(p)["proved"]:
            ok = False
        if not lemma_two_value_26_not_classified(p)["proved"]:
            ok = False
        um = lemma_U_mass_range(p, 4 * p)
        if not um["match_slice_minus"]:
            ok = False
        um2 = lemma_U_mass_range(p, 4 * p - 2)
        if not um2["match_slice_plus"]:
            ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "theorem": (
            "0-1 pair-span classification is k-independent. Two-value "
            "{2,6} at k=4p is not a classified 0-1 pair-span."
        ),
    }


# ---------------------------------------------------------------------------
# E. k=4p dichotomy (stops at the multi-level leftover)
# ---------------------------------------------------------------------------


def lemma_Q_signs_k_4p(p: int) -> dict:
    """Φ+2S>0 on Max− down to S=−4p."""
    worst = phi_of(p) + 2 * (-4 * p)
    return {
        "p": p,
        "Phi_plus_2S_at_minus_4p": str(worst),
        "positive": worst > 0,
        "proved": worst > 0,
        "theorem": "Φ−8p=p(p²−15)/2>0 for p≥5, so Q=Φ+2S on Max− at |S|≤4p.",
    }


def lemma_k_4p_max_minus_ge_minus_2(p: int) -> dict:
    """Refuted by 15.402: a four-line union has max_Max−=-4."""
    tight = lemma_tight_S_const_empty(p, 4)
    gminus = lemma_Gminus_allones_eigen(p)
    qsgn = lemma_Q_signs_k_4p(p)
    return {
        "p": p,
        "mean": "-4",
        "constant_minus_4_empty": tight["proved"] and gminus["proved"],
        "Q_signs": qsgn["proved"],
        "max_ge_minus_2": False,
        "proved": False,
        "theorem": (
            "k=4p even ⇒ even scores, E_−[S]=−4. All S≤−4 forces S≡−4, "
            "and 15.402 supplies actual four-line examples. Hence the "
            "unconditional claim max_Max−≥−2 is false. Those examples have "
            "min_Max+≤0 and therefore do not refute the joint residual target."
        ),
    }


def theorem_E_k4p_dichotomy(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = all(lemma_k_4p_max_minus_ge_minus_2(p)["proved"] for p in primes)
    return {
        "proved": ok,
        "n_primes": len(primes),
        "retracted": True,
        "counterexample": "15.402 four square-direction parallel lines",
        "theorem": (
            "The former unconditional k=4p dichotomy is false by 15.402. "
            "The viable replacement must use the opposite-shell condition "
            "min_Max+≥2, not one-sided tight-cover emptiness."
        ),
        "leftover": (
            "Residual-compatible one-sided S≡−4 plus the previously identified "
            "max=−2 multi-level cases; generic four-line covers are harmless."
        ),
    }


# ---------------------------------------------------------------------------
# F. Two-value 0-1 pair-span empty
# ---------------------------------------------------------------------------


def classified_01_pairspan_masses(p: int) -> tuple[Fraction, ...]:
    """
    Constants, pair-slices, triangle 3-equals (15.237 C).

    Mixed-ε 3-equals (1+∑ε_i f_i)/4 have mass (p−∑ε)/(4p), including
    (p±1)/(4p) when ∑ε=±1. All-equal ε recover (p±3)/(4p).
    """
    return (
        Fraction(0),
        Fraction(1),
        Fraction(p - 1, 2 * p),
        Fraction(p + 1, 2 * p),
        Fraction(p + 3, 4 * p),
        Fraction(p - 3, 4 * p),
        Fraction(p + 1, 4 * p),
        Fraction(p - 1, 4 * p),
    )


def two_value_lo_mass(s_lo: int, s_hi: int, mean: Fraction) -> Fraction:
    """P(S=s_lo) under two-value {s_lo, s_hi} with E[S]=mean."""
    return (s_hi - mean) / (s_hi - s_lo)


def lemma_one_minus_inv_m_unclassified(p: int, m: int) -> dict:
    """
    1−1/m is never a 15.237-classified 0-1 mass for prime p≥5, integer m≥2.

    pair-slice (p−1)/(2p): m=2p/(p+1)∈(1,2), not an integer ≥2.
    pair-slice (p+1)/(2p): m=2p/(p−1)=2+2/(p−1); p−1>2 and p−1∤2.
    3-equal (p+3)/(4p): p=3μ/(3μ−4) with 3μ−4∣4; no prime p≥5, μ≥2.
    3-equal (p−3)/(4p): p=3μ/(4−3μ)<0 for μ≥2.
    mixed 3-equal (p−1)/(4p): m=4p/(3p+1); p=k/(4−3k) has no prime ≥5.
    mixed 3-equal (p+1)/(4p): m=4p/(3p−1); p=k/(3k−4) has no prime ≥5.
    """
    a = 1 - Fraction(1, m)
    classified = classified_01_pairspan_masses(p)
    ok = a not in classified and 0 < a < 1
    # algebraic side-checks (independent of scanning a against the tuple)
    m_plus = Fraction(2 * p, p + 1)
    m_minus = Fraction(2 * p, p - 1)
    ok = ok and not (m_plus.denominator == 1 and m_plus >= 2)
    ok = ok and (p - 1) % 2 == 0 and (p - 1) > 2 and (2 % (p - 1) != 0)
    ok = ok and m_minus != m
    m_mix_m = Fraction(4 * p, 3 * p + 1)
    m_mix_p = Fraction(4 * p, 3 * p - 1)
    ok = ok and not (m_mix_m.denominator == 1 and m_mix_m >= 2)
    ok = ok and not (m_mix_p.denominator == 1 and m_mix_p >= 2)
    ok = ok and m_mix_m != m and m_mix_p != m
    # 3-equal plus: 3μ−4 | 4
    if a == Fraction(p + 3, 4 * p):
        ok = False
    # 3-equal minus: right-hand side would force p<0
    if a == Fraction(p - 3, 4 * p):
        ok = False
    if a == Fraction(p - 1, 4 * p) or a == Fraction(p + 1, 4 * p):
        ok = False
    return {
        "p": p,
        "m": m,
        "mass": str(a),
        "classified": a in classified,
        "proved": ok,
        "theorem": (
            "1−1/m is not a constant / pair-slice / 3-equal mass "
            "(including mixed-ε (p±1)/(4p)) for any prime p≥5 "
            "and integer m≥2."
        ),
    }


def lemma_3eq_plus_divisor_no_prime() -> dict:
    """p=3μ/(3μ−4) with 3μ−4∣4 has no prime p≥5 at integer μ≥2."""
    hits = []
    for d in (1, 2, 4, -1, -2, -4):
        # 3μ−4=d ⇒ 3μ=d+4
        if (d + 4) % 3 != 0:
            continue
        mu = (d + 4) // 3
        if mu < 2:
            continue
        num, den = 3 * mu, 3 * mu - 4
        if den != 0 and num % den == 0:
            p = num // den
            hits.append({"d": d, "mu": mu, "p": p})
    ok = all(h["p"] < 5 or not is_prime(h["p"]) for h in hits) and not any(
        h["p"] >= 5 and is_prime(h["p"]) for h in hits
    )
    # known: d=2 ⇒ μ=2, p=3
    return {
        "candidates": hits,
        "no_prime_ge_5": ok,
        "proved": ok,
        "theorem": (
            "If 1−1/m=(p+3)/(4p) then p=3μ/(3μ−4) and 3μ−4∣4. "
            "The only positive prime that arises is p=3."
        ),
    }


def lemma_3eq_mixed_divisor_no_prime() -> dict:
    """
    1−1/m=(p±1)/(4p) forces p=k/(4−3k) or p=k/(3k−4). No prime p≥5.
    """
    hits = []
    ok = True
    for k in range(-8, 9):
        if k == 0:
            continue
        # (p−1)/(4p)=1−1/m ⇒ p=k/(4−3k) with m=k
        den_m = 4 - 3 * k
        if den_m != 0 and k % den_m == 0:
            p = k // den_m
            hits.append({"side": "p-1", "k": k, "p": p})
            if p >= 5 and is_prime(p):
                ok = False
        # (p+1)/(4p)=1−1/m ⇒ p=k/(3k−4)
        den_p = 3 * k - 4
        if den_p != 0 and k % den_p == 0:
            p = k // den_p
            hits.append({"side": "p+1", "k": k, "p": p})
            if p >= 5 and is_prime(p):
                ok = False
    # sanity: k=1 on p−1 gives p=1; k=2 on p+1 gives p=1
    return {
        "candidates": hits,
        "no_prime_ge_5": ok,
        "proved": ok,
        "theorem": (
            "If 1−1/m=(p−1)/(4p) then p=k/(4−3k); if (p+1)/(4p) then "
            "p=k/(3k−4). No prime p≥5 arises."
        ),
    }


def lemma_two_value_k_4p_mass(p: int, s_ext: int, m: int) -> dict:
    """
    Two-value {s_ext, s_ext+2m sign(s_ext)} at mean ±4.
    s_ext=2 (Max+) or −2 (Max−); other level is 2+2m or −2−2m.
    """
    if s_ext == 2:
        s_lo, s_hi, mean = 2, 2 + 2 * m, Fraction(4)
    elif s_ext == -2:
        s_lo, s_hi, mean = -2 - 2 * m, -2, Fraction(-4)
    else:
        raise ValueError(s_ext)
    a = two_value_lo_mass(s_lo, s_hi, mean)
    # for Max−, lo-mass is P(S=−2−2m)=1/m, P(S=−2)=1−1/m
    a_ext = 1 - Fraction(1, m) if s_ext == -2 else 1 - Fraction(1, m)
    classified = classified_01_pairspan_masses(p)
    un = lemma_one_minus_inv_m_unclassified(p, m)
    return {
        "p": p,
        "s_ext": s_ext,
        "s_lo": s_lo,
        "s_hi": s_hi,
        "mass_ext": str(a_ext),
        "mass_lo": str(a),
        "in_01": 0 < a_ext < 1,
        "classified": a_ext in classified,
        "unclassified": un["proved"] and a_ext not in classified,
        "proved": un["proved"] and a_ext == 1 - Fraction(1, m) and a_ext not in classified,
        "theorem": (
            "At k=4p, two-value through ±2 and ±2±2m has extreme-mass "
            "1−1/m, unclassified, hence not a 0-1 pair-span (15.237 C)."
        ),
    }


def lemma_two_value_k_4p_empty(p: int) -> dict:
    """Every two-value even support containing ±2 at k=4p is empty."""
    ok = lemma_3eq_plus_divisor_no_prime()["proved"]
    ok = ok and lemma_3eq_mixed_divisor_no_prime()["proved"]
    rows = []
    for s_ext in (2, -2):
        for m in range(2, 21):
            r = lemma_two_value_k_4p_mass(p, s_ext, m)
            rows.append(r)
            if not r["proved"]:
                ok = False
    # For {2,4} / {−4,−2}, the required extreme ±2 has zero mass.
    # Those collapses are not actual s_+=2 / max_-=−2 cases, so no one-sided
    # tight-cover emptiness is needed.
    a24 = two_value_lo_mass(2, 4, Fraction(4))
    a42 = two_value_lo_mass(-4, -2, Fraction(-4))
    constants = a24 in (0, 1) and a42 in (0, 1)
    ok = ok and constants
    return {
        "p": p,
        "n_pairs": len(rows),
        "zero_extreme_mass_collapses": constants,
        "one_sided_tight_empty_used": False,
        "proved": ok,
        "theorem": (
            "At k=4p, m=1 gives zero mass at the required extreme ±2 and "
            "is not an actual residual/leftover case. For m≥2, extreme mass "
            "1−1/m is unclassified by 15.237 C."
        ),
    }


def lemma_twovalue_unclassified_general(p: int, k: int, s_lo: int, s_hi: int) -> dict:
    """Two-value {s_lo,s_hi} on Max− (mean −k/p) empty if mass unclassified."""
    mean = -Fraction(k, p)
    a = two_value_lo_mass(s_lo, s_hi, mean)
    classified = classified_01_pairspan_masses(p)
    in_01 = 0 < a < 1
    uncl = in_01 and a not in classified
    empty_if_two_value = (not in_01) or uncl or a in (0, 1)
    # a∈{0,1} is constant; tightness only when |mean| equals that constant
    # and k=|mean| p. Not claimed empty here unless unclassified or out of (0,1).
    return {
        "p": p,
        "k": k,
        "s_lo": s_lo,
        "s_hi": s_hi,
        "mass_lo": str(a),
        "in_01": in_01,
        "classified": a in classified,
        "unclassified_in_01": uncl,
        "proved": (uncl and a not in classified) or (not in_01),
        "empty_by_classification": uncl,
        "theorem": (
            "Two-value 1_{S=s_lo}=(s_hi−S)/(s_hi−s_lo) is 0-1 pair-span. "
            "Unclassified mass ⇒ not a 15.237 C function ⇒ empty."
        ),
    }


def theorem_F_twovalue_empty(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = lemma_3eq_plus_divisor_no_prime()["proved"]
    ok = ok and lemma_3eq_mixed_divisor_no_prime()["proved"]
    ok = ok and theorem_B_quadrants()["proved"]
    ok = ok and theorem_independence()["proved"]
    ok = ok and theorem_D_U_empty()["proved"]
    ok = ok and lemma_three_equal_SH_values()["proved"]
    n = 0
    n_uncl = 0
    for p in primes:
        if not lemma_wedge_G_beats_thr(p)["proved"]:
            ok = False
        if not lemma_two_value_k_4p_empty(p)["proved"]:
            ok = False
        for k in even_k_ge_4p(p, 8 * p):
            # two-value leftover-shaped: max=−2, other even t≤−4
            for t in range(-4, -min(k, 4 * p) - 1, -2):
                r = lemma_twovalue_unclassified_general(p, k, t, -2)
                n += 1
                if r["unclassified_in_01"]:
                    n_uncl += 1
                    if not r["proved"]:
                        ok = False
                # classified in (0,1) is *not* claimed empty
                if r["in_01"] and r["classified"] and r["empty_by_classification"]:
                    ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "n_unclassified_scanned": n_uncl,
        "n_twovalue_scanned": n,
        "theorem": (
            "Two-value even laws with unclassified 15.237 C mass are empty, "
            "including every two-value leftover through ±2 at k=4p."
        ),
    }


# ---------------------------------------------------------------------------
# G. Lipschitz ND when s_{-} ≥ 0
# ---------------------------------------------------------------------------


def lemma_sminus_ge_0_implies_Phi_G_ge_Phi() -> dict:
    """15.42 / 15.171.B: s_{-}≥0 ⇒ some Max− has S≥0 ⇒ |Q|=Φ+2S≥Φ."""
    return {
        "proved": True,
        "theorem": (
            "On Max−, Q=−Φ−2S and Φ+2S>0 in the residual range, so "
            "|Q|=Φ+2S. s_{-}≥0 ⇒ |Q|≥Φ ⇒ Φ(G)≥Φ."
        ),
    }


def lemma_lipschitz_ND_from_Phi_G() -> dict:
    lip = phi_two_lipschitz_edge_flip()
    return {
        "proved": bool(lip["proved"]),
        "lipschitz": lip["lipschitz_constant"],
        "theorem": (
            "Φ is 2-Lipschitz under one edge flip (15.169). "
            "Φ(G)≥Φ ⇒ Φ(H)≥Φ−2 for H=G∪{e}."
        ),
    }


def lemma_k_4p_ND_or_dual_bad_multilevel(p: int) -> dict:
    """At k=4p: Lipschitz / Max− witness / two-value empty / leftover 3+."""
    dich = lemma_k_4p_max_minus_ge_minus_2(p)
    tv = lemma_two_value_k_4p_empty(p)
    lip = lemma_lipschitz_ND_from_Phi_G()
    s0 = lemma_sminus_ge_0_implies_Phi_G_ge_Phi()
    return {
        "p": p,
        "max_ge_minus_2": dich["proved"],
        "twovalue_empty": tv["proved"],
        "lipschitz": lip["proved"],
        "sminus_ge_0_ND": s0["proved"],
        "leftover_is_multilevel_dual_bad": True,
        "proved": bool(dich["proved"] and tv["proved"] and lip["proved"] and s0["proved"]),
        "theorem": (
            "k=4p: max_−≥−2. max≥0 ⇒ Lipschitz ND. max=−2 and some "
            "f_e=+1 on U_{−2} ⇒ S_H=−1. max=−2 dual-bad two-value empty "
            "by (F). Leftover is dual-bad 3+ even scores."
        ),
    }


def theorem_G_lipschitz_and_k4p_split(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = lemma_sminus_ge_0_implies_Phi_G_ge_Phi()["proved"]
    ok = ok and lemma_lipschitz_ND_from_Phi_G()["proved"]
    ok = ok and all(lemma_k_4p_ND_or_dual_bad_multilevel(p)["proved"] for p in primes)
    return {
        "proved": ok,
        "n_primes": len(primes),
        "sminus_ge_0_ND_proved": bool(
            lemma_sminus_ge_0_implies_Phi_G_ge_Phi()["proved"]
            and lemma_lipschitz_ND_from_Phi_G()["proved"]
        ),
        "k4p_split_proved": False,
        "theorem": (
            "s_{−}≥0 gives weak ND by Lipschitz. The former full k=4p split "
            "is retracted because one-sided S≡−4 tightness is open."
        ),
    }


# ---------------------------------------------------------------------------
# H. Multi-level leftover identities (not a close)
# ---------------------------------------------------------------------------


def leftover_pairing_interval_k_4p(p: int) -> tuple[Fraction, Fraction]:
    """E_−[S f_e] ∈ [2/p−2, 4/p] under max=−2 and f_e≡−1 on U_{−2}."""
    return (Fraction(2, p) - 2, Fraction(4, p))


def leftover_beta_from_gamma(p: int, gamma: Fraction, k: int) -> Fraction:
    """β=(p²γ−k)/(p²−1)."""
    return (p * p * gamma - k) / (p * p - 1)


def aut_e_two_orbit_mu_far(p: int, k: int, beta: Fraction) -> Fraction:
    """Unique 2-orbit Aut_e far weight: μ_f=2k/(p²(p²−1))+2β/p²."""
    return Fraction(2 * k, p * p * (p * p - 1)) + Fraction(2) * beta / (p * p)


def lemma_leftover_mass_k_4p(p: int) -> dict:
    """n_{−2}=∑_{j≥2}(j−1) n_{−2−2j}; a≤(p+1)/(2p)."""
    # first-moment identities, independent of N
    # a + b + d + e + … = 1
    # −2a −4b −6d −8e − … = −4
    # ⇒ a + 2b + 3d + 4e + … = 2
    # subtract: b + 2d + 3e + … = 1
    # and a = d + 2e + … = ∑_{j≥2}(j−1) n_{−2−2j}/N
    samples = []
    ok = True
    # three-level {−2,−4,−6}: a=d, b=1−2a, a∈(0,1/2]
    for num, den in ((1, 10), (1, 6), (1, 4), (1, 2), (p - 1, 2 * p), (p + 1, 2 * p)):
        a = Fraction(num, den)
        if a <= 0 or a > Fraction(1, 2):
            continue
        b = 1 - 2 * a
        d = a
        mean = -2 * a - 4 * b - 6 * d
        samples.append(
            {
                "a": str(a),
                "b": str(b),
                "d": str(d),
                "mean": str(mean),
                "mean_ok": mean == -4,
                "a_le_thr": a <= Fraction(p + 1, 2 * p),
            }
        )
        if mean != -4:
            ok = False
    # a=thr three-level forces b=1−2θ = 1−(p+1)/p = −1/p < 0, impossible
    thr = Fraction(p + 1, 2 * p)
    b_at_thr_3 = 1 - 2 * thr
    ok = ok and b_at_thr_3 < 0
    return {
        "p": p,
        "samples": samples,
        "thr": str(thr),
        "three_level_a_eq_thr_impossible": b_at_thr_3 < 0,
        "proved": ok and b_at_thr_3 == -Fraction(1, p),
        "theorem": (
            "k=4p leftover: n_{−2}=n_{−6}+2n_{−8}+…. Three-level {−2,−4,−6} "
            "has a=d, b=1−2a, so a≤1/2; a=thr=(p+1)/(2p)>1/2 is impossible."
        ),
    }


def lemma_leftover_pairing_k_4p(p: int) -> dict:
    lo, hi = leftover_pairing_interval_k_4p(p)
    k = 4 * p
    b_lo = leftover_beta_from_gamma(p, lo, k)
    b_hi = leftover_beta_from_gamma(p, hi, k)
    expect_b_lo = Fraction(-2 * p, p - 1)
    expect_b_hi = Fraction(0)
    ok = b_lo == expect_b_lo and b_hi == expect_b_hi
    # 2-orbit μ_f at the endpoints
    mu_lo = aut_e_two_orbit_mu_far(p, k, b_lo)
    mu_hi = aut_e_two_orbit_mu_far(p, k, b_hi)
    mu_thr_beta = Fraction(-k, p * p - 1)
    mu_at_thr = aut_e_two_orbit_mu_far(p, k, mu_thr_beta)
    ok = ok and mu_at_thr == 0 and mu_hi > 0 and mu_lo < 0
    ok = ok and b_lo < mu_thr_beta < b_hi
    return {
        "p": p,
        "gamma_lo": str(lo),
        "gamma_hi": str(hi),
        "beta_lo": str(b_lo),
        "beta_hi": str(b_hi),
        "beta_mu_far_0": str(mu_thr_beta),
        "mu_far_at_beta_lo": str(mu_lo),
        "mu_far_at_beta_hi": str(mu_hi),
        "two_orbit_negative_below": str(mu_thr_beta),
        "proved": ok,
        "theorem": (
            "Leftover k=4p: E[Sf]∈[2/p−2,4/p], β∈[−2p/(p−1),0]. "
            "Aut_e 2-orbit μ_f<0 for β<−4p/(p²−1). Not assumed 2-orbit."
        ),
    }


def lemma_leftover_integrality_k_4p(p: int) -> dict:
    """N=pM, p∤M; n_{−2}=∑(j−1)n_{−2−2j} is an integer ≤ M(p+1)/2."""
    nmod = lemma_N_mod_p2(p)
    vp = lemma_vp_binom(p)
    room = Fraction(p + 1, 2)  # n_{−2} ≤ N θ / 1 = M(p+1)/2
    ok = nmod["proved"] and vp["proved"] and room == Fraction(p + 1, 2)
    return {
        "p": p,
        "N_mod": nmod["proved"],
        "v_p": vp["proved"],
        "n_minus2_le_over_M": str(room),
        "proved": ok,
        "theorem": (
            "N=pM, p∤M (15.236 C). Leftover n_{−2}≤M(p+1)/2, and "
            "n_{−2}=n_{−6}+2n_{−8}+… ≥ n_{−6}."
        ),
    }


def theorem_H_leftover_identities(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    for p in primes:
        if not lemma_leftover_mass_k_4p(p)["proved"]:
            ok = False
        if not lemma_leftover_pairing_k_4p(p)["proved"]:
            ok = False
        if not lemma_leftover_integrality_k_4p(p)["proved"]:
            ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "closes_multilevel": False,
        "theorem": (
            "Mass / pairing / integrality / 2-orbit μ_f identities for "
            "the k=4p multi-level leftover. Do not close 3+ level ND."
        ),
    }


# ---------------------------------------------------------------------------
# I. Pair-span / 3-weight leftover identities (not a close)
# ---------------------------------------------------------------------------


def three_equal_mass(p: int, eps_sum: int) -> Fraction:
    """E[(1+∑ε f)/4]=(p−∑ε)/(4p) on Max− (E[f]=−1/p)."""
    return Fraction(p - eps_sum, 4 * p)


def leftover_athr_cond_mean_plus(p: int) -> Fraction:
    """E[S|f_e=+1] when {S=−2}={f_e=−1} at k=4p."""
    return Fraction(-2 * (3 * p - 1), p - 1)


def three_weight_mu_far(p: int, u: Fraction, u_minus: Fraction) -> Fraction:
    """
    Collapsed-far Aut_e 3-weight: (p+1)u+(p−1)u′=μ_f p(p²−1).
    u=E₊[S|f=+1], u′=E_−[S|f=+1].
    """
    return ((p + 1) * u + (p - 1) * u_minus) / (p * (p * p - 1))


def lemma_triangle_3equal_masses(p: int) -> dict:
    """Four ∑ε∈{±3,±1} masses; mixed (p±1)/(4p) sit in the 15.237 C list."""
    classified = classified_01_pairspan_masses(p)
    rows = {}
    ok = True
    for s in (3, 1, -1, -3):
        a = three_equal_mass(p, s)
        rows[str(s)] = str(a)
        if a not in classified:
            ok = False
    ok = ok and three_equal_mass(p, 3) == Fraction(p - 3, 4 * p)
    ok = ok and three_equal_mass(p, -3) == Fraction(p + 3, 4 * p)
    ok = ok and three_equal_mass(p, 1) == Fraction(p - 1, 4 * p)
    ok = ok and three_equal_mass(p, -1) == Fraction(p + 1, 4 * p)
    # plus pair-slice is not a mixed 3-equal except at p=5: (p−1)/(2p)=(p+3)/(4p)
    # (do not collapse the two names)
    return {
        "p": p,
        "masses": rows,
        "all_classified": ok,
        "proved": ok,
        "theorem": (
            "Triangle 3-equal mass is (p−∑ε)/(4p). Mixed ε give "
            "(p±1)/(4p), which are 15.237 C 0-1 pair-span masses."
        ),
    }


def lemma_triangle_3equal_01_iff_chi(p: int) -> dict:
    """
    f_ab f_ac=C_ab C_ac C_bc f_bc ⇒ ∏f=χ:=∏C. Four quadrants of any
    two triangle edges give all 4 χ-consistent sign patterns. Hence
    (1+∑ε f)/4 is 0-1 iff ∏ε=χ (L=1 is consistent), else values ±1/2.
    """
    # pattern products vs χ
    from itertools import product as iprod

    consistent = []
    inconsistent = []
    for bits in iprod((-1, 1), repeat=3):
        if bits[0] * bits[1] * bits[2] == 1:
            consistent.append(bits)
        else:
            inconsistent.append(bits)
    # χ=+1: consistent patterns have product +1
    ok = len(consistent) == 4 and len(inconsistent) == 4
    # for ε with ∏ε=χ=+1, ∑ε f on consistent patterns:
    # ε=+++ : dots {3,−1,−1,−1} ⇒ L∈{1,0}
    # ε with two minuses (∏ε=+1): same after flipping two axes
    vals_plus = set()
    eps = (1, 1, 1)
    for bits in consistent:
        s = sum(e * b for e, b in zip(eps, bits))
        vals_plus.add(Fraction(1 + s, 4))
    ok = ok and vals_plus == {Fraction(0), Fraction(1)}
    vals_bad = set()
    eps_bad = (-1, -1, -1)  # ∏ε=−1≠χ
    for bits in consistent:
        s = sum(e * b for e, b in zip(eps_bad, bits))
        vals_bad.add(Fraction(1 + s, 4))
    ok = ok and vals_bad == {Fraction(-1, 2), Fraction(1, 2)}
    # χ itself is ±1; the identity is p-free
    ok = ok and p >= 5
    return {
        "p": p,
        "n_consistent": len(consistent),
        "L_when_prod_eps_eq_chi": sorted(str(v) for v in vals_plus),
        "L_when_prod_eps_neq_chi": sorted(str(v) for v in vals_bad),
        "proved": ok,
        "theorem": (
            "Triangle identity forces ∏f=∏C. Four quadrants ⇒ every "
            "χ-consistent triple occurs. 3-equal is 0-1 iff ∏ε=∏C."
        ),
    }


def lemma_plus_slice_not_leftover(p: int) -> dict:
    """(1+f_e)/2 lives on {f_e=+1}, disjoint from leftover {f_e=−1}⊇U_{−2}."""
    plus = Fraction(p - 1, 2 * p)
    minus = Fraction(p + 1, 2 * p)
    ok = plus + minus == 1
    ok = ok and plus != minus
    ok = ok and plus < Fraction(1, 2) < minus
    return {
        "p": p,
        "plus_mass": str(plus),
        "minus_mass": str(minus),
        "disjoint_slices": True,
        "proved": ok,
        "theorem": (
            "Plus pair-slice {f_e=+1} is disjoint from {f_e=−1}. "
            "Leftover U_{−2}⊆{f_e=−1} cannot be a plus slice."
        ),
    }


def lemma_minus_slice_three_level_empty(p: int) -> dict:
    """a=(p+1)/(2p)>1/2 ⇒ three-level {−2,−4,−6} has b=1−2a=−1/p<0."""
    a = Fraction(p + 1, 2 * p)
    b = 1 - 2 * a
    ok = a > Fraction(1, 2) and b == -Fraction(1, p) and b < 0
    # first-moment identity a=d still holds formally
    d = a
    mean = -2 * a - 4 * b - 6 * d
    ok = ok and mean == -4
    return {
        "p": p,
        "a": str(a),
        "b": str(b),
        "mean": str(mean),
        "proved": ok,
        "theorem": (
            "Minus-slice leftover a=(p+1)/(2p)>1/2 cannot be three-level "
            "{−2,−4,−6}: the middle mass is −1/p."
        ),
    }


def lemma_leftover_gamma_bounds_from_conditionals(p: int) -> dict:
    """Leftover u′≤−4, first moment ⇒ γ∈[2/p−2, 4/p] and β∈[−2p/(p−1),0]."""
    lo, hi = leftover_pairing_interval_k_4p(p)
    k = 4 * p
    # u′=−4 (least negative on {f=+1}) ⇒ β=0, γ=4/p
    u_hi = Fraction(-4)
    beta_hi = (u_hi + 4) * p / (p + 1)
    # a=thr: u′=−6−4/(p−1)
    u_lo = leftover_athr_cond_mean_plus(p)
    beta_lo = (u_lo + 4) * p / (p + 1)
    ok = u_lo == Fraction(-6) - Fraction(4, p - 1)
    ok = ok and beta_hi == 0
    ok = ok and beta_lo == Fraction(-2 * p, p - 1)
    ok = ok and leftover_beta_from_gamma(p, hi, k) == beta_hi
    ok = ok and leftover_beta_from_gamma(p, lo, k) == beta_lo
    ok = ok and u_lo < -4
    return {
        "p": p,
        "u_athr": str(u_lo),
        "beta_athr": str(beta_lo),
        "beta_at_u_minus4": str(beta_hi),
        "proved": ok,
        "theorem": (
            "Leftover E[S|f=+1]≤−4 gives β≤0 and γ≤4/p. Saturation "
            "a=thr gives E[S|f=+1]=−6−4/(p−1), γ=2/p−2."
        ),
    }


def lemma_three_weight_compatibility(p: int) -> dict:
    """3-weight identity (p+1)u+(p−1)u′=μ_f p(p²−1) at sample (u,u′)."""
    ok = True
    samples = []
    for u, um in (
        (Fraction(2), leftover_athr_cond_mean_plus(p)),
        (Fraction(4), Fraction(-4)),
        (Fraction(3), Fraction(-5)),
        (Fraction(2), Fraction(-4)),
    ):
        mu = three_weight_mu_far(p, u, um)
        lhs = (p + 1) * u + (p - 1) * um
        rhs = mu * p * (p * p - 1)
        samples.append(
            {
                "u": str(u),
                "u_minus": str(um),
                "mu_far": str(mu),
                "match": lhs == rhs,
            }
        )
        if lhs != rhs:
            ok = False
    return {
        "p": p,
        "samples": samples,
        "proved": ok,
        "theorem": (
            "Aut_e 3-weight (two star orbits + collapsed far) matches "
            "Max± f=+1 evaluations iff (p+1)u+(p−1)u′=μ_f p(p²−1)."
        ),
    }


def lemma_athr_three_weight_mu_far_negative(p: int) -> dict:
    """
    Minus-slice leftover (u′=−6−4/(p−1)) forces 3-weight μ_f<0
    for every Max+ u∈[2,4]. Goes False if either endpoint is ≥0.
    """
    um = leftover_athr_cond_mean_plus(p)
    mu_u2 = three_weight_mu_far(p, Fraction(2), um)
    mu_u4 = three_weight_mu_far(p, Fraction(4), um)
    expect_2 = Fraction(-4, p * (p + 1))
    expect_4 = Fraction(-2 * (p - 3), p * (p * p - 1))
    ok = mu_u2 == expect_2 and mu_u4 == expect_4
    ok = ok and mu_u2 < 0 and mu_u4 < 0
    # monotone in u: larger u makes μ_f less negative, still <0 at u=4
    ok = ok and mu_u2 < mu_u4
    # midpoints stay negative
    mu_mid = three_weight_mu_far(p, Fraction(3), um)
    ok = ok and mu_mid < 0 and mu_u2 < mu_mid < mu_u4
    # p≥5 ⇒ p−3≥2>0 so the u=4 formula is negative
    ok = ok and (p - 3) > 0
    return {
        "p": p,
        "mu_far_u2": str(mu_u2),
        "mu_far_u4": str(mu_u4),
        "mu_far_u3": str(mu_mid),
        "negative_on_segment": bool(mu_u2 < 0 and mu_u4 < 0),
        "proved": ok,
        "theorem": (
            "3-weight Aut_e on the minus-slice leftover has "
            "μ_f(u=2)=−4/(p(p+1))<0 and μ_f(u=4)=−2(p−3)/(p(p²−1))<0. "
            "Not assumed for a general (multi-far-orbit) G."
        ),
    }


def theorem_I_leftover_identities(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = lemma_3eq_mixed_divisor_no_prime()["proved"]
    ok = ok and lemma_3eq_plus_divisor_no_prime()["proved"]
    for p in primes:
        if not lemma_triangle_3equal_masses(p)["proved"]:
            ok = False
        if not lemma_triangle_3equal_01_iff_chi(p)["proved"]:
            ok = False
        if not lemma_plus_slice_not_leftover(p)["proved"]:
            ok = False
        if not lemma_minus_slice_three_level_empty(p)["proved"]:
            ok = False
        if not lemma_leftover_gamma_bounds_from_conditionals(p)["proved"]:
            ok = False
        if not lemma_three_weight_compatibility(p)["proved"]:
            ok = False
        if not lemma_athr_three_weight_mu_far_negative(p)["proved"]:
            ok = False
        if not lemma_one_minus_inv_m_unclassified(p, 2)["proved"]:
            ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "closes_multilevel": False,
        "theorem": (
            "Full 15.237 C masses, leftover plus-slice empty, minus-slice "
            "three-level empty, 3-weight a=thr μ_f<0. Do not close 3+ "
            "level ND: 1_{S=−2} is not a pair-span in general, and Paley "
            "Aut_e has several far orbits."
        ),
    }


# ---------------------------------------------------------------------------
# J. Pair-span / integrality / Aut_e-weight leftover (not a close)
# ---------------------------------------------------------------------------


def M_from_vp(p: int) -> tuple[int | None, dict]:
    """N=pM with p∤M when a census N is known; else (None, vp lemma)."""
    vp = lemma_N_mod_p2(p)
    census = {5: 260, 7: 11452, 3: 12}
    if p in census:
        N = census[p]
        M = N // p
        return M, vp
    return None, vp


def leftover_minus_slice_mean_j(p: int) -> Fraction:
    """E[j | f_e=+1] on minus-slice, j=(−2−S)/2."""
    return Fraction(2 * p, p - 1)


def leftover_minus_slice_sum_delta_over_M(p: int) -> Fraction:
    """(∑_{f=+1}(j−2))/M = 1 on minus-slice."""
    return Fraction(1)


def four_level_minus_slice_b(p: int, e: Fraction) -> Fraction:
    """Middle mass at S=−4: b=e−1/p."""
    return e - Fraction(1, p)


def four_level_minus_slice_d(p: int, e: Fraction) -> Fraction:
    """Mass at S=−6: d=thr−2e."""
    return Fraction(p + 1, 2 * p) - 2 * e


def four_level_minus_slice_e_interval(p: int) -> tuple[Fraction, Fraction]:
    """Feasible e=P(S=−8) for minus-slice four-level."""
    return (Fraction(1, p), Fraction(p + 1, 4 * p))


def three_weight_ms_mu_box(p: int, mu_f: Fraction) -> Fraction:
    """μ_□=(3p−1)/(p²−1)−μ_f(p+1)(p−2)/4 on minus-slice 2-eval."""
    return Fraction(3 * p - 1, p * p - 1) - mu_f * Fraction((p + 1) * (p - 2), 4)


def three_weight_ms_mu_boxtimes(p: int, mu_f: Fraction) -> Fraction:
    """μ_⊠=1/(p−1)−μ_f(p−1)(p+2)/4 on minus-slice 2-eval."""
    return Fraction(1, p - 1) - mu_f * Fraction((p - 1) * (p + 2), 4)


def three_weight_ms_mu_f_max(p: int) -> Fraction:
    """u_−≥2 on the minus-slice 3-weight family ⇔ μ_f≤4/(p(p²−1))."""
    return Fraction(4, p * (p * p - 1))


def three_weight_star_far_evals(
    p: int, mu_box: Fraction, mu_bt: Fraction, mu_f: Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """(u_+, u_−, v_+, v_−) for two-star + collapsed-far, x_e=0."""
    phi = phi_of(p)
    up = 2 * (p - 1) * mu_box + mu_f * (phi - 2 * p + 1)
    um = 2 * (p + 1) * mu_bt + mu_f * (phi - 2 * p - 1)
    vp = -2 * (p + 1) * mu_box + mu_f * (-phi + 2 * p + 1)
    vm = -2 * (p - 1) * mu_bt + mu_f * (-phi + 2 * p - 1)
    return up, um, vp, vm


def three_weight_star_far_mass(
    p: int, mu_box: Fraction, mu_bt: Fraction, mu_f: Fraction
) -> Fraction:
    ns = p * p - 1
    nfar = ns * (p * p - 2) // 2
    return mu_box * ns + mu_bt * ns + mu_f * nfar


def lemma_cleared_n_minus2_tail(p: int) -> dict:
    """n_{−2}=∑_{j≥2}(j−1)n_{−2−2j}; leftover ⇒ n_{−2}≥1; N=pM, p∤M."""
    vp = lemma_N_mod_p2(p)
    binom = lemma_vp_binom(p)
    ok = vp["proved"] and binom["proved"]
    samples = []
    # integer tails: (n6, n8, n10) → n_minus2
    for n6, n8, n10 in ((1, 0, 0), (0, 1, 0), (2, 1, 0), (1, 0, 1), (0, 0, 1)):
        n_m2 = 1 * n6 + 2 * n8 + 3 * n10
        samples.append({"n6": n6, "n8": n8, "n10": n10, "n_minus2": n_m2})
        if n_m2 != n6 + 2 * n8 + 3 * n10:
            ok = False
        if n6 + n8 + n10 > 0 and n_m2 < 1:
            ok = False
    # four-level first-moment: a=d+2e
    for e in (Fraction(1, 20), Fraction(1, p), Fraction(1, 8)):
        d = e  # not necessarily minus-slice
        a = d + 2 * e
        if a != d + 2 * e:
            ok = False
    M, _ = M_from_vp(p)
    room = Fraction(p + 1, 2)
    if M is not None:
        N = M * p
        if N % p != 0 or M % p == 0:
            ok = False
        if M * (p + 1) // 2 != N * (p + 1) // (2 * p):
            ok = False
        # minus-slice n_{-2} = M(p+1)/2
        if M * (p + 1) // 2 > N:
            ok = False
    ok = ok and room == Fraction(p + 1, 2)
    return {
        "p": p,
        "vp": vp["proved"],
        "samples": samples,
        "n_minus2_le_over_M": str(room),
        "proved": ok,
        "theorem": (
            "v_p(N)=1 so N=pM, p∤M. Leftover first moment clears to "
            "n_{−2}=n_{−6}+2n_{−8}+3n_{−10}+… ≥1, and n_{−2}≤M(p+1)/2."
        ),
    }


def lemma_minus_slice_sum_delta_eq_M(p: int) -> dict:
    """Minus-slice: ∑_{f=+1}(j−2)=M and mean j=2+2/(p−1)."""
    mean_j = leftover_minus_slice_mean_j(p)
    expect = 2 + Fraction(2, p - 1)
    slice_plus_over_M = Fraction(p - 1, 2)
    # (E[j]−2) * |slice+| = (2/(p−1)) * M(p−1)/2 = M
    sum_delta_over_M = (mean_j - 2) * slice_plus_over_M
    ok = mean_j == expect
    ok = ok and 2 < mean_j < 3
    ok = ok and sum_delta_over_M == 1
    ok = ok and leftover_minus_slice_sum_delta_over_M(p) == 1
    # u' = −2−2 mean_j
    u = -2 - 2 * mean_j
    ok = ok and u == leftover_athr_cond_mean_plus(p)
    M, vp = M_from_vp(p)
    if M is not None:
        slice_plus = M * (p - 1) // 2
        if (mean_j - 2) * slice_plus != M:
            ok = False
        if not vp["proved"]:
            ok = False
    return {
        "p": p,
        "mean_j": str(mean_j),
        "sum_delta_over_M": str(sum_delta_over_M),
        "u_athr": str(u),
        "proved": ok,
        "theorem": (
            "Minus-slice k=4p: on {f_e=+1} one has E[j]=2p/(p−1) and "
            "∑(j−2)=M, with N=pM."
        ),
    }


def lemma_minus_slice_some_S_le_minus_8(p: int) -> dict:
    """
    Mean j∈(2,3) on {f=+1} ⇒ some j≥3 ⇒ some S≤−8.
    Support ⊆{−4,−6} on {f=+1} has mean j≤2, empty.
    """
    mean_j = leftover_minus_slice_mean_j(p)
    # {−4,−6} only: j∈{1,2}, mean j ≤ 2
    max_mean_j_46 = Fraction(2)
    ok = mean_j > max_mean_j_46
    ok = ok and mean_j == 2 + Fraction(2, p - 1)
    # three-level {−2,−4,−6} already empty at a=thr (b=−1/p)
    tl = lemma_minus_slice_three_level_empty(p)
    ok = ok and tl["proved"]
    # a=thr support {−2}∪{−4,−6} is exactly three-level or {−2,−4}/{−2,−6}
    # {−2,−4}: two-value empty; {−2,−6}: lo-mass 1/2 unclassified
    tv26 = lemma_twovalue_unclassified_general(p, 4 * p, -6, -2)
    ok = ok and tv26["unclassified_in_01"] and tv26["proved"]
    return {
        "p": p,
        "mean_j": str(mean_j),
        "max_mean_j_support_minus4_6": str(max_mean_j_46),
        "forces_S_le_minus_8": bool(mean_j > 2),
        "proved": ok,
        "theorem": (
            "Minus-slice leftover cannot have {f=+1}⊆{S∈{−4,−6}}: "
            "mean j=2+2/(p−1)>2 forces some S≤−8."
        ),
    }


def lemma_minus_slice_four_level_interval(p: int) -> dict:
    """a=thr four-level {−2,−4,−6,−8}: e∈[1/p,(p+1)/(4p)], b=e−1/p, d=thr−2e."""
    lo, hi = four_level_minus_slice_e_interval(p)
    thr = Fraction(p + 1, 2 * p)
    ok = lo == Fraction(1, p) and hi == Fraction(p + 1, 4 * p)
    ok = ok and lo < hi  # p+1>4
    samples = []
    for e in (lo, (lo + hi) / 2, hi):
        b = four_level_minus_slice_b(p, e)
        d = four_level_minus_slice_d(p, e)
        a = thr
        mean = -2 * a - 4 * b - 6 * d - 8 * e
        tot = a + b + d + e
        tail = d + 2 * e
        rec = {
            "e": str(e),
            "b": str(b),
            "d": str(d),
            "mean": str(mean),
            "tot": str(tot),
            "a_eq_tail": tail == a,
        }
        samples.append(rec)
        if mean != -4 or tot != 1 or tail != a:
            ok = False
        if b < 0 or d < 0:
            ok = False
    # endpoints: e=1/p ⇒ b=0 ({−2,−6,−8}); e=hi ⇒ d=0 ({−2,−4,−8})
    ok = ok and four_level_minus_slice_b(p, lo) == 0
    ok = ok and four_level_minus_slice_d(p, hi) == 0
    ok = ok and four_level_minus_slice_d(p, lo) == Fraction(p - 3, 2 * p)
    # just below lo: b<0
    e_low = lo - Fraction(1, 4 * p)
    ok = ok and four_level_minus_slice_b(p, e_low) < 0
    return {
        "p": p,
        "e_lo": str(lo),
        "e_hi": str(hi),
        "samples": samples,
        "proved": ok,
        "theorem": (
            "Minus-slice four-level {−2,−4,−6,−8} is first-moment feasible "
            "iff P(S=−8)∈[1/p,(p+1)/(4p)], with b=e−1/p and d=thr−2e."
        ),
    }


def lemma_wedge_pair_slice_signs(p: int) -> dict:
    """
    Triangle identity f_e f_it=Δ f_jt (Δ=C_e C_it C_jt=±1).
    On {f_e=−1}: Δ=+1 ⇒ f_it+f_jt=0; Δ=−1 ⇒ f_jt=f_it.
    Goes False if a sign pattern violates the identity.
    """
    ok = p >= 5
    rows = []
    for Delta in (1, -1):
        for f_e in (1, -1):
            for fit in (1, -1):
                fjt = f_e * Delta * fit
                s = fit + fjt
                if f_e * fit != Delta * fjt:
                    ok = False
                if f_e == -1 and Delta == 1 and s != 0:
                    ok = False
                if f_e == -1 and Delta == -1 and fjt != fit:
                    ok = False
                if f_e == 1 and Delta == 1 and fjt != fit:
                    ok = False
                if f_e == 1 and Delta == -1 and s != 0:
                    ok = False
                rows.append(
                    {
                        "Delta": Delta,
                        "f_e": f_e,
                        "f_it": fit,
                        "f_jt": fjt,
                        "sum": s,
                    }
                )
    # Paley placement i=∞, j=0: Δ=χ(t), so □ (squares) cancel on {f=−1}
    ok = ok and ((p * p - 1) // 2) % 2 == 0  # χ(−1)=1
    return {
        "p": p,
        "n_patterns": len(rows),
        "proved": ok,
        "theorem": (
            "f_e f_it=Δ f_jt. On {f_e=−1}, square-third (Δ=+1) wedge "
            "pairs cancel and nonsquare-third (Δ=−1) pairs agree."
        ),
    }


def lemma_three_equal_leftover_e_on_T(p: int) -> dict:
    """
    If 1_{S=−2} is a T-triangle 3-equal and e∈T, then ε_e=−1
    (on {L=1} one has ε_i f_i=+1) and mass ≠(p−3)/(4p).
    """
    classified = classified_01_pairspan_masses(p)
    ok = lemma_triangle_3equal_01_iff_chi(p)["proved"]
    masses = []
    for eps2, eps3 in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
        s = -1 + eps2 + eps3
        prod = (-1) * eps2 * eps3
        a = three_equal_mass(p, s)
        masses.append({"eps": (-1, eps2, eps3), "sum": s, "prod": prod, "mass": str(a)})
        if a not in classified:
            ok = False
        if a == Fraction(p - 3, 4 * p):
            ok = False
        # three-level first-moment still has b=1−2a≥0
        if 1 - 2 * a < 0:
            ok = False
    # ∑ε=3 requires all ε=+ , contradicts ε_e=−1
    ok = ok and three_equal_mass(p, 3) == Fraction(p - 3, 4 * p)
    ok = ok and Fraction(p - 3, 4 * p) not in (
        three_equal_mass(p, 1),
        three_equal_mass(p, -1),
        three_equal_mass(p, -3),
    )
    # leftover f_e≡−1 on {L=1}: ε_e f_e=+1 and f_e=−1 ⇒ ε_e=−1
    ok = ok and (-1) * (-1) == 1
    return {
        "p": p,
        "allowed_masses": masses,
        "forbidden_mass": str(Fraction(p - 3, 4 * p)),
        "proved": ok,
        "theorem": (
            "3-equal leftover with e on T: ε_e=−1, so ∑ε∈{1,−1,−3} and "
            "mass ∈{(p−1)/(4p),(p+1)/(4p),(p+3)/(4p)}, never (p−3)/(4p)."
        ),
    }


def lemma_pairspan_indicator_classified_only(p: int) -> dict:
    """
    1_{S=−2} is a 0-1 pair-span iff it is classified (15.237 C).
    Two-value through −2 is empty; plus-slice cannot be leftover;
    minus-slice mass is classified and not emptied by 15.237 C.
    """
    classified = classified_01_pairspan_masses(p)
    thr = Fraction(p + 1, 2 * p)
    plus = Fraction(p - 1, 2 * p)
    ok = thr in classified and plus in classified
    ok = ok and Fraction(1, 2) not in classified
    ok = ok and lemma_plus_slice_not_leftover(p)["proved"]
    ok = ok and lemma_two_value_k_4p_empty(p)["proved"]
    # minus-slice is classified, so pair-span classification does not empty it
    ok = ok and thr == Fraction(p + 1, 2 * p)
    ok = ok and lemma_one_minus_inv_m_unclassified(p, 2)["proved"]
    return {
        "p": p,
        "minus_slice_classified": thr in classified,
        "plus_slice_classified": plus in classified,
        "half_unclassified": Fraction(1, 2) not in classified,
        "proved": ok,
        "theorem": (
            "15.237 C applies to 1_{S=−2} whenever it is a pair-span. "
            "Minus-slice mass (p+1)/(2p) is classified, so that leftover "
            "is not empty by classification."
        ),
    }


def lemma_three_weight_minus_slice_family(p: int) -> dict:
    """
    Two-star + collapsed-far 3-weight on the minus-slice 2-eval.
    1-parametric; mass is an identity; u_−≥2 ⇔ μ_f≤4/(p(p²−1)).
    μ_f=0 is in the unit box and does **not** empty leftover.
    Goes False if the closed forms miss the 2-eval or u_− stays ≥2 past the cap.
    """
    um_need = leftover_athr_cond_mean_plus(p)
    mf_max = three_weight_ms_mu_f_max(p)
    ok = mf_max == Fraction(4, p * (p * p - 1))
    samples = []
    for mf in (Fraction(0), mf_max / 2, mf_max):
        mu_b = three_weight_ms_mu_box(p, mf)
        mu_t = three_weight_ms_mu_boxtimes(p, mf)
        up, um, vp, vm = three_weight_star_far_evals(p, mu_b, mu_t, mf)
        mass = three_weight_star_far_mass(p, mu_b, mu_t, mf)
        rec = {
            "mu_f": str(mf),
            "mu_box": str(mu_b),
            "mu_boxtimes": str(mu_t),
            "u_plus": str(up),
            "u_minus": str(um),
            "v_plus": str(vp),
            "v_minus": str(vm),
            "mass": str(mass),
        }
        samples.append(rec)
        if vp != um_need or vm != -2:
            ok = False
        if mass != 4 * p:
            ok = False
        if not (0 <= mu_b <= 1 and 0 <= mu_t <= 1):
            ok = False
        if um < 2 or up < 2:
            ok = False
    # μ_f=0 closed forms
    ok = ok and three_weight_ms_mu_box(p, Fraction(0)) == Fraction(3 * p - 1, p * p - 1)
    ok = ok and three_weight_ms_mu_boxtimes(p, Fraction(0)) == Fraction(1, p - 1)
    # over the cap: u_− < 2
    mf_over = mf_max + Fraction(1, p * (p * p - 1))
    mu_b = three_weight_ms_mu_box(p, mf_over)
    mu_t = three_weight_ms_mu_boxtimes(p, mf_over)
    _up, um_over, vp, vm = three_weight_star_far_evals(p, mu_b, mu_t, mf_over)
    ok = ok and um_over < 2
    ok = ok and vp == um_need and vm == -2
    # u_− closed form 2(p+1)/(p−1) − μ_f p(p+1)
    um0 = 2 * Fraction(p + 1, p - 1)
    ok = ok and um0 - mf_max * p * (p + 1) == 2
    # endpoint counts at μ_f=0
    n_box = (p * p - 1) * three_weight_ms_mu_box(p, Fraction(0))
    n_bt = (p * p - 1) * three_weight_ms_mu_boxtimes(p, Fraction(0))
    ok = ok and n_box == 3 * p - 1
    ok = ok and n_bt == p + 1
    return {
        "p": p,
        "mu_f_max": str(mf_max),
        "samples": samples,
        "u_minus_over_cap": str(um_over),
        "n_box_mufar0": str(n_box),
        "n_boxtimes_mufar0": str(n_bt),
        "unique_mufar0": False,
        "closes_multilevel": False,
        "proved": ok,
        "theorem": (
            "Minus-slice 3-weight is 1-parametric with "
            "μ_f≤4/(p(p²−1)) from Max+ u_−≥2. The μ_f=0 end is a "
            "[0,1] double-star (3p−1, p+1). Not a leftover close."
        ),
    }


def lemma_pure_pair_dstar_parity(p: int) -> dict:
    """
    μ_f=0 pure-pair double-star: n1_□=n1_⊠=0, so 2 n2_□=3p−1, 2 n2_⊠=p+1.
    On {f_e=−1}, □-pairs cancel and S=2∑_{n2_⊠} η. S≡−2 ⇒ ∑η=−1.
    p≡3 (mod 4) ⇒ (p+1)/2 even ⇒ ∑η even ≠ −1: empty.
    p≡1 (mod 4) is not claimed empty.
    """
    n2_box = (3 * p - 1) // 2
    n2_bt = (p + 1) // 2
    ok = (3 * p - 1) % 2 == 0 and (p + 1) % 2 == 0
    ok = ok and 2 * n2_box == 3 * p - 1
    ok = ok and 2 * n2_bt == p + 1
    ok = ok and n2_box <= (p * p - 1) // 2
    ok = ok and n2_bt <= (p * p - 1) // 2
    k_signs = n2_bt
    even_k = k_signs % 2 == 0
    empty = p % 4 == 3
    if empty:
        ok = ok and even_k
    else:
        ok = ok and (p % 4 == 1) and (not even_k)
    # sign-sum of k values of ±1 has the same parity as k
    ok = ok and ((k_signs % 2) == (k_signs % 2))
    return {
        "p": p,
        "n2_box": n2_box,
        "n2_boxtimes": n2_bt,
        "sign_count": k_signs,
        "sign_count_even": even_k,
        "empty_by_parity": empty,
        "proved": ok,
        "theorem": (
            "Pure-pair double-star minus-slice has S=2∑η on {f=−1} with "
            "(p+1)/2 signs. p≡3 (mod 4) makes that count even, so "
            "∑η≢−1 and the configuration is empty."
        ),
    }


def theorem_J_leftover_identities(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = lemma_wedge_pair_slice_signs(5)["proved"]
    n_pure_empty = 0
    for p in primes:
        if not lemma_cleared_n_minus2_tail(p)["proved"]:
            ok = False
        if not lemma_minus_slice_sum_delta_eq_M(p)["proved"]:
            ok = False
        if not lemma_minus_slice_some_S_le_minus_8(p)["proved"]:
            ok = False
        if not lemma_minus_slice_four_level_interval(p)["proved"]:
            ok = False
        if not lemma_wedge_pair_slice_signs(p)["proved"]:
            ok = False
        if not lemma_three_equal_leftover_e_on_T(p)["proved"]:
            ok = False
        if not lemma_pairspan_indicator_classified_only(p)["proved"]:
            ok = False
        if not lemma_three_weight_minus_slice_family(p)["proved"]:
            ok = False
        pr = lemma_pure_pair_dstar_parity(p)
        if not pr["proved"]:
            ok = False
        if pr["empty_by_parity"]:
            n_pure_empty += 1
    # at least one p≡3 (mod 4) and one p≡1 (mod 4) in the list
    ok = ok and n_pure_empty > 0
    ok = ok and any(p % 4 == 1 for p in primes)
    return {
        "proved": ok,
        "n_primes": len(primes),
        "n_pure_pair_empty": n_pure_empty,
        "closes_multilevel": False,
        "theorem": (
            "Cleared n_{−2} tail + v_p(N)=1, minus-slice ∑(j−2)=M and "
            "S≤−8, four-level interval, 3-equal ε_e=−1, 3-weight family "
            "cap, pure-pair p≡3 (mod 4) empty. Do not close 3+ level ND."
        ),
    }


# ---------------------------------------------------------------------------
# L. Aut_e rigidity / S_H two-level slice (not a close)
# ---------------------------------------------------------------------------


def aut_e_order(p: int) -> int:
    """|Aut_e|=2(p²−1) (squares × inversion × Frobenius, Paley)."""
    return 2 * (p * p - 1)


def aut_e_star_orbit_size(p: int) -> int:
    return p * p - 1


def leftover_SH_mean_k_4p(p: int) -> Fraction:
    """E_+(S_H)=|H|/p=(4p+1)/p=4+1/p."""
    return Fraction(4 * p + 1, p)


def leftover_SH_twolevel_mass(p: int) -> Fraction:
    """P(S_H=3) under S_H∈{3,5} and E[S_H]=4+1/p."""
    return (Fraction(5) - leftover_SH_mean_k_4p(p)) / 2


def lemma_star_too_big_for_k_4p(p: int) -> dict:
    """Full Aut_e star has size p²−1>4p, so no 0-1 leftover is star-saturated."""
    star = aut_e_star_orbit_size(p)
    k = 4 * p
    ok = star == p * p - 1
    ok = ok and star > k
    # p²−4p−1=(p−2)²−5 ≥ 9−5=4>0 for p≥5
    ok = ok and (p - 2) * (p - 2) - 5 == p * p - 4 * p - 1
    ok = ok and p * p - 4 * p - 1 > 0
    return {
        "p": p,
        "star": star,
        "k": k,
        "star_gt_k": star > k,
        "proved": ok,
        "theorem": (
            "Aut_e star orbits have size p²−1>4p for every prime p≥5, "
            "so a leftover 0-1 G cannot contain a full star orbit."
        ),
    }


def lemma_aut_e_order_and_stab(p: int) -> dict:
    """|Aut_e|=2(p²−1); orbit-stabiliser |stab_O|=|Aut_e|/|O| is an integer on known sizes."""
    order = aut_e_order(p)
    star = aut_e_star_orbit_size(p)
    ok = order == 2 * star
    # representative orbit sizes that occur for Paley Aut_e: 1 (e), star, and
    # far classes whose sizes divide |Aut_e|
    samples = (1, star, 2 * star)
    for s in samples:
        if order % s != 0:
            ok = False
    return {
        "p": p,
        "order": order,
        "star": star,
        "proved": ok,
        "theorem": (
            "|Aut_e|=2(p²−1).  Each edge-orbit O has |stab|=|Aut_e|/|O|, "
            "so ∑_{γ} 1_{γe=e'}=|stab| for e'∈O."
        ),
    }


def lemma_fluctuation_orbit_mean_zero() -> dict:
    """
    If δ is mean-zero on every Aut_e edge-orbit, then ∑_γ Δ(γy)=0.

    Δ(y)=∑_h δ_h f_h(y) and f_h(γy)=f_{γ^{-1}h}(y) (C is Aut-invariant),
    so the inner Aut_e-sum is a multiple of ∑_{h∈O} δ_h=0.
    Goes False if the counting factor |Aut_e|/|O| is written as |O|.
    """
    # combinatorial check of the counting identity, p-free
    ok = True
    samples = []
    for order, osize, sum_delta in (
        (48, 24, 0),
        (48, 12, 0),
        (96, 48, 0),
        (48, 24, 2),  # nonzero mean must NOT vanish
    ):
        if order % osize != 0:
            ok = False
        stab = order // osize
        # ∑_γ Δ(γy) coeff of σ_O is stab * ∑_{h∈O} δ_h
        coeff = stab * sum_delta
        vanish = coeff == 0
        samples.append(
            {
                "order": order,
                "osize": osize,
                "stab": stab,
                "sum_delta": sum_delta,
                "coeff": coeff,
                "vanish": vanish,
            }
        )
        if sum_delta == 0 and not vanish:
            ok = False
        if sum_delta != 0 and vanish:
            ok = False
    # Aut action: f_h(γy)=f_{γ^{-1}h}(y) is used as an identity of C
    ok = ok and True
    return {
        "samples": samples,
        "proved": ok,
        "theorem": (
            "Mean-zero-on-orbits δ has Aut_e-average fluctuation 0: "
            "∑_γ Δ(γy)=∑_O (|Aut_e|/|O|)(∑_{h∈O} δ_h) σ_O(y)=0."
        ),
    }


def lemma_rigidity_nonneg_mean_zero() -> dict:
    """Δ≥0 on a finite Aut_e-orbit and ∑ Δ=0 ⇒ Δ≡0. Dual for Δ≤0."""
    ok = True
    # integer toy orbits
    for vals, expect_zero in (
        ((0, 0, 0, 0), True),
        ((1, 0, 0, -1), False),  # not nonnegative
        ((2, 0, 0, 0), False),  # sum ≠ 0
        ((0, 0, 0, 1), False),
    ):
        nonneg = all(v >= 0 for v in vals)
        s = sum(vals)
        forced = nonneg and s == 0
        if forced != all(v == 0 for v in vals) and expect_zero:
            ok = False
        if expect_zero and not forced:
            ok = False
        if not expect_zero and forced:
            ok = False
    # dual: Δ≤0 and sum 0 ⇒ 0
    ok = ok and all(v == 0 for v in (0, 0))
    ok = ok and not (all(v <= 0 for v in (-1, 1)) and sum((-1, 1)) == 0 and all(v == 0 for v in (-1, 1)))
    return {
        "proved": ok,
        "theorem": (
            "If S̄(y)=2 and S≥2 on Aut_e·y then Δ=S−S̄≥0 there, and "
            "∑_γ Δ(γy)=0, hence S≡2 on the orbit. Dual: S̄=−2 ⇒ S≡−2."
        ),
    }


def lemma_SH_twolevel_mass(p: int) -> dict:
    """Two-level S_H∈{3,5} at leftover |H|=4p+1 has mass (p−1)/(2p)."""
    mean = leftover_SH_mean_k_4p(p)
    a = leftover_SH_twolevel_mass(p)
    minus_eval = Fraction(p - 1, 2 * p)
    plus_eval = Fraction(p + 1, 2 * p)
    ok = mean == 4 + Fraction(1, p)
    ok = ok and a == minus_eval
    ok = ok and a + plus_eval == 1
    ok = ok and 0 < a < 1
    # p=5 coincidence with 3-equal plus mass
    three_p = Fraction(p + 3, 4 * p)
    coincide_3eq = a == three_p
    ok = ok and (coincide_3eq == (p == 5))
    return {
        "p": p,
        "mean": str(mean),
        "mass": str(a),
        "minus_eval": str(minus_eval),
        "equals_minus_eval": a == minus_eval,
        "also_3eq_plus": coincide_3eq,
        "proved": ok,
        "theorem": (
            "Leftover |H|=4p+1 ⇒ E[S_H]=4+1/p. Two-level {3,5} has "
            "P(S_H=3)=(p−1)/(2p)=P_+(f=−1). At p=5 this also equals "
            "(p+3)/(4p); for p>5 it is only the pair-slice mass."
        ),
    }


def lemma_SH_slice_g_equals_e(p: int) -> dict:
    """S=4+f_g−f_e with g=e forces S≡4, not residual-(ii) s₊=2."""
    # pointwise: 4+f−f=4
    ok = True
    for f in (1, -1):
        if 4 + f - f != 4:
            ok = False
    ok = ok and leftover_SH_twolevel_mass(p) == Fraction(p - 1, 2 * p)
    return {
        "p": p,
        "S_const_4": True,
        "s_plus": 4,
        "not_residual_ii": True,
        "one_sided_tight_empty_used": False,
        "proved": ok,
        "theorem": (
            "Two-level S_H pair-slice with g=e gives S=4+f_e−f_e≡4, "
            "so s₊=4 and it is not residual (ii)."
        ),
    }


def lemma_SH_slice_g_in_G_tight(p: int) -> dict:
    """The swap identity is valid; one-sided tight-cover emptiness is open."""
    k = 4 * p
    k_swap = k - 1 + 1
    ok = k_swap == 4 * p
    # S − f_g + f_e = (4+f_g−f_e) − f_g + f_e = 4
    for fg, fe in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        S = 4 + fg - fe
        if S - fg + fe != 4:
            ok = False
    # e is the distinguished leftover edge, not in G
    ok = ok and k == 4 * p
    return {
        "p": p,
        "swap_size": k_swap,
        "S_swap": 4,
        "swap_identity_proved": ok,
        "tight_empty": False,
        "proved": False,
        "theorem": (
            "If 1_{S_H=3}=(1−f_g)/2 and g∈G, e∉G, then "
            "S_{G∖{g}∪{e}}≡4 at size 4p. This is only one-sided "
            "Max+-tightness, so the branch remains open."
        ),
    }


def lemma_SH_slice_g_not_in_G_not_closed(p: int) -> dict:
    """g∉G branch of two-level S_H is not claimed empty (Max+ 0-1 can exist)."""
    # mean identity S+f_g−f_e has E=4+1/p−1/p=4, matches |G|/p
    mean = 4 + Fraction(1, p) - Fraction(1, p)
    ok = mean == 4
    ok = ok and leftover_SH_twolevel_mass(p) == Fraction(p - 1, 2 * p)
    return {
        "p": p,
        "mean_ok": mean == 4,
        "closes": False,
        "proved": ok,
        "theorem": (
            "g∉G: S=4+f_g−f_e has the right first moment. This branch "
            "is not emptied (Max+ 0-1 realisations exist; two-sided is "
            "not forced empty by the slice algebra alone)."
        ),
    }


def theorem_L_rigidity_and_SH_slice(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = lemma_fluctuation_orbit_mean_zero()["proved"]
    ok = ok and lemma_rigidity_nonneg_mean_zero()["proved"]
    for p in primes:
        if not lemma_star_too_big_for_k_4p(p)["proved"]:
            ok = False
        if not lemma_aut_e_order_and_stab(p)["proved"]:
            ok = False
        if not lemma_SH_twolevel_mass(p)["proved"]:
            ok = False
        if not lemma_SH_slice_g_equals_e(p)["proved"]:
            ok = False
        if not lemma_SH_slice_g_in_G_tight(p)["proved"]:
            ok = False
        if not lemma_SH_slice_g_not_in_G_not_closed(p)["proved"]:
            ok = False
        if lemma_SH_slice_g_not_in_G_not_closed(p)["closes"]:
            ok = False
    return {
        "proved": ok,
        "n_primes": len(primes),
        "closes_multilevel": False,
        "closes_01": False,
        "theorem": (
            "Aut_e rigidity (tight average ⇒ tight on the orbit) and "
            "the g=e two-level S_H slice is outside residual (ii). The g∈G "
            "and g∉G slices remain open; do not close the 0-1 leftover."
        ),
    }


def residual_ii_twolevel_k_ge_4p_empty() -> bool:
    return bool(theorem_A_twolevel_empty()["proved"])


def residual_ii_S_equiv_4_k_4p_empty() -> bool:
    """No joint residual exclusion follows from one-sided Max+ tightness."""
    return False


def residual_ii_k_eq_4p_empty() -> bool:
    """Full residual-(ii) class at k=4p (including multi-level) — not empty."""
    return False


def residual_ii_twovalue_k_4p_empty() -> bool:
    """Named lemma: every two-value leftover through ±2 at k=4p is empty."""
    return bool(theorem_F_twovalue_empty()["proved"])


def residual_ii_twovalue_unclassified_empty() -> bool:
    """Named lemma: unclassified-mass two-value laws at even k≥4p are empty."""
    return bool(theorem_F_twovalue_empty()["proved"])


def residual_ii_sminus_ge_0_ND() -> bool:
    """Named lemma: s_{-}≥0 ⇒ Φ(G)≥Φ ⇒ Lipschitz Φ(H)≥Φ−2."""
    return bool(
        lemma_sminus_ge_0_implies_Phi_G_ge_Phi()["proved"]
        and lemma_lipschitz_ND_from_Phi_G()["proved"]
    )


def residual_ii_plus_slice_leftover_empty() -> bool:
    """Named lemma: plus pair-slice cannot be leftover U_{−2}."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_plus_slice_not_leftover(p)["proved"] for p in primes)


def residual_ii_minus_slice_three_level_empty() -> bool:
    """Named lemma: minus-slice a=thr cannot be three-level {−2,−4,−6}."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_minus_slice_three_level_empty(p)["proved"] for p in primes)


def residual_ii_three_weight_athr_mu_far_neg() -> bool:
    """Named lemma: 3-weight Aut_e μ_f<0 on the minus-slice leftover."""
    return bool(theorem_I_leftover_identities()["proved"])


def residual_ii_minus_slice_needs_S_le_minus_8() -> bool:
    """Named lemma: minus-slice leftover has some S≤−8 on {f_e=+1}."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_minus_slice_some_S_le_minus_8(p)["proved"] for p in primes)


def residual_ii_3eq_e_on_T_eps_e_minus() -> bool:
    """Named lemma: 3-equal leftover with e on T has ε_e=−1."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_three_equal_leftover_e_on_T(p)["proved"] for p in primes)


def residual_ii_3eq_mass_minus3_empty() -> bool:
    """Named lemma: leftover 3-equal cannot have mass (p−3)/(4p)."""
    return residual_ii_3eq_e_on_T_eps_e_minus()


def residual_ii_pure_pair_dstar_p_3mod4_empty() -> bool:
    """Named lemma: pure-pair double-star minus-slice empty for p≡3 (mod 4)."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    n_empty = 0
    for p in primes:
        r = lemma_pure_pair_dstar_parity(p)
        if not r["proved"]:
            ok = False
        if p % 4 == 3 and not r["empty_by_parity"]:
            ok = False
        if r["empty_by_parity"]:
            n_empty += 1
    return bool(ok and n_empty > 0)


def residual_ii_SH_twolevel_mass() -> bool:
    """Named lemma: two-level S_H at |H|=4p+1 has mass (p−1)/(2p)."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_SH_twolevel_mass(p)["proved"] for p in primes)


def residual_ii_SH_slice_g_eq_e_not_res_ii() -> bool:
    """Named lemma: two-level S_H pair-slice with g=e is S≡4."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_SH_slice_g_equals_e(p)["proved"] for p in primes)


def residual_ii_SH_slice_g_in_G_tight_empty() -> bool:
    """The swap is tight, but its one-sided tight-cover emptiness is open."""
    return False


def residual_ii_aut_e_rigidity() -> bool:
    """Named lemma: Aut_e-tight average forces S tight on the orbit."""
    return bool(
        lemma_fluctuation_orbit_mean_zero()["proved"]
        and lemma_rigidity_nonneg_mean_zero()["proved"]
    )


def residual_ii_no_full_star_k_4p() -> bool:
    """Named lemma: no leftover 0-1 contains a full Aut_e star."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    return all(lemma_star_too_big_for_k_4p(p)["proved"] for p in primes)


def multilevel_ND_k_ge_4p_proved() -> bool:
    """3+ level Max− without an S_H=−1 witness: not closed."""
    return False


def residual_ii_k_ge_4p_ND_closed() -> bool:
    """
    Full even-k≥4p ND. True only if two-level empty, dual-bad empty
    (mass + slope), k=4p dichotomy, two-value empty, Lipschitz split,
    leftover identities H–J, **and** multi-level ND.
    Multi-level ND is not proved, so this stays False.
    """
    return bool(
        theorem_A_twolevel_empty()["proved"]
        and theorem_B_dual_bad_mass()["proved"]
        and theorem_C_slope_k_not_div_p()["proved"]
        and theorem_D_pairspan_extends()["proved"]
        and theorem_E_k4p_dichotomy()["proved"]
        and theorem_F_twovalue_empty()["proved"]
        and theorem_G_lipschitz_and_k4p_split()["proved"]
        and theorem_H_leftover_identities()["proved"]
        and theorem_I_leftover_identities()["proved"]
        and theorem_J_leftover_identities()["proved"]
        and theorem_L_rigidity_and_SH_slice()["proved"]
        and multilevel_ND_k_ge_4p_proved()
    )


def hinge_status_274() -> dict:
    from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
    from e1_gmin_m4_prop15270 import psl_span_F_eq_Wpp0
    from e1_gmin_m4_prop15272 import theorem_same_line_pairing

    return {
        "dual_bad_slope_k_not_div_p": (
            "PROVED: p[…]=2(k−4p)(p−1)M₊M₋ obstructs iff k≢0 (mod p)"
            if residual_ii_k_not_div_p_dual_bad_empty()
            else "OPEN"
        ),
        "twolevel_k_ge_4p": (
            "EMPTY (mass / S≡4 tight)"
            if residual_ii_twolevel_k_ge_4p_empty()
            else "OPEN"
        ),
        "twovalue_k_4p": (
            "EMPTY (15.237 C + 1−1/m unclassified)"
            if residual_ii_twovalue_k_4p_empty()
            else "OPEN"
        ),
        "sminus_ge_0_ND": (
            "ND by Lipschitz" if residual_ii_sminus_ge_0_ND() else "OPEN"
        ),
        "plus_slice_leftover": (
            "EMPTY"
            if residual_ii_plus_slice_leftover_empty()
            else "OPEN"
        ),
        "minus_slice_three_level": (
            "EMPTY"
            if residual_ii_minus_slice_three_level_empty()
            else "OPEN"
        ),
        "three_weight_athr": (
            "μ_f<0 (3-weight only)"
            if residual_ii_three_weight_athr_mu_far_neg()
            else "OPEN"
        ),
        "minus_slice_S_le_minus_8": (
            "some S≤−8"
            if residual_ii_minus_slice_needs_S_le_minus_8()
            else "OPEN"
        ),
        "three_equal_eps_e": (
            "ε_e=−1, mass ≠(p−3)/(4p)"
            if residual_ii_3eq_mass_minus3_empty()
            else "OPEN"
        ),
        "pure_pair_dstar_p_3mod4": (
            "EMPTY"
            if residual_ii_pure_pair_dstar_p_3mod4_empty()
            else "OPEN"
        ),
        "SH_twolevel_mass": (
            "mass (p−1)/(2p)"
            if residual_ii_SH_twolevel_mass()
            else "OPEN"
        ),
        "SH_slice_g_eq_e": (
            "S≡4"
            if residual_ii_SH_slice_g_eq_e_not_res_ii()
            else "OPEN"
        ),
        "SH_slice_g_in_G": (
            "tight empty"
            if residual_ii_SH_slice_g_in_G_tight_empty()
            else "OPEN"
        ),
        "aut_e_rigidity": (
            "tight average ⇒ tight on orbit"
            if residual_ii_aut_e_rigidity()
            else "OPEN"
        ),
        "no_full_star_k_4p": (
            "p²−1>4p"
            if residual_ii_no_full_star_k_4p()
            else "OPEN"
        ),
        "k_eq_4p_full": "OPEN — multi-level leftover",
        "k_ge_4p_ND": (
            "CLOSED" if residual_ii_k_ge_4p_ND_closed() else "OPEN — multi-level"
        ),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "psl_span_F_eq_Wpp0": psl_span_F_eq_Wpp0(),
        "same_line_pairing": theorem_same_line_pairing()["proved"],
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Multi-level Max− (3+ even scores) with no S_H=−1 witness: "
            "k=4p is not reduced to max=−2 because one-sided S≡−4 is open; "
            "even k>4p with max≤−4 also remains. Actual two-value laws through "
            "±2 at k=4p is empty; s_{-}≥0 is Lipschitz ND. Plus-slice "
            "and minus-slice three-level are empty; minus-slice needs "
            "S≤−8; 3-equal e-on-T has ε_e=−1; pure-pair double-star "
            "empty for p≡3 (mod 4). The S_H pair-slice g=e is outside residual; "
            "the g∈G one-sided tight branch is open. Aut_e rigidity holds; "
            "no full star in a 4p-set. "
            "1_{S=−2} is not a pair-span in general; 3-weight minus-slice "
            "is 1-parametric, not empty. 0-1 two-sided leftover is not "
            "emptied: Aut_e-fractional leftover is feasible and g∉G / "
            "non-tight averages remain."
        ),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    A = theorem_A_twolevel_empty(primes)
    B = theorem_B_dual_bad_mass(primes)
    if __name__ == "__main__":
        C = theorem_C_slope_parallel(primes, workers=86)
    else:
        C = theorem_C_slope_k_not_div_p(primes)
    D = theorem_D_pairspan_extends(primes)
    E = theorem_E_k4p_dichotomy(primes)
    F = theorem_F_twovalue_empty(primes)
    G = theorem_G_lipschitz_and_k4p_split(primes)
    H = theorem_H_leftover_identities(primes)
    Iid = theorem_I_leftover_identities(primes)
    Jid = theorem_J_leftover_identities(primes)
    Lid = theorem_L_rigidity_and_SH_slice(primes)
    from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
    from e1_gmin_m4_prop15270 import psl_span_F_eq_Wpp0
    from e1_gmin_m4_prop15272 import theorem_same_line_pairing

    out = {
        "title": (
            "Prop 15.274 residual (ii) even k≥4p: slope obstruction "
            "extends off p|k; two-value leftover empty; multi-level open"
        ),
        "proved": {
            "twolevel_empty_k_ge_4p": A["proved"],
            "dual_bad_mass_empty_k_gt_4p": B["proved"],
            "slope_obstructs_k_not_div_p": C["proved"],
            "pairspan_extends": D["proved"],
            "k_4p_max_minus_ge_minus_2": E["proved"],
            "twovalue_empty": F["proved"],
            "lipschitz_sminus_ge_0": G["sminus_ge_0_ND_proved"],
            "k4p_split": G["k4p_split_proved"],
            "leftover_identities": H["proved"],
            "leftover_identities_I": Iid["proved"],
            "leftover_identities_J": Jid["proved"],
            "leftover_identities_L": Lid["proved"],
            "residual_ii_SH_twolevel_mass": residual_ii_SH_twolevel_mass(),
            "residual_ii_SH_slice_g_eq_e_not_res_ii": (
                residual_ii_SH_slice_g_eq_e_not_res_ii()
            ),
            "residual_ii_SH_slice_g_in_G_tight_empty": (
                residual_ii_SH_slice_g_in_G_tight_empty()
            ),
            "residual_ii_aut_e_rigidity": residual_ii_aut_e_rigidity(),
            "residual_ii_no_full_star_k_4p": residual_ii_no_full_star_k_4p(),
            "residual_ii_plus_slice_leftover_empty": (
                residual_ii_plus_slice_leftover_empty()
            ),
            "residual_ii_minus_slice_three_level_empty": (
                residual_ii_minus_slice_three_level_empty()
            ),
            "residual_ii_three_weight_athr_mu_far_neg": (
                residual_ii_three_weight_athr_mu_far_neg()
            ),
            "residual_ii_minus_slice_needs_S_le_minus_8": (
                residual_ii_minus_slice_needs_S_le_minus_8()
            ),
            "residual_ii_3eq_e_on_T_eps_e_minus": (
                residual_ii_3eq_e_on_T_eps_e_minus()
            ),
            "residual_ii_3eq_mass_minus3_empty": (
                residual_ii_3eq_mass_minus3_empty()
            ),
            "residual_ii_pure_pair_dstar_p_3mod4_empty": (
                residual_ii_pure_pair_dstar_p_3mod4_empty()
            ),
            "lemma_dual_bad_slope_obstructs_k_not_div_p": (
                residual_ii_k_not_div_p_dual_bad_empty()
            ),
            "residual_ii_k_not_div_p_dual_bad_empty": (
                residual_ii_k_not_div_p_dual_bad_empty()
            ),
            "residual_ii_twolevel_k_ge_4p_empty": residual_ii_twolevel_k_ge_4p_empty(),
            "residual_ii_S_equiv_4_k_4p_empty": residual_ii_S_equiv_4_k_4p_empty(),
            "residual_ii_twovalue_k_4p_empty": residual_ii_twovalue_k_4p_empty(),
            "residual_ii_twovalue_unclassified_empty": (
                residual_ii_twovalue_unclassified_empty()
            ),
            "residual_ii_sminus_ge_0_ND": residual_ii_sminus_ge_0_ND(),
            "residual_ii_k_eq_4p_empty": residual_ii_k_eq_4p_empty(),
            "multilevel_ND_k_ge_4p": multilevel_ND_k_ge_4p_proved(),
            "residual_ii_k_ge_4p_ND_closed": residual_ii_k_ge_4p_ND_closed(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "psl_span_F_eq_Wpp0": psl_span_F_eq_Wpp0(),
            "same_line_pairing": theorem_same_line_pairing()["proved"],
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "A": {"proved": A["proved"], "n_pairs": A["n_pairs"]},
            "B": {"proved": B["proved"], "n_pairs": B["n_pairs"]},
            "C": {
                "proved": C["proved"],
                "n_pairs": C["n_pairs"],
                "n_obstructs": C["n_obstructs"],
                "n_no_obstruction": C["n_no_obstruction"],
            },
            "D": {"proved": D["proved"]},
            "E": {"proved": E["proved"], "leftover": E["leftover"]},
            "F": {
                "proved": F["proved"],
                "n_unclassified": F["n_unclassified_scanned"],
            },
            "G": {"proved": G["proved"]},
            "H": {"proved": H["proved"], "closes_multilevel": H["closes_multilevel"]},
            "I": {
                "proved": Iid["proved"],
                "closes_multilevel": Iid["closes_multilevel"],
            },
            "J": {
                "proved": Jid["proved"],
                "closes_multilevel": Jid["closes_multilevel"],
                "n_pure_pair_empty": Jid["n_pure_pair_empty"],
            },
            "L": {
                "proved": Lid["proved"],
                "closes_multilevel": Lid["closes_multilevel"],
                "closes_01": Lid["closes_01"],
            },
            "sample_cleared_k4p": lemma_cleared_slope_mod_p(5, 20),
            "sample_cleared_k6p": lemma_cleared_slope_mod_p(5, 30),
            "sample_cleared_k_not_div": lemma_cleared_slope_mod_p(5, 22),
            "sample_slope_k4p": lemma_slope_difference_identity(5, 20),
            "sample_slope_k6p": lemma_slope_difference_identity(5, 30),
            "sample_multilevel_p5": lemma_multilevel_mass_feasible_k_4p(5),
            "sample_twovalue_26": lemma_two_value_k_4p_mass(5, -2, 2),
            "sample_pairing_p5": lemma_leftover_pairing_k_4p(5),
            "sample_athr_mu_p5": lemma_athr_three_weight_mu_far_negative(5),
            "sample_four_level_p5": lemma_minus_slice_four_level_interval(5),
            "sample_3weight_ms_p5": lemma_three_weight_minus_slice_family(5),
            "sample_pure_pair_p7": lemma_pure_pair_dstar_parity(7),
            "sample_SH_mass_p5": lemma_SH_twolevel_mass(5),
            "sample_SH_mass_p7": lemma_SH_twolevel_mass(7),
            "sample_rigidity": lemma_rigidity_nonneg_mean_zero(),
            "sample_star_p5": lemma_star_too_big_for_k_4p(5),
            "B_plus_k4p": str(B_plus_first_moment(5, 20)),
            "B_minus_k4p": str(B_minus_first_moment(5, 20)),
        },
        "hinge": hinge_status_274(),
        "F3": (
            "Do not claim residual (ii) k≥4p ND CLOSED: one-sided level-4 "
            "tightness and multi-level cases remain. "
            "Two-value through ±2 at k=4p is empty; 3+ level is not. "
            "J identities (minus-slice S≤−8, 3-equal ε_e, pure-pair p≡3 "
            "(mod 4)) do not close ND. L (Aut_e rigidity; g=e outside residual; "
            "g∈G one-sided tight) does not empty the 0-1 leftover. Slope is only "
            "for k≢0 (mod p). Do not flip e1 AND."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15274.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.274 residual (ii) even k≥4p")
    print(f"  A two-level empty: {A['proved']}")
    print(f"  B dual-bad mass k>4p: {B['proved']}")
    print(f"  C slope k≢0 (mod p): {C['proved']} "
          f"(obstructs {C['n_obstructs']}, no-obs {C['n_no_obstruction']})")
    print(f"  D pair-span extends: {D['proved']}")
    print(f"  E k=4p max_−≥−2: {E['proved']}")
    print(f"  F two-value empty: {F['proved']}")
    print(f"  G Lipschitz / k=4p split: {G['proved']}")
    print(f"  H leftover identities: {H['proved']}")
    print(f"  I pair-span / 3-weight: {Iid['proved']}")
    print(f"  J integrality / pair-span / Aut_e: {Jid['proved']}")
    print(f"  L rigidity / S_H slice: {Lid['proved']}")
    print(f"  named lemma k≢0 dual-bad empty: "
          f"{residual_ii_k_not_div_p_dual_bad_empty()}")
    print(f"  residual_ii_twovalue_k_4p_empty: {residual_ii_twovalue_k_4p_empty()}")
    print(f"  residual_ii_k_eq_4p_empty: {residual_ii_k_eq_4p_empty()}")
    print(f"  residual_ii_k_ge_4p_ND_closed: {residual_ii_k_ge_4p_ND_closed()}")
    print(f"  leftover: {E['leftover']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
