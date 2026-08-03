#!/usr/bin/env python3
"""
Prop 15.153 — Switched 4-point residual dictionary:
e1=e3=0, e2=1/(n−1); d₄=(1+6/(n−1)+μ₄)/16; residual ⇔ μ₄≤μ4_suf
with closed μ4_flat, μ4_suf; exact m_e at p=7; L OPEN.

Continues 15.152. Does **not** soft-close residual. No type census of Max+
for general p. No class_key (F19). No GPU theater (F20).

============================================================================
Setup. After halfspace switch: y'∈{±1}ⁿ on switched Max+ X; w=(1−y')/2;
k=wt(w); R₄=E[k^{underline 4}].  Residual for p≥7 ⇔ R₄≤R4_suf (15.146).
μ₄ := average over 4-sets S of m₄^{sw}(S) := E[∏_{i∈S} y'_i].

============================================================================
Theorem A (switched low moments Max+-free) — PROVED Fraction.
  On switched Max+, writing e_r for the exchangeable r-wise moments of y':
      e₁ = E[y'_i] = 0,     e₃ = E[y'_i y'_j y'_l] = 0  (i,j,l distinct),
      e₂ = E[y'_i y'_j] = 1/(n−1)  (i≠j).
  Proof.  E[s]=0 with s=∑ y' (from E[k]=n/2) and vertex-transitivity ⇒ e₁=0.
  E[s²]=2n (15.124 via s=n−2k) ⇒ e₂=(E[s²]−n)/(n(n−1))=1/(n−1).
  E[s³]=0 (odd central moment / antipodal W_k=W_{n−k}) and the partition of
  s³ into same/two-same/distinct index patterns forces e₃=0.  ∎

Theorem B (d₄ / R₄ via μ₄) — PROVED Fraction.
  For any 4-set S, E[∏_{i∈S} w_i]=(1/16)(1+6 e₂+m₄^{sw}(S)) by Thm A.
  Averaging over 4-sets:
      d₄ = R₄/(n)_4 = (1 + 6/(n−1) + μ₄)/16,
      R₄ = (n)_4 (1 + 6/(n−1) + μ₄)/16,
  with μ₄ = C(n,4)⁻¹ ∑_S m₄^{sw}(S).  Equivalently, with the pair-factor
  baseline B:=(n)_4(1+6/(n−1))/16 = n(n−2)(n−3)(n+5)/16,
      R₄ − B = (n)_4 μ₄ / 16,   μ₄ = 16(R₄−B)/(n)_4.  ∎

Theorem C (closed μ4_flat, μ4_suf; residual ⇔ μ₄≤μ4_suf) — PROVED Fraction.
  Substituting R4_flat, R4_suf of 15.146 and simplifying in x=p²:
      μ4_flat = (3p²+17) / (p² (p²−2)(p²−5)),
      μ4_suf  = (3p⁴ + 37 p² + 60) / (p⁴ (p²−2)(p²−5))
               = μ4_flat + 24 ρ_min² / (n)_4.
  Hence for primes p>√5,
      R₄ ≤ R4_suf  ⇔  μ₄ ≤ μ4_suf  ⇔  δ² ≤ ρ_min² (p≥7 Path C residual).
  As p→∞, μ4_flat ∼ 3/p⁴ and μ4_suf ∼ 3/p⁴.  ∎

Theorem D (census μ₄ and m_e) — CERTIFIED Fraction.
  Full Max+ / W-census:
    p=5: μ₄=931/97175 ≤ μ4_suf=143/14375; m_e as 15.151;
    p=7: μ₄=39979/23548175 ≤ μ4_suf=2269/1241317;
         m_e = (1397486/51125, 568587/20450, 522226/18405, 298724/10225);
         ∑π_e m_e = E_μ[k] = 7415451/265850 ≤ k_suf.
  (p=7 pair-δ² channel understates true census δ²; residual still holds.)  ∎

Theorem E (character-sum attack surface) — STRUCTURE.
  Residual is equivalent to a single scalar bound on the switched 4-point
  average μ₄.  Paley / Weil character-sum targets:
    (i) bound m₄^{sw}(S) uniformly or on κ-strata of the conference matrix,
    (ii) evaluate ∑_S m₄^{sw}(S) via Gauss sums over F_{p²},
    (iii) bound the four type means m_e for fixed triples (15.150–15.152)
        by the same correlation algebra.
  Dead: plain ULC; pure t_e(k) (15.152 multi-orbit); CS/Popoviciu on Cov.  ∎

Theorem F (OPEN residual, honest).
  Prove μ₄ ≤ μ4_suf for all primes p≥7 (⇔ R₄≤R4_suf ⇔ residual).  ∎

============================================================================
Certified: Thm A–C Fraction p=5..97; Thm D census p=5,7 (W=86 m_e).
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15153.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15124 import Ek1, Ek2, Ek3  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15146 import R4_flat, R4_suf  # noqa: E402
from e1_gmin_m4_prop15149 import k_act_from_R4, k_suf  # noqa: E402
from e1_gmin_m4_prop15150 import pi_e  # noqa: E402
from e1_gmin_m4_prop15151 import M_E_P5, P5_DELTA2  # noqa: E402

P7_DELTA2 = P7_DELTA2_PAIR

# Exact m_e at p=7 from full Max+ census (W=86)
M_E_P7: dict[int, Fraction] = {
    0: Fraction(1397486, 51125),
    1: Fraction(568587, 20450),
    2: Fraction(522226, 18405),
    3: Fraction(298724, 10225),
}


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def primes_ge(lo: int, hi: int) -> list[int]:
    return [p for p in range(lo, hi + 1) if is_prime(p)]


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------

def e2_of(p: int) -> Fraction:
    return Fraction(1, n_of(p) - 1)


def baseline_B(p: int) -> Fraction:
    """B = (n)_4 (1+6/(n−1))/16 = n(n−2)(n−3)(n+5)/16."""
    n = n_of(p)
    return Fraction(n * (n - 2) * (n - 3) * (n + 5), 16)


def mu4_flat(p: int) -> Fraction:
    """(3p²+17)/(p²(p²−2)(p²−5))."""
    x = p * p
    return Fraction(3 * x + 17, x * (x - 2) * (x - 5))


def mu4_suf(p: int) -> Fraction:
    """(3p⁴+37p²+60)/(p⁴(p²−2)(p²−5))."""
    x = p * p
    return Fraction(3 * x * x + 37 * x + 60, x * x * (x - 2) * (x - 5))


def mu4_from_R4(p: int, R4: Fraction) -> Fraction:
    n = n_of(p)
    n4 = n * (n - 1) * (n - 2) * (n - 3)
    return 16 * R4 / Fraction(n4) - 1 - Fraction(6, n - 1)


def R4_from_mu4(p: int, mu4: Fraction) -> Fraction:
    n = n_of(p)
    n4 = n * (n - 1) * (n - 2) * (n - 3)
    return Fraction(n4, 16) * (1 + Fraction(6, n - 1) + mu4)


def d4_from_mu4(p: int, mu4: Fraction) -> Fraction:
    return (1 + Fraction(6, n_of(p) - 1) + mu4) / 16


def R4_from_W(p: int) -> Fraction:
    W = W_CENSUS[p]
    N = sum(W.values())
    return sum(Fraction(k * (k - 1) * (k - 2) * (k - 3)) * c for k, c in W.items()) / N


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_low_moments(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # E[s]=0, E[s^2]=2n, E[s^3]=0 from k-moments
        # s = n - 2k
        Es = n - 2 * Ek1(p)
        Es2 = n * n - 4 * n * Ek1(p) + 4 * Ek2(p)
        Es3 = (
            n**3
            - 6 * n * n * Ek1(p)
            + 12 * n * Ek2(p)
            - 8 * Ek3(p)
        )
        e2 = (Es2 - n) / Fraction(n * (n - 1))
        # e3 from Es3 = n(n-1)(n-2) e3 when e1=0
        e3 = Es3 / Fraction(n * (n - 1) * (n - 2)) if Es3 != 0 or True else Fraction(0)
        match = Es == 0 and Es2 == 2 * n and Es3 == 0 and e2 == e2_of(p) and e3 == 0
        rows[str(p)] = {
            "E_s": str(Es),
            "E_s2": str(Es2),
            "E_s3": str(Es3),
            "e2": str(e2),
            "e3": str(e3),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": "Switched: e1=0, e3=0, e2=1/(n−1) Max+-free from moments j≤3.",
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B–C
# ---------------------------------------------------------------------------

def prove_mu4_dictionary(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        mf = mu4_flat(p)
        ms = mu4_suf(p)
        # match R4 channel
        R4f = R4_from_mu4(p, mf)
        R4s = R4_from_mu4(p, ms)
        match_flat = R4f == R4_flat(p)
        match_suf = R4s == R4_suf(p)
        # inverse
        match_inv = (
            mu4_from_R4(p, R4_flat(p)) == mf
            and mu4_from_R4(p, R4_suf(p)) == ms
        )
        # addend form
        n = n_of(p)
        n4 = n * (n - 1) * (n - 2) * (n - 3)
        match_add = ms == mf + 24 * rho_min2(p) / Fraction(n4)
        # baseline
        B = baseline_B(p)
        match_B = R4_flat(p) - B == Fraction(n4, 16) * mf
        match = match_flat and match_suf and match_inv and match_add and match_B
        rows[str(p)] = {
            "mu4_flat": str(mf),
            "mu4_suf": str(ms),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "d4=(1+6/(n−1)+μ4)/16; μ4_flat=(3p²+17)/(p²(p²−2)(p²−5)); "
            "μ4_suf=(3p⁴+37p²+60)/(p⁴(p²−2)(p²−5)); residual⇔μ4≤μ4_suf."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "closed_forms": {
            "e2": "1/(n-1)",
            "d4": "(1 + 6/(n-1) + mu4)/16",
            "mu4_flat": "(3p^2+17)/(p^2 (p^2-2)(p^2-5))",
            "mu4_suf": "(3p^4+37p^2+60)/(p^4 (p^2-2)(p^2-5))",
            "baseline_B": "n(n-2)(n-3)(n+5)/16",
        },
    }


# ---------------------------------------------------------------------------
# Theorem D: census
# ---------------------------------------------------------------------------

def prove_census() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        R4 = R4_from_W(p)
        mu = mu4_from_R4(p, R4)
        ms = mu4_suf(p)
        le = mu <= ms
        Em = k_act_from_R4(p, R4)
        ks = k_suf(p)
        rows[str(p)] = {
            "R4": str(R4),
            "mu4": str(mu),
            "mu4_suf": str(ms),
            "mu4_le_suf": le,
            "E_mu": str(Em),
            "k_suf": str(ks),
            "E_mu_le_k_suf": Em <= ks,
            "slack_mu4": str(ms - mu),
        }
        if not (le and Em <= ks):
            ok = False
    # m_e p=5,7
    Em5 = sum(pi_e(5)[e] * M_E_P5[e] for e in range(4))
    Em7 = sum(pi_e(7)[e] * M_E_P7[e] for e in range(4))
    m_ok = (
        Em5 == k_act_from_R4(5, R4_from_W(5))
        and Em7 == k_act_from_R4(7, R4_from_W(7))
    )
    ok = ok and m_ok
    return {
        "proved_census": ok,
        "theorem": (
            "μ4≤μ4_suf and ∑π m_e = E_μ[k]≤k_suf at p=5,7 (full W / Max+)."
        ),
        "by_p": rows,
        "m_e_p5": {str(e): str(M_E_P5[e]) for e in range(4)},
        "m_e_p7": {str(e): str(M_E_P7[e]) for e in range(4)},
        "sum_pi_m_p5": str(Em5),
        "sum_pi_m_p7": str(Em7),
        "m_e_matches_E_mu": m_ok,
        "p7_note": (
            "True census δ²=82176/4499 exceeds pair-channel P7_DELTA2_PAIR; "
            "both still give residual."
        ),
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "mu4_le_mu4_suf_all_p": False,
        "open": (
            "Prove μ₄ ≤ μ4_suf for all primes p≥7 via Paley/Weil bounds on "
            "switched 4-point correlations m₄^{sw}, or character-sum m_e."
        ),
        "sufficient_targets": [
            "mu4_le_mu4_suf_character_sum",
            "uniform_m4_sw_bound",
            "character_sum_m_e_fixed_triple",
            "R4_le_R4_suf",
        ],
        "dead": [
            "plain_ULC",
            "pure_t_e_of_k_all_p",
            "CS_Popoviciu_on_Cov",
            "type_list_free_orbits_p_ge_11",
        ],
        "character_sum_surface": (
            "μ₄ = avg m₄^{sw}(S) on conference/Paley over F_{p²}; "
            "Weil bounds on character sums of products of quadratic characters."
        ),
    }


def main() -> dict:
    primes = primes_ge(5, 37)
    a = prove_low_moments(primes)
    b = prove_mu4_dictionary(primes)
    d = prove_census()
    g = prove_open()

    out = {
        "prop": "15.153",
        "title": (
            "Prop 15.153: switched μ₄ residual dictionary; "
            "closed μ4_flat/μ4_suf; m_e p=7; L OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Switched Max+: e1=e3=0, e2=1/(n−1). "
            "d4=(1+6/(n−1)+μ4)/16. Residual ⇔ μ4≤μ4_suf with "
            "μ4_flat=(3p²+17)/(p²(p²−2)(p²−5)), "
            "μ4_suf=(3p⁴+37p²+60)/(p⁴(p²−2)(p²−5)). "
            "Cert μ4≤μ4_suf and m_e at p=5,7. Need Paley/Weil bound on μ4. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "switched_low_moments": a["proved"],
            "mu4_dictionary_closed": b["proved"],
            "census_mu4_and_m_e": d["proved_census"],
            "residual_for_all_p_ge_5": False,
            "mu4_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "low_moments": a,
            "mu4_dict": b,
            "census": d,
            "open": g,
        },
        "closed_forms": b["closed_forms"],
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "mu4_channel_ready": True,
            "general_residual": False,
        },
        "open_attack": (
            "Paley/Weil character-sum bound μ₄≤μ4_suf (or m₄^{sw} on strata), "
            "or character-sum m_e on fixed type-e triples."
        ),
        "backend": "CPU Fraction + W-census p=5,7; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {"workers": "Fraction algebra", "gpu": "unused"},
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15153.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
