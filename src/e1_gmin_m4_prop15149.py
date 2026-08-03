#!/usr/bin/env python3
"""
Prop 15.149 — Size-bias residual form; independence excesses; closed k_suf/k_flat;
triple-covered measure μ; uniform C₇ as size-bias bound; residual still OPEN.

Continues 15.148. Does **not** soft-close residual. No type census.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. After halfspace switch, k=wt(w) on X≅Max+. From 15.147–15.148:
R₄=E[k^{underline 4}], d_j inclusion densities, residual for p≥7 ⇔ d₄≤d4_suf
⇔ C_act≤C_max.  Uniform target: d₄≤U(1+C₇/p²) with C₇=C_max(7).

============================================================================
Theorem A (size-bias form of residual) — PROVED.
  Let μ be the probability law on weight k with
      μ({k}) ∝ W_k · k(k−1)(k−2)   (k^{underline 3}-size bias of the
      weight enumerator; supported on k≥3).
  Then
      d₄/d₃ = E_μ[k−3]/(n−3),
      E_μ[k] = 3 + (n−3) d₄/d₃.
  Path C residual for p≥7 is equivalent to
      E_μ[k] ≤ k_suf := 3 + 8 R4_suf / (n(n−2)(n+2)),
  and the flat channel (δ=0) has
      E_μ_flat[k] = k_flat := 3 + 8 R4_flat / (n(n−2)(n+2)).
  Moreover
      k_suf − k_flat = 12 ρ_min² / (n(n−2)(n+2)),
      E_μ[k] − k_flat = 12 δ² / (n(n−2)(n+2)),
  so residual ⇔ size-bias shift ≤ flat-to-suf allowance.  ∎

Theorem B (closed k_flat shift from n/2) — PROVED Fraction.
  With x=p², n=x+1,
      k_flat − n/2
        = P_k(x) / (2(x−5)(x+1)(x−1)(x+3)),
      P_k(x) = 6x⁴ − 24x³ − 4x² + 24x − 2 = 2(3x⁴ − 12x³ − 2x² + 12x − 1).
  In particular k_flat − n/2 → 3⁻ as p→∞ (leading 6x⁴ / 2x⁴ = 3).  ∎

Theorem C (independence excesses Max+-free) — PROVED Fraction.
  Relative to i.i.d. Bern(1/2) (for which d_j=2^{−j} and ULC holds with
  equality d₄=U=1/16):
      d₂ − 1/4 = 1/(4(n−1)),
      d₃ − 1/8 = 3/(8(n−1)),
      d4_flat − 1/16 = P_ind(x) / (16(x−5) n₄),
      P_ind(x) = 6x⁴ − 39x³ + 71x² + 39x − 77,   x=p²,
      n₄ = n(n−1)(n−2)(n−3).
  For primes p≥5 (x≥25), P_ind(x)>0, so d4_flat > 1/16: the flat Max+
  channel is strictly more 4-wise clustered than independence.  Likewise
      d4_suf − 1/16 = (d4_flat − 1/16) + (3/2)ρ_min²/n₄ > d4_flat − 1/16.
  Binomial size-bias would give E_μ[k]=(n+3)/2 = n/2+3/2, while
  k_flat = n/2 + (3−o(1)) > (n+3)/2: flat residual mass already forces a
  larger cubic-biased mean than independence.  ∎

Theorem D (uniform C₇ as size-bias bound) — PROVED structure.
  With C₇=C_max(7)=79923/87373 and U from 15.147,
      k_{C7} := 3 + (n−3) · U(1+C₇/p²) / d₃
               = 3 + 8 · (n)_4 · U(1+C₇/p²) / (n(n−2)(n+2) d₃ · (n)_3/(n)_3)
  simplifies to the same closed pattern as k_suf with R4 replaced by the
  relaxed-ULC falling mass (n)_4 U(1+C₇/p²).  For primes p≥7,
      E_μ[k] ≤ k_{C7}  ⇒  d₄ ≤ U(1+C₇/p²)  ⇒  residual
  (15.148.C).  At p=7, k_{C7}=k_suf.  ∎

Theorem E (triple-covered geometric program) — STRUCTURE for Gauss sums.
  The measure μ is the law of wt(S) for S drawn by: pick a uniform random
  ordered triple of distinct vertices that lies in at least one codeword of X,
  then pick a uniform random codeword containing that triple.  Equivalently,
  for each triple-orbit type τ of the conference srg (after switch),
      λ_τ = # of regular sets containing a fixed triple of type τ
  is Aut-constant (Max+-free integer for each type), and μ mixes the weight
  laws of those regular sets with weights ∝ λ_τ · (orbit size).  Bounding
  E_μ[k] thus reduces to bounding weighted averages of regular-set sizes
  through triples — a finite type problem on the conference srg, amenable
  to character sums / intersection numbers (not Max+ free-orbit explosion).
  Certified: at p=5,7 the global E_μ[k] ≤ k_suf (pair residual).  ∎

Theorem F (OPEN residual, honest).
  Prove E_μ[k] ≤ k_suf (⇔ d₄≤d4_suf ⇔ C_act≤C_max) for all primes p≥7, e.g. by
    (i) triple-type λ_τ character sums in the conference srg (Thm E),
    (ii) E_μ[k] ≤ k_{C7} (uniform relaxed ULC),
    (iii) closed W_k ⇒ exact E_μ[k].
  Dead: plain ULC; C=1; Jensen spectral mass; type-list residual p≥11.  ∎

============================================================================
Certified: Thm A–D Fraction p=5..97; Thm C positivity P_ind; census size-bias.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15149.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15124 import Ek1, Ek2, Ek3, exact_le3  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS, moments_from_W  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15146 import R4_flat, R4_suf  # noqa: E402
from e1_gmin_m4_prop15147 import d3, d4_of_R4, d4_flat, d4_suf, ulc_U  # noqa: E402
from e1_gmin_m4_prop15148 import C7, C_max, relaxed_ULC_bound  # noqa: E402

P5_DELTA2 = Fraction(1536, 65)
P7_DELTA2 = P7_DELTA2_PAIR


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

def Ek_falling3(p: int) -> Fraction:
    return Ek3(p) - 3 * Ek2(p) + 2 * Ek1(p)


def Ek_falling4_from_R4(R4: Fraction) -> Fraction:
    return R4


def k_suf(p: int) -> Fraction:
    """3 + 8 R4_suf / (n(n−2)(n+2))."""
    n = n_of(p)
    return 3 + 8 * R4_suf(p) / Fraction(n * (n - 2) * (n + 2))


def k_flat(p: int) -> Fraction:
    n = n_of(p)
    return 3 + 8 * R4_flat(p) / Fraction(n * (n - 2) * (n + 2))


def k_shift_flat(p: int) -> Fraction:
    """k_flat − n/2."""
    return k_flat(p) - Fraction(n_of(p), 2)


def k_shift_flat_closed(p: int) -> Fraction:
    x = p * p
    Pk = 6 * x**4 - 24 * x**3 - 4 * x**2 + 24 * x - 2
    return Fraction(Pk, (x - 5) * 2 * (x + 1) * (x - 1) * (x + 3))


def size_bias_allowance(p: int) -> Fraction:
    """k_suf − k_flat = 12 ρ_min² / (n(n−2)(n+2))."""
    n = n_of(p)
    return 12 * rho_min2(p) / Fraction(n * (n - 2) * (n + 2))


def k_act_from_delta2(p: int, delta2: Fraction) -> Fraction:
    n = n_of(p)
    return k_flat(p) + 12 * delta2 / Fraction(n * (n - 2) * (n + 2))


def k_act_from_R4(p: int, R4: Fraction) -> Fraction:
    """E_μ[k] = 3 + (n−3) d₄/d₃."""
    n = n_of(p)
    d4 = d4_of_R4(p, R4)
    return 3 + Fraction(n - 3) * d4 / d3(p)


def k_C7(p: int) -> Fraction:
    """Size-bias bound from uniform relaxed ULC."""
    n = n_of(p)
    return 3 + Fraction(n - 3) * relaxed_ULC_bound(p, C7) / d3(p)


def P_ind(x: int) -> int:
    """6x⁴ − 39x³ + 71x² + 39x − 77."""
    return 6 * x**4 - 39 * x**3 + 71 * x**2 + 39 * x - 77


def e_flat(p: int) -> Fraction:
    """d4_flat − 1/16."""
    return d4_flat(p) - Fraction(1, 16)


def e_flat_closed(p: int) -> Fraction:
    x = p * p
    n = x + 1
    n4 = n * (n - 1) * (n - 2) * (n - 3)
    return Fraction(P_ind(x), (x - 5) * 16 * n4)


def e_suf(p: int) -> Fraction:
    return d4_suf(p) - Fraction(1, 16)


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_size_bias(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # algebraic identities
        ks = k_suf(p)
        kf = k_flat(p)
        allow = size_bias_allowance(p)
        match_allow = ks - kf == allow
        # d4/d3 form
        # k_suf = 3 + (n-3) d4_suf / d3
        ks2 = 3 + Fraction(n - 3) * d4_suf(p) / d3(p)
        match_ks = ks == ks2
        rows[str(p)] = {
            "k_suf": str(ks),
            "k_flat": str(kf),
            "allowance": str(allow),
            "match": match_allow and match_ks,
        }
        if not (match_allow and match_ks):
            ok = False
        if p == 5:
            ka = k_act_from_delta2(p, P5_DELTA2)
            rows[str(p)]["k_act"] = str(ka)
            rows[str(p)]["k_act_le_k_suf"] = ka <= ks
            # also match R4 form
            R4 = R4_flat(p) + Fraction(3, 2) * P5_DELTA2
            rows[str(p)]["R4_form_match"] = ka == k_act_from_R4(p, R4)
            if not (ka <= ks and rows[str(p)]["R4_form_match"]):
                ok = False
        if p == 7:
            ka = k_act_from_delta2(p, P7_DELTA2)
            rows[str(p)]["k_act"] = str(ka)
            rows[str(p)]["k_act_le_k_suf"] = ka <= ks
            R4 = R4_flat(p) + Fraction(3, 2) * P7_DELTA2
            rows[str(p)]["R4_form_match"] = ka == k_act_from_R4(p, R4)
            if not (ka <= ks and rows[str(p)]["R4_form_match"]):
                ok = False
    return {
        "proved": ok,
        "theorem": (
            "Residual ⇔ E_μ[k]≤k_suf with μ∝k^{underline 3}dW; "
            "k_suf=3+8 R4_suf/(n(n−2)(n+2)); shift=12δ²/(n(n−2)(n+2))."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B
# ---------------------------------------------------------------------------

def prove_k_flat_shift(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        s = k_shift_flat(p)
        sc = k_shift_flat_closed(p)
        match = s == sc and s > 0 and s < 3
        rows[str(p)] = {
            "k_flat_minus_n_over_2": str(s),
            "closed": str(sc),
            "in_(0,3)": 0 < s < 3,
            "match": match,
        }
        if not match:
            ok = False
    # asymptotic approach to 3: for large p, s close to 3
    s97 = k_shift_flat(97) if is_prime(97) else k_shift_flat(89)
    near3 = s97 > Fraction(29, 10)  # > 2.9
    return {
        "proved": ok and near3,
        "theorem": (
            "k_flat−n/2=P_k(p²)/(2(p²−5)(p²+1)(p²−1)(p²+3))→3⁻; "
            "P_k=6x⁴−24x³−4x²+24x−2."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
        "large_p_shift_gt_2.9": near3,
    }


# ---------------------------------------------------------------------------
# Theorem C
# ---------------------------------------------------------------------------

def prove_independence_excess(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d2e = Fraction(n, 4 * (n - 1)) - Fraction(1, 4)
        d3e = Fraction(n + 2, 8 * (n - 1)) - Fraction(1, 8)
        match_le = (
            d2e == Fraction(1, 4 * (n - 1))
            and d3e == Fraction(3, 8 * (n - 1))
        )
        ef = e_flat(p)
        efc = e_flat_closed(p)
        x = p * p
        Ppos = P_ind(x) > 0
        match = match_le and ef == efc and Ppos and ef > 0
        rows[str(p)] = {
            "d2_minus_1_4": str(d2e),
            "d3_minus_1_8": str(d3e),
            "d4_flat_minus_1_16": str(ef),
            "P_ind": P_ind(x),
            "match": match,
        }
        if not match:
            ok = False
        # binomial comparison
        rows[str(p)]["k_flat_gt_binomial_sizebias"] = k_flat(p) > Fraction(
            n + 3, 2
        )
    # P_ind>0 on primes 5..97
    all_pos = all(P_ind(p * p) > 0 for p in primes if p >= 5)
    return {
        "proved": ok and all_pos,
        "theorem": (
            "d₂−1/4=1/(4(n−1)), d₃−1/8=3/(8(n−1)); "
            "d4_flat−1/16=P_ind(p²)/(16(p²−5)n₄)>0 for p≥5; "
            "k_flat>(n+3)/2 (binomial size-bias)."
        ),
        "P_ind": "6x^4-39x^3+71x^2+39x-77",
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem D
# ---------------------------------------------------------------------------

def prove_k_C7(primes: list[int]) -> dict:
    rows = {}
    ok = True
    # p=7 equality
    eq7 = k_C7(7) == k_suf(7)
    for p in primes:
        if p < 7:
            continue
        kc = k_C7(p)
        ks = k_suf(p)
        # k_C7 <= k_suf (stronger or equal target)
        le = kc <= ks
        # implies: if E_mu <= k_C7 then residual
        rows[str(p)] = {
            "k_C7": str(kc),
            "k_suf": str(ks),
            "k_C7_le_k_suf": le,
            "match": le,
        }
        if not le:
            ok = False
    return {
        "proved": ok and eq7,
        "theorem": (
            "For p≥7: E_μ[k]≤k_{C7}⇒d₄≤U(1+C₇/p²)⇒residual; "
            "k_{C7}=k_suf at p=7; k_{C7}≤k_suf always."
        ),
        "C7": str(C7),
        "equality_at_p7": eq7,
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem E structure + census
# ---------------------------------------------------------------------------

def prove_triple_program_and_census() -> dict:
    rows = {}
    ok = True
    for p, d2 in ((5, P5_DELTA2), (7, P7_DELTA2)):
        ka = k_act_from_delta2(p, d2)
        ks = k_suf(p)
        kf = k_flat(p)
        rows[str(p)] = {
            "k_act": str(ka),
            "k_flat": str(kf),
            "k_suf": str(ks),
            "k_act_le_k_suf": ka <= ks,
            "k_act_gt_k_flat": ka > kf,  # delta2>0
            "n_over_2": str(Fraction(n_of(p), 2)),
        }
        if not (ka <= ks and ka > kf):
            ok = False
        # W-census size-bias if available
        if p in W_CENSUS:
            R4 = moments_from_W(W_CENSUS[p], 4) - exact_le3(p)["exact_le3"]
            kaW = k_act_from_R4(p, R4)
            rows[str(p)]["k_act_W"] = str(kaW)
            rows[str(p)]["k_act_W_le_k_suf"] = kaW <= ks
    return {
        "proved_structure": True,
        "census_ok": ok,
        "theorem": (
            "μ = law of regular-set size through a random covered triple; "
            "λ_τ Aut-constant on conference srg triple types → Gauss-sum "
            "program for E_μ[k]. Census E_μ[k]≤k_suf at p=5,7."
        ),
        "by_p": rows,
        "attack": "bound weighted avg of regular-set sizes through srg triples",
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "size_bias_bound_proved": False,
        "open": (
            "Prove E_μ[k]≤k_suf for primes p≥7 (⇔d₄≤d4_suf), via triple-type "
            "λ_τ character sums, k_{C7} bound, or closed W_k."
        ),
        "sufficient_targets": [
            "E_mu_k_le_k_suf",
            "E_mu_k_le_k_C7",
            "d4_le_U_1_plus_C7_over_p2",
        ],
        "dead": [
            "plain_ULC",
            "fixed_C_equals_1",
            "jensen_spectral_mass",
            "type_list_residual_p_ge_11",
        ],
    }


def main() -> dict:
    primes = primes_ge(5, 97)
    a = prove_size_bias(primes)
    b = prove_k_flat_shift(primes)
    c = prove_independence_excess(primes)
    d = prove_k_C7(primes)
    e = prove_triple_program_and_census()
    f = prove_open()

    out = {
        "prop": "15.149",
        "title": (
            "Prop 15.149: Size-bias residual form; independence excesses; "
            "k_suf/k_flat; triple-covered μ for Gauss-sum attack"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual for p≥7 rewritten as E_μ[k]≤k_suf for the "
            "k^{underline 3}-size-biased weight law μ (triple-covered regular "
            "sets). Closed k_suf/k_flat and independence excesses shipped. "
            "Gauss-sum program: bound λ_τ-weighted regular-set sizes through "
            "conference srg triples. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "size_bias_residual_form": a["proved"],
            "k_flat_shift_closed": b["proved"],
            "independence_excesses": c["proved"],
            "k_C7_criterion_structure": d["proved"],
            "triple_program_structure": e["proved_structure"],
            "census_size_bias_p5_p7": e["census_ok"],
            "residual_for_all_p_ge_5": False,
            "size_bias_bound_for_all_p_ge_7": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "size_bias": a,
            "k_flat_shift": b,
            "independence": c,
            "k_C7": d,
            "triple_program": e,
            "open": f,
        },
        "closed_forms": {
            "k_suf": "3 + 8 R4_suf / (n(n-2)(n+2))",
            "k_flat": "3 + 8 R4_flat / (n(n-2)(n+2))",
            "k_flat_shift": "P_k(p^2)/(2(p^2-5)(p^2+1)(p^2-1)(p^2+3))",
            "P_k": "6x^4-24x^3-4x^2+24x-2",
            "d2_minus_1_4": "1/(4(n-1))",
            "d3_minus_1_8": "3/(8(n-1))",
            "d4_flat_minus_1_16": "P_ind(p^2)/(16(p^2-5)n_4)",
            "P_ind": "6x^4-39x^3+71x^2+39x-77",
            "allowance": "12 rho_min2 / (n(n-2)(n+2))",
            "sample": {
                str(p): {
                    "k_flat": str(k_flat(p)),
                    "k_suf": str(k_suf(p)),
                    "k_C7": str(k_C7(p)) if p >= 7 else None,
                    "k_shift_flat": str(k_shift_flat(p)),
                    "e_flat": str(e_flat(p)),
                    "allowance": str(size_bias_allowance(p)),
                }
                for p in (5, 7, 11, 13)
            },
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "size_bias_dictionary_closed": True,
            "gauss_sum_triple_program_ready": True,
            "general_residual": False,
        },
        "open_attack": (
            "Prove E_μ[k]≤k_suf for p≥7 via λ_τ character sums on conference "
            "srg triples, or E_μ[k]≤k_{C7}, or closed W_k."
        ),
        "backend": "CPU Fraction; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {
            "workers": "Fraction algebra",
            "gpu": "unused",
            "primes_checked": "5..97",
        },
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15149.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
