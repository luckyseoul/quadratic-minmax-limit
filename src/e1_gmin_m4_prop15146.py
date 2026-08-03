#!/usr/bin/env python3
"""
Prop 15.146 — Type-free R₄ / central-moment residual channel; closed
ED4_from_exact3; δ²=(2/3)(R₄−R4_flat); spectral mass at 4p too weak for
Jensen; residual still OPEN for general p≥7.

Continues 15.145. Does **not** soft-close residual. No type census.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. After halfspace Seidel switch (15.123–15.124): k=wt(w) on
X=V₊∩{0,1}ⁿ ≅ Max+, D=n−2k, ED4=E[D⁴].
E[k],E[k²],E[k³] Max+-free (15.124.A). R₄:=E[k⁴]−exact_≤3 with
exact_≤3 the ≤3-collision part of E[k⁴]. Path C residual uses pair-average
forms; for p≥7, δ²≤ρ_min² suffices (15.130/15.145).

============================================================================
Theorem A (closed ED4_from_exact3) — PROVED Fraction, all n=p²+1.
  Substituting E[k]=n/2, E[k²]=n(n+2)/4, E[k³]=n²(n+6)/8 and
  E[k⁴]=exact_≤3 into the binomial expansion of E[(n−2k)⁴] yields
      ED4_from_exact3 = −n⁴ + 28 n² − 40 n
                       = −n(n³ − 28 n + 40).
  (Independent of the residual 4-wise mass R₄.)  ∎

Theorem B (R₄ residual dictionary, type-free) — PROVED Fraction.
  Define Max+-free
      R4_flat := (ED4_flat − ED4_from_exact3)/16,
      R4_suf  := (ED4_suf  − ED4_from_exact3)/16
               = R4_flat + (3/2) ρ_min²,
      R4_bud  := (ED4_bud  − ED4_from_exact3)/16
               = R4_flat + (3/2)·(room_hyp/24).
  For the pair-average residual R₄ of Max+,
      ED4 = ED4_from_exact3 + 16 R₄,
      δ²  = (2/3)(R₄ − R4_flat),
  hence
      δ² ≤ ρ_min²  ⇔  R₄ ≤ R4_suf,
      δ² ≤ room_hyp/24  ⇔  R₄ ≤ R4_bud.
  For primes p≥7, R₄ ≤ R4_suf ⇒ Path C residual.  Also
      R₄ = n(n−1)(n−2)(n−3)·Ē[w_i w_j w_l w_m]
  (ordered distinct 4-wise density of regular-set indicators), so residual
  ⇔ Ē[∏₄ w] ≤ R4_suf / (n)₄.  ∎

Theorem C (central weight moments) — PROVED Fraction.
  With μ:=E[k]=n/2:
      μ₂ := Var(k) = n/2,     μ₃ := E[(k−μ)³] = 0
  (Max+-free: from E[kʲ] j≤3 and antipodal symmetry W_k=W_{n−k}).
  The fourth central moment satisfies μ₄ = E[k⁴] − 4μ E[k³] + 6μ² E[k²] − 3μ⁴,
  and since E[k⁴]=exact_≤3+R₄ with exact_≤3 Max+-free,
      μ₄ − μ4_flat = R₄ − R4_flat = (3/2) δ²,
  so δ²≤ρ_min² ⇔ μ₄ ≤ μ4_suf := μ4_flat + (3/2)ρ_min².  ∎

Theorem D (spectral mass of f_y at eigenvalue 4p) — PROVED Fraction.
  For every y∈Max+, the spectral measure of f_y relative to T has
  mean m₁=4p−12/p and variance
      v = 24(p²−3)(p²−4)/(p²(p²−2))
  (Prop 15.115.E, Max+-free).  If λ_max(T)=4p (certified p=5,7; structure
  when 4p is an eigenvalue of T), then for any probability measure on
  (−∞,4p] with this mean and variance, the mass at 4p is at most
      w* = v / (v + (12/p)²) = (p²−3)(p²−4) / (p²(p²−1)).
  Hence
      ‖P_{E_{4p}} f_y‖₂² ≤ w* · \\binom{n}{4}
        = (p²+1)(p²−2)(p²−3)(p²−4)/24.
  This upper bound is Θ(p⁸) while ρ_min²=Θ(p²), so
      E‖P f_y‖₂² ≤ w* C(n,4)  ≫  ρ_min²
  for all primes p≥5: the Jensen path δ²≤E‖P f‖² cannot close residual
  from this spectral-moment bound alone (confirms 15.115.G quantitatively
  with a closed Max+-free UB).  ∎

Theorem E (OPEN residual, honest).
  Prove R₄ ≤ R4_suf (⇔ Ē[∏₄ w] ≤ R4_suf/(n)₄ ⇔ μ₄ ≤ μ4_suf ⇔ δ²≤ρ_min²)
  for all primes p≥7 by a type-free 4-wise bound on regular sets in V₊
  (weight enumerator / association scheme / Gauss sums).  Dead for residual:
  Jensen + spectral mass w* (Thm D); type-list residual p≥11 (15.144);
  crude envelopes of 15.120.  Census: R₄≤R4_suf at p=5,7 (prior W / pair).  ∎

============================================================================
Certified: Thm A–C Fraction all primes 5..97; Thm D algebra + weakness;
  census R₄ channel at p=5,7.
Backend: CPU Fraction; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15146.json
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
from e1_gmin_m4_prop15112 import ed4_suf  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15124 import (  # noqa: E402
    Ek1,
    Ek2,
    Ek3,
    exact_le3,
    ed4_from_exact3 as ed4_from_exact3_ref,
)
from e1_gmin_m4_prop15128 import W_CENSUS, moments_from_W  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15145 import ed4_suf_closed  # noqa: E402

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
# Theorem A: closed ED4_from_exact3
# ---------------------------------------------------------------------------

def ed4_from_exact3_closed(p: int) -> Fraction:
    """−n⁴ + 28 n² − 40 n."""
    n = n_of(p)
    return -Fraction(n**4) + 28 * Fraction(n * n) - 40 * Fraction(n)


def prove_ed4_exact3(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        c = ed4_from_exact3_closed(p)
        ref = ed4_from_exact3_ref(p)
        # direct expansion check
        n = n_of(p)
        e1, e2, e3 = Ek1(p), Ek2(p), Ek3(p)
        e4 = exact_le3(p)["exact_le3"]
        direct = (
            Fraction(n**4)
            - 8 * Fraction(n**3) * e1
            + 24 * Fraction(n**2) * e2
            - 32 * Fraction(n) * e3
            + 16 * e4
        )
        match = c == ref == direct
        rows[str(p)] = {
            "ED4_from_exact3": str(c),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": "ED4_from_exact3 = −n⁴+28n²−40n (Max+-free).",
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem B: R4 dictionary
# ---------------------------------------------------------------------------

def R4_flat(p: int) -> Fraction:
    return (ed4_flat_closed(p) - ed4_from_exact3_closed(p)) / 16


def R4_suf(p: int) -> Fraction:
    return (ed4_suf_closed(p) - ed4_from_exact3_closed(p)) / 16


def R4_bud(p: int) -> Fraction:
    return (ed4_bud_closed(p) - ed4_from_exact3_closed(p)) / 16


def bar_w4_of_R4(p: int, R4: Fraction) -> Fraction:
    n = n_of(p)
    return R4 / Fraction(n * (n - 1) * (n - 2) * (n - 3))


def prove_R4_dictionary(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        rf, rs, rb = R4_flat(p), R4_suf(p), R4_bud(p)
        # slack identities
        slack_suf = rs - rf == Fraction(3, 2) * rho_min2(p)
        slack_bud = rb - rf == Fraction(3, 2) * room_hyp_over_24(p)
        # consistency with ED4 forms
        cons = (
            ed4_from_exact3_closed(p) + 16 * rf == ed4_flat_closed(p)
            and ed4_from_exact3_closed(p) + 16 * rs == ed4_suf_closed(p)
            and ed4_from_exact3_closed(p) + 16 * rb == ed4_bud_closed(p)
            and ed4_suf_closed(p) == ed4_suf(p)
        )
        rows[str(p)] = {
            "R4_flat": str(rf),
            "R4_suf": str(rs),
            "R4_bud": str(rb),
            "bar_w4_suf": str(bar_w4_of_R4(p, rs)),
            "slack_suf_is_3_2_rho": slack_suf,
            "slack_bud_is_3_2_room": slack_bud,
            "match": slack_suf and slack_bud and cons,
        }
        if not (slack_suf and slack_bud and cons):
            ok = False
        if p >= 7 and not (rho_min2(p) < room_hyp_over_24(p)):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "δ²=(2/3)(R₄−R4_flat); δ²≤ρ_min²⇔R₄≤R4_suf=R4_flat+(3/2)ρ_min²; "
            "R₄=n(n−1)(n−2)(n−3) Ē[∏₄ w]."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem C: central moments
# ---------------------------------------------------------------------------

def mu2(p: int) -> Fraction:
    return Fraction(n_of(p), 2)


def mu3(p: int) -> Fraction:
    return Fraction(0)


def mu4_from_e4(p: int, e4: Fraction) -> Fraction:
    mu = Ek1(p)
    return e4 - 4 * mu * Ek3(p) + 6 * mu * mu * Ek2(p) - 3 * mu**4


def mu4_flat(p: int) -> Fraction:
    e4 = exact_le3(p)["exact_le3"] + R4_flat(p)
    return mu4_from_e4(p, e4)


def mu4_suf(p: int) -> Fraction:
    e4 = exact_le3(p)["exact_le3"] + R4_suf(p)
    return mu4_from_e4(p, e4)


def prove_central_moments(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # μ2 = Ek2 - Ek1^2 = n(n+2)/4 - n^2/4 = n/2
        v = Ek2(p) - Ek1(p) ** 2
        # μ3 = Ek3 - 3μ Ek2 + 2 μ^3
        mu = Ek1(p)
        m3 = Ek3(p) - 3 * mu * Ek2(p) + 2 * mu**3
        m2_ok = v == Fraction(n, 2)
        m3_ok = m3 == 0
        slack = mu4_suf(p) - mu4_flat(p) == Fraction(3, 2) * rho_min2(p)
        rows[str(p)] = {
            "mu2": str(v),
            "mu3": str(m3),
            "mu4_flat": str(mu4_flat(p)),
            "mu4_suf": str(mu4_suf(p)),
            "match": m2_ok and m3_ok and slack,
        }
        if not (m2_ok and m3_ok and slack):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "μ₂=n/2, μ₃=0 Max+-free; μ₄−μ4_flat=(3/2)δ²; "
            "δ²≤ρ_min²⇔μ₄≤μ4_suf."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem D: spectral mass bound (too weak)
# ---------------------------------------------------------------------------

def w_star(p: int) -> Fraction:
    """Max mass at λ=4p given spectral mean/var of f_y, support ≤4p."""
    return Fraction((p * p - 3) * (p * p - 4), p * p * (p * p - 1))


def Pf_energy_ub(p: int) -> Fraction:
    """w* · C(n,4) = (p²+1)(p²−2)(p²−3)(p²−4)/24."""
    return Fraction(
        (p * p + 1) * (p * p - 2) * (p * p - 3) * (p * p - 4),
        24,
    )


def var_spec(p: int) -> Fraction:
    return Fraction(24 * (p * p - 3) * (p * p - 4), p * p * (p * p - 2))


def prove_spectral_mass(primes: list[int]) -> dict:
    rows = {}
    ok = True
    weak_all = True
    for p in primes:
        if p < 5:
            continue
        v = var_spec(p)
        gap = Fraction(12, p)  # 4p - m1
        ws = v / (v + gap * gap)
        ws_c = w_star(p)
        ub = Pf_energy_ub(p)
        Cn4 = comb(n_of(p), 4)
        match = ws == ws_c and ub == ws_c * Cn4
        # weakness: ub >> rho_min2 and room
        weak = ub > rho_min2(p) and ub > room_hyp_over_24(p)
        rows[str(p)] = {
            "w_star": str(ws_c),
            "Pf_energy_ub": str(ub),
            "rho_min2": str(rho_min2(p)),
            "room": str(room_hyp_over_24(p)),
            "ub_over_rho": float(ub / rho_min2(p)),
            "too_weak_for_jensen_residual": weak,
            "match": match,
        }
        if not match:
            ok = False
        if not weak:
            weak_all = False
    return {
        "proved": ok and weak_all,
        "theorem": (
            "Under λ_max(T)=4p: ‖P f_y‖²≤(p²+1)(p²−2)(p²−3)(p²−4)/24; "
            "≪̸ residual budgets (Jensen path dead from spectral moments alone)."
        ),
        "closed_w_star": "(p^2-3)(p^2-4)/(p^2(p^2-1))",
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
        "dead_for_residual": True,
    }


# ---------------------------------------------------------------------------
# Theorem E: census + OPEN
# ---------------------------------------------------------------------------

def prove_census() -> dict:
    """W-census R4 channel at p=5,7 (basepoint hs enumerator)."""
    rows = {}
    ok = True
    for p in (5, 7):
        e4 = moments_from_W(W_CENSUS[p], 4)
        R4 = e4 - exact_le3(p)["exact_le3"]
        le_suf = R4 <= R4_suf(p)
        le_bud = R4 <= R4_bud(p)
        # pair-average δ² channel (independent of W basepoint)
        d2 = P5_DELTA2 if p == 5 else P7_DELTA2
        pair_R4 = R4_flat(p) + Fraction(3, 2) * d2
        pair_le = pair_R4 <= R4_suf(p)
        rows[str(p)] = {
            "R4_from_W": str(R4),
            "R4_W_le_suf": le_suf,
            "R4_W_le_bud": le_bud,
            "delta2_pair": str(d2),
            "R4_from_pair_delta2": str(pair_R4),
            "pair_R4_le_suf": pair_le,
        }
        if not (le_suf and le_bud and pair_le):
            ok = False
    return {
        "proved": ok,
        "theorem": "R₄≤R4_suf at p=5,7 (W-census and pair δ² channel).",
        "by_p": rows,
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "type_free_R4_le_R4_suf_all_p": False,
        "open": (
            "Prove R₄≤R4_suf (⇔ Ē[∏₄ w]≤R4_suf/(n)₄ ⇔ μ₄≤μ4_suf) for all "
            "primes p≥7 by type-free 4-wise regular-set bounds."
        ),
        "dead": [
            "jensen_plus_spectral_mass_w_star",
            "type_list_residual_p_ge_11",
            "crude_majorization_and_min_distance_envelopes",
        ],
        "preferred_attacks": [
            "closed Max+ weight enumerator W_k",
            "4-wise density bound for V+ ∩ {0,1}^n regular sets",
            "Gauss sums for m4 / character-sum c_j",
        ],
    }


def main() -> dict:
    primes = primes_ge(5, 97)
    a = prove_ed4_exact3(primes)
    b = prove_R4_dictionary(primes)
    c = prove_central_moments(primes)
    d = prove_spectral_mass(primes)
    e = prove_census()
    f = prove_open()

    out = {
        "prop": "15.146",
        "title": (
            "Prop 15.146: Type-free R₄/μ₄ residual channel; closed "
            "ED4_from_exact3; spectral mass w* too weak for Jensen"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Type-free residual rewritten as R₄≤R4_suf "
            "(⇔μ₄≤μ4_suf⇔Ē[∏₄w] bound) with closed Max+-free budgets. "
            "Spectral mass bound on ‖P f_y‖ is Max+-free but far too weak "
            "for Jensen residual. Residual closed only p=5,7. L OPEN."
        ),
        "proved": {
            "ed4_from_exact3_closed": a["proved"],
            "R4_residual_dictionary": b["proved"],
            "central_moments_mu2_mu3": c["proved"],
            "spectral_mass_w_star_too_weak": d["proved"],
            "census_R4_p5_p7": e["proved"],
            "residual_for_all_p_ge_5": False,
            "type_free_R4_le_R4_suf_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
        },
        "algebra": {
            "ed4_exact3": a,
            "R4_dictionary": b,
            "central_moments": c,
            "spectral_mass": d,
            "census": e,
            "open": f,
        },
        "closed_forms": {
            "ED4_from_exact3": "-n^4 + 28 n^2 - 40 n",
            "R4_suf": "R4_flat + (3/2) rho_min2",
            "R4_bud": "R4_flat + (3/2) room_hyp/24",
            "delta2": "(2/3)(R4 - R4_flat)",
            "w_star": "(p^2-3)(p^2-4)/(p^2(p^2-1))",
            "Pf_energy_ub": "(p^2+1)(p^2-2)(p^2-3)(p^2-4)/24",
            "sample": {
                str(p): {
                    "ED4_from_exact3": str(ed4_from_exact3_closed(p)),
                    "R4_flat": str(R4_flat(p)),
                    "R4_suf": str(R4_suf(p)),
                    "R4_bud": str(R4_bud(p)),
                    "mu4_suf": str(mu4_suf(p)),
                    "w_star": str(w_star(p)),
                    "Pf_ub": str(Pf_energy_ub(p)),
                }
                for p in (5, 7, 11, 13)
            },
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "jensen_spectral_mass_viable": False,
            "general_residual": False,
        },
        "open_attack": (
            "Prove R₄≤R4_suf for p≥7 via 4-wise regular-set density / "
            "weight enumerator / Gauss sums (type-free)."
        ),
        "backend": "CPU Fraction; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {
            "workers": "Fraction algebra (no multi-minute CPU)",
            "gpu": "unused",
            "primes_checked": "5..97",
        },
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15146.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
