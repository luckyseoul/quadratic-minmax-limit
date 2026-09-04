#!/usr/bin/env python3
"""
Prop 15.168 — E(1) structure with the corrected 15.720 bi-tight close.

REAL (checkable Fraction / prior-prop predicates):
  A. The required bi-tight levels 2 and 3 are empty by the 15.720 degree
     congruence plus ker(Gsum)=scheme+cross from 15.272/15.207.
  B. Deep tight empty for p≥5 by 15.44.3 plus the level-2 obstruction.
  C. Type I freeness ND — restated from 15.43.1 (prior proved).
  D. Type I freeness-fail k=2p−1 → tight 2p → ND when bi-tight empty
     (15.43.3 + 15.44 master lemma + 15.720).
  E. Deep non-tight freeness ⇒ ND (score eval; prior style 15.43).
  F. Auto-freeness for deep s₊=2 when k≤3p−2 (Fraction inequality).
  G. Deep freeness-fail at k=3p−1 ⇒ S≡3 tight size 3p ⇒ empty when
     its bi-tight alternative is impossible (15.720, level 3).

CURRENT GATE (no soft-close):
  (i) Multi-level Type I is closed independently by Proposition 15.750.
      The older 15.275 `3A+B` split remains incomplete as a mechanism but is
      not a live global gate.
  (ii) Non-Walsh residual (ii) for even k≥4p remains open. The bounded even
       range through 4p−2 is closed; the multi-level remainder is not.
  (iii) Proposition 15.764 exposes a separate minimal-four-gap implication
        bridge outside the historical four-unit ledger.
  Full m_n≥Φ−2 / E(1) / L remain OPEN until residual (ii) and this
  bridge close (the required bi-tight levels and Lemma D are already closed).

Does NOT set L closed. residual_closed_general=false (16N optional open).
Writes evidence/e1_gmin_m4_prop15168.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15720 import (  # noqa: E402
    bitight_level_obstruction,
    required_bitight_levels_empty_all_primes,
)


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


def freeness_threshold(p: int) -> Fraction:
    """N₁/N bound: if f_e≡1 on U then |U|/N ≤ (p+1)/(2p) (15.42.2)."""
    return Fraction(p + 1, 2 * p)


def deep_s2_freeness_lb(p: int, k: int) -> Fraction:
    """
    For even scores S≥2 on Max+: E[S]≥4−2a with a=N₂/N ⇒ a≥2−k/(2p).
    """
    return Fraction(2) - Fraction(k, 2 * p)


def k_max_auto_freeness_s2(p: int) -> int:
    """Largest k with deep s₊=2 auto-freeness: 3p−2."""
    return 3 * p - 2


def phi_of(p: int) -> Fraction:
    return Fraction(n_of(p) * p, 2)


def tight_cover_size(s: int, p: int) -> int:
    """Size of a Max+-tight level-s cover: s·p (from E[S]=s and |F|=s p)."""
    return s * p


def tight_cover_obstruction_applicable(p: int, s: int) -> dict:
    """
    Corrected compatibility wrapper: only the *bi-tight alternative* is
    excluded. Levels 2 and 3 are the required E(1) alternatives; level 4 is
    a bi-tight corollary only. A generic one-sided Max+-tight cover is not
    claimed empty.
    """
    n = n_of(p)
    E = comb(n, 2)
    size = tight_cover_size(s, p)
    bt = bitight_level_obstruction(p, s)
    # This isotropy identity is valid; only the former kernel conclusion was
    # false. Keep the fields for downstream arithmetic predicates.
    allones_quad = Fraction(size * size, p * p)
    g_perp_quad = Fraction(s * s) - allones_quad
    return {
        "p": p,
        "s": s,
        "size": size,
        "E_edges": E,
        "size_lt_E": size < E,
        "required_level": s in (2, 3),
        "level_4_bitight_corollary": s == 4,
        "bitight_empty_15720": bt["bi_tight_empty"],
        "g_perp_quad_for_tight": str(g_perp_quad),
        "g_perp_vanishes": g_perp_quad == 0,
        "generic_tight_cover_empty": False,
        "obstruction_fires": bool(bt["bi_tight_empty"] and size < E),
        "theorem": (
            "If the relevant alternative is bi-tight at level 2 or 3, 15.720 "
            "excludes it by the degree congruence. The same arithmetic excludes "
            "bi-tight level 4, but not a generic one-sided tight cover."
        ),
    }


def type_I_fail_k_2p_minus_1_ND(p: int) -> dict:
    """
    Type I freeness-fail equality k=2p−1 / S∈{1,3} produces H=G∪e tight
    size 2p (15.43.3). Master lemma (15.44): Φ(H)≥Φ or bi-tight.
    With bi-tight empty (15.720), Φ(H)≥Φ ⇒ ND for that e.

    Checkable part shipped here: bi-tight empty + size-2p obstruction fires.
    """
    obs = tight_cover_obstruction_applicable(p, s=2)
    k = 2 * p - 1
    return {
        "p": p,
        "k": k,
        "prior_structure": "15.43.3 freeness-fail → H tight size 2p; 15.44 master lemma",
        "bitight_empty": obs["obstruction_fires"],
        "tight_2p_obstruction": obs["obstruction_fires"],
        # ND for this class when bi-tight empty (prior structure + 15.720)
        "ND_for_this_class": bool(obs["obstruction_fires"]),
    }


def deep_tight_empty(p: int) -> dict:
    """Deep tight ⇒ bi-tight level 2 (15.44.3) ⇒ empty by 15.720."""
    obs = tight_cover_obstruction_applicable(p, s=2)
    return {
        "p": p,
        "prior": "15.44.3 deep tight ⇒ bi-tight level 2",
        "bitight_empty": obs["obstruction_fires"],
        "empty": bool(obs["obstruction_fires"]),
    }


def deep_s2_auto_freeness(p: int, k: int) -> dict:
    """Fraction: for s₊=2 even scores, k≤3p−2 ⇒ N₂/N lb > freeness thr ⇒ freeness."""
    thr = freeness_threshold(p)
    lb = deep_s2_freeness_lb(p, k)
    return {
        "p": p,
        "k": k,
        "N2_N_lb": str(lb),
        "freeness_thr": str(thr),
        "auto_freeness": lb > thr,
        # freeness ⇒ ND is prior score-eval (15.43 style); auto_freeness is the new algebra
        "implies_freeness_ND_for_all_e": lb > thr,
    }


def deep_fail_k_3p_minus_1_impossible(p: int) -> dict:
    """
    Freeness-fail equality + S∈{2,4} only ⇒ H=G∪e tight S≡3 size 3p.
    Impossible when Thm A fires for s=3.
    """
    obs3 = tight_cover_obstruction_applicable(p, s=3)
    k_eq = 3 * p - 1
    return {
        "p": p,
        "k_equality": k_eq,
        "H_size_if_S_in_2_4": 3 * p,
        "tight_L3_obstruction": obs3["obstruction_fires"],
        "bi_tight_alternative_impossible": obs3["obstruction_fires"],
        "impossible_when_bitight_empty": obs3["obstruction_fires"],
        "generic_tight_cover_impossible": False,
    }


def e1_open_residuals() -> list[str]:
    """Current residual and implication remainders, not the obsolete slice."""
    open_: list[str] = []
    try:
        from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
        from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
        from e1_gmin_m4_prop15276 import (
            lemma_D_2plane_amplitudes_proved,
            lemma_D_existence_written,
        )
        from e1_gmin_m4_prop15764 import minimal_gap4_shell_bridge_closed_general

        if not residual_ii_k_ge_4p_ND_closed():
            open_.append("non-Walsh residual (ii), even k≥4p")
        if not type_I_multilevel_bad_case_ND_closed():
            open_.append("Type I multi-level Max− bad case")
        if not (lemma_D_existence_written() and lemma_D_2plane_amplitudes_proved()):
            open_.append("Lemma D existence/two-plane implication")
        if not minimal_gap4_shell_bridge_closed_general():
            open_.append("minimal four-gap path outside the historical E1 units")
    except Exception:
        return ["current E1 acceptance predicates unavailable"]
    return open_


def e1_closed_general() -> bool:
    """Full corrected E(1) gate; the narrower historical AND is not used."""
    return bool(required_bitight_levels_empty_all_primes() and not e1_open_residuals())


def m_n_ge_phi_minus_2_all_p() -> bool:
    """Same as e1_closed_general — real predicate chain."""
    return e1_closed_general()


def prove_theorem_A(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        for s in (2, 3):
            r = tight_cover_obstruction_applicable(p, s)
            if not r["obstruction_fires"]:
                ok = False
        rows[str(p)] = tight_cover_obstruction_applicable(p, 2)
    return {
        "proved": ok,
        "level_4_bitight_corollary": all(
            tight_cover_obstruction_applicable(p, 4)["obstruction_fires"]
            for p in primes
        ),
        "level_4_one_sided_tight_closed": False,
        "theorem": (
            "For every prime p≥5 the required level-2 and level-3 bi-tight "
            "alternatives are empty by 15.720. Bi-tight level 4 is a corollary; "
            "no generic one-sided tight-cover emptiness is claimed."
        ),
        "n_checked": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
    }


def prove_theorem_B(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {str(p): deep_tight_empty(p) for p in primes}
    ok = all(r["empty"] for r in rows.values())
    return {
        "proved": ok,
        "theorem": "Deep tight ⇒ bi-tight level 2 (15.44.3) ⇒ empty by 15.720.",
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
        "n_checked": len(primes),
    }


def prove_theorem_C_D_type_I(primes: list[int] | None = None) -> dict:
    """
    C: freeness ND is prior 15.43.1 (not re-proved; flagged prior_proved).
    D: k=2p−1 fail class ND when 15.720 makes bi-tight empty (checkable).
    The historical k=3p−2/two-level boundary is closed downstream. This
    proposition does not close its multi-level mechanism; Proposition 15.750
    closes the current global Type-I gate independently.
    """
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {str(p): type_I_fail_k_2p_minus_1_ND(p) for p in primes}
    ok_k2pm1 = all(r["ND_for_this_class"] for r in rows.values())
    try:
        from e1_gmin_m4_prop15170 import type_I_k_3p_minus_2_closed_general as _c170

        k3pm2_closed = bool(_c170())
    except Exception:
        k3pm2_closed = False
    historical_two_level = bool(ok_k2pm1 and k3pm2_closed)
    try:
        from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed

        multilevel_closed = bool(type_I_multilevel_bad_case_ND_closed())
    except Exception:
        multilevel_closed = False
    return {
        "proved_freeness_ND_prior": True,  # 15.43.1 already in solution/repo
        "proved_k_2p_minus_1_fail_ND": ok_k2pm1,
        "k_3p_minus_2_boundary_open": not k3pm2_closed,
        "type_I_historical_two_level_classes_closed": historical_two_level,
        "type_I_closed_by_this_proposition": False,
        "type_I_all_classes_closed": multilevel_closed,
        "type_I_multilevel_bad_case_closed": multilevel_closed,
        "theorem": (
            "Type I freeness ND = 15.43.1 (prior). Type I fail k=2p−1 ND when "
            "the level-2 bi-tight alternative is empty (15.43.3+15.44+15.720). "
            "The historical k=3p−2/two-level slice is closed downstream. "
            "Proposition 15.750 independently closes the multi-level bad case."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
        "n_checked": len(primes),
    }


def prove_theorem_E_F_G_deep(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok_auto = True
    ok_k3pm1 = True
    for p in primes:
        auto = deep_s2_auto_freeness(p, k_max_auto_freeness_s2(p))
        # one past boundary: not auto
        past = deep_s2_auto_freeness(p, 3 * p - 1)
        fail = deep_fail_k_3p_minus_1_impossible(p)
        rows[str(p)] = {"auto_at_3p_minus_2": auto, "at_3p_minus_1": past, "fail_eq": fail}
        if not auto["auto_freeness"]:
            ok_auto = False
        if past["auto_freeness"]:
            ok_auto = False  # must NOT auto-free at k=3p-1
        if not fail["impossible_when_bitight_empty"]:
            ok_k3pm1 = False
    try:
        from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed as _d171

        deep_k_closed = bool(_d171())
    except Exception:
        deep_k_closed = False
    return {
        "proved_auto_freeness_k_le_3p_minus_2": ok_auto,
        "proved_fail_eq_k_3p_minus_1_impossible": ok_k3pm1,
        "deep_freeness_fail_k_ge_3p_open": not deep_k_closed,
        "deep_all_ND_closed": bool(ok_auto and ok_k3pm1 and deep_k_closed),
        "theorem": (
            "Auto-freeness for s₊=2, k≤3p−2 (Fraction). Fail-eq k=3p−1 ⇒ tight L3 "
            "has its bi-tight level-3 alternative empty under 15.720. "
            "Freeness⇒ND prior. The affine and bounded even residual-(ii) range "
            "through 4p−2 is closed downstream; the non-Walsh multi-level "
            "even-k≥4p remainder stays open."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
        "n_checked": len(primes),
    }


def prove_status() -> dict:
    open_res = e1_open_residuals()
    e1 = e1_closed_general()
    bt = required_bitight_levels_empty_all_primes()
    L = main_L_from_e1(e1, bt)
    return {
        "bi_tight_empty_for_all_p_ge_5": bt,
        "E1_closed_general": e1,
        "m_n_ge_Phi_minus_2_all_p": m_n_ge_phi_minus_2_all_p(),
        "L_status": L["L_status"],
        "residual_closed_general": False,  # 16N optional path still open
        "open_residuals": open_res,
        "note": (
            "15.720 closes required bi-tight levels. Historical two-level Type I "
            "and affine/bounded residual-(ii) slices are closed; Proposition 15.750 "
            "closes multi-level Type I. Even-k>=4p residual (ii) remains open. "
            "E1/L use that live predicate; the optional 16N route is also open."
        ),
    }


def main_L_from_e1(e1: bool, bitight: bool) -> dict:
    """L closed only if BOTH e1 and bitight — e1 is currently false."""
    L_closed = bool(e1 and bitight)
    return {
        "bi_tight": bitight,
        "E1": e1,
        "L_closed": L_closed,
        "L_status": "CLOSED" if L_closed else "OPEN",
        "rule": "L closed iff bi-tight ∧ E(1); denseness Prop 6.2",
    }


def e1_residual_open() -> dict:
    return {
        "open": e1_open_residuals(),
        "E1_closed": e1_closed_general(),
    }


def main() -> dict:
    A = prove_theorem_A()
    B = prove_theorem_B()
    CD = prove_theorem_C_D_type_I()
    EFG = prove_theorem_E_F_G_deep()
    st = prove_status()
    bt = required_bitight_levels_empty_all_primes()
    e1 = e1_closed_general()
    Lwire = main_L_from_e1(e1, bt)
    out = {
        "title": (
            "Prop 15.168 E(1) structure after corrected 15.720 bi-tight gate "
            + ("(E1 CLOSED)" if e1 else "(E1 OPEN)")
        ),
        "L_status": Lwire["L_status"],
        "proved": {
            "tight_cover_obstruction_when_bitight_empty": A["proved"],
            "deep_tight_empty_p_ge_5": B["proved"],
            "type_I_freeness_ND_prior": CD["proved_freeness_ND_prior"],
            "type_I_fail_k_2p_minus_1_ND": CD["proved_k_2p_minus_1_fail_ND"],
            "type_I_historical_two_level_classes_closed": CD[
                "type_I_historical_two_level_classes_closed"
            ],
            "type_I_all_classes_closed": CD["type_I_all_classes_closed"],
            "type_I_multilevel_bad_case_closed": CD[
                "type_I_multilevel_bad_case_closed"
            ],
            "deep_auto_freeness_k_le_3p_minus_2": EFG["proved_auto_freeness_k_le_3p_minus_2"],
            "deep_fail_eq_k_3p_minus_1_impossible": EFG["proved_fail_eq_k_3p_minus_1_impossible"],
            "deep_all_ND_closed": EFG["deep_all_ND_closed"],
            "E1_closed_general": e1,
            "m_n_ge_Phi_minus_2": e1,
            "L_closed": Lwire["L_closed"],
            "bi_tight_empty_for_all_p_ge_5": bt,
            "residual_closed_general": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {"A": A, "B": B, "CD": CD, "EFG": EFG, "status": st},
        "L_wire": Lwire,
        "open_residual": st["open_residuals"],
        "F3": "E1/L only from real bi-tight ∧ residual(i) ∧ residual(ii)",
        "F13": "15.168–171 E1 ND structure; residual/16N still open",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15168.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.168 E(1) structure")
    print(f"  A tight obstruction: {A['proved']}")
    print(f"  B deep tight empty: {B['proved']}")
    print(f"  CD type I k=2p-1 fail ND: {CD['proved_k_2p_minus_1_fail_ND']}")
    print(
        "  CD historical Type-I two-level slices closed: "
        f"{CD['type_I_historical_two_level_classes_closed']}"
    )
    print(f"  CD Type-I all classes closed: {CD['type_I_all_classes_closed']}")
    print(f"  EFG auto-freeness: {EFG['proved_auto_freeness_k_le_3p_minus_2']}")
    print(f"  EFG k=3p-1 fail impossible: {EFG['proved_fail_eq_k_3p_minus_1_impossible']}")
    print(f"  EFG deep all ND: {EFG['deep_all_ND_closed']}")
    print(f"  E1_closed_general={e1}")
    print(f"  open: {st['open_residuals']}")
    print(f"  L_status={Lwire['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
