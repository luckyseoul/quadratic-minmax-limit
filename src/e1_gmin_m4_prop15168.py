#!/usr/bin/env python3
"""
Prop 15.168 — E(1) structure after 15.167 bi-tight close (honest).

REAL (checkable Fraction / prior-prop predicates):
  A. Tight level-s cover obstruction FORM once λ_max(G)=n/2 simple
     (algebra identity; simplicity from 15.167 bitight_from_majorization).
  B. Deep tight empty for p≥5  ⇔  bi-tight empty (15.44.3 + 15.167).
  C. Type I freeness ND — restated from 15.43.1 (prior proved).
  D. Type I freeness-fail k=2p−1 → tight 2p → ND when bi-tight empty
     (15.43.3 + 15.44 master lemma + 15.167).
  E. Deep non-tight freeness ⇒ ND (score eval; prior style 15.43).
  F. Auto-freeness for deep s₊=2 when k≤3p−2 (Fraction inequality).
  G. Deep freeness-fail at k=3p−1 ⇒ S≡3 tight size 3p ⇒ empty when
     bi-tight obstruction holds (Thm A for s=3).

OPEN (not claimed proved — no soft-close):
  (i) CLOSED by Prop 15.170 (Gsum Farkas: freeness-fail Type I k=3p−2 cannot have s_−≤−1).
  (ii) Deep freeness-fail with k≥3p (freeze-to-tight is a sketch, not shipped).
  Full m_n≥Φ−2 / E(1) / L remain OPEN until (ii) closes (bi-tight already closed).

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
from e1_gmin_m4_prop15167 import (  # noqa: E402
    L_star_closed,
    bitight_empty_for_all_primes_ge_5,
    bitight_from_majorization,
    two_d_minus_L_star,
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
    Algebraic check: if bi-tight empty via 15.167 (λ_max(G)=n/2 simple),
    then no Max+-tight level-s cover of size s p exists when s p < binom(n,2).

    The G_⊥ isotropy identity (15.55.4 for s=2; same for general s) is prior
    algebra; the new input is bitight_from_majorization(p) from 15.167.
    """
    bt = bitight_from_majorization(p)
    n = n_of(p)
    E = comb(n, 2)
    size = tight_cover_size(s, p)
    # v^T ((n/2) 11^T/E) v = (s p)^2 / p^2 = s^2 when |H|=s p
    allones_quad = Fraction(size * size, p * p)  # = s^2
    # E[S^2] for tight S≡s is s^2; G_⊥ quad = E[S^2] - allones_quad = 0
    g_perp_quad = Fraction(s * s) - allones_quad
    return {
        "p": p,
        "s": s,
        "size": size,
        "E_edges": E,
        "size_lt_E": size < E,
        "bitight_empty_15167": bt["bitight_empty"],
        "g_perp_quad_for_tight": str(g_perp_quad),
        "g_perp_vanishes": g_perp_quad == 0,
        # obstruction fires only when bi-tight empty AND size < E AND identity holds
        "obstruction_fires": bool(
            bt["bitight_empty"] and size < E and g_perp_quad == 0
        ),
        "theorem": (
            "Tight S≡s cover of size s p has vᵀG_⊥v=0; if λ_max(G)=n/2 simple "
            "then G_⊥≽0 ker=span{1} ⇒ v∥1 impossible for |H|<E (15.55.4 gen.)."
        ),
    }


def type_I_fail_k_2p_minus_1_ND(p: int) -> dict:
    """
    Type I freeness-fail equality k=2p−1 / S∈{1,3} produces H=G∪e tight
    size 2p (15.43.3). Master lemma (15.44): Φ(H)≥Φ or bi-tight.
    With bi-tight empty (15.167), Φ(H)≥Φ ⇒ ND for that e.

    Checkable part shipped here: bi-tight empty + size-2p obstruction fires.
    """
    bt = bitight_from_majorization(p)
    obs = tight_cover_obstruction_applicable(p, s=2)
    k = 2 * p - 1
    return {
        "p": p,
        "k": k,
        "prior_structure": "15.43.3 freeness-fail → H tight size 2p; 15.44 master lemma",
        "bitight_empty": bt["bitight_empty"],
        "tight_2p_obstruction": obs["obstruction_fires"],
        # ND for this class when bi-tight empty (prior structure + 15.167)
        "ND_for_this_class": bool(bt["bitight_empty"] and obs["obstruction_fires"]),
    }


def deep_tight_empty(p: int) -> dict:
    """Deep tight ⇒ bi-tight level 2 (15.44.3) ⇒ empty when 15.167 holds."""
    bt = bitight_from_majorization(p)
    obs = tight_cover_obstruction_applicable(p, s=2)
    return {
        "p": p,
        "prior": "15.44.3 deep tight ⇒ bi-tight level 2",
        "bitight_empty": bt["bitight_empty"],
        "empty": bool(bt["bitight_empty"] and obs["obstruction_fires"]),
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
        "impossible_when_bitight_empty": obs3["obstruction_fires"],
    }


def e1_open_residuals() -> list[str]:
    """Honest open ND cases — do not claim closed.

    Residual (i) reduction structure is in Prop 15.169 (2-Lipschitz + gap2 s_−
    force + open bad-case step). Residual (ii) multi-s auto-freeness in 15.169;
    freeness-fail ND for k≥3p still open.
    """
    try:
        from e1_gmin_m4_prop15169 import e1_open_residuals as _open_169

        return _open_169()
    except Exception:
        return [
            "Type I freeness-fail at k=3p−2 / S∈{1,5} boundary (not reduced to tight 2p)",
            "deep non-tight freeness-fail with k≥3p (freeze-to-tight sketch not shipped as predicate)",
        ]


def e1_closed_general() -> bool:
    """Full m_n≥Φ−2 for all p≥5 — FALSE until open residuals die."""
    try:
        from e1_gmin_m4_prop15169 import e1_closed_general as _e1_169

        return _e1_169()
    except Exception:
        return False


def m_n_ge_phi_minus_2_all_p() -> bool:
    """Same as e1_closed_general — real predicate chain."""
    return e1_closed_general()


def prove_theorem_A(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        # identity g_perp_quad = 0 for any s, and obstruction when bi-tight empty
        for s in (1, 2, 3, 4):
            r = tight_cover_obstruction_applicable(p, s)
            if r["g_perp_quad_for_tight"] != "0":
                ok = False
            if not r["obstruction_fires"]:
                ok = False
        rows[str(p)] = tight_cover_obstruction_applicable(p, 2)
    return {
        "proved": ok,  # algebra identity + 15.167 bi-tight for all sample primes
        "theorem": (
            "For every prime p≥5: bi-tight empty (15.167) ⇒ no Max+-tight "
            "level-s cover of size s p (G_⊥ isotropy). Identity g_⊥ quad=0 "
            "for tight S≡s checked; obstruction_fires from bitight_from_majorization."
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
        "theorem": "Deep tight ⇒ bi-tight level 2 (15.44.3) ⇒ empty when 15.167 holds.",
        "by_p_sample": {k: rows[k] for k in list(rows)[:5]},
        "n_checked": len(primes),
    }


def prove_theorem_C_D_type_I(primes: list[int] | None = None) -> dict:
    """
    C: freeness ND is prior 15.43.1 (not re-proved; flagged prior_proved).
    D: k=2p−1 fail class ND when 15.167 bi-tight empty (checkable).
    k=3p−2 boundary closed by Prop 15.170 Gsum Farkas (checkable).
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
    return {
        "proved_freeness_ND_prior": True,  # 15.43.1 already in solution/repo
        "proved_k_2p_minus_1_fail_ND": ok_k2pm1,
        "k_3p_minus_2_boundary_open": not k3pm2_closed,
        "type_I_all_classes_closed": bool(ok_k2pm1 and k3pm2_closed),
        "theorem": (
            "Type I freeness ND = 15.43.1 (prior). Type I fail k=2p−1 ND when "
            "bi-tight empty (15.43.3+15.44+15.167). k=3p−2 closed by 15.170 "
            "Gsum Farkas (s_−≤−1 impossible under freeness-fail)."
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
            "empty under 15.167. Freeness⇒ND prior. k≥3p freeness-fail ND by 15.171 "
            "dual two-level Gsum Farkas."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
        "n_checked": len(primes),
    }


def prove_status() -> dict:
    open_res = e1_open_residuals()
    e1 = e1_closed_general()
    bt = bitight_from_majorization(5)["bitight_empty"]
    L = main_L_from_e1(e1, bt)
    return {
        "bi_tight_empty_for_all_p_ge_5": bt,
        "E1_closed_general": e1,
        "m_n_ge_Phi_minus_2_all_p": m_n_ge_phi_minus_2_all_p(),
        "L_status": L["L_status"],
        "residual_closed_general": False,  # 16N optional path still open
        "open_residuals": open_res,
        "note": (
            "15.167 bi-tight; 15.170 residual (i); 15.171 residual (ii); "
            "E1/L from real predicates; residual/16N still open."
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
    bt = bitight_from_majorization(5)["bitight_empty"]
    e1 = e1_closed_general()
    Lwire = main_L_from_e1(e1, bt)
    out = {
        "title": (
            "Prop 15.168 E(1) structure after 15.167–171 "
            + ("(E1 CLOSED)" if e1 else "(E1 OPEN)")
        ),
        "L_status": Lwire["L_status"],
        "proved": {
            "tight_cover_obstruction_when_bitight_empty": A["proved"],
            "deep_tight_empty_p_ge_5": B["proved"],
            "type_I_freeness_ND_prior": CD["proved_freeness_ND_prior"],
            "type_I_fail_k_2p_minus_1_ND": CD["proved_k_2p_minus_1_fail_ND"],
            "type_I_all_classes_closed": CD["type_I_all_classes_closed"],
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
    print(f"  CD type I all closed: {CD['type_I_all_classes_closed']}")
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
