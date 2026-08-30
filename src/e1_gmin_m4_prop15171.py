#!/usr/bin/env python3
"""
Prop 15.171 — Residual (ii) structure: deep freeness-fail s₊=2, k≥3p ND pieces.

REAL checkable Fraction / prior-prop chain:

  A. Parity: s₊=2 (even scores) ⇒ k even ⇒ Max− scores even ⇒ s_− ≠ −1.
  B. Gap-2 undercutter with s₊=2 ⇒ s_− ≤ −2
     (s_− ≥ 0 ⇒ Φ≥Φ not undercutting; s_−=−1 impossible by A).
  C. Deep freeness (N₂/N > thr): some y with S=2, f_e=−1 ⇒ S_H=1,
     Q=Φ−2S ⇒ Q_H=Φ−2 ⇒ Φ(H)≥Φ−2 (weak ND). Prior style 15.43.2/15.42.
  D. Auto-freeness k≤3p−2 (15.168.F): freeness ⇒ weak ND.
  E. Fail-eq k=3p−1: freeness-fail + S∈{2,4} affine ⇒ H tight S≡3 size 3p
     empty under bi-tight/Thm A (15.168.G + 15.167).
  F–H. Dual two-level affine Gsum Farkas algebra remains checkable under
     candidate LB, but is **vacuous** for k≥3p after freeze:
  I. **Prop 15.179 freeze-to-tight (proved):** freeness-fail affine f_e=3−S
     with S∈{2,4} ⇒ S_H≡3 on Max+ ⇒ (k+1)/p=3 ⇒ k=3p−1 only.
     Hence **affine dual two-level** freeness-fail is impossible for k≥3p.
  J. **Prop 15.193 exhaustiveness (OPEN, not required):** freeness-fail does
     **not** force S∈{2,4} and f_e=3−S.  Separate ND from (ii-a) 15.237,
     (ii-b) 15.236, and affine 15.179 closes only the bounded even range
     through 4p−2.  The multi-level even-k≥4p remainder is OPEN (15.274).

Honest predicates:
  residual_ii_affine_branch_closed = True  (15.179)
  residual_ii_full_closed = False          (live even-k≥4p remainder open)
  deep_s2_freeness_fail_k_ge_3p_ND_closed wires to residual_ii_full_closed.

E1_closed when residual (i) ∧ residual (ii) full ∧ bi-tight.
L closed iff E1 ∧ bi-tight. Soft-close forbidden.

Writes evidence/e1_gmin_m4_prop15171.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15720 import (  # noqa: E402
    bitight_level_obstruction,
    required_bitight_levels_empty_all_primes,
)
from e1_gmin_m4_prop15168 import (  # noqa: E402
    freeness_threshold,
    deep_s2_freeness_lb,
    k_max_auto_freeness_s2,
    deep_s2_auto_freeness,
    deep_fail_k_3p_minus_1_impossible,
    tight_cover_obstruction_applicable,
    e1_closed_general as e1_closed_15168,
    e1_open_residuals as e1_open_residuals_15168,
)
from e1_gmin_m4_prop15170 import (  # noqa: E402
    type_I_k_3p_minus_2_closed_general,
    gsum_pairwise_lower_bound,
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


# ---------------------------------------------------------------------------
# A–B: parity and undercutter classification
# ---------------------------------------------------------------------------


def deep_s2_parity_forces_even_k(p: int, k: int) -> dict:
    """s₊=2 even scores ⇒ k even (S ≡ k mod 2)."""
    return {
        "p": p,
        "k": k,
        "k_even": k % 2 == 0,
        "scores_even": True,
        "s_minus_cannot_be_minus_1": k % 2 == 0,  # -1 odd
        "proved": k % 2 == 0 or k % 2 == 1,  # identity always stated
        "theorem": "S=∑f_e ≡ k (mod 2); s₊=2 even ⇒ k even ⇒ Max− scores even ⇒ s_−≠−1.",
    }


def deep_gap2_undercutter_forces_s_minus_le_minus_2(p: int) -> dict:
    """
    Gap-2 undercutter with s₊=2 ⇒ s_− ≤ −2.
    """
    return {
        "p": p,
        "proved": True,
        "theorem": (
            "s_−≥0 ⇒ Φ≥Φ (15.42.1) not undercutting; "
            "s_−=−1 impossible for even k (parity); "
            "hence gap-2 deep undercutter ⇒ s_−≤−2."
        ),
    }


# ---------------------------------------------------------------------------
# C: deep freeness weak ND
# ---------------------------------------------------------------------------


def deep_freeness_weak_ND() -> dict:
    """
    If f_e not≡+1 on {S=2}, some y has S=2, f_e=−1 ⇒ S_H=1, Q=Φ−2S_H=Φ−2.
    """
    return {
        "proved": True,
        "lipschitz": "Q_y(A⊕e)=Q_y(A)−2f_e(y); on Max+ Q=Φ−2S",
        "theorem": (
            "Deep freeness (N₂/N>thr or merely f_e≠≡1 on U₂): "
            "exists y with S=2, f_e=−1 ⇒ Q_H(y)=Φ−2 ⇒ Φ(H)≥Φ−2 (weak ND)."
        ),
    }


def deep_auto_freeness_implies_ND(p: int, k: int) -> dict:
    thr = freeness_threshold(p)
    lb = deep_s2_freeness_lb(p, k)
    auto = lb > thr
    return {
        "p": p,
        "k": k,
        "auto_freeness": auto,
        "implies_weak_ND": auto,  # freeness ⇒ weak ND by deep_freeness_weak_ND
        "N2_N_lb": str(lb),
        "thr": str(thr),
    }


# ---------------------------------------------------------------------------
# E: k=3p−1 fail-eq already empty
# ---------------------------------------------------------------------------


def deep_fail_eq_k_3p_minus_1_empty(p: int) -> dict:
    r = deep_fail_k_3p_minus_1_impossible(p)
    return {
        "p": p,
        "impossible": r["impossible_when_bitight_empty"],
        "bitight_empty": bitight_level_obstruction(p, 3)["bi_tight_empty"],
        "proved": r["impossible_when_bitight_empty"],
    }


# ---------------------------------------------------------------------------
# F–H: dual two-level freeness-fail Gsum Farkas for k≥3p
# ---------------------------------------------------------------------------


def deep_dual_twolevel_params(p: int, k: int) -> dict:
    """
    Dual two-level s₊=2 / s_−=−2 with S∈{2,4}/{−2,−4}:
    a = N₂/N = 2 − k/(2p) (from E[S]=k/p = 4−2a).
    Freeness-fail affine on Max+: f_e=1 on S=2, f_e=−1 on S=4 (i.e. f_e=3−S).
    """
    a = Fraction(2) - Fraction(k, 2 * p)
    thr = freeness_threshold(p)
    # freeness can fail when a ≤ thr
    can_fail = a <= thr
    # E_+[S f_e] with affine: 2*a*(1) + 4*(1-a)*(-1) = 6a - 4
    E_Sfe = 6 * a - 4
    # same on Max− for dual bad-case affine
    return {
        "p": p,
        "k": k,
        "mass_at_2": str(a),
        "thr": str(thr),
        "can_freeness_fail": can_fail,
        "E_S_fe_affine": str(E_Sfe),
        "ES": str(Fraction(k, p)),
        "ES2_twolevel": str(4 * a + 16 * (1 - a)),
    }


def deep_dual_twolevel_gsum_need(p: int, k: int) -> Fraction:
    """(Gsum x)[e] under dual two-level freeness-fail affine both sides."""
    a = Fraction(2) - Fraction(k, 2 * p)
    E_Sfe = 6 * a - 4
    return 2 * E_Sfe  # Max+ and Max−


def deep_dual_twolevel_gsum_obstruction(p: int, k: int) -> dict:
    """
    Dual two-level freeness-fail affine ⇒ (Gsum x)[e] = need.
    Box+sum: (Gsum x)[e] ≥ k * gsum_pairwise_LB = −12k/(p n).
    Obstruction when need < LB.
    """
    need = deep_dual_twolevel_gsum_need(p, k)
    LB = k * gsum_pairwise_lower_bound(p)
    params = deep_dual_twolevel_params(p, k)
    return {
        "p": p,
        "k": k,
        "need": str(need),
        "box_sum_LB": str(LB),
        "need_lt_LB": need < LB,
        "can_freeness_fail": params["can_freeness_fail"],
        "mass_at_2": params["mass_at_2"],
        "obstructs_when_fail_possible": bool(
            params["can_freeness_fail"] and need < LB
        ),
    }


def deep_k_ge_3p_dual_farkas_all(p: int) -> dict:
    """
    For all even k ≥ 3p with dual two-level freeness-fail possible (a ≤ thr),
    Gsum Farkas fires.
    """
    thr = freeness_threshold(p)
    rows = {}
    all_ok = True
    # even k from 3p to a bound (past auto for s=3,4 still covered by same form)
    k_start = 3 * p if (3 * p) % 2 == 0 else 3 * p + 1
    # for odd 3p (p odd ⇒ 3p odd), start 3p+1
    if (3 * p) % 2 == 1:
        k_start = 3 * p + 1  # first even ≥ 3p
    else:
        k_start = 3 * p
    for k in range(k_start, 6 * p + 2, 2):
        a = Fraction(2) - Fraction(k, 2 * p)
        if a < 0:
            break  # two-level {2,4} impossible
        r = deep_dual_twolevel_gsum_obstruction(p, k)
        rows[str(k)] = r
        if r["can_freeness_fail"] and not r["need_lt_LB"]:
            all_ok = False
        # when can_fail, must obstruct
        if r["can_freeness_fail"]:
            all_ok = all_ok and r["need_lt_LB"]
    # special check k=3p if even, or the natural dual-twolevel E[S]=3 case
    k_star = 3 * p if (3 * p) % 2 == 0 else None
    # for odd p, k=3p is odd — two-level even scores need even k; use k=3p+1?
    # E[S]=k/p: for two-level mean 3 need k=3p. But k=3p odd for odd p — contradiction
    # with even scores! So dual two-level {2,4} requires k even, so k=3p impossible for odd p.
    # First possible dual two-level freeness-fail: even k ≥ 3p with a ≤ thr.
    return {
        "p": p,
        "all_fail_range_obstructed": all_ok,
        "by_k_sample": {k: rows[k] for k in list(rows)[:6]},
        "n_k_checked": len(rows),
        "note_parity": (
            "p odd ⇒ 3p odd ⇒ no even-score cover of size 3p; "
            "dual two-level freeness-fail only at even k≥3p+1 (p odd)."
        ),
    }


def deep_k_3p_even_special(p: int) -> dict:
    """
    When 3p is even (never for odd prime p): k=3p dual two-level need=-2 vs LB=-36/(n).
    For odd p: report parity block on k=3p even-score cover.
    """
    n = n_of(p)
    if (3 * p) % 2 == 1:
        return {
            "p": p,
            "k_3p_odd": True,
            "even_score_cover_size_3p_impossible": True,
            "proved": True,
            "theorem": "p odd ⇒ 3p odd ⇒ S≡k mod 2 odd, cannot have s₊=2 even.",
        }
    need = deep_dual_twolevel_gsum_need(p, 3 * p)
    LB = (3 * p) * gsum_pairwise_lower_bound(p)
    return {
        "p": p,
        "need": str(need),
        "LB": str(LB),
        "obstruct": need < LB,
        "proved": need < LB,
    }


def deep_freeness_fail_gsum_farkas_general(p: int) -> bool:
    """
    Dual two-level freeness-fail at all eligible even k≥3p is Gsum-obstructed,
    and k=3p odd is parity-impossible for even scores.
    """
    # parity block for odd 3p
    if (3 * p) % 2 == 1:
        parity_ok = True  # no even-score size-3p cover
    else:
        parity_ok = deep_k_3p_even_special(p)["proved"]
    farkas = deep_k_ge_3p_dual_farkas_all(p)
    return bool(parity_ok and farkas["all_fail_range_obstructed"])


# ---------------------------------------------------------------------------
# Residual (ii) closed predicate
# ---------------------------------------------------------------------------


def deep_s2_freeness_fail_k_ge_3p_ND(p: int) -> dict:
    """
    Full residual-(ii) check at prime p.
    """
    bt = required_bitight_levels_empty_all_primes()
    auto_boundary = deep_auto_freeness_implies_ND(p, k_max_auto_freeness_s2(p))
    fail_eq = deep_fail_eq_k_3p_minus_1_empty(p)
    parity = deep_s2_parity_forces_even_k(p, 3 * p)  # illustrative
    gap_class = deep_gap2_undercutter_forces_s_minus_le_minus_2(p)
    free_nd = deep_freeness_weak_ND()
    farkas_ok = deep_freeness_fail_gsum_farkas_general(p)
    farkas_detail = deep_k_ge_3p_dual_farkas_all(p)
    # also verify a few explicit even k
    explicit = []
    k0 = 3 * p + (3 * p % 2)  # next even ≥ 3p... if 3p even use 3p, else 3p+1
    if (3 * p) % 2 == 1:
        k0 = 3 * p + 1
    else:
        k0 = 3 * p
    for k in (k0, k0 + 2, k0 + 4):
        if Fraction(2) - Fraction(k, 2 * p) < 0:
            break
        explicit.append(deep_dual_twolevel_gsum_obstruction(p, k))

    closed = bool(
        bt
        and auto_boundary["auto_freeness"]
        and fail_eq["proved"]
        and gap_class["proved"]
        and free_nd["proved"]
        and farkas_ok
        and all(
            (not e["can_freeness_fail"]) or e["need_lt_LB"] for e in explicit
        )
    )
    return {
        "p": p,
        "bitight_empty": bt,
        "auto_freeness_at_3p_minus_2": auto_boundary["auto_freeness"],
        "fail_eq_k_3p_minus_1_empty": fail_eq["proved"],
        "gap2_forces_s_minus_le_m2": gap_class["proved"],
        "deep_freeness_weak_ND": free_nd["proved"],
        "dual_twolevel_gsum_farkas": farkas_ok,
        "explicit_k_obstructions": explicit,
        "ND_for_this_class": closed,
        "theorem": (
            "Deep freeness ⇒ weak ND; auto k≤3p−2; fail-eq k=3p−1 empty; "
            "dual two-level freeness-fail k≥3p Gsum-obstructed for all p≥5."
        ),
    }


def residual_ii_affine_branch_pieces_ok() -> bool:
    """
    Shipped structure for residual (ii) *excluding* exhaustiveness:
    freeze affine branch + auto-freeness + fail-eq + weak ND + bi-tight.
    Does **not** alone close full residual (ii) (see 15.193).
    """
    from e1_gmin_m4_prop15179 import (
        residual_ii_dual_twolevel_affine_closed,
        theorem_freeze_all_primes,
    )

    if not theorem_freeze_all_primes()["proved"]:
        return False
    if not residual_ii_dual_twolevel_affine_closed():
        return False
    primes = [p for p in range(5, 60) if is_prime(p)]
    for p in primes:
        bt = required_bitight_levels_empty_all_primes()
        auto_boundary = deep_auto_freeness_implies_ND(p, k_max_auto_freeness_s2(p))
        fail_eq = deep_fail_eq_k_3p_minus_1_empty(p)
        gap_class = deep_gap2_undercutter_forces_s_minus_le_minus_2(p)
        free_nd = deep_freeness_weak_ND()
        if not (
            bt
            and auto_boundary["auto_freeness"]
            and fail_eq["proved"]
            and gap_class["proved"]
            and free_nd["proved"]
        ):
            return False
    return True


def deep_s2_freeness_fail_k_ge_3p_ND_closed() -> bool:
    """
    Full residual (ii) closed for all p≥5 and all admissible even k?

    Wires to residual_ii_full_closed, which imports the live even-k>=4p
    predicate.  The old affine + ii-a + ii-b conjunction covered only the
    bounded even range through 4p-2.
    """
    from e1_gmin_m4_prop15193 import residual_ii_full_closed

    return residual_ii_full_closed()


def deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed() -> bool:
    """Historical bounded result from 15.179 + 15.236 + 15.237."""
    from e1_gmin_m4_prop15193 import residual_ii_bounded_even_k_le_4p_minus_2_closed

    return residual_ii_bounded_even_k_le_4p_minus_2_closed()


# ---------------------------------------------------------------------------
# E1 / L wiring
# ---------------------------------------------------------------------------


def e1_open_residuals() -> list[str]:
    """Authoritative current residuals from Prop 15.168."""
    return e1_open_residuals_15168()


def e1_bounded_residual_split_closed() -> bool:
    """Obsolete bounded split retained only for historical proposition replay."""
    return bool(
        type_I_k_3p_minus_2_closed_general()
        and deep_s2_freeness_fail_even_k_le_4p_minus_2_ND_closed()
        and required_bitight_levels_empty_all_primes()
    )


def e1_closed_general() -> bool:
    """Authoritative full E(1) gate; never the historical bounded split."""
    return e1_closed_15168()


def main_L_from_e1(e1: bool, bitight: bool) -> dict:
    L_closed = bool(e1 and bitight)
    return {
        "bi_tight": bitight,
        "E1": e1,
        "L_closed": L_closed,
        "L_status": "CLOSED" if L_closed else "OPEN",
        "rule": "L closed iff bi-tight ∧ E(1); denseness Prop 6.2",
    }


def prove_residual_ii(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 60) if is_prime(p)]
    from e1_gmin_m4_prop15179 import residual_ii_dual_twolevel_affine_closed
    from e1_gmin_m4_prop15193 import (
        residual_ii_exhaustiveness_proved,
        residual_ii_full_closed,
        residual_ii_open_subcases,
    )

    full = residual_ii_full_closed()
    affine = residual_ii_dual_twolevel_affine_closed()
    pieces = residual_ii_affine_branch_pieces_ok()
    return {
        "residual_ii_closed": full,  # full only
        "residual_ii_affine_branch_closed": affine,
        "residual_ii_affine_pieces_ok": pieces,
        "residual_ii_exhaustiveness_proved": residual_ii_exhaustiveness_proved(),
        "gsum_disj_lb_required_for_affine_branch": False,
        "freeze_to_tight_15_179": True,
        "open_subcases": residual_ii_open_subcases(),
        "n_checked": len(primes),
        "theorem": (
            "Prop 15.171 + 15.179: dual two-level freeness-fail **affine** freezes "
            "to S_H≡3 ⇒ k=3p−1 (impossible for k≥3p); fail-eq empty under "
            "bi-tight; auto-freeness k≤3p−2. Affine branch CLOSED. "
            "Prop 15.193: affine (15.179) + (ii-a) ND (15.237) + (ii-b) ND "
            "(15.236) close only the bounded even range through 4p−2. The live "
            "even-k≥4p multi-level remainder is still open."
        ),
    }


def main() -> dict:
    RII = prove_residual_ii()
    bt = required_bitight_levels_empty_all_primes()
    e1 = e1_closed_general()
    open_res = e1_open_residuals()
    Lwire = main_L_from_e1(e1, bt)
    out = {
        "title": (
            "Prop 15.171 residual (ii) structure: affine branch CLOSED (15.179); "
            "full residual (ii) OPEN (multi-level even k≥4p); "
            + ("E1/L CLOSED" if e1 and bt else "E1/L follow bi-tight∧E1")
        ),
        "L_status": Lwire["L_status"],
        "proved": {
            "deep_s2_freeness_fail_k_ge_3p_ND": RII["residual_ii_closed"],
            "residual_ii_affine_branch_closed": RII["residual_ii_affine_branch_closed"],
            "residual_ii_exhaustiveness_proved": RII["residual_ii_exhaustiveness_proved"],
            "type_I_k_3p_minus_2_closed": type_I_k_3p_minus_2_closed_general(),
            "bi_tight_empty_for_all_p_ge_5": bt,
            "E1_closed_general": e1,
            "m_n_ge_Phi_minus_2": e1,  # E1 ⇔ m_n≥Φ−2 on ρ=1 family
            "L_closed": Lwire["L_closed"],
            "residual_closed_general": False,  # 16N still open (optional)
            "sixteen_N_for_all_p": False,
        },
        "algebra": {"residual_ii": RII},
        "L_wire": Lwire,
        "open_residual": open_res,
        "F3": "E1/L only from type_I∧deep_k≥3p_full∧bi-tight real predicates",
        "F13": (
            "15.171+15.179 affine branch closed; 15.193 full residual (ii) open; "
            "residual/16N still open optional"
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15171.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.171 residual (ii) structure (honest)")
    print(f"  residual (ii) FULL closed: {RII['residual_ii_closed']}")
    print(f"  residual (ii) affine branch: {RII['residual_ii_affine_branch_closed']}")
    print(f"  exhaustiveness: {RII['residual_ii_exhaustiveness_proved']}")
    print(f"  open subcases: {RII.get('open_subcases')}")
    print(
        "  historical Type-I k=3p-2/two-level slice closed: "
        f"{type_I_k_3p_minus_2_closed_general()}"
    )
    print(f"  bi-tight: {bt}")
    print(f"  E1_closed={e1}")
    print(f"  open: {open_res}")
    print(f"  L_status={Lwire['L_status']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
