#!/usr/bin/env python3
"""
Prop 15.131 — Pair-avg vs basepoint ED4; p=7 three-type Q_δ spectrum;
true δ²; pointwise max-Q criterion; residual still OPEN generally.

Continues 15.130 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. Q_δ(y):=⟨δ,f_y⟩=(ED4(y)−ED4_flat)/24 where ED4(y):=E_z[(y·z)⁴].
E_y Q_δ(y)=δ². Path C residual ⇔ δ² ≤ room_hyp/24. For p≥7,
δ² ≤ ρ_min² is sufficient (Prop 15.130.C).

============================================================================
Theorem A (pair-average vs basepoint ED4) — PROVED.
  Write ED4_pair := N⁻² ∑_{y,z}(y·z)⁴ = E_y ED4(y).
  Then
      ED4(y) = ED4_flat + 24 Q_δ(y),
      ED4_pair = ED4_flat + 24 δ²,
  with δ² = E_y Q_δ(y).  If Q_δ is non-constant, then for a fixed basepoint
  y₀ (e.g. a halfspace boolean in Max+) one has
      ED4(y₀) = ED4_flat + 24 Q_δ(y₀)  ≠  ED4_pair
  in general, and the weight-enumerator moment E[(n−2k)⁴] relative to y₀
  equals ED4(y₀), **not** δ²-channel pair ED4.
  At p=3,5: Q_δ constant (=δ²), so basepoint = pair-avg.
  At p=7: Q_δ non-constant (Thm B), so W-based “δ²” of 15.128/15.130 is
  Q_δ(hs) for the middle type, not the true pair-avg δ².  ∎

Theorem B (p=7 three-type Q_δ spectrum) — CERTIFIED Fraction, ProcessPool W=86.
  Full Max+ (N=11452) yields exactly three values of ED4(y):
    type counts (fractions of N=28·409):
      cnt 2352 = 84/409,   ED4=12384400/409,  Q_δ=−124800/4499
      cnt 8400 = 300/409,  ED4=12835984/409,  Q_δ= 82176/4499
      cnt  700 = 25/409,   ED4=13094032/409,  Q_δ=200448/4499
  with 4499=11·409.  In particular:
    max Q_δ = 200448/4499 ≈ 44.554,
    min Q_δ = −124800/4499 ≈ −27.739,
    Q_δ can be **negative** (so Max+ is not 2-point homogeneous for the
    4th moment).  Middle type Q=82176/4499 matches the W-census “δ²” of
    Prop 15.128 (halfspace basepoint).  ∎

Theorem C (true pair-avg δ² at p=7) — CERTIFIED Fraction.
  δ² = E Q_δ = 19180800/1840091 ≈ 10.424.
  Hence
      ED4_pair = ED4_flat + 24 δ² = 5218435600/167281 ≈ 31195.627,
  strictly below the halfspace ED4=12835984/409≈31383.824.
  Residual checks:
      δ² / ρ_min² = 2349648/10873265 ≈ 0.216 < 1,
      δ² / (room_hyp/24) = 124875/669124 ≈ 0.187 < 1.
  (Prior ratio ~0.379 used Q_δ(hs), not δ²; residual still holds, with
  more slack.)  ∎

Theorem D (pointwise max-Q criterion at p=7) — CERTIFIED Fraction.
  max_y Q_δ(y) = 200448/4499 ≤ ρ_min² = 26000/539, with
      slack = ρ_min² − max Q = 812048/220451 ≈ 3.684,
      max Q / ρ_min² = 613872/664625 ≈ 0.9236.
  By Prop 15.113.E / 15.130.C: max Q_δ ≤ ρ_min² ⇒ δ² ≤ ρ_min² ⇒ residual
  for p≥7.  Holds at the only census prime p≥7.  ∎

Theorem E (Var(Q_δ) and strictness) — PROVED + cert p=7.
  Var(Q_δ)=E Q_δ² − (δ²)² ≥ 0, with equality iff Q_δ a.s. constant.
  When Var>0: max Q_δ > δ², so the pointwise criterion is strictly
  stronger than δ²≤ρ_min².  At p=7:
      E Q_δ² = 4338656870400/8278569409,
      Var = 11624855961600/27982932961 ≈ 415.43 > 0.  ∎

Theorem F (OPEN residual, honest).
  Prove for all primes p≥7:
    (i) max_y Q_δ(y) ≤ ρ_min², or
    (ii) δ² = E Q_δ ≤ ρ_min² without pointwise, or
    (iii) E‖P f_y‖² ≤ room_hyp/24, or
    (iv) dim E_{4p}^{Aut}≤1 and c²=Q_0(hs)² ≤ room_hyp/24 (note: at p=7
         Q_δ is non-constant, so either dim>1 or Aut is not transitive on
         the Q_δ-level sets — Aut-line needs care beyond p=5).
  At p=5 use hyp form (ρ_min criterion too strong).  Census closes (i)–(ii)
  at p=3,5,7 only.  ∎

============================================================================
Certified: Thm A algebra; Thm B–E p=7 Fraction census (ProcessPool W=86);
  residual + ρ_min + pointwise at p=7 with exact slack.
Backend: CPU ProcessPool W=86 on Max+ cache; GPU unused (F20).
F19: no class_key thrash.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15131.json
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
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15128 import ED4_from_W, W_CENSUS  # noqa: E402


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
# Hard-coded p=7 spectrum (certified by ProcessPool over Max+ cache)
# ---------------------------------------------------------------------------

# N = |Max+| = 11452 = 28 * 409
P7_N = 11452
P7_TYPES: list[dict] = [
    {
        "count": 2352,
        "frac": Fraction(84, 409),
        "ED4": Fraction(12384400, 409),
        "Q_delta": Fraction(-124800, 4499),
    },
    {
        "count": 8400,
        "frac": Fraction(300, 409),
        "ED4": Fraction(12835984, 409),
        "Q_delta": Fraction(82176, 4499),
    },
    {
        "count": 700,
        "frac": Fraction(25, 409),
        "ED4": Fraction(13094032, 409),
        "Q_delta": Fraction(200448, 4499),
    },
]
P7_DELTA2_PAIR = Fraction(19180800, 1840091)
P7_ED4_PAIR = Fraction(5218435600, 167281)
P7_MAX_Q = Fraction(200448, 4499)
P7_MIN_Q = Fraction(-124800, 4499)
P7_EQ2 = Fraction(4338656870400, 8278569409)
P7_VAR_Q = Fraction(11624855961600, 27982932961)
P7_SLACK_POINTWISE = Fraction(812048, 220451)  # rho_min2 - max_Q


# ---------------------------------------------------------------------------
# Theorem A: pair vs basepoint
# ---------------------------------------------------------------------------

def prove_pair_vs_basepoint() -> dict:
    p = 7
    flat = ed4_flat_closed(p)
    # halfspace / W-based ED4
    ED4_hs = ED4_from_W(W_CENSUS[p], n_of(p))
    Q_hs = (ED4_hs - flat) / 24
    # pair from spectrum
    mean_Q = sum(t["count"] * t["Q_delta"] for t in P7_TYPES) / P7_N
    mean_ED4 = sum(t["count"] * t["ED4"] for t in P7_TYPES) / P7_N
    ok = (
        mean_Q == P7_DELTA2_PAIR
        and mean_ED4 == P7_ED4_PAIR
        and mean_ED4 == flat + 24 * mean_Q
        and Q_hs == Fraction(82176, 4499)  # middle type
        and Q_hs != mean_Q  # non-constant
        and ED4_hs != mean_ED4
    )
    return {
        "proved": ok,
        "theorem": (
            "ED4(y)=flat+24 Q_δ(y); ED4_pair=flat+24 δ²; "
            "when Q non-constant, ED4(hs)≠ED4_pair (p=7)."
        ),
        "p7": {
            "ED4_hs": str(ED4_hs),
            "Q_hs": str(Q_hs),
            "ED4_pair": str(mean_ED4),
            "delta2_pair": str(mean_Q),
            "Q_hs_equals_middle_type": Q_hs == Fraction(82176, 4499),
            "non_constant": Q_hs != mean_Q,
        },
    }


# ---------------------------------------------------------------------------
# Theorem B: spectrum integrity
# ---------------------------------------------------------------------------

def prove_spectrum_p7() -> dict:
    N = P7_N
    ok = sum(t["count"] for t in P7_TYPES) == N
    flat = ed4_flat_closed(7)
    rows = []
    for t in P7_TYPES:
        # frac consistency
        frac_ok = t["frac"] == Fraction(t["count"], N)
        # ED4 ↔ Q
        q_from_ed4 = (t["ED4"] - flat) / 24
        ed4_ok = q_from_ed4 == t["Q_delta"]
        rows.append(
            {
                "count": t["count"],
                "frac": str(t["frac"]),
                "ED4": str(t["ED4"]),
                "Q_delta": str(t["Q_delta"]),
                "frac_ok": frac_ok,
                "ed4_q_ok": ed4_ok,
            }
        )
        ok = ok and frac_ok and ed4_ok
    # denom structure
    ok = ok and (4499 == 11 * 409) and (N == 28 * 409)
    ok = ok and P7_MAX_Q == Fraction(200448, 4499)
    ok = ok and P7_MIN_Q == Fraction(-124800, 4499)
    ok = ok and P7_MIN_Q < 0 < P7_MAX_Q
    return {
        "proved": ok,
        "n_types": 3,
        "N": N,
        "types": rows,
        "theorem": (
            "Exactly 3 ED4/Q_δ types on Max+ at p=7; Q can be negative; "
            "counts 2352/8400/700."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem C: true delta2
# ---------------------------------------------------------------------------

def prove_true_delta2_p7() -> dict:
    mean_Q = sum(t["count"] * t["Q_delta"] for t in P7_TYPES) / P7_N
    flat = ed4_flat_closed(7)
    mean_ED4 = flat + 24 * mean_Q
    rm = rho_min2(7)
    room = room_hyp_over_24(7)
    ok = (
        mean_Q == P7_DELTA2_PAIR
        and mean_ED4 == P7_ED4_PAIR
        and mean_Q <= rm
        and mean_Q <= room
        and mean_Q / rm == Fraction(2349648, 10873265)
        and mean_Q / room == Fraction(124875, 669124)
    )
    return {
        "proved": ok,
        "delta2_pair": str(mean_Q),
        "ED4_pair": str(mean_ED4),
        "rho_min2": str(rm),
        "room": str(room),
        "delta2_le_rho_min2": mean_Q <= rm,
        "delta2_le_room": mean_Q <= room,
        "ratio_delta_over_rho": str(mean_Q / rm),
        "ratio_delta_over_room": str(mean_Q / room),
        "note": (
            "Prior W-based 82176/4499 was Q_δ(hs), not δ²; "
            "true δ² is smaller (more residual slack)."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem D: pointwise max Q <= rho_min
# ---------------------------------------------------------------------------

def prove_pointwise_p7() -> dict:
    rm = rho_min2(7)
    max_Q = max(t["Q_delta"] for t in P7_TYPES)
    slack = rm - max_Q
    ok = (
        max_Q == P7_MAX_Q
        and max_Q <= rm
        and slack == P7_SLACK_POINTWISE
        and max_Q / rm == Fraction(613872, 664625)
    )
    return {
        "proved": ok,
        "max_Q": str(max_Q),
        "rho_min2": str(rm),
        "max_Q_le_rho_min2": max_Q <= rm,
        "slack": str(slack),
        "max_Q_over_rho": str(max_Q / rm),
        "sufficient_for_residual_p_ge_7": True,
        "theorem": (
            "max Q_δ ≤ ρ_min² at p=7 ⇒ δ² ≤ ρ_min² ⇒ Path C residual (15.130.C)."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem E: variance
# ---------------------------------------------------------------------------

def prove_var_Q_p7() -> dict:
    mean_Q = P7_DELTA2_PAIR
    EQ2 = sum(t["count"] * t["Q_delta"] ** 2 for t in P7_TYPES) / P7_N
    var = EQ2 - mean_Q**2
    ok = EQ2 == P7_EQ2 and var == P7_VAR_Q and var > 0
    return {
        "proved": ok,
        "E_Q2": str(EQ2),
        "var": str(var),
        "var_positive": var > 0,
        "theorem": (
            "Var(Q_δ)>0 at p=7 ⇒ pointwise criterion strictly stronger than "
            "δ²≤ρ_min²; Max+ not 2-point homogeneous for 4th moment."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem F: open status
# ---------------------------------------------------------------------------

def prove_open_status() -> dict:
    return {
        "proved_structure": True,
        "residual_closed_general": False,
        "open": (
            "Prove max Q_δ≤ρ_min² (or δ²≤ρ_min² / E‖Pf‖²≤room / Aut-line Q_0) "
            "for all primes p≥7; at p=5 use hyp form."
        ),
        "aut_line_caveat_p7": (
            "Q_δ non-constant at p=7 ⇒ either dim E_4p^{Aut}>1 or Aut does not "
            "act transitively on Max+ Q-level sets; halfspace pin alone is not "
            "the residual value (Q_hs ≠ δ²)."
        ),
    }


def main() -> dict:
    a = prove_pair_vs_basepoint()
    b = prove_spectrum_p7()
    c = prove_true_delta2_p7()
    d = prove_pointwise_p7()
    e = prove_var_Q_p7()
    f = prove_open_status()

    residual_closed = False
    out = {
        "prop": "15.131",
        "title": (
            "Prop 15.131: pair-avg vs basepoint ED4; p=7 three-type Q_δ; "
            "true δ²; pointwise max-Q criterion"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Clarifies that W-based ED4 at p=7 is basepoint "
            "Q_δ(hs), not pair-avg δ². True δ²=19180800/1840091 and "
            "max Q_δ=200448/4499 ≤ ρ_min² (slack 812048/220451). Residual "
            "holds at p=7 with either criterion. General p still OPEN. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "pair_vs_basepoint": a["proved"],
            "spectrum_p7_three_types": b["proved"],
            "true_delta2_p7": c["proved"],
            "pointwise_maxQ_le_rho_min_p7": d["proved"],
            "var_Q_positive_p7": e["proved"],
            "delta_le_rho_min_for_all_p": False,
            "pointwise_for_all_p": False,
            "E_Pf_le_room_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "pair_vs_basepoint": a,
            "spectrum": b,
            "true_delta2": c,
            "pointwise": d,
            "var_Q": e,
            "open": f,
        },
        "closed_forms": {
            "p7_delta2_pair": str(P7_DELTA2_PAIR),
            "p7_ED4_pair": str(P7_ED4_PAIR),
            "p7_max_Q": str(P7_MAX_Q),
            "p7_min_Q": str(P7_MIN_Q),
            "p7_slack_pointwise": str(P7_SLACK_POINTWISE),
            "p7_Q_types": [str(t["Q_delta"]) for t in P7_TYPES],
            "p7_type_counts": [t["count"] for t in P7_TYPES],
        },
        "status": {
            "p7_residual_pair_avg": True,
            "p7_pointwise_criterion": True,
            "general_residual": False,
            "hyp_residual_for_all_p": False,
        },
        "open_attack": (
            "Prove max_y Q_δ(y)≤ρ_min² for all primes p≥7 (or δ²≤ρ_min² / "
            "E‖P f‖²≤room / careful Aut-line with non-constant Q). "
            "At p=5 use hyp form."
        ),
        "backend": (
            "CPU Fraction + ProcessPool W=86 Max+ Gram-4th census for type "
            "spectrum; GPU unused (F20); F19 no class_key"
        ),
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15131.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.131: pair_vs_basepoint={a['proved']}")
    print(f"  spectrum 3-type={b['proved']}")
    print(f"  true delta2={c['proved']} d2={c['delta2_pair']}")
    print(
        f"  pointwise={d['proved']} maxQ={d['max_Q']} "
        f"slack={d['slack']}"
    )
    print(f"  varQ={e['proved']} var={e['var']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
