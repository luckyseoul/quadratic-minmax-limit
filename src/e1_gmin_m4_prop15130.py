#!/usr/bin/env python3
"""
Prop 15.130 — Identity P m₄=δ; Jensen refined; ρ_min-sufficient residual for
p≥7 with closed gap; census δ²≤ρ_min²; Aut-line program; residual OPEN generally.

Continues 15.129 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. ρ = m₄ − κ/p² = ρ_min + δ, δ ∈ E_{4p}=ker(4pI−T), ρ_min ⊥ E_{4p},
κ ⊥ E_{4p} (Prop 15.102/15.115).  f_y(S)=∏_{i∈S} y_i, m₄=E_y f_y.
P := P_{E_{4p}}.  Path C residual ⇔ δ² ≤ room_hyp/24.

============================================================================
Theorem A (P m₄ = δ) — PROVED.
  m₄ = κ/p² + ρ_min + δ with κ/p² ⊥ E_{4p} and ρ_min ⊥ E_{4p}, hence
      P m₄ = δ.
  In particular ‖P m₄‖₂² = δ².  Combined with m₄ = E_y f_y,
      δ = P E_y f_y = E_y[P f_y]
  (P linear continuous).  ∎

Theorem B (Jensen, refined) — PROVED.
  δ² = ‖E[P f_y]‖₂² ≤ E[‖P f_y‖₂²], equality iff y ↦ P f_y is a.s. constant
  on Max+.  Therefore
      E[‖P f_y‖₂²] ≤ room_hyp/24  ⇒  Path C residual.
  At p=3: E_4p={0}, both sides 0.  At p=5: ED4 saturates, single residual type,
  equality expected.  At p=7: multiple ED4 types ⇒ strict inequality expected.  ∎

Theorem C (ρ_min-sufficient residual for all primes p≥7) — PROVED Fraction.
  From Prop 15.110/15.102:
      ρ_min² = 5 n (p²−1)(p²+3) / (6 p² (p²−5)),   n=p²+1,
      room_hyp/24 = 4(p²−9)(p²−1)² / (3(p²−5)(p²+1)).
  Their difference is the Max+-free rational
      room_hyp/24 − ρ_min²
      = (p−1)(p+1)(3p⁶ − 105p⁴ + 37p² − 15)
        / (6 p² (p²−5)(p²+1)).
  For primes p≥7: numerator and denominator are positive
  (3p⁶−105p⁴+37p²−15 > 0 for p²≥49, as the cubic in x=p² has unique real
  root ≈34.65).  Hence ρ_min² < room_hyp/24 for all primes p≥7, and
      δ² ≤ ρ_min²  ⇒  δ² < room_hyp/24  ⇒  Path C residual.
  (At p=5: ρ_min²=728/25 > room=1536/65, so the ρ_min criterion is strictly
  stronger than residual and not necessary; residual uses the tighter hyp form.)  ∎

Theorem D (census: δ² ≤ ρ_min² at p=3,5,7) — CERTIFIED.
  From full W census (15.128):
    p=3: δ²=0 ≤ ρ_min²;
    p=5: δ²=1536/65 < 728/25=ρ_min² (and = room);
    p=7: δ²=82176/4499 < 26000/539=ρ_min² (ratio ≈0.379).
  Thus the sufficient criterion of Theorem C holds at the only prime p=7 among
  {5,7} where it applies, and residual holds at all three census primes.  ∎

Theorem E (Aut-line reduction, restated) — PROVED structure (15.109.F).
  If dim E_{4p}^{Aut} ≤ 1 then δ=0 or δ=c v_0 (‖v_0‖=1 Aut-invariant), and
  residual ⇔ c² ≤ room_hyp/24.  Moreover c = ⟨m₄, v_0⟩ = E_y[⟨f_y, v_0⟩]
  = Q_0(y), and if Q_0 is constant on Max+ then c = Q_0(x_hs) for the
  halfspace boolean evec (Max+-free evaluation of v_0 against hs).  ∎

Theorem F (OPEN residual, honest).
  Prove for all primes p≥5 either
    (i)  δ² ≤ ρ_min²  (sufficient for p≥7 by Thm C; at p=5 use hyp form), or
    (ii) E[‖P f_y‖₂²] ≤ room_hyp/24, or
    (iii) dim E_{4p}^{Aut}≤1 and c²=Q_0(hs)² ≤ room_hyp/24 via Gauss sums, or
    (iv) closed W_k ⇒ ED4≤ED4_bud.
  Census closes (i)–(ii) at p=3,5,7 only.  ∎

============================================================================
Certified: Thm A–C algebra (Fraction); Thm D census p=3,5,7; gap formula.
Backend: CPU Fraction; GPU unused (F20).  F19: no class_key thrash
(Aut-line only as structural reduction, not orbit refinement).

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15130.json
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
from e1_gmin_m4_prop15119 import ed4_flat_closed, m4f_flat  # noqa: E402
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
# Theorem A: P m4 = delta
# ---------------------------------------------------------------------------

def prove_P_m4_is_delta() -> dict:
    return {
        "proved": True,
        "theorem": (
            "m₄=κ/p²+ρ_min+δ with κ/p²,ρ_min ⊥ E_{4p} ⇒ P m₄=δ.  "
            "Since m₄=E f_y, δ=E[P f_y]."
        ),
        "references": "Props 15.102, 15.115, 15.129",
    }


# ---------------------------------------------------------------------------
# Theorem B: Jensen
# ---------------------------------------------------------------------------

def prove_jensen_refined() -> dict:
    return {
        "proved": True,
        "theorem": (
            "δ²=‖E[P f_y]‖² ≤ E‖P f_y‖² by convexity; "
            "E‖P f_y‖² ≤ room_hyp/24 ⇒ residual."
        ),
        "equality": "P f_y a.s. constant on Max+",
    }


# ---------------------------------------------------------------------------
# Theorem C: gap algebra rho_min < room for p>=7
# ---------------------------------------------------------------------------

def gap_room_minus_rho_min(p: int) -> Fraction:
    """room_hyp/24 − ρ_min² closed form."""
    # (p-1)(p+1)(3p^6-105p^4+37p^2-15) / (6 p^2 (p^2-5)(p^2+1))
    num = (
        (p - 1)
        * (p + 1)
        * (3 * p**6 - 105 * p**4 + 37 * p * p - 15)
    )
    den = 6 * p * p * (p * p - 5) * (p * p + 1)
    return Fraction(num, den)


def prove_gap_algebra(primes: list[int]) -> dict:
    rows = {}
    ok = True
    # cubic positivity for x>=49
    # f(x)=3x^3-105x^2+37x-15; f'(x)=9x^2-210x+37; check f(49)>0
    f49 = 3 * 49**3 - 105 * 49**2 + 37 * 49 - 15
    cubic_pos = f49 > 0
    for p in primes:
        if p * p == 5:
            continue
        gap = room_hyp_over_24(p) - rho_min2(p)
        gap_c = gap_room_minus_rho_min(p)
        match = gap == gap_c
        pos = gap > 0 if p >= 7 else gap < 0  # p=5: room < rho_min
        rows[str(p)] = {
            "gap": str(gap),
            "gap_closed": str(gap_c),
            "match": match,
            "rho_min2": str(rho_min2(p)),
            "room": str(room_hyp_over_24(p)),
            "rho_min_lt_room": rho_min2(p) < room_hyp_over_24(p),
            "sign_ok": (gap > 0) if p >= 7 else (gap < 0 or p == 3),
        }
        if not match:
            ok = False
        if p >= 7 and not (gap > 0):
            ok = False
    # all primes 7..97
    all_ge7 = True
    for p in range(7, 100):
        if is_prime(p) and gap_room_minus_rho_min(p) <= 0:
            all_ge7 = False
            break
    return {
        "proved": ok and cubic_pos and all_ge7,
        "cubic_f49_positive": cubic_pos,
        "all_primes_7_to_97_gap_positive": all_ge7,
        "theorem": (
            "room_hyp/24 − ρ_min² = (p−1)(p+1)(3p⁶−105p⁴+37p²−15)"
            "/(6p²(p²−5)(p²+1)) > 0 for all primes p≥7; "
            "hence δ²≤ρ_min² ⇒ residual for p≥7."
        ),
        "closed_form": (
            "(p−1)(p+1)(3p⁶−105p⁴+37p²−15)/(6p²(p²−5)(p²+1))"
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D: census delta2 <= rho_min2
# ---------------------------------------------------------------------------

def delta2_from_W(p: int) -> Fraction:
    ED4 = ED4_from_W(W_CENSUS[p], n_of(p))
    return (ED4 - ed4_flat_closed(p)) / 24


def prove_census_delta_le_rho_min() -> dict:
    rows = {}
    ok = True
    for p in (3, 5, 7):
        d2 = delta2_from_W(p)
        rm = rho_min2(p)
        room = room_hyp_over_24(p)
        rows[str(p)] = {
            "delta2": str(d2),
            "rho_min2": str(rm),
            "room": str(room),
            "delta2_le_rho_min2": d2 <= rm,
            "delta2_le_room": d2 <= room,
            "ratio_delta_over_rho": str(d2 / rm) if rm else None,
            "ratio_delta_over_room": str(d2 / room) if room else None,
        }
        if not (d2 <= rm and d2 <= room):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Census δ²≤ρ_min² and δ²≤room at p=3,5,7 "
            "(p=7 ratio δ²/ρ_min²=82176/4499 ÷ 26000/539 = 251664/664625)."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E: Aut-line structure
# ---------------------------------------------------------------------------

def prove_aut_line_program() -> dict:
    return {
        "proved_structure": True,
        "theorem": (
            "If dim E_{4p}^{Aut}≤1 then residual ⇔ c²≤room_hyp/24 with "
            "c=⟨m₄,v_0⟩=Q_0(hs) when Q_0 constant (15.109.F / 15.110.F)."
        ),
        "open": (
            "Prove dim E_{4p}^{Aut}≤1 for all p≥5 and evaluate Q_0(hs) "
            "by Gauss sums / halfspace character sums."
        ),
        "caveat_F19": (
            "Do not thrash class_key refinements; use Aut-invariant algebra "
            "of bounded type count or character sums on PGL."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 80) if is_prime(p)]
    a = prove_P_m4_is_delta()
    b = prove_jensen_refined()
    c = prove_gap_algebra(primes)
    d = prove_census_delta_le_rho_min()
    e = prove_aut_line_program()

    residual_closed = False
    out = {
        "prop": "15.130",
        "title": (
            "Prop 15.130: P m₄=δ; Jensen; ρ_min-sufficient residual for p≥7; "
            "census δ²≤ρ_min²; Aut-line program"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. For p≥7, residual follows from δ²≤ρ_min² "
            "(gap algebra proved). Census confirms at p=7. General δ²≤ρ_min² "
            "or E‖P f‖²≤room still open. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "P_m4_equals_delta": a["proved"],
            "jensen_refined": b["proved"],
            "gap_rho_min_lt_room_p_ge_7": c["proved"],
            "census_delta_le_rho_min": d["proved"],
            "aut_line_structure": e["proved_structure"],
            "delta_le_rho_min_for_all_p": False,
            "E_Pf_le_room_for_all_p": False,
            "weight_enumerator_closed_general": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "P_m4": a,
            "jensen": b,
            "gap": c,
            "census": d,
            "aut_line": e,
        },
        "status": {
            "rho_min_criterion_p_ge_7": True,
            "general_residual": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "gap_room_minus_rho_min": c["closed_form"],
            "rho_min2": "5n(p²−1)(p²+3)/(6p²(p²−5))",
            "room_hyp_over_24": "4(p²−9)(p²−1)²/(3(p²−5)(p²+1))",
            "P_m4": "δ",
            "jensen": "δ² ≤ E‖P f_y‖²",
        },
        "open_attack": (
            "Prove δ²≤ρ_min² for all primes p≥7 (sufficient by Thm C), "
            "or E‖P_{E_4p} f_y‖²≤room_hyp/24, or dim E_{4p}^{Aut}≤1 + "
            "Gauss-sum evaluation of Q_0(hs)²≤room_hyp/24, or closed W_k."
        ),
        "backend": (
            "CPU Fraction algebra + census from 15.128; GPU unused (F20); "
            "F19 no class_key thrash"
        ),
        "F19": "no class_key thrash; Aut-line as structure only",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15130.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.130: P m4=delta proved={a['proved']}")
    print(f"  Jensen proved={b['proved']}")
    print(f"  gap algebra proved={c['proved']} all_7_97={c['all_primes_7_to_97_gap_positive']}")
    print(f"  census delta<=rho_min proved={d['proved']}")
    for p in (5, 7, 11):
        print(
            f"  p={p}: gap={c['by_p'][str(p)]['gap']} "
            f"rho_lt_room={c['by_p'][str(p)]['rho_min_lt_room']}"
        )
    for p in (3, 5, 7):
        print(
            f"  census p={p}: d2<=rm={d['by_p'][str(p)]['delta2_le_rho_min2']} "
            f"d2<=room={d['by_p'][str(p)]['delta2_le_room']}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
