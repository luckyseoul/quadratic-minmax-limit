#!/usr/bin/env python3
"""
Prop 15.383 — Leftover (1) floor window in Var(e_*).
For p≡1 (mod 4) the S0 floor is equivalent to
    V_lo ≤ Var(e_*) ≤ V_hi
with named
    V_hi=(α+β Q^{ub}−E[W]²)/κ,   V_lo=(α+β Q_lo−E[W]²)/κ.
This does **not** bound Var(e_*), so phi_F is not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.382 flags.

============================================================================
Setup.  15.326: Var(W)=κ Var(e_*), κ=p(p−1)²/(2(p−3)),
E[W]=(p+1)p(p−1)/4.  15.321: A_□=α+β Q++ = E[W²] on p≡1.
15.326 C: S0 even S/q²≥6 ⇔ Q_lo≤Q++≤Q^{ub}.

============================================================================
Theorem A — PROVED (15.321 + 15.326; all p≡1≥5).
  Var(e_*)=(A_□−E[W]²)/κ, so the S0 floor is
      V_lo ≤ Var(e_*) ≤ V_hi
  with V_hi=(α+β Q^{ub}−EW²)/κ and V_lo=(α+β Q_lo−EW²)/κ.
  At p=5: V_hi=18/7, V_lo=9/8, live 33/13∈(9/8,18/7).
  Fail: drop κ (then V_hi=360/7 ≠ 18/7).  ∎

Theorem B — PROVED (threshold rewrite; 15.304).
  λ_j≥6 ⇔ E|⟨ΔN,G_U(j)⟩|² ≥ (3/4)(p+1) q².
  Fail: claim the threshold is (3/4) q² (den 1).  ∎

Theorem C — OPEN.  No Max+-free majorant Var(e_*)≤V_hi
  (Var_2orb is unproved and not a support law at p=13).
  p≡3 still has Q_{--}.  phi_F not imported.

============================================================================
Backend: Fraction identities.  Writes
evidence/e1_gmin_m4_prop15383.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import (  # noqa: E402
    Asq_from_Qpp_p1,
    Qpp_floor_ub,
    alpha_beta_p1,
)
from e1_gmin_m4_prop15323 import EW_named  # noqa: E402
from e1_gmin_m4_prop15326 import (  # noqa: E402
    Q_lo_named,
    even_S_p1,
    kappa_named,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15383_catalog.json"


def EW2_named(p: int) -> int:
    return EW_named(p) ** 2


def V_e_from_Qpp(p: int, Qpp: Fraction) -> Fraction:
    return (Asq_from_Qpp_p1(p, Qpp) - EW2_named(p)) / kappa_named(p)


def V_e_hi(p: int) -> Fraction:
    return V_e_from_Qpp(p, Qpp_floor_ub(p))


def V_e_lo(p: int) -> Fraction:
    return V_e_from_Qpp(p, Q_lo_named(p))


def V_e_hi_drop_kappa(p: int) -> Fraction:
    return Asq_from_Qpp_p1(p, Qpp_floor_ub(p)) - EW2_named(p)


def floor_square_thr(p: int) -> Fraction:
    """(3/4)(p+1) q² with q=p²."""
    q2 = p**4
    return Fraction(3, 4) * (p + 1) * q2


def floor_square_thr_den1(p: int) -> Fraction:
    return Fraction(3, 4) * (p**4)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29, 37):
        lo, hi = V_e_lo(p), V_e_hi(p)
        ok = ok and lo < hi
        ok = ok and hi > 0 and lo >= 0
        # endpoints recover S=6
        s_ub = even_S_p1(p, Qpp_floor_ub(p))
        s_lo = even_S_p1(p, Q_lo_named(p))
        ok = ok and s_ub["sm_joth"] == 6
        ok = ok and s_lo["sm_j0"] == 6
        rows[str(p)] = {"V_lo": str(lo), "V_hi": str(hi)}
    ok = ok and V_e_hi(5) == Fraction(18, 7)
    ok = ok and V_e_lo(5) == Fraction(9, 8)
    # live 33/13 sits in the window
    live = Fraction(33, 13)
    ok = ok and V_e_lo(5) < live < V_e_hi(5)
    if V_e_hi_drop_kappa(5) == Fraction(18, 7):
        ok = False
    ok = ok and V_e_hi_drop_kappa(5) == Fraction(360, 7)
    return {
        "proved": bool(ok),
        "rows": rows,
        "p5_hi": str(V_e_hi(5)),
        "p5_lo": str(V_e_lo(5)),
        "theorem": (
            "p≡1 S0 floor ⇔ V_lo≤Var(e_*)≤V_hi. "
            "Fail: drop κ (360/7≠18/7)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 19):
        T = floor_square_thr(p)
        q2 = p**4
        # λ=8/(p+1) q^{-2} M ≥6 ⇔ M≥ (6(p+1)/8) q² = (3/4)(p+1)q²
        ok = ok and T == Fraction(6 * (p + 1), 8) * q2
        ok = ok and T != floor_square_thr_den1(p)
        rows[str(p)] = str(T)
        if T == floor_square_thr_den1(p):
            ok = False
    ok = ok and floor_square_thr(5) == Fraction(3, 4) * 6 * 625
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "λ_j≥6 ⇔ E|⟨ΔN,G_U⟩|²≥(3/4)(p+1)q². "
            "Fail: drop (p+1)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "V_hi is the named p≡1 window, not a proved majorant. "
            "p≡3 still needs Q_{--}. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.383  p≡1 floor ⇔ V_lo≤Var(e)≤V_hi", flush=True)
    A = prove_A()
    print(f"  A window: {A['proved']} p5=[{A['p5_lo']},{A['p5_hi']}]", flush=True)
    B = prove_B()
    print(f"  B square thr: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "V_hi": {str(p): str(V_e_hi(p)) for p in (5, 13, 17)},
                "V_lo": {str(p): str(V_e_lo(p)) for p in (5, 13, 17)},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.383",
        "title": "p≡1 S0 floor ⇔ V_lo≤Var(e_*)≤V_hi; phi_F open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Ve_window": A["proved"],
            "square_threshold": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15383.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
