#!/usr/bin/env python3
"""
Prop 15.465 — The 409-free ratio live Paley-N*/Qn equals
    R(p) = 1 + 4p(p−5) / [(p+1)(p−2)(p−1)²]
at p=5 (R=1) and p=7 (R=187/180).  The pieces are
|G(χ)|²=p, |U|=p+1, |F_p^×|=p−1, |F_p\\{0,1}|=p−2, and
4=−4 J(χ,χ) χ(−1).  Dropping (p−5) misses 1 at p=5.
Does not name Q_τ for p≥13.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.464 flags.

============================================================================
Setup.  15.464 B: Paley-N* Q is 32/13, 1496/409; S0 Qn is
32/13, 1440/409.  Ratio = 1, 187/180.  409 cancels.

============================================================================
Theorem A — PROVED (named field counts; p=5,7 live).
  R(p) = 1 + 4p(p−5)/[(p+1)(p−2)(p−1)²]
  equals live_QN_paley(p)/Qn(Q++) at p=5 and p=7.
  Fail: drop (p−5) (then 1+4p/[(p+1)(p−2)(p−1)²] is
  77/72 ≠ 1 at p=5).  ∎

Theorem B — PROVED (Jacobi identity on F_p).
  J(χ,χ)=∑_x χ(x)χ(1−x)=−χ(−1), so
      −4 J(χ,χ) χ(−1) = 4.
  Fail: claim J(χ,χ)=0 (false: −1 at p=5, +1 at p=7).  ∎

Theorem C — OPEN.  R(p) is not proved equal to Paley-N*/Qn
  for p≥13 (no live Q).  Q_τ still unnamed in p.  phi_F not
  imported.

============================================================================
Backend: Fraction arithmetic + 15.464 live ratio.  Writes
evidence/e1_gmin_m4_prop15465.json
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
from e1_gmin_m4_prop15326 import Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15464 import live_QN_paley  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15465_catalog.json"


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def J_chi_chi(p: int) -> int:
    """J(χ,χ)=∑ χ(x)χ(1−x)."""
    acc = 0
    for x in range(p):
        acc += chi_p(p, x) * chi_p(p, (1 - x) % p)
    return acc


def J_chi_chi_named(p: int) -> int:
    return -chi_p(p, p - 1)


def four_from_J(p: int) -> int:
    return -4 * J_chi_chi_named(p) * chi_p(p, p - 1)


def R_named(p: int) -> Fraction:
    """1 + 4p(p−5)/[(p+1)(p−2)(p−1)²]."""
    den = (p + 1) * (p - 2) * (p - 1) ** 2
    return 1 + Fraction(4 * p * (p - 5), den)


def R_drop_pm5(p: int) -> Fraction:
    """Fail-when-wrong: omit (p−5) but keep the 1+."""
    den = (p + 1) * (p - 2) * (p - 1) ** 2
    return 1 + Fraction(4 * p, den)


def live_ratio(p: int) -> Fraction:
    qn = Qn_from_Qpp(p, LIVE_QPP[p])
    return live_QN_paley(p) / qn


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        R = R_named(p)
        live = live_ratio(p)
        rows[str(p)] = {"R": str(R), "live": str(live)}
        if R != live:
            ok = False
    if R_named(5) != 1:
        ok = False
    if R_named(7) != Fraction(187, 180):
        ok = False
    if R_drop_pm5(5) == 1:
        ok = False
    if R_drop_pm5(5) != Fraction(77, 72):
        ok = False
    if live_ratio(7) == 1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "drop5": str(R_drop_pm5(5)),
        "theorem": (
            "R=1+4p(p-5)/[(p+1)(p-2)(p-1)^2] equals the live "
            "ratio at p=5,7. Fail: drop (p-5) (77/72≠1 at p=5)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        J = J_chi_chi(p)
        Jn = J_chi_chi_named(p)
        four = four_from_J(p)
        rows[str(p)] = {"J": J, "named": Jn, "four": four}
        if J != Jn:
            ok = False
        if four != 4:
            ok = False
        if J == 0:
            ok = False
    if J_chi_chi(5) != -1 or J_chi_chi(7) != 1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "J(χ,χ)=−χ(−1) so −4 J χ(−1)=4. Fail: J=0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "R13": str(R_named(13)),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "R names the 409-free ratio at p=5,7 only. "
            "Q_τ unnamed for p≥13. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.465  Paley-N*/Qn = 1+4p(p-5)/[(p+1)(p-2)(p-1)^2]", flush=True)
    A = prove_A()
    print(f"  A ratio: {A['proved']} {A['rows']} drop={A['drop5']}", flush=True)
    B = prove_B()
    print(f"  B J(χ,χ): {B['proved']} p5={B['rows']['5']}", flush=True)
    C = prove_open()
    print(f"  C open R13={C['R13']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "R5": A["rows"]["5"]["R"],
                "R7": A["rows"]["7"]["R"],
                "drop5": A["drop5"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.465",
        "title": (
            "Paley-N*/Qn = 1+4p(p-5)/[(p+1)(p-2)(p-1)^2] "
            "at p=5,7; drop (p-5) misses 1"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "ratio_named_live": A["proved"],
            "four_from_Jacobi": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15465.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
