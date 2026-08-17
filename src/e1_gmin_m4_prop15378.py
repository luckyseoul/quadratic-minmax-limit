#!/usr/bin/env python3
"""
Prop 15.378 — QR0⊙NC Q_τ interpolants die at p=19.  The unique
deg≤2 / k fit of Q_N* through p=7,11,19 is not a name.  Named
catalog has no survivor.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.377 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  QR0⊙NC live Q from 15.377 first_nc:
  p=7  (40/9, 24/7, 76/21)
  p=11 (3976/495, 128/55, 144/55)
  p=19 (15064/855, 256/171, 272/171)
k=p(p−1)/2.

============================================================================
Theorem A — PROVED (15.377 constructor; certified p=7,11,19).
  Q-- ≟ 2(p+53)/(5p) hits p=7,11 and dies at p=19
  (144/95 ≠ 256/171).  Q++ ≟ Q_AP dies at p=11 and p=19.
  Fail: claim the interpolant at p=19.  ∎

Theorem B — PROVED (unique interpolating quadratic).
  Q_N* · k equals 76, 144, 272 at p=7,11,19.  The unique
  polynomial of degree ≤2 through those three points is
      (−p²+222p−593)/12,
  hence Q_N* = (−p²+222p−593)/(6p(p−1)).  This is the
  interpolant, not a Gauss/Jacobi name: any 3 points determine
  a unique deg≤2.  Fail: claim the formula is independent of
  the p=19 sample (drop p=19 and the 2-point line
  (17p−43)/k already dies at p=19: 280/171 ≠ 272/171).  ∎

Theorem C — OPEN.  No catalog name hits all three primes.
  Ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Fraction arithmetic on 15.377 constructor values.
Writes evidence/e1_gmin_m4_prop15378.json
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
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named  # noqa: E402
from e1_gmin_m4_prop15367 import k_named  # noqa: E402
from e1_gmin_m4_prop15377 import Qmm_interp_p7p11, qr0_nc_Q  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15378_catalog.json"


def Qn_3pt_fit(p: int) -> Fraction:
    """Unique deg≤2 / k interpolant of QR0⊙NC Q_N* at p=7,11,19."""
    return Fraction(-p * p + 222 * p - 593, 6 * p * (p - 1))


def Qn_2pt_fit(p: int) -> Fraction:
    """2-point (p=7,11) line (17p−43)/k. Dies at p=19."""
    return Fraction(17 * p - 43, k_named(p))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (7, 11, 19):
        Q = qr0_nc_Q(p)
        interp = Qmm_interp_p7p11(p)
        qap = Q_AP_named(p)
        rows[str(p)] = {
            "Qmm": str(Q["--N1"]),
            "interp": str(interp),
            "eq_interp": Q["--N1"] == interp,
            "Qpp": str(Q["++"]),
            "Q_AP": str(qap),
            "eq_QAP": Q["++"] == qap,
        }
        if p in (7,):
            ok = ok and Q["--N1"] == interp
            ok = ok and Q["++"] == qap
        else:
            ok = ok and Q["++"] != qap
        if p == 19:
            ok = ok and Q["--N1"] != interp
            ok = ok and Q["--N1"] == Fraction(256, 171)
            ok = ok and interp == Fraction(144, 95)
        if p == 19 and Q["--N1"] == interp:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "2(p+53)/(5p) and Q++=Q_AP die at p=19. "
            "Fail: claim the interpolant at p=19."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p, want in ((7, 76), (11, 144), (19, 272)):
        Q = qr0_nc_Q(p)
        kQ = Q["N*"] * k_named(p)
        ok = ok and kQ == want
        ok = ok and Qn_3pt_fit(p) == Q["N*"]
        rows[str(p)] = {
            "kQn": str(kQ),
            "fit": str(Qn_3pt_fit(p)),
            "two_pt": str(Qn_2pt_fit(p)),
            "eq_3pt": Qn_3pt_fit(p) == Q["N*"],
            "eq_2pt": Qn_2pt_fit(p) == Q["N*"],
        }
    # 2-point dies at p=19
    ok = ok and Qn_2pt_fit(19) != Fraction(272, 171)
    ok = ok and Qn_2pt_fit(19) == Fraction(280, 171)
    if Qn_2pt_fit(19) == Fraction(272, 171):
        ok = False
    # 3-point is the interpolant: it cannot fail the three samples
    ok = ok and all(Qn_3pt_fit(p) == qr0_nc_Q(p)["N*"] for p in (7, 11, 19))
    return {
        "proved": bool(ok),
        "rows": rows,
        "two_pt_p19": str(Qn_2pt_fit(19)),
        "live_qn_p19": "272/171",
        "theorem": (
            "3-point deg≤2 / k fit of Q_N* is the interpolant. "
            "2-point (17p−43)/k dies at p=19 (280/171≠272/171)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "note": (
            "No catalog name survives p=19. 3-point fit is not a name. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.378  QR0⊙NC interpolants die at p=19", flush=True)
    A = prove_A()
    print(f"  A interp: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B 3pt: {B['proved']} 2pt19={B['two_pt_p19']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.378",
        "title": "QR0⊙NC interpolants die at p=19; 3-point QN fit is not a name",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Qmm_interp_dies_p19": A["proved"],
            "Qn_3pt_is_interpolant": B["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15378.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
