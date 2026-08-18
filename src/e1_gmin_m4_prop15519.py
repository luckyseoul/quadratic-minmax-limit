#!/usr/bin/env python3
"""
Prop 15.519 — QR0 halfspace Q is 2(p+1) on F_p^× for p≡3 (mod 4).

  1_{QR0}=(1+χ)/2 on F_p^×, plus 1 at 0.  Additive Fourier:
      1̂(α)=1/2+(1/2)χ(α)G(χ)   (α≠0),
  G(χ)²=χ(−1)p.  For p≡3, χ(−1)=−1 so G is purely imaginary and
      |1̂(α)|²=(p+1)/4.
  Hence the 15.292/15.310 Fejer formula gives
      Q_{QR0}(r)/q² = 2(p+1)  for r∈F_p^×,
      Q_{QR0}(r)=0            for r∉F_p  (15.292 A).
  Fail: claim |1̂|²=(p+1)/4 at p=13 (p≡1; live hat2 depends on χ).
  Fail: claim Q=2(p+1) on the p=5 AP/QR0 orbit (16/3≠12).
  Does not name ensemble Q_τ (NL / |H_+| still unnamed).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.518 flags.
"""
from __future__ import annotations

import cmath
import json
import math
import os
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, is_prime
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_sub_over_q2
from e1_gmin_m4_prop15309 import qr0
from e1_gmin_m4_prop15316 import chi_p_pm

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15519.json"

CERT = (7, 11, 19, 23, 31)


def hat2_int(A, alpha: int, p: int) -> int:
    s = 0j
    for x in A:
        s += cmath.exp(2j * math.pi * alpha * x / p)
    return int(round(abs(s) ** 2))


def hat2_QR0_named(p: int) -> Fraction:
    """(p+1)/4.  Valid for p≡3 (mod 4) only."""
    return Fraction(p + 1, 4)


def Q_QR0_Fp_named(p: int) -> Fraction:
    """2(p+1).  Valid for p≡3 (mod 4) only."""
    return Fraction(2 * (p + 1))


def Q_QR0_Fp_wrong_p5_ap(p: int) -> Fraction:
    """Fail-when-wrong: the p=5 AP value 16/3, or 2(p+1) at p=5."""
    return Fraction(2 * (p + 1))


def prove_A() -> dict:
    """G(χ)²=χ(−1)p ⇒ |1̂_QR0|²=(p+1)/4 for p≡3."""
    ok = True
    rows = {}
    for p in CERT:
        S = qr0(p)
        want = int(hat2_QR0_named(p))
        vals = {a: hat2_int(S, a, p) for a in range(1, p)}
        hit = all(v == want for v in vals.values())
        if not hit or p % 4 != 3 or chi_p_pm(p, p - 1) != -1:
            ok = False
        rows[str(p)] = {"want": want, "vals": vals, "ok": hit}
    # fail-when-wrong: p=13 ≡1
    p = 13
    S = qr0(p)
    vals13 = {a: hat2_int(S, a, p) for a in range(1, p)}
    named13 = float(hat2_QR0_named(p))
    if all(abs(v - named13) < 1e-6 for v in vals13.values()):
        ok = False
    rows["13"] = {
        "want_wrong": named13,
        "vals": vals13,
        "constant": len(set(vals13.values())) == 1,
        "equals_named": all(abs(v - named13) < 1e-6 for v in vals13.values()),
    }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "|1̂_QR0(α)|²=(p+1)/4 on F_p^× for p≡3. Fail p=13.",
    }


def prove_B() -> dict:
    """Fejer 32/(p+1)·((p+1)/4)²=2(p+1) on F_p^×."""
    ok = True
    rows = {}
    for p in CERT:
        S = qr0(p)
        want = Q_QR0_Fp_named(p)
        qs = {}
        for r in range(1, p):
            h1 = hat2_int(S, 1, p)
            hr = hat2_int(S, r, p)
            qs[str(r)] = str(Fraction(32, p + 1) * h1 * hr)
        hit = all(Fraction(qs[str(r)]) == want for r in range(1, p))
        if not hit:
            ok = False
        rows[str(p)] = {"want": str(want), "ok": hit}
    # fail p=5: 2(p+1)=12 ≠ AP Q_sub=16/3
    if Q_QR0_Fp_wrong_p5_ap(5) == Q_1d_sub_over_q2(5):
        ok = False
    rows["5"] = {
        "two_p1": str(Q_QR0_Fp_wrong_p5_ap(5)),
        "Q_AP_sub": str(Q_1d_sub_over_q2(5)),
        "neq": Q_QR0_Fp_wrong_p5_ap(5) != Q_1d_sub_over_q2(5),
    }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_QR0=2(p+1) on F_p^× for p≡3. Fail 12=16/3 at p=5.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "QR0 Type+ Q is named for p≡3. Ensemble mix still has NL "
            "and |H_+|. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.519 QR0 Q=2(p+1) on F_p^× for p≡3", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.519",
        "title": "QR0 Q=2(p+1) on F_p^× for p≡3 (mod 4)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
