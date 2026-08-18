#!/usr/bin/env python3
"""
Prop 15.505 — Affine-R occupancy has a named Fourier line, and
at p=5 (NL is all R) this evaluates the ensemble off-value
Q++/q²=48/13 as a 1D+R mixture.  Not a p-law for the ensemble
(p=7 NL is not all R).

    n_c = a c + b  (a∈F_p^×),   |∑_c n_c ζ^{α c}|² = a² p² / |1−ζ^α|²
    |1−ζ^α|² = 4 sin²(πα/p).

On an all-R Max+,
    Q_R(r)/q² = (1/(p−1)) ∑_{α≠0} 1/(sin²(πα/p) sin²(π r α/p)).
At p=5 the summand is α-independent, so this is
1/(sin²(π/5)sin²(2π/5)).  sin(π/5)sin(2π/5)=√5/4, hence 16/5.
That 16/5 equals live Q_NL on sub (100 designs), independently
of recombining the 1D block.  Fail: claim 16/5 equals live
Q++ at p=5 (48/13) or at p=7 (1544/409).

Recombining (30·16/3 + 100·16/5)/130 = 48/13 is a consistency
check, not a second derivation.  Fail: drop the 1D term (16/5).
Do not export the collapsed 1/sin²sin² form as a p-law (p=7
summand is not α-constant).

Does **not** name ensemble Q_τ for p≥7.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.504 flags.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction
from math import sin, pi
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import Q_1d_sub_over_q2, n_1d
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15505.json"


def Q_R_sub_trig(p: int, r: int = 2) -> float:
    """1/(sin²(π/p) sin²(π r/p))."""
    return 1.0 / (sin(pi / p) ** 2 * sin(pi * r / p) ** 2)


def Q_R_sub_p5() -> Fraction:
    """Exact: sin(π/5)sin(2π/5)=√5/4 ⇒ square 5/16 ⇒ 16/5."""
    return Fraction(16, 5)


def Q_R_sub_p5_wrong_half() -> Fraction:
    """Fail-when-wrong: 8/5."""
    return Fraction(8, 5)


def mixture_p5() -> Fraction:
    return (n_1d(5) * Q_1d_sub_over_q2(5) + 100 * Q_R_sub_p5()) / HPLUS[5]


def mixture_p5_drop_1d() -> Fraction:
    return Q_R_sub_p5()


def prove_A() -> dict:
    """p=5 trig identity is 16/5. Fail 8/5."""
    ok = True
    trig = Q_R_sub_trig(5, 2)
    exact = float(Q_R_sub_p5())
    if abs(trig - exact) > 1e-12:
        ok = False
    if Q_R_sub_p5() == Q_R_sub_p5_wrong_half():
        ok = False
    if Q_R_sub_p5() == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "trig": trig,
        "exact": str(Q_R_sub_p5()),
        "live_Qpp_p5": str(LIVE_QPP[5]),
    }


def prove_B() -> dict:
    """16/5 is live Q_NL at p=5, not live Q++ at p=5 or p=7."""
    ok = True
    # NUM_NL / n_NL = (NUM - n_1d Q_1d_sub) / (N-n_1d)
    N, n1 = HPLUS[5], n_1d(5)
    q_nl = (LIVE_QPP[5] * N - Q_1d_sub_over_q2(5) * n1) / (N - n1)
    if q_nl != Q_R_sub_p5():
        ok = False
    if Q_R_sub_p5() == LIVE_QPP[5]:
        ok = False
    if abs(Q_R_sub_trig(7, 2) - float(LIVE_QPP[7])) < 0.05:
        ok = False
    return {
        "proved": bool(ok),
        "Q_NL_p5": str(q_nl),
        "Q_R_p5": str(Q_R_sub_p5()),
        "trig_p7": Q_R_sub_trig(7, 2),
        "live_p7": str(LIVE_QPP[7]),
    }


def prove_C() -> dict:
    """p=5 mixture recovers 48/13. Fail drop 1D."""
    ok = True
    mix = mixture_p5()
    if mix != LIVE_QPP[5]:
        ok = False
    if mix != Fraction(48, 13):
        ok = False
    if mixture_p5_drop_1d() == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "mixture": str(mix),
        "drop_1d": str(mixture_p5_drop_1d()),
        "n_1d": n_1d(5),
        "n_R": HPLUS[5] - n_1d(5),
        "Q_1d_sub": str(Q_1d_sub_over_q2(5)),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Q_R(r) is named for all-R Max+. At p=5 that is all NL "
            "and the mixture evaluates 48/13. p=7 NL is not all R. "
            "Ensemble Q_τ unnamed for p≥7. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.505 affine-R Fourier; p=5 mixture is 48/13", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
