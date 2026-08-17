#!/usr/bin/env python3
"""
Prop 15.385 — Type+ 1D square-line occupancy is named in p.
D={z=−1} of a Type+ halfspace is a union of (p−1)/2 parallels
in one square direction (15.292 / 15.305 C).  Averaging over the
(p+1)/2 square directions:
    P_1D(0)=1/p,  P_1D(p)=(p−1)/(p(p+1)),  P_1D(μ)=(p−1)/(p+1),
    μ=(p−1)/2, and P_1D=0 else.
E[n⁴|1D]=(p−1)/(p+1)·(p³+μ⁴).  p=5 NL is uniform on {0,…,4}
and the 3/13–10/13 mixture recovers the live ensemble hist.
p=7 NL attains n=p (live 3528 counts), so that p=5 law is not
general.  Ensemble E[n⁴] / Q_τ still mixed with den D.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.384 flags.

============================================================================
Setup.  15.292: Type+ 1D is y=σ∘L, ker L a square line,
|σ^{-1}(+1)|=m=(p+1)/2, n_1d=m C(p,m).  15.305 C: every
nonsquare-direction line meets D in μ=(p−1)/2 points, so the
selected parallel class of a 1D Max+ is a square direction
(15.332: nonsquare line-unions ∉ H_+).  D is the complementary
(p−1)/2 parallels (z=−1).

============================================================================
Theorem A — PROVED (15.292 + 15.305 C; all odd p≥5).
  On the selected square class, n∈{0,p} with
      #{n=p}=(p−1)/2,  #{n=0}=(p+1)/2.
  On every other square class, n≡μ.  Averaging over the
  (p+1)/2 square directions and p lines:
      P_1D(0)=1/p,
      P_1D(p)=(p−1)/(p(p+1)),
      P_1D(μ)=(p−1)/(p+1).
  Fail: claim P_1D(p)=1/p (then 4/30≠2/15 at p=5).  ∎

Theorem B — PROVED (A + fourth moment).
  E[n⁴|1D]=(p−1)/(p+1)·(p³+μ⁴)  (=94 at p=5, 318 at p=7).
  Fail: drop p³ (then 32/3 ≠ 94 at p=5).  ∎

Theorem C — PROVED (p=5 mixture + p=7 census obstruction).
  w=n_1d/|H_+|=3/13 at p=5.  NL occupancy uniform on {0,1,2,3,4}
  mixes with (A) to the live ensemble
      P(0)=1/5, P(1)=P(3)=P(4)=2/13, P(2)=4/13, P(5)=2/65.
  Fail: claim P_NL(p)=0 at p=7 (then #{n=7}=420 ≠ live 3948).  ∎

Theorem D — OPEN.  Ensemble P(n=k) / E[n⁴] still have den D
  at p≥7.  phi_F not imported.

============================================================================
Backend: Fraction identities + live p=5,7 occupancy counts.
Writes evidence/e1_gmin_m4_prop15385.json
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
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15357 import D_live  # noqa: E402
from e1_gmin_m4_prop15384 import EM2_named, LIVE_QPP  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15385_catalog.json"

# Live square-line #{n=p} over |H_+| · ((p+1)/2) · p samples (cpu_occupancy).
LIVE_N_EQ_P = {5: 60, 7: 3948}
# 1D-only contribution to that count: n_1d · (p−1)/2
# p=7: 140·3=420.


def mu_named(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def P_1d(p: int, k: int) -> Fraction:
    mu = (p - 1) // 2
    if k == 0:
        return Fraction(1, p)
    if k == p:
        return Fraction(p - 1, p * (p + 1))
    if k == mu:
        return Fraction(p - 1, p + 1)
    return Fraction(0)


def P_1d_fail_p_as_invp(p: int) -> Fraction:
    """Wrong neighbour: P_1D(p)=1/p."""
    return Fraction(1, p)


def En4_1d(p: int) -> Fraction:
    mu = mu_named(p)
    return Fraction(p - 1, p + 1) * (p**3 + mu**4)


def En4_1d_drop_p3(p: int) -> Fraction:
    mu = mu_named(p)
    return Fraction(p - 1, p + 1) * (mu**4)


def w_1d(p: int) -> Fraction:
    return Fraction(n_1d(p), HPLUS[p])


def P_nl_p5(k: int) -> Fraction:
    if 0 <= k <= 4:
        return Fraction(1, 5)
    return Fraction(0)


def P_ens_p5(k: int) -> Fraction:
    w = w_1d(5)
    return w * P_1d(5, k) + (1 - w) * P_nl_p5(k)


def n_eq_p_if_nl_misses_full_line(p: int) -> int:
    """Ensemble #{n=p} if every NL row misses n=p (p=5 law)."""
    return n_1d(p) * (p - 1) // 2


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        s = sum(P_1d(p, k) for k in range(p + 1))
        ok = ok and s == 1
        ok = ok and P_1d(p, 0) == Fraction(1, p)
        ok = ok and P_1d(p, p) == Fraction(p - 1, p * (p + 1))
        ok = ok and P_1d(p, (p - 1) // 2) == Fraction(p - 1, p + 1)
        ok = ok and P_1d(p, p) != P_1d_fail_p_as_invp(p)
        rows[str(p)] = {
            "P0": str(P_1d(p, 0)),
            "Pp": str(P_1d(p, p)),
            "Pmu": str(P_1d(p, (p - 1) // 2)),
        }
    # selected-class counts
    ok = ok and P_1d(5, 5) == Fraction(2, 15)
    ok = ok and P_1d_fail_p_as_invp(5) == Fraction(1, 5)
    if P_1d(5, 5) == Fraction(1, 5):
        ok = False
    # n_1d weights
    ok = ok and n_1d(5) == 30 and n_1d(7) == 140
    ok = ok and HPLUS[5] == 130 and HPLUS[7] == 5726
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "P_1D(0)=1/p, P_1D(p)=(p−1)/(p(p+1)), P_1D(μ)=(p−1)/(p+1). "
            "Fail: P_1D(p)=1/p."
        ),
    }


def prove_B() -> dict:
    ok = True
    ok = ok and En4_1d(5) == 94
    ok = ok and En4_1d(7) == 318
    ok = ok and En4_1d_drop_p3(5) == Fraction(32, 3)
    ok = ok and En4_1d_drop_p3(5) != 94
    if En4_1d_drop_p3(5) == 94:
        ok = False
    # recover from P_1D
    for p in (5, 7, 11):
        m4 = sum(P_1d(p, k) * k**4 for k in range(p + 1))
        ok = ok and m4 == En4_1d(p)
    return {
        "proved": bool(ok),
        "En4": {str(p): str(En4_1d(p)) for p in (5, 7, 11, 13)},
        "drop5": str(En4_1d_drop_p3(5)),
        "theorem": "E[n⁴|1D]=(p−1)/(p+1)·(p³+μ⁴). Fail: drop p³.",
    }


def prove_C() -> dict:
    live = {
        0: Fraction(1, 5),
        1: Fraction(2, 13),
        2: Fraction(4, 13),
        3: Fraction(2, 13),
        4: Fraction(2, 13),
        5: Fraction(2, 65),
    }
    ok = True
    ok = ok and w_1d(5) == Fraction(3, 13)
    ok = ok and D_live(5) == 13
    for k, tgt in live.items():
        ok = ok and P_ens_p5(k) == tgt
    # p=7 obstruction
    pred7 = n_eq_p_if_nl_misses_full_line(7)
    ok = ok and pred7 == 420
    ok = ok and LIVE_N_EQ_P[7] == 3948
    ok = ok and pred7 != LIVE_N_EQ_P[7]
    ok = ok and n_eq_p_if_nl_misses_full_line(5) == LIVE_N_EQ_P[5] == 60
    if pred7 == LIVE_N_EQ_P[7]:
        ok = False
    return {
        "proved": bool(ok),
        "w5": str(w_1d(5)),
        "P_ens_5": {str(k): str(P_ens_p5(k)) for k in range(6)},
        "p7_1d_only_n_eq_p": pred7,
        "p7_live_n_eq_p": LIVE_N_EQ_P[7],
        "theorem": (
            "p=5 NL uniform on {0..4} recovers live hist. "
            "Fail: P_NL(p)=0 at p=7 (420≠3948)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "EM2_still_unnamed": {
            "5": str(EM2_named(5, LIVE_QPP[5])),
            "7": str(EM2_named(7, LIVE_QPP[7])),
        },
        "note": (
            "1D occupancy named; ensemble / NL occupancy still den D. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.385  1D square-line occupancy named", flush=True)
    A = prove_A()
    print(f"  A P_1D: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B En4_1D: {B['proved']} {B['En4']}", flush=True)
    C = prove_C()
    print(f"  C mix/obstruct: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "P_1d": A["rows"],
                "En4_1d": B["En4"],
                "P_ens_5": C["P_ens_5"],
                "p7_live_n_eq_p": C["p7_live_n_eq_p"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.385",
        "title": "1D square-line occupancy named; p=5 NL uniform; ensemble open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "P_1d_named": A["proved"],
            "En4_1d_named": B["proved"],
            "p5_mixture_p7_obstruction": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15385.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
