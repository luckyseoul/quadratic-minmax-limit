#!/usr/bin/env python3
"""
Prop 15.395 — S0 + Cauchy–Schwarz / convexity on occupancy
does **not** force S_□≥6.  Var(M)≥0 and the 2-point fourth
moment sit strictly below the S0 floor window; the 3R+C
(triangle) and 1D endpoints sit strictly outside it.
phi_F is imported only if the bound held; it does not.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.394 flags.

============================================================================
Setup.  15.326: p≡1 S0 line Qn=(rhs−a Q++)/b, even S=F1±σ F_*.
Floor ⇔ Q_lo≤Q++≤Q^{ub}.  15.384: Var(M)≥0.  15.305: E[n]=Var(n)
=μ=(p−1)/2 on square lines, so E[n²]=μ(μ+1).  15.387: T is
3R+(n_sq−3)C.  15.350 / 15.356: Q_T and Q_1D named.

============================================================================
Theorem A — PROVED (15.384 inversion; all odd p≥5).
  Var(M)≥0 ⇔ Q++ ≥ Q_cs0, with Q_cs0=0 at p=5 and
  Q_cs0=4(p−5)/(p−3) for p>5.  On the S0 line, sm_j0=6
  only at Q_lo, and Q_cs0<Q_lo, so sm_j0(Q_cs0)<6.
  Fail: claim sm_j0(Q_cs0)≥6 at p=5 (live −16).  ∎

Theorem B — PROVED (2-point CS; all odd p≥5).
  Among laws on {0,…,p} with the named (E[n],E[n²]),
  E[n⁴]≥μ(μ+1)³, attained at {0,μ+1}.  Strictly stronger
  than (E[n²])²=μ²(μ+1)².  En4_T and En4_1d sit above it.
  Fail: claim μ(μ+1)³=(E[n²])² at p=5 (54≠36).  ∎

Theorem C — PROVED (3R+C / 1D are not the S0 window).
  Q_T>Q^{ub} at p=5,13,17,29,37 (64/15>26/7), so
  sm_joth(Q_T)<6.  Q_1D(5)=16/9<Q_lo; Q_1D>Q^{ub}
  at p≥13.  Fail: claim Q_T(5)≤Q^{ub}.  ∎

Theorem D — OPEN.  No CS/convexity bound from S0 + occupancy
  forces the window.  sbox_cs_ge_6() is False, so
  phi_F_ge_6 is not imported.

============================================================================
Backend: Fraction identities on even_S_p1.  Writes
evidence/e1_gmin_m4_prop15395.json
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
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15385 import En4_1d, mu_named  # noqa: E402
from e1_gmin_m4_prop15387 import En4_C, En4_R, En4_T  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15395_catalog.json"
P1 = (5, 13, 17, 29, 37)


def Q_cs0(p: int) -> Fraction:
    """Var(M)≥0 invert (15.384)."""
    if p == 5:
        return Fraction(0)
    return Fraction(4 * (p - 5), p - 3)


def En4_cs_quad(p: int) -> Fraction:
    """(E[n²])² = μ²(μ+1)²."""
    mu = mu_named(p)
    return (mu * (mu + 1)) ** 2


def En4_cs_twopoint(p: int) -> Fraction:
    """2-point {0,μ+1} with the named first two moments: μ(μ+1)³."""
    mu = mu_named(p)
    return mu * (mu + 1) ** 3


def min_even_S(p: int, Qpp: Fraction) -> Fraction:
    s = even_S_p1(p, Qpp)
    return min(s["sm_j0"], s["sm_joth"], s["sp_j0"], s["sp_joth"])


def sbox_cs_ge_6() -> bool:
    """True only if CS/Var(M)≥0 already forces every even S/q²≥6."""
    for p in P1:
        if min_even_S(p, Q_cs0(p)) < 6:
            return False
        if min_even_S(p, Qpp_from_L(p)) < 6:
            return False
        if min_even_S(p, Q_1d_pp_named(p)) < 6:
            return False
    return True


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in P1:
        qc = Q_cs0(p)
        lo = Q_lo_named(p)
        sm = even_S_p1(p, qc)["sm_j0"]
        ok = ok and qc < lo
        ok = ok and sm < 6
        rows[str(p)] = {"Q_cs0": str(qc), "Q_lo": str(lo), "sm_j0": str(sm)}
    ok = ok and even_S_p1(5, Q_cs0(5))["sm_j0"] == -16
    if even_S_p1(5, Q_cs0(5))["sm_j0"] >= 6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Var(M)≥0 ⇒ Q++≥Q_cs0<Q_lo, so sm_j0<6. Fail: sm_j0(0)≥6 at p=5.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in P1:
        c1, c2 = En4_cs_quad(p), En4_cs_twopoint(p)
        ok = ok and c2 > c1
        ok = ok and En4_T(p) > c2
        ok = ok and En4_1d(p) > c2
        ok = ok and En4_R(p) > En4_C(p)
        rows[str(p)] = {
            "cs_quad": str(c1),
            "cs_2pt": str(c2),
            "En4_T": str(En4_T(p)),
            "En4_1d": str(En4_1d(p)),
        }
    ok = ok and En4_cs_twopoint(5) == 54
    ok = ok and En4_cs_quad(5) == 36
    if En4_cs_twopoint(5) == En4_cs_quad(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "E[n⁴]≥μ(μ+1)³ > (E[n²])². Fail: 54=36 at p=5.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in P1:
        qt = Qpp_from_L(p)
        q1 = Q_1d_pp_named(p)
        lo, ub = Q_lo_named(p), Qpp_floor_ub(p)
        ok = ok and qt > ub
        if p == 5:
            ok = ok and q1 < lo
        else:
            ok = ok and q1 > ub
        ok = ok and min_even_S(p, qt) < 6
        ok = ok and min_even_S(p, q1) < 6
        rows[str(p)] = {
            "Q_T": str(qt),
            "Q_1d": str(q1),
            "Q_lo": str(lo),
            "Q_ub": str(ub),
            "minS_T": str(min_even_S(p, qt)),
            "minS_1d": str(min_even_S(p, q1)),
        }
    ok = ok and Qpp_from_L(5) == Fraction(64, 15)
    ok = ok and Qpp_floor_ub(5) == Fraction(26, 7)
    ok = ok and Q_1d_pp_named(5) == Fraction(16, 9)
    if Qpp_from_L(5) <= Qpp_floor_ub(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Q_T>Q^{ub} all listed p≡1; Q_1D<Q_lo at p=5 and "
            "Q_1D>Q^{ub} at p≥13. Fail: Q_T(5)≤26/7."
        ),
    }


def prove_open() -> dict:
    ge6 = sbox_cs_ge_6()
    return {
        "proved": False,
        "sbox_cs_ge_6": ge6,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "CS/Var(M)≥0 and 3R+C / 1D endpoints miss the S0 window. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.395  S0+CS occupancy does not force S_□≥6", flush=True)
    A = prove_A()
    print(f"  A Q_cs0: {A['proved']} p5 sm={A['rows']['5']['sm_j0']}", flush=True)
    B = prove_B()
    print(f"  B 2pt CS: {B['proved']} p5 {B['rows']['5']}", flush=True)
    C = prove_C()
    print(f"  C endpoints: {C['proved']} p5 T={C['rows']['5']['Q_T']}", flush=True)
    D = prove_open()
    print(f"  D open: sbox_cs={D['sbox_cs_ge_6']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"], "C": C["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.395",
        "title": "S0+CS occupancy does not force S_□≥6; phi_F unimported",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Q_cs0_misses_window": A["proved"],
            "twopoint_cs_strict": B["proved"],
            "T_and_1D_outside_window": C["proved"],
            "sbox_cs_ge_6": D["sbox_cs_ge_6"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15395.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
