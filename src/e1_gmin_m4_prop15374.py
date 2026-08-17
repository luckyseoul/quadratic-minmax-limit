#!/usr/bin/env python3
"""
Prop 15.374 — p≡3 S0 residual at Q++=Q_form.  S0 solves for Q_N*
in terms of Q_{--}.  The floor needs Q_{--}≥Q_m^{min} with
    Q_m^{min}=(p−3)(8+12/(p−1)−Q_form)/(p+1).
At p=7 this is 1273/409.  The named comparison Q++≥Q_N* only
gives 1232/409, which makes sm0<6.  Q_m^{min}≤Q_lo ⇔ Q_form≥Q_lo
(15.373).  Q_{--} is not proved to lie in the floor window.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.373 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  p≡3 (mod 4), p≥7.  |++|=3(p−3)/2, |N*|=(p−1)(p−3)/2,
|--N1|=(p+1)/2, Wick Q(±1)=8.  S0:
    |++|Q++ + |N*|Q_N* + |--|Q_{--} = 2p²−18.
Even S = F_1 ± σ F_* as in 15.314 / 15.302.

============================================================================
Theorem A — PROVED (S0 algebra; certified p=7).
  Pin Q++=Q_form.  Then
      Q_N* = (2p²−18 − |++| Q_form − |--| Q_{--}) / |N*|.
  At p=7 and live Q_{--}=1376/409 this returns 1496/409.
  Fail: use the S0-coarse 1440/409 as Q_N*.  ∎

Theorem B — PROVED (sm0=6; certified p=7).
  sm0≥6 ⇔ Q_{--} ≥ Q_m^{min}(p) with the named formula above.
  At p=7, Q_m^{min}=1273/409.  The comparison Q++≥Q_N* rearranges
  to Q_{--}≥1232/409.  At that point sm0=5832/409<6.
  Fail: claim Q++≥Q_N* forces sm0≥6.  ∎

Theorem C — PROVED (15.373 + algebra).
  Q_m^{min}≤Q_lo ⇔ Q_form≥Q_lo, which holds for every odd p≥5.
  Thus Q_{--}≥Q_lo already implies sm0≥6 at Q++=Q_form.
  Fail: claim Q_lo < Q_m^{min} at p=7.  ∎

Theorem D — OPEN.  The residual Q_{--} is not proved to lie in
  [Q_m^{min}, t_{smh}].  [Q_lo, Q^{ub}] is not a safe upper window
  for smh when p≥11.  phi_F not imported.  D_form untested at p≥11.

============================================================================
Backend: Fraction S0 / F_1 / F_* algebra.  Writes
evidence/e1_gmin_m4_prop15374.json
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
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15372 import Qpp_form  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15374_catalog.json"

LIVE_QMM = {7: Fraction(1376, 409)}
LIVE_QN = {7: Fraction(1496, 409)}
S0_COARSE_QN = {7: Fraction(1440, 409)}


def type_sizes_p3(p: int) -> tuple[int, int, int]:
    if p < 7 or p % 4 != 3:
        raise ValueError("need p≡3, p≥7")
    npp = 3 * (p - 3) // 2
    nN = (p - 1) * (p - 3) // 2
    nmm = (p + 1) // 2
    return npp, nN, nmm


def Qn_from_S0(p: int, Qpp: Fraction, Qm: Fraction) -> Fraction:
    npp, nN, nmm = type_sizes_p3(p)
    rhs = 2 * p * p - 18
    return (rhs - npp * Qpp - nmm * Qm) / nN


def Qm_from_Qpp_ge_Qn(p: int, Qpp: Fraction) -> Fraction:
    """Lower bound on Qm implied by Q++ ≥ Q_N* on the S0 line."""
    npp, nN, nmm = type_sizes_p3(p)
    rhs = 2 * p * p - 18
    # nN Qpp ≥ rhs - npp Qpp - nmm Qm  ⇒  nmm Qm ≥ rhs - (npp+nN) Qpp
    return (rhs - (npp + nN) * Qpp) / nmm


def Qm_min_sm0(p: int, Qpp: Fraction) -> Fraction:
    """Q_{--} at which sm0=6, for this Q++."""
    return Fraction(p - 3, 1) * (8 + Fraction(12, p - 1) - Qpp) / (p + 1)


def even_S_p3(p: int, Qpp: Fraction, Qm: Fraction) -> dict[str, Fraction]:
    Qn = Qn_from_S0(p, Qpp, Qm)
    nppU = Fraction(p - 3, 2)
    nmmU = Fraction(p + 1, 2)
    F1_0 = 16 + nppU * Qpp + nmmU * Qm
    F1_h = 16 + nppU * Qpp - nmmU * Qm
    F1_o = 16 - 2 * Qpp
    Fs_0 = 2 * Qpp + (p - 1) * Qn
    Fs_o = 2 * (Qpp - Qn)
    sig = Fraction(p - 3, 2)
    return {
        "sm0": F1_0 - Fs_0,
        "smh": F1_h - Fs_o,
        "smo": F1_o - Fs_o,
        "sp0": F1_0 + sig * Fs_0,
        "sph": F1_h + sig * Fs_o,
        "spo": F1_o + sig * Fs_o,
        "Qn": Qn,
    }


def prove_A() -> dict:
    p = 7
    Qf = Qpp_form(p)
    Qm = LIVE_QMM[p]
    Qn = Qn_from_S0(p, Qf, Qm)
    ok = Qn == LIVE_QN[p]
    ok = ok and Qn != S0_COARSE_QN[p]
    if Qn == S0_COARSE_QN[p]:
        ok = False
    return {
        "proved": bool(ok),
        "Qn": str(Qn),
        "live_QN": str(LIVE_QN[p]),
        "coarse_QN": str(S0_COARSE_QN[p]),
        "Qpp": str(Qf),
        "theorem": "S0 at Q++=Q_form recovers live Q_N*=1496/409, not 1440/409.",
    }


def prove_B() -> dict:
    p = 7
    Qf = Qpp_form(p)
    tmin = Qm_min_sm0(p, Qf)
    tpair = Qm_from_Qpp_ge_Qn(p, Qf)
    S_pair = even_S_p3(p, Qf, tpair)
    S_min = even_S_p3(p, Qf, tmin)
    S_live = even_S_p3(p, Qf, LIVE_QMM[p])
    ok = tmin == Fraction(1273, 409)
    ok = ok and tpair == Fraction(1232, 409)
    ok = ok and tpair < tmin
    ok = ok and S_pair["sm0"] < 6
    ok = ok and S_min["sm0"] == 6
    ok = ok and S_live["sm0"] == Fraction(3072, 409)
    ok = ok and S_live["sm0"] >= 6
    if S_pair["sm0"] >= 6 or tpair >= tmin:
        ok = False
    return {
        "proved": bool(ok),
        "Qm_min": str(tmin),
        "Qm_Qpp_ge_Qn": str(tpair),
        "gap": str(tmin - tpair),
        "sm0_at_pair": str(S_pair["sm0"]),
        "sm0_at_min": str(S_min["sm0"]),
        "sm0_live": str(S_live["sm0"]),
        "theorem": "Q++≥Q_N* misses the sm0=6 cut by 41/409.",
    }


def prove_C() -> dict:
    # Qm_min ≤ Q_lo ⇔ Q_form ≥ Q_lo (15.373)
    ok = True
    rows = []
    for p in (7, 11, 19, 23, 31, 43):
        Qf = Qpp_form(p)
        tmin = Qm_min_sm0(p, Qf)
        lo = Q_lo_named(p)
        ok = ok and tmin <= lo
        ok = ok and Qf >= lo
        rows.append((p, str(tmin), str(lo)))
    # p=7 explicit
    ok = ok and Qm_min_sm0(7, Qpp_form(7)) <= Q_lo_named(7)
    if Qm_min_sm0(7, Qpp_form(7)) > Q_lo_named(7):
        ok = False
    return {
        "proved": bool(ok),
        "p7_Qmmin_le_Qlo": Qm_min_sm0(7, Qpp_form(7)) <= Q_lo_named(7),
        "rows": rows,
        "theorem": "Qm_min≤Q_lo ⇔ Q_form≥Q_lo, already proved.",
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
            "Q_{--} is not proved to lie in [Qm_min, t_smh]. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.374  p≡3 residual at Q++=Q_form", flush=True)
    A = prove_A()
    print(f"  A S0: {A['proved']} Qn={A['Qn']}", flush=True)
    B = prove_B()
    print(f"  B cut: {B['proved']} Qmmin={B['Qm_min']} pair={B['Qm_Qpp_ge_Qn']}", flush=True)
    C = prove_C()
    print(f"  C Qlo: {C['proved']} p7={C['p7_Qmmin_le_Qlo']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "Qm_min_7": str(Qm_min_sm0(7, Qpp_form(7))),
                "pair_7": str(Qm_from_Qpp_ge_Qn(7, Qpp_form(7))),
                "live_Qmm": str(LIVE_QMM[7]),
                "live_QN": str(LIVE_QN[7]),
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.374",
        "title": "p≡3 S0 residual: Qm_min named; Q++≥QN* misses sm0=6",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S0_recovers_live_QN": A["proved"],
            "Qm_min_named_pair_misses": B["proved"],
            "Qm_min_le_Qlo": C["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15374.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
