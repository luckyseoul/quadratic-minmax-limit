#!/usr/bin/env python3
"""
Prop 15.375 — Affine constructors vanish on --, and
    Q_{--}^{form} := 8(A−k)/D_form
equals live 1376/409 at p=7 and lies in [Q_m^{min}, t_{smh}]
for every prime p≡3 (mod 4), p≥7.  All six even S0 families
are then ≥6 at (Q++,Q--)=(Q_form, Q_{--}^{form}).  The form
is not identified with live Q_{--} at p≥11.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.374 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  p≡3, p≥7.  -- is entirely --∩N1 (15.355), and
N1∩F_p={±1}.  15.292 A: Q_1D(r)=0 for r∉F_p.  A, D_form, k
from 15.367 / 15.371.  Q_lo=2(2p+1)(p−3)/(p−1)² (15.326).
Qm_min≤Q_lo (15.374 C).  C=C(p,(p−1)/2) ≥ C(p,4) (p=7 equal).

============================================================================
Theorem A — PROVED (15.292 A + 15.355 torus; certified p=7).
  Q_1D^{--}=Q_AP^{--}=Q_QR0^{--}=0.  Fail: claim
  Q_QR0^{--}=4(p+1)/3 (32/3≠0 at p=7).  Hence
      Q_{--}=(1−n_1d/N) Q_NL^{--}.  ∎

Theorem B — PROVED (named form; certified p=7).
  8(A−k)/D_form = 1376/409.  Drop-k is 1544/409≠live.
  8(A−n++)/D_form = 1496/409 is live Q_N*, not Q_{--}.
  At p=3 both sides of the form are 0 (S0 degeneracy).
  Fail: claim the drop-k value, or claim 1496/409 is Q_{--}.  ∎

Theorem C — PROVED (C-coefficients; all p≡3≥7).
  Q_{--}^{form}≥Q_lo ⇔ C·(17p+13)+N≥0 with
      N=(−4p⁴+12p³−69p²−38p+351)/2.
  dC>0, −N_num(7)=8784, and (−N_num)′>0.  The C(p,4)
  bound is the integer
      f_lo=34p⁵−274p⁴+506p³−1574p²−1068p+8424,
  with every derivative positive at p=7, hence f_lo>0
  for all real p≥7.  Fail: N_num(7)>0 (−8784); f_lo(7)=10944.
  Q_{--}^{form}≤Q++_form identically (k>0).  With 15.374 C,
  Qm_min≤Q_lo≤Q_{--}^{form} so sm0≥6 at Q++=Q_form.  ∎

Theorem D — PROVED (C-coefficients; all p≡3≥7).
  Q++_form≤t_smh ⇔ (9p+7)(p−3)D≥16A(p²−2p−1).
  dC_tsh=2p³+3p²−54p−47, value 408 at p=7, derivative
  6(p²+p−9)>0.  The C(p,4) bound f_tsh has every
  derivative positive at p=7.  Fail: dC_tsh(7)=0 (408).
  Then Q_{--}^{form}≤Q++_form≤t_smh ⇒ smh≥6, and
  smo>smh(Q++,Q++)≥6 (Qn increases when Qm drops).
  sph and spo have dC=84p+108 and 42p+54; their C(p,4)
  bounds have every derivative positive at p=7.  Fail:
  dC_spo(7)=0 (348).  All six even families ≥6.  ∎

Theorem E — OPEN.  The p=7 identity does not name live
  Q_{--} at p≥11.  D_form untested at p≥11.  phi_F not
  imported.

============================================================================
Backend: integer C-coefficient polynomials.  Writes
evidence/e1_gmin_m4_prop15375.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15355 import Q_QR0_p3_named, Qpp_of_y  # noqa: E402
from e1_gmin_m4_prop15367 import A_num, k_named  # noqa: E402
from e1_gmin_m4_prop15371 import D_form  # noqa: E402
from e1_gmin_m4_prop15372 import Qpp_form  # noqa: E402
from e1_gmin_m4_prop15374 import (  # noqa: E402
    LIVE_QMM,
    LIVE_QN,
    Qm_min_sm0,
    even_S_p3,
    type_sizes_p3,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15375_catalog.json"


def Qmm_form(p: int) -> Fraction:
    return Fraction(8 * (A_num(p) - k_named(p)), D_form(p))


def Qn_form_hit(p: int) -> Fraction:
    """p=7 interpolant 8(A−n++)/D; not a general name."""
    from e1_gmin_m4_prop15355 import n_pp_named

    return Fraction(8 * (A_num(p) - n_pp_named(p)), D_form(p))


def t_smh(p: int, Qpp: Fraction) -> Fraction:
    npp, nN, nmm = type_sizes_p3(p)
    nppU = Fraction(p - 3, 2)
    nmmU = Fraction(p + 1, 2)
    rhs = 2 * p * p - 18
    num = (
        10
        + (nppU - 2) * Qpp
        + Fraction(2 * rhs, nN)
        - Fraction(2 * npp, nN) * Qpp
    )
    den = nmmU + Fraction(2 * nmm, nN)
    return num / den


def dC_lo_mm(p: int) -> int:
    return 17 * p + 13


def N_num_lo(p: int) -> int:
    return -4 * p**4 + 12 * p**3 - 69 * p**2 - 38 * p + 351


def f_lo(p: int) -> int:
    return 34 * p**5 - 274 * p**4 + 506 * p**3 - 1574 * p**2 - 1068 * p + 8424


def dC_tsh(p: int) -> int:
    return 2 * p**3 + 3 * p**2 - 54 * p - 47


def dC_spo(p: int) -> int:
    return 42 * p + 54


def dC_sph(p: int) -> int:
    return 84 * p + 108


def f_tsh(p: int) -> int:
    # C4n * dC_tsh + 24 N_tsh; monomial verified by Newton
    return (
        2 * p**7
        - 9 * p**6
        - 50 * p**5
        + 226 * p**4
        - 582 * p**3
        + 2867 * p**2
        + 1110 * p
        - 10476
    )


def f_spo(p: int) -> int:
    return 42 * p**5 - 198 * p**4 + 42 * p**3 - 1314 * p**2 - 2820 * p + 11160


def f_sph(p: int) -> int:
    return 84 * p**5 - 300 * p**4 - 12 * p**3 - 2724 * p**2 - 5544 * p + 22320


def _derivs_at(mon: list[int], x: int) -> list[int]:
    out = []
    for k in range(len(mon)):
        s = 0
        for j in range(k, len(mon)):
            coef = 1
            for t in range(k):
                coef *= j - t
            s += coef * mon[j] * x ** (j - k)
        out.append(s)
    return out


def _all_derivs_pos(mon: list[int], x: int = 7) -> bool:
    return all(v > 0 for v in _derivs_at(mon, x))


def _constructor_Q_cls(p: int, which: str, cls: str) -> Fraction:
    F = _kernel_field(p)
    Omega = omega_set(F)
    rs = [r for r in F["squares"] if merge_cls(F, r) == cls]
    if which == "AP":
        y = affine_halfspace(p, (0, 1), set(range((p + 1) // 2)))
    else:
        y = affine_halfspace(p, (0, 1), qr0(p))
    if not rs:
        return Fraction(0)
    return Qpp_of_y(p, y, F, Omega, rs)


def prove_A() -> dict:
    p = 7
    q_ap = _constructor_Q_cls(p, "AP", "--N1")
    q_qr = _constructor_Q_cls(p, "QR0", "--N1")
    q_qr_pp = _constructor_Q_cls(p, "QR0", "++")
    ok = q_ap == 0 and q_qr == 0
    ok = ok and q_qr_pp == Q_QR0_p3_named(p)
    if q_qr == Q_QR0_p3_named(p) or q_qr != 0:
        ok = False
    w = Fraction(n_1d(p), 2 * p * D_form(p))
    qnl = LIVE_QMM[p] / (1 - w)
    return {
        "proved": bool(ok),
        "AP_mm": str(q_ap),
        "QR0_mm": str(q_qr),
        "QR0_pp": str(q_qr_pp),
        "QR0_pp_named": str(Q_QR0_p3_named(p)),
        "w_1d": str(w),
        "Q_NL_mm_p7": str(qnl),
        "theorem": (
            "Affine Type+ constructors vanish on --. "
            "Fail: Q_QR0^{--}=4(p+1)/3."
        ),
    }


def prove_B() -> dict:
    p = 7
    form = Qmm_form(p)
    dropk = Qpp_form(p)
    qn_hit = Qn_form_hit(p)
    ok = form == LIVE_QMM[p] == Fraction(1376, 409)
    ok = ok and dropk != LIVE_QMM[p]
    ok = ok and qn_hit == LIVE_QN[p]
    ok = ok and qn_hit != LIVE_QMM[p]
    # p=3: A-k=0 so form=0
    A3k = A_num(3) - k_named(3)
    ok = ok and A3k == 0
    if dropk == LIVE_QMM[p] or form != LIVE_QMM[p]:
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "dropk": str(dropk),
        "Qn_hit": str(qn_hit),
        "A_minus_k_p3": A3k,
        "theorem": (
            "8(A−k)/D_form hits live Q_{--}. "
            "Drop-k and 8(A−n++)/D are not Q_{--}."
        ),
    }


def prove_C() -> dict:
    ok = True
    ok = ok and dC_lo_mm(7) == 132
    ok = ok and N_num_lo(7) == -8784
    ok = ok and f_lo(7) == 10944
    if N_num_lo(7) > 0 or f_lo(7) <= 0:
        ok = False
    # k = |regular set| = p(p-1)/2 ≥ 0; A-k = (p+1)C - 2p²-4p+18
    for p in (7, 11, 19):
        ok = ok and k_named(p) == p * (p - 1) // 2
        C = comb(p, (p - 1) // 2)
        Ak_const = -2 * p * p - 4 * p + 18
        ok = ok and A_num(p) - k_named(p) == (p + 1) * C + Ak_const
        if k_named(p) < 0:
            ok = False
    mon = [8424, -1068, -1574, 506, -274, 34]
    ok = ok and _all_derivs_pos(mon, 7)
    # recon: gap = C dC + N_num/2 equals 4(A-k)(p-1)^2 - D(2p+1)(p-3)
    rows = []
    for p in (7, 11, 19, 23, 31, 43):
        C = comb(p, (p - 1) // 2)
        N = Fraction(N_num_lo(p), 2)
        gap_C = C * dC_lo_mm(p) + N
        Ak = A_num(p) - k_named(p)
        gap_Q = Fraction(
            4 * Ak * (p - 1) ** 2 - D_form(p) * (2 * p + 1) * (p - 3)
        )
        ok = ok and gap_C == gap_Q
        ok = ok and gap_Q >= 0
        ok = ok and Qmm_form(p) >= Q_lo_named(p)
        ok = ok and Qmm_form(p) <= Qpp_form(p)
        ok = ok and Qm_min_sm0(p, Qpp_form(p)) <= Q_lo_named(p)
        S = even_S_p3(p, Qpp_form(p), Qmm_form(p))
        ok = ok and S["sm0"] >= 6
        rows.append((p, str(gap_Q), str(S["sm0"])))
        if gap_Q < 0:
            ok = False
    return {
        "proved": bool(ok),
        "dC": "17p+13",
        "N_num_7": N_num_lo(7),
        "f_lo_7": f_lo(7),
        "derivs_pos": _all_derivs_pos(mon, 7),
        "rows": rows,
        "theorem": (
            "Q_{--}^{form}≥Q_lo for all p≡3≥7 by C(p,4); "
            "sm0≥6 at Q++=Q_form."
        ),
    }


def prove_D() -> dict:
    ok = True
    ok = ok and dC_tsh(7) == 408
    ok = ok and dC_spo(7) == 348
    ok = ok and dC_sph(7) == 696
    if dC_tsh(7) == 0 or dC_spo(7) == 0:
        ok = False
    mon_t = [-10476, 1110, 2867, -582, 226, -50, -9, 2]
    mon_o = [11160, -2820, -1314, 42, -198, 42]
    mon_h = [22320, -5544, -2724, -12, -300, 84]
    ok = ok and _all_derivs_pos(mon_t, 7)
    ok = ok and _all_derivs_pos(mon_o, 7)
    ok = ok and _all_derivs_pos(mon_h, 7)
    ok = ok and f_tsh(7) == 228672
    ok = ok and f_spo(7) == 171936
    ok = ok and f_sph(7) == 537408
    rows = []
    for p in (7, 11, 19, 23, 31, 43, 79):
        Qf = Qpp_form(p)
        Qm = Qmm_form(p)
        tsh = t_smh(p, Qf)
        ok = ok and Qf <= tsh
        ok = ok and Qm <= Qf
        S = even_S_p3(p, Qf, Qm)
        mins = min(S[k] for k in S if k != "Qn")
        ok = ok and mins >= 6
        ok = ok and S["smh"] >= 6 and S["smo"] >= 6
        ok = ok and S["sph"] >= 6 and S["spo"] >= 6
        ok = ok and S["sp0"] == 2 * (p * p - 1)
        # Qn = (rhs - npp Qpp - nmm Qm)/nN decreases in Qm (nmm, nN > 0)
        S_eq = even_S_p3(p, Qf, Qf)
        ok = ok and S["Qn"] > S_eq["Qn"]
        ok = ok and S["smo"] > S_eq["smh"]
        # tsh cross-multiply
        A, D = A_num(p), D_form(p)
        cross = (9 * p + 7) * (p - 3) * D - 16 * A * (p * p - 2 * p - 1)
        ok = ok and cross >= 0
        rows.append(
            (
                p,
                str(tsh),
                str(S["smh"]),
                str(S["spo"]),
                str(mins),
            )
        )
        if mins < 6 or Qf > tsh:
            ok = False
    return {
        "proved": bool(ok),
        "dC_tsh_7": dC_tsh(7),
        "dC_spo_7": dC_spo(7),
        "dC_sph_7": dC_sph(7),
        "f_tsh_7": f_tsh(7),
        "f_spo_7": f_spo(7),
        "f_sph_7": f_sph(7),
        "derivs_pos": True,
        "rows": rows,
        "theorem": (
            "Q++_form≤t_smh and sph,spo≥6 for all p≡3≥7 "
            "by C(p,4) bounds with positive jets at p=7."
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
        "Qmm_named_in_p": False,
        "note": (
            "Q_{--}^{form} sits in the floor window if D_form=D "
            "and live Q_{--} equals the form.  Neither is known "
            "at p≥11.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.375  Q--_form in the p≡3 floor window", flush=True)
    A = prove_A()
    print(f"  A affine 0: {A['proved']} QR0_mm={A['QR0_mm']}", flush=True)
    B = prove_B()
    print(f"  B form: {B['proved']} {B['form']}", flush=True)
    C = prove_C()
    print(f"  C Qlo: {C['proved']} Nnum7={C['N_num_7']}", flush=True)
    D = prove_D()
    print(f"  D tsh/S: {D['proved']} dCtsh7={D['dC_tsh_7']}", flush=True)
    E = prove_open()
    print(f"  E open: phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "Qmm_form": {str(p): str(Qmm_form(p)) for p in (7, 11, 19, 23)},
                "Qpp_form": {str(p): str(Qpp_form(p)) for p in (7, 11, 19, 23)},
                "dC_lo": "17p+13",
                "dC_tsh": "2p^3+3p^2-54p-47",
                "dC_spo": "42p+54",
                "dC_sph": "84p+108",
                "f_lo_7": f_lo(7),
                "f_tsh_7": f_tsh(7),
                "live_Qmm": str(LIVE_QMM[7]),
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.375",
        "title": (
            "Q--_form=8(A-k)/D in [Qm_min,t_smh]; "
            "p≡3 even S≥6 at the form pair"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "affine_Qmm_zero": A["proved"],
            "form_hits_live_p7": B["proved"],
            "Qmm_form_ge_Qlo": C["proved"],
            "even_S_ge_6_at_form_pair": D["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
            "Qmm_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15375.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
