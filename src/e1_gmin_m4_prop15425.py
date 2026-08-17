#!/usr/bin/env python3
"""
Prop 15.425 — N* vertical φ is ×k for k∈F_7^×\\{3,4}={1,2,5,6},
the complement of the ++ maps {×3,×4}.  Opposite-orient R×R
under N* φ equals L_par(R)=133.  R×C with relative ×2 equals
238/3.  RRRC same-pair N* is (133+238/3)/2=637/6.
The full NL mix 12349/114 still needs O orientations.
Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.424 flags.

============================================================================
Setup.  15.420: ++ extras {21,28} are vertical with φ=×3,×4.
N* extras include vertical {7,14,35,42} and slope-3/4 eights.
Field mul by those vertical N* sends slope-0 param a to
slope-∞ param k a, k∈{1,2,5,6}.  R flags occur in three
even orientations J∘×1,×2,×3; live RRRC axes are opposite.

============================================================================
Theorem A — PROVED (field mul; Max+-free).
  Vertical N* φ = ×{1,2,5,6} = F_p^× \\ {3,4}.
  Fail: claim N* uses ++ φ {×3,×4}.  ∎

Theorem B — PROVED (4-point; certified live 133).
  Opposite-orient R×R (reverse × reverse∘×2) under N* φ
  equals L_par(R)=133.  Fail: same-orient (413/3);
  fail: ++ φ on opposite orients (427/3).  ∎

Theorem C — PROVED (4-point; certified live 238/3).
  R×C under N* φ with C relabeled ×2 equals 238/3.
  Fail: relative id (245/3); relative ×3 (84).  ∎

Theorem D — PROVED (tag mix; certified 637/6).
  RRRC same-pair N* = (133+238/3)/2 = 637/6.
  Fail: ignore C tags (133).  ∎

Theorem E — OPEN.  Full NL L_N* same-pair 12349/114
  and cross 6269/57 need O-type orientations.
  --N1 and Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction 4-point + 15.416 / 15.422 J.  Writes
evidence/e1_gmin_m4_prop15425.json
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
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15416 import Lpar_R_flag  # noqa: E402
from e1_gmin_m4_prop15420 import extra_pp  # noqa: E402
from e1_gmin_m4_prop15422 import J_C_11222, J_R_reverse  # noqa: E402
from e1_gmin_m4_prop15424 import LIVE_LN_SAME  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15425_catalog.json"

P = 7
KS_NSTAR = (1, 2, 5, 6)
KS_PP = (3, 4)


def slope(x: int, p: int = P) -> str:
    x0 = x % p
    if x == 0:
        return "0"
    if x0 == 0:
        return "inf"
    return str(((x // p) * pow(x0, -1, p)) % p)


def nstar_vertical(p: int = P) -> list[int]:
    F = _kernel_field(p)
    return [
        s
        for s in range(1, p * p)
        if F["chi"][s] == 1
        and coarse_Q_class(F, s) == "N*"
        and slope(s, p) == "inf"
    ]


def phi_k_vertical(p: int = P) -> list[int]:
    """k such that vertical N* sends horiz a to vert k a."""
    F = _kernel_field(p)
    ks = []
    for s in nstar_vertical(p):
        # a=1 → param of s*1
        sd = F["mul"][s][1]
        ks.append(sd // p)
    return sorted(ks)


def L_pair_ks(JA, JB, ks, p: int = P) -> Fraction:
    acc = 0
    c = 0
    for a in range(1, p):
        for k in ks:
            b = (k * a) % p
            if b == 0:
                raise RuntimeError("k a = 0")
            acc += JA(a) * JB(b)
            c += 1
    return Fraction(acc, c)


def J_R_ferrers_via(a: int) -> int:
    return J_R_reverse((2 * a) % P)


def J_C_x2(a: int) -> int:
    return J_C_11222((2 * a) % P)


def J_C_x3(a: int) -> int:
    return J_C_11222((3 * a) % P)


def L_RR_Nstar_opp() -> Fraction:
    return L_pair_ks(J_R_reverse, J_R_ferrers_via, KS_NSTAR)


def L_RR_Nstar_same_wrong() -> Fraction:
    return L_pair_ks(J_R_reverse, J_R_reverse, KS_NSTAR)


def L_RR_opp_pp_wrong() -> Fraction:
    return L_pair_ks(J_R_reverse, J_R_ferrers_via, KS_PP)


def L_RC_Nstar() -> Fraction:
    return L_pair_ks(J_R_reverse, J_C_x2, KS_NSTAR)


def L_RC_Nstar_id_wrong() -> Fraction:
    return L_pair_ks(J_R_reverse, J_C_11222, KS_NSTAR)


def L_RC_Nstar_x3_wrong() -> Fraction:
    return L_pair_ks(J_R_reverse, J_C_x3, KS_NSTAR)


def L_RRRC_Nstar_same() -> Fraction:
    return (L_RR_Nstar_opp() + L_RC_Nstar()) / 2


def L_RRRC_Nstar_noC_wrong() -> Fraction:
    return L_RR_Nstar_opp()


def prove_A() -> dict:
    ks = phi_k_vertical()
    pp = extra_pp(7)
    ok = ks == list(KS_NSTAR)
    ok = ok and set(ks).isdisjoint(KS_PP)
    ok = ok and set(ks) | set(KS_PP) == set(range(1, P))
    ok = ok and set(pp) == {21, 28}
    if ks == list(KS_PP):
        ok = False
    return {
        "proved": bool(ok),
        "Nstar_k": ks,
        "pp_k": list(KS_PP),
        "theorem": (
            "Vertical N* φ=×{1,2,5,6}=F_7^×\\{3,4}. "
            "Fail: N* uses ++ {×3,×4}."
        ),
    }


def prove_B() -> dict:
    form = L_RR_Nstar_opp()
    same = L_RR_Nstar_same_wrong()
    ppopp = L_RR_opp_pp_wrong()
    ok = form == Lpar_R_flag(7) == 133
    ok = ok and same == Fraction(413, 3) and same != form
    ok = ok and ppopp == Fraction(427, 3) and ppopp != form
    if same == 133 or ppopp == 133 and form != 133:
        pass
    if same == form:
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "same_orient": str(same),
        "pp_opp": str(ppopp),
        "theorem": (
            "Opposite R×R under N* φ is L_par(R)=133. "
            "Fail: same-orient 413/3; ++ φ 427/3."
        ),
    }


def prove_C() -> dict:
    form = L_RC_Nstar()
    idn = L_RC_Nstar_id_wrong()
    x3 = L_RC_Nstar_x3_wrong()
    ok = form == Fraction(238, 3)
    ok = ok and idn == Fraction(245, 3) and idn != form
    ok = ok and x3 == 84 and x3 != form
    if idn == form:
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "id": str(idn),
        "x3": str(x3),
        "theorem": (
            "R×C N* relative ×2 is 238/3. Fail: id 245/3; ×3 84."
        ),
    }


def prove_D() -> dict:
    mix = L_RRRC_Nstar_same()
    noc = L_RRRC_Nstar_noC_wrong()
    ok = mix == Fraction(637, 6)
    ok = ok and noc == 133 and noc != mix
    if noc == mix:
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "noC": str(noc),
        "theorem": (
            "RRRC same-pair N*=(133+238/3)/2=637/6. Fail: 133."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "LN_same": str(LIVE_LN_SAME),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "N* φ and R/C slots named. Full mix 12349/114 needs O "
            "orientations. Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.425  N* φ=×{1,2,5,6}; R/C slots named", flush=True)
    A = prove_A()
    print(f"  A phi: {A['proved']} k={A['Nstar_k']}", flush=True)
    B = prove_B()
    print(f"  B RR: {B['proved']} {B['form']} same={B['same_orient']}", flush=True)
    C = prove_C()
    print(f"  C RC: {C['proved']} {C['form']}", flush=True)
    D = prove_D()
    print(f"  D RRRC: {D['proved']} {D['mix']}", flush=True)
    E = prove_open()
    print(f"  E open: QNL>Q1d={E['Q_NL_gt_Q1d_p5']} phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": {"k": A["Nstar_k"]}, "D": {"mix": D["mix"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.425",
        "title": (
            "N* φ=×{1,2,5,6}; opposite R×R is 133; R×C ×2 is 238/3; "
            "RRRC same-pair 637/6; Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Nstar_phi_1256": A["proved"],
            "RR_opp_133": B["proved"],
            "RC_x2_238_3": C["proved"],
            "RRRC_637_6": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": E["Q_NL_is_short_rational"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15425.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
