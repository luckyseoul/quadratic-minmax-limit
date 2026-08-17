#!/usr/bin/env python3
"""
Prop 15.427 — L_--N1 is the 15.426 CROSS Hoffman mix under
++ φ=×{3,4}.  Live 6131/57 at p=7; empty at p=5.
Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.426 flags.

============================================================================
Setup.  15.314: --N1 ⊂ U, empty if p≡1 (mod 4), |--N1|=(p+1)/2
if p≡3.  Live extras {9,12,44,47} on slopes {3,4}.  15.426
named N* CROSS 4-points under N* φ=×{1,2,5,6} (mix 6269/57).
++ extra φ is ×{3,4}.  --N1 extras send slope 0 to {3,4},
i.e. they are the complementary-axis pairing, not 3↔4.

============================================================================
Theorem A — PROVED (field geometry; Max+-free χ / Paley type).
  extra --N1 = {s∈T: coarse=--N1}.  Empty at p=5.
  {9,12,44,47} at p=7, slopes {3,4}, size (p+1)/2.
  Fail: claim extras are vertical p F_p^× (those are ++/N*);
  fail: nonempty at p=5.  ∎

Theorem B — PROVED (field mul; Max+-free).
  Each extra sends slope 0 to {3,4}, not to {inf}.
  Fail: claim --N1 extras are vertical (dest {inf}).  ∎

Theorem C — PROVED (Hoffman 4+4+4+4+2+1 over 19).
  The 15.426 CROSS type 4-points under ++ φ=×{3,4}:
      219=671=325/3, RRRC=105, 665=1309/12,
      344=213/2, 695=323/3.
  Mix = 6131/57 = L_--N1(7).
  Fail: N* φ on the same slots (6269/57);
  fail: ++ φ on SAME slots (673/6 = L_spl).  ∎

Theorem D — OPEN.  Q_NL unnamed in p.  phi_F not imported.

============================================================================
Backend: Fraction 4-point.  Writes
evidence/e1_gmin_m4_prop15427.json
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
from e1_gmin_m4_prop15314 import named_sizes  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL  # noqa: E402
from e1_gmin_m4_prop15420 import extra_pp, vert_units  # noqa: E402
from e1_gmin_m4_prop15422 import (  # noqa: E402
    J_0022566,
    J_0034455,
    J_0112557,
    J_0223347,
    J_1122366ap,
    J_1122366sg,
    J_1122447,
    J_C_11222,
    J_F1,
    J_F3,
    J_F4,
    J_G1,
    J_G3,
    J_G4,
    J_R_reverse,
)
from e1_gmin_m4_prop15423 import J_G2  # noqa: E402
from e1_gmin_m4_prop15424 import LIVE_LN_CROSS  # noqa: E402
from e1_gmin_m4_prop15425 import KS_NSTAR, KS_PP, L_pair_ks  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15427_catalog.json"

P = 7
LIVE_LM1 = {5: None, 7: Fraction(6131, 57)}
EXTRAS7 = (9, 12, 44, 47)


def slope(x: int, p: int = P) -> str:
    x0 = x % p
    if x == 0:
        return "0"
    if x0 == 0:
        return "inf"
    return str(((x // p) * pow(x0, -1, p)) % p)


def extra_mm(p: int) -> list[int]:
    F = _kernel_field(p)
    q = p * p
    return [
        s
        for s in range(1, q)
        if F["chi"][s] == 1 and coarse_Q_class(F, s) == "--N1"
    ]


def extra_is_vert_wrong(p: int) -> bool:
    return set(extra_mm(p)) == set(vert_units(p))


def dest_slopes_from_0(s: int, p: int = P) -> list[str]:
    F = _kernel_field(p)
    dest = {slope(F["mul"][s][a], p) for a in range(1, p)}
    return sorted(dest)


def sc(J, k: int):
    return lambda a, J=J, k=k: J((k * a) % P)


def Lp(JA, JB, ks) -> Fraction:
    return L_pair_ks(JA, JB, ks)


def Lcross_219(ks) -> Fraction:
    return (Lp(J_0112557, sc(J_F1, 2), ks) + Lp(J_G1, sc(J_F1, 3), ks)) / 2


def Lcross_671(ks) -> Fraction:
    return (Lp(J_0223347, sc(J_F1, 3), ks) + Lp(J_1122366sg, sc(J_F1, 2), ks)) / 2


def Lcross_RRRC(ks) -> Fraction:
    return (Lp(J_R_reverse, J_R_reverse, ks) + Lp(J_R_reverse, sc(J_C_11222, 3), ks)) / 2


def Lcross_665(ks) -> Fraction:
    return (
        Lp(J_1122447, J_G2, ks)
        + Lp(J_1122447, J_1122366ap, ks)
        + Lp(J_F3, J_G2, ks)
        + Lp(J_F3, J_1122366ap, ks)
    ) / 4


def Lcross_344(ks) -> Fraction:
    return (Lp(J_0022566, sc(J_G3, 2), ks) + Lp(J_F4, J_G3, ks)) / 2


def Lcross_695(ks) -> Fraction:
    return Lp(J_0034455, J_G4, ks)


def Lsame_219(ks) -> Fraction:
    return (Lp(J_G1, sc(J_0112557, 3), ks) + Lp(J_F1, sc(J_F1, 2), ks)) / 2


def Lsame_671(ks) -> Fraction:
    return (Lp(J_1122366sg, J_0223347, ks) + Lp(J_F1, sc(J_F1, 2), ks)) / 2


def Lsame_RRRC(ks) -> Fraction:
    return (
        Lp(J_R_reverse, sc(J_R_reverse, 2), ks)
        + Lp(J_R_reverse, sc(J_C_11222, 2), ks)
    ) / 2


def Lsame_665(ks) -> Fraction:
    return (
        Lp(J_F3, sc(J_1122447, 2), ks) + Lp(J_G2, sc(J_1122366ap, 2), ks)
    ) / 2


def Lsame_344(ks) -> Fraction:
    return (Lp(J_F4, J_0022566, ks) + Lp(J_G3, sc(J_G3, 2), ks)) / 2


def Lsame_695(ks) -> Fraction:
    return (Lp(J_0034455, sc(J_0034455, 2), ks) + Lp(J_G4, sc(J_G4, 2), ks)) / 2


def LM1_mix(ks=KS_PP) -> Fraction:
    return (
        4
        * (
            Lcross_219(ks)
            + Lcross_671(ks)
            + Lcross_RRRC(ks)
            + Lcross_665(ks)
        )
        + 2 * Lcross_344(ks)
        + Lcross_695(ks)
    ) / 19


def LM1_nstar_phi_wrong() -> Fraction:
    return LM1_mix(KS_NSTAR)


def Lspl_same_pp_wrong() -> Fraction:
    return (
        4
        * (
            Lsame_219(KS_PP)
            + Lsame_671(KS_PP)
            + Lsame_RRRC(KS_PP)
            + Lsame_665(KS_PP)
        )
        + 2 * Lsame_344(KS_PP)
        + Lsame_695(KS_PP)
    ) / 19


def prove_A() -> dict:
    e5 = extra_mm(5)
    e7 = extra_mm(7)
    sl = [slope(s) for s in e7]
    ok = e5 == []
    ok = ok and named_sizes(5)["--N1"] == 0
    ok = ok and set(e7) == set(EXTRAS7)
    ok = ok and len(e7) == named_sizes(7)["--N1"] == 4
    ok = ok and set(sl) == {"3", "4"}
    ok = ok and extra_pp(5) == []
    ok = ok and not extra_is_vert_wrong(7)
    if extra_is_vert_wrong(7) or e5:
        ok = False
    return {
        "proved": bool(ok),
        "p5": e5,
        "p7": e7,
        "slopes": sl,
        "theorem": (
            "extra --N1 empty at p=5; {9,12,44,47} slopes {3,4} at p=7. "
            "Fail: extras vertical; nonempty at p=5."
        ),
    }


def prove_B() -> dict:
    dests = {s: dest_slopes_from_0(s) for s in extra_mm(7)}
    ok = all(set(v) <= {"3", "4"} and "inf" not in v for v in dests.values())
    ok = ok and set().union(*[set(v) for v in dests.values()]) == {"3", "4"}
    if any("inf" in v for v in dests.values()):
        ok = False
    return {
        "proved": bool(ok),
        "dest0": {str(k): v for k, v in dests.items()},
        "theorem": (
            "--N1 extras send slope 0 to {3,4}, not {inf}. "
            "Fail: vertical dest."
        ),
    }


def prove_C() -> dict:
    rows = {
        "219": Lcross_219(KS_PP),
        "671": Lcross_671(KS_PP),
        "RRRC": Lcross_RRRC(KS_PP),
        "665": Lcross_665(KS_PP),
        "344": Lcross_344(KS_PP),
        "695": Lcross_695(KS_PP),
    }
    tgt = {
        "219": Fraction(325, 3),
        "671": Fraction(325, 3),
        "RRRC": 105,
        "665": Fraction(1309, 12),
        "344": Fraction(213, 2),
        "695": Fraction(323, 3),
    }
    mix = LM1_mix(KS_PP)
    nstar = LM1_nstar_phi_wrong()
    lspl = Lspl_same_pp_wrong()
    ok = all(rows[k] == tgt[k] for k in tgt)
    ok = ok and mix == LIVE_LM1[7] == Fraction(6131, 57)
    ok = ok and nstar == LIVE_LN_CROSS == Fraction(6269, 57)
    ok = ok and nstar != mix
    ok = ok and lspl == LIVE_LSPL_NL[7] == Fraction(673, 6)
    ok = ok and lspl != mix
    if nstar == mix or lspl == mix:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {k: str(v) for k, v in rows.items()},
        "mix": str(mix),
        "Nstar_phi": str(nstar),
        "same_pp": str(lspl),
        "theorem": (
            "CROSS slots under ++ φ mix to 6131/57. "
            "Fail: N* φ 6269/57; same-pair ++ 673/6."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_--N1 named at p=5,7 as CROSS mix under ++ φ. "
            "Q_NL unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.427  --N1 CROSS mix under ++ φ named", flush=True)
    A = prove_A()
    print(f"  A extras: {A['proved']} p7={A['p7']}", flush=True)
    B = prove_B()
    print(f"  B dest0: {B['proved']} {B['dest0']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} {C['mix']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"C": {"mix": C["mix"], "Nstar_phi": C["Nstar_phi"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.427",
        "title": "L_--N1 is CROSS mix under ++ φ; Q_NL still open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "LM1_extras": A["proved"],
            "LM1_dest_34": B["proved"],
            "LM1_mix_6131_57": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15427.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
