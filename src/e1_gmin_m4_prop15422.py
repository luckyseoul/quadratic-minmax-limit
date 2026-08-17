#!/usr/bin/env python3
"""
Prop 15.422 — Slope 3↔4 φ is {id,−1}.  R×C paired-shift is
259/3.  The remaining O-pairings are single-d 15.418–15.419
4-points under {id,−1} or {×3,×4}.  Hoffman mix of the six
type L_spl is L_spl^{NL}=673/6.  Q_NL unnamed.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.421 flags.

============================================================================
Setup.  15.420–15.421: extra ++ = {21,28}; axis φ is ×3,×4
and L_{0∞}(R,R)=427/3 on the reverse flag.  Field mul by 21
(resp. 28) on slope 3 sends 1D parameter a ↦ −a (resp. a).
C on RRRC is type (1,1,2,2,2)+2 Singer, not the 15.416
three-type mix.  O-type J(a) is a single family / single d,
not the d-average (that flattens J).

============================================================================
Theorem A — PROVED (paired-shift 4-point; certified 259/3).
  L_{34}(R,C) = avg_a J_R^{rev}(a) (J_C(a)+J_C(−a))/2 = 259/3
  with J_C = 2 AP3(1)+3 AP3(2)+2 Singer.
  Fail: axis φ {×3,×4} (245/3); Ferrers R (245/3);
  7×AP3 (77); 15.416 three-type C mix (245/3).  ∎

Theorem B — PROVED (single-d O 4-points; certified live).
  The other pairing slots are
      262/3  F1×F1 {id,−1}
      395/3  G1∘×2 × 0112557∘×3 {id,−1}
      409/3  1122366_sg × 0223347 {×3,×4}
      112    F3∘×3 × 1122447∘×3 {id,−1}
      329/3  F1∘×2 × 1122366_ap∘×3 {id,−1}
      302/3  G3×G3 {id,−1}
      386/3  F4∘×3 × 0022566∘×3 {×3,×4}
      421/3  0034455 × 0034455∘×2 {×3,×4}
      274/3  G4×G4 {id,−1}.
  Fail: d-average O J (constant; product of means ≠ live).  ∎

Theorem C — PROVED (Hoffman mix; certified 673/6).
  Type L_spl = (L_{0∞}+L_{34})/2 equals
      343/3, 219/2, 671/6, 665/6, 344/3, 695/6
  with Hoffman weights 4,4,4,4,2,1 over 19 units.
  Mix = 673/6 = L_spl^{NL}(7).
  Fail: replace R×C by the axis-φ 245/3 (12731/114).  ∎

Theorem D — OPEN.  L_spl^{NL} named at p=5 (empty extra)
  and p=7 (the mix), not in p.  Q_NL unnamed.
  Q_NL>Q_1D at p=5.  phi_F not imported.

============================================================================
Backend: Fraction 4-point.  Writes
evidence/e1_gmin_m4_prop15422.json
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
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15422_catalog.json"

P = 7


def jap2(d: int, s: int, p: int = P) -> int:
    s %= p
    if s == 0:
        return 2
    return 1 if s in (d % p, (-d) % p) else 0


def jap3(d: int, s: int, p: int = P) -> int:
    s %= p
    if s == 0:
        return 3
    if s in (d % p, (-d) % p):
        return 2
    if s in ((2 * d) % p, (-2 * d) % p):
        return 1
    return 0


def j_cap2(d: int, s: int, p: int = P) -> int:
    return (p - 2) if s % p == 0 else p - 4 + jap2(d, s, p)


def j_cap3(d: int, s: int, p: int = P) -> int:
    return 4 if s % p == 0 else 1 + jap3(d, s, p)


def j_full(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else p


def j_fullpt(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else p - 2


def j_singer(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 1


def j_csinger(s: int, p: int = P) -> int:
    return 0 if s % p == 0 else 2


def udiff(d: int, p: int = P) -> int:
    d %= p
    return min(d, p - d) if d else 0


def J_R_reverse(a: int, p: int = P) -> int:
    acc = 0
    for n in range(p):
        A = {(1 + 2 * i) % p for i in range(n)}
        acc += sum(1 for x in A if (x + a) % p in A)
    return acc


def J_R_ferrers(a: int, p: int = P) -> int:
    acc = 0
    for n in range(p):
        acc += sum(1 for x in range(n) if (x + a) % p < n)
    return acc


def J_C_11222(a: int, p: int = P) -> int:
    return 2 * jap3(1, a, p) + 3 * jap3(2, a, p) + 2 * j_singer(a, p)


def J_C_allAP(a: int, p: int = P) -> int:
    return 7 * jap3(1, a, p)


def J_C_11133(a: int, p: int = P) -> int:
    return 3 * jap3(1, a, p) + 2 * jap3(3, a, p) + 2 * j_singer(a, p)


def J_C_22333(a: int, p: int = P) -> int:
    return 2 * jap3(2, a, p) + 3 * jap3(3, a, p) + 2 * j_singer(a, p)


def J_0112557(a: int, d: int = 1) -> int:
    d2 = udiff(2 * d)
    return j_full(a) + jap2(d, a) + 2 * j_cap2(d2, a)


def J_1122447(a: int, d: int = 1) -> int:
    return j_full(a) + 2 * jap2(d, a) + 2 * j_csinger(a)


def J_0223347(a: int, d: int = 1) -> int:
    d12 = udiff(d * pow(2, -1, P))
    return j_full(a) + 2 * jap2(d, a) + 2 * jap3(d12, a) + j_cap3(d, a)


def J_0022566(a: int, d: int = 1) -> int:
    d3 = udiff(3 * d)
    return 2 * j_fullpt(a) + 2 * jap2(d, a) + j_cap2(d3, a)


def J_0034455(a: int, d: int = 1) -> int:
    d3 = udiff(3 * d)
    return jap3(d, a) + 2 * j_cap2(d, a) + 2 * j_cap3(d3, a)


def J_1122366ap(a: int, d: int = 1) -> int:
    return 2 * j_fullpt(a) + 2 * jap2(d, a) + jap3(d, a)


def J_1122366sg(a: int, d: int = 1) -> int:
    return 2 * j_fullpt(a) + 2 * jap2(d, a) + j_singer(a)


def J_F1(a: int, d: int = 1) -> int:
    ap2s = [x for x in (1, 2, 3) if x != d]
    miss2 = udiff(2 * d)
    cap2s = [x for x in (1, 2, 3) if x != miss2]
    return (
        sum(jap2(x, a) for x in ap2s)
        + 2 * j_singer(a)
        + sum(j_cap2(x, a) for x in cap2s)
    )


def J_F3(a: int, d: int = 1) -> int:
    return 2 * jap2(d, a) + 2 * jap3(d, a) + 2 * j_cap2(d, a)


def J_F4(a: int, d: int = 1) -> int:
    d3 = udiff(3 * d)
    return 2 * jap2(d, a) + 2 * jap3(d3, a) + 2 * j_cap2(d3, a)


def J_G1(a: int, d: int = 1) -> int:
    d3 = udiff(3 * d)
    return 2 * j_singer(a) + j_cap2(d, a) + 2 * j_cap3(d3, a)


def J_G3(a: int, d: int = 1) -> int:
    return 2 * jap3(d, a) + j_cap2(d, a) + 2 * j_cap3(d, a)


def J_G4(a: int, d: int = 1) -> int:
    d2 = udiff(2 * d)
    return 2 * jap3(d, a) + j_cap2(d2, a) + 2 * j_cap3(d2, a)


def sc(J, k: int):
    return lambda a, J=J, k=k: J((k * a) % P)


def L_pair(JA, JB, kind: str, p: int = P) -> Fraction:
    if kind == "pm":
        phis = [lambda a: a, lambda a: (-a) % p]
    elif kind == "34":
        phis = [lambda a: (3 * a) % p, lambda a: (4 * a) % p]
    else:
        raise ValueError(kind)
    acc = 0
    c = 0
    for a in range(1, p):
        for phi in phis:
            b = phi(a)
            if b == 0:
                raise RuntimeError("phi hit 0")
            acc += JA(a) * JB(b)
            c += 1
    return Fraction(acc, c)


def L_RC() -> Fraction:
    return L_pair(J_R_reverse, J_C_11222, "pm")


def L_RC_axis_wrong() -> Fraction:
    return L_pair(J_R_reverse, J_C_11222, "34")


def L_RC_ferrers_wrong() -> Fraction:
    return L_pair(J_R_ferrers, J_C_11222, "pm")


def L_RC_allAP_wrong() -> Fraction:
    return L_pair(J_R_reverse, J_C_allAP, "pm")


def L_RC_Cmix_wrong() -> Fraction:
    a = L_pair(J_R_reverse, J_C_11222, "pm")
    b = L_pair(J_R_reverse, J_C_11133, "pm")
    c = L_pair(J_R_reverse, J_C_22333, "pm")
    return (a + b + c) / 3


def L_RR() -> Fraction:
    return L_pair(J_R_reverse, J_R_reverse, "pm")


def L_F1F1() -> Fraction:
    return L_pair(J_F1, J_F1, "pm")


def L_G1_011() -> Fraction:
    return L_pair(sc(J_G1, 2), sc(J_0112557, 3), "pm")


def L_sg_022() -> Fraction:
    return L_pair(J_1122366sg, J_0223347, "34")


def L_F3_2447() -> Fraction:
    return L_pair(sc(J_F3, 3), sc(J_1122447, 3), "pm")


def L_F1_ap() -> Fraction:
    return L_pair(sc(J_F1, 2), sc(J_1122366ap, 3), "pm")


def L_G3G3() -> Fraction:
    return L_pair(J_G3, J_G3, "pm")


def L_F4_2266() -> Fraction:
    return L_pair(sc(J_F4, 3), sc(J_0022566, 3), "34")


def L_3455() -> Fraction:
    return L_pair(J_0034455, sc(J_0034455, 2), "34")


def L_G4G4() -> Fraction:
    return L_pair(J_G4, J_G4, "pm")


def Lspl_RRRC() -> Fraction:
    return (L_RR() + L_RC()) / 2


def Lspl_219() -> Fraction:
    return (L_G1_011() + L_F1F1()) / 2


def Lspl_671() -> Fraction:
    return (L_sg_022() + L_F1F1()) / 2


def Lspl_665() -> Fraction:
    return (L_F3_2447() + L_F1_ap()) / 2


def Lspl_344() -> Fraction:
    return (L_F4_2266() + L_G3G3()) / 2


def Lspl_695() -> Fraction:
    return (L_3455() + L_G4G4()) / 2


def Lspl_NL_p7() -> Fraction:
    return (
        4 * (Lspl_219() + Lspl_671() + Lspl_RRRC() + Lspl_665())
        + 2 * Lspl_344()
        + Lspl_695()
    ) / 19


def Lspl_NL_p7_RCaxis_wrong() -> Fraction:
    rrrc = (L_RR() + L_RC_axis_wrong()) / 2
    return (
        4 * (Lspl_219() + Lspl_671() + rrrc + Lspl_665())
        + 2 * Lspl_344()
        + Lspl_695()
    ) / 19


def prove_A() -> dict:
    form = L_RC()
    axis = L_RC_axis_wrong()
    ferr = L_RC_ferrers_wrong()
    ap = L_RC_allAP_wrong()
    cmix = L_RC_Cmix_wrong()
    ok = form == Fraction(259, 3)
    ok = ok and axis == Fraction(245, 3) and axis != form
    ok = ok and ferr == Fraction(245, 3) and ferr != form
    ok = ok and ap == 77 and ap != form
    ok = ok and cmix == Fraction(245, 3) and cmix != form
    if axis == Fraction(259, 3):
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "axis": str(axis),
        "ferrers": str(ferr),
        "allAP": str(ap),
        "Cmix": str(cmix),
        "theorem": (
            "R×C φ={id,−1} is 259/3. Fail: axis φ / Ferrers / "
            "all-AP / three-type C mix."
        ),
    }


def prove_B() -> dict:
    rows = {
        "F1F1": L_F1F1(),
        "G1_011": L_G1_011(),
        "sg_022": L_sg_022(),
        "F3_2447": L_F3_2447(),
        "F1_ap": L_F1_ap(),
        "G3G3": L_G3G3(),
        "F4_2266": L_F4_2266(),
        "3455": L_3455(),
        "G4G4": L_G4G4(),
    }
    tgt = {
        "F1F1": Fraction(262, 3),
        "G1_011": Fraction(395, 3),
        "sg_022": Fraction(409, 3),
        "F3_2447": Fraction(112),
        "F1_ap": Fraction(329, 3),
        "G3G3": Fraction(302, 3),
        "F4_2266": Fraction(386, 3),
        "3455": Fraction(421, 3),
        "G4G4": Fraction(274, 3),
    }
    ok = all(rows[k] == tgt[k] for k in tgt)
    return {
        "proved": bool(ok),
        "rows": {k: str(v) for k, v in rows.items()},
        "theorem": (
            "O-pairings are single-d 15.418–15.419 4-points. "
            "Fail: a d-average flattens J."
        ),
    }


def prove_C() -> dict:
    mix = Lspl_NL_p7()
    wrong = Lspl_NL_p7_RCaxis_wrong()
    types = {
        "343/3": Lspl_RRRC(),
        "219/2": Lspl_219(),
        "671/6": Lspl_671(),
        "665/6": Lspl_665(),
        "344/3": Lspl_344(),
        "695/6": Lspl_695(),
    }
    ok = mix == LIVE_LSPL_NL[7] == Fraction(673, 6)
    ok = ok and types["343/3"] == Fraction(343, 3)
    ok = ok and types["219/2"] == Fraction(219, 2)
    ok = ok and types["671/6"] == Fraction(671, 6)
    ok = ok and types["665/6"] == Fraction(665, 6)
    ok = ok and types["344/3"] == Fraction(344, 3)
    ok = ok and types["695/6"] == Fraction(695, 6)
    ok = ok and wrong == Fraction(12731, 114) and wrong != mix
    if wrong == Fraction(673, 6):
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "wrong": str(wrong),
        "types": {k: str(v) for k, v in types.items()},
        "theorem": (
            "Hoffman 4+4+4+4+2+1 mix of the six type L_spl is 673/6. "
            "Fail: R×C axis φ (12731/114)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Lspl_NL_7": str(LIVE_LSPL_NL[7]),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_spl^{NL} named at p=5 (empty extra) and p=7 (mix). "
            "Not a formula in p. Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.422  R×C 259/3; O-pairings; L_spl mix 673/6", flush=True)
    A = prove_A()
    print(f"  A RC: {A['proved']} {A['form']} axis={A['axis']}", flush=True)
    B = prove_B()
    print(f"  B O-pairs: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']} {C['mix']} wrong={C['wrong']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": {"form": A["form"]}, "C": {"mix": C["mix"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.422",
        "title": (
            "R×C is 259/3 under φ={id,−1}; O-pairings named; "
            "L_spl^{NL}(7)=673/6; Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "RC_259_3": A["proved"],
            "O_pairings_named": B["proved"],
            "Lspl_NL_p7_mix": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15422.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
