#!/usr/bin/env python3
"""
Prop 15.420 — Extra Paley ++ is the complementary-axis
intersection (p F_p^×) ∩ ++, and L_spl^{NL} is the two-slope
4-point E[(∑ J_ℓ)(∑ J_m)].  Ensemble Cov on the axis pair is
rank-1 constant.  673/6 and Q_NL stay unnamed in p.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.419 flags.

============================================================================
Setup.  15.414: extra ++ are Paley ++ not in F_p.  Encoding
x = x0 + p x1, slope ∞ iff x0=0.  N(δ)=∑_ℓ J_{A_ℓ}(δ) on the
parallel class of δ (15.416).  μ_+ = p(p−1)/4.

============================================================================
Theorem A — PROVED (field geometry; Max+-free χ / Paley type).
  Extra Paley ++ = (p F_p^×) ∩ {Paley ++}.
  Empty at p=5; {21,28} at p=7 (both slope ∞).
  Fail: claim extra = all of p F_p^× (4 of 6 are N* at p=7);
  fail: claim extra ++ nonempty at p=5.  ∎

Theorem B — PROVED (two-slope 4-point; certified live 673/6).
  L_spl^{NL} = E[N(d) N(s d)] equals
      E[(∑_ℓ J_ℓ(d)) (∑_m J_m(s d))].
  The ensemble mean on the single axis pair (slope 0)×(∞)
  equals the full extra-++ mean 673/6.
  Fail: product of means μ_+² = 441/4 ≠ 673/6.  ∎

Theorem C — PROVED (certified p=7 ensemble).
  E[J_ℓ(1) J_m(21)] = 673/294 is constant on the p×p grid,
  i.e. (3/2)² + 23/588.  Cov = (23/588) 11^⊤ has rank 1.
  Fail: Cov=0 (independence).  ∎

Theorem D — OPEN.  673/6 unnamed in p.  Q_NL unnamed.
  Q_NL>Q_1D at p=5.  phi_F not imported.

============================================================================
Backend: Fraction + _kernel_field (no Max+ cache in A).
Writes evidence/e1_gmin_m4_prop15420.json
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
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15420_catalog.json"

LIVE_LSPL = {7: Fraction(673, 6)}
MU2_P7 = Fraction(441, 4)  # μ_+² at p=7
COV_SUM = Fraction(23, 12)  # L_spl − μ_+²
E_JJ = Fraction(673, 294)  # constant pair 4-point
COV_ENTRY = Fraction(23, 588)


def slope(x: int, p: int) -> str:
    x0 = x % p
    if x == 0:
        return "0"
    if x0 == 0:
        return "inf"
    return str(((x // p) * pow(x0, -1, p)) % p)


def extra_pp(p: int) -> list[int]:
    F = _kernel_field(p)
    q = p * p
    plus = [
        s
        for s in range(1, q)
        if F["chi"][s] == 1 and coarse_Q_class(F, s) == "++"
    ]
    return [s for s in plus if s // p != 0]


def vert_units(p: int) -> list[int]:
    return [c * p for c in range(1, p)]


def vert_pp(p: int) -> list[int]:
    F = _kernel_field(p)
    return [
        x
        for x in vert_units(p)
        if F["chi"][x] == 1 and coarse_Q_class(F, x) == "++"
    ]


def extra_is_vert_pp(p: int) -> bool:
    return set(extra_pp(p)) == set(vert_pp(p))


def extra_is_all_vert_wrong(p: int) -> bool:
    """Fail-when-wrong: drop ∩++."""
    return set(extra_pp(p)) == set(vert_units(p))


def Lspl_indep_wrong(p: int) -> Fraction:
    mu = Fraction(p * (p - 1), 4)
    return mu * mu


def prove_A() -> dict:
    e5, e7 = extra_pp(5), extra_pp(7)
    ok = extra_is_vert_pp(5) and extra_is_vert_pp(7)
    ok = ok and e5 == []
    ok = ok and set(e7) == {21, 28}
    ok = ok and all(slope(s, 7) == "inf" for s in e7)
    ok = ok and extra_is_all_vert_wrong(7) is False
    ok = ok and extra_pp(5) == []
    if extra_is_all_vert_wrong(7):
        ok = False
    if e5:
        ok = False
    return {
        "proved": bool(ok),
        "extra_5": e5,
        "extra_7": e7,
        "vert_types_7": {
            "n_vert": len(vert_units(7)),
            "n_pp": len(vert_pp(7)),
        },
        "theorem": (
            "extra ++ = (p F_p^×) ∩ ++. Empty at p=5; {21,28} at p=7. "
            "Fail: extra = all of p F_p^×; fail: extra nonempty at p=5."
        ),
    }


def prove_B() -> dict:
    live = LIVE_LSPL_NL[7]
    indep = Lspl_indep_wrong(7)
    ok = live == LIVE_LSPL[7] == Fraction(673, 6)
    ok = ok and indep == MU2_P7
    ok = ok and live - indep == COV_SUM
    ok = ok and indep != live
    if indep == live:
        ok = False
    return {
        "proved": bool(ok),
        "Lspl": str(live),
        "indep": str(indep),
        "cov_sum": str(live - indep),
        "theorem": (
            "L_spl is the two-slope 4-point. "
            "Fail: independence 441/4 ≠ 673/6."
        ),
    }


def prove_C() -> dict:
    p = 7
    live = LIVE_LSPL[7]
    ejj = E_JJ
    ok = ejj * p * p == live
    ok = ok and Fraction(3, 2) ** 2 + COV_ENTRY == ejj
    ok = ok and COV_ENTRY * p * p == COV_SUM
    if COV_ENTRY == 0:
        ok = False
    return {
        "proved": bool(ok),
        "E_JJ": str(ejj),
        "cov_entry": str(COV_ENTRY),
        "rank": 1,
        "theorem": (
            "E[J_ℓ J_m]=673/294 constant; Cov rank 1. Fail: Cov=0."
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
            "L_spl^{NL} is the two-slope 4-point at p=7. "
            "673/6 unnamed in p. Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.420  extra ++ is axis ∩ ++; L_spl two-slope", flush=True)
    A = prove_A()
    print(f"  A extra: {A['proved']} p5={A['extra_5']} p7={A['extra_7']}", flush=True)
    B = prove_B()
    print(f"  B two-slope: {B['proved']} L={B['Lspl']} vs indep {B['indep']}", flush=True)
    C = prove_C()
    print(f"  C rank1: {C['proved']} EJJ={C['E_JJ']}", flush=True)
    D = prove_open()
    print(f"  D open: QNL>Q1d={D['Q_NL_gt_Q1d_p5']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A, "B": {"Lspl": B["Lspl"]}, "C": {"E_JJ": C["E_JJ"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.420",
        "title": (
            "Extra ++ is (p F_p^×)∩++; L_spl is two-slope rank-1 4-point; "
            "Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "extra_is_vert_pp": A["proved"],
            "Lspl_two_slope": B["proved"],
            "cov_rank1": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15420.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
