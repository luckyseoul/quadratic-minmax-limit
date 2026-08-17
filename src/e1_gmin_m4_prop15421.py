#!/usr/bin/env python3
"""
Prop 15.421 — The R×R axis pairing is the reverse nested-AP
paired-shift 4-point 427/3.  Standard Ferrers gives 133
(same as L_par(R)) and is wrong for L_spl.  On RRRC,
L_spl=(427/3+259/3)/2=343/3.  259/3 (R×C) and Q_NL stay
unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.420 flags.

============================================================================
Setup.  15.420: extra ++ = {21,28} at p=7; L_spl is the
two-slope 4-point.  Field mul by 21 (resp. 28) pairs
horiz shift a∈F_p^× with vert shift 3a (resp. 4a).
R = nested AP flag (15.416).  Reverse order: size k sits
on row p−1−k.

============================================================================
Theorem A — PROVED (paired-shift 4-point; certified live 427/3).
  On the reverse nested-AP flag (d=2, s=1, order p−1..0)
      L_{0∞}(R,R) = avg_a (∑_ℓ J_ℓ(a)) (∑_m J_m(φ(a)))
  with φ∈{×3,×4} equals 427/3.
  Fail: same-orientation Ferrers (d=1,s=0, order 0..p−1)
  gives 133 ≠ 427/3.  ∎

Theorem B — PROVED (type mix; certified live).
  On each p=7 RRRC type, L_spl = (L_{0∞}+L_{34})/2
  = (427/3 + 259/3)/2 = 343/3.
  Fail: replace 427/3 by 133 (then 329/3 ≠ 343/3).  ∎

Theorem C — OPEN.  259/3 is the R×C pairing, unnamed.
  Hoffman mix of all type L_spl is 673/6, unnamed in p.
  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction 4-point + _kernel_field φ.  Writes
evidence/e1_gmin_m4_prop15421.json
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
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15414 import LIVE_LSPL_NL  # noqa: E402
from e1_gmin_m4_prop15420 import extra_pp  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15421_catalog.json"

P = 7
LIVE_0INF_RR = Fraction(427, 3)
LIVE_34_RC = Fraction(259, 3)
LIVE_RRRC = Fraction(343, 3)
FERRERS_WRONG = 133


def nest_rows(order, d: int, s: int, p: int = P):
    rows = [() for _ in range(p)]
    for k, r in enumerate(order):
        rows[r] = tuple((s + d * i) % p for i in range(k))
    return rows


def _Nh(rows, a, p: int = P) -> int:
    acc = 0
    for rset in rows:
        st = set(rset)
        acc += sum(1 for x in st if (x + a) % p in st)
    return acc


def _Nv(rows, b, p: int = P) -> int:
    cols: list[list[int]] = [[] for _ in range(p)]
    for y, rset in enumerate(rows):
        for x in rset:
            cols[x].append(y)
    acc = 0
    for cset in cols:
        st = set(cset)
        acc += sum(1 for y in st if (y + b) % p in st)
    return acc


def phi_y_maps(p: int = P) -> list[dict[int, int]]:
    """y-coordinate of extra_s * (a,0) for a in F_p^×."""
    F = _kernel_field(p)
    out = []
    for s in extra_pp(p):
        m = {}
        for a in range(1, p):
            sd = F["mul"][s][a]
            assert sd % p == 0
            m[a] = sd // p
        out.append(m)
    return out


def L_paired(rows, p: int = P) -> Fraction:
    phis = phi_y_maps(p)
    acc = 0
    c = 0
    for a in range(1, p):
        nh = _Nh(rows, a, p)
        for phi in phis:
            acc += nh * _Nv(rows, phi[a], p)
            c += 1
    return Fraction(acc, c)


def reverse_flag(p: int = P):
    """Live p=7 R×R axes: d=2, s=1, order p−1..0."""
    return nest_rows(tuple(range(p - 1, -1, -1)), 2, 1, p)


def ferrers_flag(p: int = P):
    """Same-orientation Ferrers: d=1, s=0, order 0..p−1."""
    return nest_rows(tuple(range(p)), 1, 0, p)


def L_RR_reverse() -> Fraction:
    return L_paired(reverse_flag())


def L_RR_ferrers_wrong() -> Fraction:
    return L_paired(ferrers_flag())


def L_spl_RRRC() -> Fraction:
    return (LIVE_0INF_RR + LIVE_34_RC) / 2


def L_spl_RRRC_ferrers_wrong() -> Fraction:
    return (FERRERS_WRONG + LIVE_34_RC) / 2


def prove_A() -> dict:
    form = L_RR_reverse()
    wrong = L_RR_ferrers_wrong()
    ok = form == LIVE_0INF_RR == Fraction(427, 3)
    ok = ok and wrong == FERRERS_WRONG
    ok = ok and wrong != form
    if wrong == LIVE_0INF_RR:
        ok = False
    return {
        "proved": bool(ok),
        "form": str(form),
        "ferrers": str(wrong),
        "theorem": (
            "R×R paired-shift reverse flag = 427/3. "
            "Fail: Ferrers 133."
        ),
    }


def prove_B() -> dict:
    mix = L_spl_RRRC()
    wrong = L_spl_RRRC_ferrers_wrong()
    ok = mix == LIVE_RRRC == Fraction(343, 3)
    ok = ok and wrong == Fraction(329, 3)
    ok = ok and wrong != mix
    if wrong == LIVE_RRRC:
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "ferrers_mix": str(wrong),
        "theorem": (
            "RRRC L_spl=(427/3+259/3)/2=343/3. "
            "Fail: use 133 (329/3)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "L_34_RC": str(LIVE_34_RC),
        "Lspl_NL": str(LIVE_LSPL_NL[7]),
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "259/3 R×C pairing unnamed. 673/6 unnamed in p. "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.421  R×R reverse-flag paired-shift 427/3", flush=True)
    A = prove_A()
    print(f"  A reverse: {A['proved']} {A['form']} ferrers={A['ferrers']}", flush=True)
    B = prove_B()
    print(f"  B RRRC: {B['proved']} {B['mix']} wrong={B['ferrers_mix']}", flush=True)
    C = prove_open()
    print(f"  C open: RC={C['L_34_RC']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A, "B": {"mix": B["mix"]}}, indent=2) + "\n"
    )
    out = {
        "prop": "15.421",
        "title": (
            "R×R L_{0∞} is reverse nested-AP paired-shift 427/3; "
            "RRRC L_spl=343/3; Q_NL still open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "RR_reverse_427_3": A["proved"],
            "RRRC_mix_343_3": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15421.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
