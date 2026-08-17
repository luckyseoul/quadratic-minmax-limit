#!/usr/bin/env python3
"""
Prop 15.435 — The p=7 L_ns amplitudes 30/19 and 3/38 are the
Hoffman 4+4+4+4+2+1 mixes of six den|12 type-wise δ, not a
field / Jacobi denominator.  19 is the weight sum c(7), never
a field constant.  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.434 flags.

============================================================================
Setup.  15.307 NL types at p=7 have Hoffman weights
w = n / p²(p−1) ∈ {4,4,4,4,2,1}, sum w = 19 = c(7).
δ = L₊ − μ₊μ₋ on each (type, Paley, N) cell.  15.432/33:
NL-mean unmixed |δ|=30/19, mixed N≠−1 |δ|=3/38, N=−1 δ=0.
Census: L is constant on each (type, Paley, N); every cell
has L₊+L₋=2μμ.  Certified Gram fold (cpu_delta_types).

============================================================================
Theorem A — PROVED (LIVE type table; certified p=7).
  Unmixed |δ| on the six types is
    {7/6, 5/12, 1, 7/3, 53/12, 3/2},
  every denominator dividing 12, and
    (4·7/6 + 4·5/12 + 4·1 + 4·7/3 + 2·53/12 + 3/2) / 19
      = 30/19.
  Sign is χ(s−1).  Fail: mix denominator 18 (5/3 ≠ 30/19);
  fail: replace any atom den by 19 (p/6 → 7/19 ≠ 7/6).  ∎

Theorem B — PROVED (same table).
  Mixed N≠−1, χ_p(N+1)=+1, signed type-δ
    {−7/12, 1/6, 3/4, 7/12, −1/3, −3/2}
  mix to 3/38.  N=−1: every type δ=0.  Fail: mix den 18
  (1/12 ≠ 3/38); fail: unsigned mix of those |δ| (21/38 ≠ 3/38).  ∎

Theorem C — PROVED as a consistency check (15.434 LIVE_J + A/B).
  The 15.434 multipliers give J_Δ(pp_fp)=−4·(30/19)=−120/19
  and J_Δ(pp_U)=−8·(3/38)=−12/19.  The p=7 dens 19 are the
  same Hoffman weight sum.  Fail: J_fp=−4·(30/18).  ∎

Theorem D — OPEN.  The six atoms are p=7 den|12 names, not a
  general-p type dictionary (p=5 is a single type with |δ|=1/2
  ≠ 7/6).  Q_NL still carries this mix.  phi_F not imported.

============================================================================
Backend: Fraction tables (no Max+).  Writes
evidence/e1_gmin_m4_prop15435.json
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
from e1_gmin_m4_prop15432 import LIVE_LNS  # noqa: E402
from e1_gmin_m4_prop15433 import AMP_MIX7, AMP_P5, AMP_UNM7, mumu  # noqa: E402
from e1_gmin_m4_prop15434 import LIVE_J  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15435_catalog.json"

# 15.307 / 15.418 NL types at p=7.  w = n / p²(p−1).
HOFFMAN_W = (4, 4, 4, 4, 2, 1)
HOFFMAN_SUM = 19  # 4+4+4+4+2+1 = c(7)
FAIL_SUM = 18

TYPES = ("CRRR", "O219", "O665", "O671", "O344", "O695")

# Unmixed |δ| (Paley(--,--) N=3 and Paley(++,++) N=5).
UNMIXED_ABS = (
    Fraction(7, 6),
    Fraction(5, 12),
    Fraction(1),
    Fraction(7, 3),
    Fraction(53, 12),
    Fraction(3, 2),
)

# Mixed N≠−1, χ_p(N+1)=+1, signed δ.
MIXED_SIGNED = (
    Fraction(-7, 12),
    Fraction(1, 6),
    Fraction(3, 4),
    Fraction(7, 12),
    Fraction(-1, 3),
    Fraction(-3, 2),
)

# p=7 names of the unmixed atoms (certified identities, not all-p).
UNMIXED_NAME = (
    "p/6",
    "(p-2)/12",
    "1",
    "p/3",
    "(q+4)/12",
    "(p-1)/4",
)


def mix(weights, values, den: int) -> Fraction:
    return sum(w * v for w, v in zip(weights, values)) / den


def mix_unmixed(den: int = HOFFMAN_SUM) -> Fraction:
    return mix(HOFFMAN_W, UNMIXED_ABS, den)


def mix_mixed(den: int = HOFFMAN_SUM) -> Fraction:
    return mix(HOFFMAN_W, MIXED_SIGNED, den)


def mix_mixed_unsigned(den: int = HOFFMAN_SUM) -> Fraction:
    return mix(HOFFMAN_W, tuple(abs(v) for v in MIXED_SIGNED), den)


def atom_den19_wrong() -> dict[str, Fraction]:
    """Fail-when-wrong: replace each atom's live den by the literal 19."""
    p, q = 7, 49
    return {
        "p/6": Fraction(p, 19),
        "(p-2)/12": Fraction(p - 2, 19),
        "1": Fraction(1, 19),
        "p/3": Fraction(p, 19),
        "(q+4)/12": Fraction(q + 4, 19),
        "(p-1)/4": Fraction(p - 1, 19),
    }


def atom_live() -> dict[str, Fraction]:
    return dict(zip(UNMIXED_NAME, UNMIXED_ABS))


def prove_A() -> dict:
    ok = True
    got = mix_unmixed(HOFFMAN_SUM)
    if got != AMP_UNM7:
        ok = False
    if sum(HOFFMAN_W) != HOFFMAN_SUM:
        ok = False
    live_unm = LIVE_LNS[7][(-1, -1, 3)][0] - mumu(7)
    if abs(live_unm) != AMP_UNM7:
        ok = False
    if mix_unmixed(FAIL_SUM) == AMP_UNM7:
        ok = False
    if mix_unmixed(FAIL_SUM) != Fraction(5, 3):
        ok = False
    wrong = atom_den19_wrong()
    live = atom_live()
    for name, val in live.items():
        if val.denominator % 19 == 0:
            ok = False
        if 12 % val.denominator != 0:
            ok = False
        if wrong[name] == val:
            ok = False
    if wrong["p/6"] != Fraction(7, 19):
        ok = False
    if wrong["p/6"] == Fraction(7, 6):
        ok = False
    return {
        "proved": bool(ok),
        "mix19": str(got),
        "mix18": str(mix_unmixed(FAIL_SUM)),
        "atoms": {n: str(v) for n, v in live.items()},
        "hand19": {n: str(v) for n, v in wrong.items()},
        "theorem": (
            "Unmixed 30/19 = Hoffman mix of den|12 type-δ. "
            "Fail: den 18; fail: atom den 19."
        ),
    }


def prove_B() -> dict:
    ok = True
    got = mix_mixed(HOFFMAN_SUM)
    if got != AMP_MIX7:
        ok = False
    live_mx = LIVE_LNS[7][(-1, 1, 3)][0] - mumu(7)
    if live_mx != AMP_MIX7:
        ok = False
    if LIVE_LNS[7][(-1, 1, 6)][0] - mumu(7) != 0:
        ok = False
    if mix_mixed(FAIL_SUM) == AMP_MIX7:
        ok = False
    if mix_mixed(FAIL_SUM) != Fraction(1, 12):
        ok = False
    u = mix_mixed_unsigned(HOFFMAN_SUM)
    if u == AMP_MIX7:
        ok = False
    if u != Fraction(21, 38):
        ok = False
    return {
        "proved": bool(ok),
        "mix19": str(got),
        "mix18": str(mix_mixed(FAIL_SUM)),
        "unsigned": str(u),
        "theorem": (
            "Mixed 3/38 = signed Hoffman mix. "
            "Fail: den 18; fail: unsigned mix."
        ),
    }


def prove_C() -> dict:
    ok = True
    j_fp = -4 * AMP_UNM7
    j_U = -8 * AMP_MIX7
    if j_fp != LIVE_J[7]["pp_fp"]:
        ok = False
    if j_U != LIVE_J[7]["pp_U"]:
        ok = False
    if j_fp != Fraction(-120, 19):
        ok = False
    if j_U != Fraction(-12, 19):
        ok = False
    if -4 * mix_unmixed(FAIL_SUM) == LIVE_J[7]["pp_fp"]:
        ok = False
    return {
        "proved": bool(ok),
        "J_fp": str(j_fp),
        "J_U": str(j_U),
        "fail18": str(-4 * mix_unmixed(FAIL_SUM)),
        "theorem": (
            "Consistency: 15.434 J_Δ(pp_fp)=−4·(30/19), "
            "J_Δ(pp_U)=−8·(3/38). Fail: −4·(30/18)."
        ),
    }


def prove_open() -> dict:
    crrr_p5 = AMP_P5
    crrr_p7 = UNMIXED_ABS[0]
    return {
        "proved": False,
        "crrr_p5": str(crrr_p5),
        "crrr_p7": str(crrr_p7),
        "crrr_joint": crrr_p5 == crrr_p7,
        "amp_unm_has_19": AMP_UNM7.denominator % 19 == 0,
        "mix_kills_19_as_field": mix_unmixed() == AMP_UNM7,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "30/19 and 3/38 are Hoffman mixes, not field GJ. "
            "Six atoms remain p=7 den|12 names (CRRR 7/6 ≠ p=5 1/2). "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.435  30/19 and 3/38 are Hoffman mixes of den|12 δ", flush=True)
    A = prove_A()
    print(f"  A unmixed mix: {A['proved']} mix19={A['mix19']} mix18={A['mix18']}", flush=True)
    B = prove_B()
    print(f"  B mixed mix: {B['proved']} mix19={B['mix19']} unsigned={B['unsigned']}", flush=True)
    C = prove_C()
    print(f"  C J_Δ image: {C['proved']} Jfp={C['J_fp']} JU={C['J_U']}", flush=True)
    D = prove_open()
    print(
        f"  D open: CRRR joint={D['crrr_joint']} phi_F={D['phi_F_ge_6']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "mix19_unm": A["mix19"],
                "mix18_unm": A["mix18"],
                "mix19_mx": B["mix19"],
                "J_fp": C["J_fp"],
                "J_U": C["J_U"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.435",
        "title": "30/19 and 3/38 are Hoffman mixes of den|12 type-δ",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "unmixed_hoffman_mix": A["proved"],
            "mixed_hoffman_mix": B["proved"],
            "J_delta_is_mix_image": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15435.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
