#!/usr/bin/env python3
"""
Prop 15.433 — Sign of NL L_ns deviation is field-named:
    p=5: δ=χ_p(N+1)/2
    p=7 mixed, N≠−1: δ=χ_p(N+1)·3/38
    p=7 mixed, N=−1: δ=0
    p=7 unmixed: δ=χ(s−1)·30/19
The amplitudes 3/38 and 30/19 are not named in p (19 is the
Hoffman mix, not a field GJ).  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.432 flags.

============================================================================
Setup.  δ=L₊−μ₊μ₋.  15.432: L₊+L₋=2μμ, Paley×N table.
χ_p is the F_p Legendre symbol.  Unmixed means Paley=(±1,±1)
same sign; mixed means opposite signs.

============================================================================
Theorem A — PROVED (LIVE table + χ_p; certified p=5).
  δ=χ_p(N+1)/2 on every p=5 class.  Fail: L≡μμ;
  fail: replace the field 2 by the literal 19
  (χ_p(N+1)/19 ≠ ±1/2).  ∎

Theorem B — PROVED (LIVE table; certified p=7).
  Mixed N=−1: δ=0.  Mixed N≠−1: sign(δ)=χ_p(N+1) and
  |δ|=3/38.  Unmixed: sign(δ)=χ(s−1) and |δ|=30/19.
  Fail: the p=5 amplitude 1/2 at p=7 mixed (1/2≠3/38);
  fail: N-only (N=3 has two |δ|).  ∎

Theorem C — OPEN.  30/19 and 3/38 are Hoffman-mix numbers,
  not field GJ in p.  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction + 15.432 LIVE_LNS.  Writes
evidence/e1_gmin_m4_prop15433.json
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
from e1_gmin_m4_prop15298 import mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15314 import chi_p_pm  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15432 import LIVE_LNS  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15433_catalog.json"

AMP_P5 = Fraction(1, 2)
AMP_MIX7 = Fraction(3, 38)
AMP_UNM7 = Fraction(30, 19)


def mumu(p: int) -> Fraction:
    return mu_plus(p) * mu_minus(p)


def delta_plus(p: int, key: tuple[int, int, int]) -> Fraction:
    lp, _lm = LIVE_LNS[p][key]
    return lp - mumu(p)


def p5_law(N: int) -> Fraction:
    return Fraction(chi_p_pm(5, N + 1), 2)


def p5_law_19_wrong(N: int) -> Fraction:
    """Fail-when-wrong: insert 19 as the p=5 denominator."""
    return Fraction(chi_p_pm(5, N + 1), 19)


def prove_A() -> dict:
    ok = True
    rows = {}
    for key, (lp, lm) in LIVE_LNS[5].items():
        dlt = lp - mumu(5)
        law = p5_law(key[2])
        rows[str(key)] = {"δ": str(dlt), "law": str(law)}
        if dlt != law:
            ok = False
        if p5_law_19_wrong(key[2]) == dlt:
            ok = False
    if p5_law_19_wrong(2) == AMP_P5 or p5_law_19_wrong(2) == -AMP_P5:
        ok = False
    if p5_law_19_wrong(2) != Fraction(-1, 19):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "hand19": str(p5_law_19_wrong(2)),
        "theorem": (
            "p=5 δ=χ_p(N+1)/2. Fail: L≡μμ; fail: den 19 "
            "(χ_p(N+1)/19≠±1/2)."
        ),
    }


def prove_B() -> dict:
    ok = True
    mix_amp = set()
    unm_amp = set()
    n3_amp = set()
    for (a, b, N), (lp, lm) in LIVE_LNS[7].items():
        dlt = lp - mumu(7)
        mixed = a != b
        if N % 7 == 6:
            if dlt != 0:
                ok = False
            continue
        sgn = chi_p_pm(7, N + 1) if mixed else a
        if mixed:
            mix_amp.add(abs(dlt))
            if dlt != sgn * AMP_MIX7:
                ok = False
        else:
            unm_amp.add(abs(dlt))
            if dlt != sgn * AMP_UNM7:
                ok = False
        if N == 3:
            n3_amp.add(abs(dlt))
    ok = ok and mix_amp == {AMP_MIX7}
    ok = ok and unm_amp == {AMP_UNM7}
    ok = ok and n3_amp == {AMP_MIX7, AMP_UNM7}
    if AMP_P5 == AMP_MIX7:
        ok = False
    return {
        "proved": bool(ok),
        "mix_amp": str(AMP_MIX7),
        "unm_amp": str(AMP_UNM7),
        "p5_amp": str(AMP_P5),
        "N3_amps": sorted(str(x) for x in n3_amp),
        "theorem": (
            "p=7 mixed δ=χ_p(N+1)·3/38 (0 if N=−1); "
            "unmixed δ=χ(s−1)·30/19. Fail: amp 1/2 at p=7; "
            "N-only."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "amp_unm_has_19": AMP_UNM7.denominator % 19 == 0,
        "amp_mix_has_19": AMP_MIX7.denominator % 19 == 0,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "Signs named. Amplitudes 30/19 and 3/38 not field GJ. "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.433  L_ns δ signs field-named", flush=True)
    A = prove_A()
    print(f"  A p5 χ_p(N+1)/2: {A['proved']} hand19={A['hand19']}", flush=True)
    B = prove_B()
    print(f"  B p7 signs: {B['proved']} mix={B['mix_amp']} unm={B['unm_amp']}", flush=True)
    C = prove_open()
    print(f"  C open: 19|amp={C['amp_unm_has_19']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "hand19": A["hand19"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.433",
        "title": "L_ns δ signs are χ_p(N+1) / χ(s−1); amps not in p",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "p5_chi_N1": A["proved"],
            "p7_sign_law": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15433.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
