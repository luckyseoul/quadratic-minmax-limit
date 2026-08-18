#!/usr/bin/env python3
"""
Prop 15.545 — NUM_SUM splits as n_1d Q_1d^{++} + M_NL, and
    M_NL = 16 p A − n_1d Q_1d^{++}
equals the live NL mass at p=5,7.  Fail drop the 1D term;
fail Q_1d^{sub} in place of Q_1d^{++}; fail drop 3k in A.
Not proved for p≥11.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.544 flags.  Does **not** overwrite 15.522–15.544.

============================================================================
Setup.  15.366: NUM_SUM = Σ_{y∈H_+} Q_y^{++} = N·Q++.
15.356: every Type+ 1D has the same ++-average Q_1d^{++}, and
Q_1d^{++}=(p−3) Q_1d^{sub}/n_pp (zeros on ++ \\ F_p).
NL = H_+ \\ {Type+ 1D}.  15.367: NUM_SUM=16p A at p=5,7 only,
A=2 n_1d−4 n_pp−3k.

============================================================================
Theorem A — PROVED (partition + 15.356; Max+-free algebra).
  NUM_SUM = n_1d Q_1d^{++} + M_NL,   M_NL := Σ_{y∈NL} Q_y^{++}.
  Fail: write the 1D sum as n_1d Q_1d^{sub} (p=5: 160≠160/3).  ∎

Theorem B — PROVED (A + 15.367 live; p=5,7 caches).
  M_NL = 16p A − n_1d Q_1d^{++} equals (N−n_1d) Q_NL
  (1280/3, 61936/3).  Fail: drop the 1D term (480≠1280/3 at p=5).
  Fail: drop 3k in A (p=5: 3680/3≠1280/3).  ∎

Theorem C — OPEN.  16p A is not a theorem for p≥11, so the
  named remainder is only a candidate there.  Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: Fraction mix.  Writes
evidence/e1_gmin_m4_prop15545.json
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d, Q_1d_sub_over_q2
from e1_gmin_m4_prop15330 import HPLUS, LIVE_QPP
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15356 import Q_1d_pp_named, Q_NL_from_live
from e1_gmin_m4_prop15367 import A_num, NUM_SUM_named, k_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15545.json"


def one_d_mass(p: int) -> Fraction:
    """Σ_{1D} Q_y^{++} = n_1d Q_1d^{++}."""
    return n_1d(p) * Q_1d_pp_named(p)


def one_d_mass_wrong_sub(p: int) -> Fraction:
    """Fail-when-wrong: forget (p−3)/n_pp zeros on ++ \\ F_p."""
    return n_1d(p) * Q_1d_sub_over_q2(p)


def M_NL_named(p: int) -> Fraction:
    """16p A − n_1d Q_1d^{++}."""
    return NUM_SUM_named(p) - one_d_mass(p)


def M_NL_drop_1d(p: int) -> Fraction:
    """Fail-when-wrong: claim the NL mass is the whole 16p A."""
    return Fraction(NUM_SUM_named(p))


def M_NL_drop_k(p: int) -> Fraction:
    """Fail-when-wrong: drop 3k from A before subtracting 1D."""
    A0 = 2 * n_1d(p) - 4 * n_pp_named(p)
    return 16 * p * A0 - one_d_mass(p)


def M_NL_live(p: int) -> Fraction:
    n_nl = HPLUS[p] - n_1d(p)
    return n_nl * Q_NL_from_live(p)


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13):
        right = one_d_mass(p)
        wrong = one_d_mass_wrong_sub(p)
        row_ok = wrong != right
        if p == 5:
            row_ok = row_ok and right == Fraction(160, 3) and wrong == 160
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_1d_Q1d_pp": str(right),
            "n_1d_Q1d_sub": str(wrong),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "NUM_SUM = n_1d Q_1d^{++} + M_NL. Fail n_1d Q_1d^{sub} "
            f"(p=5: {rows['5']['n_1d_Q1d_sub']}≠{rows['5']['n_1d_Q1d_pp']})."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        named = M_NL_named(p)
        live = M_NL_live(p)
        drop1 = M_NL_drop_1d(p)
        dropk = M_NL_drop_k(p)
        row_ok = named == live and drop1 != live and dropk != live
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "named": str(named),
            "live": str(live),
            "drop_1d": str(drop1),
            "drop_k": str(dropk),
            "n_NL": HPLUS[p] - n_1d(p),
            "Q_NL": str(Q_NL_from_live(p)),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "M_NL=16pA−n_1d Q_1d^{++} hits live 1280/3, 61936/3. "
            f"Fail drop 1D (p=5: {rows['5']['drop_1d']}≠{rows['5']['live']}). "
            f"Fail drop 3k (p=5: {rows['5']['drop_k']}≠{rows['5']['live']})."
        ),
    }


def prove_open() -> dict:
    pred = {str(p): str(M_NL_named(p)) for p in (11, 13, 17)}
    return {
        "proved": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "M_NL_pred_no_cache": pred,
        "note": (
            "Named remainder assumes 15.367 NUM_SUM=16pA, which is "
            "only live at p=5,7. Predictions at p≥11 are untested. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.545 NUM_SUM = n_1d Q_1d++ + M_NL; M_NL=16pA-n_1d Q_1d++", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.545",
        "title": "NL mass M_NL=16pA−n_1d Q_1d^{++} at p=5,7",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "k_named": {str(p): k_named(p) for p in (5, 7)},
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
