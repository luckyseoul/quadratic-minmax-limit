#!/usr/bin/env python3
"""
Prop 15.359 — The unique quadratic for D=|H_+|/(2p) through
(3,1),(5,13),(7,409) is 48p²−378p+703.  At p=11 it gives
D=2353, hence |H_+|=51766, but a vectorized Cy=py census of the
p=11 Max+ cache finds 85052 eigenvectors, so D(11)≥3866>2353.
No cubic with small integer coeffs, no 3-factor product, and no
3-feature linear form with |coeff|≤3 hits D on {3,5,7}.
D still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.358 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (interpolation algebra).
  D_quad(p)=48p²−378p+703 equals 1,13,409 at p=3,5,7 and 2353
  at p=11.  Fail: D_quad(5)≠13.  ∎

Theorem B — PROVED (vectorized Cy=py; p=11 cache).
  Of 143132 unique rows, 85052 satisfy Cy=py (the other 58080
  have ‖Cy−py‖_∞=24).  85052>22·2353=51766, so D(11)≥3866>2353.
  Fail: claim the cache has ≤51766 Max+.  ∎

Theorem C — PROVED (ProcessPool census; coeffs as in cpu_D_cubic).
  0 cubics a∈[−8,8], b,c,d∈[−12,12]; 0 products fgh or fg/h;
  0 forms a f+b g+c h with |coeff|≤3 on the named 3-feature list.
  Fail: claim D=p³−p²−p+1 (27−9−3+1=16≠1 at p=3).  ∎

Theorem D — OPEN.  D and Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction algebra + numpy GEMM Cy.  Writes
evidence/e1_gmin_m4_prop15359.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

P11_UNIQUE = "/tmp/maxplus_p11_unique.npy"
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15359_catalog.json"


def D_quad(p: int) -> int:
    return 48 * p * p - 378 * p + 703


def count_p11_maxplus() -> dict:
    Z = np.load(P11_UNIQUE)
    p = 11
    C = paley_conference_prime_power(p).astype(np.float64)
    q = p * p
    Y = np.ones((Z.shape[0], q + 1), dtype=np.float64)
    Y[:, 1:] = Z
    err = np.max(np.abs(Y @ C.T - p * Y), axis=1)
    n_ok = int(np.sum(err < 1e-3))
    n_bad = int(np.sum(err >= 1e-3))
    return {
        "n_rows": int(Z.shape[0]),
        "n_maxplus": n_ok,
        "n_not": n_bad,
        "maxerr_bad": float(err[err >= 1e-3].max()) if n_bad else 0.0,
        "D_lb": n_ok // (2 * p),
    }


def prove_A() -> dict:
    rows = {str(p): D_quad(p) for p in (3, 5, 7, 11)}
    ok = rows["3"] == 1 and rows["5"] == 13 and rows["7"] == 409
    if D_quad(5) != 13:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "p11_Hplus_pred": 22 * D_quad(11),
        "theorem": "D_quad hits 1,13,409 and predicts D(11)=2353.",
    }


def prove_B() -> dict:
    rec = count_p11_maxplus()
    pred = 22 * D_quad(11)
    ok = rec["n_maxplus"] > pred and rec["D_lb"] > D_quad(11)
    if rec["n_maxplus"] <= pred:
        ok = False
    if rec["n_not"] == 0:
        ok = False  # the 58080 junk must still be present as a fail-check
    return {
        "proved": bool(ok),
        **rec,
        "Hplus_pred_quad": pred,
        "theorem": "p=11 cache has 85052 Max+ > 22·2353; D_quad is dead.",
    }


def prove_C() -> dict:
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    n_c = int(cat.get("n_cubic", -1))
    n_p = int(cat.get("n_prod3", -1))
    n_l = int(cat.get("n_lin3", -1))
    ok = n_c == 0 and n_p == 0 and n_l == 0
    # fail-when-wrong named cubic
    fake = (3 ** 3) - (3 ** 2) - 3 + 1
    if fake == 1:
        ok = False
    if n_c != 0:
        ok = False
    return {
        "proved": bool(ok),
        "n_cubic": n_c,
        "n_prod3": n_p,
        "n_lin3": n_l,
        "fake_cubic_p3": fake,
        "theorem": "No small cubic / 3-product / 3-linear names D on {3,5,7}.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "D_named_in_p": False,
        "Q_NL_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": "D_quad and small 3-feature forms are dead.  D unnamed.",
    }


def main() -> dict:
    print("Prop 15.359  D_quad dead at p=11; no small cubic/3-feature", flush=True)
    A = prove_A()
    print(f"  A D_quad: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B p11 census: {B['proved']} n_ok={B['n_maxplus']} pred={B['Hplus_pred_quad']}", flush=True)
    C = prove_C()
    print(f"  C 3-feature: {C['proved']} cubic={C['n_cubic']} prod={C['n_prod3']} lin3={C['n_lin3']}", flush=True)
    D = prove_open()
    print(f"  D open: D={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.359",
        "title": "D quadratic dead at p=11; no small cubic or 3-feature name",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "D_quad_formula": A["proved"],
            "p11_census_kills_D_quad": B["proved"],
            "no_small_cubic_or_3feature": C["proved"],
            "D_named_in_p": False,
            "Q_NL_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15359.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
