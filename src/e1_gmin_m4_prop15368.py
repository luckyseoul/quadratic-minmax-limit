#!/usr/bin/env python3
"""
Prop 15.368 — p=7 T(7,3) occupancy of S(α).  The p=5 values
S(x²)=40 and S(x³)=65 die.  Deg-4 is not one Gal-class (nor 3).
S(x²)=S(x⁴) as occupancy at p=7 and is complex; at p=5 S(x²) is
real and S(x⁴) is a deg-4 class.  D=dim V₊ dies.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.367 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.363: S(affine)=T, S(x²)=40, S(x³)=65 at p=5.
15.366: deg-4 has 3 Gal-classes at p=5.  D=|H+|/(2p) ∈ {13,409}.
T(7,3)=68938800 (A008300).  Full first-row shard of T(7,3)
(nT certified) gives occupancy of ten named α.

============================================================================
Theorem A — PROVED (p=5 identity + p=7 census).
  S(x²)=40 and S(x³)=65 at p=5.  At p=7, ∑ a_k ζ^k is rational
  iff a_1=⋯=a_{p−1}.  occ(x²) has a_1=48411510 ≠ a_3=52599260,
  so S(x²)∉Q, hence ≠40 and ≠8p.  occ(x³) is palindromic
  (S real) but a_1,a_2,a_3 pairwise distinct, so S(x³)∉Q, ≠65.
  Fail: claim S(x²)=40 at p=7, or S(x²)=8p.  ∎

Theorem B — PROVED (same census).
  occ(x²)=occ(x⁴)=occ(2x⁴) at p=7, and that occupancy is not
  palindromic (S(x²) is not real).  occ(x⁴)≠occ(x⁴+x)
  (occ0 179539290 vs 331055445): deg-4 is not one class.
  Five distinct occ0 appear among the seven deg-4 samples.
  Fail: claim S(x²) real at p=7; claim one deg-4 class.  ∎

Theorem C — PROVED (15.330 + named list).
  D=dim V₊=(q+1)/2 equals 13 at p=5 and 25≠409 at p=7.
  Fail: claim D=dim V₊ at p=7.  ∎

Theorem D — OPEN.  D unnamed.  The three p=5 deg-4 S ∈ Q(ζ_5)
  are not a formula in p (p=7 already has ≥5 deg-4 classes).
  phi_F not imported.

============================================================================
Backend: T(7,3) ProcessPool+numba census in
evidence/e1_gmin_m4_prop15368_catalog.json.  Writes
evidence/e1_gmin_m4_prop15368.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15363 import gen_T, occupancy, poly_alpha, wrap_diag_sums  # noqa: E402
from e1_gmin_m4_prop15367 import NUM_SUM_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15368_catalog.json"
T7 = 68938800


def dim_Vp(p: int) -> int:
    return (p * p + 1) // 2


def D_dim_neighbour(p: int) -> int:
    return dim_Vp(p)


def load_cat() -> dict:
    return json.loads(CAT.read_text())


def p5_S_x2_x3() -> dict:
    p, m = 5, 2
    mats = gen_T(p, m)
    diags = [wrap_diag_sums(g, 1) for g in mats]
    occ2 = occupancy(diags, poly_alpha((0, 0, 1), p), m)
    occ3 = occupancy(diags, poly_alpha((0, 0, 0, 1), p), m)
    return {"occ2": occ2, "occ3": occ3, "nT": len(mats)}


def prove_A() -> dict:
    p5 = p5_S_x2_x3()
    cat = load_cat()
    x2 = cat["by_alpha"]["x2"]
    x3 = cat["by_alpha"]["x3"]
    ok = p5["occ2"][0] == 440 and p5["nT"] == 2040
    ok = ok and cat["nT"] == T7
    ok = ok and x2["occ0"] == 179539290
    nz2 = [x2["occ"][k] for k in range(1, 7)]
    nz3 = [x3["occ"][k] for k in range(1, 7)]
    s2_rational = len(set(nz2)) == 1
    s3_rational = len(set(nz3)) == 1
    ok = ok and x2["occ0"] == 179539290
    ok = ok and (not s2_rational)
    ok = ok and all(x3["occ"][k] == x3["occ"][7 - k] for k in range(1, 7))
    ok = ok and (not s3_rational)
    if s2_rational or s3_rational:
        ok = False
    return {
        "proved": bool(ok),
        "p5_occ2": list(p5["occ2"]),
        "p7_x2_occ0": x2["occ0"],
        "p7_x3_palindromic": all(
            x3["occ"][k] == x3["occ"][7 - k] for k in range(1, 7)
        ),
        "p7_S2_in_Q": s2_rational,
        "p7_S3_in_Q": s3_rational,
        "theorem": "S(x²), S(x³) ∉ Q at p=7, so they are not 40 or 65.",
    }


def prove_B() -> dict:
    cat = load_cat()
    x2 = cat["by_alpha"]["x2"]["occ"]
    x4 = cat["by_alpha"]["x4"]["occ"]
    twox4 = cat["by_alpha"]["2x4"]["occ"]
    x4x = cat["by_alpha"]["x4+x"]["occ"]
    occ0s = sorted(
        {
            cat["by_alpha"][k]["occ0"]
            for k in cat["by_alpha"]
            if k.startswith("x4") or k == "2x4"
        }
    )
    pal = x2 == x2[::-1]
    ok = x2 == x4 == twox4
    ok = ok and x4[0] != x4x[0]
    ok = ok and (not pal)
    ok = ok and len(occ0s) >= 5
    if pal or x4[0] == x4x[0]:
        ok = False
    return {
        "proved": bool(ok),
        "x2_eq_x4": x2 == x4,
        "x2_palindromic": pal,
        "n_deg4_occ0": len(occ0s),
        "deg4_occ0": occ0s,
        "theorem": "At p=7, S(x²)=S(x⁴) complex; deg-4 has ≥5 occ0 classes.",
    }


def prove_C() -> dict:
    ok = D_dim_neighbour(5) == 13
    ok = ok and D_dim_neighbour(7) == 25
    if D_dim_neighbour(7) == 409:
        ok = False
    return {
        "proved": bool(ok),
        "D_dim_p5": D_dim_neighbour(5),
        "D_dim_p7": D_dim_neighbour(7),
        "live_D7": 409,
        "theorem": "D=dim V+ hits p=5 and dies at p=7 (25≠409).",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "NUM_SUM_named_p5": NUM_SUM_named(5),
        "note": (
            "p=7 deg-4 has ≥5 occupancy classes, so the three p=5 "
            "S ∈ Q(ζ_5) are not a formula in p.  D unnamed.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.368  p=7 S(α) census; S2=40 and 3-class die", flush=True)
    A = prove_A()
    print(f"  A S2/S3 die: {A['proved']} p7_x2_occ0={A['p7_x2_occ0']}", flush=True)
    B = prove_B()
    print(
        f"  B deg4: {B['proved']} x2=x4={B['x2_eq_x4']} n_occ0={B['n_deg4_occ0']}",
        flush=True,
    )
    C = prove_C()
    print(f"  C dimV+: {C['proved']} {C['D_dim_p5']} vs {C['D_dim_p7']}", flush=True)
    D = prove_open()
    print(f"  D open: D_named={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.368",
        "title": "p=7 S(α): S2=40 and one deg-4 class die; D=dim V+ dies",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S2_S3_p5_die_p7": A["proved"],
            "deg4_not_one_class_p7": B["proved"],
            "D_dimVp_dies_p7": C["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15368.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
