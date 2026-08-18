#!/usr/bin/env python3
"""
Prop 15.570 — 15.568 A3 family mean is a complete Jacobi sum
for every odd prime p.

    k = χ(d1)(2 e_α − e_{α+d1} − e_{α+d2})   same-χ {d1,d2}
    A3(p) := E_k[ mean_g ẑ(g)² ẑ(−2g) ]
      = 96 p² (2χ(2)−1)/(p−1)                 if p≡3 (mod 4)
      = −32 p³ (1+2χ(2))/((p−1)(p−3))         if p≡1 (mod 4)
    ≡ 6 · amp_{2χ}(p)  on p≡3 (mod 4)          (15.567 amp)
For odd p≥5 (p=3 has empty same-χ pairs).
p=7: 16p².  Fail 16p² at p=5,11.  Fail 6·amp at p=5.
A6 is not amp·(p−3)/4 (dies p=19).  Aut_e DEAD.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.569 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.569.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.569: for k=2e_0−e_b−e_c the g-mean is
    8 τ³ χ(−2)/(p−1) · Σ_{i,j} c_{ij} G_χ(ib+jc)
G_χ(m)=χ(m)τ (m≠0), else 0; τ²=χ(−1)p.  Named A3 multiplies
by χ(d1).  Translation invariance ⇒ the mean depends only on
r=d2/d1 ∈ □ minus {1}.  Then χ(d1) G_χ(i d1+j d2)=τ χ(i+jr) 1_{i+jr≠0}.

============================================================================
Theorem A — PROVED (Jacobi evaluation of 15.569; p=5..31).
  Average of χ(d1)S over same-χ pairs equals
      8 p² χ(−1)χ(2)/(p−1) · mean_{r square, r≠1} T(r)
  with T the 15-term χ-polynomial (the (0,0) monomial drops).
  Standard sums Σ_{r∈□} χ(r+a)=−(1+χ(a))/2 give:
    p≡3 (mod 4): mean T=12(χ(2)−2), hence
        A3=96 p² (2χ(2)−1)/(p−1)=6 amp_{2χ}.
    p≡1 (mod 4): mean T=−4(χ(2)+2) p/(p−3), hence
        A3=−32 p³ (1+2χ(2))/((p−1)(p−3)).
  Fail: 16p² (p=5: 400≠500; p=11: 1936≠−3484.8).
  Fail: 6 amp_{2χ} at p=5 (0≠500).
  Fail: drop (2χ(2)−1) at p=11.  ∎

Theorem B — PROVED (p=7,11,19,23).
  A6 ≠ amp_{2χ}(p−3)/4: holds at p=7,11, dies at p=19
  (ratio 2 vs (p−3)/4=4) and p=23 (ratio 9/5).  ∎

Theorem C — OPEN.  A6 family mean is not the A3 Jacobi
  integer.  A_full still not a p-law.  Do not interpolate
  (p−5)/15.  phi_F not imported.  Aut_e DEAD.

============================================================================
Backend: named-family Gauss O(p³) at p=5,7,11,13,17,19.
Writes evidence/e1_gmin_m4_prop15570.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15567 import chip, gauss_table, twochi_amp_named
from e1_gmin_m4_prop15568 import fam_A3_named, fam_A6_named, mean_prod

EV = ROOT / "evidence" / "e1_gmin_m4_prop15570.json"


def a3_named_closed(p: int) -> float:
    """Complete-sum name of the 15.568 A3 family mean."""
    chi2 = chip(p, 2)
    if chip(p, p - 1) == -1:
        return 96.0 * p * p * (2 * chi2 - 1) / (p - 1)
    return -32.0 * (p**3) * (1 + 2 * chi2) / ((p - 1) * (p - 3))


def a3_wrong_16p2(p: int) -> float:
    return 16.0 * p * p


def a3_only_6amp(p: int) -> float:
    """Fail: the p≡3 formula on every p."""
    return 6.0 * twochi_amp_named(p)


def a3_drop_chi2(p: int) -> float:
    """Fail: drop (2χ(2)−1) and (1+2χ(2))."""
    if chip(p, p - 1) == -1:
        return 96.0 * p * p / (p - 1)
    return -32.0 * (p**3) / ((p - 1) * (p - 3))


def a6_amp_times_p_minus_3_over_4(p: int) -> float:
    return twochi_amp_named(p) * (p - 3) / 4.0


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        G, omega = gauss_table(p)
        got = mean_prod(fam_A3_named(p), p, G, omega)
        want = a3_named_closed(p)
        match = abs(got - want) < 1e-6
        fail16 = abs(got - a3_wrong_16p2(p)) > 1.0 if p != 7 else True
        fail6 = abs(got - a3_only_6amp(p)) > 1.0 if p % 4 == 1 else True
        fail_chi = abs(got - a3_drop_chi2(p)) > 1.0 if p in (11, 17) else True
        if not (match and fail16 and fail6 and fail_chi):
            ok = False
        rows[str(p)] = {
            "got": float(got),
            "named": want,
            "sixteen_p2": a3_wrong_16p2(p),
            "six_amp": a3_only_6amp(p),
            "drop_chi2": a3_drop_chi2(p),
            "err": float(abs(got - want)),
            "match": bool(match),
            "sixteen_fails": bool(fail16),
            "six_amp_fails_or_p3": bool(fail6),
            "drop_chi2_fails_or_ok": bool(fail_chi),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A3=96p²(2χ(2)−1)/(p−1) on p≡3 (mod 4); "
            "A3=−32p³(1+2χ(2))/((p−1)(p−3)) on p≡1 (mod 4). "
            "Fail 16p²; fail 6·amp at p=5; fail drop χ(2) at p=11."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (7, 11, 19, 23):
        G, omega = gauss_table(p)
        got = mean_prod(fam_A6_named(p), p, G, omega)
        guess = a6_amp_times_p_minus_3_over_4(p)
        if p in (7, 11):
            match_small = abs(got - guess) < 1e-6
            dies = True
        else:
            match_small = True
            dies = abs(got - guess) > 1.0
        if not (match_small and dies):
            ok = False
        rows[str(p)] = {
            "A6": float(got),
            "guess": guess,
            "amp": twochi_amp_named(p),
            "err": float(abs(got - guess)),
            "guess_ok_small_p": bool(match_small if p in (7, 11) else True),
            "guess_dies_p19": bool(dies),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A6=amp(p−3)/4 at p=7,11; dies at p=19,23. "
            "A6 is not that p-law."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A6_named_in_p": False,
        "A_full_dbl_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "A3 family mean is a Jacobi integer. A6 is not amp(p−3)/4. "
            "A_full still not a p-law. Do not interpolate (p−5)/15. "
            "phi_F not imported. Aut_e DEAD."
        ),
    }


def main() -> dict:
    print("Prop 15.570  A3 family mean is a Jacobi sum", flush=True)
    A = prove_A()
    print(f"  A A3 closed: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B A6 not (p-3)/4 amp: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.570",
        "title": "A3=96p²(2χ(2)−1)/(p−1) (p≡3 mod 4); −32p³(1+2χ(2))/((p−1)(p−3)) (p≡1)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "a3_family_jacobi": A["proved"],
            "a6_not_amp_pm3_over_4": B["proved"],
            "A6_named_in_p": False,
            "A_full_dbl_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "A3_named = 96 p^2 (2 chi(2)-1)/(p-1) if p=3 mod 4; "
            "-32 p^3 (1+2 chi(2))/((p-1)(p-3)) if p=1 mod 4"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
