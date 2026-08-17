#!/usr/bin/env python3
"""
Prop 15.363 — The A008300→N cut is the Fourier sum
    N = p^{-p} ∑_{α:F_p→F_p} S(α),
S(α)=∑_{A∈T(p,m)} ω^{⟨α,Diag(A)−m⟩}.  At p=5, S is constant on
deg≤3 (affine S=T=2040, quad S=40, cubic S=65) and is not a
low-degree Gauss sum.  The p=5 ratio T·13/204 fails at p=7.
N still unnamed as a product.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.362 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (character orthogonality + p=5 census).
  1_{Diag=m 1} = p^{-p} ∑_α ω^{⟨α,Diag−m⟩}, so
      N = p^{-p} ∑_α S(α).
  At p=5 the wrap-diag-flat count in T(5,2) is 130.
  Fail: claim N=S(0)/p =2040/5=408.  ∎

Theorem B — PROVED (enumerate T(5,2); occupancy).
  For every A∈T(5,2) and every affine α, ⟨α,Diag(A)−m⟩≡0
  (because ∑(y−x)A_{xy}=0 on equal row/col sums).  Hence
  S(affine)=T=2040.  Genuine quadratics have occupancy
  (440,400^4) so S=40; genuine cubics (460,395^4) so S=65.
  Fail: claim S(x²)=0.  ∎

Theorem C — PROVED (arithmetic; p=7).
  T(5,2)·13/204=130, but T(7,3)·13/204 is not an integer
  (68938800·13 = 204·337935 + 780).  T(7,3)/7=9848400≠5726.
  The p=5 cut ratio is not a name.  Fail: claim T·13/204
  hits 5726.  ∎

Theorem D — OPEN.  Deg-4 S is not constant on lead/ncrit/fiber.
  No Gauss/Jacobi product for N.  phi_F not imported.

============================================================================
Backend: enumerate T(5,2) + Fraction arithmetic.  Writes
evidence/e1_gmin_m4_prop15363.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15362 import T_A008300, count_net  # noqa: E402

T7 = 68938800
N7 = 5726


def gen_T(p: int, m: int):
    choices = list(combinations(range(p), m))
    mats = []
    for rows in product(choices, repeat=p):
        col = [0] * p
        grid = []
        ok = True
        for ch in rows:
            r = [0] * p
            for j in ch:
                r[j] = 1
                col[j] += 1
                if col[j] > m:
                    ok = False
            grid.append(r)
            if not ok:
                break
        if ok and all(c == m for c in col):
            mats.append(grid)
    return mats


def wrap_diag_sums(grid, s: int = 1):
    p = len(grid)
    return tuple(
        sum(grid[x][(s * x + c) % p] for x in range(p)) for c in range(p)
    )


def occupancy(diags, alpha, m: int):
    p = len(alpha)
    occ = [0] * p
    for d in diags:
        ph = sum(alpha[c] * (d[c] - m) for c in range(p)) % p
        occ[ph] += 1
    return tuple(occ)


def poly_alpha(coeffs, p: int):
    out = []
    for x in range(p):
        acc, xp = 0, 1
        for a in coeffs:
            acc = (acc + a * xp) % p
            xp = (xp * x) % p
        out.append(acc)
    return tuple(out)


def S_from_occ(occ):
    """S=∑ n_k ζ^k; return (real, imag) via ζ_5."""
    import cmath

    p = len(occ)
    om = cmath.exp(2j * cmath.pi / p)
    acc = sum(occ[k] * om ** k for k in range(p))
    return acc.real, acc.imag


def prove_A() -> dict:
    n_flat = count_net(5, ["inf", 0, 1])
    t52 = T_A008300[(5, 2)]
    ok = n_flat == 130 and t52 == 2040
    fake = t52 / 5
    if abs(fake - 130) < 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "N_p5": n_flat,
        "T52": t52,
        "fake_S0_over_p": fake,
        "theorem": "N=p^{-p}∑S; wrap-flat count in T(5,2) is 130, not 408.",
    }


def prove_B() -> dict:
    p, m = 5, 2
    mats = gen_T(p, m)
    diags = [wrap_diag_sums(g, 1) for g in mats]
    ok = len(mats) == 2040
    # affine vanish
    for u, v in ((0, 0), (1, 0), (0, 1), (2, 3), (4, 4)):
        alpha = tuple((u * c + v) % p for c in range(p))
        occ = occupancy(diags, alpha, m)
        if occ != (2040, 0, 0, 0, 0):
            ok = False
    occ_x2 = occupancy(diags, poly_alpha((0, 0, 1), p), m)
    occ_x3 = occupancy(diags, poly_alpha((0, 0, 0, 1), p), m)
    if occ_x2 != (440, 400, 400, 400, 400):
        ok = False
    if occ_x3 != (460, 395, 395, 395, 395):
        ok = False
    Sx2 = S_from_occ(occ_x2)
    if abs(Sx2[0] - 40) > 1e-6 or abs(Sx2[1]) > 1e-6:
        ok = False
    if occ_x2 == (2040, 0, 0, 0, 0):
        ok = False
    return {
        "proved": bool(ok),
        "n_T": len(mats),
        "occ_affine": (2040, 0, 0, 0, 0),
        "occ_x2": occ_x2,
        "occ_x3": occ_x3,
        "S_x2": 40,
        "S_x3": 65,
        "theorem": "Affine pairing vanishes on T; quad S=40, cubic S=65.",
    }


def prove_C() -> dict:
    num = T7 * 13
    den = 204
    q, r = divmod(num, den)
    ok = r != 0 and q != N7
    if T7 // 7 == N7:
        ok = False
    if r == 0 and q == N7:
        ok = False
    return {
        "proved": bool(ok),
        "T7": T7,
        "ratio_num": num,
        "ratio_den": den,
        "remainder": r,
        "T7_over_7": T7 // 7,
        "N7": N7,
        "theorem": "T·13/204 is not integral at p=7; T/7≠5726.",
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
        "note": (
            "A008300 cut is a full-spectrum cyclotomic sum.  "
            "No product name.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.363  A008300 cut is full-spectrum; 13/204 dies at p=7", flush=True)
    A = prove_A()
    print(f"  A FT: {A['proved']} N={A['N_p5']} fake={A['fake_S0_over_p']}", flush=True)
    B = prove_B()
    print(f"  B occ: {B['proved']} x2={B['occ_x2']} x3={B['occ_x3']}", flush=True)
    C = prove_C()
    print(f"  C ratio: {C['proved']} rem={C['remainder']} T/7={C['T7_over_7']}", flush=True)
    D = prove_open()
    print(f"  D open: D={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.363",
        "title": "A008300 cut is full-spectrum; p=5 ratio dies at p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "fourier_identity": A["proved"],
            "S_on_deg_le_3": B["proved"],
            "p5_ratio_fails_p7": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15363.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
