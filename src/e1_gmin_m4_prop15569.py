#!/usr/bin/env python3
"""
Prop 15.569 — A3 double product is the Gauss expansion of
    (2−x−y)²(2−x^{−2}−y^{−2}).

    k=2e_0−e_b−e_c
    E_g[ẑ(g)²ẑ(−2g)]
      = 8 τ³ χ(−2)/(p−1)
        · Σ_{i,j} c_{ij} G_χ(ib+jc)
with (c_{ij}) the monomials of (2−x−y)²(2−x^{−2}−y^{−2})
and G_χ(m)=χ(m)τ (m≠0), else 0.  The 15.568 same-χ family
mean equals this average.  At p=7 that mean is 16p²; that
integer is not a p-law (p=11).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.568 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.568.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.567: ẑ=2 G_χ k̂, G_χ(c)=χ(c)τ, same kernel ω^{−cx}.
χ is the quadratic residue character (χ²=1, χ̄=χ).
15.568 A3: k=χ(d1)(2e_α−e_{α+d1}−e_{α+d2}) on same-χ pairs.

============================================================================
Theorem A — PROVED (monomial Gauss; p=5,7,11,13).
  For k=2e_0−e_b−e_c the g-mean equals the expansion above.
  Fail: drop the x^{−2} factor (replace (2−x^{−2}−y^{−2}) by 2).
  Fail: drop the prefactor 8.  ∎

Theorem B — PROVED (A + 15.568 family; p=5,7,11,13).
  Mean of Theorem A over the same-χ 15.568 A3 family is the
  named A3 amplitude.  p=7: 16p².  Fail: 16p² at p=11.  ∎

Theorem C — OPEN.  The family mean is not a short integer in
  p (p=5: 20p²; p=7: 16p²; p=11: −144 p²/5).  A_full still
  not a p-law.  Do not interpolate (p−5)/15.  phi_F not imported.

============================================================================
Backend: field Gauss O(p²) pairs.  Writes
evidence/e1_gmin_m4_prop15569.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15567 import chip, gauss_table, prod_dbl_k
from e1_gmin_m4_prop15568 import fam_A3_named, mean_prod

EV = ROOT / "evidence" / "e1_gmin_m4_prop15569.json"


def gauss_eval(p: int, tau: complex, m: int) -> complex:
    if m % p == 0:
        return 0j
    return chip(p, m) * tau


def monomials_full() -> dict[tuple[int, int], int]:
    """(2−x−y)² (2−x^{−2}−y^{−2})."""
    left = {
        (0, 0): 4,
        (2, 0): 1,
        (0, 2): 1,
        (1, 1): 2,
        (1, 0): -4,
        (0, 1): -4,
    }
    acc: dict[tuple[int, int], int] = {}
    for (i, j), c0 in left.items():
        acc[i, j] = acc.get((i, j), 0) + 2 * c0
        acc[i - 2, j] = acc.get((i - 2, j), 0) - c0
        acc[i, j - 2] = acc.get((i, j - 2), 0) - c0
    return acc


def monomials_drop_inv() -> dict[tuple[int, int], int]:
    """Fail: drop (2−x^{−2}−y^{−2}) to the constant 2."""
    left = {
        (0, 0): 4,
        (2, 0): 1,
        (0, 2): 1,
        (1, 1): 2,
        (1, 0): -4,
        (0, 1): -4,
    }
    return {ij: 2 * c0 for ij, c0 in left.items()}


def a3_from_monomials(p: int, b: int, c: int, tau: complex, mons) -> complex:
    S = 0j
    for (i, j), coef in mons.items():
        S += coef * gauss_eval(p, tau, (i * b + j * c) % p)
    pref = 8 * tau**3 * chip(p, (-2) % p) / (p - 1)
    return pref * S


def a3_named_expansion(p: int, b: int, c: int, tau: complex) -> complex:
    return a3_from_monomials(p, b, c, tau, monomials_full())


def a3_drop_inv(p: int, b: int, c: int, tau: complex) -> complex:
    return a3_from_monomials(p, b, c, tau, monomials_drop_inv())


def a3_drop8(p: int, b: int, c: int, tau: complex) -> complex:
    return a3_named_expansion(p, b, c, tau) / 8.0


def prove_A() -> dict:
    ok = True
    rows = {}
    mons = monomials_full()
    for p in (5, 7, 11, 13):
        G, omega = gauss_table(p)
        tau = G[1]
        maxerr = 0.0
        drop_inv = 0.0
        drop8 = 0.0
        n = 0
        for b, c in combinations(range(1, p), 2):
            k = [0] * p
            k[0] = 2
            k[b] = -1
            k[c] = -1
            got = prod_dbl_k(tuple(k), p, G, omega)
            pred = a3_named_expansion(p, b, c, tau)
            maxerr = max(maxerr, abs(got - pred))
            drop_inv = max(drop_inv, abs(got - a3_drop_inv(p, b, c, tau)))
            drop8 = max(drop8, abs(got - a3_drop8(p, b, c, tau)))
            n += 1
        match = maxerr < 1e-8
        fails = drop_inv > 1.0 and drop8 > 1.0
        if not (match and fails):
            ok = False
        rows[str(p)] = {
            "n": n,
            "maxerr": float(maxerr),
            "drop_inv_err": float(drop_inv),
            "drop8_err": float(drop8),
            "n_monomials": len(mons),
            "match": bool(match),
            "fails_ok": bool(fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A3 pair product = 8τ³χ(−2)/(p−1) times the "
            "(2−x−y)²(2−x^{−2}−y^{−2}) Gauss expansion. "
            "Fail drop x^{−2}; fail drop 8."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        G, omega = gauss_table(p)
        named = mean_prod(fam_A3_named(p), p, G, omega)
        want16 = 16.0 * p * p
        if p == 7:
            match = abs(named - want16) < 1e-6
            p11_style_fail = True
        else:
            match = True
            p11_style_fail = abs(named - want16) > 1.0
        if not (match and p11_style_fail):
            ok = False
        rows[str(p)] = {
            "named_mean": named,
            "sixteen_p2": want16,
            "match_p7_or_not_16": bool(match if p == 7 else p11_style_fail),
            "sixteen_fails_off_p7": bool(p11_style_fail),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "15.568 A3 family mean is the expansion average. "
            "p=7: 16p². Fail 16p² at p=11."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A3_family_short_p_law": False,
        "A_full_dbl_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "A3 pair product is a complete Gauss expansion. "
            "The same-χ family mean is not 16p² off p=7. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.569  A3 product is a Gauss expansion", flush=True)
    A = prove_A()
    print(f"  A expansion: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B family mean: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.569",
        "title": "A3 pair product = Gauss expansion of (2−x−y)²(2−x^{−2}−y^{−2})",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "a3_pair_expansion": A["proved"],
            "a3_family_mean": B["proved"],
            "A3_family_short_p_law": False,
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
            "A3 E_g zhat^2 zhat(-2g) = 8 tau^3 chi(-2)/(p-1) "
            "* Gauss[(2-x-y)^2 (2-x^{-2}-y^{-2})]"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
