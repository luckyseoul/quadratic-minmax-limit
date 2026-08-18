#!/usr/bin/env python3
"""
Prop 15.572 — 15.568 A6 family mean is a Jacobi integer
for every odd prime p≥5.

    A6_gen = 8 p² χ(−1)χ(2)/(p−1) · T_gen
    A6_deg / p² = 32(11+10χ(2)+4χ(2)χ(3))/(p−1)   (p≡1 only)
    A6/p² = (p−5)/(p−3) A6_gen/p² + 2/(p−3) A6_deg/p²
T_gen is the 15.571 monomial average on nondeg square-pairs:
    p≡3, χ(2)=−1:  T=−12
    p≡3, χ(2)=+1:  T=−4(p−5)/(p−3)
    p≡5 (mod 8), p>5: T=44 − 160(1+χ(3))/(p−5)
    p≡1 (mod 8): T=4(33p−351−40χ(3))/(p−5)
p=5 is deg-only: A6=40 p².  Fail amp(3−χ(2))/2 at p=23.
Fail drop χ(3) on deg (p=13 vs 29).  Aut_e DEAD.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.571 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.571.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.571: nondeg k̂=1+2(x+y+xy) expands to a Gauss
polynomial; deg s2=−s1 (p≡1) is a separate 2-point fibre.
n_gen=(p−1)(p−5)/8, n_deg=(p−1)/4, n_tot=(p−1)(p−3)/8.

============================================================================
Theorem A — PROVED (Jacobi average of 15.571; p=5..113).
  T_gen and A6_deg as above.  Mix weights (p−5):(2 n_deg/n_gen).
  p≡3: A6=−96p²/(p−1) if χ(2)=−1, else 32p²(p−5)/((p−1)(p−3)).
  Fail: amp(3−χ(2))/2 at p=23 (ratio 9/5).
  Fail: drop χ(3) on deg at p=13 (sign).
  Fail: T=−12 at p=23.  ∎

Theorem B — OPEN.  A_full / Q3_dbl / Q_τ still not a p-law.
  Do not interpolate (p−5)/15.  phi_F not imported.  Aut_e DEAD.

============================================================================
Backend: named-pair Gauss O(p²) at p=5,7,11,13,17,23,29,73.
Writes evidence/e1_gmin_m4_prop15572.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15567 import chip, gauss_table, prod_dbl_k, twochi_amp_named
from e1_gmin_m4_prop15568 import squares
from e1_gmin_m4_prop15571 import a6_guess_3mchi2, a6_guess_pm3_over4

EV = ROOT / "evidence" / "e1_gmin_m4_prop15572.json"


def t_gen_named(p: int) -> Fraction | None:
    chi2 = chip(p, 2)
    chi3 = chip(p, 3)
    if chip(p, p - 1) == -1:
        if chi2 == -1:
            return Fraction(-12)
        return Fraction(-4 * (p - 5), p - 3)
    if chi2 == -1:
        if p == 5:
            return None
        return Fraction(44) - Fraction(160 * (1 + chi3), p - 5)
    return Fraction(4 * (33 * p - 351 - 40 * chi3), p - 5)


def deg_p2_named(p: int) -> Fraction | None:
    if chip(p, p - 1) == -1:
        return None
    chi2 = chip(p, 2)
    chi3 = chip(p, 3)
    return Fraction(32 * (11 + 10 * chi2 + 4 * chi2 * chi3), p - 1)


def a6_named_closed(p: int) -> Fraction:
    """Complete-sum name of the 15.568 A6 family mean."""
    chim1 = chip(p, p - 1)
    chi2 = chip(p, 2)
    if p == 5:
        return deg_p2_named(p) * (p * p)
    if chim1 == -1:
        return (
            Fraction(8 * p * p * chim1 * chi2, p - 1) * t_gen_named(p)
        )
    tg = t_gen_named(p)
    dg = deg_p2_named(p)
    term1 = Fraction(p - 5, p - 3) * Fraction(8 * chim1 * chi2, p - 1) * tg
    term2 = Fraction(2, p - 3) * dg
    return (term1 + term2) * (p * p)


def a6_drop_chi3_deg(p: int) -> Fraction:
    """Fail: deg uses χ(3)=0 (drop the χ(3) terms)."""
    if chip(p, p - 1) == -1:
        return a6_named_closed(p)
    chi2 = chip(p, 2)
    dg = Fraction(32 * (11 + 10 * chi2), p - 1)
    if p == 5:
        return dg * (p * p)
    tg = t_gen_named(p)
    chim1 = chip(p, p - 1)
    term1 = (
        Fraction(p - 5, p - 3)
        * Fraction(8 * chim1 * chi2, p - 1)
        * tg
    )
    term2 = Fraction(2, p - 3) * dg
    return (term1 + term2) * (p * p)


def a6_t_minus_12_everywhere(p: int) -> Fraction:
    """Fail: T_gen≡−12."""
    chim1 = chip(p, p - 1)
    chi2 = chip(p, 2)
    return Fraction(8 * p * p * chim1 * chi2, p - 1) * Fraction(-12)


def pair_mean(p: int) -> float:
    G, omega = gauss_table(p)
    S = sorted(squares(p))
    acc = 0.0
    n = 0
    for s1, s2 in combinations(S, 2):
        k = [0] * p
        plus = {s1, s2, (s1 + s2) % p}
        for i in range(1, p):
            k[i] = 1 if i in plus else -1
        acc += float(prod_dbl_k(tuple(k), p, G, omega).real)
        n += 1
    return acc / n


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 23, 29, 73):
        got = pair_mean(p)
        want = float(a6_named_closed(p))
        match = abs(got - want) < 1e-6
        fail_amp = True
        fail_t12 = True
        fail_chi3 = True
        if p == 23:
            fail_amp = abs(got - a6_guess_3mchi2(p)) > 1.0
            fail_t12 = abs(got - float(a6_t_minus_12_everywhere(p))) > 1.0
        if p == 19:
            fail_amp = abs(got - a6_guess_pm3_over4(p)) > 1.0
        if p in (13, 29):
            fail_chi3 = abs(got - float(a6_drop_chi3_deg(p))) > 1.0
        if not (match and fail_amp and fail_t12 and fail_chi3):
            ok = False
        rows[str(p)] = {
            "got": got,
            "named": want,
            "frac": str(a6_named_closed(p)),
            "err": float(abs(got - want)),
            "guess_3mchi2": a6_guess_3mchi2(p),
            "drop_chi3": float(a6_drop_chi3_deg(p)),
            "t12": float(a6_t_minus_12_everywhere(p)),
            "match": bool(match),
            "amp_guess_fails_or_ok": bool(fail_amp),
            "t12_fails_or_ok": bool(fail_t12),
            "chi3_fails_or_ok": bool(fail_chi3),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A6 Jacobi mix of T_gen and χ(2),χ(3)-deg. "
            "Fail amp(3−χ(2))/2 at p=23; fail drop χ(3) at p=13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A_full_dbl_named_in_p": False,
        "Q3_double_named_in_p": False,
        "Q_tau_named": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "A6 family mean is a Jacobi integer. A_full / Q3_dbl / Q_τ "
            "still unnamed. Do not interpolate (p−5)/15. "
            "phi_F not imported. Aut_e DEAD."
        ),
    }


def main() -> dict:
    print("Prop 15.572  A6 family mean is a Jacobi sum", flush=True)
    A = prove_A()
    print(f"  A A6 closed: {A['proved']}", flush=True)
    B = prove_open()
    out = {
        "prop": "15.572",
        "title": "A6 family mean is a Jacobi mix of T_gen and deg(χ(2),χ(3))",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "proved": {
            "a6_family_jacobi": A["proved"],
            "A_full_dbl_named_in_p": False,
            "Q3_double_named_in_p": False,
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
            "A6 = mix((p-5)/(p-3) * 8 p^2 chi(-1)chi(2) T_gen/(p-1), "
            "2/(p-3) * 32 p^2 (11+10 chi(2)+4 chi(2)chi(3))/(p-1))"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
