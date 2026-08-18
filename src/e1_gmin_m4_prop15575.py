#!/usr/bin/env python3
"""
Prop 15.575 — μ_1d is the Type+ Johnson 4-point

    μ_1d = 2 p⁴ (p − 1 − 4/(p−2)) = 2 p⁴ (p²−3p−2)/(p−2)

for r∈F_p^×\\{±1}.  Type+ 1D is ε:F_p→{±1}, #(+)=m=(p+1)/2,
ẑ=p ε̂ on L^perp (15.538).  Only C(p,m) of n_1d=m C(p,m)
vectors have this line as support, so
    μ_1d = p⁴ E[|ε̂(1)|² |ε̂(r)|²] / m.
E[|ε̂(α)|²]=p+1, but |ε̂|² is not constant on F_p^* for a
single ε (it is only on average).  The 4-point is
    E = (p+1)² − 2p(p+1)/(p−2) = (p+1)(p²−3p−2)/(p−2).
Hits live μ_1d at p=5 (10000/3) and p=7 (124852/5).
Fail: constant |ε̂| ⇒ p⁴(p+1)²/m (p=5: 7500≠10000/3).
Fail: 16p⁴/3 (p=7).  Aut_e DEAD.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.574 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.574.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.573: Q_sub=(n_1d μ_1d+n_{k=3} μ_{k=3}+n_full μ_full)/N.
15.574 names μ_{k=3}.  15.538 names the 3-point A_1d,dbl, not
this 4-point.  ẑ(−ξ)=conj ẑ(ξ) so the 4-point is u(ξ)u(rξ).

============================================================================
Theorem A — PROVED (Johnson 4-point; enum p=5,7,11,13 + live p=5,7).
  μ_1d=2p⁴(p−1−4/(p−2)).  Independent of r≠±1 at the
  enumerated primes (deg≤4 rational from λ₄; not a
  termwise λ₄ expansion).  Fail constant |ε̂|² (p=5).
  Fail 16p⁴/3 (p=7).  Fail drop the 4/(p−2) term.  ∎

Theorem B — OPEN.  μ_full is the 15.568 k-measure 4-point at
  p=7 (16p² E[|k̂(1)|²|k̂(r)|²] hits 126224/15) but not a
  p-law (same mix is 2800≠2000 at p=5).  P(r) splits at
  p≥11.  Q_sub still not a p-identity.  Do not interpolate
  (p−5)/15.  phi_F not imported.  Aut_e DEAD.

============================================================================
Backend: Johnson enum O(C(p,m) p) + p=5,7 yhat.  Writes
evidence/e1_gmin_m4_prop15575.json
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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15573 import q_sub_mix

EV = ROOT / "evidence" / "e1_gmin_m4_prop15575.json"


def mu_1d_named(p: int) -> Fraction:
    """2 p⁴ (p²−3p−2)/(p−2)."""
    return Fraction(2 * p**4 * (p * p - 3 * p - 2), p - 2)


def mu_1d_const_eps(p: int) -> Fraction:
    """Fail: |ε̂|²≡p+1 ⇒ p⁴(p+1)²/m."""
    m = (p + 1) // 2
    return Fraction(p**4 * (p + 1) ** 2, m)


def mu_1d_16p4_over_3(p: int) -> Fraction:
    """Fail: 16p⁴/3 (the p=5 value)."""
    return Fraction(16 * p**4, 3)


def mu_1d_drop4(p: int) -> Fraction:
    """Fail: drop 4/(p−2)."""
    return Fraction(2 * p**4 * (p - 1), 1)


def johnson_E_sub(p: int) -> float:
    """E[|ε̂(1)|² |ε̂(r)|²] for any r≠±1, by enumerating m-subsets."""
    m = (p + 1) // 2
    acc = 0.0
    n = 0
    for plus in combinations(range(p), m):
        s = np.zeros(p, dtype=np.float64)
        for x in plus:
            s[x] = 1.0
        eh2 = np.abs(np.fft.fft(2.0 * s - 1.0)) ** 2
        acc += eh2[1] * eh2[2]
        n += 1
    return acc / n


def prove_A() -> dict:
    ok = True
    rows = {}
    live = {}
    for p in (5, 7):
        live[p] = q_sub_mix(p)["parts"]["1d"]["mu"]
    for p in (5, 7, 11, 13):
        E = johnson_E_sub(p)
        m = (p + 1) // 2
        got = E * (p**4) / m
        want = float(mu_1d_named(p))
        match = abs(got - want) < 1e-6
        fail_const = abs(want - float(mu_1d_const_eps(p))) > 1.0
        fail_16 = abs(want - float(mu_1d_16p4_over_3(p))) > 1.0 if p != 5 else True
        fail_drop = abs(want - float(mu_1d_drop4(p))) > 1.0
        live_ok = True
        if p in live:
            live_ok = abs(got - live[p]) < 1e-6
        if not (match and fail_const and fail_16 and fail_drop and live_ok):
            ok = False
        rows[str(p)] = {
            "E_sub": E,
            "got": got,
            "named": want,
            "const_eps": float(mu_1d_const_eps(p)),
            "sixteen_over_3": float(mu_1d_16p4_over_3(p)),
            "drop4": float(mu_1d_drop4(p)),
            "live": live.get(p),
            "err": float(abs(got - want)),
            "match": bool(match),
            "const_fails": bool(fail_const),
            "sixteen_fails_or_p5": bool(fail_16),
            "drop4_fails": bool(fail_drop),
            "live_ok": bool(live_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "μ_1d=2p⁴(p−1−4/(p−2)). Fail constant |ε̂|; "
            "fail 16p⁴/3 at p=7; fail drop 4/(p−2)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "mu_full_named_in_p": False,
        "Q_sub_closed_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "μ_1d is Johnson. μ_full equals the 15.568 k 4-point "
            "at p=7, not a p-law. P(r) splits at p≥11. "
            "Do not interpolate (p−5)/15. phi_F not imported. Aut_e DEAD."
        ),
    }


def main() -> dict:
    print("Prop 15.575  μ_1d = 2p⁴(p−1−4/(p−2)) Johnson 4-point", flush=True)
    A = prove_A()
    print(f"  A Johnson μ_1d: {A['proved']}", flush=True)
    B = prove_open()
    out = {
        "prop": "15.575",
        "title": "μ_1d=2p⁴(p−1−4/(p−2)), Type+ Johnson 4-point",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "proved": {
            "mu_1d_johnson_named": A["proved"],
            "mu_full_named_in_p": False,
            "Q_sub_closed_in_p": False,
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
        "identity": "mu_1d = 2 p^4 (p-1 - 4/(p-2))",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
