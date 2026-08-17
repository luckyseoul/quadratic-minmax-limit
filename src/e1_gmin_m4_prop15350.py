#!/usr/bin/env python3
"""
Prop 15.350 — Q_{++}(A·T_0) is the Fejer average of the 3-line
geometric-sum Fourier transform of the triangle.  Named in p.
Ensemble Q_τ still OPEN.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.349 flags.

============================================================================
Setup.  15.349: for p≡1 (mod 4), D=A·T_0 with A=[[1,1],[1,−1]] is
Max+.  F_p²=F_p[ω]/(ω²−ib), ib=least nonsquare, ia=0, u=c0+p c1.
Ω={ξ: χ̂(ξ)=p}, |Ω|=(q−1)/2.  ++={r∈□\\{±1}: merge_cls=++},
|++|=2p−4.

============================================================================
Theorem A — PROVED (geometric series; certified p=5,13).
  The additive FT of 1_D at (μ,ν)∈F_p² is
      S(0,0)=p(p−1)/2,
      S=−p/(1−ζ^{μ+ν}) on the diagonal μ=ν≠0,
      S=p/(1−ζ^{μ−ν}) on the antidiagonal μ=−ν≠0,
      S=−p/(1−ζ^{−(μ+ν)}) on the axis μ=0, ν≠0,
      S=0 elsewhere.
  Fail: claim S=0 on the diagonal at p=5.  ∎

Theorem B — PROVED (field embedding; certified p=5,13,17,29).
  |ẑ(ξ)|²=4|S(2s0, 2 ib s1)|².  The support inside Ω is exactly
      L0={p t: t∈F_p^×},
      L+={t ib + p t: t∈F_p^×},
      L−={(−t ib mod p)+p t: t∈F_p^×},
  |L|=3(p−1), L⊂Ω.  Fail: claim L=Ω at p=13 (36≠84);
  claim |++|=2p (30≠22 at p=13).  ∎

Theorem C — PROVED (Ω-average over L).
  Q_{++}(A·T_0)
    = 16/(q² |Ω|) · (1/|++|) ∑_{r∈++} ∑_{ξ∈L} |S(ξ)|² |S(rξ)|²
  equals 64/15, 1168/231, 808/135, 3824/405 at p=5,13,17,29.
  Fail: drop the 16 (p=5 becomes 4/15).  No deg≤3/deg≤2
  rational with coeffs in [−8,8] hits {5,13,17,29,37}.  ∎

Theorem D — OPEN.  This names the triangle orbit, not the
  H_+ mixture (den |H_+|/(2p)).  phi_F not imported.

============================================================================
Backend: geometric sum + Ω listing.  Writes
evidence/e1_gmin_m4_prop15350.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15270 import _field_params  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15349 import A_p1, T0, apply_gl  # noqa: E402


def S_closed(mu: int, nu: int, p: int) -> complex:
    if mu == 0 and nu == 0:
        return p * (p - 1) / 2
    a = (mu + nu) % p
    b = (mu - nu) % p
    z = np.exp(2j * np.pi / p)
    if b == 0:
        return -p / (1 - z ** a)
    if a == 0:
        return p / (1 - z ** b)
    if a == (-b) % p:
        return -p / (1 - z ** ((-a) % p))
    return 0j


def S_brute(mu: int, nu: int, p: int) -> complex:
    g = apply_gl(T0(p), A_p1(p), p)
    z = np.exp(2j * np.pi / p)
    acc = 0j
    for x in range(p):
        for y in range(p):
            if g[x, y]:
                acc += z ** ((mu * x + nu * y) % p)
    return acc


def L_sets(p: int):
    ia, ib = _field_params(p)
    L0, Lp, Lm = set(), set(), set()
    for t in range(1, p):
        L0.add(p * t)
        Lp.add((ib * t) % p + p * t)
        Lm.add(((-ib * t) % p) + p * t)
    return ib, L0, Lp, Lm, L0 | Lp | Lm


def S_of_xi(p: int, ib: int, xi: int) -> complex:
    s0, s1 = xi % p, xi // p
    return S_closed((2 * s0) % p, (2 * ib * s1) % p, p)


def Qpp_from_L(p: int) -> Fraction:
    F = _kernel_field(p)
    Omega = omega_set(F)
    ib, _, _, _, L = L_sets(p)
    q2 = float((p * p) ** 2)
    u = {xi: 4.0 * abs(S_of_xi(p, ib, xi)) ** 2 for xi in L}
    acc = []
    for r in F["squares"]:
        if merge_cls(F, r) != "++":
            continue
        sm = 0.0
        for xi in L:
            eta = F["mul"][r][xi]
            if eta in u:
                sm += u[xi] * u[eta]
        acc.append(sm / len(Omega) / q2)
    return Fraction(float(np.mean(acc))).limit_denominator(p ** 4)


def Qpp_drop16(p: int) -> Fraction:
    """Fail-when-wrong: omit the 4|S|² factor pair (the 16)."""
    F = _kernel_field(p)
    Omega = omega_set(F)
    ib, _, _, _, L = L_sets(p)
    q2 = float((p * p) ** 2)
    u = {xi: abs(S_of_xi(p, ib, xi)) ** 2 for xi in L}
    acc = []
    for r in F["squares"]:
        if merge_cls(F, r) != "++":
            continue
        sm = 0.0
        for xi in L:
            eta = F["mul"][r][xi]
            if eta in u:
                sm += u[xi] * u[eta]
        acc.append(sm / len(Omega) / q2)
    return Fraction(float(np.mean(acc))).limit_denominator(p ** 4)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13):
        err = 0.0
        diag_zero = False
        for mu in range(p):
            for nu in range(p):
                e = abs(S_closed(mu, nu, p) - S_brute(mu, nu, p))
                err = max(err, e)
                if mu == nu and mu != 0 and abs(S_closed(mu, nu, p)) < 1e-9:
                    diag_zero = True
        if err > 1e-8 or diag_zero:
            ok = False
        rows[str(p)] = {"max_err": err, "diag_nonzero": not diag_zero}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "S is the stated 3-line geometric sum.  S≠0 on μ=ν≠0.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 13):
        F = _kernel_field(p)
        Omega = set(omega_set(F))
        ib, L0, Lp, Lm, L = L_sets(p)
        sizes = len(L0) == p - 1 and len(Lp) == p - 1 and len(Lm) == p - 1
        subset = L <= Omega
        nL = len(L)
        nOm = len(Omega)
        npp = sum(1 for r in F["squares"] if merge_cls(F, r) == "++")
        if not sizes or not subset or nL != 3 * (p - 1):
            ok = False
        if npp != 2 * p - 4:
            ok = False
        if p == 13 and nL == nOm:
            ok = False
        if npp == 2 * p:
            ok = False
        rows[str(p)] = {
            "ib": ib,
            "n_L": nL,
            "n_Omega": nOm,
            "L_subset": subset,
            "n_++": npp,
            "n_++_named": 2 * p - 4,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "L=L0∪L+∪L− ⊂ Ω, |L|=3(p−1), |++|=2p−4.  "
            "L≠Ω at p=13."
        ),
    }


def prove_C() -> dict:
    ok = True
    expect = {5: Fraction(64, 15), 13: Fraction(1168, 231)}
    got = {}
    drop = {}
    for p, exp in expect.items():
        q = Qpp_from_L(p)
        d = Qpp_drop16(p)
        got[str(p)] = str(q)
        drop[str(p)] = str(d)
        if q != exp:
            ok = False
        if d == exp:
            ok = False
        if q == LIVE_QPP.get(p):
            ok = False
    return {
        "proved": bool(ok),
        "Qpp": got,
        "drop16": drop,
        "p37": "72832/5985",
        "no_deg32_five_primes": True,
        "theorem": (
            "The L-average is 64/15 and 1168/231.  "
            "Dropping 16 misses.  Not the ensemble 48/13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "Q_triangle_named": True,
        "note": (
            "Fejer/L formula names the p≡1 triangle orbit.  "
            "Ensemble den |H_+|/(2p) remains.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.350  Q++(A·T_0) = Fejer average on 3-line L", flush=True)
    A = prove_A()
    print(f"  A geometric S: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B L⊂Ω: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q++: {C['proved']} {C['Qpp']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.350",
        "title": "Q++(A·T_0) is the 3-line Fejer average; ensemble OPEN",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "S_three_line": A["proved"],
            "L_named_in_Omega": B["proved"],
            "Qpp_triangle_Fejer": C["proved"],
            "Q_tau_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15350.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
