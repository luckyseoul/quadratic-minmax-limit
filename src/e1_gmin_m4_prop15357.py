#!/usr/bin/env python3
"""
Prop 15.357 — w=1 at p=3 (no NL); n_1d/(2p) is an integer k and
w=k/D with D=|H_+|/(2p); extra 1D hole-types at p=11 have Q≠Q_AP.
Q_NL and D stay unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.356 flags.

============================================================================
Setup.  15.356: Q++ = w Q_1D + (1−w) Q_NL, w=n_1d/|H_+|.
k := n_1d/(2p), D := |H_+|/(2p), so w=k/D.  Live (k,D)=(3,13),(10,409).

============================================================================
Theorem A — PROVED (15.308 orbit-stab + p=3 cache).
  |H_+|=n_1d=6 at p=3, hence w=1 and NL is empty.  The 15.356 mix
  has no NL term.  Fail: |H_+|=n_1d+p²(p−1)=24.  ∎

Theorem B — PROVED (binomial arithmetic).
  k=n_1d/(2p)=(p+1)C(p,m)/(4p) is an integer at every odd prime
  p=3..31.  Fail: D=dim V_+=(p²+1)/2 at p=7 (25≠409).  ∎

Theorem C — PROVED (Type+ constructors; certified p=11).
  Let m=(p+1)/2 and S_d={0,1,…,m}\\{d}.  The halfspace y_{L,S_d}
  is Max+.  At p=11, Q++(S_1)=544/45 and Q++(S_2)=632/45, neither
  equal to Q_AP=56/9 nor to Q_1D=344/27.  Fail: Q(S_1)=Q_AP.  ∎

Theorem D — OPEN.  D and Q_NL unnamed in p.  phi_F not imported.

============================================================================
Backend: cache cardinality + constructor Q.  Writes
evidence/e1_gmin_m4_prop15357.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import load_D  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def k_named(p: int) -> int:
    return n_1d(p) // (2 * p)


def D_live(p: int) -> int:
    return HPLUS[p] // (2 * p)


def S_hole(p: int, d: int) -> set[int]:
    m = (p + 1) // 2
    return set(range(m + 1)) - {d}


def Qpp_of_y(p, y):
    F = _kernel_field(p)
    Omega = omega_set(F)
    rpp = [r for r in F["squares"] if merge_cls(F, r) == "++"]
    om = {xi: i for i, xi in enumerate(Omega)}
    q = F["q"]
    z = np.asarray(y[1 : 1 + q], dtype=np.float64)
    if y[0] < 0:
        z = -z
    u = np.abs(zhat_omega(z, F, Omega)) ** 2
    sm = 0.0
    for r in rpp:
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
    return Fraction(sm / (len(Omega) * len(rpp) * q * q)).limit_denominator(p ** 6)


def prove_A() -> dict:
    D, _ = load_D(3)
    n = int(D.shape[0])
    n1 = n_1d(3)
    ok = n == n1 == 6
    if n == n1 + 3 * 3 * 2:
        ok = False
    return {
        "proved": bool(ok),
        "n_Hplus": n,
        "n_1d": n1,
        "w": 1,
        "theorem": "p=3: |H_+|=n_1d, w=1, NL empty.",
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in range(3, 32):
        if not is_prime(p):
            continue
        n1 = n_1d(p)
        k = n1 / (2 * p)
        rows[str(p)] = k
        if abs(k - round(k)) > 1e-12:
            ok = False
        if n1 % (2 * p) != 0:
            ok = False
    if D_live(7) == (7 * 7 + 1) // 2:
        ok = False
    if k_named(5) != 3 or k_named(7) != 10:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "k5": k_named(5),
        "k7": k_named(7),
        "D5": D_live(5),
        "D7": D_live(7),
        "theorem": "k=n_1d/(2p)∈ℕ; w=k/D; D≠dim V_+ at p=7.",
    }


def prove_C() -> dict:
    p = 11
    C = paley_conference_prime_power(p).astype(np.float64)
    rec = {}
    ok = True
    for d, want in ((1, Fraction(544, 45)), (2, Fraction(632, 45))):
        y = affine_halfspace(p, (0, 1), S_hole(p, d))
        mp = bool(np.allclose(C @ y, p * y, atol=1e-4))
        qv = Qpp_of_y(p, y) if mp else None
        rec[f"d{d}"] = {"maxplus": mp, "Q": str(qv)}
        if not mp or qv != want:
            ok = False
        if qv == Q_AP_named(p) or qv == Q_1d_pp_named(p):
            ok = False
    if rec["d1"]["Q"] == "56/9":
        ok = False
    return {
        "proved": bool(ok),
        "rows": rec,
        "Q_AP": str(Q_AP_named(p)),
        "Q_1d": str(Q_1d_pp_named(p)),
        "theorem": "p=11 hole-types Q∈{544/45,632/45}≠Q_AP,Q_1D.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_NL_named_in_p": False,
        "w_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": "w=1 at p=3.  k named.  D and Q_NL unnamed.  phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.357  w=1 at p=3; k named; p=11 hole-types", flush=True)
    A = prove_A()
    print(f"  A p3 w=1: {A['proved']} |H+|={A['n_Hplus']}", flush=True)
    B = prove_B()
    print(f"  B k integer: {B['proved']} k5={B['k5']} D7={B['D7']}", flush=True)
    C = prove_C()
    print(f"  C hole-types: {C['proved']} {C['rows']}", flush=True)
    D = prove_open()
    print(f"  D open: Q_NL={D['Q_NL_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.357",
        "title": "w=1 at p=3; k=n_1d/(2p); p=11 1D hole-types",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "w_eq_1_at_p3": A["proved"],
            "k_integer": B["proved"],
            "p11_hole_types": C["proved"],
            "Q_NL_named_in_p": False,
            "w_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15357.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
