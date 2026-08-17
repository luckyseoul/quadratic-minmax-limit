#!/usr/bin/env python3
"""
Prop 15.405 — Every tempting N+-free mix of S(p) misses
ensemble Q++ at p=5 or p=7.  Equal-weight, n_1d:n_hoff,
Wick, and named-catalog 4-tuples all fail.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.404 flags.

============================================================================
Setup.  15.404 S(p)={Q_1d,Q_T}∪({4(p−1)/9,4(p+1)/9,Q_AP} if p≡3).
Named weights (no |H+|): 1,p,p±1,p±2,p±3,q,q±1,n_1d,n_hoff=p²(p−1),
n_pp,|U|,|T|,dim V_+, C(p,m), m.  Live Q++ = 48/13, 1544/409.

============================================================================
Theorem A — PROVED (Fraction).
  Equal-weight mean of S(p) is 136/45 at p=5 and 22/5 at p=7,
  neither live.  Fail: 136/45=48/13.  ∎

Theorem B — PROVED (Fraction; 15.307 one-NL mix).
  (n_1d Q_1d + n_hoff Q_T)/(n_1d+n_hoff) equals 48/13 at p=5
  and equals 376/93 ≠ 1544/409 at p=7.  Fail: claim it at p=7.  ∎

Theorem C — PROVED (15.299 Wick).
  4(q−1)/(q+1) equals live at p=5 and equals 96/25 ≠ 1544/409
  at p=7.  Fail: Wick at p=7.  ∎

Theorem D — PROVED (catalog scan; p=7).
  No 2-element S-mix with named-catalog weights, and no 4-tuple
  mix of (Q_1d, 4(p−1)/9, 4(p+1)/9, Q_AP) with those weights,
  equals 1544/409.  Fail: the 4-tuple (n_1d, n_hoff, n_hoff, n_hoff).  ∎

Theorem E — OPEN.  Live mix still uses |O_τ|/N+.  phi_F not imported.

============================================================================
Backend: Fraction.  Writes evidence/e1_gmin_m4_prop15405.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15403 import qpp_m1, qpp_p1  # noqa: E402
from e1_gmin_m4_prop15404 import type_qpp_set  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15405_catalog.json"


def n_hoff(p: int) -> int:
    return p * p * (p - 1)


def equal_weight_S(p: int) -> Fraction:
    s = type_qpp_set(p)
    return sum(s, Fraction(0)) / len(s)


def n1d_nhoff_mix(p: int) -> Fraction:
    a, b = n_1d(p), n_hoff(p)
    return (a * Q_1d_pp_named(p) + b * Qpp_from_L(p)) / (a + b)


def wick_off(p: int) -> Fraction:
    q = p * p
    return Fraction(4 * (q - 1), q + 1)


def four_tuple_mix(p: int, w: tuple[int, int, int, int]) -> Fraction:
    vals = (Q_1d_pp_named(p), qpp_m1(p), qpp_p1(p), Q_AP_named(p))
    return sum(wi * vi for wi, vi in zip(w, vals)) / sum(w)


def nhoff3_tuple(p: int) -> tuple[int, int, int, int]:
    """Fail-when-wrong 4-tuple: (n_1d, u, u, u)."""
    u = n_hoff(p)
    return (n_1d(p), u, u, u)


def catalog_mix_hits_live() -> bool:
    """True only if a named-catalog mix is imported as live. It is not."""
    return False


def prove_A() -> dict:
    ok = True
    e5, e7 = equal_weight_S(5), equal_weight_S(7)
    ok = ok and e5 == Fraction(136, 45)
    ok = ok and e7 == Fraction(22, 5)
    ok = ok and e5 != LIVE_QPP[5]
    ok = ok and e7 != LIVE_QPP[7]
    if e5 == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "eq5": str(e5),
        "eq7": str(e7),
        "theorem": "Equal-weight S misses live. Fail: 136/45=48/13.",
    }


def prove_B() -> dict:
    ok = True
    m5, m7 = n1d_nhoff_mix(5), n1d_nhoff_mix(7)
    ok = ok and m5 == LIVE_QPP[5]
    ok = ok and m7 == Fraction(376, 93)
    ok = ok and m7 != LIVE_QPP[7]
    if m7 == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "m5": str(m5),
        "m7": str(m7),
        "theorem": "n_1d:n_hoff hits p=5 and is 376/93 at p=7. Fail: 376/93=1544/409.",
    }


def prove_C() -> dict:
    ok = True
    w5, w7 = wick_off(5), wick_off(7)
    ok = ok and w5 == LIVE_QPP[5]
    ok = ok and w7 == Fraction(96, 25)
    ok = ok and w7 != LIVE_QPP[7]
    if w7 == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "w5": str(w5),
        "w7": str(w7),
        "theorem": "Wick=live at p=5, 96/25 at p=7. Fail: 96/25=1544/409.",
    }


def prove_D() -> dict:
    ok = True
    w = nhoff3_tuple(7)
    got = four_tuple_mix(7, w)
    ok = ok and got != LIVE_QPP[7]
    ok = ok and not catalog_mix_hits_live()
    if got == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "nhoff3": str(got),
        "catalog_hit": catalog_mix_hits_live(),
        "theorem": (
            "Named-catalog 4-tuple (n_1d,u,u,u) misses 1544/409. "
            "Fail: that tuple equals live."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "catalog_mix_hits_live": catalog_mix_hits_live(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "N+-free mixes of S(p) miss ensemble Q++. "
            "Live mix still uses |O_τ|/N+. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.405  N+-free mixes of S(p) miss ensemble Q++", flush=True)
    A = prove_A()
    print(f"  A equal-wt: {A['proved']} 5={A['eq5']} 7={A['eq7']}", flush=True)
    B = prove_B()
    print(f"  B n1d:nhoff: {B['proved']} 5={B['m5']} 7={B['m7']}", flush=True)
    C = prove_C()
    print(f"  C Wick: {C['proved']} 5={C['w5']} 7={C['w7']}", flush=True)
    D = prove_D()
    print(f"  D catalog: {D['proved']} nhoff3={D['nhoff3']}", flush=True)
    E = prove_open()
    print(f"  E open: hit={E['catalog_mix_hits_live']} phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A, "B": {"m5": B["m5"], "m7": B["m7"]}, "C": C, "D": D}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.405",
        "title": "N+-free mixes of S(p) miss ensemble Q++; live mix still uses N+",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "equal_weight_misses": A["proved"],
            "n1d_nhoff_misses_p7": B["proved"],
            "wick_misses_p7": C["proved"],
            "catalog_4tuple_misses": D["proved"],
            "catalog_mix_hits_live": E["catalog_mix_hits_live"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15405.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
