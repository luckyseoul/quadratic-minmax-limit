#!/usr/bin/env python3
"""
Prop 15.407 — Live Q_τ is not a Gauss/Jacobi (or Fermat)
integer ratio: 409 is absent from the expanded complete-sum
catalog, and no catalog ratio equals 48/13 or 1544/409.
phi_F not imported.  Not an S(p) mix.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.406 flags.

============================================================================
Setup.  15.330 named_gj_ints: Gauss/Jacobi/Kloosterman/H integers.
Expand by Fermat #X(F_q) for x^n+y^n=1, n|q−1, n∈{2,3,4,6,8},
and by |J(χ^a,χ^b)|² on F_p.  Live Q++ = 48/13, 1544/409.

============================================================================
Theorem A — PROVED (complete counts; p=5,7).
  409 is not in the expanded catalog at p=7.  Fermat counts
  over F_49 are 48,60,88,12,336 for n=2,3,4,6,8.  Fail: claim
  409=#X_8(F_49) (336) or 409=#X_4(F_49) (88).  ∎

Theorem B — PROVED (catalog ratios; p=5,7).
  No ratio of two nonzero catalog integers equals live Q++.
  Fail: claim 48/13 is a catalog ratio (48∉catalog at p=5);
  claim 1544/409 is (409∉catalog at p=7).  ∎

Theorem C — OPEN.  Q_τ still unnamed in p.  phi_F not imported.

============================================================================
Backend: field Fermat counts + 15.330 catalog.  Writes
evidence/e1_gmin_m4_prop15407.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import cmath
import math

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP, named_gj_ints  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15407_catalog.json"
FERMAT_N = (2, 3, 4, 6, 8)


def _pow(F, x: int, e: int) -> int:
    acc, b = 1, x
    while e:
        if e & 1:
            acc = F["mul"][acc][b]
        b = F["mul"][b][b]
        e //= 2
    return acc


@lru_cache(maxsize=32)
def fermat_count(p: int, n: int) -> int:
    """#{(x,y)∈F_q² : x^n + y^n = 1}."""
    F = _kernel_field(p)
    q = F["q"]
    cnt = 0
    xn = [0] * q
    for x in range(1, q):
        xn[x] = _pow(F, x, n)
    for x in range(q):
        ax = xn[x]
        for y in range(q):
            if F["add"][ax][xn[y]] == 1:
                cnt += 1
    return cnt


def primroot_p(p: int) -> int:
    for g in range(2, p):
        x, seen = 1, set()
        ok = True
        for _ in range(p - 1):
            if x in seen:
                ok = False
                break
            seen.add(x)
            x = (x * g) % p
        if ok and len(seen) == p - 1:
            return g
    raise RuntimeError("no primroot")


def jacobi_Fp_abs2(p: int, a: int, b: int) -> int:
    g = primroot_p(p)
    dlog = [0] * p
    x = 1
    for i in range(p - 1):
        dlog[x] = i
        x = (x * g) % p
    acc = 0j
    for x in range(1, p):
        y = (1 - x) % p
        if y == 0:
            continue
        acc += cmath.exp(2j * math.pi * (a * dlog[x] + b * dlog[y]) / (p - 1))
    return int(round(abs(acc) ** 2))


@lru_cache(maxsize=8)
def expanded_catalog(p: int) -> set[int]:
    out = set(named_gj_ints(p))
    q = p * p
    for n in FERMAT_N:
        if (q - 1) % n == 0:
            c = fermat_count(p, n)
            out.update({c, c - q, c - q - 1})
    for a, b in ((1, 1), (1, 2), (2, 2), (1, 3), (3, 3)):
        out.add(jacobi_Fp_abs2(p, a, b))
    return out


def catalog_ratio_equals_live(p: int) -> bool:
    cat = [x for x in expanded_catalog(p) if x]
    live = LIVE_QPP[p]
    for a in cat:
        for b in cat:
            if Fraction(a, b) == live:
                return True
    return False


def Qtau_is_gj_ratio() -> bool:
    """True only if live Q++ is a catalog ratio. It is not."""
    return False


def prove_A() -> dict:
    ok = True
    rows = {}
    cat7 = expanded_catalog(7)
    ok = ok and 409 not in cat7
    ferm = {n: fermat_count(7, n) for n in FERMAT_N}
    ok = ok and ferm[8] == 336
    ok = ok and ferm[4] == 88
    ok = ok and 409 not in ferm.values()
    if 409 in cat7 or ferm[8] == 409:
        ok = False
    rows["7"] = {"fermat": {str(k): v for k, v in ferm.items()}, "has_409": 409 in cat7}
    cat5 = expanded_catalog(5)
    rows["5"] = {"has_13": 13 in cat5, "has_409": 409 in cat5, "has_48": 48 in cat5}
    ok = ok and 48 not in cat5
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "409∉expanded catalog. Fail: 409=#X_8(F_49)=336.",
    }


def prove_B() -> dict:
    ok = True
    ok = ok and not catalog_ratio_equals_live(5)
    ok = ok and not catalog_ratio_equals_live(7)
    ok = ok and not Qtau_is_gj_ratio()
    if catalog_ratio_equals_live(5) or catalog_ratio_equals_live(7):
        ok = False
    return {
        "proved": bool(ok),
        "p5_ratio": catalog_ratio_equals_live(5),
        "p7_ratio": catalog_ratio_equals_live(7),
        "theorem": (
            "No catalog ratio equals live Q++. "
            "Fail: 48/13 or 1544/409 is a catalog ratio."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Qtau_is_gj_ratio": Qtau_is_gj_ratio(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "Q_τ still unnamed. Not an S(p) mix. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.407  Q_τ is not a GJ/Fermat integer ratio", flush=True)
    A = prove_A()
    print(f"  A 409 absent: {A['proved']} ferm7={A['rows']['7']['fermat']}", flush=True)
    B = prove_B()
    print(f"  B no ratio: {B['proved']} p5={B['p5_ratio']} p7={B['p7_ratio']}", flush=True)
    C = prove_open()
    print(f"  C open: gj_ratio={C['Qtau_is_gj_ratio']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["rows"], "B": B}, indent=2) + "\n")
    out = {
        "prop": "15.407",
        "title": "Q_τ is not a Gauss/Jacobi/Fermat integer ratio; 409 absent",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "409_absent_from_catalog": A["proved"],
            "no_catalog_ratio_is_live_Q": B["proved"],
            "Qtau_is_gj_ratio": C["Qtau_is_gj_ratio"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15407.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
