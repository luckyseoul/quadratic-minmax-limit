#!/usr/bin/env python3
"""
Prop 15.461 — The three p=7 remainder families of 15.456
do not extend by the same complementary shapes.
    J6 × J2^{r−1} and J4² × J2^{r−2} have ∑k = r+2,
    equal to C(r,2) iff r=4.
    aff+J2 × J2^{r−1} needs k_aff = C(r,2)−r+1 = 15
    at p=13, but the aff+J2 k-set is {4,…,12}.
C(r,3)+15 = 50 ∉ [551,617].  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.460 flags.

============================================================================
Setup.  15.456: n_NL(7)−4u = 15u splits as aff+J2×J2³,
J6×J2³, J4²J2².  15.460: CRRR = C(r,3) u extends, with
k_R=(p²−1)/24.  Jacobi J_{2m} (distinct poles, Σε=0, 2m<p)
has k=m.  aff+J2(b)=(αb+β)+χ(b−a)−χ(b−c).

============================================================================
Theorem A — PROVED (k-sum; all r).
  J6 (k=3) plus r−1 copies of J2 (k=1) sums to r+2.
  Two J4 plus r−2 copies of J2 also sum to r+2.
  r+2 = C(r,2) ⇔ r=4.  Fail: claim either shape is legal
  at r=7 (r+2=9≠21).  ∎

Theorem B — PROVED (full aff+J2 k-set; p=13).
  aff+J2 × J2^{r−1} needs k_aff = C(r,2)−(r−1).
  At p=13 this is 15.  The k-set of all aff+J2 signatures
  is {4,5,6,7,8,9,10,11,12}, so 15 is absent.
  Fail: claim 15∈k(aff+J2) at p=13.  ∎

Theorem C — PROVED (A+B + 15.460).
  Naive C(r,3)+15 equals 19 at p=7 and 50 at p=13,
  and 50∉[551,617].  Fail: claim 50 names c at p=13.  ∎

Theorem D — OPEN.  Mixed types (R_aff with J8+, or
high-k aff+J2 with R_aff) are not named.  Q_τ unnamed.
phi_F not imported.

============================================================================
Backend: integer k-sums + p=13 aff+J2 enum (24 336 sigs).
Writes evidence/e1_gmin_m4_prop15461.json
"""
from __future__ import annotations

import json
import os
import sys
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15457 import c_eq_named  # noqa: E402
from e1_gmin_m4_prop15460 import c_CRRR, k_R_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15461_catalog.json"


def r_of(p: int) -> int:
    return (p + 1) // 2


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def k_from_n(p: int, n: list[int]) -> int:
    mu = (p - 1) // 2
    M = sum(x * x for x in n)
    num = M - p * mu * mu
    if num % p:
        return -1
    lam = num // p
    if lam % 2:
        return -1
    return lam // 2


def k_J6_plus_J2(r: int) -> int:
    return 3 + (r - 1) * 1


def k_J4sq_plus_J2(r: int) -> int:
    return 2 + 2 + (r - 2) * 1


def k_aff_needed(r: int) -> int:
    return comb(r, 2) - (r - 1)


def affJ2_k_set(p: int) -> set[int]:
    out: set[int] = set()
    for al in range(1, p):
        for be in range(p):
            for a in range(p):
                for c in range(p):
                    if a == c:
                        continue
                    n = [
                        ((al * b + be) % p)
                        + chi_p(p, b - a)
                        - chi_p(p, b - c)
                        for b in range(p)
                    ]
                    out.add(k_from_n(p, n))
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for r in (3, 4, 5, 7, 9):
        s1 = k_J6_plus_J2(r)
        s2 = k_J4sq_plus_J2(r)
        tgt = comb(r, 2)
        rows[str(r)] = {"J6J2": s1, "J4J2": s2, "C_r_2": tgt}
        if s1 != r + 2 or s2 != r + 2:
            ok = False
        if r == 4 and s1 != tgt:
            ok = False
        if r != 4 and s1 == tgt:
            ok = False
    if k_J6_plus_J2(7) == 21:
        ok = False
    if k_J6_plus_J2(7) != 9:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "J6×J2^{r-1} and J4²×J2^{r-2} have sum k=r+2, "
            "legal iff r=4. Fail: legal at r=7 (9≠21)."
        ),
    }


def prove_B() -> dict:
    ks = affJ2_k_set(13)
    need = k_aff_needed(7)
    need7 = k_aff_needed(4)
    ks7 = affJ2_k_set(7)
    ok = need == 15 and 15 not in ks
    ok = ok and ks == set(range(4, 13))
    ok = ok and need7 == 3 and 3 in ks7
    if 15 in ks:
        ok = False
    return {
        "proved": bool(ok),
        "k_set13": sorted(ks),
        "need13": need,
        "k_set7": sorted(ks7),
        "need7": need7,
        "theorem": (
            "aff+J2×J2^6 needs k_aff=15∉{4..12} at p=13. "
            "Fail: 15 in the k-set."
        ),
    }


def prove_C() -> dict:
    lo, hi, nwin = c_eq_ends(13)
    naive7 = c_CRRR(7) + 15
    naive13 = c_CRRR(13) + 15
    ok = naive7 == 19 and naive13 == 50
    if lo <= naive13 <= hi:
        ok = False
    if naive13 == c_eq_named(13):
        ok = False
    if k_R_named(13) != 7:
        ok = False
    return {
        "proved": bool(ok),
        "naive7": naive7,
        "naive13": naive13,
        "window13": [lo, hi, nwin],
        "theorem": (
            "C(r,3)+15 is 19 at p=7 and 50∉[551,617] at p=13. "
            "Fail: 50 names c."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "p=7 remainder shapes are r=4-only. Mixed R_aff/J8+ "
            "types unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.461  p=7 remainder shapes do not extend", flush=True)
    A = prove_A()
    print(f"  A shapes: {A['proved']} r7={A['rows']['7']}", flush=True)
    B = prove_B()
    print(f"  B affk: {B['proved']} need={B['need13']} set={B['k_set13']}", flush=True)
    C = prove_C()
    print(f"  C naive: {C['proved']} 50 win={C['window13']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "need13": B["need13"],
                "naive13": C["naive13"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.461",
        "title": (
            "p=7 remainder families J6×J2^{r-1}, J4²×J2^{r-2}, "
            "aff+J2×J2^{r-1} do not extend; C(r,3)+15=50 misses window"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "J_shapes_r4_only": A["proved"],
            "affJ2_k15_absent": B["proved"],
            "naive_50_misses": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15461.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
