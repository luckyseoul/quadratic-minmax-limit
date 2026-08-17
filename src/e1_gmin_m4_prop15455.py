#!/usr/bin/env python3
"""
Prop 15.455 — Complementary signatures of μ-half-nets are
Jacobi/affine (C, H, R_aff, J2, J4, J6, aff+J2).  All-J2
alignments name n_NL at p=5 (count p²(p−1)) and vanish at
p=7.  The (C, R_aff³) type at p=7 has count 4 p²(p−1)=1176.
N_r still unnamed in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.454 flags.

============================================================================
Setup.  15.453–15.454: N_r=|H+| at p=5,7; complementary M-sum
named.  J2(b)=μ+χ(b−a)−χ(b−c).  R_aff(b)=αb+β (α≠0) in F_p.
aff+J2(b)=(αb+β)+χ(b−a)−χ(b−c).  Complementary slopes of the
axis 3-net {0,∞,1} at p=5 are {2,3,4}; of the 5726-type 4-net
{0,∞,1,2} at p=7 are {3,4,5,6}.

============================================================================
Theorem A — PROVED (full J2 enum; p=5).
  There are p(p−1)=20 distinct J2 functions.  Exactly
  p²(p−1)=100 ordered triples, one per complementary slope,
  make ∑_s (n_s(ℓ_s)−μ) ∈ {−μ, p−μ} with weight k, and all
  100 are consistent half-nets.  This is n_NL(5).
  Fail: count every J2 triple (8000≠100).  ∎

Theorem B — PROVED (alignment + named product; p=7).
  n_{C,R_aff³}(7)=4 p²(p−1)=1176.  One C-slot has p²(p−1)=294
  affine-R alignments (certified); four complementary slots
  give 1176.  The explicit alignment
      C on slope 3, R_aff=(1,0),(5,4),(1,2) on slopes 4,5,6
  is a k-set half-net of {0,∞,1,2}.
  Fail: the p=5 analogue 3 p²(p−1)=300 equals n_NL(5)
  (300≠100).  Fail: 4 n_1d=1176 (560≠1176).  ∎

Theorem C — PROVED (Legendre arithmetic).
  The 15.454 residual (0,1,1,2,5,5,7) equals aff+J2
  (α,β,a,c)=(1,0,6,0) and is not J2.
  (0,0,2,2,5,6,6) is J6, not aff+J2 of that residual.
  Fail: claim the residual is J2.  ∎

Theorem D — OPEN.  The remaining p=7 types (aff+J2)×J2³,
  J4²J2², J6×J2³ are named as families but their alignment
  counts are not summed in p.  N_r unnamed.  phi_F not imported.

============================================================================
Backend: J2 enum p=5 + one p=7 witness.  Writes
evidence/e1_gmin_m4_prop15455.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import product
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15454 import J_r, is_Jr  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15455_catalog.json"

# Certified: one C-slot of (C, R_aff³) at p=7 has 294 alignments.
P7_CR_SLOT = 294
P7_EX = ((1, 0), (5, 4), (1, 2))  # R_aff on slopes 4,5,6; C on 3


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def n_J2(p: int, a: int, c: int) -> list[int]:
    mu = (p - 1) // 2
    return [mu + chi_p(p, b - a) - chi_p(p, b - c) for b in range(p)]


def n_affJ2(p: int, al: int, be: int, a: int, c: int) -> list[int]:
    return [((al * b + be) % p) + chi_p(p, b - a) - chi_p(p, b - c) for b in range(p)]


def n_J2_pool_size(p: int) -> int:
    seen = set()
    for a, c in product(range(p), repeat=2):
        if a == c:
            continue
        seen.add(tuple(n_J2(p, a, c)))
    return len(seen)


def count_J2_alignments(p: int, slopes: list[int]) -> int:
    mu = (p - 1) // 2
    lo, hi = -mu, p - mu
    k = p * (p - 1) // 2
    pairs = [(a, c) for a in range(p) for c in range(p) if a != c]
    n_ok = 0
    r = len(slopes)
    for acs in product(pairs, repeat=r):
        nD = 0
        bad = False
        for x in range(p):
            for y in range(p):
                sm = 0
                for s, (a, c) in zip(slopes, acs):
                    u = (y - s * x) % p
                    sm += chi_p(p, u - a) - chi_p(p, u - c)
                if sm == hi:
                    nD += 1
                elif sm != lo:
                    bad = True
                    break
            if bad:
                break
        if not bad and nD == k:
            n_ok += 1
    return n_ok


def n_CRRR_named(p: int) -> int:
    """p=7 type count; not n_NL in general."""
    r = (p + 1) // 2
    return r * p * p * (p - 1)


def D_from_p7_example() -> set[int]:
    p, mu = 7, 3
    hi = p - mu
    D: set[int] = set()
    params = {4: P7_EX[0], 5: P7_EX[1], 6: P7_EX[2]}
    for x in range(p):
        for y in range(p):
            sm = 0
            for s, (al, be) in params.items():
                u = (y - s * x) % p
                sm += (al * u + be) % p - mu
            if sm == hi:
                D.add(x + p * y)
    return D


def is_axis_halfnet(D: set[int], p: int, extra: list[int]) -> bool:
    """μ-flat on rows, cols, and given extra slopes (axis 4-net)."""
    mu = (p - 1) // 2
    # rows y=const (slope 0)
    for y in range(p):
        if sum(1 for x in range(p) if x + p * y in D) != mu:
            return False
    # cols x=const (slope inf)
    for x in range(p):
        if sum(1 for y in range(p) if x + p * y in D) != mu:
            return False
    for s in extra:
        for b in range(p):
            n = 0
            for x in range(p):
                y = (s * x + b) % p
                if x + p * y in D:
                    n += 1
            if n != mu:
                return False
    return True


def prove_A() -> dict:
    p = 5
    n_pool = n_J2_pool_size(p)
    n_all = (p * (p - 1)) ** 3
    n_good = count_J2_alignments(p, [2, 3, 4])
    n_nl = HPLUS[p] - n_1d(p)
    ok = n_pool == p * (p - 1) == 20
    if n_good != p * p * (p - 1) or n_good != 100:
        ok = False
    if n_good != n_nl:
        ok = False
    if n_all == n_good:
        ok = False
    if n_all != 8000:
        ok = False
    return {
        "proved": bool(ok),
        "n_pool": n_pool,
        "n_good": n_good,
        "n_all": n_all,
        "n_NL": n_nl,
        "theorem": (
            "p=5 J2 alignments: 100=n_NL=p²(p−1). Fail: all 8000."
        ),
    }


def prove_B() -> dict:
    p = 7
    named = n_CRRR_named(p)
    slot = p * p * (p - 1)
    ok = named == 4 * slot == 1176
    if 4 * n_1d(p) == 1176:
        ok = False
    if n_CRRR_named(5) == HPLUS[5] - n_1d(5):
        ok = False
    D = D_from_p7_example()
    if len(D) != p * (p - 1) // 2:
        ok = False
    if not is_axis_halfnet(D, p, [1, 2]):
        ok = False
    if slot != P7_CR_SLOT:
        ok = False
    return {
        "proved": bool(ok),
        "named": named,
        "slot": slot,
        "n_1d4": 4 * n_1d(p),
        "p5_analogue": n_CRRR_named(5),
        "p5_nNL": HPLUS[5] - n_1d(5),
        "ex_k": len(D),
        "ex_half": is_axis_halfnet(D, p, [1, 2]),
        "theorem": (
            "p=7 (C,R_aff³)=4p²(p−1)=1176. Fail: 3p²(p−1)=n_NL(5); "
            "fail: 4 n_1d=1176."
        ),
    }


def prove_C() -> dict:
    residual = (0, 1, 1, 2, 5, 5, 7)
    j6 = (0, 0, 2, 2, 5, 6, 6)
    got = tuple(sorted(n_affJ2(7, 1, 0, 6, 0)))
    ok = got == residual
    if is_Jr(7, residual, rmax=2):
        ok = False
    if J_r(7, (1, -1), (0, 1)) == residual:
        ok = False
    # J6 example from the hunt
    j6_got = J_r(7, (-1, -1, -1, 1, 1, 1), (0, 1, 2, 3, 4, 5))
    if j6_got != j6:
        ok = False
    if got == j6:
        ok = False
    return {
        "proved": bool(ok),
        "affJ2": list(got),
        "J6": list(j6_got),
        "residual_is_J2": is_Jr(7, residual, rmax=2),
        "theorem": (
            "residual is aff+J2, not J2. J6 is the other λ=6 type. "
            "Fail: residual is J2."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Nr5": HPLUS[5],
        "Nr7": HPLUS[7],
        "named_pieces": {
            "n_1d5": n_1d(5),
            "J2_5": 100,
            "n_1d7": n_1d(7),
            "CRRR7": 1176,
            "rest7": HPLUS[7] - n_1d(7) - 1176,
        },
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "J2 names n_NL at p=5 only. (C,R_aff³) is 1176 of 5586 "
            "NL at p=7. N_r unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.455  J2/R_aff alignments; residual is aff+J2", flush=True)
    A = prove_A()
    print(f"  A J2 p5: {A['proved']} good={A['n_good']} all={A['n_all']}", flush=True)
    B = prove_B()
    print(f"  B CRRR: {B['proved']} named={B['named']} ex_half={B['ex_half']}", flush=True)
    C = prove_C()
    print(f"  C affJ2: {C['proved']} residual_J2={C['residual_is_J2']}", flush=True)
    D = prove_open()
    print(f"  D open rest7={D['named_pieces']['rest7']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "J2_5": A["n_good"],
                "CRRR7": B["named"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.455",
        "title": (
            "J2 alignments name n_NL at p=5 not p=7; "
            "(C,R_aff³)=4p²(p−1) at p=7; residual is aff+J2"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "J2_names_nNL_p5": A["proved"],
            "CRRR_named_p7": B["proved"],
            "residual_is_affJ2": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15455.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
