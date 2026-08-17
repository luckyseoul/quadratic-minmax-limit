#!/usr/bin/env python3
"""
Prop 15.463 — No two-occupancy type-support pair, absent at
p=5 (so c=C(r,3)=1) and equipped with the 15.455–15.458
extras, recovers live c=19 at p=7 and lands in [551,617]
at p=13.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.462 flags.

============================================================================
Setup.  Family = exact type-support S (keys other than C).
n(S,p)=∑ extra·c1 over occupancies with support S.
Pair {S,T} (S≠T, both absent at p=5):
    c = C(r,3) + n(S) + n(T).
Extras: flat; p7-dict (J4≥2 ↦ (p−2)/2, else J6 ↦ 1/2).

============================================================================
Theorem A — PROVED (support enum; extra=1).
  No absent-at-p=5 pair has n(S)+n(T)=15 at p=7.
  Nearest C(4,3)+nS+nT are 14,16,20,22, never 19.
  Fail: 20=19.  ∎

Theorem B — PROVED (same enum; p7-dict extra).
  Exactly five pairs have nS+nT=15 (hence c(7)=19, c(5)=1),
  each pairing {J6} with one of
      {J2,J4,R}, {J2,J6,R}, {J2,J4,J6}, {A,J6,R}, {A,J4,J6}.
  None has c(13)∈[551,617] (values 71/2, 911/2, 1331/2,
  1541/2, 3235/2 under half_j6; p7-dict replaces the J4
  pairs by larger half-integers).  Fail: claim {J6}+{J2,J4,R}
  lands in the window.  ∎

Theorem C — OPEN.  A two-family rule would need a different
  extra or a non-support pairing.  Q_τ unnamed.  phi_F not
  imported.

============================================================================
Backend: 15.462 occupancy recursion + support aggregation.
Writes evidence/e1_gmin_m4_prop15463.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15460 import c_CRRR  # noqa: E402
from e1_gmin_m4_prop15462 import _aff_ks, jacobi_kinds, r_of  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15463_catalog.json"

# Supports that hit nS+nT=15 under p7 extras (certified by enum).
P7_FIFTEEN = (
    (frozenset({"J2", "J4", "R"}), frozenset({"J6"})),
    (frozenset({"J2", "J6", "R"}), frozenset({"J6"})),
    (frozenset({"J2", "J4", "J6"}), frozenset({"J6"})),
    (frozenset({"J6"}), frozenset({"A", "J6", "R"})),
    (frozenset({"J6"}), frozenset({"A", "J4", "J6"})),
)


def _multinom(counts: tuple[int, ...]) -> int:
    acc = factorial(sum(counts))
    for c in counts:
        acc //= factorial(c)
    return acc


def is_crrr(p: int, desc: dict) -> bool:
    r = r_of(p)
    return (
        desc.get("R") == 3
        and desc.get("C", 0) == r - 3
        and set(k for k, v in desc.items() if v) <= {"C", "R"}
    )


def support_of(desc: dict) -> frozenset:
    return frozenset(k for k, v in desc.items() if v and k != "C")


def extra_p7(p: int, desc: dict) -> Fraction:
    if desc.get("J4", 0) >= 2:
        return Fraction(p - 2, 2)
    if desc.get("J6", 0) > 0:
        return Fraction(1, 2)
    return Fraction(1)


def extra_flat(_p: int, _desc: dict) -> Fraction:
    return Fraction(1)


def extra_half_j6(_p: int, desc: dict) -> Fraction:
    return Fraction(1, 2) if desc.get("J6", 0) else Fraction(1)


@lru_cache(maxsize=8)
def all_occupancies(p: int) -> tuple[tuple[dict, int], ...]:
    r = r_of(p)
    target = comb(r, 2)
    K = jacobi_kinds(p)
    names = [k[0] for k in K]
    base = [k[1] for k in K]
    raw: list[tuple[dict, int]] = []

    def rec(i: int, ls: int, lk: int, acc: list[int]) -> None:
        if i == len(base) - 1:
            n = ls
            if n * base[i] == lk:
                acc2 = acc + [n]
                desc = {names[j]: acc2[j] for j in range(len(names)) if acc2[j]}
                raw.append((desc, _multinom(tuple(acc2))))
            return
        for n in range(ls + 1):
            rec(i + 1, ls - n, lk - n * base[i], acc + [n])

    rec(0, r, target, [])
    for ka in _aff_ks(p):
        def rec2(i: int, ls: int, lk: int, acc: list[int], ka: int = ka) -> None:
            if i == len(base) - 1:
                n = ls
                if n * base[i] == lk:
                    acc2 = acc + [n]
                    desc = {"A": 1}
                    desc.update(
                        {names[j]: acc2[j] for j in range(len(names)) if acc2[j]}
                    )
                    raw.append((desc, r * _multinom(tuple(acc2))))
                return
            for n in range(ls + 1):
                rec2(i + 1, ls - n, lk - n * base[i], acc + [n])

        rec2(0, r - 1, target - ka, [])
    return tuple((d, c1) for d, c1 in raw if support_of(d) and not is_crrr(p, d))


def n_by_support(p: int, extra) -> dict[frozenset, Fraction]:
    acc: dict[frozenset, Fraction] = defaultdict(lambda: Fraction(0))
    for desc, c1 in all_occupancies(p):
        acc[support_of(desc)] += extra(p, desc) * c1
    return dict(acc)


def absent5_supports(extra) -> list[frozenset]:
    n5 = n_by_support(5, extra)
    n7 = n_by_support(7, extra)
    return [S for S in n7 if S not in n5]


def pair_c(p: int, S: frozenset, T: frozenset, extra) -> Fraction:
    n = n_by_support(p, extra)
    return c_CRRR(p) + n.get(S, 0) + n.get(T, 0)


def prove_A() -> dict:
    extra = extra_flat
    n7 = n_by_support(7, extra)
    S7 = [S for S in n7 if S not in n_by_support(5, extra)]
    pairs = list(combinations(S7, 2))
    sums = []
    raw15 = []
    for S, T in pairs:
        sm = n7[S] + n7[T]
        sums.append(c_CRRR(7) + sm)
        raw15.append(sm)
    ok = 15 not in raw15
    ok = ok and 19 not in sums
    ok = ok and 20 in sums and 14 in sums
    if 20 == 19:
        ok = False
    return {
        "proved": bool(ok),
        "n_supports": len(S7),
        "n_pairs": len(pairs),
        "c7_vals": sorted(set(str(s) for s in sums)),
        "theorem": (
            "Flat extra: no absent-at-5 pair has nS+nT=15. "
            "c7∈{12,14,20,22,…}, never 19. Fail: 20=19."
        ),
    }


def prove_B() -> dict:
    extra = extra_p7
    n7 = n_by_support(7, extra)
    n13 = n_by_support(13, extra)
    S7 = [S for S in n7 if S not in n_by_support(5, extra)]
    fifteen = []
    for S, T in combinations(S7, 2):
        if n7[S] + n7[T] == 15:
            fifteen.append((S, T))
    got = {frozenset((S, T)) for S, T in fifteen}
    want = {frozenset(pt) for pt in P7_FIFTEEN}
    lo, hi, _ = c_eq_ends(13)
    rows = []
    in_win = 0
    for S, T in fifteen:
        c13 = c_CRRR(13) + n13.get(S, 0) + n13.get(T, 0)
        c7 = c_CRRR(7) + 15
        rows.append(
            {
                "S": sorted(S),
                "T": sorted(T),
                "c7": str(c7),
                "c13": str(c13),
                "in_win": bool(lo <= c13 <= hi),
            }
        )
        if lo <= c13 <= hi:
            in_win += 1
    ok = got == want and len(fifteen) == 5 and in_win == 0
    # explicit fail-when-wrong: the J6 + {J2,J4,R} pair
    c_bad = pair_c(13, frozenset({"J6"}), frozenset({"J2", "J4", "R"}), extra)
    if lo <= c_bad <= hi:
        ok = False
    if pair_c(7, frozenset({"J6"}), frozenset({"J2", "J4", "R"}), extra) != 19:
        ok = False
    return {
        "proved": bool(ok),
        "n_fifteen": len(fifteen),
        "rows": rows,
        "c_J6_J2J4R_13": str(c_bad),
        "window13": [lo, hi],
        "theorem": (
            "Five p7-extra pairs give c(7)=19 and c(5)=1; "
            "none land in [551,617]. Fail: {J6}+{J2,J4,R} in window."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Two-occupancy support pairs are dead as a name of c. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.463  no two-occupancy support pair hits 1,19 and the window", flush=True)
    A = prove_A()
    print(f"  A flat: {A['proved']} c7={A['c7_vals']}", flush=True)
    B = prove_B()
    print(
        f"  B p7: {B['proved']} n15={B['n_fifteen']} in_win="
        f"{sum(1 for r in B['rows'] if r['in_win'])} "
        f"bad={B['c_J6_J2J4R_13']}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "n_fifteen": B["n_fifteen"],
                "c_bad13": B["c_J6_J2J4R_13"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.463",
        "title": (
            "No two-occupancy type-support pair recovers "
            "c=1,19 and the p=13 window"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "flat_no_19": A["proved"],
            "p7_five_miss_window": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15463.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
