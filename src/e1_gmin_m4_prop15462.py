#!/usr/bin/env python3
"""
Prop 15.462 — No single mixed complementary occupancy
(R_aff, J_{2m} with m≥4, or one aff+J2) has a named
slot-count that, added to C(r,3), lands in [c_min,c_eq]
at p=13 or p=17.  The p=7 extras {1, 1/2, (p−2)/2, 3/2}
give either a miss or the half-integer 1225/2 (not a
value of c).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.461 flags.

============================================================================
Setup.  15.460: CRRR contributes C(r,3).  15.461: the p=7
remainder shapes do not extend.  Jacobi J_{2m} (distinct
poles, Σε=0, 2m<p) has k=m.  R_aff has k=(p²−1)/24.
aff+J2 k-set from 15.461.  Named extras from 15.455–15.458:
1 (flat), 1/2 (J6), (p−2)/2 (J4-pair), 3/2 (full support).

============================================================================
Theorem A — PROVED (full occupancy enum; p=13,17).
  Mixed slot-counts c1 at p=13 are
      {7,42,105,140,210,420,630,840,1260,2520,5040}.
  None lie in [516,582] = [c_min,c_eq]−C(7,3).
  C(r,3)+c1 nearest 455=35+420 and 665=35+630.
  At p=17, none of the c1 lie in [5700,6387]; nearest
  5124=84+5040 and 7644=84+7560.
  Fail: 665∈[551,617].  Fail: 5124∈[5784,6471].  ∎

Theorem B — PROVED (same enum; p=7 extras).
  extra=1/2 or 3/2: 0 hits in either window.
  extra=(p−2)/2 at p=13 on the 10 occupancies with
  c1=105 and n_{J4}≥2 gives 1225/2=612.5, which lies
  in the real interval [551,617] but is not an integer,
  hence not c=Σ(p+1)/s_τ (stabs divide p+1).
  Fail: claim 1225/2 names c.  ∎

Theorem C — PROVED (same extras; p=7).
  C(4,3)+c1 ∈ {8,10,16,28} under extra=1, never {18,19}.
  So no single-mix+CRRR rule of this dictionary is the
  live law.  Fail: 16=19.  ∎

Theorem D — OPEN.  Unnamed two-occupancy sums can land
in the window (census).  No named filter.  Q_τ unnamed.
phi_F not imported.

============================================================================
Backend: occupancy recursion + 15.461 aff+J2 k-sets.
Writes evidence/e1_gmin_m4_prop15462.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15461 import affJ2_k_set  # noqa: E402
from e1_gmin_m4_prop15460 import c_CRRR, k_R_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15462_catalog.json"


def r_of(p: int) -> int:
    return (p + 1) // 2


def jacobi_kinds(p: int) -> list[tuple[str, int]]:
    """C, J_{2m} for 2m<p, R_aff."""
    ks: list[tuple[str, int]] = [("C", 0)]
    for m in range(1, (p - 1) // 2 + 1):
        if 2 * m >= p:
            break
        ks.append((f"J{2 * m}", m))
    ks.append(("R", k_R_named(p)))
    return ks


def _multinom(counts: tuple[int, ...]) -> int:
    acc = factorial(sum(counts))
    for c in counts:
        acc //= factorial(c)
    return acc


def is_crrr(p: int, desc: dict) -> bool:
    r = r_of(p)
    return (
        desc.get("R") == 3
        and desc.get("C") == r - 3
        and set(desc) <= {"C", "R"}
    )


def is_mixed(desc: dict) -> bool:
    if desc.get("R", 0) or desc.get("A"):
        return True
    for key, n in desc.items():
        if key.startswith("J") and n and int(key[1:]) >= 8:
            return True
    return False


@lru_cache(maxsize=8)
def _aff_ks(p: int) -> tuple[int, ...]:
    return tuple(sorted(affJ2_k_set(p)))


@lru_cache(maxsize=8)
def mixed_rows(p: int) -> tuple[dict, ...]:
    """Mixed occupancies other than CRRR, with slot-count c1."""
    r = r_of(p)
    target = comb(r, 2)
    K = jacobi_kinds(p)
    names = [k[0] for k in K]
    base = [k[1] for k in K]
    raw: list[tuple[dict, int]] = []

    def rec(i: int, left_slots: int, left_k: int, acc: list[int]) -> None:
        if i == len(base) - 1:
            n = left_slots
            if n * base[i] == left_k:
                acc2 = acc + [n]
                desc = {names[j]: acc2[j] for j in range(len(names)) if acc2[j]}
                raw.append((desc, _multinom(tuple(acc2))))
            return
        for n in range(left_slots + 1):
            rec(i + 1, left_slots - n, left_k - n * base[i], acc + [n])

    rec(0, r, target, [])
    for ka in _aff_ks(p):
        def rec2(i: int, left_slots: int, left_k: int, acc: list[int], ka: int = ka) -> None:
            if i == len(base) - 1:
                n = left_slots
                if n * base[i] == left_k:
                    acc2 = acc + [n]
                    desc = {"A": 1}
                    desc.update(
                        {names[j]: acc2[j] for j in range(len(names)) if acc2[j]}
                    )
                    raw.append((desc, r * _multinom(tuple(acc2))))
                return
            for n in range(left_slots + 1):
                rec2(i + 1, left_slots - n, left_k - n * base[i], acc + [n])

        rec2(0, r - 1, target - ka, [])

    out = []
    for desc, c1 in raw:
        if is_mixed(desc) and not is_crrr(p, desc):
            out.append({"occ": desc, "c1": c1})
    return tuple(out)


def extra_named(p: int, desc: dict, kind: str) -> Fraction:
    if kind == "flat":
        return Fraction(1)
    if kind == "half":
        return Fraction(1, 2) if desc.get("J6", 0) else Fraction(1)
    if kind == "p2":
        return Fraction(p - 2, 2) if desc.get("J4", 0) >= 2 else Fraction(1)
    if kind == "three_halves":
        npos = sum(v for k, v in desc.items() if k != "C" and v)
        # full complementary support: no C slots
        if desc.get("C", 0) == 0:
            return Fraction(3, 2)
        return Fraction(1)
    raise ValueError(kind)


def singles_in_window(p: int, kind: str) -> list[dict]:
    lo, hi, _ = c_eq_ends(p)
    cr = c_CRRR(p)
    hits = []
    for rw in mixed_rows(p):
        e = extra_named(p, rw["occ"], kind)
        c = cr + e * rw["c1"]
        if lo <= c <= hi:
            hits.append({"occ": rw["occ"], "c": str(c), "c1": rw["c1"], "e": str(e)})
    return hits


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, n_lo, n_hi, below, above in (
        (13, 516, 582, 455, 665),
        (17, 5700, 6387, 5124, 7644),
    ):
        cr = c_CRRR(p)
        lo, hi, _ = c_eq_ends(p)
        c1s = sorted({rw["c1"] for rw in mixed_rows(p)})
        gap = [c1 for c1 in c1s if n_lo <= c1 <= n_hi]
        vals = sorted({cr + c1 for c1 in c1s})
        bel = [v for v in vals if v < lo]
        abv = [v for v in vals if v > hi]
        nb = bel[-1] if bel else None
        na = abv[0] if abv else None
        rows[str(p)] = {
            "c1s": c1s,
            "gap": gap,
            "nearest_below": nb,
            "nearest_above": na,
            "window": [lo, hi],
            "c_CRRR": cr,
        }
        if gap:
            ok = False
        if nb != below or na != above:
            ok = False
        if lo <= above <= hi or lo <= below <= hi:
            ok = False
    if 551 <= 665 <= 617:
        ok = False
    if 5784 <= 5124 <= 6471:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "No mixed c1 in the needed gap. Nearest 455,665 at p=13 "
            "and 5124,7644 at p=17. Fail: 665∈window."
        ),
    }


def prove_B() -> dict:
    ok = True
    rec = {}
    for p in (13, 17):
        rec[str(p)] = {}
        for kind in ("flat", "half", "p2", "three_halves"):
            hits = singles_in_window(p, kind)
            rec[str(p)][kind] = {
                "n": len(hits),
                "cs": sorted({h["c"] for h in hits}),
            }
            if kind == "p2" and p == 13:
                if not hits:
                    ok = False
                if any(h["c"] != "1225/2" for h in hits):
                    ok = False
                if Fraction(1225, 2).denominator == 1:
                    ok = False
                if not (551 <= Fraction(1225, 2) <= 617):
                    ok = False
            elif hits and not (kind == "p2" and p == 13):
                # any integer hit is a fail of "no named mix"
                if any(Fraction(h["c"]).denominator == 1 for h in hits):
                    ok = False
    if Fraction(1225, 2) == 612:
        ok = False
    return {
        "proved": bool(ok),
        "rec": rec,
        "theorem": (
            "Named extras: 0 integer window hits. p=13 (p-2)/2 on "
            "c1=105 is 1225/2, not an integer. Fail: 1225/2 names c."
        ),
    }


def prove_C() -> dict:
    p = 7
    cr = c_CRRR(p)
    c1s = sorted({rw["c1"] for rw in mixed_rows(p)})
    vals = sorted({cr + c1 for c1 in c1s})
    lo, hi, _ = c_eq_ends(p)
    in_win = [v for v in vals if lo <= v <= hi]
    ok = in_win == [] and 16 in vals and 19 not in vals
    if 16 == 19:
        ok = False
    hits_p7 = {kind: singles_in_window(p, kind) for kind in ("flat", "half", "p2", "three_halves")}
    if any(hits_p7[k] for k in hits_p7):
        ok = False
    return {
        "proved": bool(ok),
        "c1s": c1s,
        "vals": vals,
        "window": [lo, hi],
        "hits": {k: len(v) for k, v in hits_p7.items()},
        "theorem": (
            "C(4,3)+c1 never in {18,19}. Nearest 16. Fail: 16=19."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Single mixed + CRRR is dead as a name. Unnamed "
            "two-occupancy sums can land; no filter. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.462  no single mixed + C(r,3) lands in the window", flush=True)
    A = prove_A()
    print(
        f"  A gap: {A['proved']} 13={A['rows']['13']['nearest_below']},"
        f"{A['rows']['13']['nearest_above']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B extras: {B['proved']} p13={B['rec']['13']}", flush=True)
    C = prove_C()
    print(f"  C p7: {C['proved']} vals={C['vals']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "near13": [
                    A["rows"]["13"]["nearest_below"],
                    A["rows"]["13"]["nearest_above"],
                ],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.462",
        "title": (
            "No single mixed occupancy + C(r,3) lands in "
            "[c_min,c_eq] at p=13 or p=17"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "c1_gap": A["proved"],
            "extras_no_integer_hit": B["proved"],
            "not_live_at_p7": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15462.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
