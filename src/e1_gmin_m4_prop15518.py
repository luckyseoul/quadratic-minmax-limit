#!/usr/bin/env python3
"""
Prop 15.518 — p=11 L=(0,1) affine halfspaces have four G-orbits;
the leftover (not AP, not AGL·QR0) includes a named stab=2p class.

  p=7: AP ∪ AGL·QR0 = all C(7,4) (15.309 B).
  p=11: C(11,6)=462 splits AP:55 (stab 4p), QR0:22 (stab p(p−1)),
  other:55 (stab 4p, not the AP orbit) + 330 (stab 2p).
  Fail: claim leftover S are a single stab; claim non-AP stab=2p at
  p=7 (live p(p−1)=42≠14).  Fail: these four orbits exhaust n_1d
  (1452≠2772).  Does not name NL / |H_+| / Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.517 flags.
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15139 import affine_halfspace
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15308 import G_order, frob_table, is_AP
from e1_gmin_m4_prop15309 import QR0_stab, qr0

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15518.json"

W = 40
_GENS: dict[int, list] = {}
_FOLD: dict[int, dict] = {}


def AP_stab(p: int) -> int:
    return 4 * p


def leftover_stab_named(p: int) -> int:
    """Named third-class stab at p=11. Fail at p=7."""
    return 2 * p


def _gens(p: int) -> list:
    if p in _GENS:
        return _GENS[p]
    F = _kernel_field(p)
    q = F["q"]
    n = q + 1
    mul, add = F["mul"], F["add"]
    fr = frob_table(F)

    def apply(a: int, b: int, k: int) -> np.ndarray:
        perm = np.empty(n, dtype=np.int32)
        perm[0] = 0
        for u in range(q):
            w = fr[u] if k else u
            w = add[mul[a][w]][b]
            perm[1 + u] = 1 + w
        return perm

    g2 = int(mul[F["g"]][F["g"]])
    gens = [apply(1, 1, 0), apply(1, p, 0), apply(g2, 0, 0), apply(1, 0, 1)]
    _GENS[p] = gens
    return gens


def orb_stab(p: int, y0: np.ndarray) -> tuple[int, int]:
    gens = _gens(p)
    if y0[0] < 0:
        y0 = -y0
    seen = {y0.tobytes()}
    stack = [y0]
    while stack:
        y = stack.pop()
        for g in gens:
            yn = np.empty_like(y)
            yn[g] = y
            if yn[0] < 0:
                yn = -yn
            key = yn.tobytes()
            if key not in seen:
                seen.add(key)
                stack.append(yn)
    orb = len(seen)
    return orb, G_order(p) // orb


def tag_S(p: int, S) -> str:
    S = frozenset(S)
    if is_AP(S, p):
        return "AP"
    Q = qr0(p)
    for a in range(1, p):
        for b in range(p):
            if frozenset((a * x + b) % p for x in Q) == S:
                return "QR0"
    return "other"


def _job(spec):
    p, S = spec
    y = np.sign(affine_halfspace(p, (0, 1), set(S))).astype(np.int8)
    orb, stab = orb_stab(p, y)
    return {"tag": tag_S(p, S), "orbit": orb, "stab": stab}


def fold(p: int) -> dict:
    if p in _FOLD:
        return _FOLD[p]
    m = (p + 1) // 2
    jobs = [(p, S) for S in combinations(range(p), m)]
    ws = min(W, len(jobs))
    with ProcessPoolExecutor(max_workers=ws) as ex:
        rows = list(ex.map(_job, jobs, chunksize=max(1, len(jobs) // ws)))
    by = defaultdict(list)
    for r in rows:
        by[r["tag"]].append(r)
    summary = {}
    for tag, rs in by.items():
        summary[tag] = {
            "nS": len(rs),
            "orbits": sorted(set(r["orbit"] for r in rs)),
            "stabs": sorted(set(r["stab"] for r in rs)),
            "orbit_hist": {str(k): v for k, v in sorted(Counter(r["orbit"] for r in rs).items())},
            "stab_hist": {str(k): v for k, v in sorted(Counter(r["stab"] for r in rs).items())},
        }
    rec = {"nS": len(rows), "summary": summary, "W": ws}
    _FOLD[p] = rec
    return rec


def _y(p: int, S) -> np.ndarray:
    return np.sign(affine_halfspace(p, (0, 1), set(S))).astype(np.int8)


def _in_orbit(p: int, y0: np.ndarray, y1: np.ndarray) -> bool:
    gens = _gens(p)
    if y0[0] < 0:
        y0 = -y0
    seen = {y0.tobytes()}
    stack = [y0]
    target = {y1.tobytes(), (-y1).tobytes()}
    if seen & target:
        return True
    while stack:
        y = stack.pop()
        for g in gens:
            yn = np.empty_like(y)
            yn[g] = y
            if yn[0] < 0:
                yn = -yn
            key = yn.tobytes()
            if key in target:
                return True
            if key not in seen:
                seen.add(key)
                stack.append(yn)
    return False


def prove_A() -> dict:
    rec = fold(7)
    s = rec["summary"]
    ok = set(s) == {"AP", "QR0"}
    ok = ok and s["AP"]["stabs"] == [AP_stab(7)]
    ok = ok and s["QR0"]["stabs"] == [QR0_stab(7)]
    if leftover_stab_named(7) in s.get("QR0", {}).get("stabs", []):
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "theorem": "p=7 L=(0,1) S are AP (stab 4p) or QR0 (stab p(p−1)). Not 2p.",
    }


def prove_B() -> dict:
    rec = fold(11)
    s = rec["summary"]
    other = s.get("other", {})
    ok = set(s) == {"AP", "QR0", "other"}
    ok = ok and s["AP"]["stabs"] == [AP_stab(11)]
    ok = ok and s["QR0"]["stabs"] == [QR0_stab(11)]
    ok = ok and leftover_stab_named(11) in other.get("stabs", [])
    ok = ok and AP_stab(11) in other.get("stabs", [])
    ok = ok and len(other.get("stabs", [])) == 2
    # fail-when-wrong: single leftover stab
    if other.get("stabs") == [leftover_stab_named(11)]:
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "theorem": "p=11 leftover S have stabs {2p,4p}; 2p is new. Not a single stab.",
    }


def prove_C() -> dict:
    rec = fold(11)
    y_ap = _y(11, (0, 1, 2, 3, 4, 5))
    y_o330 = _y(11, (0, 1, 2, 3, 5, 9))
    y_o660 = _y(11, (0, 1, 2, 3, 4, 6))
    y_qr = _y(11, (0, 1, 2, 4, 5, 7))
    distinct = (
        not _in_orbit(11, y_ap, y_o330)
        and not _in_orbit(11, y_ap, y_o660)
        and not _in_orbit(11, y_qr, y_o330)
        and not _in_orbit(11, y_qr, y_o660)
        and not _in_orbit(11, y_o330, y_o660)
    )
    affine_sum = 330 + 132 + 330 + 660
    ok = distinct and affine_sum != n_1d(11) and n_1d(11) == 2772
    return {
        "proved": bool(ok),
        "distinct_reps": bool(distinct),
        "affine_sum_if_four": affine_sum,
        "n_1d": n_1d(11),
        "fold": rec,
        "theorem": "Four distinct L=(0,1) orbits; sum 1452≠n_1d=2772.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Affine leftover stab 2p is named at p=11 only as a cache of "
            "L=(0,1) S. NL / |H_+| / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.518 p=11 affine leftover stab=2p; p=7 has no 2p class", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "prop": "15.518",
        "title": "p=11 affine leftover includes stab=2p; p=7 dichotomy has no 2p",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
