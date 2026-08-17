#!/usr/bin/env python3
"""
Prop 15.337 — Isolated 1176-orbit = centrally symmetric 21-sets
with k-vector (3,2,1⁵,0) on lines through the centre.  The four
special directions are a harmonic quadruple; (full, empty) is an
opposite pair.  Still 15.x.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.336 flags.  Floor stays OPEN.

============================================================================
Setup.  Isolated type as in 15.336.  After translating the unique
involution to (x,y)↦(−x,−y), D is a union of antipodal pairs plus
the origin.  On each of the p+1=8 lines through 0 there are
(p−1)/2=3 pairs; k(ℓ)=#pairs of D on ℓ.

============================================================================
Theorem A — PROVED (1176-fold; certified p=7 cache).
  Every isolated D is centrally symmetric about a unique centre
  c∈D, and the k-vector of lines through c is {3,2,1,1,1,1,1,0}
  (one full line, one empty, one double, five singles).  Fail:
  some isolated D has k-hist ≠ (0¹,1⁵,2¹,3¹).  ∎

Theorem B — PROVED (complete sweep of 244944 candidates).
  Conversely, among all origin-centred c.s. 21-sets with that
  k-vector, exactly 24 are Max+.  Fail-when-wrong: claim 0 hits
  or that a hit has type ≠ isolated.  Those 24 are the isolated
  D with centre 0.  24 × 49 field-translates = 1176.  ∎

Theorem C — PROVED (the 24 hits).
  The four special directions (full, empty, double, leftover
  single) are always the harmonic quadruple
      H = {slope 0, ∞, 3, −3=4}
  on P¹(F_7): cr(0,∞; 3,4)=−1.  (full, empty) is an opposite
  pair {0,∞} or {3,−3}.  Fail: a Max+ hit with full-slope 1.  ∎

Theorem D — OPEN.  Pair-selection on the five singles is a
  3-element scale orbit (F_p^×/{±1}) for each of the 8 opposite
  role assignments; not yet a one-line Gauss/Jacobi formula in
  general p.  Size+parity: a c.s. |D|=p(p−1)/2 through 0 exists
  only for p≡3 (mod 4).  Ensemble Q_τ unnamed.  phi_F not
  imported.  ∎

============================================================================
Backend: ProcessPool sweep of pair-choices; cache fold for A.
Writes evidence/e1_gmin_m4_prop15337.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, product
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15336 import ISOL_TUPLE, dtype_of  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from workers import default_workers  # noqa: E402

P = 7
HARMONIC = {("s", 0), ("s", 3), ("s", 4), ("inf", None)}
OPPOSITE = [({("s", 0), ("inf", None)}), ({("s", 3), ("s", 4)})]

_G: dict = {}


def ag_mask(Drow, p: int = P) -> np.ndarray:
    m = np.zeros((p, p), dtype=bool)
    for u in np.flatnonzero(Drow):
        m[int(u % p), int(u // p)] = True
    return m


def center_of(mask, p: int = P):
    xs, ys = np.where(mask)
    S = set(zip((int(x) for x in xs), (int(y) for y in ys)))
    for c in S:
        if all(((2 * c[0] - x) % p, (2 * c[1] - y) % p) in S for x, y in S):
            return c
    return None


def translate(mask, c, p: int = P) -> np.ndarray:
    out = np.zeros_like(mask)
    xs, ys = np.where(mask)
    for x, y in zip(xs, ys):
        out[(int(x) - c[0]) % p, (int(y) - c[1]) % p] = True
    return out


def k_vector(mask, p: int = P) -> dict:
    rec = {}
    for s in range(p):
        rec[("s", s)] = _k_line(mask, ("s", s), p)
    rec[("inf", None)] = _k_line(mask, ("inf", None), p)
    return rec


def _k_line(mask, dir_, p: int) -> int:
    if dir_[0] == "s":
        pts = [((t) % p, (dir_[1] * t) % p) for t in range(1, p)]
    else:
        pts = [(0, t) for t in range(1, p)]
    used = set()
    taken = 0
    for x, y in pts:
        if (x, y) in used:
            continue
        mx, my = (-x) % p, (-y) % p
        used.add((x, y))
        used.add((mx, my))
        if mask[x, y] and mask[mx, my]:
            taken += 1
    return taken


def k_hist(kv: dict) -> tuple:
    c = Counter(kv.values())
    return tuple(sorted(c.items()))


K_HIST_ISOL = ((0, 1), (1, 5), (2, 1), (3, 1))


def prove_A() -> dict:
    D, F = load_D(7)
    classes = square_classes(F)
    hists = []
    n = 0
    n_cs = 0
    for row in D:
        if dtype_of(row, classes, 7) != ISOL_TUPLE:
            continue
        n += 1
        m = ag_mask(row)
        c = center_of(m)
        if c is None:
            hists.append("no_center")
            continue
        n_cs += 1
        m0 = translate(m, c)
        hists.append(k_hist(k_vector(m0)))
    cnt = Counter(hists)
    ok = n == 1176 and n_cs == 1176 and cnt == {K_HIST_ISOL: 1176}
    return {
        "proved": bool(ok),
        "n": n,
        "n_cs": n_cs,
        "hists": {str(k): int(v) for k, v in cnt.items()},
        "theorem": "Every isolated D is c.s. with k-hist (0¹,1⁵,2¹,3¹).",
    }


# ----- sweep (Theorem B, C) -----

DIRS = [("s", s) for s in range(P)] + [("inf", None)]


def _pair_data(dir_, p=P):
    if dir_[0] == "s":
        pts = [(t % p, (dir_[1] * t) % p) for t in range(1, p)]
    else:
        pts = [(0, t) for t in range(1, p)]
    used, pairs = set(), []
    for x, y in pts:
        if (x, y) in used:
            continue
        mx, my = (-x) % p, (-y) % p
        used.add((x, y))
        used.add((mx, my))
        pairs.append((x + p * y, mx + p * my))
    return pairs


PAIR = {d: _pair_data(d) for d in DIRS}


def _init_sweep():
    _G["C"] = paley_conference_prime_power(P).astype(np.float64)


def _mask_from_spec(spec):
    mask = np.zeros(P * P, dtype=bool)
    mask[0] = True
    for d, idxs in spec:
        d = (d[0], d[1])
        for i in idxs:
            u, v = PAIR[d][i]
            mask[u] = True
            mask[v] = True
    return mask


def _sweep_work(specs):
    C = _G["C"]
    hits = []
    for spec in specs:
        y = np.ones(P * P + 1, dtype=np.float64)
        y[1:] = np.where(_mask_from_spec(spec), -1.0, 1.0)
        if np.max(np.abs(C @ y - P * y)) < 1e-5:
            hits.append(spec)
    return hits


def all_kvec_specs():
    dirs = list(DIRS)
    out = []
    for full in dirs:
        rest = [d for d in dirs if d != full]
        for empty in rest:
            rest2 = [d for d in rest if d != empty]
            for double in rest2:
                singles = [d for d in rest2 if d != double]
                for dch in combinations(range(3), 2):
                    for sch in product(range(3), repeat=5):
                        spec = [
                            (full, [0, 1, 2]),
                            (empty, []),
                            (double, list(dch)),
                        ]
                        for sdir, si in zip(singles, sch):
                            spec.append((sdir, [si]))
                        out.append(spec)
    return out


def run_sweep(W: int | None = None):
    if W is None:
        W = default_workers()
    specs = all_kvec_specs()
    shards = np.array_split(np.array(specs, dtype=object), min(W, len(specs)))
    shards = [list(s) for s in shards]
    hits = []
    with ProcessPoolExecutor(max_workers=min(W, len(shards)), initializer=_init_sweep) as ex:
        for part in ex.map(_sweep_work, shards):
            hits.extend(part)
    return hits


def _roles(spec):
    full = empty = double = None
    for d, idxs in spec:
        d = (d[0], d[1])
        k = len(idxs)
        if k == 3:
            full = d
        elif k == 0:
            empty = d
        elif k == 2:
            double = d
    return full, empty, double


def prove_B() -> dict:
    hits = run_sweep()
    D, F = load_D(7)
    classes = square_classes(F)
    isol = {
        tuple(row.astype(bool).tolist())
        for row in D
        if dtype_of(row, classes, 7) == ISOL_TUPLE
    }
    n_isol_type = 0
    for spec in hits:
        mask = _mask_from_spec(spec)
        # Drow in AG encoding
        if tuple(mask.tolist()) in isol:
            n_isol_type += 1
    ok = len(hits) == 24 and n_isol_type == 24
    ok = ok and len(hits) != 0
    return {
        "proved": bool(ok),
        "n_hits": len(hits),
        "n_isol_type": n_isol_type,
        "n_specs": 244944,
        "n_centered_times_q": 24 * 49,
        "theorem": "Exactly 24 origin-c.s. Max+ of this k-type; all isolated.",
        "_hits": hits,
    }


def prove_C(hits=None) -> dict:
    if hits is None:
        hits = run_sweep()
    ok = True
    n_full_out = 0
    n_opp = 0
    for spec in hits:
        full, empty, double = _roles(spec)
        if full not in HARMONIC or empty not in HARMONIC or double not in HARMONIC:
            ok = False
        if full == ("s", 1):
            n_full_out += 1
            ok = False
        pair = {full, empty}
        if pair == {("s", 0), ("inf", None)} or pair == {("s", 3), ("s", 4)}:
            n_opp += 1
        else:
            ok = False
    # cr(0,∞;3,4) = 3/4 = −1 in F_7
    cr = (3 * pow(4, 5, 7)) % 7  # 4^{5}=4^{-1} since 7-2=5
    ok = ok and cr == 6  # −1
    ok = ok and n_opp == len(hits) == 24
    return {
        "proved": bool(ok),
        "n_hits": len(hits),
        "n_opposite_FE": n_opp,
        "n_full_slope1": n_full_out,
        "cross_ratio_0inf_34": cr,
        "harmonic": [str(d) for d in sorted(HARMONIC, key=str)],
        "theorem": "Special dirs = harmonic {0,∞,3,−3}; (F,E) opposite.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "isolated_named_in_p": False,
        "cs_size_only_p_3mod4": True,
        "note": (
            "p=7 isolated orbit identified as the harmonic c.s. family. "
            "Pair-selection on singles is a 3-scale orbit, not a Gauss "
            "formula in general p.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.337  isolated = harmonic c.s. (3,2,1^5,0)", flush=True)
    A = prove_A()
    print(f"  A k-hist: {A['proved']} {A['hists']}", flush=True)
    B = prove_B()
    print(f"  B sweep 24: {B['proved']} hits={B['n_hits']}", flush=True)
    C = prove_C(B.get("_hits"))
    print(f"  C harmonic: {C['proved']} cr={C['cross_ratio_0inf_34']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    B_out = {k: v for k, v in B.items() if k != "_hits"}
    out = {
        "prop": "15.337",
        "title": "Isolated 1176-orbit is the harmonic (3,2,1^5,0) family",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "isol_khist_3210": A["proved"],
            "sweep_24_all_isol": B["proved"],
            "harmonic_opposite_FE": C["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B_out, "C": C, "D": D},
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15337.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
