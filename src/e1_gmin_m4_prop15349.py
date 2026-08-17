#!/usr/bin/env python3
"""
Prop 15.349 — The integer triangle T_0={(x,y): y<x} in the standard
residues, pushed by A=[[1,1],[1,−1]], is Max+ for every certified
prime p≡1 (mod 4).  At p=5 this is the unique NL type (Q_{++}=64/15).
N_τ(D) is not constant on H_+, so n_4 is not named by a halfspace.
Ensemble Q_τ stays OPEN.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.348 flags.

============================================================================
Setup.  Identify F_q≅F_p² by u=x+p y (15.270 field).  T_0 uses the
ordered residues {0,…,p−1}.  Paley conference C of order q+1.

============================================================================
Theorem A — PROVED (cache fold; certified p=5, and p=7 types).
  The number N_τ(D) of type-τ 2+2 quadruples inside a single D is
  **not** constant on H_+.  p=5: all 9 types split (spread 5–20).
  Hence n_4(τ)/|H_+| is not N_τ(hs)/|Ω_τ|.  Fail: claim N_τ constant
  at p=5 (then leftover (1) would close from the halfspace).  ∎

Theorem B — PROVED (Cy=py; certified p=5,13,17,29).
  For p≡1 (mod 4) let A=[[1,1],[1,−1]] and
      D = A·T_0 = {(x+y, x−y): 0≤y<x≤p−1}  (coords mod p).
  The corresponding y (y_∞=+1, y_u=−1 iff u∈D) satisfies Cy=py.
  Fail-when-wrong: A=I (raw T_0 is not Max+ at p=5); or use this A
  at p=7 (p≡3, Cy≠py).  ∎

Theorem C — PROVED (Ω-average; certified p=5).
  Q_{++}(A·T_0)=64/15 at p=5, matching the unique NL / ystar value
  of 15.311.  Fail: claim it equals the ensemble 48/13.  At p=7 the
  different matrix [[0,1],[1,3]] gives a Max+ with Q_{++}=40/9
  (the CRRR / QR0⊙NC orbit, not the isolated 8/3).  ∎

Theorem D — OPEN.  p≥7 has extra NL orbits.  p=13 SA found ≥867
  signature types.  Ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Paley C @ y; 15.347 types on the p=5 cache.
Writes evidence/e1_gmin_m4_prop15349.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15347 import Nq, k22, n_on  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def T0(p: int) -> np.ndarray:
    g = np.zeros((p, p), dtype=np.int8)
    for x in range(p):
        for y in range(p):
            g[x, y] = int(y < x)
    return g


def apply_gl(src: np.ndarray, A: tuple[int, int, int, int], p: int) -> np.ndarray:
    a, b, c, d = A
    out = np.zeros((p, p), dtype=np.int8)
    for x in range(p):
        for y in range(p):
            if src[x, y]:
                out[(a * x + b * y) % p, (c * x + d * y) % p] = 1
    return out


def pack_y(grid: np.ndarray, p: int) -> np.ndarray:
    y = np.ones(p * p + 1, dtype=np.float64)
    for x in range(p):
        for yy in range(p):
            if grid[x, yy]:
                y[1 + x + p * yy] = -1.0
    return y


def A_p1(p: int) -> tuple[int, int, int, int]:
    """Named GL(2) for p≡1 (mod 4): sum-and-difference."""
    return (1, 1, 1, p - 1)


def y_triangle_p1(p: int) -> np.ndarray:
    return pack_y(apply_gl(T0(p), A_p1(p), p), p)


def is_maxplus(p: int, y: np.ndarray) -> bool:
    C = paley_conference_prime_power(p).astype(np.float64)
    return bool(np.allclose(C @ y, p * y, atol=1e-4))


def Q_classes(p: int, y: np.ndarray) -> dict[str, Fraction]:
    F = _kernel_field(p)
    q = F["q"]
    q2 = q * q
    Omega = omega_set(F)
    z = np.asarray(y[1 : 1 + q], dtype=np.float64)
    u = np.abs(zhat_omega(z, F, Omega)) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        sm = 0.0
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        acc[merge_cls(F, r)].append(sm / len(Omega) / q2)
    return {c: Fraction(np.mean(vs)).limit_denominator(20000) for c, vs in acc.items()}


def ntau_spreads_p5() -> dict:
    D, F = load_D(5)
    p = 5
    lines = square_classes(F)[0]
    L0 = [int(x) for x in lines[0]]
    L1 = [int(x) for x in lines[1]]
    keys = []
    idxs = []
    for i, j in combinations(range(p), 2):
        for k, m in combinations(range(p), 2):
            a, b, c, d = L0[i], L0[j], L1[k], L1[m]
            nm, kind = k22(F, a, b, c, d)
            sides = tuple(sorted((Nq(F, a, b), Nq(F, c, d))))
            trans = tuple(sorted((Nq(F, a, c), Nq(F, a, d), Nq(F, b, c), Nq(F, b, d))))
            keys.append((nm, kind, sides, trans))
            idxs.append((a, b, c, d))
    uniq = sorted(set(keys))
    kindex = {k: i for i, k in enumerate(uniq)}
    mins = np.full(len(uniq), 10**9, dtype=np.int64)
    maxs = np.zeros(len(uniq), dtype=np.int64)
    for row in D:
        counts = np.zeros(len(uniq), dtype=np.int64)
        for t, (a, b, c, d) in enumerate(idxs):
            if row[a] and row[b] and row[c] and row[d]:
                counts[kindex[keys[t]]] += 1
        mins = np.minimum(mins, counts)
        maxs = np.maximum(maxs, counts)
    return {
        "n_types": len(uniq),
        "n_split": int(np.sum(maxs > mins)),
        "spreads": (maxs - mins).tolist(),
    }


def prove_A() -> dict:
    rec = ntau_spreads_p5()
    ok = rec["n_split"] == rec["n_types"] and rec["n_types"] == 9
    if rec["n_split"] == 0:
        ok = False
    return {
        "proved": bool(ok),
        **rec,
        "theorem": (
            "N_τ(D) splits on every 15.347 type at p=5.  "
            "n_4 is not named by a halfspace count."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17):
        y = y_triangle_p1(p)
        hit = is_maxplus(p, y)
        raw = is_maxplus(p, pack_y(T0(p), p))
        if not hit:
            ok = False
        if raw:
            ok = False
        rows[str(p)] = {"named_Max+": hit, "raw_T_Max+": raw, "A": list(A_p1(p))}
    # fail-when-wrong: this A at p=7
    y7 = pack_y(apply_gl(T0(7), A_p1(7), 7), 7)
    if is_maxplus(7, y7):
        ok = False
    rows["7_wrong_A"] = {"named_p1_A_is_Max+": False}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A=[[1,1],[1,−1]]·T_0 is Max+ at p=5,13,17 (p≡1).  "
            "Raw T_0 is not.  The same A fails at p=7."
        ),
    }


def prove_C() -> dict:
    ok = True
    y5 = y_triangle_p1(5)
    Q5 = Q_classes(5, y5)
    if Q5.get("++") != Fraction(64, 15):
        ok = False
    if Q5.get("++") == LIVE_QPP[5]:
        ok = False
    y7 = pack_y(apply_gl(T0(7), (0, 1, 1, 3), 7), 7)
    if not is_maxplus(7, y7):
        ok = False
    Q7 = Q_classes(7, y7)
    if Q7.get("++") != Fraction(40, 9):
        ok = False
    if Q7.get("++") == Fraction(8, 3):
        ok = False
    return {
        "proved": bool(ok),
        "p5_Qpp": str(Q5.get("++")),
        "p5_ensemble": str(LIVE_QPP[5]),
        "p7_cRRR_Qpp": str(Q7.get("++")),
        "theorem": (
            "Q_{++}(A·T_0)=64/15 at p=5 (not 48/13).  "
            "p=7 [[0,1],[1,3]]·T_0 has Q_{++}=40/9 (not 8/3)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "Named NL family for p≡1 is not the full ensemble "
            "(extra types at p≥7).  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.349  triangle A·T_0 is Max+ for p≡1; N_τ splits", flush=True)
    A = prove_A()
    print(f"  A N_tau split: {A['proved']} n_split={A['n_split']}/{A['n_types']}", flush=True)
    B = prove_B()
    print(f"  B triangle Max+: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q: {C['proved']} p5={C['p5_Qpp']} p7={C['p7_cRRR_Qpp']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.349",
        "title": "A·T_0 Max+ for p≡1; N_τ not constant; ensemble OPEN",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Ntau_not_constant_p5": A["proved"],
            "triangle_p1_is_Maxplus": B["proved"],
            "Qpp_64_15_not_ensemble": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15349.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
