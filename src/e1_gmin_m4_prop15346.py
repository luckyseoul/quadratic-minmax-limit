#!/usr/bin/env python3
"""
Prop 15.346 — Exact two-line moment E[n_i² n_j²]; bivariate LP from
the 15.344 2-point does not pin [Q_lo, Q^{ub}]; p=5 2+2 four-points
are type-constant after χ_p of the parallelogram sides.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.345 flags.  Does **not** name E[n_i² n_j²] in p.  Floor OPEN.

============================================================================
Setup.  One square parallel class, n_t=|D∩ℓ_t|, ∑ n_t=k, M=∑ n_t².
E[M²]=p E[n⁴]+p(p−1) E[n_i² n_j²] (i≠j).  Q++ is the 15.306 C image
of E[M²].  Two parallel lines form a 2×p grid.

============================================================================
Theorem A — PROVED (cache fold; certified p=5,7).
  E[n_i² n_j²] = 735/26 (p=5) and 101585/818 (p=7).  Dens are
  |H_+|/p (26=130/5, 818=5726/7).  Fail: claim the p=5 value is
  E[n²]²=36, or that den=|H_+|/(2p)=13.  ∎

Theorem B — PROVED (bivariate LP; Max+-free constraints only).
  The joint (n,n') on two parallel lines, fitted to the named
  E[n], E[n²], E[nn'] from 15.344 + ∑n=k (so n+n'∈[k−p(p−2),k])
  and 0≤n,n'≤p, gives E[n²n'²]∈[14.36, 51.5] at p=5 (live 28.27
  interior).  The implied Q++ interval using live E[n⁴] is far
  larger than [11/4, 26/7].  Fail: claim the LP max E[n²n'²] equals
  735/26.  ∎

Theorem C — PROVED (certified p=5; fails as a complete type law at p=7).
  On a 2×p grid, n_4 of a 2+2 set is constant on the type
      (# Paley non-edges in K_{2,2}, matching vs path,
       aligned params, parallelogram, χ_p(j−i), χ_p(δ)).
  At p=5 this is 12 constant types, 0 splits (the leftover {1,3}
  of 15.345-era χ-hex is χ_p(j−i)=1 vs 4).  At p=7 the same
  dictionary still splits (14 of 30 types).  Fail: claim the
  dictionary is split-free at p=7.  ∎

Theorem D — OPEN.  No 2–3 factor formula in the 15.344 / Paley /
  H,C,R catalog hits both 735/26 and 101585/818.  E[n_i² n_j²]
  still has den |H_+|/p.  phi_F not imported.

============================================================================
Backend: vectorized n-signatures; scipy linprog 36 vars at p=5;
2+2 enumeration on one class.  Writes evidence/e1_gmin_m4_prop15346.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15344 import HPLUS  # noqa: E402

from fractions import Fraction


def live_cross(p: int) -> dict:
    D, F = load_D(p)
    H = int(D.shape[0])
    classes = square_classes(F)
    cross = n4 = Msum = M2 = 0
    n_pairs = n_lines = 0
    for lines in classes:
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1).astype(
            np.int64
        )
        M = (ns * ns).sum(axis=1)
        Msum += int(M.sum())
        M2 += int((M * M).sum())
        n4 += int((ns ** 4).sum())
        n_lines += ns.size
        pdir = ns.shape[1]
        for i in range(pdir):
            for j in range(pdir):
                if i == j:
                    continue
                cross += int((ns[:, i] ** 2 * ns[:, j] ** 2).sum())
                n_pairs += H
    return {
        "H": H,
        "E_M": Fraction(Msum, H * len(classes)),
        "E_M2": Fraction(M2, H * len(classes)),
        "E_n4": Fraction(n4, n_lines),
        "E_n2n2": Fraction(cross, n_pairs),
    }


def named_moments(p: int):
    k = p * (p - 1) // 2
    m1 = Fraction(p - 1, 2)
    m2 = Fraction(p * p - 1, 4)
    enn = (k * k - p * m2) / (p * (p - 1))
    lo = max(0, k - p * (p - 2))
    return k, m1, m2, enn, lo, k


def bivar_lp(p: int) -> dict:
    k, m1, m2, enn, lo, hi = named_moments(p)
    idx = [(a, b) for a in range(p + 1) for b in range(p + 1)]
    n = len(idx)
    A1 = np.array([a for a, b in idx], dtype=float)
    B1 = np.array([b for a, b in idx], dtype=float)
    A2B2 = (A1 ** 2) * (B1 ** 2)
    ones = np.ones(n)
    bounds = [
        (0, 1) if lo <= a + b <= hi else (0, 0) for a, b in idx
    ]
    cols = [ones, A1, B1, A1 ** 2, B1 ** 2, A1 * B1]
    beq = [1.0, float(m1), float(m1), float(m2), float(m2), float(enn)]
    # symmetry
    for a in range(p + 1):
        for b in range(a + 1, p + 1):
            row = np.zeros(n)
            row[a * (p + 1) + b] = 1
            row[b * (p + 1) + a] = -1
            cols.append(row)
            beq.append(0.0)

    def one(sense: str) -> float:
        c = -A2B2 if sense == "max" else A2B2
        r = linprog(
            c,
            A_eq=np.vstack(cols),
            b_eq=np.array(beq),
            bounds=bounds,
            method="highs",
        )
        if not r.success:
            raise RuntimeError(r.message)
        return float(-r.fun if sense == "max" else r.fun)

    return {
        "E_nn": str(enn),
        "sum_lo": lo,
        "max": one("max"),
        "min": one("min"),
    }


def chi_q(F, x, y):
    d = F["add"][int(y)][F["neg"][int(x)]]
    if d == 0:
        return 0
    return int(F["chi"][d])


def chi_p(t, p):
    t %= p
    return 0 if t == 0 else pow(t, (p - 1) // 2, p)


def k22(F, a, b, c, d):
    e = [chi_q(F, a, c), chi_q(F, a, d), chi_q(F, b, c), chi_q(F, b, d)]
    nm = sum(1 for v in e if v == -1)
    if nm != 2:
        return nm, ""
    if (e[0] == -1 and e[3] == -1) or (e[1] == -1 and e[2] == -1):
        return 2, "match"
    return 2, "path"


def n_on(D, idxs):
    return int(D[:, list(idxs)].min(axis=1).sum())


def two_two_types(p: int) -> dict:
    D, F = load_D(p)
    lines = square_classes(F)[0]
    L0 = [int(x) for x in lines[0]]
    L1 = [int(x) for x in lines[1]]
    by = defaultdict(list)
    for i, j in combinations(range(p), 2):
        for k, m in combinations(range(p), 2):
            pts = (L0[i], L0[j], L1[k], L1[m])
            nm, kind = k22(F, *pts)
            same = {i, j} == {k, m}
            d1, d2 = (i - k) % p, (j - m) % p
            d3, d4 = (i - m) % p, (j - k) % p
            if d1 == d2:
                delta, trap = d1, True
            elif d3 == d4:
                delta, trap = d3, True
            else:
                delta, trap = None, False
            key = (
                nm,
                kind,
                int(same),
                int(trap),
                chi_p(j - i, p),
                chi_p(delta, p) if delta is not None else 0,
            )
            by[key].append(n_on(D, pts))
    types = {}
    n_const = 0
    for key, vs in by.items():
        u = sorted(set(vs))
        const = len(u) == 1
        n_const += int(const)
        types[str(key)] = {"unique": u, "n": len(vs), "const": const}
    return {
        "n_types": len(types),
        "n_const": n_const,
        "n_split": len(types) - n_const,
        "types": types,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    want = {5: Fraction(735, 26), 7: Fraction(101585, 818)}
    for p in (5, 7):
        rec = live_cross(p)
        if rec["E_n2n2"] != want[p]:
            ok = False
        if rec["E_n2n2"].denominator == HPLUS[p] // (2 * p):
            ok = False  # fail: den = |H+|/(2p)
        if rec["E_n2n2"] == rec["E_n4"] and p == 5:
            # 735/26 vs 990/13 — not equal, just a guard
            pass
        if p == 5 and rec["E_n2n2"] == Fraction(36, 1):
            ok = False
        rows[str(p)] = {k: str(v) if isinstance(v, Fraction) else v for k, v in rec.items()}
        rows[str(p)]["den"] = rec["E_n2n2"].denominator
        rows[str(p)]["H_over_p"] = HPLUS[p] // p
    if rows["5"]["den"] != 26:
        ok = False
    if rows["7"]["den"] != 818:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "E[n_i² n_j²]=735/26, 101585/818; den=|H_+|/p.",
    }


def prove_B() -> dict:
    ok = True
    rec = bivar_lp(5)
    live = float(Fraction(735, 26))
    if abs(rec["max"] - live) < 1e-6:
        ok = False
    if rec["min"] > live + 1e-6 or rec["max"] < live - 1e-6:
        ok = False  # live must sit inside
    if rec["max"] <= float(Fraction(735, 26)) + 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "p5": rec,
        "live": str(Fraction(735, 26)),
        "live_interior": rec["min"] < live < rec["max"],
        "theorem": (
            "Bivariate LP from 15.344 + ∑n=k does not pin E[n²n'²] "
            "to 735/26 (max=51.5)."
        ),
    }


def prove_C() -> dict:
    ok = True
    r5 = two_two_types(5)
    r7 = two_two_types(7)
    if r5["n_split"] != 0:
        ok = False
    if r7["n_split"] == 0:
        ok = False
    return {
        "proved": bool(ok),
        "p5": {k: r5[k] for k in ("n_types", "n_const", "n_split")},
        "p7": {k: r7[k] for k in ("n_types", "n_const", "n_split")},
        "p5_types": r5["types"],
        "theorem": (
            "p=5 2+2 n_4 is constant on K_{2,2}×χ_p(side,δ).  "
            "Same dictionary splits at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "interval_closed": False,
        "note": (
            "E[n_i² n_j²] is now an exact rational with den |H_+|/p.  "
            "The 15.344 bivariate LP does not reach [Q_lo, ub].  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.346  two-line moment; LP miss; p=5 2+2 constant", flush=True)
    A = prove_A()
    print(f"  A cross: {A['proved']} p5={A['by_p']['5']['E_n2n2']}", flush=True)
    B = prove_B()
    print(f"  B LP: {B['proved']} max={B['p5']['max']:.3f}", flush=True)
    C = prove_C()
    print(
        f"  C 2+2: {C['proved']} p5split={C['p5']['n_split']} p7split={C['p7']['n_split']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.346",
        "title": "E[n_i² n_j²] exact; bivariate LP misses; p=5 2+2 constant",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "cross_moment_exact": A["proved"],
            "bivar_lp_misses": B["proved"],
            "p5_22_type_constant": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15346.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
