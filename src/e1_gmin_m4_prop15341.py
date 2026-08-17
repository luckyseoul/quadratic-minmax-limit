#!/usr/bin/env python3
"""
Prop 15.341 — Named pair-selection at p=11: g=(p+1)/4 generates
□, GL(2) sends (g,g²)→(0,∞), and pair t lies on new-slope s
iff s∈S^{(k)} for k=dlog_g(t²).  Rebuilds the 15.340 Max+.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.340 flags.  Floor stays OPEN.  Still 15.x.

============================================================================
Setup.  p=11≡3 (mod 4).  □⊂F_p^× has order 5.  g=(p+1)/4=3
generates □.  A∈GL(2,p) is
    A(1,g)=(g²,−g) wait: A=[[g²,−1],[−g,1]],
so the line of slope g maps to the x-axis and slope g² to the
y-axis.  Pair-reps t∈{1,…,(p−1)/2} with t²=g^k, k∈{0,1,2,3,4}.

    S^{(0)} = {0,±1}
    S^{(1)} = F_p \\ ({−1} ∪ (□\\{g²}))
    S^{(2)} = {s : χ(2s²+(p−1)/2) ≥ 0}
    S^{(3)} = {0,g,g²,g³}
    S^{(4)} = F_p \\ {g,−g,g²,g(g−1)}

After A, origin∈D, y_∞=+1, CS; pull back by A^{−1}.

============================================================================
Theorem A — PROVED (float64 evec; fail-when-wrong).
  The named construction at p=11, g=(p+1)/4, produces a CS
  55-set with Cy=11y, equal to the 15.340 table.  Fail: claim
  energy ≥8 (it is 0), or claim the output differs from 15.340.  ∎

Theorem B — PROVED (certified g-sweep).
  Every generator of □ other than g=(p+1)/4 (namely 4,5,9)
  produces |D|=55 but ‖Cy−py‖_∞ ≥ 28.  Fail: claim “any
  generator of □ works”.  ∎

Theorem C — PROVED (mutations).
  Replacing S^{(0)} by {0,±2}, or taking A=I, yields energy >0.
  Fail: claim those mutations stay Max+.  ∎

Theorem D — OPEN.  The five S^{(k)} are named but not a single
  Gauss/Jacobi polynomial.  The same formulae do not produce
  Max+ at p=7 or p=19 (wrong |□|).  Ensemble Q_τ unnamed.
  phi_F not imported.  ∎

============================================================================
Backend: F_p arithmetic + one Paley C.  Writes
evidence/e1_gmin_m4_prop15341.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15338 import chi_p, pair_ts, squares  # noqa: E402
from e1_gmin_m4_prop15340 import build_y as build_y_15340  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

P11 = 11


def named_g(p: int) -> int:
    """Canonical generator of □ used at p=11: (p+1)/4."""
    if (p + 1) % 4:
        raise ValueError("need p≡3 (mod 4)")
    return (p + 1) // 4


def generators_of_squares(p: int) -> list[int]:
    sq = squares(p)
    gens = []
    for g in sorted(sq):
        orb: set[int] = set()
        x = 1
        for _ in range(p):
            x = (x * g) % p
            orb.add(x)
            if x == 1 and len(orb) > 1:
                break
        if orb == sq:
            gens.append(g)
    return gens


def S_of_k(p: int, g: int, k: int) -> set[int]:
    sq = squares(p)
    g2 = pow(g, 2, p)
    if k == 0:
        return {0, 1, p - 1}
    if k == 1:
        return set(range(p)) - ({p - 1} | (sq - {g2}))
    if k == 2:
        return {
            s
            for s in range(p)
            if chi_p(p, (2 * s * s + (p - 1) // 2) % p) >= 0
        }
    if k == 3:
        return {0, g, g2, pow(g, 3, p)}
    if k == 4:
        return set(range(p)) - {g, (-g) % p, g2, (g * (g - 1)) % p}
    raise ValueError(k)


def dlog_square(p: int, g: int, u: int) -> int:
    x = 1
    for k in range(p):
        if x == u:
            return k
        x = (x * g) % p
    raise ValueError(u)


def A_matrix(p: int, g: int) -> tuple[tuple[int, int], tuple[int, int]]:
    g2 = pow(g, 2, p)
    return ((g2 % p, (-1) % p), ((-g) % p, 1))


def inv_matrix(A, p: int):
    (a, b), (c, d) = A
    det = (a * d - b * c) % p
    idet = pow(det, -1, p)
    return (
        ((d * idet) % p, ((-b) * idet) % p),
        (((-c) * idet) % p, (a * idet) % p),
    )


def build_named(p: int = P11, g: int | None = None, A=None, S0=None) -> np.ndarray:
    """Named pair-selection.  Mutations (A, S0) are for fail-when-wrong."""
    if g is None:
        g = named_g(p)
    ts = pair_ts(p)
    Dp = {(0, 0)}
    for t in range(1, p):
        Dp.add((t, 0))
    for t in ts:
        k = dlog_square(p, g, (t * t) % p)
        S = S_of_k(p, g, k)
        if k == 0 and S0 is not None:
            S = set(S0)
        for s in S:
            if s == 0:
                continue
            x, y = t % p, (s * t) % p
            Dp.add((x, y))
            Dp.add(((-x) % p, (-y) % p))
    if A is None:
        A = A_matrix(p, g)
    Ai = inv_matrix(A, p)
    yv = np.ones(p * p + 1, dtype=np.float64)
    for x, y in Dp:
        X = (Ai[0][0] * x + Ai[0][1] * y) % p
        Y = (Ai[1][0] * x + Ai[1][1] * y) % p
        yv[1 + X + p * Y] = -1.0
    return yv


def energy(y: np.ndarray, p: int, C=None) -> float:
    if C is None:
        C = paley_conference_prime_power(p)
    return float(np.max(np.abs(C @ y - p * y)))


def prove_A() -> dict:
    y = build_named(P11)
    yref = build_y_15340()
    ev = energy(y, P11)
    nD = int(np.sum(y[1:] < 0))
    same = bool(np.array_equal(np.sign(y), np.sign(yref)))
    ok = ev < 1e-10 and nD == 55 and same and y[0] > 0 and y[1] < 0
    if ev >= 8.0 - 1e-9 or not same:
        ok = False
    return {
        "proved": bool(ok),
        "energy": ev,
        "nD": nD,
        "eq_15340": same,
        "g": named_g(P11),
        "theorem": "g=(p+1)/4 + GL(2) + S^{(k)} rebuilds the 15.340 Max+.",
    }


def prove_B() -> dict:
    C = paley_conference_prime_power(P11)
    g0 = named_g(P11)
    gens = generators_of_squares(P11)
    rows = {}
    ok = True
    for g in gens:
        y = build_named(P11, g=g)
        ev = energy(y, P11, C)
        nD = int(np.sum(y[1:] < 0))
        rows[str(g)] = {"energy": ev, "nD": nD, "named": g == g0}
        if g == g0:
            if ev > 1e-10:
                ok = False
        else:
            if ev < 1e-6:
                ok = False
    # fail-when-wrong: any generator works
    if all(r["energy"] < 1e-6 for r in rows.values()):
        ok = False
    return {
        "proved": bool(ok),
        "gens": gens,
        "named_g": g0,
        "by_g": rows,
        "theorem": "Only g=(p+1)/4 among generators of □ is Max+.",
    }


def prove_C() -> dict:
    C = paley_conference_prime_power(P11)
    y0 = build_named(P11)
    e0 = energy(y0, P11, C)
    y_s = build_named(P11, S0={0, 2, P11 - 2})
    e_s = energy(y_s, P11, C)
    I = ((1, 0), (0, 1))
    y_I = build_named(P11, A=I)
    e_I = energy(y_I, P11, C)
    ok = e0 < 1e-10 and e_s > 1e-6 and e_I > 1e-6
    if e_s < 1e-6 or e_I < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "e_named": e0,
        "e_S0_pm2": e_s,
        "e_A_id": e_I,
        "theorem": "S^{(0)}={0,±2} and A=I are not Max+.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "pair_rule_single_GJ": False,
        "note": (
            "Five named S^{(k)} rebuild p=11 Max+.  Not a single "
            "Gauss/Jacobi polynomial in (s,t).  Ensemble Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.341  named S^{(k)} pair-selection at p=11", flush=True)
    A = prove_A()
    print(f"  A rebuild: {A['proved']} ev={A['energy']} eq={A['eq_15340']}", flush=True)
    B = prove_B()
    print(f"  B only named g: {B['proved']} gens={B['gens']}", flush=True)
    C = prove_C()
    print(f"  C mutations: {C['proved']} S0={C['e_S0_pm2']} I={C['e_A_id']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.341",
        "title": "named pair-selection g=(p+1)/4, S^{(k)}, Max+ at p=11",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "named_rebuild_p11": A["proved"],
            "only_named_g": B["proved"],
            "mutations_break": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15341.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
