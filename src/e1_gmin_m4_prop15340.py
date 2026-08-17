#!/usr/bin/env python3
"""
Prop 15.340 — The 15.338 k-vector is realised by a Max+ at p=11.
Pair-selection is not the p=7 interpolant and not opposite-high.
NC-degree ≥ 1 (not the isolated p=7 type).  Still 15.x.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.339 flags.  Floor stays OPEN.  Ensemble Q_τ unnamed.

============================================================================
Setup.  p=11 ≡ 3 (mod 4).  15.338 k-vector on the p+1 AG-lines
through 0: one full k=5, one empty k=0, n_high=2 with k=3,
n_low=8 with k=2.  |D|=55.  y_∞=+1, origin ∈ D, D centrally
symmetric.

============================================================================
Theorem A — PROVED (float64 conference evec; fail-when-wrong).
  The pair-table
      inf:{1,2},  s0:{1,3,4},  s1:{1,2},  s2:{1,2},  s3:{1..5},
      s4:{1,3},   s5:{2,4},    s6:{3,4},  s7:{1,3},  s8:{2,3,4},
      s9:∅,       s10:{2,3}
  rebuilds a unique CS y with Cy=11y and k-hist (5,3²,2⁸,0).
  Roles: full=s3, empty=s9, high={s0,s8}, low=the rest
  (vertical is low; (full,empty) is not an opposite pair).
  Fail: claim energy ≥8 on every 15.338 k-vector at p=11
  (this table has energy 0), or claim (full,empty)={0,∞}.  ∎

Theorem B — PROVED (certified identities).
  The table is not the 15.338 interpolant α (that D has size 59).
  The y is not a Type+ halfspace, not the QR0⊙NC seed, not ystar.
  first_nc(y) is nonempty: NC-degree ≥ 1, so this is not the
  isolated p=7 type.  Fail: claim the hit equals the interpolant
  or has NC-degree 0.  ∎

Theorem C — PROVED as a negative (complete deg≤2 census).
  No bivariate deg≤2 f∈F_11[s,t] has
      (t,st)∈D  ⇔  χ(f(s,t))≥0
  on the 50 affine pairs.  Fail: claim a deg≤2 square predicate
  reproduces the table.  ∎

Theorem D — OPEN.  Pair-selection is not yet a Gauss/Jacobi formula
  in p.  Ensemble Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: one Paley C at p=11 + deg≤2 coefficient scan.
Writes evidence/e1_gmin_m4_prop15340.json
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

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc  # noqa: E402
from e1_gmin_m4_prop15334 import construct_qr0nc, is_halfspace_y  # noqa: E402
from e1_gmin_m4_prop15338 import (  # noqa: E402
    build_D_interpolant,
    chi_p,
    k_full,
    k_high,
    k_low,
    n_high,
    n_low,
    pair_ts,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

P = 11
TS = (1, 2, 3, 4, 5)
# slope s=0..10 and "inf"; values are pair-reps t
PAIR = {
    "inf": frozenset({1, 2}),
    0: frozenset({1, 3, 4}),
    1: frozenset({1, 2}),
    2: frozenset({1, 2}),
    3: frozenset({1, 2, 3, 4, 5}),
    4: frozenset({1, 3}),
    5: frozenset({2, 4}),
    6: frozenset({3, 4}),
    7: frozenset({1, 3}),
    8: frozenset({2, 3, 4}),
    9: frozenset(),
    10: frozenset({2, 3}),
}


def build_y(pair=None) -> np.ndarray:
    pair = PAIR if pair is None else pair
    p = P
    y = np.ones(p * p + 1, dtype=np.float64)
    y[1] = -1.0
    for t in pair["inf"]:
        y[1 + 0 + p * t] = -1.0
        y[1 + 0 + p * ((-t) % p)] = -1.0
    for s in range(p):
        for t in pair[s]:
            x, yy = t % p, (s * t) % p
            y[1 + x + p * yy] = -1.0
            y[1 + ((-x) % p) + p * ((-yy) % p)] = -1.0
    return y


def energy(y: np.ndarray, C=None) -> float:
    if C is None:
        C = paley_conference_prime_power(P)
    return float(np.max(np.abs(C @ y - P * y)))


def kvec_of(y: np.ndarray) -> dict[str, int]:
    p = P
    D = {((i - 1) % p, (i - 1) // p) for i in range(1, p * p + 1) if y[i] < 0}
    out = {}
    half = (p - 1) // 2
    # vertical
    out["inf"] = sum((0, t) in D for t in range(1, half + 1))
    out["s0"] = sum((t, 0) in D for t in range(1, half + 1))
    for s in range(1, p):
        out[f"s{s}"] = sum((t, (s * t) % p) in D for t in range(1, half + 1))
    return out


def roles_of(kv: dict[str, int]) -> dict[str, list[str]]:
    roles = {"full": [], "empty": [], "high": [], "low": []}
    kf, kh, kl = k_full(P), k_high(P), k_low(P)
    for name, k in kv.items():
        if k == kf:
            roles["full"].append(name)
        elif k == 0:
            roles["empty"].append(name)
        elif k == kh:
            roles["high"].append(name)
        elif k == kl:
            roles["low"].append(name)
    return roles


def prove_A() -> dict:
    y = build_y()
    C = paley_conference_prime_power(P)
    ev = energy(y, C)
    kv = kvec_of(y)
    roles = roles_of(kv)
    hist = {}
    for k in kv.values():
        hist[k] = hist.get(k, 0) + 1
    nD = int(np.sum(y[1:] < 0))
    ok = ev < 1e-10
    ok = ok and nD == P * (P - 1) // 2
    ok = ok and y[0] > 0 and y[1] < 0
    ok = ok and hist.get(k_full(P), 0) == 1
    ok = ok and hist.get(0, 0) == 1
    ok = ok and hist.get(k_high(P), 0) == n_high(P)
    ok = ok and hist.get(k_low(P), 0) == n_low(P)
    ok = ok and roles["full"] == ["s3"]
    ok = ok and roles["empty"] == ["s9"]
    ok = ok and set(roles["high"]) == {"s0", "s8"}
    # fail-when-wrong: (full, empty) = {0, ∞}
    if set(roles["full"] + roles["empty"]) == {"s0", "inf"}:
        ok = False
    # fail-when-wrong: energy ≥ 8
    if ev >= 8.0 - 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "energy": ev,
        "nD": nD,
        "k_hist": {str(k): v for k, v in sorted(hist.items())},
        "roles": roles,
        "kvec": kv,
        "theorem": (
            "Named 15.338 k-vector is Max+ at p=11 for this pair-table. "
            "(full,empty)=(s3,s9), not {0,∞}."
        ),
    }


def prove_B() -> dict:
    y = build_y()
    inter = build_D_interpolant(P)
    D = {((i - 1) % P, (i - 1) // P) for i in range(1, P * P + 1) if y[i] < 0}
    D_inter = {((i) % P, (i) // P) for i in range(P * P) if inter[i]}
    hs = bool(is_halfspace_y(P, y))
    qr = construct_qr0nc(P)
    ys = construct_ystar(P)
    same_qr = bool(
        qr.get("found") and np.array_equal(np.sign(qr["y"]), np.sign(y))
    )
    same_ys = bool(
        ys.get("found") and np.array_equal(np.sign(ys["y"]), np.sign(y))
    )
    nc = first_nc(P, y)
    ok = D != D_inter
    ok = ok and len(D_inter) == 59
    ok = ok and not hs and not same_qr and not same_ys
    ok = ok and nc is not None
    if D == D_inter or nc is None:
        ok = False
    return {
        "proved": bool(ok),
        "eq_interpolant": D == D_inter,
        "interpolant_size": int(len(D_inter)),
        "halfspace": hs,
        "same_qr0": same_qr,
        "same_ystar": same_ys,
        "has_nc": nc is not None,
        "theorem": "Not interpolant / 1D / QR0NC / ystar; NC-degree ≥ 1.",
    }


def prove_C() -> dict:
    """Complete scan of deg≤2 square predicates on the 50 affine pairs."""
    rows = []
    for s in range(1, P):
        for t in TS:
            rows.append((s, t, 1 if t in PAIR[s] else 0))
    n_hit = 0
    sample = None
    # 11^6 = 1771561; serial is <2s
    for a00 in range(P):
        for a10 in range(P):
            for a01 in range(P):
                for a20 in range(P):
                    for a11 in range(P):
                        for a02 in range(P):
                            good = True
                            for s, t, bit in rows:
                                val = (
                                    a00
                                    + a10 * s
                                    + a01 * t
                                    + a20 * ((s * s) % P)
                                    + a11 * ((s * t) % P)
                                    + a02 * ((t * t) % P)
                                ) % P
                                pred = 1 if chi_p(P, val) >= 0 else 0
                                if pred != bit:
                                    good = False
                                    break
                            if good:
                                n_hit += 1
                                sample = [a00, a10, a01, a20, a11, a02]
    ok = n_hit == 0
    if n_hit:
        ok = False
    return {
        "proved": bool(ok),
        "n_deg2": n_hit,
        "sample": sample,
        "n_bits": len(rows),
        "theorem": "No deg≤2 χ-predicate reproduces the p=11 pair-table.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "pair_rule_named_in_p": False,
        "note": (
            "Existence of a 15.338 Max+ at p=11 is closed.  "
            "The pair-selection is not a deg≤2 square predicate and "
            "is not yet named in p.  Ensemble Q_τ unnamed.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.340  p=11 k-vector Max+ exists; pair-rule unnamed", flush=True)
    A = prove_A()
    print(f"  A exists: {A['proved']} energy={A['energy']} roles={A['roles']}", flush=True)
    B = prove_B()
    print(f"  B not-isolated: {B['proved']} has_nc={B['has_nc']}", flush=True)
    C = prove_C()
    print(f"  C deg2 empty: {C['proved']} n={C['n_deg2']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.340",
        "title": "15.338 k-vector is Max+ at p=11; pair-rule unnamed in p",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "p11_kvec_maxplus_exists": A["proved"],
            "not_interpolant_not_isolated": B["proved"],
            "no_deg2_square_predicate": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15340.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
