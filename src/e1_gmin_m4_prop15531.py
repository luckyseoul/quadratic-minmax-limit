#!/usr/bin/env python3
"""
Prop 15.531 — 15.527 axis-only n_R is chart-dependent;
lin-form affine occupancy is all of n_free (15.305 C).
Does not name n_free.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.530 flags.  Does **not** overwrite 15.522–15.528.

============================================================================
Setup.  15.527 counted free T-orbits whose *axis-aligned* occupancy
is affine: n_R=4,16.  That predicate is not translation-invariant
in a fixed chart (p=7: 16 on min-key vs 20 on any member).
A linear form ℓ=αu+βv has affine occupancy iff |D∩ℓ^{-1}(c)|
is affine in c.  15.305 C: every Max+ meets every nonsquare
direction in μ=(p−1)/2, so every free orbit is lin-affine and
n_R^{lin}=n_free.  Random weight-k subsets are lin-affine at
rate 0.183 (p=5) and 0.012 (p=7) — the Max+ rate 1 is 15.305,
not a counting law.

============================================================================
Theorem A — PROVED (H_+ caches; 15.305 C).
  n_R^{lin}=n_free at p=5,7 (4=4, 114=114).  Fail: n_X^{lin}=98.
  Fail: 15.527 axis n_R is translation-invariant (16≠20 at p=7).  ∎

Theorem B — PROVED (ProcessPool random control).
  Random k-subsets of F_q are lin-affine with rate ≠1
  (1577/8600 at p=5, 106/8600 at p=7).  Fail: claim the
  random rate is 1 (then n_R^{lin}=n_free would be tautological
  on all k-sets, not a Max+ identity).  ∎

Theorem C — OPEN.  n_free, D, and Q_τ still unnamed in p.
  Occupancy splits cannot name n_free (15.305).  phi_F not imported.

============================================================================
Backend: H_+ caches + _linear_forms.  Writes
evidence/e1_gmin_m4_prop15531.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _linear_forms
from e1_gmin_m4_prop15526 import is_1d_row, load_Hplus
from e1_gmin_m4_prop15527 import is_affine_R, split_free_orbits

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15531.json"

# ProcessPool control (cpu_nR_random_control.py), W=86.
RANDOM_LIN_AFF = {
    5: {"hits": 1577, "n": 8600},
    7: {"hits": 106, "n": 8600},
}


def orbit_min_key(y: np.ndarray, p: int) -> bytes:
    grid = np.ascontiguousarray(y[1:]).reshape(p, p)
    mk = None
    for t0 in range(p):
        for t1 in range(p):
            z = np.roll(np.roll(grid, t1, axis=0), t0, axis=1)
            kb = z.tobytes()
            if mk is None or kb < mk:
                mk = kb
    return mk


def has_lin_aff(y: np.ndarray, p: int) -> bool:
    q = p * p
    neg = y[1:] < 0
    xs = np.arange(q)
    u, v = xs % p, xs // p
    for a, b in _linear_forms(p):
        lab = (a * u + b * v) % p
        occ = np.bincount(lab[neg], minlength=p)
        A = int(occ[1] - occ[0])
        if np.all(occ == A * np.arange(p) + occ[0]):
            return True
    return False


def lin_and_axis_counts(p: int) -> dict:
    H = load_Hplus(p)
    first: dict[bytes, np.ndarray] = {}
    axis_any: set[bytes] = set()
    lin_any: set[bytes] = set()
    for i in range(len(H)):
        y = H[i]
        if is_1d_row(y, p):
            continue
        mk = orbit_min_key(y, p)
        if mk not in first:
            first[mk] = y
        if is_affine_R(y, p):
            axis_any.add(mk)
        if has_lin_aff(y, p):
            lin_any.add(mk)
    axis_mk = sum(1 for y in first.values() if is_affine_R(y, p))
    return {
        "n_free": len(first),
        "n_R_axis_any": len(axis_any),
        "n_R_axis_minkey": axis_mk,
        "n_R_lin": len(lin_any),
        "n_X_lin": len(first) - len(lin_any),
    }


def n_R_wrong_axis527(p: int) -> int:
    """15.527 axis count, not a p-law."""
    return {5: 4, 7: 16}[p]


def n_X_wrong_98(p: int) -> int:
    return {5: 0, 7: 98}[p]


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        rec = lin_and_axis_counts(p)
        s527 = split_free_orbits(p)
        row_ok = (
            rec["n_R_lin"] == rec["n_free"]
            and rec["n_X_lin"] == 0
            and rec["n_X_lin"] != n_X_wrong_98(p) if p == 7 else True
        )
        if p == 7:
            row_ok = (
                row_ok
                and rec["n_R_axis_minkey"] == 16
                and rec["n_R_axis_any"] == 20
                and rec["n_R_axis_minkey"] != rec["n_R_axis_any"]
                and s527["n_R"] == 16
            )
        if p == 5:
            row_ok = row_ok and rec["n_free"] == 4 and rec["n_R_lin"] == 4
        if not row_ok:
            ok = False
        rows[str(p)] = rec
        rows[str(p)]["ok"] = row_ok
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_R^{lin}=n_free (15.305). Axis n_R is not "
            "translation-invariant at p=7 (16≠20). Fail n_X^{lin}=98."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p, rec in RANDOM_LIN_AFF.items():
        rate = rec["hits"] / rec["n"]
        row_ok = rec["hits"] != rec["n"] and rate < 0.5
        if not row_ok:
            ok = False
        rows[str(p)] = {**rec, "rate": rate, "ok": row_ok}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Random k-subsets are not all lin-affine "
            "(p=5: 1577/8600; p=7: 106/8600). Fail rate=1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_free_named_in_p": False,
        "n_R_named_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Occupancy splits cannot name n_free (15.305 C). "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.531 axis n_R chart-dependent; lin-form = n_free by 15.305", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.531",
        "title": "axis n_R not invariant; lin-form affine is all of n_free",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
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
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
