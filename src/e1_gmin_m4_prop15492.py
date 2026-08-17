#!/usr/bin/env python3
"""
Prop 15.492 — Per-Max+ Var(e_*) is not a single value and is not
bounded by V_floor/κ on the Type+ 1D family.

15.326 D restated the p≡1 floor as Var(e_*)≤V_floor/κ (18/7 at
p=5).  That cannot be a per-design bound: at p=5 the live 1D
lifts (n_1d=30) all have Var(e_*)=5>18/7, while the 100 NL
designs have Var=9/5<18/7.  The cited 33/13 is the mixture
mean.  Fourier energy at α≠0 is complement-invariant, so the
choice {y=+1} vs {y=−1} as D does not change Var.

Fail-when-wrong: claim every Max+ has Var≤18/7; claim a single
Var on all 130; claim the 1d count is not n_1d.

Does not name Q_τ or prove ⟨δ,ψ⟩≤2.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.491 flags.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15294 import _fp_dirs, _is_1d
from e1_gmin_m4_prop15326 import V_floor, kappa_named, theta_plus_norms

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15492.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _field_norm(F, r: int) -> int:
    p = F["p"]
    acc, b, ee = 1, int(r), p + 1
    mul = F["mul"]
    while ee:
        if ee & 1:
            acc = mul[acc][b]
        b = mul[b][b]
        ee >>= 1
    return acc


def _slices(F):
    q = F["q"]
    want = set(theta_plus_norms(F["p"]))
    out = {s: [] for s in want}
    for a in range(1, q):
        n = _field_norm(F, a)
        if n in want:
            out[n].append(a)
    return out, want


def var_e_of_y(F, y, slices, want) -> float:
    q = F["q"]
    D = [x for x in range(q) if y[x] > 0]
    e_s = []
    for s in sorted(want):
        acc = 0.0
        for a in slices[s]:
            h = 0j
            for x in D:
                h += _e_tr(F, F["mul"][a][x])
            acc += h.real * h.real + h.imag * h.imag
        e_s.append(acc / q)
    return float(np.asarray(e_s, dtype=np.float64).var(ddof=0))


def census_p5() -> dict:
    p = 5
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    dirs = _fp_dirs(p)
    slices, want = _slices(F)
    by = defaultdict(int)
    n1 = 0
    for yp in Yp:
        z = yp[1 : 1 + q]
        one = bool(_is_1d(z, dirs, F["add"], p))
        if one:
            n1 += 1
        ve = var_e_of_y(F, z, slices, want)
        by[(one, round(ve, 10))] += 1
    return {
        "n_total": int(len(Yp)),
        "n_1d_flag": n1,
        "n_1d_named": int(n_1d(p)),
        "by": {f"1d={a}|Var={v}": c for (a, v), c in sorted(by.items())},
        "need": float(V_floor(p) / kappa_named(p)),
    }


def prove_A() -> dict:
    """p=5: 1d all Var=5>18/7; NL all Var=9/5; counts match n_1d."""
    row = census_p5()
    ok = True
    if row["n_1d_flag"] != row["n_1d_named"] or row["n_1d_named"] != 30:
        ok = False
    if row["n_total"] != 130:
        ok = False
    if row["by"].get("1d=True|Var=5.0") != 30:
        ok = False
    if row["by"].get("1d=False|Var=1.8") != 100:
        ok = False
    # fail-when-wrong: a single value, or all ≤ 18/7
    if len(row["by"]) == 1:
        ok = False
    if 5.0 <= row["need"]:
        ok = False
    return {"proved": bool(ok), "row": row}


def prove_B() -> dict:
    """p=7: 1d count is n_1d=140 and is not a single Var."""
    p = 7
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    dirs = _fp_dirs(p)
    slices, want = _slices(F)
    n1 = 0
    vars_1d = defaultdict(int)
    for yp in Yp:
        z = yp[1 : 1 + q]
        if not _is_1d(z, dirs, F["add"], p):
            continue
        n1 += 1
        ve = var_e_of_y(F, z, slices, want)
        vars_1d[round(ve, 8)] += 1
    ok = n1 == n_1d(p) == 140 and len(vars_1d) >= 2
    # fail: claim 1d is a single Var
    if len(vars_1d) == 1:
        ok = False
    return {
        "proved": bool(ok),
        "n_1d_flag": n1,
        "n_1d_named": int(n_1d(p)),
        "vars_1d": {str(k): int(v) for k, v in sorted(vars_1d.items())},
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Per-D Var is not the 15.326 floor. Ensemble mean 33/13 "
            "still unnamed as a law. ⟨δ,ψ⟩≤2 open."
        ),
    }


def main() -> dict:
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "A": A,
        "B": B,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
