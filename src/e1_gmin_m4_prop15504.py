#!/usr/bin/env python3
"""
Prop 15.504 — Autocorrelation of u=|ẑ|² on F_p^×-orbits of
Ω is at most two-level at p=5,7 (C constant on {±1}, and
constant on the complement) by evenness plus |(p−1)/2|∈{2,3}.
Not necessarily two distinct values (u may be constant on an
orbit).  Not a p-law: Type+ 1D at p=13 has three off values.
Does not name Q_τ.

    C_ξ(r) := mean_{a∈F_p^×} u(aξ) u(r a ξ),   r∈F_p^×.

A. u(−ξ)=u(ξ) (z real) ⇒ C(r)=C(−r).  Reindexing a↦ra gives
   C(r)=C(r^{−1}) with no evenness.  Fail: claim u(−ξ)=−u(ξ).

B. Sufficient, not a characterization: if (p−1)/2∈{2,3} then
   F_p^×/{±1} has at most one non-identity inv-pair, so A forces
   C constant on {±1} and constant on F_p^×\\{±1} (at most two
   levels; one if u is orbit-constant).  That is p=5,7 only and
   is the reason for 15.503.  Fail: claim the same count at p=13
   ((p−1)/2=6).  C kills the unrestricted p-law.

C. A Type+ 1D lift at p=13 has a unique live F_p^×-orbit in Ω
   on which C takes three distinct values on F_p^×\\{±1}.
   Fail: claim that 1D C is 2-level at p=13.

D. OPEN.  The tautology does not determine the off value, so
   Q_τ is unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.503 flags.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import (
    _e_tr,
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15290 import omega_set

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15504.json"


def G_over_H(p: int) -> int:
    return (p - 1) // 2


def C_on_orbit(u: np.ndarray, F, om_idx: dict, xi0: int) -> dict[int, float]:
    p = F["p"]
    vals = {a: float(u[om_idx[int(F["mul"][a][xi0])]]) for a in range(1, p)}
    C = {}
    for r in range(1, p):
        acc = 0.0
        for a in range(1, p):
            acc += vals[a] * vals[(r * a) % p]
        C[r] = acc / (p - 1)
    return C


def one_typeplus_1d(p: int):
    forms = [f for f in _linear_forms(p) if is_type_plus_form(p, *f)]
    plus = set(range((p + 1) // 2))
    y = lift_sigma_vector(p, *forms[0], plus)
    return np.asarray(y[1 : 1 + p * p], dtype=np.float64)


def u_of(z: np.ndarray, F, Omega) -> np.ndarray:
    q = F["q"]
    u = np.empty(len(Omega), dtype=np.float64)
    for i, xi in enumerate(Omega):
        s = 0j
        for x in range(q):
            s += z[x] * _e_tr(F, F["mul"][xi][x])
        u[i] = float(abs(s) ** 2)
    return u


def prove_A() -> dict:
    """u(−ξ)=u(ξ) for real z. Fail: odd."""
    F = _kernel_field(5)
    Omega = omega_set(F)
    z = one_typeplus_1d(5)
    u = u_of(z, F, Omega)
    om_idx = {int(xi): i for i, xi in enumerate(Omega)}
    ok = True
    n_check = 0
    for xi in Omega:
        eta = int(F["mul"][F["neg1"]][int(xi)])
        if abs(u[om_idx[int(xi)]] - u[om_idx[eta]]) > 1e-8:
            ok = False
        n_check += 1
        # fail-when-wrong: odd
        if abs(u[om_idx[int(xi)]] + u[om_idx[eta]]) < 1e-8 and u[om_idx[int(xi)]] > 1e-8:
            ok = False
    return {"proved": bool(ok), "n_check": n_check, "neg1": int(F["neg1"])}


def prove_B() -> dict:
    """(p−1)/2 ∈ {2,3} only at p=5,7. Fail p=13."""
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17):
        g = G_over_H(p)
        rows[p] = {"G_over_H": g, "forces_2level": g in (2, 3)}
    if rows[5]["G_over_H"] != 2 or rows[7]["G_over_H"] != 3:
        ok = False
    if rows[13]["G_over_H"] != 6:
        ok = False
    if rows[13]["forces_2level"]:
        ok = False
    if rows[11]["forces_2level"]:
        ok = False
    if not (rows[5]["forces_2level"] and rows[7]["forces_2level"]):
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """Type+ 1D at p=13: live orbit has ≥3 distinct off C. Fail 2-level."""
    p = 13
    F = _kernel_field(p)
    Omega = omega_set(F)
    z = one_typeplus_1d(p)
    u = u_of(z, F, Omega)
    om_idx = {int(xi): i for i, xi in enumerate(Omega)}
    seen = set()
    n_zero = 0
    n_live = 0
    max_n_off = 0
    max_spread = 0.0
    live_Coff: dict[int, float] | None = None
    for xi in Omega:
        xi = int(xi)
        if xi in seen:
            continue
        orb = [int(F["mul"][a][xi]) for a in range(1, p)]
        seen.update(orb)
        if max(float(u[om_idx[y]]) for y in orb) < 1e-12:
            n_zero += 1
            continue
        n_live += 1
        C = C_on_orbit(u, F, om_idx, orb[0])
        off = [C[r] for r in range(2, p - 1)]
        n_off = len({round(c, 6) for c in off})
        sp = max(off) - min(off)
        if n_off > max_n_off:
            max_n_off = n_off
        if sp > max_spread:
            max_spread = sp
            live_Coff = {int(r): float(C[r]) for r in range(2, p - 1)}
    ok = max_n_off >= 3 and max_spread > 1.0 and n_live >= 1
    if max_n_off == 1:
        ok = False
    return {
        "proved": bool(ok),
        "n_zero_orbits": n_zero,
        "n_live_orbits": n_live,
        "max_n_distinct_off": max_n_off,
        "max_off_spread": float(max_spread),
        "live_C_off": live_Coff,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "2-level C at p=5,7 is evenness plus |G/H|≤3, not a name "
            "of the off value. p=13 1D kills the p-law. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.504 two-level C is tautological at p=5,7; false at p=13 1D", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
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
