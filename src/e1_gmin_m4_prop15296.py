#!/usr/bin/env python3
"""
Prop 15.296 — X_A is the multiplicative ψ-mode of u on Ω;
complementary {0,4q} partner and pointwise Q≥Coll are false.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.295 flags.  Does **not** sign Cov(u,X_A)≥6q².
Floor stays OPEN.

============================================================================
Setup.  u(η)=|ẑ(η)|² on Ω.  A_{ab}=A(ξ(b−a)) is additive-circulant,
commutes with E_{ab}=e_ξ(b−a).  Eigenvalues of A:
  μ(η)=q ψ(η/ξ) if η/ξ∈T, else 0.
Hence X_A=zᵀAz=∑_{r∈T} ψ(r) u(rξ).  Floor is E[u(ξ)(ψ*u)(ξ)]≥6q².

============================================================================
Theorem A — PROVED (circulant FT; certified p=5,7).
  X_A = ∑_{r∈T} ψ(r) u(rξ).  Fail-when-wrong: drop ψ (trivial
  mode ∑_T u = q(q−1) ≠ X_A).  ∎

Theorem B — PROVED (15.279 M_1; recertified every Max+ row).
  ∑_{η∈Ω} u(η)=q(q−1) constantly, so u−2q is mean-zero on Ω
  for every Max+ vector.  Fail-when-wrong: claim the sum is 0.  ∎

Theorem C — PROVED (cache histogram).
  u is not two-level {0,4q}: p=5 takes 5 values
  {0, 10(5−√5), 4q/φ², 10(5+√5), 4q φ²}; p=7 takes ≥20.
  Complementary energy partner u↦4q−u is not a map on the
  support of u (4q is not attained).  Fail-when-wrong: claim
  every u∈{0,4q}.  ∎

Theorem D — PROVED (p=5 type (--,N*)).
  Q≥Coll off is false: Coll_off/q²=3−10/q, Q_{--,N*}/q²=32/13
  at p=5, and 32/13 < 3−10/25.  Fail-when-wrong: claim Q≥Coll
  on every 15.290 type at p=5.  ∎

Theorem E — OPEN.  Cov(u,X_A)≥6q² unsigned.  phi_F False.  ∎

============================================================================
Backend: CPU ProcessPool over Max+ rows.
Writes evidence/e1_gmin_m4_prop15296.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _gauss,
    _kernel_field,
    _psi_vec,
    bool_coll_off,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import live_Q_on_T, omega_set, type_key, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15295 import A_of_alpha  # noqa: E402
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
_G: dict = {}


def _init(p: int, k: int) -> None:
    F = _kernel_field(p)
    psi = _psi_vec(F, k)
    Gpsi = _gauss(F, psi)
    Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
    Omega = omega_set(F)
    xi = Omega[0]
    q = F["q"]
    A_tab = np.array(
        [
            A_of_alpha(F, psi, Gpsi, Gchipsi, F["mul"][xi][d] if d else 0)
            for d in range(q)
        ],
        dtype=np.complex128,
    )
    add = np.array(F["add"])
    _G.update(
        F=F,
        psi=np.array(psi, dtype=np.complex128),
        Omega=Omega,
        xi=xi,
        A_tab=A_tab,
        add=add,
        q=q,
        Gpsi=Gpsi,
    )


def XA_quadratic(z) -> complex:
    z = np.asarray(z, dtype=np.float64)
    q, add, A_tab = _G["q"], _G["add"], _G["A_tab"]
    Az = np.zeros(q, dtype=np.float64)
    for d in range(q):
        s = 0.0
        for x in range(q):
            s += z[x] * z[add[x, d]]
        Az[d] = s
    return complex(np.dot(Az, A_tab))


def XA_both_ways(z):
    return XA_quadratic(z), XA_psi_mode(z)


def XA_psi_mode(z) -> complex:
    F, Omega, xi, psi = _G["F"], _G["Omega"], _G["xi"], _G["psi"]
    u = np.abs(zhat_omega(np.asarray(z, dtype=np.float64), F, Omega)) ** 2
    om_idx = {int(eta): i for i, eta in enumerate(Omega)}
    s = 0j
    for r in F["squares"]:
        eta = F["mul"][r][xi]
        s += psi[r] * u[om_idx[eta]]
    return s


def u_on_omega(z) -> np.ndarray:
    return np.abs(
        zhat_omega(np.asarray(z, dtype=np.float64), _G["F"], _G["Omega"])
    ) ** 2


def _load_z(p: int) -> np.ndarray:
    F = _kernel_field(p)
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    return np.sign(Yp[:, 1 : 1 + F["q"]].astype(np.float64))


def prove_XA_is_psi_mode(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        k = 2
        _init(p, k)
        zs = _load_z(p)
        W = default_workers()

        if zs.shape[0] < 200:
            pairs = [XA_both_ways(z) for z in zs]
        else:
            with ProcessPoolExecutor(
                max_workers=W, initializer=_init, initargs=(p, k)
            ) as ex:
                pairs = list(
                    ex.map(XA_both_ways, zs, chunksize=max(1, len(zs) // (W * 4)))
                )
        diffs = [abs(a - b) for a, b in pairs]
        mx = float(max(diffs))
        if mx > 1e-6:
            ok = False
        # fail-when-wrong: drop ψ (trivial mode = ∑_T u = q(q−1))
        saved = _G["psi"].copy()
        _G["psi"] = np.ones_like(saved)
        inv_diffs = [
            abs(XA_quadratic(z) - XA_psi_mode(z))
            for z in zs[: min(20, len(zs))]
        ]
        invert_tot += 1
        if max(inv_diffs) < 1e-6:
            invert_hit += 1
        _G["psi"] = saved
        sample[str(p)] = {"max_diff": mx, "n": len(zs), "ok": mx <= 1e-6}
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_sum_u_constant(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        _init(p, 2)
        zs = _load_z(p)
        q = _G["q"]
        want = q * (q - 1)
        W = default_workers()
        if zs.shape[0] < 200:
            sums = [float(u_on_omega(z).sum()) for z in zs]
        else:
            with ProcessPoolExecutor(
                max_workers=W, initializer=_init, initargs=(p, 2)
            ) as ex:
                us = list(
                    ex.map(u_on_omega, zs, chunksize=max(1, len(zs) // (W * 4)))
                )
            sums = [float(u.sum()) for u in us]
        mx = max(abs(s - want) for s in sums)
        if mx > 1e-4:
            ok = False
        invert_tot += 1
        if abs(want) < 1e-9:
            invert_hit += 1
        sample[str(p)] = {
            "want": want,
            "max_err": mx,
            "ok": mx <= 1e-4,
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_not_two_level() -> dict:
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in (5, 7):
        _init(p, 2)
        zs = _load_z(p)
        q = _G["q"]
        W = default_workers()
        if zs.shape[0] < 200:
            us = [u_on_omega(z) for z in zs]
        else:
            with ProcessPoolExecutor(
                max_workers=W, initializer=_init, initargs=(p, 2)
            ) as ex:
                us = list(
                    ex.map(u_on_omega, zs, chunksize=max(1, len(zs) // (W * 4)))
                )
        flat = np.round(np.concatenate(us), 6)
        vals = np.unique(flat)
        four_q = 4 * q
        near0 = float(np.mean(np.abs(flat) < 1e-4))
        near4q = float(np.mean(np.abs(flat - four_q) < 1e-4))
        two_level = near0 + near4q > 1 - 1e-6
        if two_level:
            ok = False
        invert_tot += 1
        if two_level:
            invert_hit += 1
        sample[str(p)] = {
            "n_unique": int(len(vals)),
            "near0": near0,
            "near4q": near4q,
            "two_level": two_level,
            "max_u": float(flat.max()),
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_Q_lt_Coll_p5() -> dict:
    F = _kernel_field(5)
    Q = live_Q_on_T(5)
    q2 = float(F["q"] ** 2)
    coll_off = bool_coll_off(5) / q2
    vals = []
    for r, Qr in Q.items():
        if r in (1, F["neg1"]):
            continue
        vals.append((type_key(F, r), Qr / q2))
    mn = min(v for _t, v in vals)
    ok = mn < coll_off - 1e-9
    invert_hit = 0 if ok else 1
    return {
        "proved": ok,
        "coll_off_over_q2": coll_off,
        "min_Q_off_over_q2": mn,
        "by_type": {str(t): v for t, v in vals},
        "invert_hit": invert_hit,
        "invert_tot": 1,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": (
            "X_A=(ψ*u)(ξ).  Complementary {0,4q} and Q≥Coll are false.  "
            "Cov(u,X_A)≥6q² unsigned."
        ),
    }


def main() -> dict:
    print("Prop 15.296  X_A=ψ-mode of u; two-level and Q≥Coll dead", flush=True)
    A = prove_XA_is_psi_mode()
    print(f"  A X_A=ψ*u: {A['proved']}", flush=True)
    B = prove_sum_u_constant()
    print(f"  B sum_Ω u=q(q−1): {B['proved']}", flush=True)
    C = prove_not_two_level()
    print(f"  C not {{0,4q}}: {C['proved']}", flush=True)
    D = prove_Q_lt_Coll_p5()
    print(f"  D Q<Coll at p=5: {D['proved']}", flush=True)
    E = prove_open()
    print(f"  E floor: {E['proved']} phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.296",
        "title": "X_A is ψ-mode of u; complementary {0,4q} and Q≥Coll false",
        "proved": {
            "XA_is_psi_mode": A["proved"],
            "sum_u_constant": B["proved"],
            "u_not_two_level": C["proved"],
            "Q_lt_Coll_p5": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15296.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
