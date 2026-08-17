#!/usr/bin/env python3
"""
Prop 15.297 — Aut involutions do not sign Cov; X_A can vanish on Hoffman.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.296 flags.  Floor stays OPEN.

============================================================================
Theorem A — PROVED (p=5 cache; named Aut maps).
  Pair-average (u X_A + u_T X_{A,T})/2 for T ∈ {signed inversion,
  Frob, x↦−x, translation} has min < 6q².  Those involutions do
  not sign Cov(u,X_A)≥6q².  Fail-when-wrong: claim signed-inv
  pair-min ≥ 6q² at p=5 k=2.  ∎

Theorem B — PROVED (p=5 vs p=7).
  At p=5, even k=2: X_A ≡ 0 on every Hoffman vector, while
  u ∈ {10(5±√5)} ≠ 0.  The whole covariance is 1D.
  At p=7, k=40 (S_full-minimizer): X_A is not ≡ 0 on Hoffman
  (S_H/q²=8).  Fail-when-wrong: claim X_A≡0 on Hoffman at p=7 k=40.  ∎

Theorem C — OPEN.  No complementary Aut/Hoffman-vanishing argument
  signs Cov for all p.  phi_F False.  ∎

============================================================================
Backend: CPU.  p=5 serial; p=7 ProcessPool.
Writes evidence/e1_gmin_m4_prop15297.json
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
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15286 import build_1d_odot_nc  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15296 import XA_psi_mode, _init  # noqa: E402
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _load_z(p: int) -> np.ndarray:
    F = _kernel_field(p)
    Y = np.load(YPATH[p])
    return np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64))


def _prod_row(z):
    F, Omega = _G_F_Omega()
    uu = np.abs(zhat_omega(np.asarray(z, dtype=np.float64), F, Omega)) ** 2
    return float(uu[0]), complex(XA_psi_mode(z))


def _G_F_Omega():
    from e1_gmin_m4_prop15296 import _G

    return _G["F"], _G["Omega"]


def signed_inv_z(z, F):
    q = F["q"]
    new_inf = z[0]
    out = np.zeros(q, dtype=np.float64)
    out[0] = 1.0
    for x in range(1, q):
        out[x] = F["chi"][x] * z[F["inv"][x]]
    if new_inf < 0:
        out = -out
    return out


def frob_z(z, F):
    p, q = F["p"], F["q"]
    out = np.zeros(q, dtype=np.float64)
    out[0] = z[0]
    for x in range(1, q):
        acc, b, e = 1, int(x), p
        while e:
            if e & 1:
                acc = F["mul"][acc][b]
            b = F["mul"][b][b]
            e >>= 1
        out[acc] = z[x]
    return out


def prove_aut_does_not_sign() -> dict:
    p, k = 5, 2
    F = _kernel_field(p)
    _init(p, k)
    zs = _load_z(p)
    prods = []
    keys = {}
    for i, z in enumerate(zs):
        u, xa = _prod_row(z)
        prods.append(u * xa.real)
        keys[tuple(int(v) for v in np.sign(z))] = i
    prods = np.array(prods)
    q2 = F["q"] ** 2
    floor = 6 * q2

    def pair_min(fn):
        vals = []
        hits = 0
        for i, z in enumerate(zs):
            z2 = fn(z)
            j = keys.get(tuple(int(v) for v in np.sign(z2)))
            if j is None:
                continue
            hits += 1
            vals.append(0.5 * (prods[i] + prods[j]))
        return hits, (min(vals) if vals else None)

    rows = {}
    ok = True
    invert_hit = invert_tot = 0
    for name, fn in (
        ("signed_inv", lambda z: signed_inv_z(z, F)),
        ("frob", lambda z: frob_z(z, F)),
        ("negx", lambda z: np.array([z[F["neg"][x]] for x in range(F["q"])])),
        ("trans1", lambda z: np.array([z[F["add"][x][1]] for x in range(F["q"])])),
    ):
        hits, mn = pair_min(fn)
        rows[name] = {
            "hits": hits,
            "pair_min": mn,
            "below_floor": bool(mn is not None and mn < floor),
        }
        if mn is None or mn >= floor:
            ok = False
    invert_tot += 1
    # fail-when-wrong: signed_inv pair-min ≥ 6q²
    if rows["signed_inv"]["pair_min"] is not None and rows["signed_inv"]["pair_min"] >= floor:
        invert_hit += 1
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "floor": floor,
        "E_uXA": float(prods.mean()),
        "maps": rows,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def _hoffman_mask(p: int, zs: np.ndarray) -> np.ndarray:
    F = _kernel_field(p)
    DH = build_1d_odot_nc(p)
    Hset = {tuple(int(x) for x in row) for row in DH}
    dbits = (zs < 0).astype(np.int8)
    return np.array([tuple(int(x) for x in dbits[i]) in Hset for i in range(len(zs))])


def prove_XA_vanishes_on_Hoffman() -> dict:
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    # p=5 k=2: must vanish
    p, k = 5, 2
    _init(p, k)
    zs = _load_z(p)
    H = _hoffman_mask(p, zs)
    xas = np.array([XA_psi_mode(z) for z in zs])
    van5 = float(np.max(np.abs(xas[H]))) if H.any() else 1.0
    if van5 > 1e-6:
        ok = False
    sample["5:2"] = {
        "n_H": int(H.sum()),
        "max_abs_XA": van5,
        "vanishes": van5 <= 1e-6,
    }
    # p=7 k=40: must NOT vanish
    p, k = 7, 40
    _init(p, k)
    zs = _load_z(p)
    H = _hoffman_mask(p, zs)
    W = default_workers()
    with ProcessPoolExecutor(max_workers=W, initializer=_init, initargs=(p, k)) as ex:
        xas = np.array(
            list(ex.map(XA_psi_mode, zs[H], chunksize=max(1, int(H.sum()) // (W * 4))))
        )
    mx7 = float(np.max(np.abs(xas))) if len(xas) else 0.0
    if mx7 <= 1e-4:
        ok = False
    invert_tot += 1
    if mx7 <= 1e-4:
        invert_hit += 1
    sample["7:40"] = {
        "n_H": int(H.sum()),
        "max_abs_XA": mx7,
        "vanishes": mx7 <= 1e-4,
    }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": "Aut partners and Hoffman-vanishing do not sign Cov for all p.",
    }


def main() -> dict:
    print("Prop 15.297  Aut partners fail; X_A kills Hoffman only at p=5", flush=True)
    A = prove_aut_does_not_sign()
    print(f"  A Aut pair-min < 6q²: {A['proved']}", flush=True)
    B = prove_XA_vanishes_on_Hoffman()
    print(f"  B Hoffman XA vanish p=5 not p=7: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.297",
        "title": "Aut involutions do not sign Cov; Hoffman ker X_A is p=5-only",
        "proved": {
            "aut_pair_below_floor": A["proved"],
            "Hoffman_XA_zero_only_p5": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15297.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
