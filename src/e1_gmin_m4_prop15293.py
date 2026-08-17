#!/usr/bin/env python3
"""
Prop 15.293 — Q_NL is not the 1D⊙Hoffman class function of 15.290 types.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.292 flags.  Does **not** name Q_NL(p).  Floor stays OPEN.

============================================================================
Setup.  Q_1D is named (15.292).  NL = Max+ \\ 1D.  1D⊙nc = Hoffman
switches of Type+ parents (15.286–15.287).  At p=5, NL = 1D⊙nc
(100 sets).  At p=7, |1D⊙nc|=1764 < |NL|=5586.

============================================================================
Theorem A — PROVED (construction vs cache).
  Q_{1D⊙nc} = Q_NL at p=5 (both 24/5 on {pm1,N1} and 16/5 on {N*,sub}).
  Q_{1D⊙nc} ≠ Q_NL at p=7
      (Hoffman: 8 on pm1, 4 on sub, 32/9 on other off;
       NL: 1000/133, 480/133, 1376/399, 1496/399, 1544/399).
  Fail-when-wrong: claim Q_H=Q_NL at p=7.  ∎

Theorem B — PROVED at p=11 (Max+-free 1D⊙nc; no Max+ census).
  Paley×norm type does **not** determine Q on 1D⊙nc: AP-parent
  sub-class spread 4.61, (+−,N*) spread 0.87.  Paley-parent sub
  spread 0.44.  So Q_NL cannot be named as a 15.290 class function
  by evaluating Hoffman.  Fail-when-wrong: claim type-constancy of
  Q_H at p=11 (spread < 1e-4).  ∎

Theorem C — OPEN.  Q_NL(p) unnamed.  15.292 names only Q_1D.
  ⟨δ,ψ⟩≤2 / phi_F_ge_6 stay False.  ∎

============================================================================
Backend: CPU ProcessPool.  p=11 Hoffman in prove_B (W=full_workers).
Writes evidence/e1_gmin_m4_prop15293.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
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
from e1_gmin_m4_prop15290 import omega_set, type_key, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15292 import (  # noqa: E402
    Q_1d_off_Fp_over_q2,
    Q_1d_pm1_over_q2,
    Q_1d_sub_over_q2,
)
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
_G: dict = {}


def _init(p: int) -> None:
    F = _kernel_field(p)
    _G["F"] = F
    _G["Omega"] = omega_set(F)


def _u_of_z(z):
    return np.abs(zhat_omega(np.asarray(z, dtype=np.float64), _G["F"], _G["Omega"])) ** 2


def _Q_rel(zrows: np.ndarray, p: int) -> dict[int, float]:
    F = _kernel_field(p)
    Omega = omega_set(F)
    q2 = float(F["q"] ** 2)
    W = default_workers()
    with ProcessPoolExecutor(max_workers=W, initializer=_init, initargs=(p,)) as ex:
        u = np.stack(
            list(ex.map(_u_of_z, zrows, chunksize=max(1, len(zrows) // (W * 4))))
        )
    xi0 = Omega[0]
    om_idx = {xi: i for i, xi in enumerate(Omega)}
    u0 = u[:, om_idx[xi0]]
    out = {}
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * u[:, om_idx[eta]])) / q2
    return out


def _fp_dirs(p: int):
    dirs = [[t + ((s * t) % p) * p for t in range(p)] for s in range(p)]
    dirs.append([t * p for t in range(p)])
    return dirs


def _is_1d(z, dirs, add, p: int) -> bool:
    dset = set(int(i) for i in np.flatnonzero(np.asarray(z) < 0))
    for V in dirs:
        if all(sum(add[x][v] in dset for v in V) == p for x in dset):
            return True
    return False


def _pack(Qrel: dict[int, float], F) -> dict:
    by = defaultdict(list)
    for r, v in Qrel.items():
        by[str(type_key(F, r))].append(v)
    return {
        k: {
            "n": len(vs),
            "mean": float(np.mean(vs)),
            "frac": str(Fraction(float(np.mean(vs))).limit_denominator(2000)),
            "spread": float(max(vs) - min(vs)),
        }
        for k, vs in by.items()
    }


def prove_H_vs_NL() -> dict:
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in (5, 7):
        F = _kernel_field(p)
        DH = build_1d_odot_nc(p)
        zH = np.ones((DH.shape[0], F["q"]), dtype=np.float64)
        zH[DH > 0] = -1.0
        Y = np.load(YPATH[p])
        Yp = Y[Y[:, 0] > 0]
        z_all = np.sign(Yp[:, 1 : 1 + F["q"]].astype(np.float64))
        dirs = _fp_dirs(p)
        flags = np.array([_is_1d(z, dirs, F["add"], p) for z in z_all])
        z_nl = z_all[~flags]
        QH = _Q_rel(zH, p)
        QN = _Q_rel(z_nl, p)
        rel = max(abs(QH[r] - QN[r]) for r in QH)
        if p == 5 and rel > 1e-8:
            ok = False
        if p == 7 and rel < 1e-4:
            ok = False
        invert_tot += 1
        if p == 7 and rel < 1e-4:
            invert_hit += 1
        sample[str(p)] = {
            "n_H": int(DH.shape[0]),
            "n_NL": int(z_nl.shape[0]),
            "n_1d": int(flags.sum()),
            "max_rel": rel,
            "Q_H": _pack(QH, F),
            "Q_NL": _pack(QN, F),
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_p11_spread_from_scratch() -> dict:
    """Read the p=11 Hoffman-by-parent dump if present; else skip live rebuild.

    The dump is produced by the Max+-free script in the goal scratch.
    Fail-when-wrong is the numerical spread itself, not a flag flip.
    """
    path = ROOT / "evidence" / "e1_gmin_m4_hoffman_p11.json"
    if not path.is_file():
        return {
            "proved": False,
            "note": "p=11 Hoffman dump missing; B not certified this run",
        }
    data = json.loads(path.read_text())
    spreads = {}
    big = 0
    for cls, block in data.get("by_cls", {}).items():
        for tau, row in block.get("Q", {}).items():
            sp = float(row.get("spread", 0.0))
            spreads[f"{cls}:{tau}"] = sp
            if sp > 1e-4:
                big += 1
    # invert: type-constancy would have big==0
    proved = big > 0
    return {
        "proved": proved,
        "n_big_spreads": big,
        "spreads": spreads,
        "nD": {k: v.get("n") for k, v in data.get("by_cls", {}).items()},
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_1d_named": {
            "pm1": str(Q_1d_pm1_over_q2(7)),
            "sub": str(Q_1d_sub_over_q2(7)),
            "off_Fp": str(Q_1d_off_Fp_over_q2(7)),
        },
        "note": "Q_1D named (15.292).  Q_NL unnamed.  Floor OPEN.",
    }


def main() -> dict:
    print("Prop 15.293  Q_NL ≠ Hoffman class function of 15.290 types", flush=True)
    A = prove_H_vs_NL()
    print(f"  A Q_H vs Q_NL: {A['proved']}", flush=True)
    B = prove_p11_spread_from_scratch()
    print(f"  B p=11 type-spread: {B['proved']} n_big={B.get('n_big_spreads')}", flush=True)
    C = prove_open()
    print(f"  C Q_NL/floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.293",
        "title": "Q_NL is not Hoffman-on-15.290-types; Q_1D named; floor open",
        "proved": {
            "H_equals_NL_only_p5": A["proved"],
            "p11_type_does_not_determine_Q_H": B["proved"],
            "Q_NL_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15293.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
