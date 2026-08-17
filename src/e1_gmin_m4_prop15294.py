#!/usr/bin/env python3
"""
Prop 15.294 — Q on rest = NL \\ Hoffman; familywise S≥6 is false.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.293 flags.  Does **not** name Q_rest as a formula in p
(only p=7 has rest).  Floor stays OPEN.

============================================================================
Setup.  Max+ (y_∞=+1) = 1D ⊔ Hoffman ⊔ rest.  Hoffman = 1D⊙nc
(15.286).  rest = NL \\ Hoffman = ync⊙nc double Seidel (15.139)
at p=7 (|rest|=3822=13·294).  Empty at p=5.

============================================================================
Theorem A — PROVED (construction vs cache).
  rest is empty at p=5.  At p=7 rest is nonempty, Q_rest is constant
  on 15.290 types, and Q_rest ≠ Q_H:
      pm1 664/91,  sub 24/7,
      (--,N1) 928/273, (++ ,N1) 1096/273, N* 1048/273.
  Fail-when-wrong: rest empty at p=7; Q_rest=Q_H at p=7.  ∎

Theorem B — PROVED (same partition; even ψ, k even, k∉{0,(q−1)/2}).
  Familywise S_□≥6q² is FALSE:
      S_H/q² hits 0 at p=5 (k=16);
      S_rest/q² hits 4.923…=64/13 at p=7.
  S_1D/q² ≥ 4p(p−1)/(p−2) > 6 (15.292).  S_full/q² ≥6 at p=5,7
  (15.290 B) only as a mixture: S_full=(n_1D S_1D+n_H S_H+n_rest S_rest)/N+.
  N+ unnamed, so the mixture does not import the floor.
  Fail-when-wrong: claim min S_H/q²≥6 at p=5.  ∎

Theorem C — OPEN.  Q_rest(p) unnamed (rest empty at p=5; one prime
  cannot interpolate).  ⟨δ,ψ⟩≤2 / phi_F stay False.  ∎

============================================================================
Backend: CPU ProcessPool W=full_workers over Max+ rows.
Writes evidence/e1_gmin_m4_prop15294.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15286 import build_1d_odot_nc  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, type_key, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15292 import (  # noqa: E402
    S_1d_over_q2_ntriv_on_Fp,
    n_1d,
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


def _Q_from_u(u: np.ndarray, F, Omega) -> dict[int, float]:
    q2 = float(F["q"] ** 2)
    xi0 = Omega[0]
    om_idx = {xi: i for i, xi in enumerate(Omega)}
    u0 = u[:, om_idx[xi0]]
    out = {}
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * u[:, om_idx[eta]])) / q2
    return out


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


def even_ks(p: int) -> list[int]:
    q = p * p
    chi_k = (q - 1) // 2
    return [k for k in range(0, q - 1) if k % 2 == 0 and k not in (0, chi_k)]


def S_over_q2(Qrel: dict[int, float], F, k: int) -> complex:
    psi = _psi_vec(F, k)
    return sum(psi[r] * val for r, val in Qrel.items())


@lru_cache(maxsize=4)
def partition_cache(p: int):
    F = _kernel_field(p)
    Omega = omega_set(F)
    DH = build_1d_odot_nc(p)
    Hset = {tuple(int(x) for x in row) for row in DH}
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    z_all = np.sign(Yp[:, 1 : 1 + F["q"]].astype(np.float64))
    dirs = _fp_dirs(p)
    flags_1d = np.array([_is_1d(z, dirs, F["add"], p) for z in z_all])
    d_bits = (z_all < 0).astype(np.int8)
    in_H = np.array(
        [tuple(int(x) for x in d_bits[i]) in Hset for i in range(len(z_all))]
    )
    W = default_workers()
    with ProcessPoolExecutor(max_workers=W, initializer=_init, initargs=(p,)) as ex:
        u_all = np.stack(
            list(ex.map(_u_of_z, z_all, chunksize=max(1, len(z_all) // (W * 4))))
        )
    idx = {
        "1d": np.flatnonzero(flags_1d),
        "H": np.flatnonzero(in_H),
        "rest": np.flatnonzero((~flags_1d) & (~in_H)),
    }
    Qs = {}
    for name, ix in idx.items():
        if len(ix) == 0:
            Qs[name] = {}
        else:
            Qs[name] = _Q_from_u(u_all[ix], F, Omega)
    Qs["full"] = _Q_from_u(u_all, F, Omega)
    return F, Omega, idx, Qs, int(DH.shape[0])


def prove_rest_Q() -> dict:
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in (5, 7):
        F, _Omega, idx, Qs, nHset = partition_cache(p)
        n_rest = int(len(idx["rest"]))
        n_H = int(len(idx["H"]))
        if p == 5 and n_rest != 0:
            ok = False
        if p == 7 and n_rest == 0:
            ok = False
        rel = None
        if Qs["rest"] and Qs["H"]:
            rel = max(abs(Qs["rest"][r] - Qs["H"][r]) for r in Qs["H"])
            if rel < 1e-4:
                ok = False
        invert_tot += 1
        if p == 7 and (n_rest == 0 or (rel is not None and rel < 1e-4)):
            invert_hit += 1
        # type-constancy of rest at p=7
        packed_rest = _pack(Qs["rest"], F) if Qs["rest"] else {}
        if p == 7:
            if not packed_rest or max(v["spread"] for v in packed_rest.values()) > 1e-8:
                ok = False
        sample[str(p)] = {
            "n_1d": int(len(idx["1d"])),
            "n_1d_formula": n_1d(p),
            "n_H": n_H,
            "n_Hset": nHset,
            "n_rest": n_rest,
            "max_rel_rest_H": rel,
            "Q_rest": packed_rest,
            "Q_H": _pack(Qs["H"], F) if Qs["H"] else {},
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_familywise_floor_false() -> dict:
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in (5, 7):
        F, _Omega, idx, Qs, _nH = partition_cache(p)
        mins = {}
        worst = {}
        for name, Q in Qs.items():
            if not Q:
                mins[name] = None
                continue
            vals = [(k, float(S_over_q2(Q, F, k).real)) for k in even_ks(p)]
            mn_k, mn = min(vals, key=lambda t: t[1])
            mins[name] = mn
            worst[name] = {"k": mn_k, "S_over_q2": mn}
        # B: familywise ≥6 is false
        if mins.get("H") is None or mins["H"] >= 6 - 1e-9:
            if p == 5:
                ok = False
        if p == 7:
            if mins.get("rest") is None or mins["rest"] >= 6 - 1e-9:
                ok = False
        if mins.get("1d") is None or mins["1d"] < 6:
            ok = False
        if mins.get("full") is None or mins["full"] < 6:
            ok = False
        invert_tot += 1
        if p == 5 and mins.get("H") is not None and mins["H"] >= 6 - 1e-9:
            invert_hit += 1
        sample[str(p)] = {
            "Smin": mins,
            "worst": worst,
            "S_1d_ntriv_formula": str(S_1d_over_q2_ntriv_on_Fp(p)),
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "theorem": (
            "S_H/q²=0 at p=5 and S_rest/q²<6 at p=7, so familywise "
            "S≥6 cannot import the floor.  Mixture needs N+."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": (
            "Q_rest type-constant at p=7, unnamed in p.  "
            "Familywise S≥6 is false.  Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.294  Q_rest type-constant at p=7; familywise S≥6 false", flush=True)
    A = prove_rest_Q()
    print(f"  A rest Q: {A['proved']}", flush=True)
    B = prove_familywise_floor_false()
    print(f"  B familywise floor false: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C general floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.294",
        "title": "Q_rest at p=7; familywise S≥6 false; floor open",
        "proved": {
            "rest_Q_type_constant_p7": A["proved"],
            "familywise_S_ge_6_false": B["proved"],
            "Q_rest_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15294.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
