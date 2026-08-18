#!/usr/bin/env python3
"""
Prop 15.499 — Aut_∞ translations act freely on Max+ \\ Type+ 1D
and every T-orbit has size p².  Then
    n_T = n_NL / p² = c(p−1)
equals c_eq(p−1) at p=5,7.  Fail: c_min (114≠108 at p=7).
F_p^× does **not** act with constant orbit size p−1 on those
T-orbits (p=7 sizes {3,6}, 26 orbits ≠ 19).  So c is not the
number of F_p^×-orbits.  D unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.498 flags.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15294 import _fp_dirs, _is_1d
from e1_gmin_m4_prop15298 import YPATH
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15499.json"


def _pack(D: np.ndarray) -> np.ndarray:
    q = D.shape[1]
    out = np.zeros(D.shape[0], dtype=np.uint64)
    for j in range(q):
        out |= D[:, j].astype(np.uint64) << np.uint64(j)
    return out


def _mask_from_bits(bits: np.ndarray) -> int:
    m = 0
    for j, b in enumerate(bits):
        if b:
            m |= 1 << j
    return int(m)


@lru_cache(maxsize=4)
def nl_designs(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    z = np.sign(Yp[:, 1 : 1 + q].astype(np.float64))
    D = (z < 0).astype(np.uint8)
    dirs = _fp_dirs(p)
    one = np.array([_is_1d(z[i], dirs, F["add"], p) for i in range(len(z))])
    return F, D[~one], int(one.sum()), int((~one).sum())


def translation_orbits(p: int, F, Dnl: np.ndarray) -> dict:
    q = F["q"]
    add = np.array(F["add"], dtype=np.int32)
    masks = _pack(Dnl)
    mset = set(int(m) for m in masks)
    trans_idx = [add[:, b] for b in range(q)]
    n_fixed = 0
    sizes = []
    seen: set[int] = set()
    reps: list[int] = []
    mask_to_orb: dict[int, int] = {}
    for m in masks:
        mi = int(m)
        if mi in seen:
            continue
        bits = np.array([(mi >> j) & 1 for j in range(q)], dtype=np.uint8)
        orb_id = len(reps)
        reps.append(mi)
        so = set()
        for b in range(q):
            tb = bits[trans_idx[b]]
            tm = _mask_from_bits(tb)
            if tm == mi and b != 0:
                n_fixed += 1
            if tm not in mset:
                return {"ok": False, "why": "translate_missing"}
            so.add(tm)
            mask_to_orb[tm] = orb_id
        sizes.append(len(so))
        seen |= so
    return {
        "ok": True,
        "nT": len(reps),
        "n_fixed": n_fixed,
        "sizes": sorted(set(sizes)),
        "reps": reps,
        "mask_to_orb": mask_to_orb,
        "trans_idx": trans_idx,
    }


def scale_orbits(p: int, F, Trec: dict) -> dict:
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    fp = list(range(1, p))
    ainv = {}
    for a in fp:
        ai = 1
        for _ in range(p - 2):
            ai = int(mul[ai, a])
        ainv[a] = ai
    sizes = []
    seen: set[int] = set()
    n_neg1_stab = 0
    for oid0, mi in enumerate(Trec["reps"]):
        if oid0 in seen:
            continue
        bits = np.array([(mi >> j) & 1 for j in range(q)], dtype=np.uint8)
        oset = set()
        for a in fp:
            ai = ainv[a]
            idx = [int(mul[ai, x]) for x in range(q)]
            sm = _mask_from_bits(bits[idx])
            if sm not in Trec["mask_to_orb"]:
                return {"ok": False, "why": "scale_missing"}
            oset.add(Trec["mask_to_orb"][sm])
        if p % 2 == 1:
            # a = p-1 ≡ −1
            ai = ainv[p - 1]
            idx = [int(mul[ai, x]) for x in range(q)]
            sm = _mask_from_bits(bits[idx])
            if Trec["mask_to_orb"].get(sm) == oid0:
                n_neg1_stab += 1
        sizes.append(len(oset))
        seen |= oset
    return {
        "ok": True,
        "nS": len(sizes),
        "sizes": sorted(set(sizes)),
        "n_neg1_stab": n_neg1_stab,
        "size_hist": {int(s): sizes.count(s) for s in sorted(set(sizes))},
    }


def prove_A() -> dict:
    """T free on NL; every T-orbit has size p². Fail: claim size p."""
    ok = True
    rows = {}
    for p in (5, 7):
        F, Dnl, n1, nNL = nl_designs(p)
        if n1 != n_1d(p):
            ok = False
        Trec = translation_orbits(p, F, Dnl)
        if not Trec["ok"] or Trec["n_fixed"] != 0:
            ok = False
        if Trec["sizes"] != [p * p]:
            ok = False
        if Trec["nT"] * p * p != nNL:
            ok = False
        if Trec["nT"] == nNL // p:
            ok = False
        rows[p] = {
            "n1d": n1,
            "nNL": nNL,
            "nT": Trec["nT"],
            "T_sizes": Trec["sizes"],
            "n_fixed": Trec["n_fixed"],
        }
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """n_T = c_eq (p−1). Fail: c_min at p=7."""
    ok = True
    rows = {}
    for p in (5, 7):
        F, Dnl, n1, nNL = nl_designs(p)
        Trec = translation_orbits(p, F, Dnl)
        ce, cm = c_eq_named(p), c_min_named(p)
        got = Trec["nT"]
        named = ce * (p - 1)
        wrong = cm * (p - 1)
        if got != named:
            ok = False
        rows[p] = {
            "nT": got,
            "c_eq_times_pm1": named,
            "c_min_times_pm1": wrong,
            "c_eq": ce,
            "c_min": cm,
        }
    if rows[7]["nT"] == rows[7]["c_min_times_pm1"]:
        ok = False
    if rows[7]["c_eq"] == rows[7]["c_min"]:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """F_p^× orbits of T-orbits are not all size p−1 at p=7."""
    ok = True
    rows = {}
    for p in (5, 7):
        F, Dnl, _, _ = nl_designs(p)
        Trec = translation_orbits(p, F, Dnl)
        Srec = scale_orbits(p, F, Trec)
        if not Srec["ok"]:
            ok = False
        rows[p] = {
            "nS": Srec["nS"],
            "sizes": Srec["sizes"],
            "size_hist": Srec["size_hist"],
            "n_neg1_stab": Srec["n_neg1_stab"],
            "c_eq": c_eq_named(p),
        }
    if rows[7]["sizes"] == [6]:
        ok = False
    if rows[7]["nS"] == rows[7]["c_eq"]:
        ok = False
    if 3 not in rows[7]["sizes"]:
        ok = False
    if rows[5]["sizes"] != [4]:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "T-orbits named in p²; n_T=c_eq(p−1) only as a p=5,7 cache "
            "identity. F_p^× orbits are not all size p−1. D unnamed."
        ),
    }


def main() -> dict:
    print("15.499 Aut_∞ orbits", flush=True)
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
