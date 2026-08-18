#!/usr/bin/env python3
"""
Prop 15.533 — Aut_∞ involutions do not freely pair T_ns-orbits.
n_inv(x↦−x)=Fix(g_0)=2·3^{C((p−1)/2,2)} (15.509).  Frob and
Frob∘(−id) also have invariant orbits.  Fail: n_inv=0.
Does not name D.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.532 flags.  Does **not** overwrite 15.522–15.532.

============================================================================
Setup.  15.532: T_v is free on H_+ iff χ_q(v)=−1, so
|H_+|=p n_orb(T_ns) and n_orb=2D live.  A free involution
normalising ⟨T_v⟩ with 0 invariant orbits would give
n_orb=2 n_orb(⟨T_v,σ⟩) and D=n_orb(⟨T_v,σ⟩).  x↦−x and
x↦−x+b normalise every line-group (conjugate T_v to T_{−v}).
Frob does too (v^p=N(v)v).

============================================================================
Theorem A — PROVED (H_+ caches; 15.509).
  On T_ns-orbits, n_inv(x↦−x)=n_inv(x↦−x+1)
  =2·3^{C((p−1)/2,2)} = Fix(g_b) (6 at p=5, 54 at p=7).
  Fail: n_inv=0; fail n_inv=n_1d/p at p=7 (20≠54).  ∎

Theorem B — PROVED (same caches).
  Frob has n_inv=2,218.  Frob∘(−id) has n_inv=22,218.
  Fail: either is 0.  ∎

Theorem C — OPEN.  No tested Aut_∞ involution is a free
  pairing, so D=n_orb/2 is not a named orbit count.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: ProcessPool over (p, involution).  Writes
evidence/e1_gmin_m4_prop15533.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15526 import load_Hplus
from e1_gmin_m4_prop15532 import DIRS

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15533.json"

W = 58
NS_DIR = {5: DIRS[5]["ns"], 7: DIRS[7]["ns"]}


def n_inv_neg_named(p: int) -> int:
    """15.509 Fix(g_b)=2·3^{C((p−1)/2,2)}."""
    return 2 * 3 ** comb((p - 1) // 2, 2)


def n_inv_wrong_zero() -> int:
    return 0


def n_inv_wrong_1d_orb(p: int) -> int:
    return n_1d(p) // p


def _key(y: np.ndarray) -> bytes:
    return np.ascontiguousarray(y.astype(np.int8)).tobytes()


def _tv_perm(p: int, da: int, db: int) -> np.ndarray:
    q = p * p
    perm = np.empty(q + 1, dtype=np.int32)
    perm[0] = 0
    for x in range(q):
        perm[1 + x] = 1 + ((x % p + da) % p) + ((x // p + db) % p) * p
    return perm


def _pow_el(F, x: int, e: int) -> int:
    if x == 0:
        return 0
    r, b = 1, int(x)
    while e:
        if e & 1:
            r = F["mul"][r][b]
        b = F["mul"][b][b]
        e >>= 1
    return r


def _sigma(F, kind: str) -> np.ndarray:
    q, p = F["q"], F["p"]
    sig = np.empty(q + 1, dtype=np.int32)
    sig[0] = 0
    if kind == "neg":
        for u in range(q):
            sig[1 + u] = 1 + F["neg"][u]
        return sig
    if kind == "g1":
        for u in range(q):
            sig[1 + u] = 1 + F["add"][F["neg"][u]][1]
        return sig
    if kind == "frob":
        for u in range(q):
            sig[1 + u] = 1 + _pow_el(F, u, p)
        return sig
    if kind == "frob_neg":
        for u in range(q):
            sig[1 + u] = 1 + _pow_el(F, F["neg"][u], p)
        return sig
    raise KeyError(kind)


def _apply(y: np.ndarray, perm: np.ndarray) -> np.ndarray:
    z = np.empty_like(y)
    z[0] = y[0]
    z[1:] = y[perm[1:]]
    return z


def _inv_job(spec) -> dict:
    p, kind = spec
    F = _kernel_field(p)
    H = load_Hplus(p)
    da, db = NS_DIR[p]
    tv = _tv_perm(p, da, db)
    sigma = _sigma(F, kind)
    table = {_key(H[i]): i for i in range(len(H))}
    oid = -np.ones(len(H), dtype=np.int32)
    n_orb = 0
    for i in range(len(H)):
        if oid[i] >= 0:
            continue
        j = i
        while oid[j] < 0:
            oid[j] = n_orb
            j = table[_key(_apply(H[j], tv))]
        n_orb += 1
    inv = set()
    for i in range(len(H)):
        j = table[_key(_apply(H[i], sigma))]
        if oid[i] == oid[j]:
            inv.add(int(oid[i]))
    return {
        "p": p,
        "kind": kind,
        "n_orb": n_orb,
        "n_inv": len(inv),
    }


def prove_A() -> dict:
    specs = [(p, k) for p in (5, 7) for k in ("neg", "g1")]
    with ProcessPoolExecutor(max_workers=min(W, len(specs))) as ex:
        rows = list(ex.map(_inv_job, specs))
    by = {(r["p"], r["kind"]): r for r in rows}
    ok = True
    rec = {}
    for p in (5, 7):
        named = n_inv_neg_named(p)
        n_neg = by[(p, "neg")]["n_inv"]
        n_g1 = by[(p, "g1")]["n_inv"]
        row_ok = (
            n_neg == named == n_g1
            and n_neg != n_inv_wrong_zero()
            and (p != 7 or n_neg != n_inv_wrong_1d_orb(p))
        )
        if not row_ok:
            ok = False
        rec[str(p)] = {
            "n_inv_neg": n_neg,
            "n_inv_g1": n_g1,
            "named": named,
            "n_1d_over_p": n_inv_wrong_1d_orb(p),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "n_inv(−id)=n_inv(g_1)=2·3^{C((p−1)/2,2)}. "
            "Fail 0; fail n_1d/p at p=7."
        ),
    }


def prove_B() -> dict:
    specs = [(p, k) for p in (5, 7) for k in ("frob", "frob_neg")]
    with ProcessPoolExecutor(max_workers=min(W, len(specs))) as ex:
        rows = list(ex.map(_inv_job, specs))
    by = {(r["p"], r["kind"]): r for r in rows}
    ok = all(by[(p, k)]["n_inv"] != 0 for p in (5, 7) for k in ("frob", "frob_neg"))
    rec = {
        str(p): {
            "frob": by[(p, "frob")]["n_inv"],
            "frob_neg": by[(p, "frob_neg")]["n_inv"],
        }
        for p in (5, 7)
    }
    if rec["5"]["frob"] != 2 or rec["7"]["frob"] != 218:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": "Frob n_inv=2,218; Frob∘(−id)=22,218. Fail 0.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "D_named_as_orbits": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "No tested Aut_∞ involution freely pairs T_ns-orbits. "
            "D and Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.533 Aut_∞ involutions do not freely pair T_ns-orbits", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.533",
        "title": "n_inv(−id)=15.509 Fix; no free Aut_∞ pairing of T_ns-orbits",
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
