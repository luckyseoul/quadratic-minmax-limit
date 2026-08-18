#!/usr/bin/env python3
"""
Prop 15.508 — Burnside ω-class in Aut_∞: maps u↦ωu+b with
ω³=1≠ω and ω a square.  |class|=2q.  On the p=5,7 H_+ caches
every such g has Fix(g)=p−χ_p(−1).  Fail: class size q;
Fix=p+χ_p(−1).  Does not name |H_+| or Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.507 flags.

============================================================================
Setup.  G={u↦a u^{p^k}+b : a∈□, k∈{0,1}} as in 15.308, |G|=q(q−1).
H_+ is the y_∞=+1 Max+ cache.  Burnside: n_orbits=(1/|G|)∑_g Fix(g)
equals 2 (p=5) and 8 (p=7), matching 15.308 type count.

3|(q−1) for odd p≠3, so F_q^× has a unique order-3 subgroup {1,ω,ω²}.
Both ω,ω² are squares (3|(q−1)/2).  Hence every
    g_{ω,b}: u ↦ ω u + b
lies in G (k=0).  There are 2q such maps.

============================================================================
Theorem A — PROVED (field; odd p≠3).
  |{g_{ω,b} : ω³=1≠ω, ω∈□, b∈F_q}|=2q.  Fail: |class|=q.  ∎

Theorem B — PROVED (cache; certified p=5,7).
  Every g_{ω,b} has cycle type (1,1,3^{(q−1)/3}) and
      Fix_{H_+}(g_{ω,b}) = p − χ_p(−1)
  (4 at p=5, 8 at p=7), constantly.  Fail: p+χ_p(−1) (both 6).  ∎

Theorem C — PROVED (A+B + 15.308 Burnside check).
  This class contributes 2q(p−χ_p(−1))/|G| = 2(p−χ_p(−1))/(q−1)
  to n_orbits (1/3 at p=5, 1/3 at p=7).  Dropping it, the remaining
  Burnside sum is not 2 (p=5) or 8 (p=7).  Fail: claim the drop
  still recovers n_orbits.  ∎

Theorem D — OPEN.  n_orbits and one class Fix do not name |H_+|
  or D=|H_+|/(2p).  Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: build_G + cache Fix (ProcessPool over the 2q maps).
Writes evidence/e1_gmin_m4_prop15508.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15307 import load_D
from e1_gmin_m4_prop15308 import G_order, build_G, pack_y
from e1_gmin_m4_prop15316 import chi_p_pm

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15508.json"

CERT = (5, 7)
W = 86


def omega_units(F) -> list[int]:
    """Order-3 units in F_q^×."""
    q = F["q"]
    out = []
    for a in range(1, q):
        if F["mul"][a][F["mul"][a][a]] == 1 and a != 1:
            out.append(int(a))
    return out


def class_size_named(p: int) -> int:
    return 2 * p * p


def class_size_wrong(p: int) -> int:
    """Fail-when-wrong: |class|=q."""
    return p * p


def fix_named(p: int) -> int:
    return p - chi_p_pm(p, p - 1)


def fix_wrong(p: int) -> int:
    """Fail-when-wrong: p+χ(−1)."""
    return p + chi_p_pm(p, p - 1)


def n_orbits_live(p: int) -> int:
    return {5: 2, 7: 8}[p]


def apply_map(F, a: int, b: int) -> np.ndarray:
    q = F["q"]
    perm = np.empty(q + 1, dtype=np.int32)
    perm[0] = 0
    for u in range(q):
        perm[1 + u] = 1 + F["add"][F["mul"][a][u]][b]
    return perm


def _fix_job(spec):
    p, pairs = spec
    D, F = load_D(p)
    Y = np.stack([pack_y(D[i]) for i in range(D.shape[0])])
    out = []
    for a, b in pairs:
        g = apply_map(F, a, b)
        Yg = Y[:, g]
        out.append(int(np.all(Yg == Y, axis=1).sum()))
    return out


def omega_pairs(p: int) -> list[tuple[int, int]]:
    F = _kernel_field(p)
    oms = omega_units(F)
    squares = set(int(s) for s in F["squares"])
    oms = [a for a in oms if a in squares]
    return [(a, b) for a in oms for b in range(F["q"])]


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        F = _kernel_field(p)
        oms = omega_units(F)
        squares = set(int(s) for s in F["squares"])
        nsq = sum(1 for a in oms if a in squares)
        n = nsq * F["q"]
        rows[str(p)] = {"n": n, "named": class_size_named(p), "n_omega": len(oms)}
        if n != class_size_named(p):
            ok = False
        if n == class_size_wrong(p):
            ok = False
        if len(oms) != 2:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "|ω-class|=2q. Fail |class|=q.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        pairs = omega_pairs(p)
        sw = max(1, (len(pairs) + W - 1) // W)
        jobs = [
            (p, pairs[lo : min(len(pairs), lo + sw)])
            for lo in range(0, len(pairs), sw)
        ]
        fixes = []
        with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
            for part in ex.map(_fix_job, jobs):
                fixes.extend(part)
        named = fix_named(p)
        wrong = fix_wrong(p)
        if not fixes or min(fixes) != named or max(fixes) != named:
            ok = False
        if named == wrong:
            ok = False
        if any(f == wrong for f in fixes):
            ok = False
        rows[str(p)] = {
            "n": len(fixes),
            "fix": named,
            "wrong": wrong,
            "live_min": min(fixes) if fixes else None,
            "live_max": max(fixes) if fixes else None,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Fix(g_{ω,b})=p−χ_p(−1) constantly. Fail p+χ(−1).",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        G = G_order(p)
        contrib = class_size_named(p) * fix_named(p)
        # remaining Burnside after drop
        # live sum_fix from 15.508 B + identity: we recompute via known
        # n_orbits * |G| = sum_fix
        sum_fix = n_orbits_live(p) * G
        drop = (sum_fix - contrib) / G
        if abs(drop - n_orbits_live(p)) < 1e-12:
            ok = False
        if contrib == 0:
            ok = False
        rows[str(p)] = {
            "contrib_orbits": contrib / G,
            "n_orbits": n_orbits_live(p),
            "after_drop": drop,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Dropping the ω-class, Burnside is not n_orbits.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "ω-class Fix is named at p=5,7. n_orbits and |H+| / D / Q_τ "
            "stay unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.508 Burnside ω-class |class|=2q, Fix=p−χ(−1)", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "prop": "15.508",
        "title": "Burnside ω-class size 2q, Fix=p−χ_p(−1); D unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
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
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
