#!/usr/bin/env python3
"""
Prop 15.512 — Frobenius ω-class: u↦ω u^p+b.

  Same multipliers as 15.508, but k=1 (Frobenius).  |class|=2q.
  On the p=5,7 H_+ caches every such g has Fix(g)=2 constantly.
  Fail: |class|=q; Fix=p−χ_p(−1) (the 15.508 k=0 value 4,8).
  Does not name |H_+| or Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.511 flags.

============================================================================
Setup.  15.508 names g_{ω,b}: u↦ωu+b (k=0).  The twin class is
    g_{ω,b}^{(p)}: u ↦ ω u^p + b
with the same ω³=1≠ω, ω∈□.  There are 2q such maps.

============================================================================
Theorem A — PROVED (field; odd p≠3).
  |{g_{ω,b}^{(p)}}|=2q.  Fail: |class|=q.  ∎

Theorem B — PROVED (cache; certified p=5,7).
  Every g_{ω,b}^{(p)} has Fix_{H_+}=2 constantly.
  Fail: p−χ_p(−1) (4 at p=5, 8 at p=7).  ∎

Theorem C — OPEN.  Cycle type is not unique at p=5 (two types).
  Remaining Frob Fix (high-order / N≠1) unnamed.
  D / Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: frob_table + cache Fix (ProcessPool over the 2q maps).
Writes evidence/e1_gmin_m4_prop15512.json
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
from e1_gmin_m4_prop15308 import frob_table, pack_y
from e1_gmin_m4_prop15508 import (
    class_size_named,
    class_size_wrong,
    fix_named as omega_k0_fix,
    omega_units,
)

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15512.json"

CERT = (5, 7)
W = 86


def fix_named(_p: int) -> int:
    return 2


def fix_wrong(p: int) -> int:
    """Fail-when-wrong: 15.508 k=0 Fix p−χ(−1)."""
    return omega_k0_fix(p)


def apply_frob(F, a: int, b: int, fr) -> np.ndarray:
    q = F["q"]
    mul, add = F["mul"], F["add"]
    perm = np.empty(q + 1, dtype=np.int32)
    perm[0] = 0
    for u in range(q):
        w = int(fr[u])
        perm[1 + u] = 1 + add[mul[a][w]][b]
    return perm


def omega_frob_pairs(p: int) -> list[tuple[int, int]]:
    F = _kernel_field(p)
    squares = set(int(s) for s in F["squares"])
    oms = [a for a in omega_units(F) if a in squares]
    return [(a, b) for a in oms for b in range(F["q"])]


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        n = len(omega_frob_pairs(p))
        rows[str(p)] = {"n": n, "named": class_size_named(p)}
        if n != class_size_named(p):
            ok = False
        if n == class_size_wrong(p):
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "|ω-Frob class|=2q. Fail |class|=q.",
    }


def _fix_job(spec):
    p, pairs = spec
    D, F = load_D(p)
    Y = np.stack([pack_y(D[i]) for i in range(D.shape[0])])
    fr = frob_table(F)
    out = []
    for a, b in pairs:
        g = apply_frob(F, a, b, fr)
        out.append(int(np.all(Y[:, g] == Y, axis=1).sum()))
    return out


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        pairs = omega_frob_pairs(p)
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
        "theorem": "Fix(g_{ω,b}^{(p)})=2 constantly. Fail p−χ(−1).",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "ω-Frob Fix is 2 at p=5,7. Remaining Frob Fix unnamed. "
            "D / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.512 Frobenius ω-class |class|=2q, Fix=2", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.512",
        "title": "Frobenius ω-class size 2q, Fix=2; D unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": D,
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
