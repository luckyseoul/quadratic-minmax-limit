#!/usr/bin/env python3
"""
Prop 15.485 — Aut_∞ orbit types at p=7 are the AG direction
classes: 2×1D (stabs 4p,6p), 4×s=2 (3 NF1 + 1 TR), 1×s=4 and
1×s=8 (OTHER0).  c=19.  Not all s=2 orbits are TR.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.484 flags.

============================================================================
Setup.  15.331: |H_+|=n_1d + c p²(p−1), c=Σ_τ (p+1)/s_τ, s_τ|(p+1)
on NL.  15.476–15.483: TR / NF1 / OTHER0 / PAR.  Live c(7)=19.

============================================================================
Theorem A — PROVED (orbit–stab fold; p=7 cache).
  Aut_∞ has 8 orbits on H_+: two 1D (sizes 84,56; stabs 28=4p
  and 42=6p), four NL of stab 2 and size 1176, one of stab 4
  size 588, one of stab 8 size 294.  Σ_NL (p+1)/s = 19.
  Fail: claim a single 1D orbit of size n_1d=140.  ∎

Theorem B — PROVED (A + 15.476 types).
  Exactly one of the four s=2 orbits is TR; the other three
  are NF1.  s=4 and s=8 are OTHER0.  Fail: claim every s=2
  orbit is TR.  ∎

Theorem C — OPEN.  Multiplicities (4,1,1) of (s=2,4,8) unnamed
  in p.  p=11 AGL(T)∩Max+ is TR-only and gives c_NL≥20≠c_eq=198.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.308 Aut_∞ + 15.476 types on the p=7 cache.  Writes
evidence/e1_gmin_m4_prop15485.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15308 import build_G, G_order  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15476 import is_triple_rainbow, live_delta2_and_sigs  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15485_catalog.json"


def geom_tag(sig, nf, p: int) -> str:
    if is_triple_rainbow(sig, p):
        return "TR"
    if int(nf) >= 2:
        return "PAR"
    if int(nf) == 1:
        return "NF1"
    return "OTHER0"


@lru_cache(None)
def fold_p7() -> dict:
    p = 7
    F = _kernel_field(p)
    q = F["q"]
    G = build_G(F, use_frob=True)
    Ds = D_rows(p, F).astype(np.int8)
    _d2, sigs, n_full = live_delta2_and_sigs(p)
    lookup = {Ds[i].tobytes(): i for i in range(Ds.shape[0])}
    remaining = set(range(Ds.shape[0]))
    recs = []
    c_nl = 0
    while remaining:
        i = remaining.pop()
        Dset = np.where(Ds[i] > 0)[0]
        members = set()
        for g in G:
            img = np.zeros(q, dtype=np.int8)
            img[g[1 + Dset] - 1] = 1
            j = lookup[img.tobytes()]
            members.add(j)
            remaining.discard(j)
        tags = Counter(geom_tag(sigs[j], n_full[j], p) for j in members)
        sz = len(members)
        s = len(G) // sz
        is_1d = tags.get("PAR", 0) == sz
        piece = 0
        if (not is_1d) and (p + 1) % s == 0:
            piece = (p + 1) // s
            c_nl += piece
        recs.append(
            {
                "sz": sz,
                "s": s,
                "tags": dict(tags),
                "is_1d": is_1d,
                "c_piece": piece,
            }
        )
    return {
        "n_orbits": len(recs),
        "c_nl": c_nl,
        "n_1d": n_1d(p),
        "G": len(G),
        "G_named": G_order(p),
        "recs": recs,
    }


def prove_A() -> dict:
    f = fold_p7()
    recs = f["recs"]
    ones = [r for r in recs if r["is_1d"]]
    nl = [r for r in recs if not r["is_1d"]]
    ok = f["n_orbits"] == 8 and f["c_nl"] == 19 and f["G"] == f["G_named"]
    ok = ok and len(ones) == 2
    ok = ok and sorted(r["sz"] for r in ones) == [56, 84]
    ok = ok and sorted(r["s"] for r in ones) == [28, 42]
    ok = ok and sum(r["sz"] for r in ones) == n_1d(7) == 140
    s2 = [r for r in nl if r["s"] == 2]
    ok = ok and len(s2) == 4 and all(r["sz"] == 1176 for r in s2)
    ok = ok and any(r["s"] == 4 and r["sz"] == 588 for r in nl)
    ok = ok and any(r["s"] == 8 and r["sz"] == 294 for r in nl)
    if len(ones) == 1 and ones[0]["sz"] == 140:
        ok = False
    return {
        "proved": bool(ok),
        "n_orbits": f["n_orbits"],
        "c_nl": f["c_nl"],
        "ones": [{"sz": r["sz"], "s": r["s"]} for r in ones],
        "nl_s": sorted(r["s"] for r in nl),
        "theorem": (
            "8 Aut_∞ orbits; 1D splits 84+56; NL stabs 2⁴·4·8; c=19. "
            "Fail: one 1D orbit of size 140."
        ),
    }


def prove_B() -> dict:
    f = fold_p7()
    s2 = [r for r in f["recs"] if r["s"] == 2]
    n_tr = sum(1 for r in s2 if r["tags"].get("TR", 0) == r["sz"])
    n_nf1 = sum(1 for r in s2 if r["tags"].get("NF1", 0) == r["sz"])
    oth = [r for r in f["recs"] if r["s"] in (4, 8)]
    ok = n_tr == 1 and n_nf1 == 3 and len(s2) == 4
    ok = ok and all(r["tags"].get("OTHER0", 0) == r["sz"] for r in oth)
    if n_tr == 4:
        ok = False
    return {
        "proved": bool(ok),
        "n_s2": len(s2),
        "n_s2_TR": n_tr,
        "n_s2_NF1": n_nf1,
        "theorem": (
            "Exactly one s=2 orbit is TR; three are NF1. "
            "Fail: every s=2 orbit is TR."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Multiplicities of NL stabs unnamed in p. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.485  Aut_∞ orbits are AG types; not all s=2 are TR", flush=True)
    A = prove_A()
    print(f"  A orbits: {A['proved']} n={A['n_orbits']} c={A['c_nl']}", flush=True)
    B = prove_B()
    print(f"  B TR/NF1: {B['proved']} TR={B['n_s2_TR']} NF1={B['n_s2_NF1']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.485",
        "title": "Aut_∞ orbits at p=7 are AG types; not all s=2 are TR",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "orbit_table_p7": A["proved"],
            "one_s2_is_TR": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15485.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
