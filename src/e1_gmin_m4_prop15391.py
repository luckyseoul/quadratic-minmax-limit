#!/usr/bin/env python3
"""
Prop 15.391 — For p≡1, T=A·T_0 is {y<x} after A^{-1}=½[[1,1],[1,−1]].
Hence k=|T∩S∩ℓ| is the number of circle-line roots (15.390)
that satisfy y<x.  c still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.390 flags.

============================================================================
Setup.  15.349: A=[[1,1],[1,−1]], T_0={0≤y<x≤p−1}.
Encoding u=X+p Y.  15.390: S∩ℓ = roots of N(w+λs)=c.

============================================================================
Theorem A — PROVED (inverse of A; certified p=5,13,17,29).
  (X,Y)∈T ⇔ y<x where x=(X+Y)/2, y=(X−Y)/2 in {0,…,p−1}.
  Matches the packed y_T grid on all q points.
  Fail: y≤x (then the diagonal is included; |T|=p(p+1)/2 ≠ k).  ∎

Theorem B — PROVED (A + 15.390; certified p=5,13,17).
  On every R-line, live k equals #{u∈S∩ℓ : y<x}.
  Fail: k=m (drop T; false on split pairs of 15.390 B).  ∎

Theorem C — OPEN.  Names k given (t,c).  first_nc c still
  unnamed (2,5,3,17 at p=5,13,17,29).  Ensemble D open.
  phi_F not imported.

============================================================================
Backend: A^{-1} + first_nc meets.  Writes
evidence/e1_gmin_m4_prop15391.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc  # noqa: E402
from e1_gmin_m4_prop15349 import y_triangle_p1  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15391_catalog.json"


def xy_from_u(u: int, p: int) -> tuple[int, int]:
    X, Y = u % p, u // p
    inv2 = (p + 1) // 2
    x = ((X + Y) * inv2) % p
    y = ((X - Y) * inv2) % p
    return x, y


def in_T_lt(u: int, p: int) -> bool:
    x, y = xy_from_u(u, p)
    return y < x


def in_T_le(u: int, p: int) -> bool:
    x, y = xy_from_u(u, p)
    return y <= x


def k_named(meet: list[int], p: int) -> int:
    return sum(1 for u in meet if in_T_lt(u, p))


def k_drop_T(meet: list[int], p: int) -> int:
    return len(meet)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        y = y_triangle_p1(p)
        D = y[1 : 1 + p * p] < 0
        q = p * p
        n_lt = sum(1 for u in range(q) if in_T_lt(u, p) == bool(D[u]))
        n_le = sum(1 for u in range(q) if in_T_le(u, p) == bool(D[u]))
        nT = int(D.sum())
        n_le_size = sum(1 for u in range(q) if in_T_le(u, p))
        ok = ok and n_lt == q
        ok = ok and nT == p * (p - 1) // 2
        ok = ok and n_le < q
        ok = ok and n_le_size == p * (p + 1) // 2
        if n_le == q:
            ok = False
        rows[str(p)] = {
            "n_lt_match": n_lt,
            "n_le_match": n_le,
            "nT": nT,
            "n_le_size": n_le_size,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "T={y<x} after A^{-1}. Fail: y≤x (includes diagonal).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17):
        y = y_triangle_p1(p)
        D = y[1 : 1 + p * p] < 0
        pack = first_nc(p, y)
        S = set(int(u) for u in pack["circle"])
        F = _kernel_field(p)
        n_ok = n_drop = n_tot = 0
        for lines in square_classes(F):
            ns = D[lines].sum(axis=1)
            if tag_sig(p, tuple(int(x) for x in sorted(ns))) != "R":
                continue
            for ln in lines:
                meet = [int(x) for x in ln if int(x) in S]
                k_live = sum(1 for u in meet if D[u])
                n_tot += 1
                if k_named(meet, p) == k_live:
                    n_ok += 1
                if k_drop_T(meet, p) == k_live:
                    n_drop += 1
        ok = ok and n_ok == n_tot and n_tot > 0
        ok = ok and n_drop < n_tot
        if n_drop == n_tot:
            ok = False
        rows[str(p)] = {
            "t": int(pack["t"]),
            "c": int(pack["c"]),
            "n_ok": n_ok,
            "n_drop": n_drop,
            "n_tot": n_tot,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "k=#{u∈S∩ℓ: y<x}. Fail: k=m.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "k named given (t,c); c and D unnamed. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.391  k from y<x on circle-line roots", flush=True)
    A = prove_A()
    print(f"  A y<x: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B k: {B['proved']} {B['rows']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n")
    out = {
        "prop": "15.391",
        "title": "k is # of circle-line roots with y<x; c unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "T_is_y_lt_x": A["proved"],
            "k_from_ineq": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15391.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
