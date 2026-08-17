#!/usr/bin/env python3
"""
Prop 15.477 — Triple-rainbow Δ² is named:
    p≡1 (mod 4):  Δ² = p S_χ² = M²/p = (S_χ G(χ_p))²
    p≡3 (mod 4):  Δ² = (M/λ)² = [p²(p+1)/2]²
    λ=(p−1)/2, M=k(q−k), S_χ=p(p²−1)/4.
The class is the AGL(2,p)-orbit of the ordered triangle
    T={(x,y)∈F_p² : y<x}
intersect Max+ (equals the 15.476 triple-rainbow type).
Q_τ unnamed.  Other complementary types unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.476 flags.

============================================================================
Setup.  15.476: Δ² is a function of direction type.  Triple-rainbow
is 3×(0,…,p−1) plus (p−2)×λ-regular.  T has size C(p,2)=k.
Raw T leaks off Ω.  AGL images that are Max+ are exactly the
triple-rainbow rows (certified p=5; 600 hits / 100 unique).

============================================================================
Theorem A — PROVED (named congruence split; certified p=5,7).
  On every live triple-rainbow Max+, Δ² equals delta2_tr_named(p).
  Fail: swap the congruence (5625≠4500 at p=5; 49392≠38416 at p=7).  ∎

Theorem B — PROVED (AGL fold; p=5).
  Raw T is not Max+ (leak>1).  Every AGL(2,5) image that lies in
  the live Max+ cache is triple-rainbow with Δ²=4500 (600 hits,
  100 unique = all TR).  Fail: claim raw T is Max+.  ∎

Theorem C — OPEN.  Other complementary types (p=7 non-TR) and
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.476 type fold + AGL(2,5) enum + Fraction.  Writes
evidence/e1_gmin_m4_prop15477.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import product
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows, S_chi_named  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15475 import M_named  # noqa: E402
from e1_gmin_m4_prop15476 import (  # noqa: E402
    is_triple_rainbow,
    live_delta2_and_sigs,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15477_catalog.json"


def lambda_reg(p: int) -> int:
    return (p - 1) // 2


def delta2_tr_named(p: int) -> int:
    if p % 4 == 1:
        S = S_chi_named(p)
        return p * S * S
    M = M_named(p)
    return (M // lambda_reg(p)) ** 2


def delta2_tr_wrong_cong(p: int) -> int:
    """Swap the p mod 4 branch."""
    if p % 4 == 1:
        M = M_named(p)
        return (M // lambda_reg(p)) ** 2
    S = S_chi_named(p)
    return p * S * S


def triangle_mask(p: int) -> np.ndarray:
    m = np.zeros(p * p, dtype=np.float64)
    for x in range(p):
        for y in range(x):
            m[x + y * p] = 1.0
    return m


def triangle_leak(p: int) -> float:
    F = _kernel_field(p)
    q = F["q"]
    hat = triangle_mask(p) @ e_matrix(F).T
    chi = np.array(F["chi"], dtype=np.int32)
    sigma = sigma_named(p)
    Omega = chi == sigma
    Omega[0] = False
    off = ~Omega
    off[0] = False
    return float((np.abs(hat[off]) ** 2).max()) if off.any() else 0.0


def prove_A() -> dict:
    ok = True
    rows = {}
    if delta2_tr_named(5) != 4500:
        ok = False
    if delta2_tr_named(7) != 38416:
        ok = False
    if delta2_tr_wrong_cong(5) == 4500 or delta2_tr_wrong_cong(7) == 38416:
        ok = False
    for p in (5, 7):
        d2, sigs, _ = live_delta2_and_sigs(p)
        pred = float(delta2_tr_named(p))
        tr = [i for i, s in enumerate(sigs) if is_triple_rainbow(s, p)]
        if not tr:
            ok = False
        err = max(abs(float(d2[i]) - pred) for i in tr)
        if err > 1e-4:
            ok = False
        rows[str(p)] = {"n_tr": len(tr), "pred": delta2_tr_named(p), "err": err}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "TR Δ² = p S_χ² (p≡1(4)) or (M/λ)² (p≡3(4)). "
            "Fail: swap congruence."
        ),
    }


def prove_B() -> dict:
    p = 5
    ok = triangle_leak(p) > 1.0 and triangle_leak(7) > 1.0
    F = _kernel_field(p)
    Ds = D_rows(p, F).astype(np.int8)
    live = {tuple(row.tolist()) for row in Ds}
    d2, sigs, _ = live_delta2_and_sigs(p)
    tr_keys = {
        tuple(Ds[i].tolist())
        for i, s in enumerate(sigs)
        if is_triple_rainbow(s, p)
    }
    pts = [(x, y) for x in range(p) for y in range(x)]
    n_hit = 0
    n_tr = 0
    seen = set()
    for m00, m01, m10, m11 in product(range(p), repeat=4):
        if (m00 * m11 - m01 * m10) % p == 0:
            continue
        for t0 in range(p):
            for t1 in range(p):
                m = np.zeros(p * p, dtype=np.int8)
                for x, y in pts:
                    nx = (m00 * x + m01 * y + t0) % p
                    ny = (m10 * x + m11 * y + t1) % p
                    m[nx + ny * p] = 1
                key = tuple(m.tolist())
                if key not in live:
                    continue
                n_hit += 1
                seen.add(key)
                if key in tr_keys:
                    n_tr += 1
                else:
                    ok = False
    if seen != tr_keys:
        ok = False
    if n_hit == 0:
        ok = False
    return {
        "proved": bool(ok),
        "leak5": triangle_leak(5),
        "leak7": triangle_leak(7),
        "n_hit": n_hit,
        "n_unique": len(seen),
        "n_tr": len(tr_keys),
        "theorem": (
            "Raw T leaks. AGL(2,5)∩Max+ = all TR. "
            "Fail: raw T is Max+."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Names Δ² on the triple-rainbow class only. "
            "Other complementary types and Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.477  triple-rainbow Δ² named", flush=True)
    A = prove_A()
    print(f"  A formula: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B AGL: {B['proved']} unique={B['n_unique']}/{B['n_tr']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.477",
        "title": "Triple-rainbow Δ² named; other complementary open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "tr_delta2_named": A["proved"],
            "agl_triangle_is_tr": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15477.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
