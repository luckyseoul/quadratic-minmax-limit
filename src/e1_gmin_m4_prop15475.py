#!/usr/bin/env python3
"""
Prop 15.475 — Δ² on the parallel-union class is named:
    p≡3 (mod 4):  Δ² = M² = [p²(p²−1)/4]²
    p≡1 (mod 4):  Δ² = p⁵(p−1)²/4
from the F_p Gauss sum on a union of (p−1)/2 parallel lines
whose annihilator lies in Ω.  Complementary Max+ have a
different Δ².  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.474 flags.

============================================================================
Setup.  15.474: hat 1_D=0 off Ω, W=τ̄Δ, Δ not constant.
M=k(q−k)=p²(p²−1)/4.  A parallel-union is D=ℓ^{-1}(A) for
F_p-linear surjective ℓ and |A|=(p−1)/2.  After dilating
into Ω when σ=−1, support is an annihilator line of size p−1.

============================================================================
Theorem A — PROVED (Gauss on F_p; certified p=3,5,7).
  On this construction, p≡3(4) ⇒ ψ≡1 on the annihilator and
  |hat|² is constant, so Δ²=M².  p≡1(4) ⇒ two energy levels
  (p²/4)(√p±1)² and Δ=± p²(p−1)√p/2, hence
      Δ² = p⁵(p−1)²/4.
  Fail: use the p≡3 formula at p=5 (22500 ≠ 12500).  ∎

Theorem B — PROVED (cache fold).
  Live parallel-unions (n_full=(p−1)/2) match A: 6, 30, 140
  rows.  Complementary rows have a different Δ²
  (4500 at p=5; 0 or 38416 at p=7).
  Fail: claim every Max+ has the parallel-union Δ² at p=5.  ∎

Theorem C — OPEN.  Complementary Δ² unnamed.  Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: construction DFT + Fraction.  Writes
evidence/e1_gmin_m4_prop15475.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows, k_named  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15475_catalog.json"


def M_named(p: int) -> int:
    k = k_named(p)
    return k * (p * p - k)


def delta2_parunion_named(p: int) -> int:
    if p % 4 == 3:
        return M_named(p) ** 2
    return (p ** 5) * (p - 1) ** 2 // 4


def delta2_parunion_wrong_cong(p: int) -> int:
    """Swap the congruence: p≡3 formula used at every p."""
    return M_named(p) ** 2


def _qr(p: int) -> set[int]:
    return {a for a in range(1, p) if pow(a, (p - 1) // 2, p) == 1}


def parunion_mask(F, p: int) -> np.ndarray:
    """Tr^{-1}(QR), dilated by a trace-zero element when p≡1 (mod 4)."""
    q = F["q"]
    trt = F["trt"]
    qr = _qr(p)
    mask = np.array([int(trt[x]) in qr for x in range(q)], dtype=np.float64)
    if p % 4 == 1:
        th = next(x for x in range(1, q) if trt[x] == 0)
        mul = F["mul"]
        out = np.zeros(q, dtype=np.float64)
        for x in range(q):
            if mask[x]:
                out[mul[th][x]] = 1.0
        return out
    return mask


def delta2_of_mask(F, mask: np.ndarray, p: int) -> tuple[float, float]:
    q = F["q"]
    hat = mask @ e_matrix(F).T
    chi = np.array(F["chi"], dtype=np.int32)
    sigma = sigma_named(p)
    Omega = chi == sigma
    Omega[0] = False
    off = ~Omega
    off[0] = False
    leak = float((np.abs(hat[off]) ** 2).max()) if off.any() else 0.0
    psi = make_psi(F)
    plus = Omega & ((psi.imag > 0) if sigma < 0 else (psi.real > 0))
    Up = float((np.abs(hat[plus]) ** 2).sum())
    Um = float((np.abs(hat[Omega]) ** 2).sum() - Up)
    return (Up - Um) ** 2, leak


def live_delta2(p: int) -> np.ndarray:
    F = _kernel_field(p)
    Ds = D_rows(p, F).astype(np.float64)
    hat = Ds @ e_matrix(F).T
    chi = np.array(F["chi"], dtype=np.int32)
    sigma = sigma_named(p)
    Omega = chi == sigma
    Omega[0] = False
    psi = make_psi(F)
    plus = Omega & ((psi.imag > 0) if sigma < 0 else (psi.real > 0))
    Up = (np.abs(hat[:, plus]) ** 2).sum(1)
    Um = (np.abs(hat[:, Omega]) ** 2).sum(1) - Up
    return (Up - Um) ** 2


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        mask = parunion_mask(F, p)
        d2, leak = delta2_of_mask(F, mask, p)
        pred = delta2_parunion_named(p)
        if leak > 1e-10 or abs(d2 - pred) > 1e-4:
            ok = False
        if int(mask.sum()) != k_named(p):
            ok = False
        rows[str(p)] = {"d2": float(d2), "pred": pred, "leak": leak}
    if abs(delta2_parunion_wrong_cong(5) - 12500) < 1:
        ok = False
    if delta2_parunion_named(5) != 12500:
        ok = False
    if delta2_parunion_named(7) != M_named(7) ** 2:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Parallel-union Δ² = M² (p≡3(4)) or p⁵(p−1)²/4 (p≡1(4)). "
            "Fail: p≡3 formula at p=5."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    expected = {3: 6, 5: 30, 7: 140}
    for p in (3, 5, 7):
        d2 = live_delta2(p)
        pred = float(delta2_parunion_named(p))
        nmatch = int(np.sum(np.abs(d2 - pred) < 1.0))
        ntot = int(d2.shape[0])
        if nmatch != expected[p]:
            ok = False
        if p == 5 and nmatch == ntot:
            ok = False
        rows[str(p)] = {"nmatch": nmatch, "ntot": ntot, "pred": pred}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Live parallel-unions match A; complementary Δ² differs. "
            "Fail: all p=5 Max+ have Δ²=12500."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Names Δ² on the parallel-union class only. "
            "Complementary Δ² and Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.475  parallel-union Δ² named", flush=True)
    A = prove_A()
    print(f"  A formula: {A['proved']} p5={A['rows']['5']['pred']}", flush=True)
    B = prove_B()
    print(f"  B live: {B['proved']} { {k: v['nmatch'] for k,v in B['rows'].items()} }", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.475",
        "title": "Δ² named on parallel-union class; complementary open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "parunion_delta2_named": A["proved"],
            "complement_differs": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15475.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
