#!/usr/bin/env python3
"""
Prop 15.601 — Paley QR indicator is in rowspan(S) or rowspan(S)+span{ℓ}.
A general-p F2 identity (Max-free).  Does **not** close Walsh / 15.406 E:
the proposed single-orbit Aut_e spanning proof is false at p=5 (Claude and
OpenAI both BLOCK).  residual_ii stays False.

============================================================================
Setup (15.598–15.600).  S = incidence of square-direction {∞}∪L.
ℓ = e_∞ + e_0.  QR = {e∈F_q^× : χ(e)=1}, 1_QR ∈ F_2^{P¹}.

============================================================================
Theorem A — PROVED (pencil; all odd p; Max-free).
  Let w be the indicator of the (p+1)/2 square-direction affine lines
  through 0.  For each point x,
      (Sᵀ w)_x = # { such lines through x }  (mod 2).
    x=0 or ∞: that count is (p+1)/2.
    x finite ≠0: the unique line 0x has square direction iff χ(x)=1,
      so the count is 1_{QR}(x).
  Hence Sᵀ w = ((p+1)/2)(e_0+e_∞) + 1_QR.
  (p+1)/2 is odd iff p≡1 (mod 4), even iff p≡3 (mod 4), therefore
      p≡1 (mod 4)  ⇒  1_QR + ℓ ∈ rowspan(S),
      p≡3 (mod 4)  ⇒  1_QR ∈ rowspan(S).
  Fail: swap the two congruences.  ∎

Theorem B — PROVED (corollary on H0).
  On H0=ker S, rowspan annihilates, so
      1_QR · x = ℓ(x)   if p≡1 (mod 4),
      1_QR · x = 0      if p≡3 (mod 4).
  In particular 1_QR is constant on each xor-slice of H0 (equal to ℓ,
  or 0).  It is not an extra dual of A∩{ℓ=c}.  Fail: claim 1_QR is
  independent of {S,ℓ} on H0.  ∎

Theorem C — OPEN.  Walsh (15.406 E) is still dir(U)=ker S ∩ ker ℓ.
  Aut_e-invariant extra duals are ruled out by B, but Aut_e is
  reducible and a single Aut_e-orbit of a U-difference does not span
  (p=5: ≤11 < 12).  residual_ii stays False.

============================================================================
Writes evidence/e1_gmin_m4_prop15601.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime, field_ctx  # noqa: E402
from walsh_linecode_rank import square_line_matrix  # noqa: E402


def qr_and_ell(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1
    qr = np.zeros(n, dtype=np.uint8)
    for e in range(1, q):
        if chi(e) == 1:
            qr[1 + e] = 1
    ell = np.zeros(n, dtype=np.uint8)
    ell[0] = 1
    ell[1] = 1
    return qr, ell, chi, q, mul, add


def pencil_through_zero(p: int) -> tuple[np.ndarray, np.ndarray]:
    """Rows of S that contain field-0, and the 0/1 selector w of length b."""
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    S = square_line_matrix(p)
    # Rebuild lines in the same order as square_line_matrix and mark those with 0.
    used = set()
    dirs = []
    for b in range(1, q):
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    w = []
    for b in dirs:
        if chi(b) != 1:
            continue
        covered = set()
        for a in range(q):
            if a in covered:
                continue
            pts = []
            for t in range(p):
                e = add(a, mul(t, b))
                pts.append(e)
                covered.add(e)
            w.append(1 if 0 in pts else 0)
    w = np.array(w, dtype=np.uint8)
    assert len(w) == S.shape[0]
    return S, w


def theorem_A_pencil(primes=None) -> dict:
    if primes is None:
        primes = (5, 7, 11, 13, 17, 19)
    ok = True
    rows = {}
    for p in primes:
        S, w = pencil_through_zero(p)
        qr, ell, chi, q, mul, add = qr_and_ell(p)
        got = (S.astype(np.int32).T @ w.astype(np.int32)) % 2
        got = got.astype(np.uint8)
        half = (p + 1) // 2
        if p % 4 == 1:
            target = (qr ^ ell) & 1
        else:
            target = qr
        match = bool(np.array_equal(got, target))
        # fail-when-wrong: swapped congruence
        swapped = qr if p % 4 == 1 else (qr ^ ell) & 1
        ok = ok and match
        ok = ok and int(w.sum()) == half
        if np.array_equal(got, swapped) and not np.array_equal(target, swapped):
            ok = False
        rows[str(p)] = {
            "wt_w": int(w.sum()),
            "half": half,
            "p_mod_4": p % 4,
            "match": match,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Sᵀw=((p+1)/2)(e_0+e_∞)+1_QR. p≡1⇒QR+ℓ∈rowspan S; "
            "p≡3⇒QR∈rowspan S. Fail: swap congruences."
        ),
    }


def theorem_B_on_H0(primes=None) -> dict:
    if primes is None:
        primes = (5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        S, w = pencil_through_zero(p)
        qr, ell, chi, q, mul, add = qr_and_ell(p)
        rS = gf2_rref(S)[2]
        rQ = gf2_rref(np.vstack([S, qr]))[2]
        rQE = gf2_rref(np.vstack([S, (qr ^ ell) & 1]))[2]
        rEll = gf2_rref(np.vstack([S, ell]))[2]
        if p % 4 == 1:
            ok = ok and rQE == rS and rQ == rS + 1
        else:
            ok = ok and rQ == rS and rQE == rS + 1
        ok = ok and rEll == rS + 1  # ℓ never in rowspan
        rows[str(p)] = {
            "rS": int(rS),
            "rQ": int(rQ),
            "rQE": int(rQE),
            "ell_extra": rEll == rS + 1,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "On H0, QR·x=ℓ(x) (p≡1) or 0 (p≡3). Fail: QR independent of {S,ℓ}.",
    }


def theorem_C_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "single_orbit_spans_W": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "15.406 E stays OPEN. Single Aut_e-orbit spanning is false at p=5. "
            "residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.601  QR indicator in rowspan(S) or S+ℓ", flush=True)
    A = theorem_A_pencil()
    print(f"  A pencil: {A['proved']}", flush=True)
    B = theorem_B_on_H0()
    print(f"  B on H0: {B['proved']}", flush=True)
    C = theorem_C_open()
    print(f"  C Walsh open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.601",
        "title": "QR indicator is in rowspan(S) or rowspan(S)+span{ℓ}",
        "proved": {
            "pencil_identity": A["proved"],
            "QR_on_H0": B["proved"],
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15601.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
