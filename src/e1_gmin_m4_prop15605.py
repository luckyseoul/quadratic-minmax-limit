#!/usr/bin/env python3
"""
Prop 15.605 — Paley adjacency is an F2-projection; H0 = ⟨1⟩ ⊕ W
with W the translate-span of the 15.604 extra vector, dim (q−1)/2.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** prove W irreducible as a G_aff^□-module, nor Walsh 15.406 E.

============================================================================
Setup.  Paley graph of order q=p² (p odd ⇒ q≡1 (mod 8)): vertices F_q,
A_{xy}=1 iff χ_q(x−y)=1.  P_{xy}=1 iff χ_q(x−y)=−1.  Degree
N=(q−1)/2 even.  extra = 1_QR if p≡1 (mod 4), else 1_QNR (15.604 A).
W ⊂ F2^{P¹} is the F2-span of additive translates of extra (∞ fixed,
extra_∞=0).  S_aff = affine square-line incidence on F_q (no ∞).

============================================================================
Theorem A — PROVED (srg parameters; all odd p).
  Paley is srg(q, (q−1)/2, (q−5)/4, (q−1)/4), so
      A² = ((q−1)/4) J + ((q−1)/4) I − A.
  q=p²≡1 (mod 8) ⇒ (q−1)/4 is even, hence A²=A over F2.
  Row-sums N even ⇒ A1=0.  P=A+I+J over F2 (off-diagonal partition
  QR / QNR, diagonal I).  Then P²=A+I+J=P, using A²=A, AJ=JA=0,
  J²=J (q odd).  Fail: the same for Paley of order 13 (q≡5 (mod 8),
  (q−1)/4 odd, A²≠A); fail: A²=0.  ∎

Theorem B — PROVED (15.600, 15.601 pencil, 15.604 A; all odd p).
  dim{x∈H0: x_∞=0} = dim H0 − 1 = N, and this slice is ker S_aff
  (S x = x_∞ + S_aff x_aff).  W ⊂ ker S_aff: extra ∈ H0 (15.604)
  with extra_∞=0, and translations permute square affine lines
  (15.602 A).  1 ∉ W (every vector of W vanishes at ∞).

  p≡1: (p+1)/2 odd.  15.601: sum of square 0-lines is 1_QR+e_0
  on F_q.  Translating, the Paley neighborhood {τ}∪(τ+QR) is the
  F2-sum of square lines through τ.  So f∈ker S_aff ⇒ Af=f ⇒
  ker S_aff ⊂ im A.  Columns of A are the affine restrictions of
  the QR-translates, so W=im A (embedded at ∞=0).  Sandwich:
  N = dim ker S_aff ≤ rank A = dim W ≤ N.

  p≡3: (p+1)/2 even.  Sum of square lines through τ is 1_{τ+QR}
  (no e_τ), so Af=0 on ker S_aff.  Even weight (parallel-class
  sums) and P=A+I+J give Pf=f, so ker S_aff ⊂ im P = W.

  Thus H0 = ⟨1⟩ ⊕ W with dim W = N.  Fail: dim W = N−1;
  fail: 1 ∈ W.  ∎

Theorem C — PROVED as a G_aff^□-splitting; irreducibility OPEN.
  W is G_aff^□-invariant: translations by construction; square
  dilation fixes extra and permutes its translates; Frobenius
  preserves QR and QNR.  extra generates W as a translation
  module.  15.602 B: the unique 1-dim G-invariant of H0 is ⟨1⟩,
  so W^G=0.  Irreducibility of W (equivalently of H0/⟨1⟩) is
  OPEN — char-2 averaging over D does not force every submodule
  to meet ⟨extra⟩ (N even).

Theorem D — OPEN.  Walsh 15.406 E still needs the xor-slice.
  residual_ii stays False.

============================================================================
Backend: identities serial (inherently); F2 matmul/rref p=3..11.
GPU unused.  Writes evidence/e1_gmin_m4_prop15605.json
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
from e1_gmin_m4_prop15598 import field_ctx, legendre  # noqa: E402
from e1_gmin_m4_prop15604 import _in_ker, _qr_qnr  # noqa: E402
from walsh_linecode_rank import square_line_matrix  # noqa: E402


def _sub(add, mul, x: int, y: int, p: int) -> int:
    if y == 0:
        return x
    return add(x, mul(y, p - 1))


def paley_AP(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    A = np.zeros((q, q), dtype=np.uint8)
    P = np.zeros((q, q), dtype=np.uint8)
    for x in range(q):
        for y in range(q):
            if x == y:
                continue
            d = _sub(add, mul, x, y, p)
            c = chi(d)
            if c == 1:
                A[x, y] = 1
            elif c == -1:
                P[x, y] = 1
    return A, P, q, mul, add, chi


def paley_prime(q: int) -> np.ndarray:
    """Paley graph of prime order q≡1 (mod 4)."""
    A = np.zeros((q, q), dtype=np.uint8)
    for x in range(q):
        for y in range(q):
            if x == y:
                continue
            if legendre(x - y, q) == 1:
                A[x, y] = 1
    return A


def theorem_A_projection(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        A, P, q, mul, add, chi = paley_AP(p)
        A2 = (A.astype(np.int32) @ A.astype(np.int32)) % 2
        P2 = (P.astype(np.int32) @ P.astype(np.int32)) % 2
        J = np.ones((q, q), dtype=np.uint8)
        I = np.eye(q, dtype=np.uint8)
        recon = (A + I + J) % 2
        ones = np.ones(q, dtype=np.int32)
        A1 = (A.astype(np.int32) @ ones) % 2
        P1 = (P.astype(np.int32) @ ones) % 2
        rA = int(gf2_rref(A)[2])
        rP = int(gf2_rref(P)[2])
        N = (q - 1) // 2
        even_q14 = ((q - 1) // 4) % 2 == 0
        ok = ok and np.array_equal(A2, A)
        ok = ok and np.array_equal(P2, P)
        ok = ok and np.array_equal(recon, P)
        ok = ok and int(A1.max()) == 0
        ok = ok and int(P1.max()) == 0
        ok = ok and even_q14
        ok = ok and rA == N == rP
        ok = ok and not np.array_equal(A2, np.zeros_like(A))
        rows[str(p)] = {
            "q": q,
            "A2_eq_A": True,
            "P2_eq_P": True,
            "rank_A": rA,
            "rank_P": rP,
            "N": N,
            "q_minus_1_over_4_even": even_q14,
        }
    # fail-when-wrong: Paley of order 13 (≡5 mod 8)
    A13 = paley_prime(13)
    A13_2 = (A13.astype(np.int32) @ A13.astype(np.int32)) % 2
    ok = ok and not np.array_equal(A13_2, A13)
    ok = ok and ((13 - 1) // 4) % 2 == 1
    return {
        "proved": bool(ok),
        "paley13_is_projection": bool(np.array_equal(A13_2, A13)),
        "rows": rows,
        "theorem": (
            "Paley A of order p² satisfies A²=A over F2; P=A+I+J "
            "likewise.  Fail: Paley of order 13; fail: A²=0."
        ),
    }


def theorem_B_splitting(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    Aproj = theorem_A_projection(primes)
    ok = Aproj["proved"]
    rows = {}
    for p in primes:
        S = square_line_matrix(p)
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        extra = qr if p % 4 == 1 else qnr
        N = (q - 1) // 2
        cols = []
        in_h0 = True
        for tau in range(q):
            v = np.zeros(q + 1, dtype=np.uint8)
            for z in range(q):
                v[1 + add(z, tau)] = extra[1 + z]
            in_h0 = in_h0 and _in_ker(S, v)
            cols.append(v)
        W = np.stack(cols, axis=1)
        dimW = int(gf2_rref(W)[2])
        one = np.ones(q + 1, dtype=np.uint8)
        dimW1 = int(gf2_rref(np.concatenate([W, one[:, None]], axis=1))[2])
        one_in_W = dimW1 == dimW
        # extra itself is the τ=0 translate
        extra_in = _in_ker(S, extra) and extra[0] == 0
        n = q + 1
        ok = ok and in_h0 and extra_in
        ok = ok and dimW == N
        ok = ok and dimW1 == n // 2
        ok = ok and not one_in_W
        ok = ok and dimW != N - 1
        rows[str(p)] = {
            "dim_W": dimW,
            "N": N,
            "dim_W_plus_1": dimW1,
            "dim_H0": n // 2,
            "translates_in_H0": in_h0,
            "one_in_W": bool(one_in_W),
            "extra_infty_zero": bool(extra[0] == 0),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "H0 = ⟨1⟩ ⊕ W, W = span of extra-translates, dim N.  "
            "Fail: dim W = N−1; fail: 1 ∈ W."
        ),
    }


def theorem_C_module_open(primes=None) -> dict:
    """W is a G_aff complement to ⟨1⟩; irreducibility not claimed."""
    B = theorem_B_splitting(primes)
    return {
        "proved": False,
        "splitting": B["proved"],
        "W_irreducible": False,
        "H0_quotient_irreducible": False,
        "rows": B["rows"],
        "note": (
            "W is G_aff^□-invariant, generated by extra, and W^G=0 "
            "by 15.602 B.  Irreducibility of W is OPEN."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "H0_quotient_irreducible": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "H0=⟨1⟩⊕W does not make W irreducible and does not span "
            "the xor-slice.  15.406 E and residual_ii stay OPEN/False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.605  Paley F2-projection; H0=⟨1⟩⊕W", flush=True)
    A = theorem_A_projection()
    print(f"  A Paley A²=A: {A['proved']} (Paley13 proj={A['paley13_is_projection']})", flush=True)
    B = theorem_B_splitting()
    print(f"  B H0=⟨1⟩⊕W dim N: {B['proved']}", flush=True)
    C = theorem_C_module_open()
    print(f"  C W irred open: splitting={C['splitting']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.605",
        "title": "Paley A²=A over F2; H0=⟨1⟩⊕ translate-span of extra",
        "proved": {
            "paley_projection": A["proved"],
            "H0_splits": B["proved"],
            "W_irreducible": False,
            "H0_quotient_irreducible": False,
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "serial F2 identities (inherently); rref p=3..11; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15605.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
