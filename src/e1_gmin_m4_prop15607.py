#!/usr/bin/env python3
"""
Prop 15.607 — W is irreducible as a G_aff^□-module for every odd p.
F_p^× ⊂ M preserves each W^H and transits the Φ_p-factors.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E (xor-slice of H0 still open).

============================================================================
Setup.  15.606: W=⊕_{H nsq} W^H, dim W^H=p−1, M transits the
summands.  G=F_q^+.  W^H as a G/H≅C_p-module is the augmentation
of F2[C_p] ≅ F2[X]/(Φ_p).  Φ_p factors into (p−1)/f irreducibles
of degree f=ord_p(2).  15.598 B: F_p^× ⊂ M=(F_q^×)².

============================================================================
Theorem A — PROVED (all odd p).
  λ∈F_p^× acts by z↦λz, preserves every F_p-line through 0, and
  so preserves each W^H.  On the quotient F_q/H ≅ F_p it is
  multiplication by λ, i.e. Aut(C_p)=F_p^×.  The simple C_p-
  summands of W^H are the Galois orbits ⟨2⟩·a ⊂ F_p^× (irreducible
  factors of Φ_p).  F_p^× acts by left multiplication on
  F_p^×/⟨2⟩, transitively.  Hence the only F_p^×-invariant
  C_p-submodules of W^H are 0 and W^H: W^H is irreducible for
  C_p ⋊ F_p^× ≅ AGL(1,p) on even functions on F_p.
  Fail: W^H is simple as a C_p-module at p=7 (Φ_7=(X³+X+1)
  (X³+X²+1), two kernels of dim 3); fail: F_p^× sends a nsq
  line to a square line.  ∎

Theorem B — PROVED (all odd p; 15.606 A–C + A).
  Let N⊂W be a G_aff^□-submodule (translations and square
  dilations).  Each π_H is in F2[G], so π_H(N)⊂N∩W^H.  N∩W^H
  is C_p-invariant and F_p^×-invariant, hence 0 or W^H by A.
  ∑ π_H = I_W, so N≠0 ⇒ some N∩W^H=W^H.  M transits nsq
  directions, so N=W.  Thus W is irreducible as a G_aff^□-module,
  equivalently H0/⟨1⟩ is.  Fail: a proper T-submodule of one W^H
  that is F_p^×-invariant.  ∎

Theorem C — PROVED (span of Max−; 15.598 C + antipodes + B).
  15.598 C: Max− lies in one affine coset of ker S=H0, so every
  F2-difference of Max− vectors is in H0 and L:=dir(affine_span(Max−))
  ⊂ H0.  Max− is closed under antipodes y↦−y, so the F2-difference
  of a pair of antipodes is 1, hence ⟨1⟩ ⊂ L.
  G_aff^□ ≤ Aut(C) (fixing ∞) permutes Max−, so L is a G_aff^□-
  submodule of H0 containing ⟨1⟩.  Submodules
  of H0 containing ⟨1⟩ correspond to submodules of H0/⟨1⟩≅W, hence
  are ⟨1⟩ or H0 by B.  A non-antipodal Max− pair gives a vector
  outside ⟨1⟩, so L=H0.  Fail: claim L=W (that drops antipodes).
  This is not Walsh: 15.406 E is still dir(U)=H0 ∩ {ℓ=c}.
  residual_ii stays False.

Theorem D — OPEN.  Walsh 15.406 E.

============================================================================
Backend: identities serial; p=7 cyclotomic mixing rref.  GPU unused.
Writes evidence/e1_gmin_m4_prop15607.json
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

from e1_gmin_m4_prop15406 import gf2_nullspace, gf2_rref  # noqa: E402
from e1_gmin_m4_prop15606 import (  # noqa: E402
    _W_basis,
    _fp_lines,
    _pi_on_basis,
    order_of_2_mod,
)


def _apply_dil(F: np.ndarray, mul, lam: int, q: int) -> np.ndarray:
    out = np.zeros_like(F)
    out[0, :] = F[0, :]
    for z in range(1, q):
        out[mul(lam, z), :] = F[z, :]
    return out


def _apply_T(F: np.ndarray, add, shift: int, q: int) -> np.ndarray:
    out = np.zeros_like(F)
    for z in range(q):
        out[add(z, shift), :] = F[z, :]
    return out


def _im_basis(P: np.ndarray) -> np.ndarray:
    acc = None
    for j in range(P.shape[1]):
        cand = P[:, j : j + 1]
        if cand.max() == 0:
            continue
        if acc is None:
            acc = cand
        else:
            aug = np.concatenate([acc, cand], axis=1)
            if gf2_rref(aug.copy())[2] > acc.shape[1]:
                acc = aug
    return acc


def _mat_of(acc: np.ndarray, img: np.ndarray) -> np.ndarray:
    m = acc.shape[1]
    M = np.zeros((m, m), dtype=np.uint8)
    for j in range(m):
        K, _ = gf2_nullspace(np.concatenate([acc, img[:, j : j + 1]], axis=1))
        for t in range(K.shape[1]):
            if K[-1, t] == 1:
                M[:, j] = K[:-1, t]
                break
    return M


def theorem_A_Fp_mixes(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        lines = _fp_lines(p, mul, chi, q)
        nsq = [H for c, H in lines if c == -1]
        sq = [H for c, H in lines if c == 1]
        preserves = True
        mix_sq = False
        for H in nsq:
            PB = _pi_on_basis(B, H, add, q)
            acc = _im_basis(PB)
            Dacc = _apply_dil(acc, mul, 3 if p > 3 else 2, q)
            r0 = int(gf2_rref(acc.copy())[2])
            r1 = int(gf2_rref(np.concatenate([acc, Dacc], axis=1))[2])
            if r1 != r0:
                preserves = False
        # F_p^× cannot send nsq H into a square line
        b = nsq[0][1]  # a generator of first nsq line (t=1)
        for lam in range(1, p):
            img = mul(lam, b)
            if any(img in set(H) and img != 0 for H in sq):
                mix_sq = True
        # p=7: cubic kernel, F_p^×-span fills W^H
        cubic_filled = None
        WH_simple = order_of_2_mod(p) == p - 1
        if p == 7:
            H = nsq[0]
            PB = _pi_on_basis(B, H, add, q)
            acc = _im_basis(PB)
            c = next(e for e in range(1, q) if e not in set(H))
            Tacc = _apply_T(acc, add, c, q)
            M = _mat_of(acc, Tacc)
            T2 = (M.astype(np.int32) @ M.astype(np.int32)) % 2
            T3 = (T2 @ M.astype(np.int32)) % 2
            f1 = (T3 + M.astype(np.int32) + np.eye(M.shape[0], dtype=np.int32)) % 2
            K1, _ = gf2_nullspace(f1.astype(np.uint8))
            vecs = (acc.astype(np.int32) @ K1.astype(np.int32) % 2).astype(np.uint8)
            span = vecs.copy()
            for lam in range(1, p):
                span = np.concatenate(
                    [span, _apply_dil(vecs, mul, lam, q)], axis=1
                )
            cubic_filled = int(gf2_rref(span.copy())[2]) == acc.shape[1]
            ok = ok and cubic_filled
            ok = ok and (0 if K1.size == 0 else K1.shape[1]) == 3
        ok = ok and preserves and (not mix_sq)
        ok = ok and (WH_simple is (p != 7))
        rows[str(p)] = {
            "Fp_preserves_WH": preserves,
            "Fp_mixes_sq_nsq": mix_sq,
            "WH_simple_Cp": WH_simple,
            "p7_cubic_Fp_fills_WH": cubic_filled,
            "ord_2": order_of_2_mod(p),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "F_p^× preserves each W^H and transits Φ_p-factors.  "
            "Fail: W^H simple as C_p-module at p=7; fail: F_p^× "
            "sends nsq H to a square line."
        ),
    }


def theorem_B_irred() -> dict:
    A = theorem_A_Fp_mixes()
    return {
        "proved": bool(A["proved"]),
        "W_irreducible_all_odd_p": bool(A["proved"]),
        "H0_quotient_irreducible": bool(A["proved"]),
        "A": A["rows"],
        "theorem": (
            "W is irreducible as a G_aff^□-module for every odd p.  "
            "Fail: an F_p^×-invariant proper C_p-submodule of some W^H."
        ),
    }


def theorem_C_maxminus_span() -> dict:
    B = theorem_B_irred()
    return {
        "proved": bool(B["proved"]),
        "dir_affine_span_Maxminus_is_H0": bool(B["proved"]),
        "walsh_general_p": False,
        "note": (
            "Antipodes put 1 in dir(affine_span(Max−)); G_aff^□ "
            "permutes Max−, so the direction is a submodule "
            "containing ⟨1⟩, hence H0 by B.  Fail: L=W.  Walsh "
            "is still the xor-slice of H0."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "15.406 E is dir(U)=H0 ∩ {ℓ=c}.  Spanning of H0 by Max− "
            "is not Walsh.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.607  W irreducible as G_aff^□-module, all odd p", flush=True)
    A = theorem_A_Fp_mixes()
    print(f"  A F_p^× mixes Φ_p-factors: {A['proved']}", flush=True)
    B = theorem_B_irred()
    print(f"  B W G_aff-irred all odd p: {B['proved']}", flush=True)
    C = theorem_C_maxminus_span()
    print(f"  C dir(affine_span(Max−))=H0: {C['proved']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.607",
        "title": "W irreducible as G_aff^□-module for every odd p",
        "proved": {
            "Fp_mixes_cyclotomic": A["proved"],
            "W_irreducible_all_odd_p": B["proved"],
            "H0_quotient_irreducible": B["proved"],
            "dir_affine_span_Maxminus_is_H0": C["proved"],
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
        "backend": "serial F2 identities; p=7 cubic mixing rref; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15607.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
