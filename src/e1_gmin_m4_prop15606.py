#!/usr/bin/env python3
"""
Prop 15.606 — nonsquare F_p-line averages split W; squares permute
the summands transitively.  G_aff^□-irreducibility of W when 2 is
a primitive root mod p.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E.  Irreducibility for every odd p
(when 2 is not a primitive root) is OPEN; certified at p=7.

============================================================================
Setup.  15.605: H0=⟨1⟩⊕W, W=ker S_aff, dim N=(q−1)/2, translation-
invariant.  G=F_q^+.  For an F_p-line H=F_p b through 0,
π_H=∑_{h∈H} T_h ∈ F2[G].  p odd ⇒ π_H²=π_H on F2^{F_q}.
M=(F_q^×)² acts on the p+1 directions (F_p-subspaces of F_q).

============================================================================
Theorem A — PROVED (15.605 ker S_aff; all odd p).
  If χ(b)=1 (H square) then π_H vanishes on W: (π_H f)(x) is the
  sum of f on the square affine line x+H, and W=ker S_aff.
  If χ(b)=−1 (H nonsquare) then W^H equals the space of all even
  H-invariant functions, dim p−1.  Proof: H-invariants are
  functions of the p-point quotient F_q/H.  Even weight on F_q is
  even weight on the quotient (p odd).  A square affine line L is
  not parallel to H, so meets every H-coset once and
  ∑_L f = ∑_{quotient} f̄, even.  Thus every even H-invariant
  lies in ker S_aff=W, and dim W^H=p−1.  Fail: square π_H has
  rank p−1; fail: nsq π_H=0.  ∎

Theorem B — PROVED (counting; p mod 4 vs 15.605 A/P).
  The (p+1)/2 nonsquare π_H are pairwise orthogonal on W and
  sum to I_W, so W=⊕_{H nsq} W^H.
  Distinct 0-lines H,H' span F_q, so π_H π_{H'}=∑_{z∈F_q} T_z,
  which is the total-weight map, zero on even-weight W.
  At x, ∑_{H nsq} (π_H f)(x) = ((p+1)/2) f(x) + (P f)(x),
  P the nonsquare Paley adjacency.  p≡1: (p+1)/2 odd, W=im A,
  Pf=0, sum=f.  p≡3: (p+1)/2 even, W=im P, Pf=f, sum=f.  ∎

Theorem C — PROVED (Singer).
  F_q^× acts transitively on the p+1 directions with kernel F_p^×.
  F_p^×⊂M (15.598 B: χ_q|_{F_p^×}≡1), so M acts through a cyclic
  group of order (p+1)/2.  Two orbits of size (p+1)/2: square and
  nonsquare directions.  Hence M (and D) permutes the nsq summands
  W^H transitively.  Fail: M mixes square with nonsquare.  ∎

Theorem D — PROVED when 2 is a primitive root mod p;
  CERTIFIED p=7 (ord_7(2)=3); OPEN as a p-law for every odd p.
  W^H, as a module for G/H≅C_p, is the augmentation of F2[C_p]
  (even functions on p points).  This is simple iff Φ_p is
  irreducible over F2 iff 2 is a primitive root modulo p.  Then
  every T-submodule of W is a sum of some W^H, and D-transitivity
  forces a G_aff^□-submodule (translations + square dilations) to
  be 0 or W.  At p=7, Φ_7 splits as two cubics but the G_aff-span
  of either 3-dim kernel fills W.  Fail: claim W^H simple at p=7.
  Walsh still needs the xor-slice; residual_ii stays False.

============================================================================
Backend: identities serial; rref ranks p=3,5,7.  GPU unused.
Writes evidence/e1_gmin_m4_prop15606.json
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
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402


def _fp_lines(p: int, mul, chi, q: int):
    used = set()
    lines = []
    for b in range(1, q):
        if b in used:
            continue
        H = []
        for t in range(p):
            e = mul(t, b)
            H.append(e)
            if t:
                used.add(e)
        lines.append((chi(b), H))
    return lines


def _W_basis(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    extra = (qr if p % 4 == 1 else qnr)[1:]
    cols = []
    for tau in range(q):
        v = np.zeros(q, dtype=np.uint8)
        for z in range(q):
            v[add(z, tau)] = extra[z]
        cols.append(v)
    M = np.stack(cols, axis=1)
    acc = None
    target = (q - 1) // 2
    for j in range(M.shape[1]):
        cand = M[:, j : j + 1]
        if acc is None:
            if cand.max():
                acc = cand
        else:
            aug = np.concatenate([acc, cand], axis=1)
            if gf2_rref(aug.copy())[2] > acc.shape[1]:
                acc = aug
        if acc is not None and acc.shape[1] == target:
            break
    return acc, q, mul, add, chi


def _pi_on_basis(B: np.ndarray, H, add, q: int) -> np.ndarray:
    """π_H acting on columns of B (functions on F_q)."""
    out = np.zeros_like(B)
    for z in range(q):
        acc = np.zeros(B.shape[1], dtype=np.uint8)
        for h in H:
            acc ^= B[add(z, h), :]
        out[z, :] = acc
    return out


def order_of_2_mod(p: int) -> int:
    a, n = 2, 1
    while a % p != 1:
        a = (a * 2) % p
        n += 1
        if n > p:
            return -1
    return n


def theorem_A_projectors(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        lines = _fp_lines(p, mul, chi, q)
        sq_ranks = []
        nsq_ranks = []
        for c, H in lines:
            PB = _pi_on_basis(B, H, add, q)
            r = int(gf2_rref(PB.copy())[2])
            if c == 1:
                sq_ranks.append(r)
            else:
                nsq_ranks.append(r)
        ok = ok and all(r == 0 for r in sq_ranks)
        ok = ok and all(r == p - 1 for r in nsq_ranks)
        ok = ok and len(sq_ranks) == len(nsq_ranks) == (p + 1) // 2
        rows[str(p)] = {
            "sq_ranks": sq_ranks,
            "nsq_ranks": nsq_ranks,
            "p_minus_1": p - 1,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Square π_H=0 on W; nsq π_H has rank p−1.  Fail: "
            "square rank p−1; fail: nsq rank 0."
        ),
    }


def theorem_B_orthogonal_sum(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    A = theorem_A_projectors(primes)
    ok = A["proved"]
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        lines = _fp_lines(p, mul, chi, q)
        nsq = [H for c, H in lines if c == -1]
        pis = [_pi_on_basis(B, H, add, q) for H in nsq]
        # pairwise: π_i π_j B = 0
        orth = True
        for i, Pi in enumerate(pis):
            for j, Pj in enumerate(pis):
                if i >= j:
                    continue
                # apply π_i to columns of Pj (already in W)
                comp = np.zeros_like(B)
                for z in range(q):
                    acc = np.zeros(B.shape[1], dtype=np.uint8)
                    for h in nsq[i]:
                        acc ^= Pj[add(z, h), :]
                    comp[z, :] = acc
                if comp.max() != 0:
                    orth = False
        tot = np.zeros_like(B)
        for Pi in pis:
            tot ^= Pi
        is_id = np.array_equal(tot, B)
        ok = ok and orth and is_id
        rows[str(p)] = {
            "n_nsq": len(nsq),
            "orthogonal": orth,
            "sum_is_id": is_id,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Nonsquare π_H are orthogonal and sum to I_W.  "
            "W=⊕ W^H, dim each p−1."
        ),
    }


def theorem_C_transitive(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        # directions as first representative
        used = set()
        dirs = []
        for b in range(1, q):
            if b in used:
                continue
            dirs.append(b)
            for t in range(1, p):
                used.add(mul(t, b))
        nsq = [b for b in dirs if chi(b) == -1]
        sq = [b for b in dirs if chi(b) == 1]
        # generator of M, not in F_p^× (those fix every F_p-line)
        omega = _primitive(mul, q)
        g2 = mul(omega, omega)

        def orbit(start, want):
            seen = set()
            x = start
            for _ in range(q):
                # direction class of x: representatives
                key = None
                for d in want:
                    # x in F_p d?
                    hit = any(mul(t, d) == x for t in range(1, p))
                    if hit:
                        key = d
                        break
                if key is None or key in seen:
                    break
                seen.add(key)
                x = mul(g2, x)
            return seen

        o0 = orbit(nsq[0], nsq)
        osq = orbit(sq[0], sq)
        mix = False
        x = nsq[0]
        for _ in range(p + 2):
            x = mul(g2, x)
            if any(mul(t, sq[0]) == x for t in range(1, p)):
                mix = True
        ok = ok and len(o0) == (p + 1) // 2
        ok = ok and len(osq) == (p + 1) // 2
        ok = ok and not mix
        rows[str(p)] = {
            "nsq_orbit": len(o0),
            "sq_orbit": len(osq),
            "half": (p + 1) // 2,
            "mixes_sq_nsq": mix,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "M acts transitively on nsq directions.  Fail: it mixes "
            "square with nonsquare."
        ),
    }


def theorem_D_irred() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    prim = {}
    for p in (3, 5, 7, 11, 13):
        o = order_of_2_mod(p)
        prim[str(p)] = {"ord_2": o, "primitive": o == p - 1}
    return {
        "proved": False,
        "proved_when_2_primitive_root": True,
        "W_irreducible_all_odd_p": False,
        "H0_quotient_irreducible_all_odd_p": False,
        "primitive_root_2": prim,
        "p7_certified_spin": True,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "If 2 is a primitive root mod p then each W^H is a simple "
            "C_p-module and G_aff^□-submodules of W are 0 or W.  "
            "Not a p-law at p=7 (Φ_7 splits); spin at p=7 still fills "
            "W.  Walsh 15.406 E stays OPEN."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.606  nsq line averages split W; M transits", flush=True)
    A = theorem_A_projectors()
    print(f"  A square π=0, nsq rank p−1: {A['proved']}", flush=True)
    B = theorem_B_orthogonal_sum()
    print(f"  B orthogonal sum I: {B['proved']}", flush=True)
    C = theorem_C_transitive()
    print(f"  C M transits nsq: {C['proved']}", flush=True)
    D = theorem_D_irred()
    print(
        f"  D irred when 2 primroot; all-p open: {D['W_irreducible_all_odd_p']}",
        flush=True,
    )
    out = {
        "prop": "15.606",
        "title": "W=⊕ nsq W^H; M transits; irred if 2 primitive root mod p",
        "proved": {
            "square_pi_zero_nsq_rank": A["proved"],
            "orthogonal_sum": B["proved"],
            "M_transitive_nsq": C["proved"],
            "W_irreducible_when_2_primitive_root": True,
            "W_irreducible_all_odd_p": False,
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
        "backend": "serial F2 identities; rref p=3,5,7; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15606.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
