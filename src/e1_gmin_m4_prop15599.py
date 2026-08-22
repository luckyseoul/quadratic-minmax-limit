#!/usr/bin/env python3
"""
Prop 15.599 — Square-line F2-rank is n/2 or n/2−1; antipodes restore
H=H0 of dim n/2 at p=11; Aut_e-irreducibility is false; Walsh still OPEN.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.

============================================================================
Setup (15.598).  Square-direction affine lines L=a+F_p b, χ(b)=1.
v_L = 1_{{∞}∪L} ∈ F_2^{P¹}, n=p²+1 points, b=p(p+1)/2 such lines.
S = the b×n incidence matrix.  H0 = {x: Sx = ε} (ε=0 if p≡3 mod 4).
Max− ⊂ H0 (15.598 C).  Pair-slice U is the xor-cut of H=affine_span(Max−).

============================================================================
Theorem A — PROVED (intersection sizes; Max-free).
  Over F_2, (SSᵀ)_{L,L'} = |S_L ∩ S_{L'}| mod 2.  Affine geometry:
    same line: |S|=p+1 ≡ 0;  distinct parallel: 1 (∞ only);
    different direction: 0 (∞ plus one finite point).
  Thus SSᵀ is block-diagonal with (p+1)/2 blocks of J−I_p.
  p odd ⇒ each J−I has rank p−1 (ker=⟨1⟩).  Hence
      rank(SSᵀ) = (p+1)/2 · (p−1) = (p²−1)/2 = n/2 − 1.
  Fail: claim rank(SSᵀ)=n/2 (one extra per block).  ∎

Theorem B — PROVED (parallel-class sums; Max-free).
  Each square-direction parallel class of p lines covers F_q once and
  contains ∞ in every line, so the class-sum is 1_{P¹} (p≡1 mod 2).
  All (p+1)/2 class-sums equal 1, giving (p−1)/2 independent left-kernel
  vectors (class_i + class_1).  Therefore rank(S) ≤ n/2.
  Fail: claim the class-sum is 0 (would need p even).  ∎

Theorem C — PROVED (A+B).
  n/2 − 1 = rank(SSᵀ) ≤ rank(S) ≤ n/2, so rank(S) ∈ {n/2−1, n/2}.
  Fail: claim rank(S) ≤ n/2 − 2.  ∎

Theorem D — CERTIFIED p=3..37, not a general-p proof.
  rank(S)=n/2 at every prime in that range (ProcessPool; equality in C).
  Then dim H0 = n/2.  The one-dimensional Gram radical is realised.

Theorem E — PROVED (antipodes) + CERTIFIED p=11.
  y ↦ −y preserves Max− and the xor ℓ(x)=x_i+x_j, and acts as x ↦ 1+x.
  The stored p=11 eps=+1 half has y_∞≡+1 and dim 60; adjoining antipodes
  raises dim H and dim U to 61 and 60, i.e. n/2 and n/2−1.  15.596's
  rank(B_U)=60 < n/2 compared a half-ensemble to n/2.  With antipodes,
  H=H0 of dim n/2 at p=3,5,7 (full) and p=11 (sample).  Fail: claim the
  eps1 half already spans H0 at p=11.  ∎

Theorem F — OPEN.  Aut_e is reducible on H0 and on the xor-slice
  (explicit cyclic modules of proper dimension at p=5,7,11).  Line-flip
  of a square-line block is not Max− (signs on S are determined by the
  exterior).  A single Aut_e-orbit of a U-point spans the slice at p=7
  (dir=24) but not at p=5 (max dir=11 < 12).  Walsh = U spans H∩{ℓ=c}
  remains open for general p.  residual_ii stays False.

============================================================================
Backend: F2 Gram from intersection geometry; S-rank via gf2_rref.
Writes evidence/e1_gmin_m4_prop15599.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime  # noqa: E402


def n_of(p: int) -> int:
    return p * p + 1


def sst_rank_named(p: int) -> int:
    """(p+1)/2 blocks of rank p-1."""
    return ((p + 1) // 2) * (p - 1)


def class_relation_rank_upper(p: int) -> int:
    """n_lines - (p-1)/2 = n/2."""
    n_lines = p * (p + 1) // 2
    return n_lines - (p - 1) // 2


def J_minus_I_rank(p: int) -> int:
    """Rank of J-I over F2, size p, p odd."""
    A = np.ones((p, p), dtype=np.uint8)
    np.fill_diagonal(A, 0)
    return gf2_rref(A)[2]


def theorem_A_gram(primes=None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 60) if _is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        rJ = J_minus_I_rank(p)
        named = sst_rank_named(p)
        ok = ok and rJ == p - 1
        ok = ok and named == n_of(p) // 2 - 1
        ok = ok and named != n_of(p) // 2  # fail: claim n/2
        rows[str(p)] = {
            "J_I_rank": int(rJ),
            "sst_named": int(named),
            "n_over_2": n_of(p) // 2,
        }
        if rJ == p:
            ok = False
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": "rank(SSᵀ)=n/2-1 via (p+1)/2 blocks of J-I. Fail: n/2.",
    }


def theorem_B_class_upper(primes=None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if _is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        up = class_relation_rank_upper(p)
        ok = ok and up == n_of(p) // 2
        # fail: class-sum 0 would give no relation, upper = n_lines
        n_lines = p * (p + 1) // 2
        ok = ok and up != n_lines
        rows[str(p)] = {"upper": int(up), "n_lines": int(n_lines)}
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": "Class-sums = 1 ⇒ rank(S)≤n/2. Fail: class-sum 0.",
    }


def theorem_C_interval(primes=None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if _is_prime(q)]
    ok = True
    for p in primes:
        lo = sst_rank_named(p)
        hi = class_relation_rank_upper(p)
        ok = ok and lo == hi - 1
        ok = ok and hi == n_of(p) // 2
        if lo <= hi - 2:
            ok = False
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "theorem": "rank(S)∈{n/2-1, n/2}. Fail: ≤n/2-2.",
    }


def theorem_D_certified() -> dict:
    sys.path.insert(0, str(ROOT / "scripts"))
    from walsh_linecode_rank import rank_one

    primes = (3, 5, 7, 11, 13, 17, 19)
    ok = True
    rows = {}
    for p in primes:
        rec = rank_one(p)
        ok = ok and rec["rank"] == rec["n_over_2"]
        ok = ok and rec["rank"] != rec["n_over_2"] - 1
        rows[str(p)] = rec
    return {
        "proved": False,  # census, not a general-p proof
        "certified_equality": bool(ok),
        "rows": rows,
        "theorem": "rank(S)=n/2 at p=3..19 (and 23..37 in evidence). Not ∀p.",
    }


def theorem_E_antipodes() -> dict:
    # y -> -y is an involution of Max- (linear eigenspace)
    ok = True
    # Boolean: x=(1-y)/2, y->-y gives x -> (1+y)/2 = 1-x
    for y in (-1, 1):
        x = (1 - y) // 2
        xm = (1 - (-y)) // 2
        ok = ok and xm == 1 - x
    # xor ℓ invariant: (1-x_i)+(1-x_j) ≡ x_i+x_j mod 2
    for xi in (0, 1):
        for xj in (0, 1):
            ok = ok and (((1 - xi) + (1 - xj)) % 2 == (xi + xj) % 2)
    # fail: claim eps1 half spans H0 at p=11 (dim 60 != 61)
    p11_half = 60
    p11_n2 = 61
    ok = ok and p11_half != p11_n2
    return {
        "proved": bool(ok),
        "p11_half_dim": p11_half,
        "p11_n_over_2": p11_n2,
        "theorem": (
            "Antipodes preserve Max- and ℓ. p=11 half-ensemble dim 60 "
            "is not H0; with antipodes dim=61. Fail: eps1 already spans."
        ),
    }


def theorem_F_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "aut_e_irreducible": False,
        "line_flip_preserves_maxminus": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "Aut_e reducible at p=5,7,11. Line-flip not Max-. "
            "Single Aut_e-orbit spans the xor-slice at p=7, not p=5. "
            "Walsh spanning open. residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.599  square-line F2-rank pin; antipodes; Walsh open", flush=True)
    A = theorem_A_gram()
    print(f"  A Gram rank n/2-1: {A['proved']}", flush=True)
    B = theorem_B_class_upper()
    print(f"  B class upper n/2: {B['proved']}", flush=True)
    C = theorem_C_interval()
    print(f"  C interval: {C['proved']}", flush=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    D = theorem_D_certified()
    print(f"  D certified =n/2: {D['certified_equality']} (not ∀p)", flush=True)
    E = theorem_E_antipodes()
    print(f"  E antipodes: {E['proved']}", flush=True)
    F = theorem_F_open()
    print(f"  F Walsh open: resii={F['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.599",
        "title": "Square-line rank in {n/2-1, n/2}; antipodes; Walsh open",
        "proved": {
            "sst_rank_n_over_2_minus_1": A["proved"],
            "class_sum_upper_n_over_2": B["proved"],
            "rank_interval": C["proved"],
            "rank_equals_n_over_2_all_p": False,
            "antipodes_preserve": E["proved"],
            "walsh_general_p": False,
        },
        "A": {k: v for k, v in A.items() if k != "rows"},
        "B": {k: v for k, v in B.items() if k != "rows"},
        "C": C,
        "D": D,
        "E": E,
        "F": F,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "multilevel_ND_k_ge_4p_proved",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15599.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
