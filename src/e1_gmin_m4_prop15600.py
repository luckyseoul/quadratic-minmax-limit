#!/usr/bin/env python3
"""
Prop 15.600 — rank(S)=n/2 for every odd prime: the Gram gap of 15.599
is filled by the radical ⟨1⟩ of ker S.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Walsh spanning still OPEN.

============================================================================
Setup (15.598–15.599).  S is the b×n F2-incidence matrix of
square-direction {∞}∪L, n=p²+1, b=p(p+1)/2.  15.599 A: rank(SSᵀ)=n/2−1.
15.599 B: rank(S)≤n/2.  15.599 C: rank(S)∈{n/2−1, n/2}.

============================================================================
Theorem A — PROVED (all odd p; Max-free).
  Let K=ker S ⊂ F_2^n.  Then 1 ∈ K ∩ K^⊥, so dim(K ∩ K^⊥) ≥ 1.
    (i)  S1=0 because each row v_L has weight p+1 even.  Thus 1∈K.
    (ii) The sum of the p rows in any parallel class equals 1
         (15.599 B), so 1 lies in the row-span of S, i.e. in
         (ker S)^⊥ = K^⊥.
  Fail: claim S1=1 (needs p+1 odd).  ∎

Theorem B — PROVED (linear algebra over F2; all odd p).
  For any matrix S over a field, im(SSᵀ)=S((ker S)^⊥), hence
      rank(SSᵀ) = rank(S) − dim(K ∩ K^⊥).
  A gives dim(K ∩ K^⊥)≥1, so rank(SSᵀ) ≤ rank(S)−1.  15.599 A,B:
      n/2 − 1 = rank(SSᵀ) ≤ rank(S)−1  ⇒  rank(S) ≥ n/2,
      and rank(S) ≤ n/2, therefore rank(S)=n/2.
  In particular dim K = n/2 and dim(K ∩ K^⊥)=1, the radical is ⟨1⟩.
  Fail: claim rank(S)=n/2−1.  ∎

Theorem C — CERTIFIED p=3,5,7 (PSL generators; not a p-law).
  Random cyclic PΓL-modules of K=H0 are the full n/2-space (8/8 at
  each of p=3,5,7).  Aut_e remains reducible (15.599 F).  15.602 B:
  unique 1-dim G_aff^□-invariant subspace is ⟨1⟩ (theorem).  If
  H0/⟨1⟩ is irreducible then dir(affine_span(Max−))=H0.  That
  quotient irreducibility is OPEN.  Walsh still needs U to span
  H∩{ℓ=c}.

Theorem D — OPEN.  Walsh (15.406 E) is affine_span(U)=H∩{ℓ=c} with
  H=affine_span(Max−) ⊂ H0, dim H0=n/2 now a theorem.  residual_ii
  stays False.

============================================================================
Backend: identities + gf2_rref cross-check at p=5,7,11.
Writes evidence/e1_gmin_m4_prop15600.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15406 import gf2_nullspace, gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import _is_prime  # noqa: E402
from e1_gmin_m4_prop15599 import (  # noqa: E402
    class_relation_rank_upper,
    n_of,
    sst_rank_named,
)


def theorem_A_radical(primes=None) -> dict:
    """1 ∈ ker S ∩ (ker S)^⊥: S1=0 and class-sum=1, all odd p."""
    if primes is None:
        primes = [q for q in range(3, 80) if _is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        # (i) |v_L|=p+1 even
        wt = p + 1
        ok = ok and wt % 2 == 0
        ok = ok and wt != p  # fail: p odd would make S1=1 if wt odd
        # (ii) p lines in a class, p odd ⇒ class-sum includes ∞ once
        n_class = p
        ok = ok and n_class % 2 == 1
        rows[str(p)] = {"row_weight": wt, "class_size": n_class, "n": n_of(p)}
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": "1∈ker S ∩ (ker S)^⊥. Fail: S1=1 (p+1 odd).",
    }


def theorem_B_rank_eq(primes=None) -> dict:
    if primes is None:
        primes = [q for q in range(3, 80) if _is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        lo = sst_rank_named(p)  # n/2-1
        hi = class_relation_rank_upper(p)  # n/2
        # rank(S) ≥ lo+1 = n/2, and ≤ hi = n/2
        ok = ok and lo + 1 == hi
        ok = ok and hi == n_of(p) // 2
        # fail: rank = n/2-1 would mean lo == hi
        ok = ok and lo != hi
        rows[str(p)] = {"sst": int(lo), "upper": int(hi), "forced": int(hi)}
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": "rank(S)=n/2. Fail: rank(S)=n/2-1.",
    }


def theorem_B_crosscheck() -> dict:
    """gf2_rref at p=5,7,11 matches the forced n/2."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from walsh_linecode_rank import rank_one

    ok = True
    rows = {}
    for p in (5, 7, 11):
        rec = rank_one(p)
        ok = ok and rec["rank"] == rec["n_over_2"]
        ok = ok and rec["rank"] != rec["n_over_2"] - 1
        rows[str(p)] = rec
    return {"proved": bool(ok), "rows": rows}


def theorem_C_psl_certified() -> dict:
    # recorded from scripts/walsh_rank_equality.py: 8/8 full at p=3,5,7
    return {
        "proved": False,
        "certified": {"3": True, "5": True, "7": True},
        "aut_e_reducible": True,
        "theorem": (
            "PΓL-cyclic modules of H0 are full at p=3,5,7 (8/8). "
            "Not a general-p irreducibility proof."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "rank_S_equals_n_over_2": True,  # this prop
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "dim H0=n/2 is a theorem. Walsh is still U spanning the "
            "xor-hyperplane of H=affine_span(Max−)⊂H0. residual_ii False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.600  rank(S)=n/2 for every odd prime", flush=True)
    A = theorem_A_radical()
    print(f"  A radical ⟨1⟩: {A['proved']}", flush=True)
    B = theorem_B_rank_eq()
    print(f"  B rank=n/2: {B['proved']}", flush=True)
    Bx = theorem_B_crosscheck()
    print(f"  B rref p=5,7,11: {Bx['proved']}", flush=True)
    C = theorem_C_psl_certified()
    print(f"  C PSL certified (not ∀p): {C['certified']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.600",
        "title": "rank(S)=n/2 for every odd prime",
        "proved": {
            "one_in_radical": A["proved"],
            "rank_S_equals_n_over_2": B["proved"] and Bx["proved"],
            "psl_irreducible_all_p": False,
            "walsh_general_p": False,
        },
        "A": {k: v for k, v in A.items() if k != "rows"},
        "B": {k: v for k, v in B.items() if k != "rows"},
        "B_rref": Bx,
        "C": C,
        "D": D,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15600.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
