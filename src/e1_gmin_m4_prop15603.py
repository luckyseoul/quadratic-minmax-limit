#!/usr/bin/env python3
"""
Prop 15.603 — square and nonsquare duals of the ∞-line code:
H0 ∩ H0' = ⟨1⟩ and H0 + H0' = even-weight code.  All odd p, Max-free.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** prove H0/⟨1⟩ irreducible or Walsh 15.406 E.

============================================================================
Setup.  S (resp. S') is the b×n F2-incidence of square-direction
(resp. nonsquare-direction) sets {∞}∪L, L an affine F_p-line in F_q,
q=p², n=q+1, b=p(p+1)/2.  H0=ker S, H0'=ker S'.  15.600: rank(S)=n/2,
dim H0=n/2, radical ⟨1⟩.  Even-weight code E={x: ⟨1,x⟩=0}, dim n−1.

============================================================================
Theorem A — PROVED (all odd p; same as 15.600 for S').
  rank(S')=n/2.  Rows have weight p+1 even, so S'1=0.  The p lines of
  a nonsquare parallel class plus ∞ sum to 1 (p odd ⇒ p·e_∞ + 1_V = 1).
  Thus 1 ∈ ker S' ∩ (ker S')^⊥, and rank(S'S'ᵀ)=n/2−1 as for S
  (Gram still J−I on each parallel class).  Fail: rank(S')=n/2−1.  ∎

Theorem B — PROVED (AG(2,p) line sums; all odd p).
  H0 ∩ H0' = ⟨1⟩.  If x lies in both kernels then
      x_∞ + ∑_{u∈L} x_u = 0
  for every affine F_p-line L (square and nonsquare).  Write σ=x_∞
  and f=x on V=AG(2,p).  Every affine line has sum σ.
  Fix a∈V.  The p+1 lines through a: left side (p+1)σ=0 (p+1 even).
  Right side: (p+1)f(a) + ∑_{b≠a} f(b) = p f(a) + ∑_V f, and p odd
  so f(a)+∑_V f = 0 for every a.  Thus f is constant on V, equal to
  Σ=∑_V f.  Line sums then force σ=Σ, so x is constant on P¹.
  Fail: claim a non-constant common kernel vector; fail: p even
  (the factor p vanishes).  ∎

Theorem C — PROVED (dimension).
  H0 + H0' = E.  Both sit in E (1 annihilates both).  Theorem B:
  dim(H0 ∩ H0')=1, so dim(H0+H0') = n/2+n/2−1 = n−1 = dim E.  Fail:
  H0=H0' (then dim intersection n/2>1).  ∎

  So the even-weight PSL-heart Q = E/⟨1⟩ splits as
      (H0/⟨1⟩) ⊕ (H0'/⟨1⟩)
  of dimensions (q−1)/2 each.  Mortimer: the F2-heart of PSL(2,q) on q+1 points is not
  absolutely simple for q odd; these are the two extra
  constituents (q=p²≡1 mod 8).  Irreducibility of each
  summand is OPEN.

Theorem D — OPEN.  H0/⟨1⟩ irreducible would give
  dir(affine_span(Max−))=H0.  Walsh 15.406 E still requires the
  xor-slice.  residual_ii stays False.

============================================================================
Backend: identities serial; rref cross-check p=3,5,7,11. GPU unused.
Writes evidence/e1_gmin_m4_prop15603.json
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
from e1_gmin_m4_prop15598 import _is_prime, field_ctx  # noqa: E402
from e1_gmin_m4_prop15599 import n_of  # noqa: E402
from walsh_linecode_rank import square_line_matrix  # noqa: E402


def direction_line_matrix(p: int, square: bool) -> np.ndarray:
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1
    used = set()
    dirs = []
    for b in range(1, q):
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    want = 1 if square else -1
    rows = []
    for b in dirs:
        if chi(b) != want:
            continue
        covered = set()
        for a in range(q):
            if a in covered:
                continue
            v = np.zeros(n, dtype=np.uint8)
            v[0] = 1
            for t in range(p):
                e = add(a, mul(t, b))
                v[1 + e] = 1
                covered.add(e)
            rows.append(v)
    return np.stack(rows, axis=0)


def theorem_A_nonsquare_rank(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        S1 = square_line_matrix(p)
        S0 = direction_line_matrix(p, square=False)
        n = S0.shape[1]
        r0 = gf2_rref(S0)[2]
        r1 = gf2_rref(S1)[2]
        one = np.ones(n, dtype=np.uint8)
        S01 = (S0.astype(np.int32) @ one.astype(np.int32)) % 2
        # class-sum: p rows of one parallel class = 1
        # (checked as 1 in rowspan)
        r_one = gf2_rref(np.vstack([S0, one]))[2]
        ok = ok and r0 == n // 2 == r1
        ok = ok and int(S01.max()) == 0
        ok = ok and r_one == r0  # 1 already in rowspan
        ok = ok and S0.shape[0] == S1.shape[0] == p * (p + 1) // 2
        rows[str(p)] = {
            "rank_S": int(r1),
            "rank_Sprime": int(r0),
            "n_over_2": n // 2,
            "one_in_ker": bool(S01.max() == 0),
            "one_in_rowspan": r_one == r0,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "rank(S')=n/2 by the same radical ⟨1⟩ as 15.600.  "
            "Fail: rank(S')=n/2−1."
        ),
    }


def theorem_B_intersection(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        S1 = square_line_matrix(p)
        S0 = direction_line_matrix(p, square=False)
        n = S1.shape[1]
        H0, _ = gf2_nullspace(S1)
        Hp, _ = gf2_nullspace(S0)
        # intersection of column spaces: vectors H0 t = Hp s
        # stack [H0 | Hp] nullspace projects
        M = np.concatenate([H0, Hp], axis=1)
        K, _ = gf2_nullspace(M)
        # K is (k0+kp) x d; H0 part is first k0 rows of... wait
        # M x = 0 with x=(t,s), H0 t + Hp s = 0 so H0 t = Hp s (F2).
        # dim of intersection = dim{H0 t : exists s, H0 t = Hp s}
        # = dim ker of the map, more carefully: rank(M)=dim(H0+H0'),
        # dim int = dim H0 + dim H0' - rank(M).
        rM = gf2_rref(M)[2]
        dim_int = H0.shape[1] + Hp.shape[1] - rM
        one = np.ones(n, dtype=np.uint8)
        # 1 in both
        in0 = bool(((S1.astype(np.int32) @ one.astype(np.int32)) % 2).max() == 0)
        in1 = bool(((S0.astype(np.int32) @ one.astype(np.int32)) % 2).max() == 0)
        ok = ok and dim_int == 1 and in0 and in1
        ok = ok and H0.shape[1] == n // 2 == Hp.shape[1]
        rows[str(p)] = {
            "dim_H0": int(H0.shape[1]),
            "dim_H0p": int(Hp.shape[1]),
            "dim_int": int(dim_int),
            "rank_sum": int(rM),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "H0 ∩ H0' = ⟨1⟩: all affine line sums equal ⇒ f constant "
            "on AG(2,p).  Fail: a non-constant common kernel vector."
        ),
    }


def theorem_C_sum_even(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    B = theorem_B_intersection(primes)
    ok = B["proved"]
    rows = {}
    for p, rec in B["rows"].items():
        n = n_of(int(p))
        dim_sum = rec["rank_sum"]
        ok = ok and dim_sum == n - 1
        ok = ok and rec["dim_int"] == 1
        rows[p] = {
            "dim_sum": dim_sum,
            "dim_even": n - 1,
            "not_equal": rec["dim_H0"] != rec["dim_int"],
        }
        ok = ok and rec["dim_H0"] != rec["dim_int"]
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "H0+H0'=even-weight code.  Fail: H0=H0'.  The heart "
            "E/⟨1⟩ = (H0/⟨1⟩)⊕(H0'/⟨1⟩)."
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
            "Irreducibility of H0/⟨1⟩ would give Max− spanning of H0, "
            "not Walsh.  15.406 E and residual_ii stay OPEN/False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.603  square/nonsquare duals of the ∞-line code", flush=True)
    A = theorem_A_nonsquare_rank()
    print(f"  A rank(S')=n/2: {A['proved']}", flush=True)
    B = theorem_B_intersection()
    print(f"  B H0∩H0'=⟨1⟩: {B['proved']}", flush=True)
    C = theorem_C_sum_even()
    print(f"  C H0+H0'=even: {C['proved']}", flush=True)
    D = theorem_D_open()
    print(f"  D irred/Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.603",
        "title": "H0 ∩ H0' = ⟨1⟩ and H0 + H0' = even-weight",
        "proved": {
            "rank_Sprime": A["proved"],
            "intersection_ones": B["proved"],
            "sum_even_weight": C["proved"],
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
        "backend": "serial F2 identities; rref p=3,5,7,11; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15603.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
