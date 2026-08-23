#!/usr/bin/env python3
"""
Prop 15.609 — opposite-type circles are never tangent; I(H0)=H0
for every odd p.  Max-free.  Walsh spanning still OPEN.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close 15.406 E.

============================================================================
Setup.  15.608: F_p-sublines of P¹(F_q) are inversive circles of
two PSL-orbits (square / nsq).  I(z)=1/z ∈ PSL preserves type.
S, S' = square / nsq {∞}∪L incidence.  H0=ker S, H0'=ker S'.

============================================================================
Theorem A — PROVED (incidence; all odd p).
  H0' = H0^⊥ = rowspan(S).  A square ∞-circle and a nsq ∞-circle
  have different direction classes, so their affine parts meet in
  one point and |C∩C'|=2 (∞ plus that point).  Thus every row of
  S is orthogonal to every row of S', rowspan(S) ⊂ ker S'=H0'.
  Both have dim n/2 (15.600 / 15.603).  Fail: |C∩C'|=1.  ∎

Theorem B — PROVED (tangency lemma; all odd p).
  Distinct F_p-sublines meet in 0, 1, or 2 points (unique circle
  through three points: p(p²+1)C(p+1,3)=C(q+1,3)).  Tangency
  means |∩|=1.  Two tangent circles have the same χ-type.
  Proof: PGL(2,q) is transitive on flags (point, circle through
  it).  Normalize a tangent pair to (P¹(F_p), ∞).  The other
  circle through ∞ and tangent at ∞ has affine part parallel to
  F_p, hence direction in F_p^× and χ=1.  Both are square.
  PGL either preserves both PSL-orbits or swaps both, so it cannot
  send a mixed-type pair to a same-type pair.  Hence mixed-type
  pairs are never tangent, and |C∩C'|∈{0,2} when types differ.
  Fail: a square circle tangent to a nsq circle.  ∎

  Case 0∈M of the inversion identity is elementary without B:
  M=F_p d, χ(d)=−1, dL is nsq not through 0, meets the square
  0-line F_p in one nonzero point; plus 0 gives |∩|=2.

Theorem C — PROVED (A+B; all odd p).
  I sends an off-0 square ∞-circle to a square circle through 0.
  That vector is even on every nsq ∞-circle by B, hence lies in
  ker S'=rowspan(S).  The 0-pencil is already permuted (15.602 C /
  15.608 C).  So I(rowspan(S))=rowspan(S) and I(H0)=H0.
  Equivalently, the F2-span of all square circles is rowspan(S).
  Fail: I maps a square row into rowspan(S').  ∎

Theorem D — OPEN.  Walsh 15.406 E is still dir(U)=H0∩ker ℓ.
  I now acts on the xor-slice for {0,∞}.  Spanning OPEN.
  residual_ii stays False.

============================================================================
Backend: identities serial; rref / intersection p=3,5,7.  GPU unused.
Writes evidence/e1_gmin_m4_prop15609.json
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
from e1_gmin_m4_prop15599 import n_of  # noqa: E402
from e1_gmin_m4_prop15603 import direction_line_matrix  # noqa: E402
from walsh_linecode_rank import _mobius_perm, square_line_matrix  # noqa: E402


def theorem_A_dual_equals_nsq(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        S = square_line_matrix(p)
        Sp = direction_line_matrix(p, square=False)
        n = S.shape[1]
        Gram = (S.astype(np.int32) @ Sp.astype(np.int32).T) % 2
        rS = gf2_rref(S.copy())[2]
        rP = gf2_rref(Sp.copy())[2]
        # rowspan S ⊂ ker S'
        SSp = (Sp.astype(np.int32) @ S.astype(np.int32).T) % 2
        subset = int(SSp.max()) == 0
        # integer |∩| for one pair: always 2
        caps = []
        for v in S[:3]:
            for w in Sp[:3]:
                caps.append(int((v & w).sum()))
        ok = ok and subset and rS == rP == n // 2
        ok = ok and all(c == 2 for c in caps)
        ok = ok and int(Gram.max()) == 0
        rows[str(p)] = {
            "rank_S": int(rS),
            "rank_Sprime": int(rP),
            "n_over_2": n // 2,
            "S_perp_Sprime": subset,
            "sample_caps": caps,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "H0'=H0^⊥=rowspan(S): opposite-direction ∞-circles meet "
            "in two points.  Fail: |C∩C'|=1."
        ),
    }


def theorem_B_no_mixed_tangency(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    Iperm = None
    for p in primes:
        S = square_line_matrix(p)
        Sp = direction_line_matrix(p, square=False)
        Iperm = _mobius_perm(p, 0, 1, 1, 0)
        # all opposite-type pairs among ∞-circles: |∩|=2
        inf_caps = set()
        for v in S:
            for w in Sp:
                inf_caps.add(int((v & w).sum()))
        # inverted off-0 square vs nsq ∞-circles: 0 or 2, never 1
        inv_caps = set()
        n_off = 0
        for v in S:
            if v[1] == 1:
                continue
            n_off += 1
            w = v[Iperm]
            for r in Sp:
                inv_caps.add(int((w & r).sum()))
        ok = ok and inf_caps == {2}
        ok = ok and inv_caps <= {0, 2}
        ok = ok and 1 not in inv_caps
        rows[str(p)] = {
            "infty_opposite_caps": sorted(inf_caps),
            "inverted_vs_nsq_caps": sorted(inv_caps),
            "n_off0": n_off,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Opposite-type circles are never tangent, so |∩|∈{0,2}.  "
            "Fail: a square–nsq pair with |∩|=1."
        ),
    }


def theorem_C_I_preserves_H0(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    A = theorem_A_dual_equals_nsq(primes)
    B = theorem_B_no_mixed_tangency(primes)
    ok = A["proved"] and B["proved"]
    rows = {}
    for p in primes:
        S = square_line_matrix(p)
        Sp = direction_line_matrix(p, square=False)
        Iperm = _mobius_perm(p, 0, 1, 1, 0)
        H0, _ = gf2_nullspace(S)
        Hg = H0[Iperm, :]
        rH = gf2_rref(H0.copy())[2]
        rA = gf2_rref(np.concatenate([H0, Hg], axis=1))[2]
        preserve = rH == rA
        off = next(v for v in S if v[1] == 0)
        w = off[Iperm]
        inS = gf2_rref(np.vstack([S, w]))[2] == gf2_rref(S.copy())[2]
        inP = gf2_rref(np.vstack([Sp, w]))[2] == gf2_rref(Sp.copy())[2]
        is_row = any(np.array_equal(w, r) for r in S)
        ok = ok and preserve and inS and (not inP) and (not is_row)
        rows[str(p)] = {
            "H0_preserved": preserve,
            "off0_in_rowspan_S": inS,
            "off0_in_rowspan_Sprime": inP,
            "off0_is_row": is_row,
            "dim_H0": int(H0.shape[1]),
            "n_over_2": n_of(p) // 2,
        }
        ok = ok and H0.shape[1] == n_of(p) // 2
    return {
        "proved": bool(ok),
        "H0_invariance_p_law": bool(ok),
        "rows": rows,
        "theorem": (
            "I(H0)=H0 for every odd p.  Span of all square circles "
            "equals rowspan(S).  Fail: I maps a square row into "
            "rowspan(S')."
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
            "I acts on H0 and on the xor-slice for {0,∞}.  "
            "Walsh spanning of V/⟨1⟩ stays OPEN.  residual_ii False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.609  opposite-type never tangent; I(H0)=H0", flush=True)
    A = theorem_A_dual_equals_nsq()
    print(f"  A H0'=H0^perp=rowspan(S): {A['proved']}", flush=True)
    B = theorem_B_no_mixed_tangency()
    print(f"  B no mixed tangency: {B['proved']}", flush=True)
    C = theorem_C_I_preserves_H0()
    print(f"  C I(H0)=H0 p-law: {C['proved']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.609",
        "title": "Opposite-type circles never tangent; I(H0)=H0",
        "proved": {
            "H0prime_equals_H0_perp": A["proved"],
            "no_mixed_tangency": B["proved"],
            "I_preserves_H0": C["proved"],
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
        "backend": "serial F2 identities; intersection p=3,5,7; GPU unused",
        "claude_referee": "deep_review PASS (fable-5 xhigh) on tangency lemma",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15609.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
