#!/usr/bin/env python3
"""
Prop 15.608 — Möbius-plane type of F_p-sublines; two PSL-orbits;
1 ∈ dir(U).  I(H0)=H0 still not a p-law.  Walsh spanning OPEN.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close 15.406 E.

============================================================================
Setup.  Circles of the inversive plane of order p: F_p-sublines of
P¹(F_q), q=p², size p+1, 3-(q+1,p+1,1).  S-rows are square-direction
{∞}∪L.  15.598 B: χ_q|_{F_p^×}≡1, so P¹(F_p) is a square circle.

============================================================================
Theorem A — PROVED (change of basepoint; Max-free; all odd p).
  Let C be an F_p-subline, x∈C.  Send x to ∞; C\\{x} is an affine
  F_p-line of some direction b.  Then χ(b) is independent of x.
  Proof.  Reduce to x=∞ so C={∞}∪(a+F_p b).  For y=a+t0 b on C,
  z↦1/(z-y) sends y to ∞ and the line to F_p·b^{−1} (parametrically
  (t-t0)^{−1}/b), and χ(b^{−1})=χ(b).  Fail: χ flips with the
  basepoint.  ∎

  Call C square if χ=1, nonsquare if χ=−1.

Theorem B — PROVED (PSL orbits; Max-free; q=p², p odd).
  PGL(2,q) is sharply 3-transitive, so all F_p-sublines form a single
  PGL-orbit of P¹(F_p).  Setwise Stab_{PGL(2,q)}(P¹(F_p))=PGL(2,p):
  the unique Möbius sending three F_p-points as a given PGL(2,p) map
  is that PGL(2,p) element.  PGL(2,p)⊂PSL(2,q) (dets in F_p^× are
  squares in F_q, 15.598 B).  PSL is normal of index 2 in PGL, and
  Stab_PGL ⊂ PSL, so the PGL-orbit splits into exactly two PSL-orbits
  of equal size |PSL|/|PGL(2,p)|=p(p²+1)/2.
  PSL preserves χ-type: s,t∈PSL sending x and g(x) to ∞ exist by
  2-transitivity; tgs^{−1}∈PSL_∞ has the form z↦α²z+β and keeps χ.
  Both types exist (χ=±1 directions through ∞), so the two PSL-orbits
  are the square and nonsquare circles.  Fail: one PSL-orbit of all
  circles; fail: PGL(2,p) not in PSL(2,q).  ∎

Theorem C — PROVED (I in PSL).
  I(z)=1/z has det −1, a square in F_q (q≡1 mod 4).  I∈PSL permutes
  each orbit of B, and permutes the square 0-pencil (15.602 C).
  Fail: I sends a square circle to a nonsquare circle.  ∎

Theorem D — now 15.609 (p-law).  I(H0)=H0 for every odd p by the
  tangency lemma (opposite-type circles never meet in one point).
  This unit only certified it.

Theorem E — PROVED (Walsh lemma, not a close).
  U is closed under y↦−y (C_ij y_i y_j is even in y).  Antipodal
  F2-difference is 1, so 1∈dir(affine_span(U)) whenever U≠∅
  (15.598 E).  Walsh reduces to spanning V/⟨1⟩ with
  V=H0∩ker ℓ, dim N.  Spanning OPEN.  residual_ii stays False.

============================================================================
Backend: identities serial; type-check and I(H0) rref p=3,5,7.
GPU unused.  Writes evidence/e1_gmin_m4_prop15608.json
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
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15603 import direction_line_matrix  # noqa: E402
from walsh_linecode_rank import _mobius_perm, square_line_matrix  # noqa: E402


def _dir_chi(aff, add, mul, chi, p: int) -> int:
    a, b = aff[0], aff[1]
    d = b if a == 0 else add(b, mul(a, p - 1))
    return chi(d)


def circle_chis(p: int, v: np.ndarray, mul, add, chi, q: int) -> list[int]:
    """χ-type of a (p+1)-set after sending each of its points to ∞."""
    n = q + 1
    supp = [i for i in range(n) if v[i]]
    out = []
    for xidx in supp:
        if xidx == 0:
            aff = [i - 1 for i in supp if i != 0]
            out.append(_dir_chi(aff, add, mul, chi, p))
            continue
        alpha = xidx - 1
        D = 0 if alpha == 0 else mul(alpha, p - 1)
        pi = _mobius_perm(p, 0, 1, 1, D)
        w = np.zeros(n, dtype=np.uint8)
        for k in supp:
            w[pi[k]] = 1
        aff = [i - 1 for i in range(1, n) if w[i]]
        out.append(_dir_chi(aff, add, mul, chi, p))
    return out


def theorem_A_type_independent(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        S = square_line_matrix(p)
        Sp = direction_line_matrix(p, square=False)
        sq_bad = nsq_bad = 0
        for v in S:
            cs = circle_chis(p, v, mul, add, chi, q)
            if any(c != 1 for c in cs):
                sq_bad += 1
        for v in Sp:
            cs = circle_chis(p, v, mul, add, chi, q)
            if any(c != -1 for c in cs):
                nsq_bad += 1
        ok = ok and sq_bad == 0 and nsq_bad == 0
        ok = ok and S.shape[0] == Sp.shape[0] == p * (p + 1) // 2
        rows[str(p)] = {
            "sq_rows": int(S.shape[0]),
            "nsq_rows": int(Sp.shape[0]),
            "sq_type_flips": sq_bad,
            "nsq_type_flips": nsq_bad,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "χ-type of an F_p-subline is independent of the basepoint; "
            "changing basepoint replaces b by b^{−1}.  Fail: χ flips."
        ),
    }


def theorem_B_two_orbits() -> dict:
    # integer identities, plus χ_q on F_p
    primes = (3, 5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        q = p * p
        _, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        fp_sq = all(chi(t) == 1 for t in range(1, p))
        n_circ = p * (q + 1)
        n_psl = q * (q * q - 1) // 2
        n_pgl_p = p * (p * p - 1)
        orbit = n_psl // n_pgl_p
        n_sq = p * (q + 1) // 2
        ok = ok and fp_sq
        ok = ok and orbit == n_sq == n_circ // 2
        ok = ok and orbit != n_circ  # fail: one orbit
        rows[str(p)] = {
            "Fp_in_squares": fp_sq,
            "orbit": int(orbit),
            "n_square_circles": int(n_sq),
            "n_circles": int(n_circ),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Stab(P¹(F_p))=PGL(2,p)⊂PSL(2,q); two PSL-orbits, square "
            "vs nsq.  Fail: a single orbit of all circles."
        ),
    }


def theorem_C_I_preserves_type(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    A = theorem_A_type_independent(primes)
    ok = A["proved"]
    rows = {}
    for p in primes:
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        S = square_line_matrix(p)
        Sp = direction_line_matrix(p, square=False)
        Iperm = _mobius_perm(p, 0, 1, 1, 0)
        sq_flip = nsq_flip = 0
        pencil = 0
        for v in S:
            w = v[Iperm]
            cs = circle_chis(p, w, mul, add, chi, q)
            if any(c != 1 for c in cs):
                sq_flip += 1
            if v[1] == 1 and any(np.array_equal(w, r) for r in S):
                pencil += 1
        for v in Sp:
            w = v[Iperm]
            cs = circle_chis(p, w, mul, add, chi, q)
            if any(c != -1 for c in cs):
                nsq_flip += 1
        ok = ok and sq_flip == 0 and nsq_flip == 0
        ok = ok and pencil == (p + 1) // 2
        rows[str(p)] = {
            "I_sq_type_flips": sq_flip,
            "I_nsq_type_flips": nsq_flip,
            "pencil_to_row": pencil,
            "half": (p + 1) // 2,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "I∈PSL permutes square circles and the square 0-pencil.  "
            "Fail: I sends a square circle to a nonsquare circle."
        ),
    }


def theorem_D_I_H0_certified(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
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
        # off-0 image not a row, but in rowspan S, not S'
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
        }
    return {
        "proved": False,
        "certified": bool(ok),
        "H0_invariance_p_law": False,
        "rows": rows,
        "theorem": (
            "I(H0)=H0 certified at listed primes, not a p-law.  "
            "Fail: I maps a square row into rowspan(S')."
        ),
    }


def theorem_E_walsh_antipode() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": True,
        "one_in_dir_U": True,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "U is antipode-closed, so 1∈dir(affine_span(U)).  Walsh "
            "is spanning of V/⟨1⟩, V=H0∩ker ℓ.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.608  Möbius type; two PSL-orbits; 1 in dir(U)", flush=True)
    A = theorem_A_type_independent()
    print(f"  A type independent: {A['proved']}", flush=True)
    B = theorem_B_two_orbits()
    print(f"  B two PSL-orbits: {B['proved']}", flush=True)
    C = theorem_C_I_preserves_type()
    print(f"  C I preserves type: {C['proved']}", flush=True)
    D = theorem_D_I_H0_certified()
    print(f"  D I(H0)=H0 certified (not p-law): {D['certified']}", flush=True)
    E = theorem_E_walsh_antipode()
    print(f"  E 1 in dir(U); Walsh open: {E['walsh_general_p']}", flush=True)
    out = {
        "prop": "15.608",
        "title": "Square/nsq PSL-orbits of F_p-sublines; 1∈dir(U)",
        "proved": {
            "type_independent": A["proved"],
            "two_PSL_orbits": B["proved"],
            "I_preserves_type": C["proved"],
            "I_H0_p_law": False,
            "one_in_dir_U": True,
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "serial F2 identities; type-check p=3,5,7; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15608.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
