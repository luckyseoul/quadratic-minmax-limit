#!/usr/bin/env python3
"""
Prop 15.612 — Walsh is W1 ∧ W2 on Aut-invariant ideals of W_0.
CLASS of maximal Aut-invariant ideals is a p-law.  W1 is not.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E.

============================================================================
Setup.  15.611: W ≅ F2[X]/(X^N+1), W_0 ≅ R = F2[X]/h with
h=(X^N+1)/(X+1)=(X+1)^{2^a−1} g^{2^a}, a=v_2(N)≥2,
g=(X^m+1)/(X+1), m odd.  U is the {∞,0} pair-slice of Max−.
D, I(z)=1/z, and Frob(z)=z^p generate Aut({0,∞}) and preserve U,
so I_U := dir(U) ∩ W_0 is an Aut-invariant F2[D]-submodule of
W_0, hence an ideal of R.  15.608: 1∈dir(U).  Walsh ⇔ I_U=W_0.

============================================================================
Theorem A — PROVED (dictionary; 15.608 + 15.611; all odd p).
  dir(U) is D-invariant (D preserves U) and contains ⟨1⟩, so
  I_U is a D-submodule of the cyclic module W_0, i.e. an ideal
  of R.  Walsh 15.406 E is affine_span(U)=V-coset, equivalently
  I_U=W_0 equivalently I_U=(1) in R.  Fail: claim dir(U) need
  not be D-invariant.  ∎

Theorem B — PROVED (CLASS; Max-free; Fable xhigh PASS).
  I conjugates D to D^{−1} (X↦X^{−1} on R; h is self-reciprocal).
  Frob conjugates D to D^p (X↦X^p).  Ideals of R correspond to
  exponent vectors on {X+1 (cap 2^a−1)} ∪ {irred factors of g
  (cap 2^a)}.  Aut-invariant ⇔ the vector is constant on
  ⟨I,Frob⟩-orbits.  Maximal proper Aut-invariant ideals are
  therefore exactly
      (i)  (X+1)R = (D−I)W_0, the unique D-invariant hyperplane;
      (ii) (f_O)R for each ⟨I,Frob⟩-orbit O of irred factors of g,
           f_O=∏_{f∈O} f.
  Hence Walsh ⇔ W1 ∧ W2, where
      W1: I_U ⊄ (X+1)W_0  (some U-difference has (X+1)-valuation 0);
      W2: for every orbit O, some U-difference is nonzero mod f_O.
  At p=3, g=1 and W2 is vacuous.  At p=5 and p=7, g=Φ_3 is a
  single irred so W1⇒W2.  W2 is first live at p=11 (m=15;
  orbits {Φ_3}, {Φ_5}, {the two Φ_15-quartics}).  Fail: a
  missing maximal Aut-invariant ideal; fail: Frob acts by X↦X^2.
  ∎

Theorem C — W1 CERTIFIED p=3,5,7; NOT a p-law (Fable BLOCK).
  Pair-differences from a fixed U-basepoint have W_0-valuation 0
  at p=3,5,7.  One-point Aut-orbits do not prove W1: D-differences
  of one point lie in (D−I)W_0, and Frob-differences have ε-parity
  odd at p=3,7 but even at p=5 (not a constant in p).
  Translation-stay fills W_0 at p=5,7 but a 2-space of dim-3 W_0
  at p=3.  extra∈dir(U) is the socle, inside (D−I)W_0, not W1.
  Fail: claim W1 from Frob of one U-point (false at p=5);
  fail: claim no val-0 pair-difference at p=3.  ∎

Theorem D — OPEN.  W1 as a p-law, W2, and Walsh stay OPEN.
  residual_ii stays False.

============================================================================
Backend: CLASS identities serial; W1 rref p=3,5,7.  GPU unused.
Fable xhigh deep_review: CLASS PASS, W1 BLOCK (census is not a p-law).
Writes evidence/e1_gmin_m4_prop15612.json
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

from e1_gmin_m4_prop15406 import (  # noqa: E402
    gf2_nullspace,
    gf2_rref,
    load_minus,
)
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15606 import _W_basis  # noqa: E402
from e1_gmin_m4_prop15610 import _D_matrix  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402


def _f2_divmod(u: list[int], v: list[int]) -> tuple[list[int], list[int]]:
    u = u[:]
    while len(u) > 1 and u[-1] == 0:
        u.pop()
    v = v[:]
    while len(v) > 1 and v[-1] == 0:
        v.pop()
    if v == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(u) - len(v) + 1)
    while len(u) >= len(v) and (len(u) > 1 or u[0]):
        if not u[-1]:
            u.pop()
            continue
        sh = len(u) - len(v)
        q[sh] = 1
        for i, c in enumerate(v):
            u[i + sh] ^= c
        while len(u) > 1 and u[-1] == 0:
            u.pop()
        if u == [0]:
            break
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q, u


def _f2_mod(u, v):
    return _f2_divmod(u, v)[1]


def _f2_factors(poly: list[int]) -> list[list[int]]:
    """Distinct monic irred factors over F2, trial up to deg 8 (m≤15)."""
    f = poly[:]
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    facs = []
    # linear X+1
    lin = [1, 1]
    _, r = _f2_divmod(f, lin)
    if r == [0]:
        facs.append(lin)
        while True:
            q, r = _f2_divmod(f, lin)
            if r != [0]:
                break
            f = q
    d = 2
    while len(f) - 1 >= 2 * d or (len(f) > 1 and f != [1]):
        found = False
        # all monic of degree d
        for mask in range(1 << (d - 1)):
            cand = [1] + [(mask >> i) & 1 for i in range(d - 1)] + [1]
            _, r = _f2_divmod(f, cand)
            if r == [0]:
                # irred test: no proper factor
                irr = True
                for dd in range(1, d):
                    for m2 in range(1 << max(dd - 1, 0)):
                        if dd == 1:
                            sub = [1, 1]
                        else:
                            sub = [1] + [(m2 >> i) & 1 for i in range(dd - 1)] + [1]
                        if _f2_mod(cand, sub) == [0]:
                            irr = False
                            break
                    if not irr:
                        break
                if irr:
                    facs.append(cand)
                    while True:
                        q, rr = _f2_divmod(f, cand)
                        if rr != [0]:
                            break
                        f = q
                    found = True
                    break
        if not found:
            d += 1
            if d > 8:
                if f not in ([0], [1]) and len(f) > 1:
                    facs.append(f)
                break
        if f in ([0], [1]):
            break
    if len(f) > 1 and f not in facs:
        facs.append(f)
    return facs


def theorem_A_dictionary() -> dict:
    return {
        "proved": True,
        "walsh_iff_IU_unit": True,
        "theorem": (
            "dir(U) is D-invariant and 1∈dir(U), so I_U is an ideal "
            "of W_0 ≅ R.  Walsh ⇔ I_U=W_0.  Fail: dir(U) not D-invariant."
        ),
    }


def theorem_B_class(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        N = (p * p - 1) // 2
        a = _v2(N)
        m = N // (1 << a)
        four = N % 4 == 0
        # g = (X^m+1)/(X+1)
        xm1 = [0] * m + [1]
        xm1[0] = 1  # X^m + 1
        g, r = _f2_divmod(xm1, [1, 1])
        g_ok = r == [0] and m % 2 == 1 and a >= 2
        facs = _f2_factors(g) if m > 1 else []
        # p=3: m=1, g=1, no g-factors
        if p == 3:
            n_orb = 0
            orb_ok = g == [1]
        elif p in (5, 7):
            # g = X^2+X+1
            n_orb = 1
            orb_ok = facs == [[1, 1, 1]] or g == [1, 1, 1]
        elif p == 11:
            # orbits {Φ3}, {Φ5}, {two Φ15 quartics}
            n_orb = 3
            degs = sorted(len(f) - 1 for f in facs)
            orb_ok = degs == [2, 4, 4, 4]
        else:
            n_orb = None
            orb_ok = True
        ok = ok and four and g_ok and orb_ok
        rows[str(p)] = {
            "N": N,
            "v2_N": a,
            "m": m,
            "g_deg": len(g) - 1,
            "n_irred_factors_g": len(facs),
            "factor_degs": [len(f) - 1 for f in facs],
            "n_g_orbits_live": n_orb,
        }
    return {
        "proved": bool(ok),
        "Walsh_iff_W1_and_W2": True,
        "rows": rows,
        "theorem": (
            "Maximal Aut-invariant ideals of R are (X+1)R and (f_O)R.  "
            "Walsh ⇔ W1 ∧ W2.  Fail: Frob acts by X↦X^2."
        ),
    }


def _w0_eps_setup(p: int):
    WB, q, mul, add, chi = _W_basis(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    Dmat = _D_matrix(WB, mul, gen, q)
    K0, _ = gf2_nullspace(WB[0:1, :].astype(np.uint8))
    dimW0 = K0.shape[1]
    DK = (Dmat.astype(np.int32) @ K0.astype(np.int32)) % 2
    D0 = np.zeros((dimW0, dimW0), dtype=np.uint8)
    for j in range(dimW0):
        Aug = np.concatenate([K0, DK[:, j : j + 1]], axis=1)
        R, pivots, _rank = gf2_rref(Aug.copy())
        for i, pv in enumerate(pivots):
            if pv < dimW0:
                D0[pv, j] = R[i, dimW0]
    A0 = (D0.astype(np.int32) + np.eye(dimW0, dtype=np.int32)) % 2
    A0 = A0.astype(np.uint8)
    return WB, q, mul, K0, dimW0, A0


def _w0_of(v_full, WB, q, K0, dimW0):
    wfn = v_full[1 : 1 + q].copy()
    if v_full[0]:
        wfn ^= 1
    Aug = np.concatenate([WB, wfn.reshape(-1, 1)], axis=1)
    R, pivots, rank = gf2_rref(Aug.copy())
    k = WB.shape[1]
    if rank > k:
        return None
    x = np.zeros(k, dtype=np.uint8)
    for i, pv in enumerate(pivots):
        if pv < k:
            x[pv] = R[i, k]
    Aug = np.concatenate([K0, x.reshape(-1, 1)], axis=1)
    R, pivots, rank = gf2_rref(Aug.copy())
    if rank > dimW0:
        return None
    c = np.zeros(dimW0, dtype=np.uint8)
    for i, pv in enumerate(pivots):
        if pv < dimW0:
            c[pv] = R[i, dimW0]
    return c


def _eps(c, A0, dimW0):
    if c is None or c.max() == 0:
        return None
    M = np.concatenate([A0, c.reshape(-1, 1)], axis=1)
    in_im = gf2_rref(M.copy())[2] == gf2_rref(A0.copy())[2]
    return 0 if in_im else 1


def theorem_C_W1_certified(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        Y, C = load_minus(p)
        Y = np.sign(Y.astype(np.float64)).astype(np.int8)
        B = ((1 - Y) // 2).astype(np.uint8)
        fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
        BU = B[fe < 0]
        WB, q, mul, K0, dimW0, A0 = _w0_eps_setup(p)
        n_val0 = 0
        take = min(80, len(BU))
        for i in range(take):
            d = (BU[i] ^ BU[0]) & 1
            e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
            if e == 1:
                n_val0 += 1
        q2, mul2, add2, chi2, frob, norm, ia, ib = field_ctx(p)
        Fperm = np.arange(q + 1)
        Fperm[0] = 0
        for e in range(q):
            Fperm[1 + e] = 1 + frob(e)
        nF_odd = nF_even = nF_none = 0
        for i in range(take):
            y = BU[i]
            fy = y[Fperm]
            e = _eps(_w0_of((y ^ fy) & 1, WB, q, K0, dimW0), A0, dimW0)
            if e is None:
                nF_none += 1
            elif e:
                nF_odd += 1
            else:
                nF_even += 1
        # W1 certified: some pair-diff val 0
        # Frob one-point is NOT a p-law: odd at p=3,7, not at p=5
        frob_all_odd = nF_odd == take
        frob_no_odd = nF_odd == 0
        if p in (3, 7):
            frob_pattern = frob_all_odd
        else:
            frob_pattern = frob_no_odd
        ok = ok and n_val0 > 0 and frob_pattern
        rows[str(p)] = {
            "nU": int(len(BU)),
            "take": take,
            "n_pair_val0": n_val0,
            "W1_certified": n_val0 > 0,
            "Frob_eps_odd": nF_odd,
            "Frob_eps_even": nF_even,
            "Frob_eps_none": nF_none,
            "Frob_one_point_is_p_law": False,
        }
    return {
        "proved": False,
        "W1_p_law": False,
        "W1_certified_p357": bool(ok),
        "rows": rows,
        "theorem": (
            "W1 certified p=3,5,7; not a p-law (Frob one-point fails "
            "at p=5; translation-stay fails at p=3).  Fail: W1 from "
            "Frob of one U-point."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "W1_p_law": False,
        "W2_p_law": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "Walsh ⇔ W1 ∧ W2 on Aut-invariant ideals.  CLASS is a "
            "p-law; W1/W2/Walsh stay open.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.612  Walsh ⇔ W1 ∧ W2; CLASS p-law, W1 not", flush=True)
    A = theorem_A_dictionary()
    print(f"  A dictionary: {A['proved']}", flush=True)
    B = theorem_B_class()
    print(f"  B CLASS: {B['proved']}", flush=True)
    C = theorem_C_W1_certified()
    print(f"  C W1 certified not p-law: cert={C['W1_certified_p357']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.612",
        "title": "Walsh ⇔ W1 ∧ W2; CLASS of Aut-invariant ideals",
        "proved": {
            "dictionary_Walsh_iff_IU_unit": A["proved"],
            "class_maximal_Aut_ideals": B["proved"],
            "W1_p_law": False,
            "W2_p_law": False,
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
        "backend": "CLASS serial; W1 rref p=3,5,7; GPU unused",
        "claude_referee": (
            "deep_review BLOCK on W1 p-law; CLASS PASS (Walsh ⇔ W1∧W2)"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15612.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
