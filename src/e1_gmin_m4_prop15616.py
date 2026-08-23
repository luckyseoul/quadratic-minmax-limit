#!/usr/bin/env python3
"""
Prop 15.616 — W2 via z+Dz; Walsh at p=11.
The old Krylov-gcd test was wrong: f(D)(z+Dz)≠0 for every
irred factor of g.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh for general p, W1 for p≡1, or leftover 2.

============================================================================
Setup.  z named halfspace-anti in U (15.613 A).  z+Dz is a
U-difference.  15.614: ε(z+Dz)=1 iff p≡3 (mod 4).  W_0 ≅ R =
F2[X]/h.  For irred f dividing g, w∈(f)R iff f(D)w=0.

============================================================================
Theorem A — PROVED (correct coprime test; p=5,7,11, the live g).
  Write w=z+Dz ∈ W_0.  For every irred f of g=(X^m+1)/(X+1),
      f(D)w ≠ 0.
  So w∉(f)R, hence I_U ⊄ (f)R, for each such f.  The earlier
  Krylov-first-dependence gcd with g was not the minpoly
  (it recovered a multiple of h).  Fail: Φ3(D)(z+Dz)=0 at p=5.  ∎

Theorem B — PROVED at p=11 (A + 15.612 CLASS + 15.614 B).
  p=11≡3 so W1 by 15.614.  A gives W2 (orbits {Φ3},{Φ5},
  {Φ15-pair}).  I_U is Aut-invariant and not contained in any
  maximal proper Aut-invariant ideal, so I_U=W_0.  Walsh 15.406 E
  at p=11.  Fail: some f(D)w=0 at p=11.  ∎

Theorem C — OPEN.  W2 as a p-law (A for every odd p); W1 for
  p≡1; Walsh for general p; residual_ii.  leftover+splus at
  p=5 k=20 remains 15.528, not a general close.

============================================================================
Backend: serial F2; apply f(D) p=5,7,11.  GPU unused.
Writes evidence/e1_gmin_m4_prop15616.json
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

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402
from e1_gmin_m4_prop15612 import _f2_divmod, _f2_factors  # noqa: E402
from e1_gmin_m4_prop15613 import _Dperm, named_z  # noqa: E402


def _apply_poly(w, poly, mul, gen, q):
    acc = np.zeros_like(w)
    cur = w.copy()
    for i, c in enumerate(poly):
        if c:
            acc ^= cur
        if i < len(poly) - 1:
            cur = _dil_fn(cur, mul, gen, q)
    return acc


def _g_factors(p: int):
    N = (p * p - 1) // 2
    m = N >> _v2(N)
    xm = [0] * m + [1]
    xm[0] = 1
    g, _ = _f2_divmod(xm, [1, 1])
    if m == 1:
        return g, []
    return g, _f2_factors(g)


def _w_zDz(p: int):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    Dperm = _Dperm(mul, gen, q)
    d = (bits ^ bits[Dperm]) & 1
    wfn = d[1 : 1 + q].copy()
    if d[0]:
        wfn ^= 1
    return wfn, q, mul, gen, eigen, inU


def theorem_A_coprime_test(primes=None) -> dict:
    if primes is None:
        primes = (5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        wfn, q, mul, gen, eigen, inU = _w_zDz(p)
        g, facs = _g_factors(p)
        hits = []
        all_nz = True
        for f in facs:
            r = _apply_poly(wfn, f, mul, gen, q)
            nz = bool(r.max())
            hits.append({"deg": len(f) - 1, "fD_nonzero": nz})
            all_nz = all_nz and nz
        ok = ok and eigen and inU and all_nz and len(facs) > 0
        rows[str(p)] = {
            "n_factors": len(facs),
            "all_fD_nonzero": all_nz,
            "hits": hits,
            "wt": int(wfn.sum()),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "f(D)(z+Dz)≠0 for every irred f of g.  Fail: Φ3(D)(z+Dz)=0 at p=5."
        ),
    }


def theorem_B_walsh_p11() -> dict:
    A = theorem_A_coprime_test((11,))
    from e1_gmin_m4_prop15614 import theorem_A_lift

    lift = theorem_A_lift((11,))
    e = lift["rows"]["11"]["eps"]
    w2 = A["proved"]
    w1 = e == 1
    return {
        "proved": bool(w1 and w2),
        "W1": w1,
        "W2": w2,
        "eps": e,
        "A11": A["rows"]["11"],
        "theorem": (
            "p=11: W1 by 15.614 and W2 by A, so I_U=W_0 (15.612 CLASS).  "
            "Walsh at p=11.  Fail: some f(D)(z+Dz)=0 at p=11."
        ),
    }


def theorem_C_open() -> dict:
    return {
        "proved": False,
        "W2_p_law": False,
        "W1_p_eq_1": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "W2 as a p-law, W1 for p≡1, Walsh ∀p, leftover 2 remain open.  "
            "Fail: Walsh from p=11 census alone; fail: residual_ii from Walsh."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.616  W2 via z+Dz; Walsh at p=11", flush=True)
    A = theorem_A_coprime_test()
    print(f"  A f(D)(z+Dz)≠0: {A['proved']}", flush=True)
    B = theorem_B_walsh_p11()
    print(f"  B Walsh p=11: {B['proved']}", flush=True)
    C = theorem_C_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.616",
        "title": "W2 via z+Dz; Walsh at p=11",
        "proved": {
            "fD_zDz_nonzero": A["proved"],
            "walsh_p11": B["proved"],
            "W2_p_law": False,
            "walsh_general_p": False,
            "W1_p_eq_1": False,
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
        "walsh_15_406_E": "OPEN (p=11 proved in this unit; not ∀p)",
        "backend": "serial F2 apply f(D); p=5,7,11; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15616.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
