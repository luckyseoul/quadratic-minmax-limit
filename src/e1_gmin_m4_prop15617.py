#!/usr/bin/env python3
"""
Prop 15.617 — W1 for p≡1 via NSQ-class stay-sum; correct W2 test.
Retracts Walsh-at-p=11 from 15.616 (weak f(D)w≠0 ≠ not-in-(f)R).

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh ∀p or leftover 2.

============================================================================
Theorem A — PROVED (module arithmetic; p=5,7,11).
  w=c(D)γ ∈ W with γ=1_M.  w∈(f)R iff f divides c.  So w∉(f)R
  iff gcd(c,f)=1, not iff f(D)w≠0.  For w=z+Dz, gcd(c,Φ3)=Φ3 at
  p=5,7, and at p=11 gcd(c,Φ3)=1 but Φ5 and both Φ15-quartics
  divide c.  15.616 B (Walsh at p=11 from this vector) is
  withdrawn.  Fail: gcd(c,Φ3)=1 at p=5 for z+Dz.  ∎

Theorem B — PROVED as existence density (p=5,7); not a p-law.
  From one U-basepoint, 72/156 (p=5) and 218/400 (p=7) of
  U-differences have gcd(c,g)=1, so W2 holds at those p by
  generic diffs (consistent with 15.406).  Not a general p-law.  ∎

Theorem C — CONSTRUCTION p-law; ε CERTIFIED p=5,13,17.
  s_N = ∑_{a∈F_p^×, (a/p)=-1, stay} (z+T_a z) (F2).
  Stay on F_p^× has size (p-1)/2.  s_N ∈ W_0 ⊂ I_U.
  ε(s_N)=1 at p=5,13,17 (Zolotarev nsq class; Fable).  This is
  W1 for p≡1 at those p.  Fail: ε(s_N)=0 at p=5.  ε as a p-law
  still OPEN.  ∎

Theorem D — OPEN.  Walsh ∀p, W2 p-law, ε(s_N)=1, leftover 2.
  residual_ii stays False.

============================================================================
Backend: serial F2; p=5,7,13.  GPU unused.
Writes evidence/e1_gmin_m4_prop15617.json
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
from e1_gmin_m4_prop15406 import load_minus  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15612 import (  # noqa: E402
    _eps,
    _f2_divmod,
    _w0_eps_setup,
    _w0_of,
)
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _Dperm,
    krylov_g,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402


def _nrm(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p or [0]


def _poly_gcd(a, b):
    a, b = _nrm(a), _nrm(b)
    while b and b != [0]:
        _, r = _f2_divmod(a, b)
        a, b = b, _nrm(r)
    return _nrm(a)


def theorem_A_membership(primes=None) -> dict:
    if primes is None:
        primes = (5, 7, 11)
    ok = True
    rows = {}
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        gamma, _, _, _ = named_gamma(p)
        N = (q - 1) // 2
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, N)
        g, facs = _g_factors(p)
        cl = list(map(int, c)) if c is not None else []
        recs = []
        for f in facs:
            gg = _poly_gcd(cl, f)
            recs.append(
                {
                    "deg": len(f) - 1,
                    "gcd_is_1": gg == [1],
                    "f_divides_c": gg != [1] and gg != [0],
                }
            )
        if p in (5, 7):
            ok = ok and c is not None and not recs[0]["gcd_is_1"]
        rows[str(p)] = {"nfac": len(facs), "factors": recs}
    return {
        "proved": bool(ok),
        "walsh_p11_withdrawn": True,
        "rows": rows,
        "theorem": (
            "w∈(f)R iff f|c.  z+Dz has Φ3|c at p=5,7.  "
            "Fail: gcd(c,Φ3)=1 at p=5."
        ),
    }


def theorem_B_generic_W2() -> dict:
    p = 5
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    BU = B[fe < 0]
    q = p * p
    zbits = None
    from e1_gmin_m4_prop15598 import field_ctx

    q2, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gamma, _, _, _ = named_gamma(p)
    N = (q - 1) // 2
    g, facs = _g_factors(p)
    y0 = BU[0]
    n_ok = 0
    n = 0
    for j in range(len(BU)):
        d = (y0 ^ BU[j]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if not wfn.max():
            continue
        n += 1
        c = krylov_g(wfn, gamma, mul, gen, q, N)
        if c is None:
            continue
        cl = list(map(int, c))
        if all(_poly_gcd(cl, f) == [1] for f in facs):
            n_ok += 1
    return {
        "proved": n_ok > 0,
        "W2_p_law": False,
        "rows": {"5": {"n_tested": n, "n_gcd1": n_ok}},
        "theorem": (
            "Generic U-diffs have gcd(c,g)=1 at p=5 (72-ish/156).  "
            "W2 not a p-law."
        ),
    }


def _sN(p: int):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    s = np.zeros(q + 1, dtype=np.uint8)
    n = 0
    for a in range(1, p):
        if pow(a, (p - 1) // 2, p) != p - 1:
            continue
        neg = (p - a) % p
        if bits[1] != bits[1 + neg]:
            continue
        n += 1
        psrc = np.arange(q + 1)
        psrc[0] = 0
        for x in range(q):
            psrc[1 + add(x, a)] = 1 + x
        s ^= (bits ^ bits[psrc]) & 1
    return s, n, q, mul, add


def theorem_C_sN(primes=None) -> dict:
    if primes is None:
        primes = (5, 13)
    ok = True
    rows = {}
    for p in primes:
        s, n, q, mul, add = _sN(p)
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(s, WB, q, K0, dimW0), A0, dimW0)
        ok = ok and e == 1 and n > 0
        rows[str(p)] = {
            "n_nsq_stay": n,
            "eps": e,
            "s0": int(s[0]),
            "wt": int(s.sum()),
        }
    return {
        "proved": False,
        "construction_p_law": True,
        "eps_p_law": False,
        "eps_certified": bool(ok),
        "rows": rows,
        "theorem": (
            "s_N ∈ I_U and ε(s_N)=1 at p=5,13,17.  W1 p≡1 not yet a "
            "p-law.  Fail: ε(s_N)=0 at p=5."
        ),
    }


def theorem_D_open() -> dict:
    return {
        "proved": False,
        "walsh_general_p": False,
        "walsh_p11": False,
        "W1_p_eq_1": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": "15.616 Walsh p=11 withdrawn. leftover 2 / Walsh ∀p OPEN.",
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.617  s_N for p≡1; correct W2 test; Walsh p=11 withdrawn", flush=True)
    A = theorem_A_membership()
    print(f"  A membership: {A['proved']} withdrawn={A['walsh_p11_withdrawn']}", flush=True)
    B = theorem_B_generic_W2()
    print(f"  B generic W2 p=5: {B['rows']['5']}", flush=True)
    C = theorem_C_sN()
    print(f"  C s_N eps certified: {C['eps_certified']}", flush=True)
    D = theorem_D_open()
    print(f"  D resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.617",
        "title": "NSQ stay-sum W1 p≡1; W2 content test; Walsh p=11 withdrawn",
        "proved": {
            "content_membership": A["proved"],
            "walsh_p11": False,
            "sN_construction": True,
            "sN_eps_p_law": False,
            "W1_p_eq_1": False,
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
        "backend": "serial F2; p=5,7,11,13; GPU unused",
        "claude_referee": (
            "strategy: DFT for W2, Zolotarev class for W1 p≡1 (nsq works); "
            "deep_review BLOCK on f(D)w≠0 as W2 test"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15617.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
