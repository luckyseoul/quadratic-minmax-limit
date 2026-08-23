#!/usr/bin/env python3
"""
Prop 15.613 — W1 reduces to one ε-bit per p mod 4, via a named
Max- point z in U.  The bit is not yet a p-law.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E or W1 as a p-law.

============================================================================
Setup.  15.612: Walsh ⇔ W1 ∧ W2.  W1 is I_U ⊄ (D−I)W_0, i.e. some
U-difference has ε=1, where ε: W_0→F2 is the unique D-invariant
functional with kernel (D−I)W_0.  Paley halfspace h (ρ=1):
h_∞=+1, h_{a+bω}=+1 iff b∈S={0,…,(p−1)/2}.  15.254: σ nsq,
z=D_monomial(h∘π_σ^{−1}) lies in Max−; z_∞=−1, z_x=h(σ^{−1}x).

============================================================================
Theorem A — PROVED (15.254 + ρ=1; all odd p).
  C_∞0=1 and h(0)=+1, so C_∞0 z_∞ z_0=−1.  Thus z∈U.  Fail:
  z_0=−1.  ∎

Theorem B — PROVED (affine; all odd p).
  U ⊂ y_*+V with V=⟨1⟩⊕W_0.  D preserves U and annihilates ⟨1⟩,
  and (D−I)W_0=ker ε, so ε(y+Dy)=ε(y_*+D y_*) is CONSTANT on U.
  Fail: ε(y+Dy) takes both values on U at p=3.  ∎

Theorem C — PROVED (15.611 F2[M]; all odd p).
  Let H be a nsq 0-line, γ=1_{H ∪ (1+H)} (even H-invariants on
  the 0- and 1-cosets).  This is 1_M under W ≅ F2[M], hence a
  cyclic generator of W.  For w∈W_0 write w=g(D)γ with g(1)=0;
  then w=(X+1)f γ and ε(w)=f(1).  Char 2: f(1)=g'(1)=∑_{k odd} g_k.
  Fail: D-orbit of γ has rank <N; fail: odd-sum ≠ ε.  ∎

Theorem D — W1-3 / W1-1 CERTIFIED, values not a p-law
  (Fable xhigh BLOCK on the ε-bit).
  p≡3 (mod 4): the U-constant ε(y+Dy) equals 1 at p=3,7,11
  (any nsq σ).  p≡1 (mod 4): let α=(p+1)/2∈F_p and σ with
  L(σ^{−1})=p−2 (exists: the off-0 square-direction line
  {c+(p−2)ω} has |QNR|=(p+1)/2 by 15.604).  Then T_α z∈U
  because −L(σ^{−1}α)=1∈S, and ε(z+T_α z)=1 at p=5,13.
  Stay/existence of this σ are p-laws; the value of ε is not.
  Fail: claim ε(y+Dy)=1 at p=5.  ∎

Theorem E — OPEN.  W1 as a p-law, W2, Walsh.  residual_ii False.

============================================================================
Backend: identities serial; rref p=3,5,7.  GPU unused.
Fable xhigh: construction/stay/constancy stand; ε-value BLOCK.
Writes evidence/e1_gmin_m4_prop15613.json
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

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def _finv(mul, q, u):
    r, base = 1, u
    e = q - 2
    while e:
        if e & 1:
            r = mul(r, base)
        base = mul(base, base)
        e >>= 1
    return r


def named_z(p: int, sig=None):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    h = np.sign(halfspace_boolean_vector(p)).astype(np.int8)
    if sig is None:
        sig = next(e for e in range(1, q) if chi(e) == -1)
    sinv = _finv(mul, q, sig)
    z = np.zeros(q + 1, dtype=np.int8)
    z[0] = np.int8(-h[0])
    for x in range(q):
        z[1 + x] = h[1 + mul(sinv, x)]
    C = paley_conference_prime_power(p)
    eigen = bool(
        np.allclose(C @ z.astype(np.float64), -p * z.astype(np.float64))
    )
    bits = ((1 - z) // 2).astype(np.uint8)
    inU = int(np.rint(C[0, 1])) * int(z[0]) * int(z[1]) == -1
    return z, bits, eigen, inU, q, mul, add, chi, sig


def named_gamma(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    b = next(e for e in range(1, q) if chi(e) == -1)
    fn = np.zeros(q, dtype=np.uint8)
    for t in range(p):
        fn[add(1, mul(t, b))] = 1
        fn[mul(t, b)] ^= 1
    return fn, q, mul, b


def _Dperm(mul, g, q):
    perm = np.arange(q + 1)
    perm[0] = 0
    for e in range(q):
        perm[1 + mul(g, e)] = 1 + e
    return perm


def krylov_g(wfn, gamma, mul, g, q, N):
    cols = [gamma.copy()]
    cur = gamma.copy()
    for _ in range(N - 1):
        cur = _dil_fn(cur, mul, g, q)
        cols.append(cur.copy())
    M = np.stack(cols, axis=1)
    Aug = np.concatenate([M, wfn.reshape(-1, 1)], axis=1)
    R, pivots, rank = gf2_rref(Aug.copy())
    if rank > N:
        return None
    c = np.zeros(N, dtype=np.uint8)
    for i, pv in enumerate(pivots):
        if pv < N:
            c[pv] = R[i, N]
    recon = (M.astype(np.int32) @ c.astype(np.int32) % 2).astype(np.uint8)
    if not np.array_equal(recon, wfn):
        return None
    return c


def theorem_A_z_in_U(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        ok = ok and eigen and inU and int(z[0]) == -1 and int(z[1]) == 1
        rows[str(p)] = {
            "eigen_minus": eigen,
            "in_U": inU,
            "z_inf": int(z[0]),
            "z_0": int(z[1]),
            "sig": int(sig),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Named halfspace-anti z lies in Max- and in U.  Fail: z_0=−1.",
    }


def theorem_B_constant_on_U(primes=None) -> dict:
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
        omega = _primitive(mul, q)
        g = mul(omega, omega)
        Dperm = _Dperm(mul, g, q)

        def epsb(d):
            return _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)

        take = min(40, len(BU))
        vals = [epsb((BU[i] ^ BU[i][Dperm]) & 1) for i in range(take)]
        constant = len(set(vals)) == 1
        expected = 1 if p % 4 == 3 else 0
        ok = ok and constant and vals[0] == expected
        rows[str(p)] = {
            "take": take,
            "constant": constant,
            "value": vals[0],
            "expected_mod4": expected,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "ε(y+Dy) is constant on U (affine).  Fail: both values at p=3."
        ),
    }


def theorem_C_odd_coeffs(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        gamma, qg, mulg, b = named_gamma(p)
        N = (q - 1) // 2
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, N)
        odd = int(c[1::2].sum() % 2) if c is not None else None
        even = int(c[0::2].sum() % 2) if c is not None else None
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        orbit = []
        cur = gamma.copy()
        for _ in range(N):
            orbit.append(cur.copy())
            cur = _dil_fn(cur, mul, gen, q)
        rank = int(gf2_rref(np.stack(orbit, axis=1).copy())[2])
        ok = ok and c is not None and odd == e and rank == N
        rows[str(p)] = {
            "gamma_cyclic": rank == N,
            "odd_sum": odd,
            "even_sum": even,
            "eps": e,
            "odd_equals_eps": odd == e,
            "even_equals_eps": even == e,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "ε(w)=∑_{k odd} g_k for w=g(D)γ, γ=1_{H∪(1+H)}.  "
            "Fail: D-orbit of γ has rank <N; fail: odd-sum ≠ ε."
        ),
    }


def theorem_D_certified(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    B = theorem_B_constant_on_U(primes)
    ok_pat = B["proved"]
    rows = {}
    for p in primes:
        rec = dict(B["rows"][str(p)])
        rec["W1_from_D"] = rec["value"] == 1
        rec["p_mod_4"] = p % 4
        # T_α at p≡1
        if p % 4 == 1:
            z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
            # pin σ: L(σ^{-1})=p-2
            q2, mul2, add2, chi2, frob, norm, ia, ib = field_ctx(p)
            target = p - 2
            sig_pin = None
            for s in range(1, q):
                if chi(s) != -1:
                    continue
                sinv = _finv(mul, q, s)
                if (sinv // p) == target:
                    sig_pin = s
                    break
            z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p, sig=sig_pin)
            a = (p + 1) // 2
            c0, c1 = a % p, a // p
            neg = ((p - c0) % p) + ((p - c1) % p) * p
            stay = bool(bits[1] == bits[1 + neg])
            WB, q2, mul3, K0, dimW0, A0 = _w0_eps_setup(p)
            eT = None
            if stay:
                psrc = np.arange(q + 1)
                psrc[0] = 0
                for x in range(q):
                    psrc[1 + add(x, a)] = 1 + x
                eT = _eps(
                    _w0_of((bits ^ bits[psrc]) & 1, WB, q, K0, dimW0),
                    A0,
                    dimW0,
                )
            rec["sig_pin"] = sig_pin
            rec["T_stay"] = stay
            rec["T_eps"] = eT
            rec["W1_from_T"] = eT == 1
            ok_pat = ok_pat and stay and eT == 1
        rows[str(p)] = rec
    return {
        "proved": False,
        "W1_p_law": False,
        "pattern_certified": bool(ok_pat),
        "rows": rows,
        "theorem": (
            "ε(y+Dy)=1 iff p≡3 (mod 4) certified p=3,5,7; T_α at p=5 "
            "gives ε=1.  Not a p-law.  Fail: ε(y+Dy)=1 at p=5."
        ),
    }


def theorem_E_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "W1_p_law": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "W1 is one ε-bit per p mod 4 on a named U-point.  "
            "Value of the bit stays open.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.613  W1 named z; ε-bit not a p-law", flush=True)
    A = theorem_A_z_in_U()
    print(f"  A z in U: {A['proved']}", flush=True)
    B = theorem_B_constant_on_U()
    print(f"  B constant on U: {B['proved']}", flush=True)
    C = theorem_C_odd_coeffs()
    print(f"  C odd-coeff formula: {C['proved']}", flush=True)
    D = theorem_D_certified()
    print(f"  D pattern certified: {D['pattern_certified']}", flush=True)
    E = theorem_E_open()
    print(f"  E Walsh open: resii={E['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.613",
        "title": "W1 named Max- in U; ε-bit per p mod 4 not a p-law",
        "proved": {
            "named_z_in_U": A["proved"],
            "eps_Dy_constant_on_U": B["proved"],
            "eps_odd_krylov": C["proved"],
            "W1_p_law": False,
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
        "backend": "serial F2; rref p=3,5,7; GPU unused",
        "claude_referee": (
            "deep_review BLOCK on ε-value; construction/stay/constancy kept"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15613.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
