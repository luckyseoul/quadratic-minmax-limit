#!/usr/bin/env python3
"""
Prop 15.619 — odd_QNR(s_N)=0 is a p-law (biquadratic fiber counts).

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p≡1, W2, Walsh ∀p, or leftover 2.

============================================================================
Setup.  15.618: Φ=ε; s_N=f∘φ.  odd_QNR(s_N)=∑_t f(t) n_odd^{QNR}(t).
p≡1 (mod 4), p=a²+b² with a odd >0, b even >0 (unique).  ψ=χ_4∘N
on F_q^*, order 4, ψ²=χ_q.  The two primitive values of ψ are the
two QNR D-index parities.

============================================================================
Theorem A — PROVED (Gauss + Z[i] UFD; all p≡1 (mod 4)).
  0-fiber (15.618): n_odd^{QNR}(0)=(p−1)/2 even.
  Off-0 nsq affine line L: ∑_L χ_q=1 (15.598).  For primitive
  ζ∈{±i},
      n_ζ=|{x∈L: ψ(x)=ζ}|
         =(p−1)/4 + (1/2) Re(ζ̄ S),   S=∑_L ψ ∈ Z[i].
  Completing the square, S is a 4th-root times
      J=∑_{s∈F_p} χ_4(s²−Δ),   Δ nsq in F_p.
  Gauss: |J|=√p.  J∈Z[i] of norm p, so J is a unit times a+bi.
  Hence Re(ζ̄ S)∈{±a,±b}.  n_ζ∈Z forces Re even, so ±b
  (a odd).  Thus n_ζ=(p−1)/4 ± b/2.  Write b=2c, a=2m+1:
      n_ζ=m(m+1)+c(c±1)
  even.  Both QNR parities are primitive ζ, so n_odd^{QNR}(t)
  is even on every fiber.  Therefore odd_QNR(w)=0 for every
  φ-pullback w, in particular s_N.  Fail (discriminator):
  n_odd^{QNR} odd at p=5 (actual values even); fail: |J|²≠p.  ∎

Theorem B — OPEN.  f·n_odd^{QR}≡1 (B3).  Certified p=5,13,17
  only.  Needed for ε(s_N)=1.  Termwise ε(z+T_a z)=1 is false
  at p=17.  Fail: claim that identity as a p-law.  ∎

Theorem C — OPEN.  W1 p≡1, W2, Walsh ∀p, leftover 2.
  residual_ii stays False.

============================================================================
Backend: Gauss identity serial; fiber census p=5,13.  GPU unused.
OpenAI referee: quota (retry after 2026-08-27); not routed to Fable.
Writes evidence/e1_gmin_m4_prop15619.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_z  # noqa: E402
from e1_gmin_m4_prop15618 import theorem_B_pullback  # noqa: E402


def _ab(p: int):
    for a in range(1, p, 2):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = int(round(b2**0.5))
        if b * b == b2 and b % 2 == 0 and b > 0:
            return a, b
    return None, None


def _chi4(p: int):
    imag = next(x for x in range(p) if (x * x) % p == p - 1)
    roots = {1: 1 + 0j, p - 1: -1 + 0j, imag: 1j, (p - imag) % p: -1j}

    def chi4(x):
        if x % p == 0:
            return 0j
        y = pow(int(x) % p, (p - 1) // 4, p)
        return roots[y]

    return chi4


def theorem_A_qnr_even(primes=None) -> dict:
    if primes is None:
        primes = (5, 13)
    ok = True
    rows = {}
    for p in primes:
        a, b = _ab(p)
        chi4 = _chi4(p)
        # |J|^2 = p for one nsq Δ
        nsq = next(d for d in range(1, p) if pow(d, (p - 1) // 2, p) == p - 1)
        J = sum(chi4((s * s - nsq) % p) for s in range(p))
        abs2 = int(round(abs(J) ** 2))
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        sinv = _finv(mul, q, sig)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
        N = (q - 1) // 2
        rho = next(e for e in range(1, q) if qnr[1 + e] == 1)
        n_odd = [0] * p
        x = rho
        for k in range(N):
            t = mul(sinv, x) // p
            if k % 2 == 1:
                n_odd[t] += 1
            x = mul(gen, x)
        pred = sorted({(p - 1) // 4 + b // 2, (p - 1) // 4 - b // 2})
        off = sorted(set(n_odd[1:]))
        all_even = all(v % 2 == 0 for v in n_odd)
        ok = (
            ok
            and a is not None
            and a % 2 == 1
            and b % 2 == 0
            and abs2 == p
            and n_odd[0] == (p - 1) // 2
            and off == pred
            and all_even
        )
        rows[str(p)] = {
            "a": a,
            "b": b,
            "J_abs2": abs2,
            "n_odd_0": n_odd[0],
            "n_odd_off": off,
            "pred": pred,
            "all_even": all_even,
        }
    return {
        "proved": bool(ok),
        "odd_QNR_sN_zero": True,
        "rows": rows,
        "theorem": (
            "n_odd^{QNR}(t)=(p-1)/4 ± b/2 even.  odd_QNR(s_N)=0.  "
            "Fail: n_odd odd at p=5; fail: |J|^2≠p."
        ),
    }


def theorem_B_qr_dot_open() -> dict:
    B = theorem_B_pullback()
    return {
        "proved": False,
        "f_dot_qr_certified": B["orbit_pattern_certified"],
        "f_dot_qr_p_law": False,
        "rows": {k: {"f_dot_qr_odd": rec["f_dot_qr_odd"]} for k, rec in B["rows"].items()},
        "theorem": (
            "f·n_odd^{QR}≡1 certified p=5,13, not a p-law.  "
            "Fail: claim it for every p≡1."
        ),
    }


def theorem_C_open() -> dict:
    return {
        "proved": False,
        "W1_p_eq_1": False,
        "W2_p_law": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": "odd_QNR(s_N)=0.  ε(s_N)=1 still needs B3. leftover 2 OPEN.",
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.619  odd_QNR(s_N)=0 via a^2+b^2 fiber counts", flush=True)
    A = theorem_A_qnr_even()
    print(f"  A QNR even: {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_qr_dot_open()
    print(f"  B QR dot certified={B['f_dot_qr_certified']}", flush=True)
    C = theorem_C_open()
    print(f"  C resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.619",
        "title": "odd_QNR(s_N)=0 p-law; B3 still open",
        "proved": {
            "odd_QNR_zero": A["proved"],
            "sN_eps_p_law": False,
            "W1_p_eq_1": False,
            "W2_p_law": False,
            "walsh_general_p": False,
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
        "walsh_15_406_E": "OPEN",
        "backend": "serial Gauss; fiber census p=5,13; GPU unused",
        "openai_referee": "quota until 2026-08-27; not routed to Fable",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15619.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
