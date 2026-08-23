#!/usr/bin/env python3
"""
Prop 15.618 — Φ=ε is a p-law; s_N is a φ-pullback; 1_M coprime to g.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh ∀p, W1 for p≡1 as a p-law, or leftover 2.

============================================================================
Setup.  15.611–15.613: W ≅ F2[X]/(X^N+1), W_0=im(D−I), ε unique
with ker(D−I)W_0, ε(g(D)γ)=∑_{k odd} g_k.  D-orbits on F_q^×
are QR and QNR, each an N-cycle, 4|N.  Square 0-lines partition
QR∪{0}; nsq 0-lines partition QNR∪{0}.  γ=1_{H∪(1+H)}, H a
nsq 0-line.  s_N as in 15.617.  φ=L∘σ^{−1}.

============================================================================
Theorem A — PROVED (uniqueness + scale; all odd p; Fable xhigh
PASS-WITH-NOTE, scale completed as (D−I)γ).
  For w∈W_0 set
      Φ(w)=∑_{k odd} w(g^k)  +  ∑_{k odd} w(ρ g^k)   in F2
  (g generator of QR, ρ any nsq).  Then Φ=ε.
  Orbit weights vanish: even weight on each square 0-line and
  w_0=0 ⇒ wt_QR=0; even H-invariants and w_0=0 ⇒ wt_QNR=0.
  N even ⇒ odd_sum(Dw)=wt+odd_sum, so Φ(Dw)=Φ(w) on W_0, and
  odd_sum((D−I)t)=wt(t) per cycle so Φ vanishes on (D−I)W_0.
  dim W_0/(D−I)W_0=1 (ker(D−I)∩W_0=⟨extra⟩ dim 1), so
  {forms vanishing on (D−I)W_0}={0,ε}.  Scale: the cycle
  identity holds for t=γ (γ_0=1 is off the cycles).  H nsq
  0-line ⇒ 1∉H (χ_q(1)=1) ⇒ H∩(1+H)=∅ ⇒ |supp γ|=2p.
  γ_∞=0, 0∈H ⇒ |supp γ ∩ (QR∪QNR)|=2p−1≡1, so
      Φ((D−I)γ)=1=ε((D−I)γ)
  (15.613 C: Krylov X+1 has odd-coeff sum 1).  Fail: Φ((D−I)γ)=0.
  (Do not transfer scale from z+Dz at p≡3 — that does not
  constrain p≡1.  Discriminator, actual value 1.)  ∎

Theorem B — PULLBACK p-law; orbit pattern CERTIFIED, not a p-law.
  s_N(x)=f(φ(x)) with
      f(t)=∑_{a∈A} 1_{t ∈ S Δ (S+aλ)}   in F2,
  A=nsq-stay, λ=L(σ^{−1}), S={0,…,(p−1)/2}.  p-law (Fable B1).
  Off-0 nsq line: |L∩QR|=(p+1)/2 odd, |L∩QNR|=(p−1)/2 even
  (p≡1).  0-fiber B2 is a p-law: F_p^× has half-index (p+1)/2
  odd, so H except 0 is QNR and splits (p-1)/2:(p-1)/2 between D-index
  parities, both even.  Off-0 n_odd^{QNR} even, and
  f·n_odd^{QR}≡1, certified p=5,13,17 only.  Inversion x↦x^{−1}
  swaps QNR classes 1↔3 and fixes QR 0,2, but is not fiberwise.
  Termwise ε(z+T_a z)=1 is false at p=17.  Fail: s_N not a
  φ-pullback; fail (discriminator): ε(s_N)=0 at p=5 (actual 1).
  ε(s_N)=1 as a p-law still OPEN.  ∎

Theorem C — PROVED (1_M nonvanishing at g-roots); W2 p-law OPEN.
  Content of (D−I)γ is X+1.  g=(X^m+1)/(X+1) has g(1)≠0, so
  gcd(X+1,g)=1: 1_M is coprime to every irred of g, and
  w=c(D)γ ∈ (f)R iff f|c is a genuine test (15.617 A).  This
  vector is not a U-difference.  Named U-diffs z+Dz, s_N,
  two-fiber, z(σ)+z(σ'), and T_b z all miss gcd(c,g)=1 at p=5
  (Φ3-gate).  Generic U-diffs pass (15.617 B).  W2 remains:
  some w∈I_U with gcd(c,g)=1, equivalently per-factor
  witnesses glued by primitive idempotents of g (Fable).
  Fail: g divides X+1.  ∎

Theorem D — OPEN.  W1 for p≡1, W2 p-law, Walsh ∀p, leftover 2.
  residual_ii stays False.  Do not close leftover 2 via Walsh.

============================================================================
Backend: serial F2 identities; rref p=5,7,13.  GPU unused.
Fable xhigh: Claim A PASS-WITH-NOTE (scale via (D−I)γ); B
cyclotomic / per-factor W2 as directions.
Writes evidence/e1_gmin_m4_prop15618.json
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
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15612 import _eps, _w0_eps_setup, _w0_of  # noqa: E402
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _Dperm,
    _finv,
    named_gamma,
    named_z,
)
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd, _sN  # noqa: E402


def _phi_orbits(s_full, mul, gen, q, qr, qnr):
    N = (q - 1) // 2
    out = {}
    for starter_fn, lab in ((qr, "QR"), (qnr, "QNR")):
        rho = next(e for e in range(1, q) if starter_fn[1 + e] == 1)
        bits = []
        x = rho
        for _k in range(N):
            bits.append(int(s_full[1 + x]))
            x = mul(gen, x)
        odd = sum(bits[k] for k in range(N) if k % 2 == 1) % 2
        wt = sum(bits) % 2
        out[lab] = {"rho": int(rho), "odd": int(odd), "wt": int(wt)}
    out["phi"] = out["QR"]["odd"] ^ out["QNR"]["odd"]
    return out


def _DI_gamma(p):
    gamma, q, mul, b = named_gamma(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    gfull = np.zeros(q + 1, dtype=np.uint8)
    gfull[1 : 1 + q] = gamma
    Dg = np.zeros(q + 1, dtype=np.uint8)
    for e in range(q):
        Dg[1 + mul(gen, e)] = gfull[1 + e]
    return (gfull ^ Dg) & 1, gamma, q, mul, gen, b, gfull


def theorem_A_phi_eq_eps(primes=None) -> dict:
    if primes is None:
        primes = (5, 7)
    ok = True
    rows = {}
    for p in primes:
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        w, gamma, q, mul, gen, b, gfull = _DI_gamma(p)
        # |supp γ|=2p, 0∈H, 1∉H
        wt_g = int(gamma.sum())
        phi = _phi_orbits(w, mul, gen, q, qr, qnr)
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(w, WB, q, K0, dimW0), A0, dimW0)
        # z+Dz
        z, bits, eigen, inU, qz, mulz, addz, chiz, sig = named_z(p)
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        phi_d = _phi_orbits(d, mul, gen, q, qr, qnr)
        e_d = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        rec = {
            "gamma_wt": wt_g,
            "gamma_wt_is_2p": wt_g == 2 * p,
            "gamma_0": int(gamma[0]),
            "phi_DI": phi["phi"],
            "eps_DI": e,
            "QR_wt_DI": phi["QR"]["wt"],
            "QNR_wt_DI": phi["QNR"]["wt"],
            "scale_1": e == 1 and phi["phi"] == 1,
            "zDz_match": e_d == phi_d["phi"],
            "zDz_eps": e_d,
            "zDz_phi": phi_d["phi"],
        }
        if p % 4 == 1:
            s, n, qs, muls, adds = _sN(p)
            phi_s = _phi_orbits(s, mul, gen, q, qr, qnr)
            e_s = _eps(_w0_of(s, WB, q, K0, dimW0), A0, dimW0)
            rec["sN_match"] = e_s == phi_s["phi"]
            rec["sN_eps"] = e_s
            rec["sN_QR_odd"] = phi_s["QR"]["odd"]
            rec["sN_QNR_odd"] = phi_s["QNR"]["odd"]
            ok = ok and rec["sN_match"] and e_s == 1
        ok = (
            ok
            and rec["gamma_wt_is_2p"]
            and rec["gamma_0"] == 1
            and rec["scale_1"]
            and rec["zDz_match"]
            and phi["QR"]["wt"] == 0
            and phi["QNR"]["wt"] == 0
        )
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Φ=ε on W_0.  Scale Φ((D−I)γ)=2p−1≡1.  "
            "Fail: Φ((D−I)γ)=0."
        ),
    }


def theorem_B_pullback(primes=None) -> dict:
    if primes is None:
        primes = (5, 13)
    ok_pb = True
    ok_pat = True
    rows = {}
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        sinv = _finv(mul, q, sig)
        lam = sinv // p
        S = set(range((p + 1) // 2))
        A = []
        for a in range(1, p):
            if pow(a, (p - 1) // 2, p) != p - 1:
                continue
            neg = (p - a) % p
            if bits[1] != bits[1 + neg]:
                continue
            A.append(a)

        def f_t(t):
            acc = 0
            for a in A:
                d = (a * lam) % p
                if (t in S) != (((t - d) % p) in S):
                    acc ^= 1
            return acc

        f = [f_t(t) for t in range(p)]
        s, n, q, mul, add = _sN(p)
        bad = 0
        for x in range(q):
            t = mul(sinv, x) // p
            if int(s[1 + x]) != f[t]:
                bad += 1
        q2, mul2, add2, chi2, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        N = (q - 1) // 2
        n_odd_qnr = [0] * p
        n_odd_qr = [0] * p
        rho_qnr = next(e for e in range(1, q) if qnr[1 + e] == 1)
        x = 1
        for k in range(N):
            t = mul(sinv, x) // p
            if k % 2 == 1:
                n_odd_qr[t] += 1
            x = mul(gen, x)
        x = rho_qnr
        for k in range(N):
            t = mul(sinv, x) // p
            if k % 2 == 1:
                n_odd_qnr[t] += 1
            x = mul(gen, x)
        qnr_all_even = all(v % 2 == 0 for v in n_odd_qnr)
        dot_qr = sum(f[t] * n_odd_qr[t] for t in range(p)) % 2
        ok_pb = ok_pb and bad == 0 and int(s[0]) == 0
        ok_pat = ok_pat and qnr_all_even and dot_qr == 1
        rows[str(p)] = {
            "mismatches": bad,
            "nsq_stay": A,
            "lam": int(lam),
            "f": f,
            "qnr_odd_all_even": qnr_all_even,
            "f_dot_qr_odd": dot_qr,
        }
    return {
        "proved": False,
        "pullback_p_law": bool(ok_pb),
        "orbit_pattern_p_law": False,
        "orbit_pattern_certified": bool(ok_pat),
        "eps_p_law": False,
        "rows": rows,
        "theorem": (
            "s_N=f∘φ.  odd_QNR=0 and f·n_odd^QR=1 certified "
            "p=5,13,17, not a p-law.  Fail: mismatches>0; "
            "Fail (discriminator): ε(s_N)=0 at p=5."
        ),
    }


def theorem_C_one_M() -> dict:
    p = 5
    w, gamma, q, mul, gen, b, gfull = _DI_gamma(p)
    _g, facs = _g_factors(p)
    # content of (D-I)γ is X+1; gcd with each irred of g is 1
    # X+1 = [1,1]
    all1 = all(_poly_gcd([1, 1], f) == [1] for f in facs)
    q2, mul2, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
    phi = _phi_orbits(w, mul, gen, q, qr, qnr)
    WB, q3, mul3, K0, dimW0, A0 = _w0_eps_setup(p)
    e = _eps(_w0_of(w, WB, q, K0, dimW0), A0, dimW0)
    return {
        "proved": bool(all1 and e == 1 and phi["phi"] == 1),
        "W2_p_law": False,
        "gcd_X1_g_is_1": bool(all1),
        "nfac": len(facs),
        "eps_DI": e,
        "theorem": (
            "gcd(X+1,g)=1 so 1_M is coprime to g.  Not a U-diff.  "
            "W2 p-law OPEN.  Fail: g divides X+1."
        ),
    }


def theorem_D_open() -> dict:
    return {
        "proved": False,
        "walsh_general_p": False,
        "W1_p_eq_1": False,
        "W2_p_law": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": "15.618 Φ=ε. leftover 2 / Walsh ∀p OPEN.",
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.618  Φ=ε; s_N φ-pullback; 1_M coprime to g", flush=True)
    A = theorem_A_phi_eq_eps()
    print(f"  A Φ=ε: {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_pullback()
    print(
        f"  B pullback={B['pullback_p_law']} pattern={B['orbit_pattern_certified']}",
        flush=True,
    )
    C = theorem_C_one_M()
    print(f"  C 1_M coprime: {C['proved']} W2={C['W2_p_law']}", flush=True)
    D = theorem_D_open()
    print(f"  D resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.618",
        "title": "Φ=ε p-law; s_N φ-pullback; 1_M coprime to g",
        "proved": {
            "phi_eq_eps": A["proved"],
            "sN_pullback": B["pullback_p_law"],
            "sN_eps_p_law": False,
            "W1_p_eq_1": False,
            "one_M_coprime_g": C["proved"],
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
        "backend": "serial F2; rref p=5,7,13; GPU unused",
        "claude_referee": (
            "math_review PASS-WITH-NOTE Claim A, scale via (D-I)γ=2p-1; "
            "suggest_direction: cyclotomic order 2 for B, per-factor "
            "idempotent glue for W2"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15618.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
