#!/usr/bin/env python3
"""
Prop 15.604 — QR / QNR indicators in H0=ker S, and the dilation
fixed space ker(D−I)∩H0.  All odd p, Max-free.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** prove H0/⟨1⟩ irreducible, Walsh 15.406 E, or that
H0→F2^{QR} is surjective.

============================================================================
Setup.  S is square-direction {∞}∪L incidence (15.598–15.600).
H0=ker S, dim n/2, radical ⟨1⟩.  QR = (F_q^×)², QNR = F_q^×\\QR,
each of size N=(q−1)/2.  Column 0 is ∞; column 1+e is field e.
D is square dilation z↦g z on P¹, g a generator of M=(F_q^×)²,
order N; D fixes ∞ and 0.  15.598 B: χ_q|_{F_p^×}≡1 and, for
j∉{∞}∪L, Σ_{i∈L} C_{ij}=−χ(b).

============================================================================
Theorem A — PROVED (line counts; 15.598 A,B; all odd p).
  1_QR ∈ H0  ⇔  p≡1 (mod 4);
  1_QNR ∈ H0  ⇔  p≡3 (mod 4).
  They are never both in H0 (their F2-sum is 1_{F_q^×}=1+e_0+e_∞,
  and e_0+e_∞ ∉ H0: a square line with 0∉L meets {0,∞} only at ∞).

  Square line L=a+F_p b, χ(b)=1.  ⟨1_QR, 1_{∞}∪L⟩ = |L∩QR| (mod 2),
  likewise QNR; ∞ is in neither.

  (i) 0∈L.  Then L=F_p b.  For t∈F_p^×, χ_q(tb)=χ_q(t)χ_q(b)=1
      (15.598 B: χ_q|_{F_p^×}≡1).  So L\\{0}⊂QR,
      |L∩QR|=p−1 even, |L∩QNR|=0 even.

  (ii) 0∉L.  15.598 B with j=0: Σ_{z∈L} χ_q(z)=−χ(b)=−1.
      |L∩QR|−|L∩QNR|=−1 and |L∩QR|+|L∩QNR|=p, so
      |L∩QR|=(p−1)/2, |L∩QNR|=(p+1)/2.
      (p−1)/2 even ⇔ p≡1 (mod 4);
      (p+1)/2 even ⇔ p≡3 (mod 4).

  Fail: swap the congruences; fail: 1_QR ∈ H0 at p≡3;
  fail: |L∩QR|=(p+1)/2 on a 0-line.  ∎

Theorem B — PROVED (orbits of D; all odd p).
  D^N=I as a permutation of P¹.  Ambient ker(D−I) is 4-dimensional,
  spanned by {e_∞, e_0, 1_QR, 1_QNR} (orbits {∞}, {0}, QR, QNR).
  ker(D−I)∩H0 has dimension 2:
      p≡1 (mod 4)  ⇒  ⟨1, 1_QR⟩,
      p≡3 (mod 4)  ⇒  ⟨1, 1_QNR⟩.
  Proof.  Write x=α e_∞+β e_0+γ 1_QR+δ 1_QNR.  Square 0-line:
  ⟨x,row⟩=α+β (p−1 even).  Need α=β.  Square line 0∉L:
  ⟨x,row⟩=α+γ(p−1)/2+δ(p+1)/2.
  p≡1 ⇒ (p−1)/2 even, (p+1)/2 odd ⇒ α=δ, so
  x=α(e_∞+e_0+1_QNR)+γ 1_QR=α·1+(α+γ)1_QR.
  p≡3 ⇒ α=γ and x ∈ ⟨1, 1_QNR⟩.
  Fail: dim=1 (only ⟨1⟩); fail: dim=4 (whole ambient);
  fail: e_∞ ∈ H0.  ∎

Theorem C — CENSUS, not a p-law.
  Restriction H0→F2^{QR} is not surjective in general (ker dims
  2,3,6,12,7 at p=3,5,7,11,13).  Do not quote a cycle-minpoly
  X^N+1 via QR coordinates from a false onto map.

Theorem D — OPEN.  H0/⟨1⟩ irreducible as G_aff^□-module, and
  Walsh 15.406 E, stay OPEN.  residual_ii stays False.

============================================================================
Backend: identities serial (inherently); line-count / ker checks
p=3,5,7,11,13.  GPU unused.
Writes evidence/e1_gmin_m4_prop15604.json
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
from walsh_linecode_rank import square_line_matrix  # noqa: E402


def _in_ker(S: np.ndarray, v: np.ndarray) -> bool:
    return int((S.astype(np.int32) @ v.astype(np.int32) % 2).max()) == 0


def _qr_qnr(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    qr = np.zeros(q + 1, dtype=np.uint8)
    qnr = np.zeros(q + 1, dtype=np.uint8)
    n_qr = n_qnr = 0
    for e in range(q):
        c = chi(e)
        if c == 1:
            qr[1 + e] = 1
            n_qr += 1
        elif c == -1:
            qnr[1 + e] = 1
            n_qnr += 1
    return q, mul, add, chi, qr, qnr, n_qr, n_qnr


def _line_counts(S: np.ndarray, chi, q: int) -> tuple[list[int], list[int], list[bool]]:
    qr_counts = []
    qnr_counts = []
    through0 = []
    for v in S:
        nqr = nqnr = 0
        for e in range(q):
            if v[1 + e] == 0:
                continue
            c = chi(e)
            if c == 1:
                nqr += 1
            elif c == -1:
                nqnr += 1
        qr_counts.append(nqr)
        qnr_counts.append(nqnr)
        through0.append(bool(v[1] == 1))
    return qr_counts, qnr_counts, through0


def _primitive(mul, q: int) -> int:
    n = q - 1
    # prime factors of n (n=(p-1)(p+1), tiny)
    fac = []
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        fac.append(m)

    def fpow(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    for e in range(2, q):
        if all(fpow(e, n // r) != 1 for r in fac):
            return e
    raise RuntimeError("no primitive root")


def _dilation_perm(p: int, mul, q: int) -> np.ndarray:
    omega = _primitive(mul, q)
    g = mul(omega, omega)  # generator of squares
    pi = np.zeros(q + 1, dtype=np.int64)
    pi[0] = 0
    pi[1] = 1
    for e in range(1, q):
        pi[1 + e] = 1 + mul(g, e)
    return pi, g


def theorem_A_qr_in_h0(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        S = square_line_matrix(p)
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        chi_fp = all(chi(t) == 1 for t in range(1, p))
        qr_in = _in_ker(S, qr)
        qnr_in = _in_ker(S, qnr)
        want_qr = p % 4 == 1
        want_qnr = p % 4 == 3
        qr_c, qnr_c, th0 = _line_counts(S, chi, q)
        thru = [i for i, t in enumerate(th0) if t]
        non = [i for i, t in enumerate(th0) if not t]
        thru_qr = {qr_c[i] for i in thru}
        thru_qnr = {qnr_c[i] for i in thru}
        non_qr = {qr_c[i] for i in non}
        non_qnr = {qnr_c[i] for i in non}
        e0inf = np.zeros(q + 1, dtype=np.uint8)
        e0inf[0] = 1
        e0inf[1] = 1
        e0inf_in = _in_ker(S, e0inf)
        einf = np.zeros(q + 1, dtype=np.uint8)
        einf[0] = 1
        einf_in = _in_ker(S, einf)
        ok = ok and chi_fp
        ok = ok and n_qr == n_qnr == (q - 1) // 2
        ok = ok and qr_in is want_qr
        ok = ok and qnr_in is want_qnr
        ok = ok and not (qr_in and qnr_in)
        ok = ok and (qr_in or qnr_in)
        ok = ok and thru_qr == {p - 1}
        ok = ok and thru_qnr == {0}
        ok = ok and non_qr == {(p - 1) // 2}
        ok = ok and non_qnr == {(p + 1) // 2}
        ok = ok and e0inf_in is False
        ok = ok and einf_in is False
        # fail-when-wrong: a 0-line is not all-QNR, and an off-0
        # line is not empty of QR.  (p=3 has p−1=(p+1)/2 so that
        # pair of numbers is not a discriminator.)
        ok = ok and thru_qr != {0}
        ok = ok and non_qr != {0}
        if p > 3:
            ok = ok and thru_qr != {(p + 1) // 2}
        rows[str(p)] = {
            "p_mod_4": p % 4,
            "qr_in_H0": qr_in,
            "qnr_in_H0": qnr_in,
            "chi_Fp_one": chi_fp,
            "thru0_L_cap_QR": p - 1,
            "off0_L_cap_QR": (p - 1) // 2,
            "e0_einf_in_H0": e0inf_in,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "1_QR ∈ H0 iff p≡1 (mod 4); 1_QNR ∈ H0 iff p≡3 (mod 4).  "
            "Fail: swap congruences; fail: 1_QR in H0 at p≡3."
        ),
    }


def theorem_B_dilation_fixed(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7, 11, 13)
    A = theorem_A_qr_in_h0(primes)
    ok = A["proved"]
    rows = {}
    for p in primes:
        S = square_line_matrix(p)
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        N = (q - 1) // 2
        pi, g = _dilation_perm(p, mul, q)
        # D^N = I
        x = np.arange(q + 1, dtype=np.int64)
        for _ in range(N):
            x = pi[x]
        dn_id = bool(np.array_equal(x, np.arange(q + 1)))
        # D ≠ I (N>1 for p≥3)
        d_id = bool(np.array_equal(pi, np.arange(q + 1)))
        einf = np.zeros(q + 1, dtype=np.uint8)
        einf[0] = 1
        e0 = np.zeros(q + 1, dtype=np.uint8)
        e0[1] = 1
        one = np.ones(q + 1, dtype=np.uint8)
        basis = [einf, e0, qr, qnr]
        inker = [v for v in basis if _in_ker(S, v)]
        # F2-span of {e_∞,e_0,1_QR,1_QNR} ∩ H0
        span_in = []
        for mask in range(16):
            v = np.zeros(q + 1, dtype=np.uint8)
            for i in range(4):
                if mask & (1 << i):
                    v ^= basis[i]
            if _in_ker(S, v):
                span_in.append(v)
        M = np.stack(span_in, axis=1) if span_in else np.zeros((q + 1, 0), dtype=np.uint8)
        dim = int(gf2_rref(M)[2]) if M.size else 0
        extra = qr if p % 4 == 1 else qnr
        extra_in = _in_ker(S, extra)
        one_in = _in_ker(S, one)
        # predicted 2-space: 0, extra, one, one+extra
        pred = [np.zeros(q + 1, dtype=np.uint8), extra, one, one ^ extra]
        pred_in = all(_in_ker(S, v) for v in pred)
        ok = ok and dn_id and (not d_id) and dim == 2 and extra_in and one_in
        ok = ok and pred_in
        ok = ok and dim != 1 and dim != 4
        ok = ok and not _in_ker(S, einf)
        rows[str(p)] = {
            "N": N,
            "D_N_is_id": dn_id,
            "dim_ker_DminusI_cap_H0": dim,
            "n_ambient_in_H0": len(inker),
            "extra_in_H0": extra_in,
            "span_card": len(span_in),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "ker(D−I)∩H0 = ⟨1, 1_QR⟩ (p≡1) or ⟨1, 1_QNR⟩ (p≡3), dim 2.  "
            "D^N=I.  Fail: dim 1 or 4; fail: e_∞ in H0."
        ),
    }


def theorem_C_restriction_census(primes=None) -> dict:
    """Not a p-law.  Restriction H0→F2^{QR} is not onto in general."""
    if primes is None:
        primes = (3, 5, 7, 11)
    rows = {}
    onto = True
    for p in primes:
        S = square_line_matrix(p)
        q, mul, add, chi, qr, qnr, n_qr, n_qnr = _qr_qnr(p)
        H0, _ = gf2_nullspace(S)
        qr_idx = [1 + e for e in range(q) if chi(e) == 1]
        qnr_idx = [1 + e for e in range(q) if chi(e) == -1]
        R = H0[qr_idx, :]
        Rn = H0[qnr_idx, :]
        r = int(gf2_rref(R)[2])
        rn = int(gf2_rref(Rn)[2])
        ker = int(H0.shape[1] - r)
        kern = int(H0.shape[1] - rn)
        onto = onto and (r == n_qr)
        rows[str(p)] = {
            "dim_H0": int(H0.shape[1]),
            "n_QR": n_qr,
            "rank_to_QR": r,
            "ker_to_QR": ker,
            "rank_to_QNR": rn,
            "ker_to_QNR": kern,
        }
    return {
        "proved": False,
        "surjective_in_general": False,
        "onto_at_listed_primes": bool(onto),
        "rows": rows,
        "note": (
            "Restriction H0→F2^{QR} is not a p-law and is not onto "
            "(ker 2,3,6,12 at p=3,5,7,11).  Do not use it to name "
            "minpoly(D)=X^N+1 by cycle projection."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "H0_quotient_irreducible": False,
        "walsh_general_p": False,
        "restriction_QR_surjective": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "15.604 A,B do not make H0/⟨1⟩ irreducible.  Any G-submodule "
            "is an M-submodule; the unique trivial M-line is the extra "
            "fixed vector, which translations mix.  That mixing argument "
            "is not yet a proof of irreducibility.  Walsh 15.406 E still "
            "needs the xor-slice.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.604  QR/QNR in H0 and ker(D−I)∩H0", flush=True)
    A = theorem_A_qr_in_h0()
    print(f"  A 1_QR/1_QNR in H0 by p mod 4: {A['proved']}", flush=True)
    B = theorem_B_dilation_fixed()
    print(f"  B ker(D−I)∩H0 dim 2: {B['proved']}", flush=True)
    C = theorem_C_restriction_census()
    print(f"  C restriction census (not a p-law): onto={C['onto_at_listed_primes']}", flush=True)
    D = theorem_D_open()
    print(f"  D irred/Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.604",
        "title": "1_QR ∈ H0 iff p≡1 (mod 4); ker(D−I)∩H0 has dim 2",
        "proved": {
            "qr_qnr_in_H0": A["proved"],
            "dilation_fixed_dim2": B["proved"],
            "restriction_QR_surjective": False,
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
        "backend": "serial F2 identities (inherently); line-count p=3..13; GPU unused",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15604.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
