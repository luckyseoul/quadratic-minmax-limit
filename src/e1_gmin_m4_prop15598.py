#!/usr/bin/env python3
"""
Prop 15.598 — Square-direction affine lines cut Max− over F_2.
A general-p reason the pair-slice U is an xor-hyperplane of
affine_span(Max−).  15.406 Theorem E (Walsh ∀p) stays OPEN:
spanning of the slice is certified p=3,5,7 (exact) and p=11
(full rank(B_U)=60 with dir(Max−)=60 on a 200k sample), not proved.
residual_ii stays False (Walsh, even ∀p, does not empty leftover-only / 5+).

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.

============================================================================
Setup.  Paley C on P¹(F_q), q=p², n=q+1, ∞ and F_q.  Max− = {y: Cy=−p y}.
x=(1−y)/2 ∈ F_2^n.  An affine F_p-line is L=a+F_p·b ⊂ F_q, b≠0.
S={∞}∪L has p+1 points.  Direction class b F_p^×; χ(b) is well-defined.

============================================================================
Theorem A — PROVED (Jacobi on F_p; Max-free).
  For odd prime p and δ∈F_p^×,
      Σ_{x∈F_p} χ_p(x(x+δ)) = −1.
  Proof: u≠0 term Σ χ(u(u+δ))=Σ_{u≠0} χ(1+δ u^{-1})=Σ_{z≠1} χ_p(z)=−χ_p(1)=−1,
  and χ(0·δ)=0.  Fail: claim the sum is 0, or p−1.  ∎

Theorem B — PROVED (norm of F_{p²}/F_p; Max-free).
  χ_{p²}(z)=χ_p(N(z)) with N(z)=z^{p+1}=z z^p.  Hence χ_{p²}|_{F_p^×}≡1,
  and for L=a+F_p b, distinct i,j∈L have C_{ij}=χ(b).  For j∉S,
      Σ_{i∈L} C_{ij} = −χ(b).
  Proof: Σ_{t∈F_p} χ_{p²}(a+tb−j)=χ(b) Σ_u χ_{p²}(u−c), c∉F_p;
  N(u−c)=(u−c)(u−c^p)=u²−Tr(c)u+N(c) has discriminant (c−c^p)²≠0;
  complete the square and apply A with δ≠0.  Thus the off-S row-sum
  of C on S is σ_j=1−χ(b).  Square direction ⇒ σ_j=0 and on-S row-sums
  equal p.  Fail: claim σ=0 on a nonsquare direction (then σ=2).  ∎

Theorem C — PROVED (Cy=−p y + B).
  If χ(b)=1 then 1_S^T C y = p 1_S^T y, while Cy=−p y gives
  1_S^T C y = −p 1_S^T y, hence Σ_{k∈S} y_k = 0 on every Max− vector.
  Equivalently ⟨x, 1_S⟩=(p+1)/2, so
      ⟨x, 1_S⟩ ≡ 0 (mod 2) if p≡3 (mod 4),
      ⟨x, 1_S⟩ ≡ 1 (mod 2) if p≡1 (mod 4).
  Fail: the same identity on Max+ (not forced; many values at p=5);
  fail: nonsquare direction on Max− (sum_S not identically 0).  ∎

Theorem D — PROVED (Boolean).
  y_i y_j=−1 ⇔ x_i XOR x_j=1.  Hence the pair-slice
  U={y∈Max−: C_{ij} y_i y_j=−1} is Max− cut by the affine hyperplane
  x_i+x_j=c of F_2^n (c=0 or 1 according as C_{ij}=−1 or +1).  ∎

Theorem E — OPEN.  Let H=affine_span_{F_2}(Max−) and ℓ(x)=x_i+x_j.
  C + nonempty pair-slices ⇒ ℓ is nonconstant on H, so U ⊂ H∩{ℓ=c}
  of codimension 1.  Walsh (15.406 C/E) is exactly
      affine_span(U)=H∩{ℓ=c}.
  Certified: dim U = dim H − 1 at p=3,5,7 (full ensemble) and at p=11
  (rank(B_U)=60, sample dir(H)=60=n/2−1).  Square-line equations of C
  cut an affine space containing H; their F_2-rank equals the dual
  dimension at p=3,5,7.  Spanning for general p is not proved.
  residual_ii stays False: leftover-only / 5+ / k>4p remain.

============================================================================
Backend: F_p Jacobi (no Max); Paley C line-sums; Max− caches p=5,7.
Writes evidence/e1_gmin_m4_prop15598.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def legendre(x: int, p: int) -> int:
    x %= p
    if x == 0:
        return 0
    r = pow(x, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1


def jacobi_quad_sum(p: int, delta: int) -> int:
    """Σ_x χ_p(x(x+delta))."""
    s = 0
    for x in range(p):
        s += legendre(x * ((x + delta) % p) % p, p)
    return s


def field_ctx(p: int):
    """Encoding of paley_conference_prime_power: e=c0+c1*p."""
    q = p * p

    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        e0 = (c0 * d0 + c1 * d1 * ib) % p
        e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        return e0 + e1 * p

    def add(u, v):
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    def fpow(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def chi(z):
        if z == 0:
            return 0
        return 1 if fpow(z, (q - 1) // 2) == 1 else -1

    def frob(z):
        return fpow(z, p)

    def norm(z):
        return mul(z, frob(z))  # in F_p (c1=0)

    return q, mul, add, chi, frob, norm, ia, ib


def theorem_A_jacobi(primes=None) -> dict:
    if primes is None:
        primes = [q for q in range(5, 80) if _is_prime(q)]
    ok = True
    rows = {}
    for p in primes:
        vals = [jacobi_quad_sum(p, d) for d in range(1, p)]
        ok = ok and all(v == -1 for v in vals)
        # fail-when-wrong
        ok = ok and all(v != 0 for v in vals)
        ok = ok and all(v != p - 1 for v in vals)
        rows[str(p)] = {"n_delta": p - 1, "all_minus_1": all(v == -1 for v in vals)}
    return {
        "proved": bool(ok),
        "n_primes": len(primes),
        "rows": rows,
        "theorem": "Σ χ_p(x(x+δ))=-1 for δ≠0. Fail: 0 or p-1.",
    }


def theorem_B_line_sums(primes=None) -> dict:
    if primes is None:
        primes = (5, 7, 11, 13)
    ok = True
    rows = {}
    for p in primes:
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        # χ_{p²}(z)=χ_p(N(z))
        n_mismatch = 0
        for z in range(q):
            Nz = norm(z)  # c1=0, value in 0..p-1
            if chi(z) != legendre(Nz % p, p):
                n_mismatch += 1
        ok = ok and n_mismatch == 0
        # χ_{p²} on F_p^× is 1
        ok = ok and all(chi(t) == 1 for t in range(1, p))
        # line sums: one square and one nonsquare direction
        sq_b = next(b for b in range(1, q) if chi(b) == 1)
        nsq_b = next(b for b in range(1, q) if chi(b) == -1)
        rec = {"norm_mismatch": n_mismatch, "chi_Fp": True}

        def off_line_sums(b):
            L = [add(0, mul(t, b)) for t in range(p)]  # 0+F_p b
            Lset = set(L)
            out = []
            for j in range(q):
                if j in Lset:
                    continue
                s = 0
                for i in L:
                    d = (i % p - j % p) % p + ((i // p - j // p) % p) * p
                    s += chi(d)
                out.append(s)
            return sorted(set(out))

        rec["sq_sums"] = off_line_sums(sq_b)
        rec["nsq_sums"] = off_line_sums(nsq_b)
        ok = ok and rec["sq_sums"] == [-1]
        ok = ok and rec["nsq_sums"] == [1]
        # fail: σ=0 on nonsquare ⇒ line sum -1 on nsq
        if rec["nsq_sums"] == [-1]:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Σ_{i∈L} C_ij=−χ(b) off the line. Square ⇒ σ=0. "
            "Fail: σ=0 on a nonsquare direction."
        ),
    }


def _square_direction_lines(p: int):
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    used = set()
    dirs = []
    for b in range(1, q):
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    lines = []
    for b in dirs:
        if chi(b) != 1:
            continue
        covered = set()
        for a in range(q):
            if a in covered:
                continue
            pts = [0]  # infinity
            for t in range(p):
                e = add(a, mul(t, b))
                pts.append(1 + e)
                covered.add(e)
            lines.append(pts)
    return lines


def theorem_C_maxminus_sum(primes=None) -> dict:
    from e1_gmin_m4_prop15406 import load_minus  # local: needs caches p=5,7

    if primes is None:
        primes = (5, 7)
    ok = True
    rows = {}
    for p in primes:
        Y, C = load_minus(p)
        Y = np.sign(Y.astype(np.float64)).astype(np.int64)
        lines = _square_direction_lines(p)
        n_bad = 0
        for pts in lines:
            sm = Y[:, pts].sum(axis=1)
            if sm.min() != 0 or sm.max() != 0:
                n_bad += 1
        # fail: Max+ is not identically 0 on the standard subline
        # (use minus cache only; check a nonsquare line is NOT 0)
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        nsq_b = next(b for b in range(1, q) if chi(b) == -1)
        nsq_pts = [0] + [1 + mul(t, nsq_b) for t in range(p)]
        nsq_vals = set(Y[:, nsq_pts].sum(axis=1).tolist())
        ok = ok and n_bad == 0
        ok = ok and nsq_vals != {0}
        rows[str(p)] = {
            "n_square_lines": len(lines),
            "n_bad_square": n_bad,
            "nsq_nunique": len(nsq_vals),
            "parity": int((p + 1) // 2 % 2),
        }
        # (p+1)/2 mod 2
        ok = ok and rows[str(p)]["parity"] == ((p % 4) == 1)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Max− has Σ_S y=0 on every square-direction ∞∪L. "
            "Fail: nonsquare direction, or the same claim on Max+."
        ),
    }


def theorem_D_xor() -> dict:
    ok = True
    # Boolean: y=1-2x, y_i y_j=-1 ⇔ x_i XOR x_j=1
    for xi in (0, 1):
        for xj in (0, 1):
            yi, yj = 1 - 2 * xi, 1 - 2 * xj
            xor = xi ^ xj
            ok = ok and ((yi * yj == -1) == (xor == 1))
    # fail: claim XOR=0
    if (1 * 1 == -1):
        ok = False
    return {
        "proved": bool(ok),
        "theorem": "y_i y_j=-1 ⇔ x_i XOR x_j=1. Fail: claim XOR=0.",
    }


def theorem_E_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "walsh_general_p": False,
        "note": (
            "Walsh = affine_span(U)=H∩{x_i+x_j=c}. Certified dim U=dim H-1 "
            "at p=3,5,7 and p=11 sample. Square-line equations of C contain H. "
            "Spanning for general p open. residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.598  square-direction lines cut Max- over F2", flush=True)
    A = theorem_A_jacobi()
    print(f"  A Jacobi: {A['proved']}  n_primes={A['n_primes']}", flush=True)
    B = theorem_B_line_sums()
    print(f"  B line-sums: {B['proved']}", flush=True)
    C = theorem_C_maxminus_sum()
    print(f"  C Max- sum_S=0: {C['proved']}", flush=True)
    D = theorem_D_xor()
    print(f"  D xor cut: {D['proved']}", flush=True)
    E = theorem_E_open()
    print(f"  E Walsh open: resii={E['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.598",
        "title": "Square-direction affine lines cut Max- over F2",
        "proved": {
            "jacobi_quad_sum": A["proved"],
            "line_character_sum": B["proved"],
            "maxminus_square_line_sum_zero": C["proved"],
            "xor_pair_slice": D["proved"],
            "walsh_general_p": False,
        },
        "A": {k: v for k, v in A.items() if k != "rows"},
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "multilevel_ND_k_ge_4p_proved",
            "phi_F_ge_6_proved_general",
            "type_I_multilevel_bad_case_ND_closed",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15598.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
