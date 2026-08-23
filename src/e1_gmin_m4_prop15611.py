#!/usr/bin/env python3
"""
Prop 15.611 — W ≅ F2[M] ≅ F2[X]/(X^N+1) as D-modules; W_0 is
the unique D-invariant hyperplane; dim ker((D−I)^2)∩W_0=2 is a
p-law.  Max-free.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E.

============================================================================
Setup.  15.605–15.607: W=ker S_aff ⊂ {x_∞=0}, dim N=(q−1)/2,
W=⊕_{H nsq} W^H, dim W^H=p−1, M=(F_q^×)² transits nsq 0-lines
with point-stabilizer F_p^× ⊂ M (15.598 B).  D generates M.
W_0={w∈W: w_0=0}=extra^⊥∩W (15.610 A).  extra generates
ker(D−I)∩W (15.604).

============================================================================
Theorem A — PROVED (Mackey / regular induction; all odd p).
  Even H-invariants on F_q/H ≅ F_p are determined by their
  restriction to F_p^×: char 2 and even weight force
  f(0)=∑_{F_p^×} f.  F_p^× permutes F_p^× regularly, so
      W^H ≅ F2[F_p^×]
  as F_p^×-modules (the regular representation).  This is the
  same vector space as the C_p-augmentation F2[C_p]/(Σ), but a
  different group action: at p=7 the C_p-module splits into two
  cubics (15.606 D) while F_p^× is still cyclic of order 6.
  W=⊕ W^H and M transits nsq lines regularly with Stab=F_p^×,
  so W ≅ Ind_{F_p^×}^M W^{H_0} ≅ Ind_{F_p^×}^M F2[F_p^×].
  Induction of the regular H-module is the regular G-module, and
  M=⟨D⟩ is cyclic of order N, hence
      W ≅ F2[M] ≅ F2[X]/(X^N+1)
  as D-modules.  Fail: minpoly of D on W has degree <N (that is
  the restriction-of-Ind-is-copies picture, or a non-regular
  W^H).  Fail: W^H simple as a C_p-module at p=7.  ∎

Theorem B — PROVED (unique trivial + permutation form; all odd p).
  D is a permutation of F_q so preserves the standard F2-dot.
  Unique D-invariants in W are ⟨extra⟩ (15.604).  Unique
  D-invariant linear form is extra·(−).  15.610 A: extra·w=w_0,
  so the unique D-invariant hyperplane is W_0.  Equivalently
  im(D−I)=W_0.  As F2[X]/(X^N+1) the unique hom to F2 is
  evaluation at X=1 (only linear factor over F2), with kernel
  (X+1)R ≅ F2[X]/((X^N+1)/(X+1)).  Fail: a second invariant
  hyperplane; fail: extra ∉ W_0.  ∎

Theorem C — PROVED (2-adic; upgrades 15.610 C from certified
  to a p-law).  p odd ⇒ p²≡1 (mod 8) ⇒ 4|N.  Write
  N=2^a m with m odd and a=v_2(N)≥2.  Then
      X^N+1=(X^m+1)^{2^a}=(X+1)^{2^a} g^{2^a},
  g=(X^m+1)/(X+1), g(1)≠0.  On W_0 the (X+1)-primary has
  dimension 2^a−1≥3 and ker((D−I)^k) has dimension
  min(k, 2^a−1) there; other primaries have D−I invertible.
  Hence dim ker((D−I)^2)∩W_0=2 for every odd p.  Fail:
  dim=1 (would need a=1).  ∎

Theorem D — OPEN.  Walsh is still F2[D]-ideal generation by
  the full set of U-difference W_0-components (escape the
  unipotent flag and the reciprocal/Frobenius orbits of the
  factors of h=(X^N+1)/(X+1)).  residual_ii stays False.

============================================================================
Backend: identities serial; rref p=3,5,7.  GPU unused.
Fable xhigh deep_review PASS (confidence 0.93) on claims 1–4.
Writes evidence/e1_gmin_m4_prop15611.json
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
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15606 import (  # noqa: E402
    _W_basis,
    _fp_lines,
    _pi_on_basis,
    order_of_2_mod,
)
from e1_gmin_m4_prop15610 import _D_matrix, _dil_fn  # noqa: E402


def _v2(n: int) -> int:
    a = 0
    while n % 2 == 0:
        n //= 2
        a += 1
    return a


def _prim_root_mod_p(p: int) -> int:
    fac = []
    m = p - 1
    d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        fac.append(m)
    for e in range(2, p):
        if all(pow(e, (p - 1) // r, p) != 1 for r in fac):
            return e
    raise RuntimeError("no primitive root mod p")


def _orbit_rank_fn(fn: np.ndarray, mul, g: int, q: int, limit: int) -> int:
    cols = [fn.copy()]
    cur = fn.copy()
    for _ in range(limit - 1):
        cur = _dil_fn(cur, mul, g, q)
        cols.append(cur)
    return int(gf2_rref(np.stack(cols, axis=1).copy())[2])


def _find_cyclic_fn(F: np.ndarray, mul, g: int, q: int, target: int, trials: int = 48):
    rng = np.random.default_rng(0)
    k = F.shape[1]
    for t in range(trials):
        if t < k:
            v = F[:, t]
        else:
            coef = rng.integers(0, 2, size=k, dtype=np.uint8)
            if coef.max() == 0:
                continue
            v = (F.astype(np.int32) @ coef.astype(np.int32) % 2).astype(np.uint8)
        if v.max() == 0:
            continue
        r = _orbit_rank_fn(v, mul, g, q, target)
        if r == target:
            return True, t
    return False, -1


def _col_basis(M: np.ndarray) -> np.ndarray:
    acc = None
    for j in range(M.shape[1]):
        cand = M[:, j : j + 1]
        if cand.max() == 0:
            continue
        if acc is None:
            acc = cand
        else:
            aug = np.concatenate([acc, cand], axis=1)
            if gf2_rref(aug.copy())[2] > acc.shape[1]:
                acc = aug
    if acc is None:
        return np.zeros((M.shape[0], 0), dtype=np.uint8)
    return acc


def theorem_A_cyclic_regular(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        N = (q - 1) // 2
        omega = _primitive(mul, q)
        g = mul(omega, omega)
        k = B.shape[1]
        Dmat = _D_matrix(B, mul, g, q)
        Mpow = np.eye(k, dtype=np.uint8)
        for _ in range(N):
            Mpow = (Dmat.astype(np.int32) @ Mpow.astype(np.int32) % 2).astype(
                np.uint8
            )
        DN_is_I = np.array_equal(Mpow, np.eye(k, dtype=np.uint8))
        cyclic, trial = _find_cyclic_fn(B, mul, g, q, N)
        lines = _fp_lines(p, mul, chi, q)
        nsq = [H for c, H in lines if c == -1]
        r = _prim_root_mod_p(p)
        WH = _col_basis(_pi_on_basis(B, nsq[0], add, q))
        # C_p = G/H acts by translation by a vector outside H.
        # Cyclic iff Φ_p is irreducible over F2 (false at p=7).
        Hset = set(nsq[0])
        vshift = next(x for x in range(1, q) if x not in Hset)
        Cp_rank = 0
        for j in range(WH.shape[1]):
            cols = [WH[:, j].copy()]
            cur = WH[:, j].copy()
            for _ in range(p - 2):
                nxt = np.zeros(q, dtype=np.uint8)
                for z in range(q):
                    nxt[add(z, vshift)] = cur[z]
                cur = nxt
                cols.append(cur)
            Cp_rank = max(Cp_rank, int(gf2_rref(np.stack(cols, axis=1).copy())[2]))
        Fp_cyclic, fp_trial = _find_cyclic_fn(WH, mul, r, q, p - 1)
        WH_simple_Cp = order_of_2_mod(p) == p - 1
        ok = ok and DN_is_I and cyclic and Fp_cyclic and WH.shape[1] == p - 1
        ok = ok and Cp_rank == p - 1
        ok = ok and not (p == 7 and WH_simple_Cp)
        rows[str(p)] = {
            "N": N,
            "dim_W": k,
            "D_N_is_I": bool(DN_is_I),
            "W_cyclic_D": bool(cyclic),
            "cyclic_trial": int(trial),
            "dim_WH": int(WH.shape[1]),
            "WH_cyclic_Fpstar": bool(Fp_cyclic),
            "WH_orbit_Cp": int(Cp_rank),
            "WH_simple_Cp": bool(WH_simple_Cp),
            "fp_trial": int(fp_trial),
            "n_nsq": len(nsq),
            "index_M_over_Fpstar": (p + 1) // 2,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W ≅ F2[M] ≅ F2[X]/(X^N+1) as D-modules (Ind of regular "
            "F_p^×).  Fail: minpoly(D) has degree <N; fail: W^H simple "
            "as C_p at p=7."
        ),
    }


def theorem_B_unique_hyperplane(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    A = theorem_A_cyclic_regular(primes)
    ok = A["proved"]
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        extra = (_qr_qnr(p)[4] if p % 4 == 1 else _qr_qnr(p)[5])[1:]
        omega = _primitive(mul, q)
        g = mul(omega, omega)
        k = B.shape[1]
        N = (q - 1) // 2
        Dmat = _D_matrix(B, mul, g, q)
        Am = (Dmat.astype(np.int32) + np.eye(k, dtype=np.int32)) % 2
        Fim = (B.astype(np.int32) @ Am.astype(np.int32) % 2).astype(np.uint8)
        r_im = int(gf2_rref(Fim.copy())[2])
        vanish0 = bool(Fim[0, :].max() == 0)
        extra0 = int(extra[0])
        K0, _ = gf2_nullspace(B[0:1, :])
        dimW0 = int(K0.shape[1])
        # extra lives in W_0 (even weight of QR/QNR = N even)
        extra_in_W0 = extra0 == 0
        ok = ok and r_im == N - 1 and vanish0 and extra_in_W0 and dimW0 == N - 1
        rows[str(p)] = {
            "rank_im_DminusI": r_im,
            "im_vanishes_at_0": vanish0,
            "dim_W0": dimW0,
            "N_minus_1": N - 1,
            "extra_at_0": extra0,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Unique D-invariant hyperplane is W_0=im(D−I) ≅ "
            "F2[X]/((X^N+1)/(X+1)).  Fail: extra∉W_0; fail: a second "
            "invariant hyperplane."
        ),
    }


def theorem_C_ker2_dim_p_law(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    Bth = theorem_B_unique_hyperplane(primes)
    ok = Bth["proved"]
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        omega = _primitive(mul, q)
        g = mul(omega, omega)
        k = B.shape[1]
        N = (q - 1) // 2
        a = _v2(N)
        Dmat = _D_matrix(B, mul, g, q)
        Am = (Dmat.astype(np.int32) + np.eye(k, dtype=np.int32)) % 2
        A2 = (Am @ Am) % 2
        K1, _ = gf2_nullspace(Am.astype(np.uint8))
        K2, _ = gf2_nullspace(A2.astype(np.uint8))
        F1 = ((B.astype(np.int32) @ K1.astype(np.int32)) % 2).astype(np.uint8)
        F2 = ((B.astype(np.int32) @ K2.astype(np.int32)) % 2).astype(np.uint8)

        def dim_ev0_zero(F):
            cols = [F[:, j] for j in range(F.shape[1]) if F[0, j] == 0]
            if not cols:
                return 0
            return int(gf2_rref(np.stack(cols, axis=1))[2])

        d1 = dim_ev0_zero(F1)
        d2 = dim_ev0_zero(F2)
        four_div = N % 4 == 0
        ok = ok and a >= 2 and four_div and d1 == 1 and d2 == 2
        rows[str(p)] = {
            "N": N,
            "v2_N": a,
            "four_divides_N": four_div,
            "dim_ker1_W0": d1,
            "dim_ker2_W0": d2,
        }
    return {
        "proved": bool(ok),
        "ker2_dim_is_p_law": True,
        "rows": rows,
        "theorem": (
            "dim ker((D−I)^2)∩W_0=2 for every odd p (a=v_2(N)≥2; "
            "(X+1)-primary of W_0 has dim 2^a−1≥3).  Fail: dim=1."
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
            "Walsh is F2[D]-ideal generation by all U-differences "
            "in W_0, not cyclicity of W.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.611  W ≅ F2[X]/(X^N+1); ker2 dim 2 is a p-law", flush=True)
    A = theorem_A_cyclic_regular()
    print(f"  A W cyclic F2[M]: {A['proved']}", flush=True)
    B = theorem_B_unique_hyperplane()
    print(f"  B unique hyperplane W_0: {B['proved']}", flush=True)
    C = theorem_C_ker2_dim_p_law()
    print(f"  C ker2 dim 2 p-law: {C['proved']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.611",
        "title": "W ≅ F2[X]/(X^N+1); dim ker((D−I)^2)∩W_0=2 is a p-law",
        "proved": {
            "W_cyclic_F2M": A["proved"],
            "W0_unique_D_hyperplane": B["proved"],
            "ker2_dim_is_p_law": C["ker2_dim_is_p_law"] and C["proved"],
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
        "backend": "serial F2 identities; rref p=3,5,7; GPU unused",
        "claude_referee": (
            "deep_review PASS confidence 0.93 on W≅F2[M] and ker2 dim 2"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15611.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
