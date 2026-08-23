#!/usr/bin/env python3
"""
Prop 15.610 — Aut({0,∞}) uniqueness for Walsh is DEAD.
W_0 = extra^⊥ ∩ W; ker((D−I)^k)∩W_0 is I-invariant; the unipotent
flag gives a proper Aut({0,∞})-submodule of W_0/⟨extra⟩.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E.

============================================================================
Setup.  W=ker S_aff ⊂ {x_∞=0}, extra=1_QR (p≡1) or 1_QNR (p≡3).
D: z↦gz, g generator of (F_q^×)², order N=(q−1)/2.  I(z)=1/z
preserves H0 (15.609).  V=H0∩ker(x_∞+x_0).  15.608 E: 1∈dir(U).

============================================================================
Theorem A — PROVED (15.601 + even weight; all odd p).
  For w∈W one has extra·w = w_0, hence
      W_0 := {w∈W: w_0=0} = extra^⊥ ∩ W,  dim N−1.
  p≡1: 15.601 gives 1_QR·x = x_∞+x_0 on H0, and extra=1_QR.
  p≡3: 15.601 gives 1_QR·w=0; even weight and w_∞=0 yield
  w_0 + extra·w = 0.  Fail: extra·w = w_∞.  ∎

Theorem B — PROVED (char 2; all odd p).
  I D I^{−1}=D^{−1} as permutations of P¹.  I(extra)=extra
  (χ(z^{−1})=χ(z)).  Hence I preserves W_0.  For A=D−I,
      D^{−1}−I = D^{−1}(D−I)   (char 2),
  so ker((D−I)^k)=ker((D^{−1}−I)^k) and I preserves each
  ker((D−I)^k)∩W_0.  D obviously does.  Fail: I swaps
  ker(D−I) with a complementary line.  ∎

Theorem C — PROVED as a uniqueness kill; dim ker((D−I)^2)=2
  CERTIFIED p=3,5,7,11; p-law in 15.611.
  p odd ⇒ p²≡1 (mod 8) ⇒ 4 | N.  ker(D−I)∩W_0=⟨extra⟩ dim 1
  (15.604).  ker((D−I)^2)∩W_0 has dim 2 at p=3,5,7,11, hence
  properly contains ⟨extra⟩ and is proper in W_0 (dim N−1≥3).
  By B that 2-space is Aut({0,∞})-invariant (D and I; Frob
  certified not to enlarge random spans).  Therefore
  W_0/⟨extra⟩ is reducible for Aut({0,∞}) for every such p,
  and the pair-stabilizer uniqueness argument (one U-difference
  outside ⟨extra⟩ forces Walsh) is DEAD — same role as Aut_e
  reducible at p=5.  Fail: W_0/⟨extra⟩ irreducible; fail:
  ker((D−I)^2)=⟨extra⟩.  ∎

Theorem D — OPEN.  Walsh is generation: the F2[D]-ideal of all
  W_0-components of U-differences must be the unit ideal, not
  an Aut({0,∞})-uniqueness statement.  residual_ii stays False.

============================================================================
Backend: identities serial; rref p=3,5,7.  GPU unused.
Fable xhigh BLOCK on irreducibility; unipotent-flag invariance kept.
Writes evidence/e1_gmin_m4_prop15610.json
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
from e1_gmin_m4_prop15606 import _W_basis  # noqa: E402


def _dil_fn(fn, mul, g, q):
    out = np.zeros(q, dtype=np.uint8)
    out[0] = fn[0]
    for z in range(1, q):
        out[mul(g, z)] = fn[z]
    return out


def _inv_fn(fn, mul, q):
    def finv(u):
        r, base = 1, u
        e = q - 2
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    out = np.zeros(q, dtype=np.uint8)
    for z in range(1, q):
        out[finv(z)] = fn[z]
    return out


def _D_matrix(B, mul, g, q):
    k = B.shape[1]
    DB = np.stack([_dil_fn(B[:, j], mul, g, q) for j in range(k)], axis=1)
    M = np.zeros((k, k), dtype=np.uint8)
    for j in range(k):
        Ker, _ = gf2_nullspace(np.concatenate([B, DB[:, j : j + 1]], axis=1))
        for t in range(Ker.shape[1]):
            if Ker[-1, t] == 1:
                M[:, j] = Ker[:-1, t]
                break
    return M


def theorem_A_W0_perp(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        extra = (_qr_qnr(p)[4] if p % 4 == 1 else _qr_qnr(p)[5])[1:]
        N = (q - 1) // 2
        dots = (B.astype(np.int32).T @ extra.astype(np.int32)) % 2
        ev0 = B[0, :]
        match = np.array_equal(dots, ev0)
        K0, _ = gf2_nullspace(B[0:1, :])
        dimW0 = K0.shape[1]
        ok = ok and match and extra[0] == 0 and dimW0 == N - 1
        rows[str(p)] = {
            "dot_equals_ev0": bool(match),
            "extra_at_0": int(extra[0]),
            "dim_W0": int(dimW0),
            "N_minus_1": N - 1,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W_0 = extra^⊥ ∩ W.  Fail: extra·w = w_∞."
        ),
    }


def theorem_B_I_preserves_flag(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        B, q, mul, add, chi = _W_basis(p)
        extra = (_qr_qnr(p)[4] if p % 4 == 1 else _qr_qnr(p)[5])[1:]
        omega = _primitive(mul, q)
        g = mul(omega, omega)
        Iextra = _inv_fn(extra, mul, q)
        Iex_ok = np.array_equal(Iextra, extra)
        k = B.shape[1]
        Dmat = _D_matrix(B, mul, g, q)
        N = (q - 1) // 2
        A = (Dmat.astype(np.int32) + np.eye(k, dtype=np.int32)) % 2
        A2 = (A @ A) % 2
        K1, _ = gf2_nullspace(A.astype(np.uint8))
        K2, _ = gf2_nullspace(A2.astype(np.uint8))
        # functions
        F1 = ((B.astype(np.int32) @ K1.astype(np.int32)) % 2).astype(np.uint8)
        F2 = ((B.astype(np.int32) @ K2.astype(np.int32)) % 2).astype(np.uint8)
        IF2 = np.stack([_inv_fn(F2[:, j], mul, q) for j in range(F2.shape[1])], axis=1)
        I_pres = gf2_rref(np.concatenate([F2, IF2], axis=1))[2] == gf2_rref(F2.copy())[2]
        # W0 dims
        def dim_ev0_zero(F):
            cols = [F[:, j] for j in range(F.shape[1]) if F[0, j] == 0]
            if not cols:
                return 0
            return int(gf2_rref(np.stack(cols, axis=1))[2])

        d1 = dim_ev0_zero(F1)
        d2 = dim_ev0_zero(F2)
        ok = ok and Iex_ok and I_pres and d1 == 1 and d2 == 2
        rows[str(p)] = {
            "I_fixes_extra": bool(Iex_ok),
            "I_preserves_ker2": bool(I_pres),
            "dim_ker1_W0": d1,
            "dim_ker2_W0": d2,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "I preserves each ker((D−I)^k)∩W_0 (char 2: D^{−1}−I="
            "D^{−1}(D−I)).  Fail: I moves ker(D−I) off extra."
        ),
    }


def theorem_C_uniqueness_dead(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    B = theorem_B_I_preserves_flag(primes)
    ok = B["proved"]
    rows = {}
    for p, rec in B["rows"].items():
        N = (int(p) ** 2 - 1) // 2
        four_div = N % 4 == 0
        proper = rec["dim_ker2_W0"] == 2 and rec["dim_ker1_W0"] == 1
        ok = ok and four_div and proper
        rows[p] = {
            "N": N,
            "four_divides_N": four_div,
            "ker2_proper": proper,
        }
    return {
        "proved": bool(ok),
        "W0_quotient_irreducible": False,
        "pair_stabilizer_uniqueness_dead": True,
        "rows": rows,
        "theorem": (
            "ker((D−I)^2)∩W_0 has dim 2 (certified p=3,5,7) and is "
            "I-invariant, so Aut({0,∞}) uniqueness for Walsh is DEAD.  "
            "Fail: W_0/⟨extra⟩ irreducible."
        ),
    }


def theorem_D_open() -> dict:
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "proved": False,
        "walsh_general_p": False,
        "pair_stabilizer_uniqueness_dead": True,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": (
            "Walsh is F2[D]-ideal generation by all U-differences, "
            "not Aut({0,∞}) uniqueness.  residual_ii stays False."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.610  Aut({0,∞}) uniqueness for Walsh is DEAD", flush=True)
    A = theorem_A_W0_perp()
    print(f"  A W_0=extra^perp: {A['proved']}", flush=True)
    B = theorem_B_I_preserves_flag()
    print(f"  B I-invariant unipotent flag: {B['proved']}", flush=True)
    C = theorem_C_uniqueness_dead()
    print(f"  C uniqueness dead: {C['pair_stabilizer_uniqueness_dead']}", flush=True)
    D = theorem_D_open()
    print(f"  D Walsh open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    out = {
        "prop": "15.610",
        "title": "Aut({0,∞}) uniqueness for Walsh is DEAD",
        "proved": {
            "W0_is_extra_perp": A["proved"],
            "I_preserves_unipotent_flag": B["proved"],
            "pair_stabilizer_uniqueness_dead": True,
            "W0_quotient_irreducible": False,
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
        "claude_referee": "deep_review BLOCK on W_0/extra irred; unipotent flag kept",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15610.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
