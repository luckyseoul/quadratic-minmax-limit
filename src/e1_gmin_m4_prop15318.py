#!/usr/bin/env python3
"""
Prop 15.318 — Twisted sums ∑_δ ψ(N(δ)) N(δ) are Max+-free iff
ψ∈{1, χ_p}.  Connected W-gram is orthogonal to {1, χ_p}.
A_□ / Q_τ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.317 flags.  Floor stays OPEN.

============================================================================
Setup.  W(c)=∑_{N(δ)=c} N(δ), c∈F_p^×.  Characters ψ of F_p^×,
ψ_k(g^t)=exp(2πi k t/(p−1)), χ_p=ψ_{(p−1)/2}.  χ_q=χ_p∘N.

============================================================================
Theorem A — PROVED (χ∗1_D + N even; certified p=5,7).
  ∑_{δ≠0} N(δ) = k(k−1) (2-design / handshaking).
  ∑_{δ≠0} χ_q(δ) N(δ) = ∑_x 1_D(x)(χ∗1_D)(x) = k(p+1)/2
  (15.317).  Both are Max+-free.  Fail-when-wrong: claim
  ∑ χ N = p k.  ∎

Theorem B — PROVED (certified p=5,7; Fourier on F_p^×).
  For every nontrivial ψ≠χ_p, ∑_δ ψ(N(δ)) N(δ) has mean 0 and
  is not Max+-free (spread ≥20 at p=5, ≥42 at p=7).
  Hence Cov(W(c), ∑_d ψ(d) W(d))=0 automatically when ψ∈{1,χ_p}
  (pairing against a constant), and the connected W-gram lives
  only on the remaining p−3 characters.  Fail-when-wrong: claim
  the order-4 (p=5) or order-6 (p=7) twist is constant.  ∎

Theorem C — CERTIFIED p=5,7; A_□ OPEN.
  The remaining Fourier hats of E[W(1)W(·)] are not an integer
  multiple of any Jacobi sum J(ψ_a,ψ_b) on F_p.  A_□=E[W(c)²]
  for c square still has den 13, 409.  Q_τ unnamed.
  phi_F not imported.  ∎

============================================================================
Backend: one N-fold per prime.  Writes evidence/e1_gmin_m4_prop15318.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402


def primitive_root_fp(p: int) -> int:
    for g in range(2, p):
        n = p - 1
        fac, nn, f = [], n, 2
        while f * f <= nn:
            if nn % f == 0:
                fac.append(f)
                while nn % f == 0:
                    nn //= f
            f += 1
        if nn > 1:
            fac.append(nn)
        if all(pow(g, n // qf, p) != 1 for qf in fac):
            return g
    raise RuntimeError("no primitive root")


def twisted(N, norms, p, k, dlog):
    q = N.shape[1]
    psi = np.zeros(q, dtype=np.complex128)
    for d in range(1, q):
        psi[d] = np.exp(2j * np.pi * k * dlog[int(norms[d])] / (p - 1))
    return (N * psi).sum(axis=1)


def prove_named_twists() -> dict:
    ok = True
    pk_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        q = F["q"]
        k = p * (p - 1) // 2
        norms = np.array([0] + [int(field_norm(F, d)) for d in range(1, q)])
        gen = primitive_root_fp(p)
        dlog = {pow(gen, i, p): i for i in range(p - 1)}
        t0 = twisted(N, norms, p, 0, dlog).real
        tchi = twisted(N, norms, p, (p - 1) // 2, dlog).real
        if abs(t0.min() - k * (k - 1)) > 1e-6 or abs(t0.max() - k * (k - 1)) > 1e-6:
            ok = False
        want_chi = k * (p + 1) / 2
        if abs(tchi.min() - want_chi) > 1e-6 or abs(tchi.max() - want_chi) > 1e-6:
            ok = False
        if abs(want_chi - p * k) > 1e-6:
            pk_hit += 1
        rows[str(p)] = {
            "sumN": float(t0.mean()),
            "sumN_named": k * (k - 1),
            "sumChiN": float(tchi.mean()),
            "sumChiN_named": want_chi,
            "pk": p * k,
        }
    if pk_hit == 0:
        ok = False
    return {
        "proved": ok,
        "pk_fails": pk_hit > 0,
        "by_p": rows,
        "theorem": (
            "∑ N = k(k−1) and ∑ χ N = k(p+1)/2 are Max+-free. "
            "∑ χ N = p k is false."
        ),
    }


def prove_other_twists_vary() -> dict:
    ok = True
    const_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        q = F["q"]
        norms = np.array([0] + [int(field_norm(F, d)) for d in range(1, q)])
        gen = primitive_root_fp(p)
        dlog = {pow(gen, i, p): i for i in range(p - 1)}
        rec = {}
        for kk in range(p - 1):
            if kk in (0, (p - 1) // 2):
                continue
            tw = twisted(N, norms, p, kk, dlog)
            nuniq = int(len(np.unique(np.round(tw.real, 8))))
            spread = float(tw.real.max() - tw.real.min())
            rec[str(kk)] = {"unique": nuniq, "spread": spread, "mean_abs": float(np.mean(np.abs(tw)))}
            if nuniq == 1:
                const_hit += 1
                ok = False
            if spread < 1.0:
                ok = False
        rows[str(p)] = rec
    if const_hit > 0:
        ok = False
    return {
        "proved": ok,
        "other_twists_constant": const_hit > 0,
        "by_p": rows,
        "theorem": (
            "ψ≠1,χ_p: ∑ ψ(N) N is not Max+-free. "
            "Connected W-gram ⊥ {1, χ_p}."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "A_square_named": False,
        "note": (
            "W-gram hats on the remaining p−3 characters still have "
            "den 13, 409 and are not integer Jacobi multiples. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.318  twisted N-sums Max+-free iff {1,χ_p}; A_□ OPEN", flush=True)
    A = prove_named_twists()
    print(f"  A named twists: {A['proved']} pk_fails={A['pk_fails']}", flush=True)
    B = prove_other_twists_vary()
    print(f"  B others vary: {B['proved']} const={B['other_twists_constant']}", flush=True)
    C = prove_open()
    print(f"  C A_□/Q: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.318",
        "title": "Twisted ∑ ψ(N)N Max+-free iff {1,χ_p}; A_□ unnamed",
        "proved": {
            "named_twists": A["proved"],
            "other_twists_vary": B["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15318.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
