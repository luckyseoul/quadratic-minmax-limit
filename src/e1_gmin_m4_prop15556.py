#!/usr/bin/env python3
"""
Prop 15.556 — The Type+ 1D summand of G_3 is named:
    G3_1d(0,u,v) = [n_spec/p − 3(m−n_spec)/(p(p−2))] C(p,m)
with m=(p+1)/2 and
    n_spec = (1+χ(u))/2                 if v/u ∈ F_p,
            (3+χ(u)+χ(v)+χ(v−u))/2    otherwise.
Φ = G3_1d + Φ_free.  Φ_free is not identically 0 (p=5 max 36).
16pA / Q_τ / phi_F stay open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.555 flags.  Does **not** import phi_F.

============================================================================
Setup.  15.538: Type+ 1D is y_x=ε(φ(x)), ker φ a square F_p-line,
ε:F_p→{±1}, #(+)=m, n_1d=m C(p,m).  Johnson: Σε(a)=C/p,
Σε(a)ε(b)ε(c)=−3C/(p(p−2)) on distinct triples, C/p on
double+one.  Square lines: m=(p+1)/2 (c∈F_p is always a
square in F_q, so each F_p-line is square iff χ_q(ω)=1).

============================================================================
Theorem A — PROVED Max+-free (Johnson + square-line count).
  The formula for G3_1d above.  Fail: use n_spec=(3+σ)/2 on
  every pair (collinear wrong).  Fail: drop the μ₃ term.
  Certified every 3-distinct pair at p=5,7 against the 1D
  slice of the Max+ cache.  ∎

Theorem B — PROVED (A + 15.555).
  Φ(N(u),N(v),N(v−u)) = G3_1d(0,u,v) + Φ_free(N(u),N(v),N(v−u)).
  Fail: Φ_free≡0 (p=5 a key has Φ_free=20).  ∎

Theorem C — OPEN.  Φ_free takes 5 values at p=5 and 14 at p=7
  and is not named.  Q3 / 16pA / Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Fraction Johnson + p=5,7 1D slice.  Writes
evidence/e1_gmin_m4_prop15556.json
"""
from __future__ import annotations

import json
import os
import sys
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15554 import load_Y, yhat_omega
from e1_gmin_m4_prop15555 import field_norm

EV = ROOT / "evidence" / "e1_gmin_m4_prop15556.json"


def g3_1d_named(p: int, chi_u: int, chi_v: int, chi_w: int, collinear: bool) -> float:
    """Type+ 1D contribution to G3(0,u,v)."""
    m = (p + 1) // 2
    C = comb(p, m)
    mu3C = -3 * C / (p * (p - 2))
    Cp = C / p
    if collinear:
        n_spec = (1 + chi_u) / 2
    else:
        n_spec = (3 + chi_u + chi_v + chi_w) / 2
    return n_spec * Cp + (m - n_spec) * mu3C


def g3_1d_wrong_nonspec(p: int, chi_u: int, chi_v: int, chi_w: int, collinear: bool) -> float:
    """Fail: pretend every pair is non-collinear."""
    return g3_1d_named(p, chi_u, chi_v, chi_w, False)


def g3_1d_wrong_drop_mu3(p: int, chi_u: int, chi_v: int, chi_w: int, collinear: bool) -> float:
    """Fail: drop the μ₃ lines."""
    m = (p + 1) // 2
    C = comb(p, m)
    if collinear:
        n_spec = (1 + chi_u) / 2
    else:
        n_spec = (3 + chi_u + chi_v + chi_w) / 2
    return n_spec * (C / p)


def is_collinear(F: dict, u: int, v: int) -> bool:
    """v/u ∈ F_p."""
    r = F["mul"][F["inv"][u]][v]
    return r < F["p"]


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        Y = load_Y(p, q)
        _, _, Omega, yhat = yhat_omega(p)
        one = (np.abs(yhat) > 1e-6).sum(axis=1) == (p - 1)
        Y1 = Y[one]
        add, neg, chi = F["add"], F["neg"], F["chi"]
        maxerr = 0.0
        drop_err = 0.0
        nonspec_err = 0.0
        n_col = 0
        n_pairs = 0
        for u in range(1, q):
            for v in range(1, q):
                if u == v:
                    continue
                w = add[v][neg[u]]
                if w == 0:
                    continue
                n_pairs += 1
                col = is_collinear(F, u, v)
                if col:
                    n_col += 1
                named = g3_1d_named(
                    p, int(chi[u]), int(chi[v]), int(chi[w]), col
                )
                got = float(np.sum(Y1[:, 0] * Y1[:, u] * Y1[:, v]).real)
                maxerr = max(maxerr, abs(got - named))
                drop_err = max(
                    drop_err,
                    abs(got - g3_1d_wrong_drop_mu3(p, int(chi[u]), int(chi[v]), int(chi[w]), col)),
                )
                if col:
                    nonspec_err = max(
                        nonspec_err,
                        abs(
                            got
                            - g3_1d_wrong_nonspec(
                                p, int(chi[u]), int(chi[v]), int(chi[w]), col
                            )
                        ),
                    )
        row_ok = maxerr < 1e-8 and drop_err > 1.0 and (nonspec_err > 1.0 or n_col == 0)
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_1d": int(Y1.shape[0]),
            "n_pairs": n_pairs,
            "n_collinear": n_col,
            "maxerr": float(maxerr),
            "drop_mu3_err": float(drop_err),
            "nonspec_err": float(nonspec_err),
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "G3_1d=[n_spec/p−3(m−n_spec)/(p(p−2))] C(p,m). "
            "Fail drop μ₃; fail n_spec=(3+σ)/2 on collinear."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        Y = load_Y(p, q)
        G3 = (Y * Y[:, 0:1]).T @ Y
        add, neg, chi = F["add"], F["neg"], F["chi"]
        Ns = [field_norm(F, x) for x in range(q)]
        free_vals = []
        max_free = 0.0
        # Φ_free constant on Norm keys
        buckets: dict[tuple[int, int, int], list[float]] = {}
        for u in range(1, q):
            for v in range(1, q):
                if u == v:
                    continue
                w = add[v][neg[u]]
                if w == 0:
                    continue
                col = is_collinear(F, u, v)
                named = g3_1d_named(
                    p, int(chi[u]), int(chi[v]), int(chi[w]), col
                )
                tot = float(np.real(G3[u, v]))
                fr = tot - named
                key = (Ns[u], Ns[v], Ns[w])
                buckets.setdefault(key, []).append(fr)
                free_vals.append(fr)
                max_free = max(max_free, abs(fr))
        n_split = 0
        for vals in buckets.values():
            arr = np.array(vals)
            if arr.max() - arr.min() > 1e-6:
                n_split += 1
        uniq = sorted(set(np.round(free_vals, 6)))
        row_ok = n_split == 0 and max_free > 1.0
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_keys": len(buckets),
            "n_split": n_split,
            "max_abs_free": float(max_free),
            "nuniq_free": len(uniq),
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Φ=G3_1d+Φ_free. Fail Φ_free≡0.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Phi_free_named_in_p": False,
        "Phi_named_in_p": False,
        "Q3_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "G3_1d is named. Φ_free takes 5 values at p=5 and 14 at p=7. "
            "16pA / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.556  G3_1d named by Johnson + square-line count", flush=True)
    A = prove_A()
    print(f"  A G3_1d: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Φ=G3_1d+Φ_free: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.556",
        "title": "G3_1d is named; Φ=G3_1d+Φ_free with Φ_free unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "G3_1d_named": A["proved"],
            "Phi_split": B["proved"],
            "Phi_named_in_p": False,
            "NUM_SUM_16pA_general": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "G3_1d=[n_spec/p−3(m−n_spec)/(p(p−2))] C(p,m), "
            "n_spec=(1+χ(u))/2 if v/u∈F_p else (3+χ(u)+χ(v)+χ(v−u))/2"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
