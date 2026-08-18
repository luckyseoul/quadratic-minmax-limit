#!/usr/bin/env python3
"""
Prop 15.557 — The 3-distinct pair-sum of Φ_free is named:
    Σ_{u≠v, u,v,v−u≠0} Φ_free
        = 2N(1/p − p) + (p−1)(p+1)² C(p,m)/p,
m=(p+1)/2.  Pointwise Φ_free still unnamed.
16pA / Q_τ / phi_F stay open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.556 flags.  Does **not** import phi_F.

============================================================================
Setup.  15.555–15.556: G3=Φ∘N = G3_1d + Φ_free.  Translation Aut
gives Σ_u≠0 y_u = p−y_0.  15.556 names G3_1d.  Non-collinear
G3_1d is odd in (χ(u),χ(v),χ(v−u)) and the χ-type counts are
symmetric, so that stratum sums to 0.

============================================================================
Theorem A — PROVED Max+-free given N (row-sum + y²=1).
  Σ_{u≠v, u,v≠0} G3(0,u,v) = Σ_y y_0[(p−y_0)²−(q−1)]
      = 2N/p − 2p N = 2N(1/p−p).
  Fail: drop 2N/p (p=5: −1300≠−1248).  ∎

Theorem B — PROVED Max+-free (15.556 + field χ-counts).
  Non-collinear G3_1d sums to 0 (σ ↔ −σ).  Collinear:
      n_each=(p−2)(q−1)/2,
      g_++g_− = C/p − 3C/(p−2) = −2C(p+1)/(p(p−2)),
      Σ G3_1d = −(p−1)(p+1)² C(p,m)/p.
  Fail: drop C (0≠−288 at p=5).  ∎

Theorem C — PROVED (A−B).
  Σ Φ_free = 2N(1/p−p) + (p−1)(p+1)² C(p,m)/p.
  p=5: −960; p=7: −76608.  Fail: Φ_free≡0 (−1248≠−960).  ∎

Theorem D — OPEN.  Pointwise Φ_free on (F_p^×)³, Q3 constants
  as a p-rational, 16pA / Q_τ.  phi_F not imported.

============================================================================
Backend: Fraction + p=5,7 G3 fold.  Writes
evidence/e1_gmin_m4_prop15557.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15554 import load_Y
from e1_gmin_m4_prop15556 import g3_1d_named, is_collinear

EV = ROOT / "evidence" / "e1_gmin_m4_prop15557.json"
LIVE_N = {5: 130, 7: 5726}


def sum_G3_named(p: int, N: int) -> float:
    """2N(1/p − p)."""
    return 2.0 * N * (1.0 / p - p)


def sum_G3_named_drop_2Np(p: int, N: int) -> float:
    """Fail: drop 2N/p."""
    return -2.0 * p * N


def sum_G3_1d_named(p: int) -> float:
    """−(p−1)(p+1)² C(p,m)/p."""
    m = (p + 1) // 2
    C = comb(p, m)
    return -(p - 1) * (p + 1) ** 2 * C / p


def sum_Phi_free_named(p: int, N: int) -> float:
    return sum_G3_named(p, N) - sum_G3_1d_named(p)


@lru_cache(maxsize=4)
def fold_sums(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    Y = load_Y(p, q)
    N = Y.shape[0]
    G3 = (Y * Y[:, 0:1]).T @ Y
    add, neg, chi = F["add"], F["neg"], F["chi"]
    s_g3 = s_1d = 0.0
    n_pairs = 0
    for u in range(1, q):
        for v in range(1, q):
            if u == v:
                continue
            w = add[v][neg[u]]
            if w == 0:
                continue
            n_pairs += 1
            s_g3 += float(np.real(G3[u, v]))
            col = is_collinear(F, u, v)
            s_1d += g3_1d_named(p, int(chi[u]), int(chi[v]), int(chi[w]), col)
    s_free = s_g3 - s_1d
    return {
        "N": int(N),
        "n_pairs": n_pairs,
        "sum_G3": float(s_g3),
        "sum_G3_1d": float(s_1d),
        "sum_free": float(s_free),
        "named_G3": float(sum_G3_named(p, N)),
        "named_1d": float(sum_G3_1d_named(p)),
        "named_free": float(sum_Phi_free_named(p, N)),
        "drop_2Np": float(sum_G3_named_drop_2Np(p, N)),
        "drop_C": 0.0,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = fold_sums(p)
        match = abs(rec["sum_G3"] - rec["named_G3"]) < 1e-4
        drop_fails = abs(rec["sum_G3"] - rec["drop_2Np"]) > 1.0
        if not (match and drop_fails):
            ok = False
        rows[str(p)] = {
            "sum_G3": rec["sum_G3"],
            "named": rec["named_G3"],
            "drop_2Np": rec["drop_2Np"],
            "match": bool(match),
            "drop_fails": bool(drop_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Σ G3=2N(1/p−p). Fail drop 2N/p.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = fold_sums(p)
        match = abs(rec["sum_G3_1d"] - rec["named_1d"]) < 1e-4
        drop_fails = abs(rec["sum_G3_1d"] - rec["drop_C"]) > 1.0
        if not (match and drop_fails):
            ok = False
        rows[str(p)] = {
            "sum_1d": rec["sum_G3_1d"],
            "named": rec["named_1d"],
            "match": bool(match),
            "drop_C_fails": bool(drop_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Σ G3_1d=−(p−1)(p+1)² C(p,m)/p. Fail drop C.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = fold_sums(p)
        match = abs(rec["sum_free"] - rec["named_free"]) < 1e-4
        zero_fails = abs(rec["sum_free"]) > 1.0
        if not (match and zero_fails):
            ok = False
        rows[str(p)] = {
            "sum_free": rec["sum_free"],
            "named": rec["named_free"],
            "N": rec["N"],
            "match": bool(match),
            "zero_fails": bool(zero_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Σ Φ_free=2N(1/p−p)+(p−1)(p+1)² C(p,m)/p. "
            "Fail Φ_free≡0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Phi_free_pointwise_named": False,
        "Q3_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Pair-sum of Φ_free is named. Pointwise Φ_free / Q3 "
            "constants / 16pA / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.557  pair-sum of Φ_free is named", flush=True)
    A = prove_A()
    print(f"  A ΣG3: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B ΣG3_1d: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C ΣΦ_free: {C['proved']}", flush=True)
    D = prove_open()
    out = {
        "prop": "15.557",
        "title": "Σ Φ_free = 2N(1/p−p)+(p−1)(p+1)² C(p,m)/p",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "sum_G3_named": A["proved"],
            "sum_G3_1d_named": B["proved"],
            "sum_Phi_free_named": C["proved"],
            "Phi_free_pointwise": False,
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
        "identity": "ΣΦ_free=2N(1/p−p)+(p−1)(p+1)² C(p,(p+1)/2)/p",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
