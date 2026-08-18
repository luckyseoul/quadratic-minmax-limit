#!/usr/bin/env python3
"""
Prop 15.574 — μ_{k=3} is the Fejér 4-point complete sum

    P(r) = Σ_{t∈F_p^×} 1 / (|ω^t−1|² |ω^{rt}−1|²)
    μ_{k=3}(r) = 96 p⁴ P(r) / (p²−1)

for r∈F_p^×\\{±1} (type sub).  15.276: on a support line,
ẑ(c)=F_{λ,s}(c)=2p ω^{−cλ^{-1}s}/(ω^{cλ^{-1}}−1), so
|ẑ(c)|²|ẑ(rc)|²=16p⁴/(|ω^c−1|²|ω^{rc}−1|²).  The line
through ξ sits in C(m−1,2) of the C(m,3) supports, hence
the factor 3/m=6/(p+1) and 16p⁴·3/m·P/(p−1)=96p⁴ P/(p²−1).
Hits live μ_{k=3} at p=5 (2000) and p=7 (9604=4p⁴).
Fail: 4p⁴ at p=5.  Fail: drop 3/m at p=7.  Fail: P(1) at p=7.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.573 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.573.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.573: Q_sub=(n_1d μ_1d+n_{k=3} μ_{k=3}+n_{k>3} μ_full)/N.
15.276 / 15.562: k=3 is locked Fejér on each support line.

============================================================================
Theorem A — PROVED (15.276 Fejér + support count; p=5,7 cache).
  μ_{k=3}(r)=96 p⁴ P(r)/(p²−1).  Independent of r∈F_p^×\\{±1}
  at p=5,7.  Fail 4p⁴ (p=5: 2500≠2000).  Fail drop 3/m
  (p=7: 16p⁴ P/(p−1)=12805.3≠9604).  Fail P↦P(1) (p=7).  ∎

Theorem B — OPEN.  μ_1d and μ_full unnamed as a p-law.
  P(r) splits by r at p≥11 (no H_+ cache).  Q_N* unnamed.
  Do not interpolate (p−5)/15.  phi_F not imported.

============================================================================
Backend: field Fejér O(p) + p=5,7 yhat μ.  Writes
evidence/e1_gmin_m4_prop15574.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15573 import q_sub_mix

EV = ROOT / "evidence" / "e1_gmin_m4_prop15574.json"


def P_fejer(p: int, r: int) -> float:
    omega = np.exp(2j * np.pi / p)
    acc = 0.0
    for t in range(1, p):
        d1 = abs(omega**t - 1.0) ** 2
        d2 = abs(omega ** ((r * t) % p) - 1.0) ** 2
        acc += 1.0 / (d1 * d2)
    return acc


def mu_k3_named(p: int, r: int) -> float:
    return 96.0 * (p**4) * P_fejer(p, r) / (p * p - 1)


def mu_k3_drop_frac(p: int, r: int) -> float:
    """Fail: drop 3/m, i.e. 16 p⁴ P/(p−1)."""
    return 16.0 * (p**4) * P_fejer(p, r) / (p - 1)


def mu_k3_four_p4(p: int) -> float:
    return 4.0 * (p**4)


def mu_k3_P_of_one(p: int) -> float:
    """Fail: replace P(r) by P(1)."""
    return 96.0 * (p**4) * P_fejer(p, 1) / (p * p - 1)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = q_sub_mix(p)["parts"]["k3"]["mu"]
        mus = [mu_k3_named(p, r) for r in range(2, p - 1)]
        match = all(abs(m - live) < 1e-6 for m in mus)
        spread = float(max(mus) - min(mus))
        # 4p⁴ coincides at p=7; it must fail at p=5
        fail4 = abs(live - mu_k3_four_p4(p)) > 1.0 if p != 7 else True
        fail_frac = abs(live - mu_k3_drop_frac(p, 2)) > 1.0 if p == 7 else True
        fail_p1 = abs(live - mu_k3_P_of_one(p)) > 1.0
        if not (match and spread < 1e-8 and fail4 and fail_frac and fail_p1):
            ok = False
        rows[str(p)] = {
            "live": live,
            "named": mus[0],
            "spread": spread,
            "four_p4": mu_k3_four_p4(p),
            "drop_frac": mu_k3_drop_frac(p, 2),
            "P_of_1": mu_k3_P_of_one(p),
            "match": bool(match),
            "indep_r": bool(spread < 1e-8),
            "four_p4_fails": bool(fail4),
            "drop_frac_fails_or_p5": bool(fail_frac),
            "P1_fails": bool(fail_p1),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "μ_{k=3}=96 p⁴ P(r)/(p²−1). Fail 4p⁴ at p=5; "
            "fail drop 3/m at p=7; fail P(1) at p=7."
        ),
    }


def prove_open() -> dict:
    # P(r) splits at p=11
    mus = [mu_k3_named(11, r) for r in range(2, 10)]
    return {
        "proved": False,
        "p11_mu_spread": float(max(mus) - min(mus)),
        "Q_sub_closed_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "μ_{k=3} is the Fejér 4-point. μ_1d / μ_full / Q_N* unnamed. "
            "P(r) splits at p≥11. Do not interpolate (p−5)/15. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.574  μ_k3 = 96 p⁴ P(r)/(p²−1) Fejér 4-point", flush=True)
    A = prove_A()
    print(f"  A Fejér μ_k3: {A['proved']}", flush=True)
    B = prove_open()
    out = {
        "prop": "15.574",
        "title": "μ_{k=3}=96 p⁴ P(r)/(p²−1), P=Σ 1/(|ω^t−1|²|ω^{rt}−1|²)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "proved": {
            "mu_k3_fejer_named": A["proved"],
            "Q_sub_closed_in_p": False,
            "Q_tau_named_in_p": False,
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
        "identity": "mu_k3 = 96 p^4 P(r) / (p^2-1)",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
