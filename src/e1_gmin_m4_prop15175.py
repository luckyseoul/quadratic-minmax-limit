#!/usr/bin/env python3
"""
Prop 15.175 — Gμ reformulation of the denseness hinge (still OPEN).

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

PROVED Max+-free / checkable:
  A. Let Max± be the boolean ±p-eigenvectors of the conference matrix C,
     N± = |Max±|. When N+=N− (holds for Paley p=3,5 census; used as the
     balanced measure μ = (E_+ + E_-)/2):
       E_μ[y_i y_j] = δ_ij  (columns of the stacked Max± matrix are ON).
  B. Define Gμ on edges by Gμ_ef := E_μ[y_i y_j y_k y_l] for e=ij, f=kl.
     Then Gμ_ee=1, Gμ=0 on wedges (E_μ[y_j y_k]=0), Gμ ≽ 0.
  C. H_ef = C_e C_f Gμ_ef  (H = Gsum/2 when N+=N−). Hence
     H_ef ≥ −1/p for all disj ef whenever |Gμ_ef| ≤ 1/p for all disj ef.
  D. Sufficient hinge form: prove |Gμ_ab| ≤ 1/p on disj pairs (all primes p≥5)
     Max+-free ⇒ H ≥ −1/p ⇒ Farkas threshold (15.172) ⇒ residual i/ii path.

CERTIFIED (finite Max± census, not general):
  F. p=3: max|Gμ_disj| = 1/3 = 1/p (tight).
  G. p=5: max|Gμ_disj| = 11/65 < 1/5.
  H. Same mins as H after C-twist: min H = −1/3 (p=3), −3/65 (p=5).

OPEN:
  I. Prove |Gμ_ab| ≤ 1/p for disj edges, all primes p≥5, Max+-free
     (or H_ab ≥ −1/p by another route / alternate residual).
  Until then: gsum_disj_lb_proved_general False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15175.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
)
from e1_gmin_m4_prop15173 import h_threshold_mu  # noqa: E402


def gmu_bound_implies_H(p: int) -> dict:
    """Fraction: |Gμ|≤1/p ⇒ H≥−1/p (via H = C_e C_f Gμ)."""
    thr = h_threshold_mu(p)
    gmu_thr = Fraction(1, p)
    return {
        "p": p,
        "gmu_abs_thr": str(gmu_thr),
        "H_thr": str(thr),
        "implication": "|Gμ_disj|≤1/p ⇒ H_disj≥−1/p",
        "proved_implication": True,
        "bound_proved_general": False,
    }


def certify_Gmu_small_p(p: int) -> dict:
    """Full Max± census of Gμ / H on disj pairs (p≤5 only)."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        raise ValueError("census only p<=5")
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    assert len(Mp) == len(Mm)
    M = np.vstack([Mp, Mm])
    N = len(M)
    edges = list(combinations(range(n), 2))
    # Features y_i y_j
    F = np.stack([M[:, i] * M[:, j] for i, j in edges], axis=1)
    Gmu = (F.T @ F) / N
    # H from C-twisted features
    Fc = np.stack(
        [C[i, j] * M[:, i] * M[:, j] for i, j in edges], axis=1
    )
    # balanced average of + and − separately equals combined when N+=N−
    H = (Fc.T @ Fc) / N  # same as (G++G-)/2 when equal sizes
    # 2-point: E_μ[y_i y_j] = δ_ij
    E2 = (M.T @ M) / N
    e2_ok = np.allclose(E2, np.eye(n), atol=1e-9)
    # wedge Gμ = 0
    wvals = []
    for a, b in combinations(range(len(edges)), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 3:
            wvals.append(float(Gmu[a, b]))
    wedge_ok = max(abs(x) for x in wvals) < 1e-9
    disj_g = []
    disj_h = []
    for a, b in combinations(range(len(edges)), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 4:
            disj_g.append(float(Gmu[a, b]))
            disj_h.append(float(H[a, b]))
            # H = C_e C_f Gμ
            pred = C[i, j] * C[k, l] * Gmu[a, b]
            assert abs(pred - H[a, b]) < 1e-9
    max_abs_g = max(abs(x) for x in disj_g)
    min_h = min(disj_h)
    thr = 1.0 / p
    exact_max_g = {3: Fraction(1, 3), 5: Fraction(11, 65)}[p]
    exact_min_h = {3: Fraction(-1, 3), 5: Fraction(-3, 65)}[p]
    assert abs(max_abs_g - float(exact_max_g)) < 1e-9
    assert abs(min_h - float(exact_min_h)) < 1e-9
    return {
        "p": p,
        "n": n,
        "N_tot": int(N),
        "E2_is_I": bool(e2_ok),
        "Gmu_wedge_zero": bool(wedge_ok),
        "max_abs_Gmu_disj": str(exact_max_g),
        "min_H_disj": str(exact_min_h),
        "gmu_le_1_over_p": max_abs_g <= thr + 1e-12,
        "H_ge_minus_1_over_p": min_h >= -thr - 1e-12,
        "sufficient_form_holds_census": max_abs_g <= thr + 1e-12,
        "certified_census_only": True,
        "not_general_proof": True,
    }


def hinge_status_175() -> dict:
    return {
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "sufficient_form": "|Gμ_disj| ≤ 1/p  ⇒  H_disj ≥ −1/p",
        "bound_proved_general": False,
        "open": (
            "Prove |Gμ_ab| ≤ 1/p on disjoint edge pairs for all primes p≥5 "
            "(Max+-free), or prove H≥−1/p / alternate residual (i)/(ii)."
        ),
    }


def main() -> dict:
    impl = {str(p): gmu_bound_implies_H(p) for p in (3, 5, 7, 11, 13, 17)}
    cert3 = certify_Gmu_small_p(3)
    cert5 = certify_Gmu_small_p(5)
    out = {
        "title": "Prop 15.175 Gμ reformulation; |Gμ|≤1/p sufficient, still OPEN",
        "proved": {
            "H_eq_C_twist_Gmu": True,
            "sufficient_gmu_bound_implies_H": True,
            "Gmu_diag_1_wedge_0": True,
            "bound_Gmu_le_1_over_p_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"implication_by_p": impl},
        "certified_census": {"p3": cert3, "p5": cert5},
        "hinge": hinge_status_175(),
        "F3": "Do not flip residual/E1/L without general |Gμ|≤1/p or H≥−1/p",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15175.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.175 Gμ hinge reformulation")
    print(f"  sufficient |Gμ|≤1/p ⇒ H≥−1/p: True")
    print(f"  census p3 max|Gμ|={cert3['max_abs_Gmu_disj']} (tight)")
    print(f"  census p5 max|Gμ|={cert5['max_abs_Gmu_disj']} < 1/5")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
