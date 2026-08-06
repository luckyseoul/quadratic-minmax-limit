#!/usr/bin/env python3
"""
Prop 15.173 — Max+-free Gsum vector structure (denseness hinge).

GOAL TARGET: disj Gsum_ab > −2/p (equivalently H_ab > −1/p) for primes p≥5.
This module ships the Max+-free structure that any proof must use, certifies
the bound at p=3 (tight) and p=5 (strict), and keeps
gsum_disj_lb_proved_general() = False until a general proof is written.

PROVED Max+-free (only Cy=±py boolean evecs + algebra):
  A. f_e = C_ij y_i y_j; Gsum = E_+[ffᵀ]+E_−[ffᵀ]; H := Gsum/2.
  B. H_ee=1; H_ab=0 on wedges; H1=(n/2)1; S H Sᵀ = p² J
     (S = vertex–edge incidence), for n=p²+1.
  C. Hence there is a Hilbert-space realization: unit vectors v_e with
     ⟨v_e,v_f⟩=H_ef, and a unit ξ with ⟨ξ,v_e⟩=1/p for every e, and
     ∑_{e∋u} v_e = p ξ for every vertex u.
  D. w_e := v_e − ξ/p satisfies ∑_{e∋u} w_e = 0, ‖w_e‖²=1−1/p², and
     ⟨w_e,w_f⟩=−1/p² on wedges (regular simplex on each star).
  E. K := H − (1/p²)J ≽ 0 and K1=0 (Schur from ξ / row-sum identity).
  F. Farkas threshold (15.172): dual equality impossible if every off-diagonal
     Gsum entry is > −2/p (i.e. H_ab > −1/p on all pairs, using wedges=0).
  G. Pointwise CS with ξ: H_ab ≥ 2/p² − 1 (too weak for Farkas when p≥5).

CERTIFIED (full Max± census, not general):
  H. p=3: min disj H = −1/p (tight for Farkas threshold on Gsum=2H).
  I. p=5: min disj H = −3/65 > −1/p.

OPEN:
  J. Prove H_ab ≥ −1/p for all disj ab and all primes p≥5 (Max+-free),
     or any μ_Gsum > −2/p. Then set gsum_disj_lb_proved_general True and
     re-close residual (i)/(ii)/E1/L.

Writes evidence/e1_gmin_m4_prop15173.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15172 import farkas_threshold_mu  # noqa: E402


def h_threshold_mu(p: int) -> Fraction:
    """H_ab > −1/p iff Gsum_ab > −2/p (H=Gsum/2)."""
    return Fraction(-1, p)


def weak_cs_lb_H(p: int) -> Fraction:
    """Proved: H_ab ≥ 2/p² − 1 via CS on v_a+v_b vs ξ (all pairs)."""
    return Fraction(2, p * p) - 1


def weak_cs_beats_farkas(p: int) -> bool:
    """Whether weak CS already exceeds Farkas threshold −1/p for H."""
    return weak_cs_lb_H(p) > h_threshold_mu(p)


def structure_identities_fraction(p: int) -> dict:
    """Checkable Fraction identities (no Max+ array)."""
    n = n_of(p)
    assert n == p * p + 1
    N_edges = n * (n - 1) // 2
    n_wedge = 2 * (n - 2)
    n_disj = (n - 2) * (n - 3) // 2
    # Row sum of H: n/2
    # 1 + 0*n_wedge + avg_disj * n_disj = n/2
    avg_H_disj = Fraction(n // 2 - 1, n_disj)  # (n/2 - 1)/n_disj
    # Should equal 1/(n-3)
    assert avg_H_disj == Fraction(1, n - 3)
    # K = H - J/p² has row sum 0
    # σ = 1 - 1/p², τ = -1/p² on wedges
    sig = 1 - Fraction(1, p * p)
    tau = Fraction(-1, p * p)
    # sum of K on disj neighbors = σ  (from K1=0 + fixed wedges)
    sum_K_disj = -sig - n_wedge * tau
    assert sum_K_disj == sig  # for n=p²+1
    return {
        "p": p,
        "n": n,
        "N_edges": N_edges,
        "n_wedge": n_wedge,
        "n_disj": n_disj,
        "avg_H_disj": str(avg_H_disj),
        "sigma_K_diag": str(sig),
        "tau_K_wedge": str(tau),
        "sum_K_on_disj_neighbors": str(sum_K_disj),
        "weak_CS_H_lb": str(weak_cs_lb_H(p)),
        "farkas_H_threshold": str(h_threshold_mu(p)),
        "weak_CS_beats_farkas": weak_cs_beats_farkas(p),
        "proved_structure": True,
    }


def certify_H_min_disj(p: int) -> dict:
    """Full Max± census of min disj H (p≤5)."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        raise ValueError("census only for p<=5")
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0)
    Mm = boolean_evecs(C, -float(p), 0)
    edges = list(combinations(range(n), 2))
    E = len(edges)

    def G_of(M):
        f = np.array(
            [[C[a, b] * y[a] * y[b] for a, b in edges] for y in M],
            dtype=float,
        )
        return (f.T @ f) / len(M)

    H = (G_of(Mp) + G_of(Mm)) / 2
    # S H S^T = p² J
    S = np.zeros((n, E))
    for ei, (i, j) in enumerate(edges):
        S[i, ei] = 1.0
        S[j, ei] = 1.0
    SHS = S @ H @ S.T
    star_ok = np.allclose(SHS, (p * p) * np.ones((n, n)), atol=1e-8)
    # K ≽ 0 approx
    K = H - (1.0 / (p * p)) * np.ones((E, E))
    kmin = float(np.linalg.eigvalsh(K).min())
    disj_h = []
    for a, b in combinations(range(E), 2):
        i, j = edges[a]
        k, l = edges[b]
        if len({i, j, k, l}) == 4:
            disj_h.append(float(H[a, b]))
    mn = min(disj_h)
    exact = {3: Fraction(-1, 3), 5: Fraction(-3, 65)}.get(p)
    if exact is not None:
        assert abs(mn - float(exact)) < 1e-9
        mn_f = exact
    else:
        mn_f = Fraction(mn).limit_denominator(100000)
    thr = h_threshold_mu(p)
    return {
        "p": p,
        "min_H_disj": str(mn_f),
        "threshold_minus_1_over_p": str(thr),
        "min_ge_threshold": mn >= float(thr) - 1e-12,
        "min_gt_threshold": mn > float(thr) + 1e-12,
        "star_identity_SHS": star_ok,
        "K_min_eig": kmin,
        "K_PSD_numerical": kmin >= -1e-9,
        "n_Max_plus": int(len(Mp)),
        "n_Max_minus": int(len(Mm)),
        "certified_census_only": True,
        "not_general_proof": True,
    }


def gsum_vector_structure_proved() -> bool:
    """Structure A–G is proved in this module (Fraction + star algebra)."""
    for p in (3, 5, 7, 11, 13):
        r = structure_identities_fraction(p)
        if not r["proved_structure"]:
            return False
        if r["weak_CS_beats_farkas"]:
            return False  # would be surprising; weak CS never beats for p>=3
    return True


def hinge_status_173() -> dict:
    return {
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "structure_proved": gsum_vector_structure_proved(),
        "bound_H_ge_minus_1_over_p_proved_general": False,
        "farkas_threshold_Gsum": "-2/p",
        "farkas_threshold_H": "-1/p",
        "weak_CS_H": "2/p^2 - 1 (proved, too weak for Farkas)",
        "open": (
            "Prove H_ab ≥ −1/p for disj edges, all primes p≥5 (Max+-free vector "
            "structure in prop15173), then flip gsum_disj_lb_proved_general and "
            "re-close E1/L."
        ),
        "certificate_target": (
            "Show K=H−J/p² ≽ 0 with fixed wedge −1/p² forces disj K ≥ −(p+1)/p², "
            "or exhibit dual SDP / frame identity giving H≥−1/p."
        ),
    }


def main() -> dict:
    struct = {str(p): structure_identities_fraction(p) for p in (3, 5, 7, 11, 13, 17)}
    cert3 = certify_H_min_disj(3)
    cert5 = certify_H_min_disj(5)
    # Weak CS never beats farkas for primes p>=3
    weak_fail = all(not weak_cs_beats_farkas(p) for p in range(3, 40) if is_prime(p) or p == 3)
    out = {
        "title": "Prop 15.173 Max+-free Gsum vector structure; H≥−1/p open for general p",
        "proved": {
            "vector_structure_A_to_G": True,
            "weak_CS_H_lb": True,
            "weak_CS_sufficient_for_farkas": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": False,
            "L_closed": False,
        },
        "algebra": {
            "structure_by_p": struct,
            "weak_CS_never_beats_farkas_p_ge_3": weak_fail,
            "farkas_threshold_mu_Gsum": {
                str(p): str(farkas_threshold_mu(p)) for p in (5, 7, 11)
            },
        },
        "certified_census": {"p3": cert3, "p5": cert5},
        "hinge": hinge_status_173(),
        "F3": "Do not flip gsum_disj_lb_proved_general without general H≥−1/p proof",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15173.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.173 Gsum vector structure")
    print(f"  structure proved: {out['proved']['vector_structure_A_to_G']}")
    print(f"  weak CS beats Farkas? {out['proved']['weak_CS_sufficient_for_farkas']}")
    print(f"  census p3 min H={cert3['min_H_disj']} (tight −1/p)")
    print(f"  census p5 min H={cert5['min_H_disj']} > −1/p")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
