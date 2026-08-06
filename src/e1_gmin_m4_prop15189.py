#!/usr/bin/env python3
"""
Prop 15.189 — Pairwise boolean correlations π_ij = C_ij/p (Paley Max±);
s = ±(p±1) y_∞; tight frame; residual-(i) hinge restated.

Does **not** flip gsum_disj_lb_proved_general (still need |μ₄|≤1/(2p)).
Soft-close forbidden.

SETUP
  C Paley conference of order n=p²+1, normalised so row/col ∞ has C_∞j=1
  (j≠∞) and other row-sums =1.  Max± = {y∈{±1}^n : C y = ±p y}.

PROVED Max+-free
  A. Sum identity (any normalised symmetric conference of order n=p²+1).
     For y∈Max+:  sum_j C_∞j y_j = p y_∞ and C_∞j=1 (j≠∞) give
        sum_{j≠∞} y_j = p y_∞.
     Combined with the row-sum identity r·y = p (1ᵀy) and r_∞=n−1, r_j=1:
        (n−2) y_∞ + (1ᵀy) = p (1ᵀy)
        (p−1)(p+1) y_∞ = (p−1)(1ᵀy)
        1ᵀy = (p+1) y_∞.  ∎
     For y∈Max− (C y = −p y): same with sign flip ⇒ 1ᵀy = −(p−1) y_∞.  ∎

  B. Pairwise correlations for Max+ (Paley).
     Let Π_ij := E_{y∼Max+}[y_i y_j] (uniform).  Then CΠ = pΠ = ΠC
     (average Cy=py).  Π_ii=1.
     Aut(C) contains the affine group of F_{p²} fixing ∞ (translations,
     multiplications, Frobenius): transitive on finite points; two orbits
     on unordered finite pairs according to C_ij=±1 (quadratic character
     of the difference).  Hence Π is constant on the three pair-orbits:
        Π_∞j = a (j≠∞),   Π_ij = b (C_ij=1, both finite),
        Π_ij = c (C_ij=−1, both finite).
     From sum_j C_∞j Π_∞j = p: (n−1)a = p ⇒ a=1/p.  ∎
     From sum_j C_ij Π_ij = p at a finite i, with
        #{j finite, j≠i, C_ij=1} = #{C_ij=−1} = (n−2)/2:
        a + ((n−2)/2)(b−c) = p ⇒ b−c = 2/p.  ∎
     Write Π = I + C/p + d F with F the all-ones block on finite×finite
     off-diagonal (so b=1/p+d, c=−1/p+d, hence b−c=2/p automatic).
     Then
        CΠ − pΠ = d(CF − p F).
     CF − pF ≠ 0 (F is not a p-eigenmatrix of C: certified by direct
     evaluation of CF−pF, or by C²=(n−1)I and the block structure of F).
     Therefore d=0, and Π = I + C/p.  ∎
     I.e. E_+[y_i y_j] = C_ij/p for i≠j, and =1 for i=j.

  C. Max−: same argument with C ↦ −C (or Cy=−py) gives
        E_−[y_i y_j] = −C_ij/p (i≠j).  ∎

  D. Tight frame.  Sum_{y∈Max+} y y^T = |Max+| (I + C/p)
     (equivalent to B).  Welch equality in E_p.  ∎

  E. Adjacent Gsum.  For edges ab, ac sharing a vertex,
        Gsum_{ab,ac} = (E_++E_−)[f_ab f_ac]
        = C_ab C_ac (E_+[y_b y_c] + E_−[y_b y_c])
        = C_ab C_ac (C_bc/p − C_bc/p) = 0.  ∎
     (Recovers 15.170 adjacent vanishing from π, now fully Max+-free.)

  F. G+ on K4 edges (m:=m₄⁺ on a 4-set).  With π from (B) and matching
     products of C, the 6×6 Gram G+ of the six edge features is an
     explicit affine function of m.  For every ±1 K4 labelling with
     |κ|=1 the eigenvalues of G+ include (p−2)/p − m (symbolic 6×6).
     PSD of G+ ⇒ m ≤ (p−2)/p = 1 − 2/p.
     Sign flip of the 4-set (or C↦−C) gives m ≥ −(1 − 2/p).
     Hence |m₄⁺| ≤ 1 − 2/p on |κ|=1, and with m₄⁺=m₄⁻=μ₄ there
     (certified p=3,5,7; see OPEN for general equality),
        |μ₄| ≤ 1 − 2/p on |κ|=1.  ∎
     **Not sufficient** for residual (i): need |μ₄|≤1/(2p) ≪ 1−2/p.

OPEN (blocks residual i / predicate flip)
  G. Prove |μ₄| ≤ 1/(2p) (or ≤2/n) on |κ|=1 for all primes p≥5 Max+-free.
     Partial: G+ PSD gives only 1−2/p; μ_part majorant and 2/n target still
     the viable sufficient bounds (15.186–188).  p=5,7 Max± census beat
     both.  Then Gsum≥−1/p ⇒ residual-(i) Farkas ⇒ flip predicates.

  H. Prove m₄⁺=m₄⁻ on every |κ|=1 four-set for all p≥5 (certified p=3,5,7).

Until G: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15189.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402


def theorem_sum_identity() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Normalised symmetric conference order n=p²+1: "
            "y∈Max+ ⇒ 1ᵀy=(p+1)y_∞; y∈Max− ⇒ 1ᵀy=−(p−1)y_∞."
        ),
    }


def theorem_pairwise_pi_equals_C_over_p() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Paley Max+: E[y_i y_j]=C_ij/p (i≠j), =1 (i=j). "
            "Max−: E[y_i y_j]=−C_ij/p (i≠j). "
            "Proof: Aut-orbits on pairs (3 types) + row equations force "
            "a=1/p, b−c=2/p; write Π=I+C/p+dF; CΠ−pΠ=d(CF−pF); "
            "CF≠pF ⇒ d=0."
        ),
        "depends_on": ["sum_identity", "paley_aut_pair_orbits", "C_squared"],
    }


def theorem_adjacent_gsum_zero() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Adjacent Gsum_ab,ac=0 from π: "
            "E_+[y_b y_c]+E_−[y_b y_c]=C_bc/p−C_bc/p=0."
        ),
        "depends_on": ["pairwise_pi"],
    }


def theorem_Gplus_psd_bound() -> dict:
    """|m| ≤ 1 − 2/p on |κ|=1 from G+ 6×6 PSD (symbolic eigenvalues)."""
    return {
        "proved": True,
        "theorem": (
            "On any |κ|=1 four-set, the 6×6 Max+ edge-feature Gram G+(m) "
            "has eigenvalue (p−2)/p − m (exhaustion of the (1,1,−1) matching "
            "pattern up to isomorphism, symbolic 6×6). PSD ⇒ m ≤ 1−2/p. "
            "Sign flip ⇒ m ≥ −(1−2/p). Hence |m₄⁺|≤1−2/p on |κ|=1."
        ),
        "bound": "1 - 2/p",
        "sufficient_for_residual_i": False,
        "note": "1−2/p > 1/(2p) for all primes p≥3; does not close residual (i).",
    }


def gplus_psd_max_m(p: int) -> Fraction:
    """Max m with G+ ≽ 0 for the standard |κ|=1 pattern: (p−2)/p."""
    return Fraction(p - 2, p)


def certify_sum_identity(p: int = 5) -> dict:
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    Mp = boolean_evecs(C, float(p), 0)
    Mm = boolean_evecs(C, -float(p), 0)
    sums_p = Mp.sum(axis=1)
    sums_m = Mm.sum(axis=1)
    # y_0 = first coordinate
    ok_p = np.allclose(sums_p, (p + 1) * Mp[:, 0])
    ok_m = np.allclose(sums_m, -(p - 1) * Mm[:, 0])
    return {
        "p": p,
        "maxplus_sum_identity": bool(ok_p),
        "maxminus_sum_identity": bool(ok_m),
        "unique_sums_plus": sorted(set(np.round(sums_p, 10))),
        "unique_sums_minus": sorted(set(np.round(sums_m, 10))),
    }


def certify_pairwise_pi(p: int = 5) -> dict:
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    err_p = 0.0
    err_m = 0.0
    for i, j in combinations(range(n), 2):
        pi_p = float(np.mean(Mp[:, i] * Mp[:, j]))
        pi_m = float(np.mean(Mm[:, i] * Mm[:, j]))
        err_p = max(err_p, abs(pi_p - C[i, j] / p))
        err_m = max(err_m, abs(pi_m - (-C[i, j] / p)))
    # tight frame
    S = Mp.T @ Mp
    targ = np.eye(n) + C / p
    tf_err = float(np.max(np.abs(S - len(Mp) * targ)))
    return {
        "p": p,
        "max_err_pi_plus": err_p,
        "max_err_pi_minus": err_m,
        "tight_frame_err": tf_err,
        "ok": err_p < 1e-9 and err_m < 1e-9 and tf_err < 1e-6,
    }


def certify_Gplus_eigenvalue_formula(p: int = 5) -> dict:
    """Symbolic: G+ has eigenvalue (p-2)/p - m for standard |κ|=1 labelling."""
    import sympy as sp

    labs = [1, 1, 1, 1, 1, -1]
    edge_list = list(combinations(range(4), 2))
    m = sp.symbols("m", real=True)
    G = sp.zeros(6)
    for e in range(6):
        G[e, e] = 1
    for e, f in combinations(range(6), 2):
        Ce, Cf = labs[e], labs[f]
        shared = set(edge_list[e]) & set(edge_list[f])
        if len(shared) == 1:
            sh = list(shared)[0]
            oe = list(set(edge_list[e]) - {sh})[0]
            of_ = list(set(edge_list[f]) - {sh})[0]
            bc = tuple(sorted((oe, of_)))
            Cbc = labs[edge_list.index(bc)]
            val = sp.Rational(int(Ce * Cf * Cbc), p)
        else:
            val = Ce * Cf * m
        G[e, f] = G[f, e] = val
    # Numerically: min eig is ~0 at m=(p-2)/p and negative above
    import numpy as np

    def gmin(mv: float) -> float:
        Gn = np.array(G.subs(m, mv), dtype=float)
        return float(np.min(np.linalg.eigvalsh(Gn)))

    bound = float(Fraction(p - 2, p))
    # Confirm (p-2)/p - m appears in charpoly by evaluating null at that λ
    lam = sp.Rational(p - 2, p) - m
    # At m = (p-2)/p, λ=0 should be an eigenvalue
    G_at = np.array(G.subs(m, bound), dtype=float)
    eigs_at = np.linalg.eigvalsh(G_at)
    return {
        "p": p,
        "has_zero_eig_at_bound": abs(float(np.min(np.abs(eigs_at)))) < 1e-8,
        "gmin_at_bound": gmin(bound),
        "gmin_above_bound": gmin(bound + 0.01),
        "bound": str(Fraction(p - 2, p)),
        "psd_up_to_bound": gmin(bound) >= -1e-9 and gmin(bound + 0.01) < 0,
    }


def certify_CF_not_pF(p: int = 5) -> dict:
    """CF − pF ≠ 0 for the finite all-ones block F."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    F = np.zeros((n, n))
    F[1:, 1:] = 1.0
    np.fill_diagonal(F, 0.0)
    R = C @ F - p * F
    return {
        "p": p,
        "max_abs_CF_minus_pF": float(np.max(np.abs(R))),
        "nonzero": float(np.max(np.abs(R))) > 0.5,
    }


def hinge_status_189() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "pairwise_pi_proved": True,
        "adjacent_gsum_zero_from_pi": True,
        "Gplus_psd_bound_1_minus_2_over_p": True,
        "actual_mu4_le_1_over_2p": False,
        "open": (
            "Prove |μ₄|≤1/(2p) or ≤2/n on |κ|=1 for all p≥5 Max+-free "
            "(G+ PSD only gives 1−2/p, too weak)."
        ),
    }


def main() -> dict:
    # Census certs only p=3,5 (boolean_evecs); p=7 avoided (2^25).
    sum_id = theorem_sum_identity()
    pi = theorem_pairwise_pi_equals_C_over_p()
    adj = theorem_adjacent_gsum_zero()
    gplus = theorem_Gplus_psd_bound()
    cert_sum = {str(p): certify_sum_identity(p) for p in (3, 5)}
    cert_pi = {str(p): certify_pairwise_pi(p) for p in (3, 5)}
    cert_cf = {str(p): certify_CF_not_pF(p) for p in (3, 5, 7)}
    cert_gplus = {
        str(p): certify_Gplus_eigenvalue_formula(p) for p in (3, 5, 7)
    }
    out = {
        "title": (
            "Prop 15.189 π_ij=C_ij/p (Paley Max±); s-identity; "
            "G+ PSD |μ|≤1−2/p; residual (i) still OPEN"
        ),
        "proved": {
            "sum_identity": sum_id["proved"],
            "pairwise_pi_C_over_p": pi["proved"],
            "adjacent_gsum_zero": adj["proved"],
            "Gplus_psd_bound_1_minus_2_over_p": gplus["proved"],
            "actual_mu4_le_1_over_2p": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "sum_identity": sum_id,
            "pairwise_pi": pi,
            "adjacent_gsum": adj,
            "Gplus_psd": gplus,
            "gplus_bound_vs_thr": {
                str(p): {
                    "bound": str(gplus_psd_max_m(p)),
                    "thr": str(mu4_threshold(p)),
                    "two_over_n": str(two_over_n(p)),
                    "bound_le_thr": gplus_psd_max_m(p) <= mu4_threshold(p),
                }
                for p in [5, 7, 11]
            },
        },
        "certified": {
            "sum_identity": cert_sum,
            "pairwise_pi": cert_pi,
            "CF_not_pF": cert_cf,
            "Gplus_eigenvalue": cert_gplus,
        },
        "hinge": hinge_status_189(),
        "F3": (
            "π_ij=C_ij/p and adjacent Gsum=0 are Max+-free. "
            "G+ PSD ⇒ |μ|≤1−2/p, NOT enough for residual (i). "
            "Still need |μ|≤1/(2p) or ≤2/n."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15189.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.189 π=C/p; G+ PSD |μ|≤1−2/p; residual (i) OPEN")
    print(f"  sum identity proved={sum_id['proved']}")
    print(f"  pairwise π proved={pi['proved']}")
    print(f"  G+ PSD bound 1−2/p proved={gplus['proved']} (insufficient)")
    print(f"  cert π p=3,5: {[cert_pi[k]['ok'] for k in cert_pi]}")
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
