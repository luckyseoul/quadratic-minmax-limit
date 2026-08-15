"""Tests for Prop 15.222 — ⟨D,TD⟩ closed form; CS δ bound; residual-i open."""
from __future__ import annotations

import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15217 import delta_room_for_R  # noqa: E402
from e1_gmin_m4_prop15221 import D_energy_closed  # noqa: E402
from e1_gmin_m4_prop15222 import (  # noqa: E402
    ADA_energy,
    T1_closed,
    T2_closed,
    T4_closed,
    dtd_closed,
    dtd_closed_via_n,
    hinge_status_222,
    jensen_m0_lower,
    mean_T_on_D,
    residual_i_closed_via_222,
    theorem_CS_delta_le_k2_over_N,
    theorem_dtd_closed_form,
    theorem_m0_k2_bound_open,
)


def test_dtd_closed_proved():
    assert theorem_dtd_closed_form()["proved"] is True


def test_CS_delta_proved():
    assert theorem_CS_delta_le_k2_over_N()["proved"] is True


def test_m0_bound_still_open():
    assert theorem_m0_k2_bound_open()["proved"] is False


def test_dtd_T_pieces_sum():
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        assert is_prime(p)
        assert dtd_closed(p) == dtd_closed_via_n(p)
        assert dtd_closed(p) == T1_closed(p) + 2 * T2_closed(p) + T4_closed(p)
        assert dtd_closed(p) > 0


def test_dtd_values():
    # Verified by full T-matvec census at p=3,5
    assert dtd_closed(3) == 10560
    assert dtd_closed(5) == 2433600
    assert dtd_closed(7) == 63907200


def test_dtd_over_D2_rational():
    assert mean_T_on_D(3) == Fraction(11, 6)
    assert mean_T_on_D(5) == Fraction(13, 2)
    assert mean_T_on_D(7) == Fraction(317, 28)


def test_ADA_positive():
    for p in (3, 5, 7, 11, 13):
        assert ADA_energy(p) > 0
        assert ADA_energy(p) == 4 * p * D_energy_closed(p) - dtd_closed(p)


def test_jensen_m0_insufficient_at_p5():
    """Jensen lower bound on m0 alone does not reach the CS room threshold."""
    n = n_of(5)
    need = Fraction(comb(n, 4)) - delta_room_for_R(5) * 260
    assert jensen_m0_lower(5) < need
    assert theorem_m0_k2_bound_open()["jensen_insufficient_p5"]["jensen_lt_need"]


def test_certified_CS_numerics_in_open_theorem():
    c = theorem_m0_k2_bound_open()
    assert c["certified_small_p"]["p5"]["CS_closes_numerically"] is True
    assert c["certified_small_p"]["p7"]["CS_closes_numerically"] is True


def test_dtd_matches_census_p3_p5():
    """Full sparse T: ⟨D,TD⟩ equals closed form on one Max+ vector."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power
    from scipy.sparse import csr_matrix

    for p in (3, 5):
        C = paley_conference_prime_power(p).astype(float)
        y = boolean_evecs(C, float(p), 0)[0].astype(float)
        n = C.shape[0]
        Ss = list(combinations(range(n), 4))
        N = len(Ss)
        S_idx = {S: i for i, S in enumerate(Ss)}
        rows, cols, vals = [], [], []
        for i, S in enumerate(Ss):
            sset = set(S)
            for v in S:
                for r in range(n):
                    if r in sset:
                        continue
                    Sp = tuple(sorted((sset - {v}) | {r}))
                    rows.append(i)
                    cols.append(S_idx[Sp])
                    vals.append(C[v, r])
        T = csr_matrix((vals, (rows, cols)), shape=(N, N))
        D = np.empty(N)
        for i, S in enumerate(Ss):
            a, b, c, d = S
            prod = y[a] * y[b] * y[c] * y[d]
            q = 2 * (
                C[a, b] * y[a] * y[b]
                + C[a, c] * y[a] * y[c]
                + C[a, d] * y[a] * y[d]
                + C[b, c] * y[b] * y[c]
                + C[b, d] * y[b] * y[d]
                + C[c, d] * y[c] * y[d]
            )
            D[i] = prod * q
        dtd = float(D @ (T @ D))
        assert abs(dtd - float(dtd_closed(p))) < 1e-4
        assert abs(float(D @ D) - float(D_energy_closed(p))) < 1e-4


def test_residual_i_still_open():
    assert residual_i_closed_via_222() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_222()
    assert h["dtd_closed"] is True
    assert h["CS_delta"] is True
    assert h["m0_k2_bound"] is False
    assert h["residual_i_closed_via_222"] is False
