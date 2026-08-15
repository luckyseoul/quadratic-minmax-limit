"""Prop 15.251: Cy-identity structure; 2/n suffices; residual_i open."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold
from e1_gmin_m4_prop15251 import (
    hinge_status_251,
    residual_i_closed_via_251,
    theorem_Cy_identity_structure,
    theorem_per_equals_derang,
    theorem_two_over_n_suffices,
    two_over_n_le_half_p,
)


def test_per_equals_derang():
    t = theorem_per_equals_derang()
    assert t["proved"] is True


def test_Cy_identity_structure():
    t = theorem_Cy_identity_structure()
    assert t["proved"] is True


def test_two_over_n_suffices():
    t = theorem_two_over_n_suffices()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert two_over_n_le_half_p(p)
        assert Fraction(2, n_of(p)) <= mu4_threshold(p)


def test_predicates_remain_false():
    assert residual_i_closed_via_251() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
    h = hinge_status_251()
    assert h["residual_i_closed"] is False
    assert h["mu4_bound_general"] is False


def test_census_mu4_if_caches():
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            pytest.skip(f"no Max± cache p={p}")
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(float)
        n = C.shape[0]
        # sample 200 |κ|=1 sets
        count = 0
        max_mu = 0.0
        max_diff = 0.0
        for a, b, c, d in combinations(range(n), 4):
            ka = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
            if abs(abs(ka) - 1) > 1e-9:
                continue
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            mm = float(np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d]))
            max_mu = max(max_mu, abs(0.5 * (mp + mm)))
            max_diff = max(max_diff, abs(mp - mm))
            count += 1
            if count >= 200:
                break
        assert count >= 200
        assert max_diff < 1e-9  # m4+=m4- on sample
        assert max_mu <= 2 / n + 1e-12
