"""Tests for Prop 15.232 — k=3 per closed form, C² pairing, unsigned layers."""
from __future__ import annotations

import sys
from fractions import Fraction
from itertools import combinations
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
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15231 import permanent_submatrix  # noqa: E402
from e1_gmin_m4_prop15232 import (  # noqa: E402
    c2_pairing_off_S,
    c2_pairing_off_S_direct,
    hinge_status_232,
    k3_layer_R_safe_via_unsigned,
    k3_permanent_closed,
    k3_slack_poly,
    k3_unsigned_majorant,
    layer_count,
    residual_i_closed_via_232,
    theorem_c2_pairing,
    theorem_intersection_decomposition,
    theorem_k3_column_linear,
    theorem_k3_layer_R_safe_large_p,
    theorem_k3_permanent_closed_form,
    theorem_unsigned_layers_012_dead,
    ultra_crude_layer,
    ultra_crude_layer_closed,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def _kappa(C, S):
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def test_intersection_decomposition_proved():
    assert theorem_intersection_decomposition()["proved"] is True


def test_k3_permanent_closed_form_proved():
    t = theorem_k3_permanent_closed_form()
    assert t["proved"] is True
    assert t["n_surviving_perms"] == 11
    assert t["abs_per_ub"] == 11


def test_k3_column_linear_proved():
    assert theorem_k3_column_linear()["proved"] is True
    assert theorem_k3_column_linear()["sufficient_for_residual_i"] is False


def test_c2_pairing_proved():
    assert theorem_c2_pairing()["proved"] is True
    assert theorem_c2_pairing()["sufficient_for_residual_i"] is False


def test_unsigned_layers_012_dead_proved():
    t = theorem_unsigned_layers_012_dead()
    assert t["proved"] is True
    assert t["sufficient_for_residual_i"] is False


def test_k3_layer_R_safe_large_p_proved_but_not_a_close():
    t = theorem_k3_layer_R_safe_large_p()
    assert t["proved"] is True
    assert t["sufficient_for_residual_i"] is False
    assert t["f_89"] == 666520


def test_residual_i_still_open():
    assert residual_i_closed_via_232() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status_open():
    h = hinge_status_232()
    assert h["k3_permanent_closed_form"] is True
    assert h["c2_pairing"] is True
    assert h["unsigned_layers_012_dead"] is True
    assert h["k3_layer_R_safe_large_p"] is True
    assert h["residual_i_closed_via_232"] is False
    assert h["residual_i_dual_eq_empty"] is False
    assert h["gsum_disj_lb"] is False
    assert h["type_I"] is False
    assert h["e1"] is False


def test_layer_counts_partition():
    from math import comb

    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        assert sum(layer_count(n, k) for k in range(4)) == comb(n, 4) - 1
        assert layer_count(n, 3) == 4 * (n - 4)


def test_ultra_crude_closed_matches_count():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        for k in (0, 1, 2, 3):
            assert ultra_crude_layer_closed(p, k) == ultra_crude_layer(p, k)


def test_unsigned_012_exceeds_B_independent():
    """Drive shipped ultra_crude_layer against the shipped B, not a re-impl."""
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 89, 101):
        assert is_prime(p)
        B = R_bar_budget_for_mu_thr(p)
        for k in (0, 1, 2):
            assert Fraction(ultra_crude_layer(p, k)) > B


def test_k3_slack_identity():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 89, 97):
        slack = R_bar_budget_for_mu_thr(p) - Fraction(k3_unsigned_majorant(p))
        assert slack == Fraction(k3_slack_poly(p), 2 * p)


def test_k3_R_safe_threshold():
    for p in range(5, 89):
        if is_prime(p):
            assert k3_layer_R_safe_via_unsigned(p) is False
    for p in (89, 97, 101, 103, 107, 109, 113):
        assert k3_layer_R_safe_via_unsigned(p) is True


def test_k3_permanent_matches_per_on_paley_p3():
    """Shipped closed form vs shipped permanent_submatrix on real Paley C."""
    C = paley_conference_prime_power(3).astype(int)
    n = C.shape[0]
    fours = list(combinations(range(n), 4))
    checked = 0
    for S in fours:
        if abs(_kappa(C, S)) != 1:
            continue
        for d in S:
            for t in range(n):
                if t in S:
                    continue
                T = tuple(sorted([x for x in S if x != d] + [t]))
                got = k3_permanent_closed(C, S, d, t)
                exp = permanent_submatrix(C, S, T)
                assert got == exp
                assert abs(got) <= 11
                checked += 1
    assert checked == 180 * 4 * 6  # |κ|=1 count 180, 4 choices of d, n-4=6


def test_k3_permanent_matches_per_on_paley_p5_sample():
    C = paley_conference_prime_power(5).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        if abs(_kappa(C, S)) != 1:
            continue
        d = S[0]
        for t in range(n):
            if t in S:
                continue
            T = tuple(sorted([x for x in S if x != d] + [t]))
            got = k3_permanent_closed(C, S, d, t)
            assert got == permanent_submatrix(C, S, T)
            assert abs(got) <= 11
            checked += 1
        if checked >= 88:
            break
    assert checked >= 88


def test_c2_pairing_matches_direct_paley():
    for p in (3, 5):
        C = paley_conference_prime_power(p).astype(int)
        n = C.shape[0]
        S = (0, 1, 2, 3)
        pairs = [(0, 0), (0, 1), (1, 2), (0, 5 if n > 5 else 4), (S[3], S[3])]
        for u, v in pairs:
            assert c2_pairing_off_S(C, S, u, v) == c2_pairing_off_S_direct(C, S, u, v)
