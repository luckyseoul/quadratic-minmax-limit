"""Tests for Prop 15.234 — signed k=1 Laplace; trilinear; residual open."""
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
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15231 import permanent_submatrix  # noqa: E402
from e1_gmin_m4_prop15232 import layer_count  # noqa: E402
from e1_gmin_m4_prop15234 import (  # noqa: E402
    hinge_status_234,
    k1_permanent_closed,
    k1_unsigned_majorant,
    residual_i_closed_via_234,
    theorem_k1_laplace,
    theorem_k1_surviving_18,
    theorem_k1_trilinear,
    theorem_k1_unsigned_still_dead,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_k1_theorems_proved():
    assert theorem_k1_surviving_18()["proved"] is True
    assert theorem_k1_surviving_18()["n_surviving_perms"] == 18
    assert theorem_k1_laplace()["proved"] is True
    assert theorem_k1_trilinear()["proved"] is True
    assert theorem_k1_trilinear()["sufficient_for_residual_i"] is False
    assert theorem_k1_unsigned_still_dead()["proved"] is True


def test_residual_i_still_open():
    assert residual_i_closed_via_234() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status_open():
    h = hinge_status_234()
    assert h["k1_laplace"] is True
    assert h["k1_trilinear"] is True
    assert h["residual_i_closed_via_234"] is False
    assert h["e1"] is False


def test_k1_majorant_matches_count_and_exceeds_B():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 89):
        assert is_prime(p)
        assert k1_unsigned_majorant(p) == 18 * layer_count(p * p + 1, 1)
        assert Fraction(k1_unsigned_majorant(p)) > R_bar_budget_for_mu_thr(p)


def test_k1_laplace_matches_per_paley_p3():
    C = paley_conference_prime_power(3).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        a = S[0]
        outside = [x for x in range(n) if x not in S]
        for r, s, t in combinations(outside, 3):
            T = tuple(sorted((a, r, s, t)))
            got = k1_permanent_closed(C, S, a, r, s, t)
            assert got == permanent_submatrix(C, S, T)
            assert abs(got) <= 18
            checked += 1
    assert checked > 0


def test_k1_laplace_matches_per_paley_p5_sample():
    C = paley_conference_prime_power(5).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        a = S[0]
        outside = [x for x in range(n) if x not in S]
        for r, s, t in combinations(outside[:7], 3):
            T = tuple(sorted((a, r, s, t)))
            got = k1_permanent_closed(C, S, a, r, s, t)
            assert got == permanent_submatrix(C, S, T)
            assert abs(got) <= 18
            checked += 1
        if checked >= 60:
            break
    assert checked >= 60
