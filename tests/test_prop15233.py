"""Tests for Prop 15.233 — signed k=2 per2 form; bilinear; residual open."""
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
from e1_gmin_m4_prop15233 import (  # noqa: E402
    hinge_status_233,
    k2_permanent_closed,
    k2_unsigned_majorant,
    residual_i_closed_via_233,
    theorem_k2_bilinear,
    theorem_k2_permanent_closed_form,
    theorem_k2_surviving_14,
    theorem_k2_unsigned_still_dead,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_k2_theorems_proved():
    assert theorem_k2_surviving_14()["proved"] is True
    assert theorem_k2_surviving_14()["n_surviving_perms"] == 14
    assert theorem_k2_permanent_closed_form()["proved"] is True
    assert theorem_k2_bilinear()["proved"] is True
    assert theorem_k2_bilinear()["sufficient_for_residual_i"] is False
    assert theorem_k2_unsigned_still_dead()["proved"] is True
    assert theorem_k2_unsigned_still_dead()["sufficient_for_residual_i"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_233() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status_open():
    h = hinge_status_233()
    assert h["k2_permanent_closed_form"] is True
    assert h["k2_bilinear"] is True
    assert h["residual_i_closed_via_233"] is False
    assert h["e1"] is False


def test_k2_majorant_matches_count_and_exceeds_B():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 89, 101):
        assert is_prime(p)
        assert k2_unsigned_majorant(p) == 14 * layer_count(p * p + 1, 2)
        assert Fraction(k2_unsigned_majorant(p)) > R_bar_budget_for_mu_thr(p)


def test_k2_formula_matches_per_paley_p3():
    C = paley_conference_prime_power(3).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        a, b = S[0], S[1]
        outside = [x for x in range(n) if x not in S]
        for s, t in combinations(outside, 2):
            T = tuple(sorted((a, b, s, t)))
            got = k2_permanent_closed(C, S, a, b, s, t)
            assert got == permanent_submatrix(C, S, T)
            assert abs(got) <= 14
            checked += 1
    assert checked > 0


def test_k2_formula_matches_per_paley_p5_sample():
    C = paley_conference_prime_power(5).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        a, b = S[0], S[1]
        outside = [x for x in range(n) if x not in S]
        for s, t in combinations(outside[:8], 2):
            T = tuple(sorted((a, b, s, t)))
            got = k2_permanent_closed(C, S, a, b, s, t)
            assert got == permanent_submatrix(C, S, T)
            assert abs(got) <= 14
            checked += 1
        if checked >= 80:
            break
    assert checked >= 80
