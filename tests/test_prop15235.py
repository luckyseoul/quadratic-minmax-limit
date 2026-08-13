"""Tests for Prop 15.235 — signed k=0 cycle-type + inverse pairing; residual open."""
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
from e1_gmin_m4_prop15235 import (  # noqa: E402
    cycle_type_permanent_parts,
    four_cycle_inverse_pair_sum,
    hinge_status_235,
    k0_unsigned_majorant,
    residual_i_closed_via_235,
    theorem_k0_cycle_type_split,
    theorem_k0_four_cycle_inverse_pairing,
    theorem_k0_full_s4,
    theorem_k0_unsigned_still_dead,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_k0_theorems_proved():
    assert theorem_k0_full_s4()["proved"] is True
    assert theorem_k0_full_s4()["n_surviving_perms"] == 24
    assert theorem_k0_cycle_type_split()["proved"] is True
    assert theorem_k0_four_cycle_inverse_pairing()["proved"] is True
    assert theorem_k0_unsigned_still_dead()["proved"] is True
    assert theorem_k0_unsigned_still_dead()["sufficient_for_residual_i"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_235() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status_open():
    h = hinge_status_235()
    assert h["k0_cycle_type_split"] is True
    assert h["k0_four_cycle_inverse_pairing"] is True
    assert h["residual_i_closed_via_235"] is False
    assert h["e1"] is False


def test_k0_majorant_matches_count_and_exceeds_B():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 89):
        assert is_prime(p)
        assert k0_unsigned_majorant(p) == 24 * layer_count(p * p + 1, 0)
        assert Fraction(k0_unsigned_majorant(p)) > R_bar_budget_for_mu_thr(p)


def test_cycle_parts_sum_to_permanent_p3():
    C = paley_conference_prime_power(3).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        outside = [x for x in range(n) if x not in S]
        if len(outside) < 4:
            continue
        for T in combinations(outside, 4):
            parts = cycle_type_permanent_parts(C, S, T)
            assert sum(parts.values()) == permanent_submatrix(C, S, T)
            assert abs(sum(parts.values())) <= 24
            assert parts["4cyc"] == four_cycle_inverse_pair_sum(C, S, T)
            assert parts["4cyc"] % 2 == 0
            assert abs(parts["4cyc"]) <= 6
            checked += 1
    assert checked > 0


def test_cycle_parts_sum_to_permanent_p5_sample():
    C = paley_conference_prime_power(5).astype(int)
    n = C.shape[0]
    checked = 0
    for S in combinations(range(n), 4):
        outside = [x for x in range(n) if x not in S]
        for T in combinations(outside[:8], 4):
            parts = cycle_type_permanent_parts(C, S, T)
            per = permanent_submatrix(C, S, T)
            assert sum(parts.values()) == per
            assert parts["4cyc"] == four_cycle_inverse_pair_sum(C, S, T)
            assert parts["4cyc"] % 2 == 0
            assert abs(parts["4cyc"]) <= 6
            checked += 1
        if checked >= 40:
            break
    assert checked >= 40
