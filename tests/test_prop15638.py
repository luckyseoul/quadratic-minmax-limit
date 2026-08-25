from itertools import product

import pytest

from scripts.r1_next_shell_half_conic_audit import audit_prime
from src.e1_gmin_m4_prop15638 import (
    candidate_common_sum_magnitudes,
    candidate_scaled_norm,
    candidate_shell_excluded,
    candidate_shell_gap_theorem,
    excess_over_balancing,
    sum_one_doubled_energy5_defect,
    sum_one_doubled_energy5_factor,
    sum_one_energy3_defect,
    sum_one_energy3_factor,
    t2_energy4_nondegenerate_cubic_defect,
    t2_energy4_nondegenerate_cubic_factor,
    t2_energy6_doubled_cubic_defect,
    t2_energy6_doubled_cubic_factor,
    t2_hasse_obstruction,
    unsigned_pair_cubic_defect,
    unsigned_pair_quartic_defect,
)


def coefficient_patterns(total: int, energy: int) -> set[tuple[int, ...]]:
    """Exhaust local patterns through energy six, modulo zeroes/permutation."""
    patterns = set()
    for values in product(range(-2, 3), repeat=6):
        if sum(values) == total and sum(x * x for x in values) == energy:
            patterns.add(tuple(sorted(x for x in values if x)))
    return patterns


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 43])
def test_balancing_and_defect_budgets(p):
    assert candidate_scaled_norm(p) == 2 * (p + 3)
    assert candidate_common_sum_magnitudes(p) == (0, 2, p - 1, p + 1)
    assert excess_over_balancing(p, 2) == 4
    assert excess_over_balancing(p, p - 1) == 4
    assert excess_over_balancing(p, p + 1) == 2


def test_all_low_energy_integer_profile_patterns_are_accounted_for():
    assert coefficient_patterns(1, 1) == {(1,)}
    assert coefficient_patterns(1, 3) == {(-1, 1, 1)}
    assert coefficient_patterns(1, 5) == {
        (-1, 2),
        (-1, -1, 1, 1, 1),
    }
    assert coefficient_patterns(2, 2) == {(1, 1)}
    assert coefficient_patterns(2, 4) == {
        (2,),
        (-1, 1, 1, 1),
    }
    assert coefficient_patterns(2, 6) == {
        (-1, 1, 2),
        (-1, -1, 1, 1, 1, 1),
    }


def test_unsigned_pair_moment_recurrences():
    for a, b in ((0, 1), (2, 7), (-3, 11), (5, -8)):
        q1 = a + b
        q2 = a * a + b * b
        q3 = a**3 + b**3
        q4 = a**4 + b**4
        assert unsigned_pair_cubic_defect(q1, q2, q3) == 0
        assert unsigned_pair_quartic_defect(q1, q2, q4) == 0


def test_t2_exception_factorizations():
    for values in ((0, 1, 2), (3, -2, 7), (11, 4, -5)):
        assert t2_energy6_doubled_cubic_defect(
            *values
        ) == t2_energy6_doubled_cubic_factor(*values)
    for values in ((0, 1, 2, 3), (1, 4, 7, 9), (-3, 2, 5, 11)):
        assert t2_energy4_nondegenerate_cubic_defect(
            *values
        ) == t2_energy4_nondegenerate_cubic_factor(*values)


def test_sum_one_exception_factorizations():
    for values in ((0, 1, 2), (3, -2, 7), (11, 4, -5)):
        assert sum_one_energy3_defect(*values) == sum_one_energy3_factor(
            *values
        )
    for values in ((0, 1), (3, -2), (11, 4)):
        assert sum_one_doubled_energy5_defect(
            *values
        ) == sum_one_doubled_energy5_factor(*values)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 43, 101])
def test_two_root_half_conic_forces_hasse_contradiction(p):
    cert = t2_hasse_obstruction(p)
    assert cert["required_abs_character_sum"] == p - 3
    assert cert["required_value_squared"] > cert["hasse_upper_bound_squared"]
    assert cert["contradiction"]


def test_p11_finite_half_conic_audit():
    row = audit_prime(11)
    assert row["quadratic_forms_checked"] == 11**3 - 1
    assert row["forbidden_two_root_nonnegative_forms"] == 0
    assert row["checks"]


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 43])
def test_complete_candidate_shell_is_excluded(p):
    assert candidate_shell_excluded(p)


def test_theorem_scope_keeps_the_actual_tail_open():
    theorem = candidate_shell_gap_theorem()
    assert theorem["proved"]
    assert "scaled norm 2(p+3)" in theorem["scope"]
    assert all(row["checks"] for row in theorem["rows"].values())
