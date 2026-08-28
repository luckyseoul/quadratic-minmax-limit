from fractions import Fraction

from e1_gmin_m4_prop15679 import (
    _relaxed_minimum_dp,
    large_prime_lift_ledger,
    next_boundary_exclusion,
    next_even_boundary,
    pair_row,
    phase_one_minimum,
    phase_zero_interior_minimum,
    phase_zero_u0_minimum,
    small_prime_lift_ledger,
    symbolic_residue_reduction,
    theorem_record,
)


def test_next_boundary_closed_forms_and_middle_range():
    expected = {43: 34, 47: 38, 53: 42, 59: 46, 73: 58}
    for p, s in expected.items():
        assert next_even_boundary(p) == s
        assert s <= p - 5


def test_phase_minima_match_independent_relaxed_dp():
    for p in (43, 47, 53, 59, 73):
        phase_zero = {
            row["u"]: row for row in _relaxed_minimum_dp(p, 0)["rows"]
        }
        phase_one = _relaxed_minimum_dp(p, 1)["rows"]
        assert phase_zero[0]["minimum_deficit"] == (
            phase_zero_u0_minimum(p)["minimum_deficit"]
        )
        for u in range(2, 8):
            assert phase_zero[u]["minimum_deficit"] == (
                phase_zero_interior_minimum(p, u)["minimum_deficit"]
            )
        assert len(phase_one) == 1
        assert phase_one[0]["minimum_deficit"] == (
            phase_one_minimum(p)["minimum_deficit"]
        )


def test_u0_and_u8_are_over_pair_budget():
    for p in (43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101):
        u0 = phase_zero_u0_minimum(p)
        d1 = phase_one_minimum(p)["minimum_deficit"]
        s = next_even_boundary(p)
        assert u0["minimum_deficit"] + d1 > s * (s - 1)
        assert pair_row(p, 8)["survives_pair_budget"] is False


def test_small_prime_pair_survivors_and_lift_floors():
    expected = {
        43: ([2, 3, 4], 12),
        47: ([2, 3, 4, 5, 6], 14),
        53: ([2, 3, 4, 5], 14),
    }
    for p, (residues, cost) in expected.items():
        row = small_prime_lift_ledger(p)
        assert row["surviving_residues"] == residues
        assert row["nonzero_quadratic_lift_floor"] == cost
        assert row["maximum_scaled_mean"] < cost
        assert row["excluded"] is True


def test_large_prime_support_floor_is_strictly_above_fourteen():
    row = large_prime_lift_ledger()
    assert Fraction(row["value_at_59"]) > 14
    assert row["excluded"] is True


def test_symbolic_reduction_and_samples_are_closed():
    reduction = symbolic_residue_reduction()
    assert reduction["only_residues_requiring_lift"] == [2, 3, 4, 5, 6, 7]
    assert reduction["proved"] is True
    for p in (43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101):
        assert next_boundary_exclusion(p)["excluded"] is True


def test_theorem_record_is_honestly_scoped():
    record = theorem_record()
    assert record["proved"] is True
    theorem = record["theorem"]
    assert theorem["smaller_endpoints_status"] == "OPEN_AT_THIS_BOUNDARY_SIZE"
    assert theorem["general_residual_ii"] is False
    assert theorem["R1"] is False
    assert theorem["limit_exists"] is False
