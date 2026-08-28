import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15669 import (
    abstract_infinity_gap,
    abstract_no_infinity_gap,
    full_symbolic_floor,
    hypergeometric_moments,
    infinity_range_exclusion,
    infinity_type_deficit_ledger,
    largest_even_in_general_range,
    middle_floor_quadrature,
    no_infinity_range_exclusion,
    phase_zero_deficit_ledger,
    range_gap_concavity_ledger,
    small_prime_extensions,
    symbolic_floor_margin_ledger,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_hypergeometric_moments_and_all_sampled_middle_quadratures():
    for p in (17, 19, 23, 29, 31):
        assert symbolic_floor_margin_ledger(p)["all_strictly_positive"]
        assert range_gap_concavity_ledger(p)["positive_on_full_interval"]
        for b in range(5, p - 4):
            for phase in (0, 1):
                row = middle_floor_quadrature(p, b, phase)
                assert row["exact_positive_quadrature_certificate"]
                assert row["candidate_expectation"] == 1
                assert row["scaled_floor"] == 2 * p
                assert sum(row["quadrature_weights"].values(), Fraction()) == 1
                mean, second, _variance = hypergeometric_moments(p, b)
                assert sum(
                    Fraction(node) * weight
                    for node, weight in row["quadrature_weights"].items()
                ) == mean
                assert sum(
                    Fraction(node * node) * weight
                    for node, weight in row["quadrature_weights"].items()
                ) == second


def test_symbolic_floor_matches_generic_exact_lp_at_route_boundaries():
    for p in (17, 19, 23):
        probes = {
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            (p - 1) // 2,
            p - 6,
            p - 5,
            p - 1,
            p,
        }
        for b in sorted(probes):
            for phase in (0, 1):
                assert full_symbolic_floor(p, b, phase) == scaled_direction_floor(
                    p, b, phase
                )


def test_general_no_infinity_range_and_first_abstract_survivors():
    expected_first = {
        17: 14,
        19: 14,
        23: 18,
        29: 22,
        31: 24,
        37: 28,
        47: 36,
    }
    for p, first in expected_first.items():
        endpoint = largest_even_in_general_range(p)
        assert first == endpoint + 2
        row = no_infinity_range_exclusion(p, endpoint)
        assert row["excluded"]
        assert row["contradiction_gap"] > 0
        assert phase_zero_deficit_ledger(p, endpoint)["proved"]
        assert abstract_no_infinity_gap(p, endpoint) > 0
        assert abstract_no_infinity_gap(p, first) <= 0


def test_general_infinity_range_and_first_abstract_survivors():
    for p in (17, 19, 23, 29, 31, 37, 47):
        endpoint = p - 4
        for phase in (0, 1):
            row = infinity_range_exclusion(p, endpoint, phase)
            assert row["excluded"]
            assert abstract_infinity_gap(p, endpoint, phase) > 0
            assert abstract_infinity_gap(p, p - 2, phase) <= 0


def test_infinity_endpoint_complement_saving_is_accounted_for_exactly():
    # At p=17, phase one and b=p-4=13 has the complemented b=4 saving 6.
    # It does not invalidate the forced m-1 directions with b=1.
    row = infinity_type_deficit_ledger(17, 13, 1)
    assert row["savings_from_middle_baseline"][13] == 6
    assert row["maximum_other_saving"] == 6
    assert row["minimum_forced_b1_directions"] == 8
    assert row["saving_shortfall"] == 6 == row["symbolic_shortfall_lower_bound"]
    assert row["required_type_deficit_lower_bound"] == 96

    for p in (17, 19, 23, 29, 31, 37, 47):
        for s in range(5, p - 3, 2):
            for phase in (0, 1):
                ledger = infinity_type_deficit_ledger(p, s, phase)
                assert ledger["proved"]
                assert ledger["maximum_other_saving"] <= 6


def test_small_prime_extensions_are_exact_and_scoped():
    row = small_prime_extensions()
    assert row["p11_infinity_plus_7_gaps"] == {0: 30, 1: 18}
    assert row["p13_8_finite_gap"] == 4
    assert row["p13_infinity_gaps"] == {
        7: {0: 42, 1: 30},
        9: {0: 40, 1: 24},
    }
    assert row["p11_infinity_plus_7_excluded_both_phases"]
    assert row["p13_8_finite_excluded"]
    assert row["p13_infinity_plus_7_or_9_excluded_both_phases"]
    assert row["first_floor_pair_survivors"] == {
        "p11_8_finite_gap": -10,
        "p11_infinity_plus_9_gaps": {0: -72, 1: -72},
        "p13_10_finite_gap": -6,
        "p13_infinity_plus_11_gaps": {0: -110, 1: -110},
    }
    assert row["first_floor_pair_survivors_verified"]


def test_prop15669_record_keeps_every_top_level_gate_open():
    row = theorem_record()
    assert row["proved"]
    theorem = row["theorem"]
    assert theorem["general_residual_ii"] is False
    assert theorem["all_non_Walsh_multilevel"] is False
    assert theorem["R1"] is False
    assert theorem["global_QVAR"] is False
    assert theorem["type_I"] is False
    assert theorem["limit_exists"] is False
    assert row["remaining"]["floor_pair_survivor_is_not_an_actual_graph"]
    committed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15669.json").read_text()
    )
    assert committed == json.loads(json.dumps(row, default=str))
