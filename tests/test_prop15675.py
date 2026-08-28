import json
from pathlib import Path

from e1_gmin_m4_prop15675 import (
    first_even_survivor,
    first_survivor_gap,
    phase_one_quantized_minimum,
    phase_zero_quantized_minimum,
    quantized_type_minimum_dp,
    symbolic_mod8_ledger,
    theorem_record,
)


def test_first_even_survivor_in_each_mod_eight_class():
    assert first_even_survivor(41) == 32
    assert first_even_survivor(19) == 14
    assert first_even_survivor(29) == 22
    assert first_even_survivor(23) == 18


def test_symbolic_type_minima_match_independent_dp():
    for p in (19, 23, 29, 31, 37, 41, 43):
        zero = phase_zero_quantized_minimum(p)
        one = phase_one_quantized_minimum(p)
        assert zero["u"] == 2
        assert zero["u_zero_strict_gap"] > 0
        assert one["u"] == (p - 1) // 2
        assert quantized_type_minimum_dp(p, 0)["minimum_deficit"] == zero[
            "minimum_deficit"
        ]
        assert quantized_type_minimum_dp(p, 1)["minimum_deficit"] == one[
            "minimum_deficit"
        ]


def test_exact_mod_eight_gap_formulas():
    expected = {
        41: -10,
        19: 10,
        29: 14,
        23: -4,
        43: 22,
        53: 26,
        47: -10,
    }
    for p, gap in expected.items():
        row = first_survivor_gap(p)
        assert row["gap"] == gap
        assert row["closed"] is (p % 8 in (3, 5))


def test_symbolic_scope_is_honest():
    row = symbolic_mod8_ledger()
    assert row["p=3 mod 8"]["closed"] is True
    assert row["p=5 mod 8"]["closed"] is True
    assert row["p=1 mod 8"]["closed"] is False
    assert row["p=7 mod 8"]["closed"] is False


def test_theorem_scope_and_generated_evidence():
    row = theorem_record()
    assert row["proved"] is True
    assert row["independent_dp_agrees"] is True
    assert row["theorem"]["p_mod_8_in_1_7"] == "OPEN_AT_THIS_RELAXATION"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15675.json").read_text()
    )
    assert stored == row
