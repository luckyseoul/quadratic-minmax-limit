import json
from pathlib import Path

from e1_gmin_m4_prop15747 import (
    EXPECTED_MODEL_SHA256,
    PARALLEL_COUNTS,
    build_mass12_height_four_model,
    mass12_boolean_second_moment_exclusion,
    mass12_height_four_arithmetic,
    proposition_15747,
)


ROOT = Path(__file__).resolve().parents[1]


def test_boolean_mass12_cut_second_moment_is_impossible_mod_seven():
    row = mass12_boolean_second_moment_exclusion()
    assert row["proved"] is True
    assert row["one_edge_cut_probability"] == "7/13"
    assert row["two_adjacent_edges_both_cut_probability"] == "7/26"
    assert row["two_disjoint_edges_both_cut_probability"] == "42/143"
    assert row["boolean_support_size"] == 396
    assert row["general_residual_is_one_mod_7"] is True
    assert row["all_integral_Q_boolean_mass12_lifts_excluded"] is True
    assert row["parallel_rows"]["3"]["required_cut_second_moment"] == "552/13"
    assert row["parallel_rows"]["5"]["required_cut_second_moment"] == "748/13"
    assert all(
        item["residual_mod_7"] == 1 and item["impossible_mod_7"]
        for item in row["parallel_rows"].values()
    )


def test_height_four_models_have_exact_projected_dimensions_and_hashes():
    expected = {
        3: {"sum": -12, "l1": 58, "cuts": (-14, -6)},
        5: {"sum": 14, "l1": 56, "cuts": (0, 8)},
    }
    for q in PARALLEL_COUNTS:
        arithmetic = mass12_height_four_arithmetic(q)
        model, weights, metadata = build_mass12_height_four_model(q)
        assert arithmetic["proved"] is True
        assert arithmetic["coefficient_sum"] == expected[q]["sum"]
        assert arithmetic["l1_budget"] == expected[q]["l1"]
        assert (
            arithmetic["cut_lower_at_C4"],
            arithmetic["cut_upper_at_C0"],
        ) == expected[q]["cuts"]
        assert arithmetic["derived_value_sum"] == 396
        assert arithmetic["height_four_anchor_is_wlog_before_field_moments"]
        assert len(weights) == 78
        assert metadata["integer_variable_count"] == 169
        assert metadata["constraint_count"] == 3526
        assert metadata["model_textproto_sha256"] == EXPECTED_MODEL_SHA256[q]
        assert model.Validate() == ""


def test_p3_branch_closed_and_p5_minimum_cells_forced_literal():
    row = proposition_15747()
    assert row["proved"] is True
    assert row["p13_t4_u4_P3_branch_closed"] is True
    assert row["P5_Q3_minimum_cells_forced_literal"] is True
    assert row["P5_minimum_literal_count_at_least"] == 2
    assert row["p13_t4_u4_closed"] is False
    for q in PARALLEL_COUNTS:
        exclusion = row["height_four_exclusions"][str(q)]
        assert exclusion["proved"] is True
        assert exclusion["solver"]["status"] == "INFEASIBLE"
        assert exclusion["solver"]["num_search_workers"] == 1
        assert exclusion["model"]["model_textproto_sha256"] == EXPECTED_MODEL_SHA256[q]

    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15747.json").read_text()
    )
    assert evidence == row
