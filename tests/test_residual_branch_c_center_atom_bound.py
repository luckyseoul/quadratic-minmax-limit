import numpy as np

from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    PHYSICAL_CENTERS,
    PHYSICAL_FIXED_POINT,
)
from scripts.residual_branch_c_center_atom_bound import (
    _normalized_rows,
    fixed_edges_in_target_direction,
    score_all_fixed_edges,
    score_centers,
)
from scripts.residual_branch_c_center_transverse_gpu import (
    _cpu_tensor,
    build_exact_tables,
)


def test_published_center_graph_has_exact_stronger_row_score() -> None:
    result = score_centers(PHYSICAL_CENTERS, PHYSICAL_FIXED_POINT)
    assert result["status"] == "EXACT_NECESSARY_SCORE"
    assert result["graph_sha256"] == (
        "c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d"
    )
    assert result["geometry"]["ternary"] is True
    assert result["geometry"]["cancellation_units"] == 1
    assert result["geometry"]["cancellation_multiplicity"] == 2
    assert result["geometry"]["cancellation_half_indices"] == (2, 13)
    assert result["geometry"]["fixed_direction_compatible"] is True
    assert result["violating_row_count"] == 32
    assert result["total_incremental_search_cost"] == 30132
    assert result["maximum_row_edit_lower_bound"] == 194
    assert result["all_496_cut_row_bounds_pass"] is False


def test_all_fifteen_fixed_edge_magnitudes_are_scored_exactly() -> None:
    edges = fixed_edges_in_target_direction()
    assert len(edges) == len(set(edges)) == 15
    result = score_all_fixed_edges(PHYSICAL_CENTERS)
    assert result["status"] == "EXACT_ALL_FIXED_EDGE_FRONTIER"
    assert result["fixed_edge_count"] == 15
    assert len(result["records"]) == 15
    assert result["best"]["total_incremental_search_cost"] == 29910
    assert result["best"]["maximum_row_edit_lower_bound"] == 194


def test_additive_center_tensor_equals_actual_graph_on_all_transverse_cells() -> None:
    tables = build_exact_tables()
    fixed_edges = tables["fixed_edges"]
    fixed = tuple(
        sorted(
            (
                PHYSICAL_FIXED_POINT,
                (-PHYSICAL_FIXED_POINT[0] % 31, -PHYSICAL_FIXED_POINT[1] % 31),
            )
        )
    )
    fixed_index = fixed_edges.index(fixed)
    choices = np.asarray([center - 1 for center in PHYSICAL_CENTERS], dtype=np.int16)
    additive = _cpu_tensor(choices, fixed_index, tables)

    _graph, rows = _normalized_rows(PHYSICAL_CENTERS, fixed)
    cells = tables["cells"]
    direct = np.zeros_like(additive)
    for direction_index, _parallel, coefficients in rows:
        direct[direction_index] = np.asarray(
            [coefficients.get(cell, 0) for cell in cells], dtype=np.int16
        )
    assert additive.shape == (32, 465)
    assert np.array_equal(additive, direct)


def test_rejected_center_tuple_is_not_scored_as_a_graph() -> None:
    result = score_centers((1,) * 16)
    assert result["status"] == "REJECTED_BY_TERNARITY_OR_TOP_CANCELLATION"
    assert result["geometry"]["admissible_top_geometry"] is False
