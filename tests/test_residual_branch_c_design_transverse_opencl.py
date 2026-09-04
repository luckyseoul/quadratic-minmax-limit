"""CPU exactness checks for the moved-design OpenCL table builder."""

from __future__ import annotations

from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    HALVES,
    direct_parallel_design_certificate,
)
from scripts.residual_branch_c_design_transverse_opencl import (
    build_design_tables,
    exact_design_graph_replay,
)


def test_generic_tables_reproduce_the_recorded_base_design_winner() -> None:
    profile = direct_parallel_design_certificate()
    design = {
        "halves": HALVES,
        "raw_parallel_profile": profile["raw_parallel_profile"],
        "fixed_direction_index": 5,
    }
    collision = {
        "half_indices": [2, 13],
        "centers": [28, 1],
        "orbit": [[2, 25], [29, 1]],
        "coefficients": [-1, 1],
        "spatial_direction_index": 5,
    }
    tables = build_design_tables(design, collision)
    result = exact_design_graph_replay(
        [8, 9, 28, 12, 7, 9, 11, 6, 3, 15, 12, 3, 15, 1, 10, 11],
        2,
        tables,
    )
    assert (
        result["total_positive_l1_excess"],
        result["maximum_l1_excess"],
    ) == (4604, 178)
    assert result["graph_sha256"] == (
        "46b4a125eee33a6b66a358b7a04bf9cb640a9f1398205e345c882f2d10a02889"
    )
    assert result["exact_integer_graph_replay"] is True
