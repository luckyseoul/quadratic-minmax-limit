from __future__ import annotations

import pytest

from e1_gmin_m4_prop15766 import (
    beta_p,
    closed_t_max,
    directional_identity_ledger,
    even_band_row,
    odd_residual_nozero_profile,
    theorem_record,
    zero_direction_rigidity,
)


@pytest.mark.parametrize(
    ("p", "beta", "t_max"),
    [
        (11, 6, 3),
        (13, 6, 3),
        (17, 8, 3),
        (19, 10, 6),
        (23, 12, 7),
        (29, 14, 9),
        (31, 16, 9),
        (37, 18, 11),
        (41, 20, 12),
        (43, 22, 12),
    ],
)
def test_closed_ranges(p: int, beta: int, t_max: int) -> None:
    assert beta_p(p) == beta
    assert closed_t_max(p) == t_max
    for t in range(1, t_max + 1):
        row = even_band_row(p, t)
        assert row["proved"]
        assert row["isolated_vertex_margin"] > 0
        assert row["maximum_possible_residue_sum_from_average"] < row[
            "if_both_types_no_zero_required_residue_sum"
        ]
        assert row["opposite_type_average"] < row["opposite_least_positive"]


def test_directional_type_budget_and_zero_rigidity() -> None:
    row = directional_identity_ledger(29, 9, 1)
    assert row["type_sum"] == 30 * 9
    assert row["type_average"] == 18
    assert row["a_d_even"]
    rigid = zero_direction_rigidity(29, 9)
    assert rigid["proved"]
    assert rigid["P_d"] == 4
    assert rigid["strict_edge_gap"] > 0


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31])
@pytest.mark.parametrize("t", [2, 3, 7])
def test_odd_residual_exact_nozero_method_barrier(p: int, t: int) -> None:
    row = odd_residual_nozero_profile(p, t)
    assert row["proved"]
    assert row["all_directions_nonzero"]
    assert not row["Paley_graph_realizability_claimed"]
    assert sum(row["plus_type"]["a_values"]) == row["target_type_sum"]
    assert sum(row["minus_type"]["a_values"]) == row["target_type_sum"]


def test_theorem_scope_stays_narrow() -> None:
    theorem = theorem_record()
    assert theorem["proved"]["even_bridge_band"]
    assert not theorem["proved"]["all_even_sizes"]
    assert not theorem["proved"]["p5_p7"]
    assert not theorem["proved"]["residual_ii"]
    assert not theorem["proved"]["minimal_gap4_shell_bridge_closed_general"]
    assert not theorem["proved"]["e1_closed_general"]
