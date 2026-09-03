import json
from pathlib import Path

from e1_gmin_m4_prop15761 import (
    compact_ray_norm_bounds,
    full_radon_spectrum_certificate,
    signed_least_norm_certificate,
    symbolic_gap_certificate,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_full_midpoint_spectrum_has_all_three_orthogonal_blocks():
    row = full_radon_spectrum_certificate()
    assert row["proved"]
    assert row["distinct_direction_row_intersections"] == {
        "parallel_parallel": 0,
        "parallel_pair": "p",
        "pair_pair": 2,
    }
    assert row["one_direction_row_norms_squared"] == {
        "parallel": "p^2*m",
        "off_diagonal_pair": "p^2",
    }
    assert set(row["orthogonal_blocks"]) == {
        "within_direction_pair_zero_sum",
        "directional_aggregate_zero_sum",
        "uniform",
    }
    assert all(row["symbolic_eigenvalue_checks"].values())


def test_signed_boolean_source_implies_the_exact_least_norm_bound():
    row = signed_least_norm_certificate()
    assert row["proved"]
    assert row["simple_source_squared_norm"] == "||z||^2=|H|"
    assert row["necessary_inequality"] == "Q<=|H|"
    assert row["strictly_finer_than_difference_aggregate_parseval"]


def test_balanced_compact_atom_bounds_are_uniform_in_labels():
    row = compact_ray_norm_bounds()
    assert row["proved"]
    assert row["uniform_in_atom_labels"]
    assert row["p_1_mod_4"]["balanced_count_bounds"] == {
        "e_L": "<=2r-3",
        "Q_L": "<=2r-1",
    }
    assert row["p_3_mod_4"]["balanced_count_bounds"] == {
        "e_L": "<=2r-2",
        "Q_L": "<=2r+2",
    }


def test_all_r_gaps_are_proved_by_positive_shifted_coefficients():
    row = symbolic_gap_certificate()
    assert row["proved"]
    assert row["no_prime_or_parameter_scan"]
    assert row["substitution"] == "u=r-7>=0"
    assert row["shifted_numerators"]["p_1_main_gap"] == [
        30706,
        29864,
        9796,
        1328,
        64,
    ]
    assert row["shifted_numerators"]["p_3_main_gap"] == [
        154867,
        85374,
        18108,
        1744,
        64,
    ]
    assert all(
        all(coefficient > 0 for coefficient in coefficients)
        for coefficients in row["shifted_numerators"].values()
    )


def test_checked_in_evidence_preserves_the_nonclosure_boundary():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15761.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["proved"]["both_prop15758_rays_pass_with_strict_room_for_r_ge_7"]
    assert not observed["proved"]["one_common_simple_graph_constructed"]
    assert not observed["proved"]["residual_ii_closed"]
    assert observed["L_status"] == "OPEN"
