import json
from pathlib import Path

from e1_gmin_m4_prop15760 import (
    difference_quotient_certificate,
    integral_image_criterion,
    midpoint_kernel_certificate,
    ordinary_lattice_certificate,
    smith_extension_certificate,
    symbolic_rank_identities,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_all_rank_identities_are_symbolic_in_m_not_finite_prime_scans():
    row = symbolic_rank_identities()
    assert row["proved"]
    assert all(row["checks"].values())
    assert row["p"] == [1, 2]
    assert row["directions_p_plus_one"] == [2, 2]
    assert row["difference_class_count"] == [0, 2, 2]
    assert row["ordinary_lattice_rank"] == [0, 2, 6, 4]
    assert row["midpoint_kernel_rank"] == [0, 0, 4, 4]
    assert row["difference_quotient_rank"] == [0, 2, 2]
    assert row["difference_quotient_mod_p_rank"] == [1, 2, 1]
    assert row["difference_quotient_defect"] == [-1, 0, 1]
    assert row["midpoint_kernel_defect"] == [0, "-1/6", "-1/2", "2/3"]
    assert row["total_defect"] == [-1, "-1/6", "1/2", "2/3"]
    assert row["total_defect"] == row["prop15759_defect"]


def test_ordinary_lattice_is_a_primitive_integral_kernel_with_a_basis():
    row = ordinary_lattice_certificate()
    assert row["proved"]
    assert row["split_surjective_constraint_map"]
    assert row["ordinary_lattice_is_primitive"]
    assert row["surjection_witnesses"] == {
        "direction_coordinate_L": "K_L(c0)",
        "last_coordinate": "P_L-K_L(c0), for any L!=L0",
    }
    assert len(row["integral_basis"]) == 2


def test_both_exact_layers_have_elementary_p_cokernel():
    quotient = difference_quotient_certificate()
    midpoint = midpoint_kernel_certificate()
    assert quotient["exact_sequence"] == "0 -> M -> A -> A_D -> 0"
    assert quotient["S_is_injective"]
    assert quotient["p_times_quotient_lattice_in_image"]
    assert quotient["cokernel"] == "(Z/pZ)^(m^2-1)"
    assert midpoint["source_exact_sequence"] == "0 -> E_0 -> E -> Z^Delta -> 0"
    assert midpoint["other_directions_cancel_exactly"]
    assert midpoint["p_times_midpoint_lattice_in_image"]
    assert midpoint["cokernel"] == "(Z/pZ)^[m(m-1)(4m+1)/6]"


def test_snake_extension_has_no_hidden_p_squared_or_other_torsion():
    row = smith_extension_certificate()
    assert row["proved"]
    assert row["snake_tail"] == "0 -> coker R_0 -> coker R -> coker S -> 0"
    assert row["no_p_squared_torsion"]
    assert row["no_other_prime_torsion"]
    assert row["full_cokernel"] == "(Z/pZ)^[(m-1)(4m^2+7m+6)/6]"


def test_integral_sufficiency_does_not_promote_to_a_simple_graph():
    row = integral_image_criterion()
    assert row["proved"]
    assert row["all_linear_congruence_obstructions_exhausted"]
    assert row["simple_graph_normalized_domain"] == "z_e in {0,tau_e}"
    assert not row["simple_nonnegative_lift_proved"]
    assert not row["compact_prop15758_rays_pass_all_moments_proved"]


def test_checked_in_evidence_is_the_symbolic_theorem_record():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15760.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["prop"] == "15.760"
    assert observed["proved"]["prop15759_moments_are_integer_image_sufficient"]
    assert not observed["proved"]["compact_aggregate_survivor_has_simple_lift"]
    assert observed["L_status"] == "OPEN"
