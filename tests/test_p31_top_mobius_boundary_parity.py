import hashlib
import json
from pathlib import Path

from e1_gmin_m4_mobius_half_symmetric import (
    _relative_coefficients,
    mobius_parameter_edges,
    paley_edge_sign,
)
from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    HALVES,
    PHYSICAL_CENTERS,
)
from e1_gmin_m4_p31_top_mobius_boundary_parity import (
    P,
    adaptive_design_boundary_parity_certificate,
    all_prime_adaptive_product_theorem,
    design_kernel_parity_certificate,
    fixed_antipodal_edges,
    half_projective_signature_certificate,
    half_kernel_parity,
    inversion_orbit_pair_signature,
    kernel_selector,
    selector_pairing,
    boundary_mask,
)
from e1_gmin_m4_inversion_antisymmetric_radon import projective_functionals


ROOT = Path(__file__).resolve().parents[1]
GLOBAL_REPLAY = (
    ROOT
    / "evidence"
    / "e1_gmin_m4_p31_top_mobius_boundary_parity_global_replay.json"
)
GLOBAL_REPLAY_SHA256 = (
    "9ce413b62c65b7a541388f49bb390c690af44389556e6e8676ec138ef2cbc533"
)
COMPONENT_REPLAY = (
    ROOT
    / "evidence"
    / "e1_gmin_m4_p31_top_mobius_boundary_parity_component_replay.json"
)
COMPONENT_REPLAY_SHA256 = (
    "1e6cb15792a681f6e9b02290e60e0e14f16225282c56258c562564d1b0556650"
)


def test_frozen_top_design_has_exact_kernel_boundary_obstruction() -> None:
    row = design_kernel_parity_certificate(
        HALVES,
        5,
        centers=PHYSICAL_CENTERS,
        cancelled_edge=((2, 25), (29, 1)),
    )
    assert row["proved"]
    assert row["selector_size"] == 16
    assert row["half_kernel_parity_bits"] == (
        1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0
    )
    assert row["half_parity_sum_mod_two"] == 0
    assert row["cancelled_orbit_removal_pairing"] == 0
    assert row["necessary_total_pairing"] == 1
    assert row["boundary_parity_obstruction"]


def test_every_center_and_fixed_edge_has_the_replayed_pairing() -> None:
    directions = projective_functionals(P)
    fixed = directions[5]
    selector = kernel_selector(fixed)
    for target, auxiliary in HALVES:
        assert {
            half_kernel_parity(target, auxiliary, fixed, center)
            for center in range(1, P)
        } == {half_kernel_parity(target, auxiliary, fixed)}
    assert {
        selector_pairing(boundary_mask((edge,)), selector)
        for edge in fixed_antipodal_edges(fixed)
    } == {1}


def test_generic_kernel_bit_matches_the_unique_endpoint_paley_sign() -> None:
    directions = projective_functionals(P)
    fixed = directions[5]
    for target, auxiliary in HALVES:
        a, b = _relative_coefficients(P, target, auxiliary, fixed)
        assert a and b and (a + b) % P
        parameter = -a * pow(a + b, -1, P) % P
        tau = paley_edge_sign(
            P, mobius_parameter_edges(P, target, auxiliary, 1)[parameter]
        )
        assert half_kernel_parity(target, auxiliary, fixed) == int(tau == -1)


def test_all_prime_adaptive_product_algebra() -> None:
    for p in (3, 7, 11, 19, 23, 31, 43):
        row = all_prime_adaptive_product_theorem(p)
        assert row["proved"]
        assert row["half_count_is_even"]
        assert row["one_half_signature_product_formula"] == "-epsilon_M"
        assert row["auxiliary_sign_product"] == 1
        assert row["all_half_signature_product"] == 1
        assert row["required_signature_product"] == -1
        assert row["some_projective_kernel_is_a_contradiction"]


def test_each_frozen_half_has_the_full_projective_product_formula() -> None:
    directions = projective_functionals(P)
    for target, auxiliary in HALVES:
        row = half_projective_signature_certificate(target, auxiliary)
        assert row["proved"]
        assert row["all_centers_have_same_signature"]
        assert row["closed_formula_matches_every_direction"]
        assert row["sigma_at_target"] == -1
        assert row["sigma_at_l_minus_m"] == row["target_sign"] == 1
        assert row["signature_product"] == -row["auxiliary_sign"]
        assert len(row["kernel_sigma"]) == len(directions) == P + 1


def test_frozen_design_has_adaptive_all_design_obstruction() -> None:
    row = adaptive_design_boundary_parity_certificate(
        HALVES,
        5,
        centers=PHYSICAL_CENTERS,
        cancelled_edge=((2, 25), (29, 1)),
    )
    assert row["proved"]
    assert row["top_profile_multisets_match"]
    assert row["top_opposite_fixed_sdr"]
    assert (row["auxiliary_hard_count"], row["auxiliary_opposite_count"]) == (
        14,
        2,
    )
    assert row["auxiliary_sign_product"] == 1
    assert row["half_kernel_signature_product"] == 1
    assert row["fixed_edge_signature_product"] == -1
    assert row["origin_orbit_count"] == 16
    assert row["origin_orbits_distinct"]
    assert row["nonorigin_orbit_pair_aggregate_parity"] == 0
    assert row["final_signature_product_after_arbitrary_nonorigin_reductions"] == -1
    assert row["adaptive_kernel_product_obstruction"]
    assert row["top_family_boundary_parity_obstruction"]
    contradiction_directions = row[
        "clean_collision_contradiction_direction_indices"
    ]
    assert len(contradiction_directions) % 2 == 1
    assert 20 in contradiction_directions


def test_nonorigin_pair_is_even_and_origin_pair_is_the_exact_caveat() -> None:
    nonorigin = inversion_orbit_pair_signature(((2, 25), (29, 1)))
    origin = inversion_orbit_pair_signature(((0, 0), (1, 0)))
    assert sum(nonorigin) == 2
    assert sum(nonorigin) % 2 == 0
    assert sum(origin) == 1
    assert sum(origin) % 2 == 1


def test_exhaustive_global_half_option_replay_is_pinned() -> None:
    raw = GLOBAL_REPLAY.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == GLOBAL_REPLAY_SHA256
    row = json.loads(raw)
    assert row["p"] == 31
    assert row["half_option_count"] == 14_880
    assert row["center_half_instance_count"] == 446_400
    assert row["selector_pairing_count"] == 14_284_800
    assert row["center_variation_count"] == 0
    assert row["centrality_failure_count"] == 0
    assert row["origin_failure_count"] == 0
    assert row["selector_partition_failure_count"] == 0
    assert row["product_failure_count"] == 0
    assert row["proved_for_enumerated_p31_half_options"]
    assert not row["residual_ii_closed"]


def test_bounded_component_adaptive_witness_histogram_is_pinned() -> None:
    raw = COMPONENT_REPLAY.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == COMPONENT_REPLAY_SHA256
    row = json.loads(raw)
    assert row["input_catalog_sha256"] == (
        "196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b"
    )
    assert row["searchable_design_count"] == 2_969
    assert row["unique_half_choice_count"] == 225
    assert row["center_invariance_failure_count"] == 0
    assert row["universal_single_direction_count"] == 0
    assert row["unwitnessed_design_count"] == 0
    assert row["minimum_witness_direction_count"] == 7
    assert sum(row["witness_direction_count_histogram"].values()) == 2_969
    assert row["all_10000_design_double_sum_parity_counts"] == [10_000, 0]
    assert not row["residual_ii_closed"]
