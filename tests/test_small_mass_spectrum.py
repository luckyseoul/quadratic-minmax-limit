"""Sharp mass threshold, both parity phases, and fail-closed dependencies."""
from fractions import Fraction
import json
from pathlib import Path

import pytest

import e1_gmin_m4_small_mass_spectrum as proof


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", params=[29, 31])
def record(request):
    return proof.small_mass_spectrum(request.param)


def test_strict_sharp_local_spectrum_in_both_congruence_classes(record):
    p = record["p"]
    assert record["strict_upper_mass"] == record["nonboolean_mass_lower_bound"] == 2 * p - 10
    assert record["allowed_positive_masses"] == [p - 3, p + 1]
    assert record["boolean_below_strict_upper_mass"]
    assert record["proved"]
    assert not record["new_catalog_or_prime_census_used"]
    assert not record["residual_ii_closed"]
    assert not record["limit_closed"]


def test_two_stage_height_bootstrap_uses_strict_endpoint(record):
    p, height = record["p"], record["height_gap"]
    assert Fraction(height["first_height_strictly_greater_than"]) == 3
    assert height["half_mean_cube_maximum"] == 3
    assert height["all_maximizing_cubes_refined_mean_at_least"] == "3/4"
    assert Fraction(height["refined_height_strictly_greater_than"]) == Fraction(p + 13, 4) > 6
    assert height["three_quarter_mean_cube_maximum"] == 6
    assert height["all_maximizing_cubes_final_mean_at_least"] == 1
    endpoint_average = Fraction(height["paired_cube_average_endpoint_upper"])
    assert endpoint_average <= Fraction(2 * p - 10, 2 * (p - 1)) < 1
    assert height["all_prime_congruence_classes"]


def test_uniform_seven_coordinate_bound_and_fixed_four_bit_spectrum(record):
    p, boolean = record["p"], record["boolean_spectrum"]
    upper = Fraction(boolean["largest_invariant_class_complement_bound"])
    assert upper == Fraction(8 * (p - 1) * (p - 2), (p + 1) * (p - 3)) < 8
    assert 8 - upper == Fraction(8 * (p - 5), (p + 1) * (p - 3))
    assert boolean["junta_coordinates_at_most"] == 7 < (p - 1) // 2
    assert boolean["all_kept_coordinate_patterns_extend_to_slice"]
    assert boolean["symmetrization_preserves_degree_at_most_two"]
    assert boolean["cube_coordinates_actually_needed_at_most"] == 4
    assert boolean["all_boolean_scaled_masses"] == [
        0, p - 3, p + 1, 2 * p - 2, 2 * p + 2, 3 * p - 1, 3 * p + 3, 4 * p]
    catalog = boolean["fixed_catalog_dependency"]
    assert catalog["valid_tables"] == 222
    assert catalog["profile_count"] == 14
    assert catalog["evidence_sha256"] == proof.CATALOG_EVIDENCE_SHA256
    assert catalog["catalog_rerun"] is False


def test_existing_catalog_enumerator_is_never_called(monkeypatch):
    import e1_gmin_m4_prop15751 as prior

    def forbidden():
        raise AssertionError("old fixed catalog must not be rerun")

    monkeypatch.setattr(prior, "exact_four_cube_catalog", forbidden)
    assert proof.small_mass_spectrum(29)["proved"]


def test_both_boolean_masses_and_nonboolean_endpoint_are_attained(record):
    p = record["p"]
    assert [row["scaled_mass"] for row in record["attaining_boolean_examples"]] == [p - 3, p + 1]
    endpoint = record["sharp_endpoint"]
    assert endpoint["strict_upper_endpoint_is_necessary"]
    assert not endpoint["classification_of_endpoint_equalities_claimed"]
    assert [row["support_size"] for row in endpoint["examples"]] == [4, 5]
    assert [row["layer_values"] for row in endpoint["examples"]] == [
        [3, 1, 0, 0, 1], [3, 1, 0, 0, 1, 3]]
    for row in endpoint["examples"]:
        assert row["height"] == 3
        assert 4 * p * Fraction(row["mean"]) == row["scaled_mass"] == 2 * p - 10


@pytest.mark.parametrize("p", [29, 31])
def test_both_phase_spectra_and_genuine_parity_baselines(p):
    row = proof.affine_parity_small_mass_spectrum(p)
    assert row["proved"]
    assert row["strict_upper_mass"] == 2 * p - 10
    assert row["union_allowed_masses"] == [0, p - 3, p - 1, p + 1]
    assert row["phase_zero_allowed_masses"] == (
        [0, p - 3, p - 1, p + 1] if p % 4 == 1 else [0, p - 3, p + 1])
    assert row["phase_one_allowed_masses"] == (
        [p - 1, p + 1] if p % 4 == 1 else [p - 1])
    assert row["only_genuine_pointwise_parity_minima_subtracted"]
    assert row["no_punctured_complement_triple_difference_used"]
    for phase, data in row["phases"].items():
        assert data["candidate_boundary_sizes"] == ([0, 2, p - 1] if phase == "0" else [2, p - 1])
        for baseline in data["pointwise_baselines"]:
            assert baseline["pointwise_parity_minimum"]
            assert baseline["C_equals_half_difference_is_nonnegative_integral_quadratic"]
            assert baseline["lift_mass_equals_a_minus_baseline_mass"]
            assert 0 < baseline["positive_lift_mass_strictly_below"] < p - 3
            assert baseline["positive_lift_excluded"]
        pair, last = data["pointwise_baselines"]
        assert pair["scaled_mass"] == (p + 1 if phase == "0" else p - 1)
        assert pair["truth_values"] == ([0, 1, 1, 0] if phase == "0" else [1, 0, 0, 1])
        expected_last = "x_j" if (int(phase) + (p + 1) // 2) % 2 == 0 else "1-x_j"
        assert last["formula"] == expected_last


@pytest.mark.parametrize("p", [29, 31])
def test_local_mass_exclusion_accepts_only_the_open_interval(p):
    for mass in (p + 3, p + 5, p + 7, p + 9, 2 * p - 11,
                 Fraction(4 * p - 21, 2)):
        row = proof.local_mass_exclusion(p, mass)
        assert row["proved"] and row["excluded"]
    for mass in (0, p - 3, p - 1, p + 1, 2 * p - 10, 2 * p - 9):
        row = proof.local_mass_exclusion(p, mass)
        assert row["proved"] is False and row["excluded"] is False


@pytest.mark.parametrize("p", [True, False, 0, 13, 23, 25, 33, 49, 29.0])
def test_public_apis_reject_primes_outside_the_claim(p):
    for function in (proof.small_mass_spectrum, proof.affine_parity_small_mass_spectrum):
        with pytest.raises(ValueError, match="prime p>=29"):
            function(p)
    with pytest.raises(ValueError, match="prime p>=29"):
        proof.local_mass_exclusion(p, 34)


@pytest.mark.parametrize("mass", [True, False, 34.0, "34", None])
def test_mass_api_requires_exact_rationals(mass):
    with pytest.raises(ValueError, match="exact integer or Fraction"):
        proof.local_mass_exclusion(29, mass)


@pytest.mark.parametrize("dependency", ["sharp_integral_quadratic_lift_floor",
                                       "cube_half_mean_height_certificate",
                                       "cube_three_quarter_height_certificate"])
def test_missing_cube_or_stabilizer_proof_fails_closed(monkeypatch, dependency):
    monkeypatch.setattr(proof, dependency, lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="cube/stabilizer proof"):
        proof.small_mass_spectrum(29)


def test_wrong_half_mean_height_cannot_certify_bootstrap(monkeypatch):
    actual = proof.cube_half_mean_height_certificate()
    monkeypatch.setattr(proof, "cube_half_mean_height_certificate",
                        lambda: {**actual, "maximum_upper_bound": 4})
    with pytest.raises(ArithmeticError, match="dimension-free cube"):
        proof.small_mass_spectrum(29)


def test_wrong_stabilizer_coefficient_fails_closed(monkeypatch):
    actual = proof.sharp_integral_quadratic_lift_floor(29)
    monkeypatch.setattr(proof, "sharp_integral_quadratic_lift_floor",
                        lambda p: {**actual, "H_at_least_two_stabilizer_coefficient": Fraction(4)})
    with pytest.raises(ArithmeticError, match="paired-cube/stabilizer statement"):
        proof.small_mass_spectrum(29)


def test_pinned_catalog_hash_is_enforced(monkeypatch):
    monkeypatch.setattr(proof, "CATALOG_EVIDENCE_SHA256", "0" * 64)
    with pytest.raises(ArithmeticError, match="evidence hash"):
        proof.small_mass_spectrum(29)


def test_missing_catalog_proof_fails_closed(monkeypatch):
    monkeypatch.setattr(proof, "_fixed_catalog_dependency", lambda: {"proved": False})
    with pytest.raises(ArithmeticError, match="fixed Boolean catalog"):
        proof.small_mass_spectrum(29)


def test_wrong_density_conversion_fails_closed(monkeypatch):
    monkeypatch.setattr(proof, "profile_density", lambda *args: Fraction(0))
    with pytest.raises(ArithmeticError, match="density spectrum"):
        proof.small_mass_spectrum(29)


@pytest.mark.parametrize("dependency", ["_height_gap_certificate", "_boolean_spectrum_certificate",
                                       "_sharp_endpoint_certificate"])
def test_top_level_theorem_checks_each_branch_proof(monkeypatch, dependency):
    monkeypatch.setattr(proof, dependency, lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="small-mass theorem dependency"):
        proof.small_mass_spectrum(29)


@pytest.mark.parametrize("consumer", ["affine_parity_small_mass_spectrum", "local_mass_exclusion"])
def test_corollaries_require_the_local_spectrum(monkeypatch, consumer):
    monkeypatch.setattr(proof, "small_mass_spectrum", lambda p: {"proved": False})
    args = (29, 34) if consumer == "local_mass_exclusion" else (29,)
    with pytest.raises(ArithmeticError, match="spectrum dependency|exclusion dependency"):
        getattr(proof, consumer)(*args)


def test_affine_parity_requires_floor_proof(monkeypatch):
    monkeypatch.setattr(proof, "residual_even_floor_table", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="affine-parity spectrum dependency"):
        proof.affine_parity_small_mass_spectrum(29)


def test_affine_parity_rejects_a_new_low_boundary(monkeypatch):
    actual = proof.residual_even_floor_table(29)
    altered = {**actual, "phase_zero_floors": {**actual["phase_zero_floors"], 4: 30}}
    monkeypatch.setattr(proof, "residual_even_floor_table", lambda p: altered)
    with pytest.raises(ArithmeticError, match="boundary partition"):
        proof.affine_parity_small_mass_spectrum(29)


def test_affine_parity_does_not_accept_a_nonpointwise_baseline(monkeypatch):
    actual = proof._parity_baseline
    monkeypatch.setattr(proof, "_parity_baseline",
                        lambda *args: {**actual(*args), "pointwise_parity_minimum": False})
    with pytest.raises(ArithmeticError, match="genuine parity-baseline lift"):
        proof.affine_parity_small_mass_spectrum(29)


def test_affine_parity_requires_baseline_proof(monkeypatch):
    monkeypatch.setattr(proof, "_parity_baseline", lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="parity-baseline proof dependency"):
        proof.affine_parity_small_mass_spectrum(29)


def test_saved_evidence_equals_the_live_theorem_and_keeps_global_gates_open():
    saved = json.loads((ROOT / "evidence" / "e1_gmin_m4_small_mass_spectrum.json").read_text())
    assert saved == proof.theorem_record()
    assert saved["proved"]
    assert saved["residual_ii_closed"] is False
    assert saved["E1_closed"] is False
    assert saved["quadratic_minmax_limit_closed"] is False
