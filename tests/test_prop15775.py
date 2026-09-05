"""Eventual first-layer theorem and its quantitative, not global, corollary."""
import json
from fractions import Fraction
from pathlib import Path

import pytest

import e1_gmin_m4_prop15775 as proof

ROOT = Path(__file__).resolve().parents[1]


def test_dimension_free_noise_and_integer_fourier_constants():
    rec = proof.cube_height_bound_certificate()
    assert rec["two_point_nonnegative_remainder"] == ["0", "0", "8/9"]
    assert rec["degree_two_L4_over_L2"] == "3"
    assert rec["L2_over_L1"] == "9"
    assert rec["Fourier_coefficient_lattice"] == "(1/4) Z"
    assert rec["height_over_mean_squared"] == 324
    assert rec["proof_is_dimension_free"] and rec["proved"]


@pytest.mark.parametrize("p", [29, 31, 259201, 524287, 6700417])
def test_paired_cube_and_middle_moments_exactly(p):
    rec = proof.paired_cube_basis_certificate(p)
    assert len(rec["monomial_orbits"]) == 6
    for row in rec["monomial_orbits"]:
        assert Fraction(row["paired_mean"]) == (
            row["at_X"] + p * Fraction(row["slice_mean"])) / (p + 1)
    mean = proof.middle_slice_mean_certificate(p)
    for a, b, c in zip(mean["slice_moments"], mean["cube_moments"], mean["at_zero"]):
        assert p * Fraction(a) == (p + 1) * Fraction(b) - c
    assert mean["full_cube_mean_lattice"] == "(1/4) Z"


@pytest.mark.parametrize("p", [29, 31, 259201, 524287, 6700417])
@pytest.mark.parametrize("s", [4, 6])
def test_height_and_strict_integer_influence_junta_bounds(p, s):
    height = proof.slice_height_bound_certificate(p, s)
    assert Fraction(height["mean"]) == Fraction(2 * p + s, 2 * p)
    assert Fraction(height["paired_cube_mean_upper"]) == Fraction(2 * p + s, p - 1)
    assert Fraction(height["height_upper_exact"]) <= Fraction(82944, 49) < 1800
    derivative = proof.derivative_support_certificate(p)
    assert Fraction(derivative["conditional_nonzero_support_lower"]) > Fraction(1, 4)
    assert Fraction(derivative["relevant_pair_influence_lower"]) == Fraction(
        p * p - 1, 32 * p * (p - 2)) > Fraction(1, 32)
    assert derivative["total_influence_variance_factor"] == p - 1
    junta = proof.bounded_junta_certificate(p, s)
    assert Fraction(junta["largest_invariant_class_complement_upper"]) < 129600
    assert junta["junta_coordinates_at_most"] == 129599
    assert junta["all_kept_patterns_extend"] == ((p - 1) // 2 >= 129599)


@pytest.mark.parametrize("p", [259201, 524287, 6700417])
@pytest.mark.parametrize("s", [4, 6])
def test_local_theorem_needs_no_prime_or_affine_phase(p, s):
    rec = proof.local_mass_exclusion(p, s)
    lo, hi = map(Fraction, rec["cube_nu_minus_one_interval"])
    assert lo == Fraction(s // 2 - 1, p + 1)
    assert hi == Fraction(1800 + s // 2 - 1, p + 1)
    assert 0 < lo <= hi < Fraction(1, 4)
    assert rec["affine_parity_hypothesis_needed"] is False
    assert rec["excluded"] and rec["proved"]


@pytest.mark.parametrize("p", [37, 43, 524287, 6700417])
def test_first_layer_exact_quota_without_pinching_P_or_T(p):
    rec = proof.first_layer_quota_certificate(p)
    m = (p + 1) // 2
    assert rec["H"] == 5 * p + 6 and rec["t_shell"] == 2 * m + 2
    assert rec["only_ordered_residue_pairs"] == [[1, 2], [2, 1]]
    assert rec["allowed_residue_intervals"] == [[0, 2], [m - 6, m - 1]]
    assert rec["u1_quotient_histogram"] == {"2": m - 1, "3": 1}
    assert rec["u2_quotient_histogram"] == {"2": m}
    assert rec["forced_low_masses"] == [2 * p + 4, 2 * p + 6]
    assert rec["P_and_Q_sum"] == 10
    assert rec["P_equals_5_or_T_equals_pm1_assumed"] is False
    assert rec["large_residue_or_boundary_scan_used"] is False


@pytest.mark.parametrize("p", proof.SAMPLE_PRIMES)
def test_whole_eventual_first_layer_is_excluded(p):
    rec = proof.eventual_first_layer_exclusion(p)
    assert rec["H"] == 5 * p + 6 and rec["k"] == 5 * p + 5
    assert rec["t_residual"] == (p - 1) // 2 + 3
    assert rec["isolation_margin"] == p * p + 1 - 2 * (5 * p + 6) > 0
    assert rec["both_signed_shell_floors"] == 3
    assert rec["all_boundary_sizes_excluded"] and rec["whole_layer_excluded"]


@pytest.mark.parametrize("r", [3, 4, 5])
def test_power_band_exact_extension_and_divisibility(r):
    p = 6700417
    h = r * p + 2
    rec = proof.power_band_exclusion(p, r, h)
    B = Fraction(h, 2 * p)
    assert 46656 * h**3 <= p**3 * (p - 1)
    assert Fraction(rec["type_average_A"]) == Fraction(1, p)
    assert Fraction(rec["junta_strict_upper"]) == 186624 * B**3 <= rec["q"]
    assert 0 < rec["positive_multiple_lower"] <= Fraction(rec["positive_multiple_upper"]) < p + 1
    assert rec["parity_required"] == "h=r mod 2"
    assert rec["full_cube_extensions_exist"] and rec["neither_phase_omitted"]
    assert rec["excluded"] and rec["proved"]


@pytest.mark.parametrize("r", [3, 4, 5])
def test_power_band_below_frame_mean(r):
    rec = proof.power_band_exclusion(29, r, r * 29 - 2)
    assert rec["case"] == "below the signed frame mean"
    assert rec["excluded"] and rec["proved"]


@pytest.mark.parametrize("p", [29, 524287, 6700417])
@pytest.mark.parametrize("B", [Fraction(9, 8), Fraction(3, 2), Fraction(5, 2)])
def test_generic_bounded_mean_constants_are_exact(p, B):
    rec = proof.generic_bounded_mean_certificate(p, B)
    assert Fraction(rec["mean_cap"]) == B
    assert Fraction(rec["height_upper"]) == 2916 * B**2
    assert Fraction(rec["junta_size_strict_upper"]) == 186624 * B**3
    assert 324 * (Fraction(2 * p, p - 1) * B)**2 <= Fraction(rec["height_upper"])
    assert rec["mean_dependency"]["full_cube_mean_lattice"] == "(1/4) Z"
    assert rec["proved"]


@pytest.mark.parametrize("B", [True, 0, -1, Fraction(-1, 2), 0.0, 1.0])
def test_generic_mean_cap_rejects_inexact_or_nonpositive_inputs(B):
    with pytest.raises(ValueError, match="positive exact mean cap"):
        proof.generic_bounded_mean_certificate(29, B)


@pytest.mark.parametrize("p,error", [
    (True, "odd p>=29"), (27, "odd p>=29"), (28, "odd p>=29"),
    (29.0, "odd p>=29"), (39, "prime p>=29"),
])
def test_power_band_rejects_unproved_orders(p, error):
    with pytest.raises(ValueError, match=error):
        proof.power_band_exclusion(p, 3, 1)


@pytest.mark.parametrize("r", [True, 2, 6, 3.0])
def test_power_band_rejects_unproved_or_inexact_shell_floors(r):
    with pytest.raises(ValueError, match="shell floor"):
        proof.power_band_exclusion(6700417, r, 1)


@pytest.mark.parametrize("h", [True, -1, 1.0, 2])
def test_power_band_rejects_negative_inexact_or_wrong_parity_support(h):
    with pytest.raises(ValueError, match="nonnegative h congruent to r modulo two"):
        proof.power_band_exclusion(6700417, 3, h)


@pytest.mark.parametrize("p", [True, 27, 28, 29.0])
def test_bad_local_orders_rejected(p):
    with pytest.raises(ValueError, match="odd p>=29"):
        proof.slice_height_bound_certificate(p, 4)


@pytest.mark.parametrize("s", [True, 2, 8, 4.0])
def test_unproved_local_shifts_rejected(s):
    with pytest.raises(ValueError, match="s=4,6"):
        proof.local_mass_exclusion(259201, s)


def test_threshold_prime_parity_and_support_domain_guards():
    with pytest.raises(ValueError, match="odd p>=259201"):
        proof.local_mass_exclusion(259199, 4)
    with pytest.raises(ValueError, match="prime p>=37"):
        proof.first_layer_quota_certificate(39)
    with pytest.raises(ValueError, match="prime p>=259201"):
        proof.eventual_first_layer_exclusion(300003)
    with pytest.raises(ValueError, match="modulo two"):
        proof.power_band_exclusion(6700417, 3, 3 * 6700417 + 1)
    with pytest.raises(ValueError, match="shell floor"):
        proof.power_band_exclusion(6700417, 6, 6 * 6700417)
    with pytest.raises(ValueError, match="cubic support band"):
        proof.power_band_exclusion(29, 3, 3 * 29)
    with pytest.raises(ValueError, match="positive exact mean cap"):
        proof.generic_bounded_mean_certificate(29, 1.0)


@pytest.mark.parametrize("constant,value,call,error", [
    ("CUBE_HEIGHT_CONSTANT", 323, lambda: proof.cube_height_bound_certificate(), "noise/Fourier"),
    ("HEIGHT_BOUND", 1799, lambda: proof.slice_height_bound_certificate(29, 6), "slice height"),
    ("MEAN_UPPER", Fraction(8, 7), lambda: proof.slice_height_bound_certificate(29, 6), "slice height"),
    ("JUNTA_STRICT_BOUND", 129599, lambda: proof.bounded_junta_certificate(29, 6), "complement bound"),
    ("PRIME_THRESHOLD", 259199, lambda: proof.local_mass_exclusion(259201, 4), "order threshold"),
])
def test_silent_constant_changes_are_not_certified(monkeypatch, constant, value, call, error):
    monkeypatch.setattr(proof, constant, value)
    with pytest.raises(ArithmeticError, match=error):
        call()


@pytest.mark.parametrize("name,field,value,call,error", [
    ("cube_height_bound_certificate", "proved", False,
     lambda: proof.slice_height_bound_certificate(29, 4), "height proof dependency"),
    ("sharp_integral_quadratic_lift_floor", "paired_cube_identity", "wrong",
     lambda: proof.slice_height_bound_certificate(29, 4), "15.688"),
    ("stabilizer_mass_certificate", "weights", (Fraction(0), Fraction(0), Fraction(0)),
     lambda: proof.slice_height_bound_certificate(29, 4), "15.642"),
    ("derivative_support_certificate", "uniform_strict_influence_lower", "1/16",
     lambda: proof.bounded_junta_certificate(29, 4), "height/influence dependency"),
    ("bounded_junta_certificate", "all_kept_patterns_extend", False,
     lambda: proof.local_mass_exclusion(259201, 4), "mean/junta/extension"),
    ("middle_slice_mean_certificate", "identity", "wrong",
     lambda: proof.local_mass_exclusion(259201, 4), "mean/junta/extension"),
    ("affine_parity_small_mass_spectrum", "proved", False,
     lambda: proof.first_layer_quota_certificate(37), "spectrum dependency"),
    ("affine_parity_small_mass_spectrum", "union_allowed_masses", [0, 2, 34, 36, 38],
     lambda: proof.first_layer_quota_certificate(37), "spectrum dependency"),
    ("local_mass_exclusion", "excluded", False,
     lambda: proof.eventual_first_layer_exclusion(524287), "first-layer closure"),
    ("first_layer_quota_certificate", "forced_low_masses", [4, 6],
     lambda: proof.eventual_first_layer_exclusion(524287), "first-layer closure"),
    ("generic_bounded_mean_certificate", "junta_coefficient", 1,
     lambda: proof.power_band_exclusion(6700417, 3, 3 * 6700417 + 2), "bounded-mean dependency"),
    ("generic_bounded_mean_certificate", "junta_size_strict_upper", "1",
     lambda: proof.power_band_exclusion(6700417, 3, 3 * 6700417 + 2), "extension/divisibility"),
])
def test_failed_child_premises_block_closure(monkeypatch, name, field, value, call, error):
    actual = getattr(proof, name)
    monkeypatch.setattr(proof, name, lambda *args: {**actual(*args), field: value})
    with pytest.raises(ArithmeticError, match=error):
        call()


@pytest.mark.parametrize("pipeline", ["local", "power_band"])
@pytest.mark.parametrize("name,field,value", [
    ("cube_height_bound_certificate", "height_over_mean_squared", 323),
    ("paired_cube_basis_certificate", "identity", "wrong"),
    ("sharp_integral_quadratic_lift_floor", "integer_quadratic_cube_mean_lattice", "(1/8) Z"),
    ("stabilizer_mass_certificate", "value", Fraction(0)),
    ("derivative_support_certificate", "uniform_strict_influence_lower", "1/16"),
    ("derivative_support_certificate", "total_influence_variance_factor", 0),
    ("middle_slice_mean_certificate", "identity", "wrong"),
    ("middle_slice_mean_certificate", "full_cube_mean_lattice", "(1/8) Z"),
])
def test_shared_dependency_corruption_blocks_both_theorems(
        monkeypatch, pipeline, name, field, value):
    actual = getattr(proof, name)
    monkeypatch.setattr(proof, name, lambda *args: {**actual(*args), field: value})
    with pytest.raises(ArithmeticError):
        if pipeline == "local":
            proof.local_mass_exclusion(524287, 4)
        else:
            proof.power_band_exclusion(6700417, 3, 3 * 6700417 + 2)


def test_weakened_quotient_floor_is_not_discarded(monkeypatch):
    actual = proof.quotient_floor
    monkeypatch.setattr(proof, "quotient_floor", lambda p, u: 1 if u == 2 else actual(p, u))
    with pytest.raises(ArithmeticError, match="quotient floor dependency"):
        proof.first_layer_quota_certificate(37)


def test_saved_artifact_and_no_global_or_smaller_prime_overclaim():
    result = proof.proposition_15775()
    assert result == json.loads((ROOT / "evidence/e1_gmin_m4_prop15775.json").read_text())
    assert result["proved"] and result["minimum_prime"] == 259201
    assert result["status"] == "PROVED_INFINITE_FAMILY"
    assert result["records_are_identity_replays_not_a_prime_census"]
    assert len(result["records"]) == 2 and len(result["superlinear_support_records"]) == 3
    for name in ("smaller_prime_frontiers_changed", "new_equality_catalog_used",
                 "residual_ii_closed_general", "minimal_four_gap_bridge_closed_general",
                 "eventual_E1_proved", "e1_closed_general", "original_MO_limit_closed"):
        assert result[name] is False
    assert (ROOT / result["proof_note"]).is_file()
