"""Prop 15.589 — PSL decomposition and exceptional-scalar reduction."""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import e1_gmin_m4_prop15589 as M
from e1_gmin_m4_prop15534 import field_tables
from minmax_quadratic import paley_conference_prime_power


@pytest.mark.parametrize("p,r", [(5, 2), (7, 5), (11, 14), (13, 20), (19, 44)])
def test_character_decomposition_dimensions(p, r):
    assert M.n_principal_constituents(p) == r
    assert M.d_of(p) + r * M.n_of(p) == M.dim_Z(p)
    assert 1 + 2 * r == M.dim_F(p)


def test_character_decomposition_theorem_and_gap_audit():
    A = M.theorem_A_character_decomposition()
    assert A["proved"], A
    assert A["gap_audit"]["25"]["principal"] == 2
    assert A["gap_audit"]["49"]["principal"] == 5
    assert A["gap_audit"]["121"]["principal"] == 14
    assert all(row["multiplicity_free"] for row in A["by_p"].values())


def test_Z_decomposition_excludes_other_families():
    B = M.theorem_B_Z_decomposition()
    assert B["proved"], B
    assert B["multiplicity_free"]
    assert B["no_trivial"] and B["no_steinberg"] and B["no_cuspidal"]


def test_phi_has_exactly_one_possible_small_block():
    C = M.theorem_C_phi_multiplicity_reduction()
    assert C["proved"], C
    assert C["exact_remaining_scalar"] == "lambda_exc >= 6"
    assert C["mult_lambda_min_ge_n_proved_unconditionally"] is False


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_variance_room_halves_for_exceptional_block(p):
    assert M.variance_room_exceptional(p) * 2 == M.variance_room_principal(p)
    assert M.delta2_room_exceptional(p) * 2 == M.delta2_room_principal(p)
    n, D = M.n_of(p), M.dim_Z(p)
    gap = M.spectral_mean(p) - 6
    assert Fraction(n, D) * gap * gap == M.variance_room_principal(p)
    assert Fraction(M.d_of(p), D) * gap * gap == M.variance_room_exceptional(p)


def test_FWW_wrong_principal_count_breaks_dimension():
    for p in (5, 7, 11):
        r = M.n_principal_constituents(p)
        assert M.d_of(p) + (r + 1) * M.n_of(p) != M.dim_Z(p)
        assert M.d_of(p) + (r - 1) * M.n_of(p) != M.dim_Z(p)


def test_FWW_wrong_U_fixed_dimensions_break_F():
    for p in (5, 7, 11):
        r = M.n_principal_constituents(p)
        assert 1 + 2 * r == M.dim_F(p)
        assert 1 + r != M.dim_F(p)


def test_exceptional_scalar_is_quartic_variance():
    E = M.theorem_E_exceptional_quartic_variance()
    assert E["proved_reduction"] and E["proved_census"], E
    assert E["proved_general_inequality"] is False
    assert M.lambda_exc_from_quartic_variance(5, Fraction(3300, 13)) == Fraction(176, 13)
    assert M.lambda_exc_from_quartic_variance(7, Fraction(317520, 409)) == Fraction(4320, 409)


@pytest.mark.parametrize("p,total", [(7, 84), (11, 330), (19, 1710)])
def test_p3mod4_profile_energy_conservation_constant(p, total):
    assert M.profile_energy_total(p) == total
    # Dropping the factor 1/4 is the natural pair-counting error.
    assert p * (p * p - 1) != total


@pytest.mark.parametrize("p", [7, 11, 19, 23, 31, 43])
def test_profile_energy_2p_divisibility_normalization_and_parity(p):
    T = M.normalized_profile_energy_total(p)
    assert M.profile_energy_total(p) == 2 * p * T
    assert M.normalized_quartic_variance_threshold(p) == Fraction(3 * T, 8)
    assert M.quartic_variance_floor_threshold(p) == (
        4 * p * p * M.normalized_quartic_variance_threshold(p)
    )
    assert M.quartic_pointwise_parity(p) == (1 if p % 8 == 3 else 0)
    assert Fraction(M.quartic_pointwise_parity(p)) < (
        M.normalized_quartic_variance_threshold(p)
    )


def test_profile_energy_arithmetic_theorem():
    FA = M.theorem_F_profile_energy_arithmetic()
    assert FA["proved_energy_divisibility_2p"], FA
    assert FA["proved_pointwise_parity"], FA
    assert all(row["ok"] for row in FA["by_p"].values())


def test_low_strata_exact_values_match_censuses():
    assert M.k1_quartic_variance(5) == 500
    assert M.k1_quartic_variance(7) == 7056
    assert M.k1_quartic_variance(11) == 108900
    assert M.k3_quartic_variance_p3mod4(7) == 784
    assert M.k3_quartic_variance_p3mod4(11) == 21780


@pytest.mark.parametrize("p", [7, 11, 19, 23, 31])
def test_p3mod4_k1_k3_clear_QVAR(p):
    threshold = M.quartic_variance_floor_threshold(p)
    assert M.k1_quartic_variance(p) >= threshold
    assert M.k3_quartic_variance_p3mod4(p) >= threshold


@pytest.mark.parametrize("p", [13, 17, 29, 37])
def test_p1mod4_euler_product_k3_bound_clears_QVAR(p):
    lower = M.k3_quartic_variance_lower_p1mod4(p)
    assert lower >= M.quartic_variance_floor_threshold(p)


def test_p1mod4_euler_constant_is_load_bearing_at_p13():
    threshold = M.quartic_variance_floor_threshold(13)
    assert Fraction(13**6, 30**2) >= threshold
    assert Fraction(13**6, 31**2) < threshold


def test_profile_energy_reduction_closes_exactly_low_strata():
    FG = M.theorem_FG_profile_energy_and_low_strata()
    assert FG["proved_profile_energy_identity_p3mod4"]
    assert FG["proved_profile_energy_conservation_p3mod4"]
    assert FG["proved_k1_k3_QVAR_all_primes"]
    assert FG["remaining_exceptional_strata"] == "k>=4"
    assert all(row["k1_clears"] and row["k3_clears"] for row in FG["by_p"].values())


@pytest.mark.parametrize(
    "p,value",
    [
        (5, Fraction(130)),
        (7, Fraction(4900, 9)),
        (11, Fraction(73810, 21)),
    ],
)
def test_odd_coset_spherical_benchmark_exact_values(p, value):
    assert M.spherical_quartic_variance(p) == value
    assert M.spherical_QVAR_gap(p) > 0


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_spherical_gap_closed_form_and_shorter_ordinary_shell(p):
    q = p * p
    assert M.spherical_QVAR_gap(p) == Fraction(
        q * (q - 1) * (q - 11), 16 * (q + 5)
    )
    assert p + 1 < M.n_of(p)


def test_odd_coset_reduction_does_not_claim_the_harmonic_bound():
    H = M.theorem_H_odd_coset_spherical_benchmark()
    assert H["proved_reduction"]
    assert H["maxplus_is_odd_coset_first_shell"]
    assert H["maxplus_is_ordinary_lattice_first_shell"] is False
    assert H["ordinary_minimum_shell_design_route_applies"] is False
    assert H["sufficient_harmonic_target"].endswith(">= 0")


@pytest.mark.parametrize("p", [7, 11, 19, 23, 31, 43])
def test_coarse_profile_counterexample_has_all_claimed_constraints(p):
    rec = M.coarse_profile_counterexample(p)
    assert rec["m"] == (p + 1) // 2
    assert rec["sum_a"] == M.profile_energy_total(p)
    assert rec["full_support"]
    assert rec["all_energies_divisible_by_2p"]
    assert rec["all_energies_individually_profile_admissible"]
    assert rec["equal_directional_means_under_cyclic_orbit"]
    assert rec["signed_energy_magnitude"] == 2 * p * rec["parity"]
    assert rec["cyclic_orbit_variance"] < rec["QVAR_threshold"]
    assert rec["violates_QVAR"]

    t = (p - 3) // 4
    assert set(rec["b"]) == {t, t + 1}
    assert rec["b"].count(t) == rec["b"].count(t + 1) == (p + 1) // 4
    for b, h in rec["line_profile_witnesses"].items():
        assert len(h) == p
        assert sum(h) == 0
        assert sum(value * value for value in h) == 2 * p * b
        assert min(h) >= -(p + 1) // 2
        assert max(h) <= (p - 1) // 2
        sigma = [1 + 2 * value for value in h]
        assert sum(sigma) == p
        assert all(-p <= value <= p and value % 2 for value in sigma)
        assert rec["line_profile_degrees_mod_p"][b] <= (
            rec["line_profile_degree_bound"]
        )


def test_coarse_profile_countermechanism_kills_only_coarse_route():
    I = M.theorem_I_coarse_profile_constraints_insufficient()
    assert I["proved_countermechanism"], I
    assert "Boolean ridge reconstruction" in I["missing_kind_of_input"]
    assert "coefficient kernels" in I["missing_kind_of_input"]
    assert all(row["ok"] for row in I["by_p"].values())
    with pytest.raises(ValueError):
        M.coarse_profile_counterexample(13)


def test_p11_k4_requires_active_subset_mixing():
    J = M.theorem_J_p11_k4_active_subset_mixing()
    assert J["proved_counterexample"], J
    assert J["normalized_QVAR_threshold"] == "45/8"
    assert J["balanced"]["E_B2"] == "5"
    assert J["balanced"]["fails_QVAR"]
    assert J["unbalanced"]["E_B2"] == "63"
    assert J["unbalanced"]["clears_QVAR"]
    assert J["n_pure_reps"] == 480
    assert J["n_full_k4_vectors"] == 58_080
    assert J["aggregate_E_B2"] == "39/2"
    assert J["aggregate_E_Z2"] == "9438"
    assert J["aggregate_clears_QVAR"]

    data = json.loads(
        (Path(__file__).parents[1] / "evidence" / "maxplus_p11" /
         "k4_active_subset_quartic_p11.json").read_text()
    )
    assert data["n_pure_reps"] == J["n_pure_reps"]
    assert data["n_full_vectors"] == J["n_full_k4_vectors"]
    assert data["aggregate_E_B2"] == J["aggregate_E_B2"]
    assert {row["E_B2"] for row in data["subsets"].values()} == {"5", "63"}


def test_full_support_requires_top_degree_mixing():
    K = M.theorem_K_full_support_top_degree_mixing()
    assert K["proved_counterexample"], K
    assert K["p7"]["top_zero_count"] == 0
    assert K["p7"]["E_B2_per_nonzero_class"] == "44/15"
    assert K["p11"]["top_zero_E_B2"] == "137/36"
    assert K["p11"]["top_zero_fails_QVAR"]
    assert K["p11"]["E_B2_per_nonzero_class"] == "111483/14039"
    assert K["p11"]["nonzero_classes_clear_QVAR"]
    assert K["p11"]["degree_drops_twice_count"] == 0
    assert K["p11"]["aggregate_E_B2"] == "114771/14903"
    assert K["p11"]["aggregate_clears_QVAR"]
    assert [
        (row["n_projective_classes"], row["vectors_per_class"], row["E_B2_per_class"])
        for row in K["p11"]["degree3_projective_orbits"]
    ] == [(6, 123_420, "151/51"), (6, 225_060, "397/93")]

    data = json.loads(
        (Path(__file__).parents[1] / "evidence" / "maxplus_p11" /
         "full_support_top_degree_p7_p11.json").read_text()
    )
    assert data["p7"]["full_support_count"] == K["p7"]["full_support_count"]
    assert data["p11"]["full_support_count"] == K["p11"]["full_support_count"]
    assert data["p11"]["top_zero_count"] == K["p11"]["top_zero_count"]
    assert data["p11"]["aggregate_E_B2"] == K["p11"]["aggregate_E_B2"]
    projective = list(data["p11"]["degree_drop_projective_classes"].values())
    assert Counter((row["count"], row["E_B2"]) for row in projective) == Counter(
        {(123_420, "151/51"): 6, (225_060, "397/93"): 6}
    )


@pytest.mark.parametrize(
    "p,minimum",
    [(7, 1), (11, 3), (19, 10), (23, 16), (31, 30),
     (41, 54), (43, 60), (47, 74), (53, 96), (59, 119), (61, 122)],
)
def test_exact_minimum_quadratic_profile_energy(p, minimum):
    assert M.quadratic_profile_min_b(p) == minimum


def test_quadratic_energy_barrier_makes_k4_empty_from_p41():
    L = M.theorem_L_k4_empty_for_p_ge_41()
    assert L["proved"], L
    assert L["scope"] == "every odd prime p>=41"
    assert L["analytic_range"]["first_possible_prime"] == 67
    assert all(
        row["ok"] and row["four_profiles_exceed_total"]
        for row in L["finite_exact_range"].values()
    )
    assert 4 * M.quadratic_profile_min_b(31) == (31 * 31 - 1) // 8
    assert 4 * M.quadratic_profile_min_b(37) < (37 * 37 - 1) // 8
    with pytest.raises(ValueError):
        M.quadratic_profile_min_b(5)


def test_general_weil_energy_barrier_excludes_low_activity():
    theorem = M.theorem_M_general_low_activity_exclusion()
    assert theorem["proved"], theorem
    assert theorem["empty_condition"] == "k>=4 and p>4k^2"
    assert not M.weil_activity_barrier_excludes(101, 3)
    assert M.weil_activity_barrier_excludes(101, 4)
    assert M.weil_activity_barrier_excludes(101, 5)
    assert not M.weil_activity_barrier_excludes(101, 6)
    assert M.weil_activity_barrier_excludes(1009, 15)
    assert not M.weil_activity_barrier_excludes(1009, 16)


def test_k4_is_completely_closed_for_p3mod4():
    N = M.theorem_N_k4_closed_p3mod4()
    assert N["proved"], N
    assert N["nonempty_primes"] == [7, 11]
    assert N["nonempty_QVAR_moments"] == {"7": "44/15", "11": "39/2"}
    assert all(
        record["k4_empty"] and record["total_coefficient_candidates"] == 0
        for record in N["finite_coefficient_sieve"].values()
    )

    data = json.loads(
        (Path(__file__).parents[1] / "evidence" /
         "k4_p3mod4_coefficient_sieve.json").read_text()
    )
    assert data["19"]["coefficient_candidate_histogram"] == {"0": 210}
    assert data["23"]["coefficient_candidate_histogram"] == {"0": 495}
    assert data["31"]["coefficient_candidate_histogram"] == {"0": 1820}
    assert all(record["k4_empty"] for record in data.values())


def test_k4_QVAR_is_closed_for_every_prime():
    O = M.theorem_O_k4_QVAR_all_primes()
    assert O["proved"], O
    assert O["remaining_exceptional_strata"] == "k>=5"
    assert O["p1mod4"]["13"]["E_abs_Zpsi_sq"] == "8788"
    assert O["p1mod4"]["17"]["E_abs_Zpsi_sq"] == "314432/3"
    assert O["p1mod4"]["29"]["coefficient_candidates"] == 0
    assert O["p1mod4"]["37"]["coefficient_candidates"] == 0

    data = json.loads(
        (Path(__file__).parents[1] / "evidence" / "k4_p1mod4_closure.json").read_text()
    )
    assert data["empty"]["29"]["coefficient_candidate_histogram"] == {"0": 1365}
    assert data["empty"]["37"]["coefficient_candidate_histogram"] == {"0": 3876}
    assert data["nonempty"]["13"]["quartic"]["E_abs_Zpsi_sq"] == "8788"
    assert data["nonempty"]["17"]["quartic"]["E_abs_Zpsi_sq"] == "314432/3"
    assert data["nonempty"]["17"]["quartic"]["abs_Zpsi_sq_histogram"] == {
        "0": 41_616,
        "314432": 20_808,
    }


@pytest.mark.parametrize("p", [5, 7])
def test_square_affine_line_is_a_shorter_ordinary_lattice_vector(p):
    C = paley_conference_prime_power(p)
    r = np.zeros(p * p + 1)
    r[0] = 1
    r[1 : p + 1] = 1  # infinity plus the square-direction line F_p
    assert np.array_equal(C @ r, p * r)
    assert r @ r == p + 1 < p * p + 1
    wrong = r.copy()
    wrong[0] = 0
    assert not np.array_equal(C @ wrong, p * wrong)


def test_spherical_benchmark_trace_and_hilbert_schmidt_inputs_at_p5():
    p = 5
    q = p * p
    F = field_tables(p)
    mul, sub = F["mul"], F["sub"]

    generator = None
    for g in range(2, q):
        x, seen = 1, set()
        for _ in range(q - 1):
            seen.add(x)
            x = int(mul[x, g])
        if len(seen) == q - 1:
            generator = g
            break
    assert generator is not None

    psi = np.zeros(q, dtype=np.complex128)
    x = 1
    for exponent in range(q - 1):
        psi[x] = (1j) ** (exponent % 4)
        x = int(mul[x, generator])
    K = psi[sub]
    np.fill_diagonal(K, 0)
    K_ext = np.zeros((q + 1, q + 1), dtype=np.complex128)
    K_ext[1:, 1:] = K
    C = paley_conference_prime_power(p)
    P = (np.eye(q + 1) + C / p) / 2
    A = P @ K_ext @ P / 4

    assert abs(np.trace(A)) < 1e-10
    assert np.vdot(A, A).real == pytest.approx(q * (q - 1) / 32)


@pytest.mark.parametrize("p", [5, 7, 11, 13])
def test_quartic_variance_threshold_is_exactly_lambda_six(p):
    threshold = M.quartic_variance_floor_threshold(p)
    assert M.lambda_exc_from_quartic_variance(p, threshold) == 6
    assert M.lambda_exc_from_quartic_variance(p, threshold - 1) < 6


def test_floor_flag_remains_open():
    assert M.leftover_flags_unchanged()
