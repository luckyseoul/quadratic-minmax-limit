"""Prop. 15.725: exact parabola-plus-internal family closure."""
from __future__ import annotations

from e1_gmin_m4_prop15725 import (
    FINITE_PRIMES,
    exact_case,
    finite_exact_certificate,
    large_prime_symbolic_branch,
    normalization_to_reference,
    normalized_fibre_label_scale,
    normalized_target_direction,
    symbolic_phase_transport,
    theorem_parabola_internal_family,
)


EXPECTED_DIRECTION_HASHES = {
    17: "bd9a07d5bdaaa71eda5033b867084409c41e6f88fe6c74500476c21e492989c7",
    19: "6e97255e918278861bf16f9c601aa4e95dbff4e3912765c4dee2f47803c860ef",
    23: "090a1a57921d2996070ddbc2eb56e6ee5c4ef3cb33084d5c9911e8af5b382d28",
    29: "6e0b822bd9fa00bfbbd88dc6bd36741753d698c0a75229b9d06ec39f88ac45e5",
    31: "10d6ab05053239f8bb681c3a1c1eb4fefd9e589d18bd1c529b74320af48fc9d3",
    37: "eea4c0063e2e6f20cad8a8d3902a67f4d02d5127af7d3472d0fde2fddc6ab453",
    41: "6d100a78e3464cb120aa7f14759e340636287b3c6c9f4d18a7ef524cbab43f49",
    43: "d60475fc746d61ea48bdaaec69d6924eb71323a24633ba95f1969f8007ec4e11",
    47: "9f0fe360e3d9f7526ee18096487eeac1d8565a9d68906fd4fadc5824eefbff79",
}

EXPECTED_CASE_HASHES = {
    17: "4ab894599d3bbe6606c04a8739a5b552a8e546bb7e95300d100b62e2a5f3127d",
    19: "f3fb6417ef393ad3e22df971c1c80f373fefbad402f0cd83e89abc35e15ef59a",
    23: "32987fbd4a46fb1b372742deabca9e2b8cfbf3553f2d2165645696efad812465",
    29: "b09777ebb615a1a857e41abd5d1542c858f91f6c4e01fcee52db744b8ca504c6",
    31: "3630e5a00803c185454d0b41a689a01d504e147012d7e4378fab404f3135bf3c",
    37: "0f2ac903211aedde5f60302f3cf6f189b7a99e3f5ee7ee3f74c3afdf3717c570",
    41: "8c8166df43ffacd24c4c5318161bbdfb4b6c40f51584a916ad9dc2c419b5c3d5",
    43: "e2d90624578736745fab3c1ea07497d61473e2fc6ec749d8a4aba064581d517e",
    47: "4fbe680d7af8b89d2b4150baffd60967f646d9aa1a4ff889ab60255bd5507706",
}

EXPECTED_NORMALIZATION_HASHES = {
    17: "e547a3f362190936760fcca3930115b2dedfcb972e135cb9271a679734d8cc25",
    19: "b00b91689bf1c414407eea278ddd1f2a08babf41257dc8c3001eac40b9f51844",
    23: "5a81beeea784d33ed682d80589be426209c676cd69d850f46078a764e2b79ce5",
    29: "9bf7a2c432b29a10e098d043ac02022d60cdcd8f4469ad1f64c25cb07d4caa7c",
    31: "d7c311024ee86ad3890d9c68fba95a9ae98eb6045efe609013f0c04eb9d98c63",
    37: "252312f8ac48ee5f3e4709d7d14d79ea67f3e17a593065ab6887b38cda676fbe",
    41: "cd710a67cf1cea923e8975b12744a76bef743d444dc07d6a2668db6b4237c0fe",
    43: "8d5d97b507619cb83b185d819fbdb7855bd84177d0a38b893cd5d9516ddf1a5e",
    47: "b5cb88af3b92ee6ac4fa18a28e35aeee9107b8be7d3d7319648f8249150874ad",
}

EXPECTED_TYPE_SUM_HISTOGRAMS = {
    17: {"282,290": 8, "294,278": 8, "306,290": 48},
    19: {"380,338": 9, "380,356": 18, "380,362": 54},
    23: {"552,530": 121},
    29: {"870,842": 196},
    31: {"992,962": 225},
    37: {"1406,1370": 324},
    41: {"1722,1682": 400},
    43: {"1892,1850": 441},
    47: {"2256,2210": 529},
}


def test_representative_exact_case_uses_actual_types_and_phase_zero_floors():
    """A transparent coordinate check; the exhaustive test is separate below."""
    case = exact_case(17, 3, 5)
    rows = {tuple(row["direction"]): row for row in case["rows"]}

    exceptional = rows[(0, 1)]
    assert exceptional["type"] == 1
    assert exceptional["b"] == 1
    assert exceptional["floor"] == 18
    assert exceptional["occupancy"] == (
        0,
        4,
        0,
        0,
        2,
        0,
        4,
        0,
        0,
        1,
        2,
        2,
        0,
        0,
        0,
        0,
        2,
    )

    finite_direction = rows[(1, 4)]
    assert finite_direction["type"] == 1
    assert finite_direction["b"] == 3
    assert finite_direction["floor"] == 28
    assert finite_direction["occupancy"] == (
        2,
        1,
        4,
        0,
        0,
        0,
        2,
        0,
        1,
        0,
        0,
        0,
        3,
        2,
        2,
        0,
        0,
    )

    assert case["source_vertical_type"] == -1
    assert case["boundary_multiplier"] == -1
    assert case["transported_c_H"] == 1
    assert case["common_phase"] == 0
    assert case["type_counts"] == {-1: 9, 1: 9}
    assert case["type_sums"] == {-1: 294, 1: 278}
    assert case["type_budget"] == 162
    assert case["both_types_strict"] is True


def test_nonsquare_symmetry_includes_direction_and_fibre_label_normalization():
    normalization = normalization_to_reference(17, 5, 3)
    assert normalization == {"reference_nu": 3, "scale": 8, "reference_a": 5}
    assert normalized_target_direction(17, 8, (1, 4)) == (1, 9)
    assert normalized_fibre_label_scale(17, 8, (1, 4)) == 8
    assert normalized_target_direction(17, 8, (0, 1)) == (0, 1)
    assert normalized_fibre_label_scale(17, 8, (0, 1)) == 1


def test_symbolic_phase_and_unproved_p_at_least_53_threshold_arithmetic():
    phase = symbolic_phase_transport()
    assert phase["involution_identity"] == "h(iota(x))=(-a/x^2)h(x)"
    assert phase["fixed_point_equation"] == "x^2=-a"
    assert phase["nonzero_factor_product"] == 1
    assert phase["boundary_multiplier"] == "chi_p(-nu)"
    assert phase["transported_c_H"] == 1
    assert phase["residual_exponent"] == "(|H|-3)/2=2p-1 is odd"
    assert phase["sign_reduction"] == "eps*(-1)*1*eps*(-1)=1"
    assert phase["phase_independent_of_direction_type"] is True
    assert phase["common_phase"] == 0
    assert phase["proved"] is True

    branch = large_prime_symbolic_branch()
    assert branch["scope"] == "every prime p>=53"
    assert branch["separate_from_finite_enumeration"] is True
    assert branch["odd_integrality_consequence"] == (
        "3 <= b_c <= p-4 for every finite c"
    )
    assert branch["exceptional_direction"] == {
        "direction": "[0:1]",
        "reason": "B(x) is even in x, leaving only x=0 unpaired",
        "b": 1,
        "type": 1,
    }
    assert branch["threshold_squared_polynomial_values"] == {
        "c_nonzero_lower_gt_1": 117,
        "c_nonzero_upper_lt_p_minus_2": 1044,
        "c_zero_lower_gt_1": 1456,
        "c_zero_upper_lt_p_minus_2": 572,
    }
    assert branch["forward_differences_at_53"] == {
        "c_nonzero_lower_gt_1": 55,
        "c_nonzero_upper_lt_p_minus_2": 70,
        "c_zero_lower_gt_1": 81,
        "c_zero_upper_lt_p_minus_2": 53,
    }
    assert branch["type_minus_exact_lower_gap"] == "m(p-7)>0"
    assert branch["type_plus_exact_lower_gap"] == "(m-1)(p-7)>0"
    assert branch["threshold_arithmetic_verified"] is True
    assert branch["character_curve_bounds_status"] == "UNPROVED"
    assert branch["admissible_degenerate_locus"] == "4*a*nu+1=0"
    assert branch["conditional_only"] is True
    assert branch["proved"] is False


def test_all_nine_primes_are_exhausted_with_exact_typed_budget_certificates():
    certificate = finite_exact_certificate()
    assert certificate["scope"] == list(FINITE_PRIMES)
    assert certificate["separate_from_large_prime_proof"] is True
    assert certificate["all_nonsquare_nu_enumerated_directly"] is True
    assert certificate["all_a_with_chi_minus_a_negative_enumerated"] is True
    assert certificate["all_projective_directions_enumerated"] is True
    assert certificate["actual_norm_character_types_used"] is True
    assert certificate["common_phase_zero_checked_per_case"] is True
    assert certificate["parameter_case_count"] == 2381
    assert certificate["symmetry_reduced_case_count"] == 139
    assert certificate["direction_record_count"] == 92664
    assert certificate["survivor_count"] == 0

    for p in FINITE_PRIMES:
        row = certificate["rows"][str(p)]
        half = (p - 1) // 2
        cases = half * half
        directions = cases * (p + 1)
        budget = (p + 1) * (p + 1) // 2
        assert row["nonsquare_count"] == half
        assert row["admissible_a_count"] == half
        assert row["parameter_case_count"] == cases
        assert row["direction_record_count"] == directions
        assert row["phase_transport_check_count"] == cases
        assert row["normalization_direction_check_count"] == directions
        assert row["canonical_a_case_count"] == half
        assert row["type_budget"] == budget
        assert row["b1_directions_per_case"] == {"1": cases}
        assert row["all_cases_have_one_b1_direction"] is True
        assert row["both_types_strict_case_count"] == cases
        assert row["all_cases_excluded_in_both_types"] is True
        assert row["survivor_count"] == 0

        histogram = row["type_sum_pair_histogram"]
        assert histogram == EXPECTED_TYPE_SUM_HISTOGRAMS[p]
        assert sum(histogram.values()) == cases
        for pair in histogram:
            minus, plus = map(int, pair.split(","))
            assert minus > budget
            assert plus > budget

        assert row["direction_audit_sha256"] == EXPECTED_DIRECTION_HASHES[p]
        assert row["case_audit_sha256"] == EXPECTED_CASE_HASHES[p]
        assert (
            row["normalization_audit_sha256"]
            == EXPECTED_NORMALIZATION_HASHES[p]
        )
        assert row["proved"] is True

    assert certificate["proved"] is True


def test_theorem_honestly_retracts_the_explicit_family_close():
    theorem = theorem_parabola_internal_family()
    assert theorem["proved"] is False
    assert "unconditional family exclusion remains open" in theorem["scope"]
    assert theorem["theorem"]["parabola_plus_internal_family"] == "OPEN"
    assert theorem["opposite_product_sign_checked"] is False
    assert theorem["retraction"] == {
        "all_prime_character_bounds": "UNPROVED",
        "opposite_product_sign": "UNCHECKED",
        "finite_phase_zero_census": "EXACT",
    }
    assert theorem["theorem"]["whole_p_plus_one_shell"] == "OPEN"
    assert theorem["theorem"]["residual_ii"] is False
    assert theorem["theorem"]["type_I"] is False
    assert theorem["theorem"]["limit_exists"] is False
    assert theorem["L_status"] == "OPEN"
