from e1_gmin_m4_prop15659 import (
    p7_negative_infinity_catalog_classification,
    p7_negative_infinity_certificate,
    p7_negative_infinity_floor_rigidity,
    theorem_p7_negative_infinity_size_six_exclusion,
)


def test_negative_infinity_floor_rigidity():
    out = p7_negative_infinity_floor_rigidity()
    assert out["proved"] is True
    assert out["scaled_floors"] == {1: 6, 3: 14, 5: 6}
    assert out["same_type_congruence_modulus"] == 8
    assert out["conclusion"] == "exactly one mean-14 direction per type"


def test_negative_infinity_odd_catalog_counts():
    out = p7_negative_infinity_catalog_classification()
    assert out["proved"] is True
    assert out["exact_johnson_catalog_counts"] == {
        "b1_phase1_mean6": 1,
        "b5_phase1_mean6": 1,
        "b1_phase1_mean14": 1764,
        "b5_phase1_mean14": 1764,
        "b3_phase1_mean14": 36,
    }


def test_complete_negative_infinity_certificate_has_no_survivor():
    out = p7_negative_infinity_certificate()
    assert out["floor_survivors"] + out["floor_rejected"] == out[
        "finite_boundaries"
    ]
    assert out["stabilizer_orbits"] == 1_750
    assert out["affine_span_rejected_cases"] + out["exact_catalog_cases"] == out[
        "elevation_cases"
    ]
    assert out["checked_exact_catalog_pairs"] == 25 * 36 * 36
    assert out["surviving_catalog_pairs"] == 0
    assert out["complete_branch_mod7_infeasible"] is True


def test_prop15659_closes_only_the_claimed_branch():
    out = theorem_p7_negative_infinity_size_six_exclusion()
    assert out["proved"] is True
    assert out["p7_cH_positive_infinity_plus_five_finite"] == "CLOSED_BY_15.658"
    assert out["p7_cH_negative_infinity_plus_five_finite"] == "CLOSED"
    assert out["p7_six_finite"] == "OPEN"
    assert out["closes_all_p7_size_six"] is False
    assert out["closes_residual_ii"] is False
    assert out["closes_R1"] is False
    assert out["L_status"] == "OPEN"
