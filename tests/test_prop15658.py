from e1_gmin_m4_prop15658 import (
    p7_positive_infinity_certificate,
    p7_positive_infinity_slack_classification,
    same_type_scaled_slack_congruence,
    theorem_p7_positive_infinity_size_six_exclusion,
)


def test_same_type_scaled_slack_congruence():
    for p in (5, 7, 11, 13):
        out = same_type_scaled_slack_congruence(p)
        assert out["proved"] is True
        assert out["modulus"] == p + 1


def test_p7_positive_infinity_slacks_are_unique_mean_eight():
    out = p7_positive_infinity_slack_classification()
    assert out["proved"] is True
    assert out["scaled_floors"] == {1: 8, 3: 8, 5: 8}
    assert out["all_scaled_means"] == 8


def test_complete_mod7_certificate_has_no_survivor():
    out = p7_positive_infinity_certificate()
    assert out["finite_boundaries"] == 1_906_884
    assert out["v100_checked"] == out["finite_boundaries"]
    assert out["nuka_cpu_checked"] == out["finite_boundaries"]
    assert out["v100_survivors"] == out["nuka_cpu_survivors"] == 0
    assert sum(out["matching_direction_mask_histogram"].values()) == 8 * out[
        "finite_boundaries"
    ]


def test_prop15658_closes_only_the_claimed_branch():
    out = theorem_p7_positive_infinity_size_six_exclusion()
    assert out["proved"] is True
    assert out["p7_cH_positive_infinity_plus_five_finite"] == "CLOSED"
    assert out["p7_cH_negative_infinity_plus_five_finite"] == "OPEN"
    assert out["p7_six_finite"] == "OPEN"
    assert out["closes_all_p7_size_six"] is False
    assert out["closes_residual_ii"] is False
    assert out["closes_R1"] is False
    assert out["L_status"] == "OPEN"
