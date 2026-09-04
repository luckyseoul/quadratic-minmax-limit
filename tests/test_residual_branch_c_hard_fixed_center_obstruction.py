"""Checks for the fixed complement-pair centre obstruction."""

from scripts.residual_branch_c_hard_fixed_center_obstruction import (
    fixed_family_center_obstruction_certificate,
)


def test_fixed_complement_pair_family_has_no_required_collision_direction() -> None:
    row = fixed_family_center_obstruction_certificate()
    assert row["exact_profile_replay"] is True
    assert row["aggregate_signature_hex"] == "00800005"
    assert row["correction_signature_support"] == (2, 23)
    assert row["pair_center_cases_checked"] == 108_000
    assert row["intersection_size_histogram"] == {0: 107_700, 1: 300}
    assert row["shared_orbit_incidence_count"] == 300
    assert row["shared_orbit_spatial_direction_histogram"] == {
        4: 60,
        9: 60,
        20: 60,
        21: 60,
        29: 60,
    }
    assert row["required_direction_shared_orbit_incidence_count"] == 0
    assert row["fixed_16_half_family_physically_excluded"] is True
