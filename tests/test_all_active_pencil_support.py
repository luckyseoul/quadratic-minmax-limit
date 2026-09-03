from __future__ import annotations

import pytest

from e1_gmin_m4_all_active_pencil_support import (
    all_active_pencil_equality_exclusion,
)


def test_p31_all_active_equality_is_excluded_symbolically():
    out = all_active_pencil_equality_exclusion(31)
    assert out["oriented_other_endpoint_set_size"] == 29
    assert out["completed_point_count"] == 31
    assert out["hard_phase_sign_changes_support"] is False
    assert out["physical_orbit_orientation_does_not_change_support"] is True
    assert out["hard_undetermined_direction_count"] == 16
    assert out["maximum_completed_direction_count"] == 16
    assert out["redei_megyesi_noncollinear_direction_minimum"] == 17
    assert out["full_ray_opposite_parallel_quota_maximum"] == 16
    assert out["parallel_quota_deficit"] == 13
    assert out["equality_c_eq_p_minus_2_excluded"] is True
    assert out["conclusion"] == "c>=p-1 when every hard center is nonzero"
    assert out["residual_ii_closed"] is False
    assert out["proved"] is True


def test_scope_rejects_nonprime_small_and_wrong_congruence_inputs():
    for p in (7, 27, 29, 35):
        with pytest.raises(ValueError):
            all_active_pencil_equality_exclusion(p)
