"""Tests for Prop 15.440 — Lemma D plane is Geo, not A5."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15276 import (
    lemma_D_2plane_amplitudes_proved,
    locked_phases,
    theorem_amplitudes_not_parallel,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15440 import (
    cleared_numerator,
    geo_amplitude_vector,
    main,
    plane_identity_proved,
    prove_B,
    prove_C,
    prove_open,
)


def test_rational_and_geo_lock():
    assert cleared_numerator(2, 3) == 0
    assert (2 - 3) - (2 - 1) != 0
    s = locked_phases(0, 1, 5)
    g = geo_amplitude_vector(5, 1, 2, *s, 1)
    assert abs(g.sum()) < 1e-8
    s_bad = (s[0], s[1], (s[2] + 2) % 5)
    g_bad = geo_amplitude_vector(5, 1, 2, *s_bad, 1)
    assert abs(abs(g_bad.sum()) - 100.0) < 1e-6
    B = prove_B()
    assert B["proved"] is True
    assert B["lock_drop_sum_abs"] > 1e-6


def test_A5_is_not_the_plane():
    C = prove_C()
    assert C["proved"] is True
    assert C["A5_sum_abs"]["5"] > 10.0
    assert plane_identity_proved() is True
    par = theorem_amplitudes_not_parallel()
    assert par["proved"] is True
    assert "prop_15_440_geo_plane" in par["depends_on"]
    assert lemma_D_2plane_amplitudes_proved() is True


def test_flags_untouched():
    O = prove_open()
    assert O["proved"] is False
    assert O["plane_identity"] is True
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.440"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["rational_identity"] is True
    assert out["proved"]["geo_on_plane"] is True
    assert out["proved"]["A5_free_not_plane"] is True
    assert out["proved"]["geo_span_plane"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
