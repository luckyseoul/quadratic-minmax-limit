"""Tests for Prop 15.417 — index ID dead; 15.298 slack too wide."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15308 import G_order
from e1_gmin_m4_prop15331 import HPLUS, c_from_stabs, nl_stabs
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15417 import (
    dQ_dL_plus,
    index_hits_den,
    main,
    prove_open,
    slack_width,
    window_width,
)


def test_index_identification_dead():
    assert LIVE_QNL[5].denominator == 15
    assert LIVE_QNL[7].denominator == 171
    assert HPLUS[5] % 15 != 0
    assert G_order(7) % 19 != 0
    assert HPLUS[7] % 19 != 0
    assert G_order(7) % 171 != 0
    assert c_from_stabs(5, nl_stabs(5)) == 1
    assert c_from_stabs(7, nl_stabs(7)) == 19
    assert c_from_stabs(7, nl_stabs(7)) != 171
    assert index_hits_den(5) is False
    assert index_hits_den(7) is False


def test_slack_exceeds_window():
    assert dQ_dL_plus(5) != 0
    assert slack_width(5) > window_width(5)
    assert slack_width(7) > window_width(7)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["index_identification_dead"] is True
    assert out["proved"]["slack_too_wide"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["rows"]["5"]["c"] == 1
    assert out["algebra"]["A"]["rows"]["7"]["c"] == 19
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.417"
