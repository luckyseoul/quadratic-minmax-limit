"""Tests for Prop 15.479 — complementary Δ² by balanced full line."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15477 import delta2_tr_named
from e1_gmin_m4_prop15479 import main, prove_A, prove_B, prove_open


def test_comp_law():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["n_bal"] == 0
    assert A["rows"]["5"]["err"] < 1e-4
    assert A["rows"]["7"]["n_bal"] == 1176
    assert A["rows"]["7"]["n_nfull1"] == 3528
    assert A["rows"]["7"]["n_nfull1"] > A["rows"]["7"]["n_bal"]
    assert A["rows"]["7"]["err"] < 1e-4


def test_nfull_is_not_the_name():
    A = prove_A()
    # the fail-when-wrong: n_full=1 is strictly larger than the zero class
    assert A["rows"]["7"]["n_nfull1"] - A["rows"]["7"]["n_bal"] == 2352


def test_p5_no_balanced():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_bal5"] == 0
    assert delta2_tr_named(5) == 4500


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.479"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["comp_delta2_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
