"""Tests for Prop 15.295 — unsplit pairing; Paley×norm not a scheme."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15295 import (
    main,
    paley_norm_scheme_spread,
    prove_E_XA_zero,
    prove_cov_is_S,
    prove_not_scheme,
    prove_open,
)


def test_paley_norm_not_scheme():
    A = prove_not_scheme([5, 7, 11])
    assert A["proved"] is True
    assert paley_norm_scheme_spread(5) > 0
    assert A["sample"]["5"]["scheme"] is False


def test_E_XA_zero_and_I_even():
    B = prove_E_XA_zero([5, 7, 11])
    assert B["proved"] is True
    assert B["invert_hit"] < B["invert_tot"]
    assert B["sample"]["5"]["I_even"] is True
    assert B["sample"]["5"]["triv_sum_sq"] != 0


def test_cov_is_S_p5():
    C = prove_cov_is_S([(5, 2), (5, 16)])
    assert C["proved"] is True
    assert C["sample"]["5:2"]["spectrum_hit"] is True
    assert C["sample"]["5:2"]["E_XA_zero"] is True
    assert C["invert_hit"] < C["invert_tot"]


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["paley_norm_not_scheme"] is True
    assert out["proved"]["E_XA_zero"] is True
    assert out["proved"]["E_u_XA_is_S"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
