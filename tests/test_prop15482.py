"""Tests for Prop 15.482 — W_c real Φ(Tr,N); mean-field pairing ≡4."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15482 import (
    Nrm,
    beta_den_candidate,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    two_sq,
    wick_Q,
)


def test_Wc_real_Phi():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["imag_max"] < 1e-8
    assert A["rows"]["5"]["collisions"] == 0
    assert A["rows"]["5"]["homog_err"] > 9.0


def test_mean_field_is_4():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["rows"]["5"]["wick"] - 4.0) < 1e-8
    assert abs(B["rows"]["5"]["wick"] - float(LIVE_QPP[5])) > 0.05
    assert abs(wick_Q(5) - 4.0) < 1e-8


def test_48_not_norm():
    C = prove_C()
    assert C["proved"] is True
    assert (2, 3) in two_sq(13)
    assert (3, 20) in two_sq(409)
    assert two_sq(48) == []


def test_den_candidate_dies_at_p3():
    D = prove_D()
    assert D["proved"] is True
    assert Nrm(beta_den_candidate(5)) == 13
    assert Nrm(beta_den_candidate(7)) == 409
    assert Nrm(beta_den_candidate(3)) == 153
    assert D["Hplus3_over_2p"] == 1


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.482"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Wc_real_Phi_TrN"] is True
    assert out["proved"]["mean_field_pairing_is_4"] is True
    assert out["proved"]["48_not_gaussian_norm"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
