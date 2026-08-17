"""Tests for Prop 15.352 — complete p=7 type dictionary."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15312 import first_nc
from e1_gmin_m4_prop15351 import sig_y, y_T
from e1_gmin_m4_prop15352 import (
    SIG_294,
    SIG_588,
    SIG_LAST,
    SIG_TNC,
    first_nc_c,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_seven_types_sum_5726_and_Tnc2_is_294():
    A = prove_A()
    assert A["proved"] is True
    assert A["sum"] == 5726
    assert A["rows"]["294_Tnc2"]["census"] == 294
    assert A["rows"]["588_ystar"]["census"] == 588
    assert A["rows"]["1176_c1"]["census"] == 1176


def test_c1_is_leftover_c2_is_Tnc():
    B = prove_B()
    assert B["proved"] is True
    y0 = y_T(7)
    raw = first_nc_c(7, y0, 1)
    c1 = first_nc_c(7, y0, 1, skip_sigs=(SIG_588, SIG_294))
    c2 = first_nc_c(7, y0, 2)
    assert sig_y(7, raw["y"]) == SIG_588
    assert sig_y(7, c1["y"]) == SIG_LAST
    assert sig_y(7, c2["y"]) == SIG_TNC
    assert sig_y(7, c1["y"]) != sig_y(7, first_nc(7, y0)["y"])


def test_Tnc2_signature_is_294_not_1176():
    y0 = y_T(7)
    tnc = first_nc(7, y0)
    tnc2 = first_nc(
        7, tnc["y"], skip_pred=lambda y: np_sign_eq(y, y0)
    )
    assert sig_y(7, tnc2["y"]) == SIG_294


def np_sign_eq(y, y0):
    import numpy as np

    return np.array_equal(np.sign(y), np.sign(y0))


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["p7_seven_types_named"] is True
    assert out["proved"]["Tnc_c1_is_leftover_1176"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
