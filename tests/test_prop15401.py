"""Tests for Prop 15.401 — mix-aware λ by D-type; NL means can be <6."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15400 import one_d_weight_lb, s1d_named
from e1_gmin_m4_prop15401 import (
    main,
    nl_type_means_ge_6,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_D_type_splits_at_p7():
    A = prove_A()
    assert A["proved"] is True
    assert A["p7_n1176_split"] == 4
    assert A["p7_n294_const0"] == 1
    # p=5 stays two-level
    ns = sorted(r["n"] for r in A["p5"])
    assert ns == [30, 100]


def test_NL_type_means_can_be_below_6():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_mean_8_3"] == 2
    assert B["rows_at_8_3"] == 2352
    assert B["nl_type_means_ge_6"] is False
    assert nl_type_means_ge_6() is False
    assert Fraction(B["mix"]) == Fraction(3360, 409)


def test_one_d_not_constant():
    C = prove_C()
    assert C["proved"] is True
    assert C["n0_AP"] == 56
    assert C["nunique"] == 4
    assert abs(C["mean"] - float(s1d_named(7))) < 1e-6
    assert C["n_1d"] == n_1d(7) == 140


def test_one_d_lb_still_misses():
    D = prove_open()
    assert D["proved"] is False
    assert one_d_weight_lb(7) < 6
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["D_type_not_lambda_constant"] is True
    assert out["proved"]["NL_type_means_can_be_lt_6"] is True
    assert out["proved"]["one_d_not_constant_p7"] is True
    assert out["proved"]["nl_type_means_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
