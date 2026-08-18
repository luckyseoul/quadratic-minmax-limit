"""Tests for Prop 15.553 — Term0 of K_all is a named Kloosterman sum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15550 import omega_set_local, pp_r
from e1_gmin_m4_prop15553 import (
    main,
    prove_A,
    prove_B,
    prove_open,
    term0_coeff,
    term0_coeff_G,
    term0_coeff_drop16,
    term0_coeff_drop_chi,
)


def test_term0_coeff_live_and_field():
    A = prove_A()
    assert A["proved"] is True
    assert abs(A["by_p"]["5"]["coeff_re"] - 5.0) < 1e-8
    assert abs(A["by_p"]["7"]["coeff_re"] - 101.0) < 1e-8
    assert abs(A["by_p"]["11"]["coeff_re"] - 197.0) < 1e-8
    assert abs(A["by_p"]["13"]["coeff_re"] - 485.0) < 1e-8
    assert abs(A["by_p"]["5"]["Term0"] - 650.0) < 1e-4
    assert abs(A["by_p"]["7"]["Term0"] - 578326.0) < 1e-4


def test_fail_when_wrong_dies():
    B = prove_B()
    assert B["proved"] is True
    F = _kernel_field(5)
    xi = omega_set_local(F)[0]
    r = pp_r(F)
    got = term0_coeff(F, xi, r)
    assert abs(got - 5) < 1e-8
    assert abs(got - term0_coeff_drop_chi(F, xi, r)) > 1.0
    assert abs(got - term0_coeff_G(F, xi, r)) > 1.0
    # 16=1 in F_5: /16 cannot fail there; it must fail at p=7.
    assert B["by_p"]["5"]["drop16_invisible_p5"] is True
    assert B["by_p"]["7"]["drop16"] > 1.0
    for p in (5, 7, 11, 13):
        row = B["by_p"][str(p)]
        assert row["ok"] is True
        assert row["drop_chi"] > 1.0
        assert row["Kl_to_G"] > 1.0
        if p != 5:
            assert row["drop16"] > 1.0


def test_flags_untouched_Kall_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["K_all_named_in_p"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    assert C["term0_is_bulk"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.553"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15553.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
