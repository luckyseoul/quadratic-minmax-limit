"""Tests for Prop 15.574 — μ_k3 Fejér 4-point."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15573 import q_sub_mix
from e1_gmin_m4_prop15574 import (
    main,
    mu_k3_P_of_one,
    mu_k3_drop_frac,
    mu_k3_four_p4,
    mu_k3_named,
    prove_A,
    prove_open,
)


def test_mu_k3_fejer_hits_live_and_fails():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        live = q_sub_mix(p)["parts"]["k3"]["mu"]
        assert abs(mu_k3_named(p, 2) - live) < 1e-6
        assert abs(mu_k3_named(p, 3) - live) < 1e-6
        assert A["by_p"][str(p)]["match"] is True
        assert A["by_p"][str(p)]["indep_r"] is True
        assert A["by_p"][str(p)]["four_p4_fails"] is True
        assert A["by_p"][str(p)]["P1_fails"] is True
    assert abs(mu_k3_named(5, 2) - mu_k3_four_p4(5)) > 1.0
    assert abs(mu_k3_named(7, 2) - mu_k3_drop_frac(7, 2)) > 1.0
    assert abs(mu_k3_named(7, 2) - mu_k3_P_of_one(7)) > 1.0


def test_flags_and_Q_tau_still_open():
    B = prove_open()
    assert B["proved"] is False
    assert B["Q_tau_named_in_p"] is False
    assert B["phi_F_imported"] is False
    assert B["p11_mu_spread"] > 1.0
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.574"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15574.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
