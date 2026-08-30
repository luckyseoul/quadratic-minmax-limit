"""Tests for Prop 15.474 — Ω-support and Plancherel split."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15474 import main, omega_pack, prove_A, prove_B, prove_C, prove_open


def test_omega_support():
    A = prove_A()
    assert A["proved"] is True
    rec = omega_pack(5)
    assert rec["leak"] < 1e-10
    assert rec["sigma"] == -1
    assert rec["live_k"] == 10
    assert rec["live_k"] not in (0, 25)
    assert omega_pack(7)["leak"] < 1e-10
    assert omega_pack(7)["sigma"] == 1


def test_plancherel_split():
    B = prove_B()
    assert B["proved"] is True
    rec = omega_pack(5)
    assert rec["sum_err"] < 1e-8
    assert rec["split_err"] < 1e-8
    assert abs(rec["kqmk"] - 150) < 1e-9
    assert abs(rec["delta_mean"]) < 1e-6


def test_delta_not_constant():
    C = prove_C()
    assert C["proved"] is True
    assert omega_pack(5)["nunique"] >= 4
    assert omega_pack(7)["nunique"] >= 5


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.474"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["omega_support"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
