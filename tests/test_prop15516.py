"""Tests for Prop 15.516 — ensemble L_ns_mix=μ₊μ₋; Q_τ unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_minus, mu_plus
from e1_gmin_m4_prop15448 import Lns_11_named, Lns_mix_named
from e1_gmin_m4_prop15516 import ensemble_slot_L, main, prove_A, prove_open


def test_ensemble_nsmix_is_mupm():
    A = prove_A()
    assert A["proved"] is True
    mean5, vals5 = ensemble_slot_L(5, "ns_mix")
    assert mean5 == Lns_mix_named(5) == mu_plus(5) * mu_minus(5)
    assert len(vals5) >= 2
    assert mean5 != Lns_11_named(5)
    mean7, vals7 = ensemble_slot_L(7, "ns_mix")
    assert mean7 == Lns_mix_named(7)
    assert len(vals7) >= 2


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["phi_F_ge_6"] is False
