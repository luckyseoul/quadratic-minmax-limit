"""Tests for Prop 15.500 — central Max+ counts; NL stab ⊆ {±1}."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15500 import (
    hyp_1d_central,
    hyp_nl_central,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_nl_central_count():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"][5]["n_nl_central"] == 0
    assert A["by_p"][7]["n_nl_central"] == 42
    assert A["by_p"][7]["n_nl_central"] == hyp_nl_central(7)
    assert A["by_p"][5]["n_nl_central"] == hyp_nl_central(5)
    # fail-when-wrong: p≡3 formula at p=5; include 1D; p=3
    assert A["by_p"][5]["n_nl_central"] != 20
    assert A["by_p"][7]["n_nl_central"] != A["by_p"][7]["n_all_central"]
    assert A["by_p"][3]["n_nl_central"] != hyp_nl_central(3)
    assert A["by_p"][3]["n_nl_central"] == 0


def test_1d_central_count():
    B = prove_B()
    assert B["proved"] is True
    for p in (3, 5, 7):
        assert B["by_p"][p]["n_1d_central"] == hyp_1d_central(p)
        assert B["by_p"][p]["n_1d_central"] != n_1d(p)


def test_nl_stab_subseteq_pm1():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"][5]["n_cand"] == 0
    assert C["by_p"]["7_ord6"]["n_cand"] == 0
    assert C["by_p"]["7_ord3"]["hits"]["without0"] == 8
    assert C["by_p"]["7_ord3"]["hits"]["without0_1d"] == 8
    assert C["by_p"]["7_ord3"]["hits"]["without0_nl"] == 0


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert D["nS_free_named"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
