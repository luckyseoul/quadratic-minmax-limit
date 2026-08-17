"""Tests for Prop 15.469 — Cov++ isolated; named atoms miss."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_plus
from e1_gmin_m4_prop15414 import Lpar_1d_named
from e1_gmin_m4_prop15469 import (
    LIVE_COV,
    atom,
    cov_pp,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_cov_values():
    A = prove_A()
    assert A["proved"] is True
    assert cov_pp(5) == LIVE_COV[5] == Fraction(10, 39)
    assert cov_pp(7) == LIVE_COV[7] == Fraction(2793, 1636)
    assert cov_pp(5) != 0
    assert cov_pp(5) != mu_plus(5) ** 2


def test_named_atoms_miss():
    B = prove_B()
    assert B["proved"] is True
    assert atom("Lpar_1d", 5) != cov_pp(5)
    assert atom("Lpar_1d", 5) - mu_plus(5) ** 2 != Fraction(10, 39)
    assert Lpar_1d_named(5) - mu_plus(5) ** 2 != cov_pp(5)
    assert atom("S_chi", 5) != cov_pp(5)
    assert atom("chi+", 7) != cov_pp(7)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.469"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["cov_values"] is True
    assert out["proved"]["named_atoms_miss_cov"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
