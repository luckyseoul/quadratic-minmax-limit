"""Tests for Prop 15.441 — Var_2orb is not floor-sufficient for p≡1≥13."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15323 import Var_2orb
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1
from e1_gmin_m4_prop15441 import (
    Q_from_V,
    V_lo,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_V2_below_Vlo_except_p5():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["V2_lt_Vlo"] is False
    assert A["by_p"]["13"]["V2_lt_Vlo"] is True
    assert A["by_p"]["17"]["V2_lt_Vlo"] is True
    assert Var_2orb(13) == Fraction(3276, 11)
    assert V_lo(13) == Fraction(819, 2)
    assert Var_2orb(13) < V_lo(13)
    assert Q_from_V(13, Var_2orb(13)) == Fraction(40, 11)
    assert Q_from_V(13, Var_2orb(13)) < Q_lo_named(13) == Fraction(15, 4)
    # fail-when-wrong neighbour
    assert not (Var_2orb(13) >= V_lo(13))


def test_mode_below_6_at_Q2():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["ge6"] is True
    assert B["by_p"]["13"]["ge6"] is False
    assert even_S_p1(13, Fraction(40, 11))["sm_j0"] == Fraction(48, 11)
    assert even_S_p1(13, Fraction(40, 11))["sm_j0"] < 6


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.441"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["V2_below_Vlo"] is True
    assert out["proved"]["mode_below_6"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["by_p"]["13"]["Q2"] == "40/11"
    assert out["algebra"]["B"]["by_p"]["13"]["sm_j0"] == "48/11"
