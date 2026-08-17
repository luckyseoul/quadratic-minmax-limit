"""Tests for Prop 15.432 — L_ns Paley×norm indexed."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_minus, mu_plus
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15432 import (
    LIVE_LNS,
    Lns_const_wrong,
    Lns_p5_law,
    main,
    prove_open,
)


def test_table_and_fail():
    mm5 = mu_plus(5) * mu_minus(5)
    assert mm5 == Fraction(25, 2)
    assert Lns_const_wrong(5) == mm5
    vals5 = {v[0] for v in LIVE_LNS[5].values()}
    assert vals5 == {12, 13}
    assert mm5 not in vals5
    n3 = {v[0] for k, v in LIVE_LNS[7].items() if k[2] == 3}
    assert len(n3) == 2
    assert Fraction(2733, 38) in n3
    assert Fraction(1398, 19) in n3


def test_sum_and_p5_law():
    for p, tab in LIVE_LNS.items():
        twomm = 2 * mu_plus(p) * mu_minus(p)
        for v in tab.values():
            assert v[0] + v[1] == twomm
    for (a, b, N), (lp, lm) in LIVE_LNS[5].items():
        assert lp == Lns_p5_law(N, 5)
    assert Lns_p5_law(3, 7) != LIVE_LNS[7][(-1, -1, 3)][0]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert LIVE_QNL[5] > Q_1d_pp_named(5)
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False


def test_main():
    out = main()
    assert out["proved"]["Lns_not_N_only"] is True
    assert out["proved"]["Lns_sum_2mumu"] is True
    assert out["proved"]["p5_law_dies_p7"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["fold_matches_LIVE"] is True
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.432"
