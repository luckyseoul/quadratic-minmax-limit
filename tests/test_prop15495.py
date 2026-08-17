"""Tests for Prop 15.495 — L_NL not named in p; 15.298 misses Q_NL."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import named_gj_ints
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15409 import D_form_on_lattice_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational
from e1_gmin_m4_prop15430 import Lpar_R_named
from e1_gmin_m4_prop15495 import (
    LIVE_LNL_PLUS,
    P7_DENS,
    fold_NL,
    gj_closure,
    main,
    p5_square_name,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_p5_law_is_n_or_qm1_and_dies_at_p7():
    assert LIVE_LNL_PLUS[5][("pm1",)] == Fraction(26) == Fraction(5 * 5 + 1)
    assert LIVE_LNL_PLUS[5][(1, 1, "sub")] == Fraction(24) == Fraction(5 * 5 - 1)
    assert p5_square_name(("pm1",), 1, 7) == 50
    assert p5_square_name((1, 1, "sub"), 1, 7) == 48
    assert LIVE_LNL_PLUS[7][("pm1",)] == Fraction(4435, 38)
    assert LIVE_LNL_PLUS[7][("pm1",)] != 50
    assert LIVE_LNL_PLUS[7][(1, 1, "sub")] == Fraction(2113, 19)
    assert LIVE_LNL_PLUS[7][(1, 1, "sub")] != 48
    dens = {v.denominator for v in LIVE_LNL_PLUS[7].values()}
    assert dens == {6, 19, 38, 57, 114}
    assert dens == P7_DENS - {2}
    assert 409 not in dens
    for d in dens:
        assert 409 % d != 0


def test_p7_slots_miss_gj_and_Lpar_R():
    S7 = named_gj_ints(7)
    clos = gj_closure(7)
    for val in LIVE_LNL_PLUS[7].values():
        assert val not in clos
    assert 19 not in S7
    assert 38 not in S7
    assert 57 not in S7
    assert 114 not in S7
    assert 2113 not in S7
    assert 673 not in S7
    assert Lpar_R_named(7) == 133
    assert LIVE_LNL_PLUS[7][(1, 1, "sub")] != Lpar_R_named(7)


def test_live_fold_matches_table():
    A = prove_A()
    assert A["proved"] is True
    pack5 = fold_NL(5)
    pack7 = fold_NL(7)
    assert pack5["n_1d"] == 30
    assert pack7["n_1d"] == 140
    assert pack5["n_NL"] == 100
    assert pack7["n_NL"] == 5586
    for p, pack in ((5, pack5), (7, pack7)):
        for tk, val in LIVE_LNL_PLUS[p].items():
            assert pack["L"][(1, tk)] == val


def test_298_live_hits_named_table_misses():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"][5]["live"] == "64/15"
    assert C["by_p"][7]["live"] == "632/171"
    assert C["by_p"][5]["p5law"] != str(LIVE_QNL[5])
    assert C["by_p"][7]["p5law"] != str(LIVE_QNL[7])
    assert C["by_p"][5]["mumu"] == "104/25"
    assert C["by_p"][7]["mumu"] == "216/49"


def test_mixture_still_has_N():
    assert Q_1d_pp_named(5) != LIVE_QNL[5]
    assert Q_1d_pp_named(7) != LIVE_QNL[7]
    assert Q_1d_pp_named(5) == Fraction(16, 9)
    assert LIVE_QNL[5] == Fraction(64, 15)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["Q_1d_eq_Q_NL"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert D_form_on_lattice_general() is False
    assert Q_NL_is_short_rational() is False
    B = prove_B()
    assert B["proved"] is True
    out = main()
    assert out["prop"] == "15.495"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert "phi_F_ge_6" in out["flags_not_flipped"]
