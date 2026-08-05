"""Tests for Prop 15.165: exact Es4, closed Es4_*/η_*, GoG↔Φ."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15165 import (
    ES4_EXACT,
    Es4_exact_census,
    Es4_from_spectrum,
    Es4_star_closed,
    eta_from_Es4,
    eta_star_closed,
    main,
    prove_theorem_B,
    prove_theorem_C,
    prove_theorem_E,
    prove_theorem_F,
    SPECTRUM_P3,
)
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7
from e1_gmin_m4_prop15164 import Es4_star, eta_star, sixteen_N_from_dminus1_Es4
from e1_gmin_m4_prop15100 import d_of


def test_closed_forms_match_api():
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert Es4_star_closed(p) == Es4_star(p)
        assert eta_star_closed(p) == eta_star(p)


def test_closed_forms_algebra():
    c = prove_theorem_C()
    assert c["proved"]
    assert c["asymp"]["Es4_star_over_n2_to"] == 12
    assert c["asymp"]["eta_star_over_n_to"] == str(Fraction(19, 6))


def test_exact_es4_from_spectrum():
    assert Es4_from_spectrum(SPECTRUM_P3, 3) == ES4_EXACT[3]
    assert Es4_from_spectrum(SPECTRUM_P5, 5) == ES4_EXACT[5]
    assert Es4_from_spectrum(SPECTRUM_P7, 7) == ES4_EXACT[7]


def test_es4_le_star_census():
    e = prove_theorem_E()
    assert e["certified"]
    assert ES4_EXACT[3] == Es4_star(3)  # saturation at p=3
    assert ES4_EXACT[5] <= Es4_star(5)
    assert ES4_EXACT[7] <= Es4_star(7)
    # p=7 W_CENSUS trap
    assert e["W_CENSUS_p7_is_not_global_Es4"]["equal"] is False


def test_sixteen_N_with_exact_es4():
    for p, mult in ((5, 13), (7, 25)):
        pred = sixteen_N_from_dminus1_Es4(
            mult, d_of(p), ES4_EXACT[p], p
        )
        assert pred["sixteen_N"] is True


def test_eta_room():
    # p=3 saturates η; p=5,7 strict room
    assert eta_from_Es4(3, ES4_EXACT[3]) == eta_star(3)
    assert eta_from_Es4(5, ES4_EXACT[5]) < eta_star(5)
    assert eta_from_Es4(7, ES4_EXACT[7]) < eta_star(7)


def test_H_saturation_p5():
    f = prove_theorem_F()
    assert f["proved_saturation_p5"]
    assert f["p7"]["strict"]


def test_GoG_Phi_link():
    b = prove_theorem_B()
    assert b["proved"]


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Es4_star_for_all_p"] is False
    assert out["proved"]["exact_Es4_census_p3_5_7"] is True
    assert out["proved"]["Es4_star_eta_star_closed_form"] is True
    assert Es4_exact_census(11) is None
