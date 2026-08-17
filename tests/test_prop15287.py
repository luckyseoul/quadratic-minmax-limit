"""Tests for Prop 15.287 — Paley-parent Hoffman ncirc.

Fail-when-wrong: 3p formula; AP=Paley at p=7 or p=13; AP=p at p=13.
Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15287 import (
    main,
    paley_ncirc_formula,
    prove_ap_not_paley,
    prove_open,
    prove_paley_ncirc,
)


def test_paley_formula_kills_threep():
    assert paley_ncirc_formula(5) == 10
    assert paley_ncirc_formula(7) == 21
    assert paley_ncirc_formula(11) == 55
    assert paley_ncirc_formula(13) == 78
    assert paley_ncirc_formula(11) != 3 * 11
    A = prove_paley_ncirc([5, 7, 11])
    assert A["proved"] is True
    b11 = next(b for b in A["blocks"] if b["p"] == 11)
    assert b11["counts"] == [55, 55, 55]
    assert b11["formula"] == 55


def test_ap_splits_except_p5():
    B = prove_ap_not_paley([5, 7])
    assert B["proved"] is True
    byp = {b["p"]: b for b in B["blocks"]}
    assert byp[5]["equal"] is True
    assert byp[5]["paley"] == 10
    assert byp[7]["equal"] is False
    assert byp[7]["ap"] == 7
    assert byp[7]["paley"] == 21


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["n_H_named_in_p"] is False
    assert op["n_NL_named_in_p"] is False
    assert op["Nplus_named"] is False
    assert op["ap_ncirc_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_includes_p13_ap_zero():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["paley_ncirc_formula"] is True
    assert out["proved"]["ap_not_paley"] is True
    assert out["proved"]["three_p_is_identity"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15287.json"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.287"
    assert data["status"]["floor"] is False
    b13 = next(b for b in data["B"]["blocks"] if b["p"] == 13)
    assert b13["ap"] == 0
    assert b13["paley"] == 78
    a13 = next(b for b in data["A"]["blocks"] if b["p"] == 13)
    assert a13["formula"] == 78
    assert all(c == 78 for c in a13["counts"])
