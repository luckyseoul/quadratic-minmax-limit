"""Tests for Prop 15.285 — n_1D Type+ fiber binomial.

Fail-when-wrong: Type− swap, inverted binomial, pal-without-dim,
n=n_1D at p=5 on a line 4-set.  Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15285 import (
    C,
    main,
    n1_formula_3,
    n1_formula_4_dim1,
    n1d_set_count,
    prove_formula_3,
    prove_formula_4_dim1,
    prove_lifts,
    prove_open,
)


def test_closed_forms_values():
    assert n1d_set_count(5) == 30
    assert n1d_set_count(7) == 140
    assert n1d_set_count(3) == 6
    # p=5 m=3
    assert n1_formula_3(5, (1, 1, 1), 1) == 4
    assert n1_formula_3(5, (1, 1, 1), 2) == 3
    assert n1_formula_3(5, (-1, -1, -1), 1) == 0
    assert n1_formula_3(5, (-1, 1, 1), 2) == 2
    assert n1_formula_3(5, (-1, -1, 1), 2) == 1
    # p=7 m=4
    assert n1_formula_3(7, (1, 1, 1), 1) == 18
    assert n1_formula_3(7, (1, 1, 1), 2) == 16
    assert n1_formula_3(7, (-1, -1, -1), 1) == 4
    assert n1_formula_4_dim1(5, (1, 1, 1, 1, 1, 1)) == 4
    assert n1_formula_4_dim1(5, (-1, -1, -1, -1, -1, -1)) == 0
    assert n1_formula_4_dim1(7, (1, 1, 1, 1, 1, 1)) == 15
    assert n1_formula_4_dim1(7, (-1, -1, -1, -1, -1, -1)) == 0
    # pal-without-dim is not a function
    assert n1_formula_3(5, (1, 1, 1), 1) != n1_formula_3(5, (1, 1, 1), 2)
    assert n1_formula_3(7, (1, 1, 1), 1) != n1_formula_3(7, (1, 1, 1), 2)


def test_lifts_and_formula_3():
    L = prove_lifts([3, 5, 7])
    assert L["proved"] is True
    a = prove_formula_3([5, 7])
    assert a["proved"] is True
    assert a["pal_without_dim_splits"] is True
    assert a["invert_tot"] > 0
    assert a["invert_eq"] == 0
    for b in a["blocks"]:
        assert b["bad"] == 0
        assert b["checked"] > 0


def test_formula_4_dim1_and_ignore_nl():
    b = prove_formula_4_dim1([5, 7])
    assert b["proved"] is True
    assert b["ignore_nl_dead"] is True
    for row in b["blocks"]:
        assert row["bad"] == 0
        assert row["typeplus_n1"] != row["typeplus_n1m"]
        if row["p"] == 5:
            assert row["nNL_pos"] > 0
            assert row["typeplus_n1"] == 4
            assert row["typeplus_n1m"] == 0


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["n_NL_named_in_p"] is False
    assert op["Nplus_named"] is False
    assert op["term4_reduced_to_1d"] is False
    assert phi_F_ge_6_proved_general() is False
    assert C(2, 3) == 0


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["typeplus_binomial"] is True
    assert out["proved"]["formula_3"] is True
    assert out["proved"]["formula_4_dim1"] is True
    assert out["proved"]["pal_without_dim_is_identity"] is False
    assert out["proved"]["invert_binom_is_identity"] is False
    assert out["proved"]["n_equals_n1D_p5_4sets"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    assert out["proved"]["n_NL_named_in_p"] is False
    assert out["proved"]["term4_reduced_to_1d"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15285.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.285"
    assert data["status"]["floor"] is False
