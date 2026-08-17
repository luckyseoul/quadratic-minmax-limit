"""Tests for Prop 15.288 — Type+ |U|-histogram in (p,e).

Fail-when-wrong: inverted binomial index; Type− forms; dim-2
formula on a dim-1 3-set.  Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15285 import n1_formula_3, n1d_set_count
from e1_gmin_m4_prop15288 import (
    main,
    nU_formula,
    nU_formula_inverted,
    prove_formula,
    prove_open,
)


def test_formula_matches_15285_n1_and_mass():
    for p, e in ((5, 0), (5, 1), (5, 2), (5, 3), (7, 0), (7, 2), (7, 3)):
        n0, n1, n2, n3 = nU_formula(p, e)
        pal = tuple(sorted((1,) * e + (-1,) * (3 - e)))
        assert n0 == n1_formula_3(p, pal, 2)
        assert n0 + n1 + n2 + n3 == n1d_set_count(p)


def test_invert_and_live_count():
    A = prove_formula([5, 7])
    assert A["proved"] is True
    assert A["invert_hit"] < A["invert_tot"]
    assert A["invert_tot"] > 0
    assert A["minus_eq"] == 0
    assert A["minus_tot"] > 0
    assert A["dim1_formula_applied_hits"] == 0
    for b in A["blocks"]:
        assert b["bad"] == 0
        assert b["checked"] > 0
    # e=3 at p=5 makes the inverted 3-fiber index collide (both 3).
    # e=0 is the fail-when-wrong witness.
    inv = nU_formula_inverted(5, 0)
    assert inv != nU_formula(5, 0)
    assert nU_formula_inverted(5, 3) == nU_formula(5, 3)


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["n_H_named_in_p"] is False
    assert op["n_NL_named_in_p"] is False
    assert op["tau_k_named"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["nU_formula"] is True
    assert out["proved"]["invert_is_identity"] is False
    assert out["proved"]["type_minus_is_identity"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15288.json"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.288"
    assert data["status"]["floor"] is False
