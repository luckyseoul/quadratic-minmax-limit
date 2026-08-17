"""Tests for Prop 15.280 — DΔT Weil expansion; circle |N̂_T|²; Hoffman.

Floor / phi_F_ge_6 / e1 / L stay open.  Fail-when-wrong is required.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general
from e1_gmin_m4_prop15280 import (
    main,
    prove_circle_nhat,
    prove_expansion,
    prove_hoffman_support,
    prove_open,
    ratio_r,
    unit_circle,
)


def test_expansion_identity_and_invert_cross():
    a = prove_expansion([5, 7])
    assert a["proved"] is True
    assert a["invert_cross_is_identity"] is False
    assert a["n_checks"] > 0


def test_circle_nhat_closed_form():
    b = prove_circle_nhat([5, 7, 11, 13])
    assert b["proved"] is True
    for p in (5, 7, 11, 13):
        row = b["by_p"][str(p)]
        assert row["r_formula"] is True
        assert row["pred"] == p * (p + 1) ** 2
        assert row["inverted_slot_vanishes"] is False
        assert row["slot_p(p+1)^2"].endswith(row["slot_p(p+1)^2"].split("/")[-1])


def test_r_vanishes_on_fp():
    F = _kernel_field(7)
    H = unit_circle(F)
    assert len(H) == 8
    r = ratio_r(F, H)
    assert int(r[0]) == 0
    assert int(r[1]) == 7
    for c in range(2, 7):
        assert int(r[c]) == 0
    for c in range(7, 49):
        assert int(r[c]) == 1


def test_hoffman_and_nonhoffman_and_p13():
    c = prove_hoffman_support()
    assert c["proved"] is True
    assert c["nonhoffman_leaks_p5"] is True
    assert c["p13_norm_circle_empty"] is True
    assert c["ystar"]["5"]["balanced"] is True
    assert c["ystar"]["7"]["pair_chi"] is True


def test_open_floor_flags_untouched():
    d = prove_open()
    assert d["floor_closed"] is False
    assert d["cross_closed"] is False
    assert d["Nplus_named"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["d_symdiff_expansion"] is True
    assert out["proved"]["circle_Nhat_eval"] is True
    assert out["proved"]["hoffman_support"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    assert out["proved"]["invert_cross_is_identity"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15280.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.280"
    assert data["status"]["floor"] is False
