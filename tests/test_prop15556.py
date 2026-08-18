"""Tests for Prop 15.556 — G3_1d named; Φ=G3_1d+Φ_free."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15556 import (
    g3_1d_named,
    g3_1d_wrong_drop_mu3,
    g3_1d_wrong_nonspec,
    is_collinear,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_G3_1d_matches_1d_slice():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["ok"] is True
        assert row["maxerr"] < 1e-8
        assert row["drop_mu3_err"] > 1.0
        assert row["n_collinear"] > 0
        assert row["nonspec_err"] > 1.0
    # explicit p=5 non-collinear all-plus: n_spec=3, m=3, C=10
    # G3_1d = 3*(10/5) + 0 = 6
    assert abs(g3_1d_named(5, 1, 1, 1, False) - 6.0) < 1e-12
    # collinear plus: n_spec=1, G3_1d=2+(2)*(-2)= -2
    assert abs(g3_1d_named(5, 1, 1, 1, True) + 2.0) < 1e-12
    F = _kernel_field(5)
    # some square-line pair is collinear
    found = False
    for u in range(1, 25):
        for v in range(1, 25):
            if u != v and is_collinear(F, u, v):
                found = True
                assert g3_1d_wrong_nonspec(5, 1, 1, 1, True) != g3_1d_named(
                    5, 1, 1, 1, True
                )
                break
        if found:
            break
    assert found is True


def test_fail_when_wrong_dies():
    F = _kernel_field(5)
    # pick a collinear plus pair
    u = v = None
    for x in range(1, 25):
        for y in range(1, 25):
            if x != y and is_collinear(F, x, y) and F["chi"][x] == 1:
                u, v = x, y
                break
        if u is not None:
            break
    assert u is not None
    w = F["add"][v][F["neg"][u]]
    got = g3_1d_named(5, 1, 1, int(F["chi"][w]), True)
    assert abs(got - g3_1d_wrong_nonspec(5, 1, 1, int(F["chi"][w]), True)) > 1.0
    assert abs(got - g3_1d_wrong_drop_mu3(5, 1, 1, int(F["chi"][w]), True)) > 1.0


def test_flags_untouched_Phi_free_open():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["max_abs_free"] > 1.0
    assert B["by_p"]["5"]["n_split"] == 0
    assert B["by_p"]["7"]["n_split"] == 0
    C = prove_open()
    assert C["proved"] is False
    assert C["Phi_named_in_p"] is False
    assert C["Phi_free_named_in_p"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.556"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15556.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
