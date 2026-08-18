"""Tests for Prop 15.569 — A3 product is a Gauss expansion."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import gauss_table, prod_dbl_k
from e1_gmin_m4_prop15569 import (
    a3_drop8,
    a3_drop_inv,
    a3_named_expansion,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_a3_pair_expansion_matches_prod():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["fails_ok"] is True
        assert row["drop_inv_err"] > 1.0
        assert row["drop8_err"] > 1.0
    G, omega = gauss_table(7)
    tau = G[1]
    k = [0] * 7
    k[0], k[1], k[2] = 2, -1, -1
    got = prod_dbl_k(tuple(k), 7, G, omega)
    assert abs(got - a3_named_expansion(7, 1, 2, tau)) < 1e-8
    assert abs(got - a3_drop_inv(7, 1, 2, tau)) > 1.0
    assert abs(got - a3_drop8(7, 1, 2, tau)) > 1.0


def test_family_mean_16p2_only_at_p7_and_flags():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["by_p"]["7"]["named_mean"] - 16.0 * 49.0) < 1e-6
    assert abs(B["by_p"]["11"]["named_mean"] - 16.0 * 121.0) > 1.0
    C = prove_open()
    assert C["proved"] is False
    assert C["A_full_dbl_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.569"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15569.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
