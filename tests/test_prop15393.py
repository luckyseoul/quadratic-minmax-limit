"""Tests for Prop 15.393 — first_nc(T) empty at p=37; Paley-half."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15393 import (
    LIVE,
    least_nsq,
    main,
    n_clique_circles,
    paley_half_hits,
    prove_A,
    prove_B,
    prove_open,
)


def test_paley_half_live_and_empty():
    A = prove_A()
    assert A["proved"] is True
    assert paley_half_hits(5) == [2, 3]
    assert paley_half_hits(13) == [5]
    assert least_nsq(13) == 2
    assert 2 not in paley_half_hits(13)
    assert paley_half_hits(37) == []
    for p, c in LIVE.items():
        assert min(paley_half_hits(p)) == c


def test_no_circle_p37():
    B = prove_B()
    assert B["proved"] is True
    rec = n_clique_circles(37)
    assert rec["Maxplus"] is True
    assert rec["n_clique"] == 0


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["paley_half_matches_or_empty"] is True
    assert out["proved"]["no_circle_at_37"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
