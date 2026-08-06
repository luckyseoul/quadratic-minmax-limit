"""Prop 15.180 — residual (i) support constraints; bound still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed
from e1_gmin_m4_prop15180 import (
    Q_pairs_dual_eq,
    hinge_status_180,
    main,
    score_parity_viable_nd,
    theorem_Q_pairs_negative,
    theorem_open_core_nd_ge_2,
)


def test_Q_pairs_negative():
    r = theorem_Q_pairs_negative([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    assert Q_pairs_dual_eq(5) == Fraction(-24, 5)
    assert Q_pairs_dual_eq(5) < 0


def test_open_core_nd_ge_2():
    r = theorem_open_core_nd_ge_2([5, 7, 11, 13])
    assert r["proved_open_core_nonempty"] is True
    for p in (5, 7, 11):
        s = score_parity_viable_nd(p)
        assert s["min_open_nd"] >= 2


def test_predicates_honest():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True  # residual ii done
    h = hinge_status_180()
    assert h["bound_proved_general"] is False
    assert h["residual_i_closed"] is False
    assert h["residual_ii_closed"] is True
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["mu4_or_gsum_bound_general"] is False
