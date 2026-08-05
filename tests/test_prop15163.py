"""Tests for Prop 15.163: Wick m4 mass, eta room, H-split census."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15163 import (
    certify_H_split,
    certify_T_counts,
    certify_eta_room_usage,
    eta_room_16N,
    main,
    mean_T2,
    prove_theorem_A,
    prove_theorem_B,
    wick_T2_sum,
    wick_m4_sumsq,
)
from e1_gmin_m4_prop15162 import m4_sumsq_budget_16N
from e1_gmin_m4_prop15100 import n_of


def test_mean_T2_identity():
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        assert mean_T2(p) == Fraction(3 * (n - 5), n - 3)
        # 3 - 6/(n-3)
        assert mean_T2(p) == 3 - Fraction(6, n - 3)


def test_wick_closed_forms():
    a = prove_theorem_A()
    assert a["proved"]
    for p in (5, 7, 11, 13):
        assert wick_m4_sumsq(p) == wick_T2_sum(p) / (p**4)
        n = n_of(p)
        assert wick_T2_sum(p) == Fraction(n * (n - 1) * (n - 2) * (n - 5), 8)


def test_eta_room_algebra():
    b = prove_theorem_B()
    assert b["proved"]
    for p in (5, 7, 11, 13, 17, 19):
        room = eta_room_16N(p)
        assert room == m4_sumsq_budget_16N(p) - wick_m4_sumsq(p)
        n = n_of(p)
        assert room == Fraction(n * (19 * p * p - 3), 6 * p * p)
        assert room > 0


def test_T_counts_census():
    c = certify_T_counts()
    assert c["certified"]
    for p in (5, 7):
        assert c["by_p"][str(p)]["match"]


def test_eta_room_usage_census():
    e = certify_eta_room_usage()
    assert e["certified"]
    # p=5 uses most of the room; p=7 less
    assert e["by_p"]["5"]["eta2_over_room"] < 1.0
    assert e["by_p"]["7"]["eta2_over_room"] < 1.0


def test_H_split_census():
    path5 = Path("/tmp/maxplus_p5.npy")
    path7 = Path("/tmp/maxplus_p7.npy")
    if not (path5.is_file() and path7.is_file()):
        pytest.skip("Max+ caches missing")
    h = certify_H_split()
    assert h["certified"]
    for p in (5, 7):
        r = h["by_p"][str(p)]
        assert r["is_1_plus_d_minus_1"]
        assert r["mult_eq_d"]


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["wick_T2_and_m4_mass"] is True
    assert out["proved"]["mult_ge_d_all_p"] is False
    assert out["proved"]["residual_for_all_p_ge_5"] is False
