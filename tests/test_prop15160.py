"""Tests for Prop 15.160 — H vs dual-gap algebra."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15160 import (
    H_of,
    certify_H_and_gap,
    prove_open,
    prove_theorem_A,
    prove_theorem_B,
    prove_theorem_C,
    thr_ray_of,
)


def test_theorem_A_H_le_thr_ray():
    r = prove_theorem_A([3, 5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    assert thr_ray_of(5) == H_of(5)
    assert thr_ray_of(7) > H_of(7)
    # identity thr_ray - H = (3p+7)(p-5)/(2d)
    for p in (5, 7, 11, 13):
        d = (p * p + 1) // 2
        assert thr_ray_of(p) - H_of(p) == Fraction((3 * p + 7) * (p - 5), 2 * d)


def test_theorem_B_and_C():
    assert prove_theorem_B()["proved"] is True
    assert prove_theorem_C()["proved"] is True
    assert H_of(3) == 5
    assert H_of(5) < 5


def test_certified_spectra_H_gap():
    c = certify_H_and_gap()
    assert c["certified"] is True
    assert c["p5"]["ray_eq_H"] is True
    assert c["p5"]["H_eq_thr_ray"] is True
    assert c["p7"]["ray_lt_H"] is True


def test_honest_open():
    o = prove_open()
    assert o["residual_closed_general"] is False
    assert o["H_proved_for_all_p"] is False
