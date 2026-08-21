"""Tests for Prop 15.591 — closed form of ν_part; strict negativity on the locus.

Fail-when-wrong: the three character-sum lemmas are verified by exhaustive
direct summation (data-free, any p); the closed form is checked against the
independently measured K(mu_part) samples; negativity and the 12(p-1)/D bound
are checked over the entire locus.
"""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15591 import (
    lemma_A_S_kappa,
    lemma_B_S_star,
    lemma_C_S_phi,
    nu_part_closed,
    theorem_D_closed_form,
    theorem_E_negativity,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


@pytest.mark.parametrize("p", [5, 7])
def test_lemmas_exhaustive_small(p):
    assert lemma_A_S_kappa(p)["proved"]
    assert lemma_B_S_star(p)["proved"]
    assert lemma_C_S_phi(p)["proved"]
    assert theorem_D_closed_form(p)["proved"]


def test_lemmas_p11_sampled():
    assert lemma_A_S_kappa(11, samples=30)["proved"]
    assert lemma_B_S_star(11, samples=30)["proved"]
    assert lemma_C_S_phi(11, samples=30)["proved"]
    assert theorem_D_closed_form(11, samples=12)["proved"]


def test_measured_samples_match():
    # measured by the independent data pipeline (§10 of the convolution note)
    assert Fraction(7) * nu_part_closed(7, 3) == Fraction(-10, 77)
    # field element 2 (harmonic) is index 2*p in this module's encoding
    assert Fraction(11) * nu_part_closed(11, 22) == Fraction(-30, 319)


@pytest.mark.parametrize("p", [5, 7, 11])
def test_negativity_and_bound_on_locus(p):
    e = theorem_E_negativity(p)
    assert e["all_negative"] is True
    assert e["within_bound"] is True
    assert e["hasse_ok"] is True


def test_bound_tight_iff_supersingular():
    # p ≡ 3 mod 4: harmonic/equianharmonic fibers are supersingular, phi=2p,
    # so max|nu_part| attains 12(p-1)/(p^2(p^2-5)) exactly.
    for p in (7, 11):
        e = theorem_E_negativity(p)
        assert Fraction(e["max_abs"]) == Fraction(12 * (p - 1), p * p * (p * p - 5))
    # p=5 ≡ 1 mod 4: strictly below
    e5 = theorem_E_negativity(5)
    assert Fraction(e5["max_abs"]) < Fraction(12 * 4, 25 * 20)
