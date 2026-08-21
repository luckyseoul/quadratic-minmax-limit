"""Tests for Prop 15.590 — μ/δ contraction closure; p=5 determination; level-4 kill.

Fail-when-wrong: every assertion is either an exact integer identity on
enumerated Max±, an exact-elimination rank, or an external census anchor
(max|μ| = 3/65 at p=5, 109/2863 at p=7).  Nothing is fitted.
"""
from __future__ import annotations

import os
from fractions import Fraction

import numpy as np
import pytest

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15590 import MuLab, psd_kill

FULL = os.environ.get("PROP15590_FULL", "") == "1"

_cache: dict = {}


def lab5() -> MuLab:
    if "lab5" not in _cache:
        _cache["lab5"] = MuLab(5, with_deg6=True)
    return _cache["lab5"]


def lab7() -> MuLab:
    if "lab7" not in _cache:
        _cache["lab7"] = MuLab(7, with_deg6=FULL)
    return _cache["lab7"]


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_p5_enumeration_and_census_anchor():
    lb = lab5()
    assert len(lb.Yp) == 260 and len(lb.Ym) == 260  # ground truth |Max±|
    k1 = np.abs(lb.kap) == 1
    mx = int(np.abs(lb.mu4[k1]).max())
    assert Fraction(mx, lb.N) == Fraction(3, 65)  # census max|μ| on |κ|=1
    assert lb.s_identity()


def test_p5_contraction_identities_exact():
    lb = lab5()
    tested, bad = lb.check_star_mu(exhaustive=True)
    assert tested == 7800 and bad == 0
    _, bado = lb.check_out_mu(samples=800)
    assert bado == 0
    _, bad6 = lb.check_deg6(samples=400)
    assert bad6 == [0, 0, 0, 0]


def test_p5_delta_dead_gives_m4plus_eq_m4minus_on_kappa1():
    assert lab5().kappa1_delta_dead()


def test_p5_deg4_kernel_dim_one_touching_kappa1():
    lb = lab5()
    sysres = lb.equivariant_system(deg6=False)
    assert not sysres["inconsistent"]
    assert sysres["kernel_dim"] == 1
    kv = sysres["kernel"][0]
    k1 = np.abs(lb.kap) == 1
    k1labs = {int(x) for x in np.unique(lb.labM4[k1])} - lb.deadM4
    assert any(kv[sysres["maps"]["M4"][l]] != 0 for l in k1labs)


def test_p5_joint_deg46_determination():
    lb = lab5()
    j = lb.equivariant_system(deg6=True)
    assert not j["inconsistent"]
    assert j["kernel_dim"] == 0
    assert lb.solution_matches_data(j)  # exact per-set equality with data


def test_p5_level4_psd_kill():
    lb = lab5()
    d4 = lb.equivariant_system(deg6=False)
    res = psd_kill(lb, d4)
    assert res["killed"] is True
    assert res["lambda_min_at_point"] > -1e-8
    assert res["feasible_mu_int"] > res["threshold_2n_int"]
    assert res["feasible_mu_int"] > res["threshold_L_int"]
    # and the bound does hold at the true point — the relaxation, not the
    # statement, is what fails
    assert res["true_max_mu_int"] < res["threshold_2n_int"]


def test_p7_enumeration_census_anchor_and_deg4_kernel():
    lb = lab7()
    assert len(lb.Yp) == 11452  # ground truth |Max+|
    k1 = np.abs(lb.kap) == 1
    mx = int(np.abs(lb.mu4[k1]).max())
    assert Fraction(mx, lb.N) == Fraction(109, 2863)  # census anchor
    assert lb.s_identity()
    assert lb.kappa1_delta_dead()
    _, bad = lb.check_star_mu(samples=300)
    assert bad == 0
    sysres = lb.equivariant_system(deg6=False)
    assert not sysres["inconsistent"]
    assert sysres["kernel_dim"] == 2


@pytest.mark.skipif(not FULL, reason="set PROP15590_FULL=1 for the slow p=7 blocks")
def test_p7_level4_psd_kill_full():
    lb = lab7()
    d4 = lb.equivariant_system(deg6=False)
    res = psd_kill(lb, d4)
    assert res["killed"] is True
    assert res["feasible_mu_int"] > res["threshold_L_int"]
    assert res["true_max_mu_int"] < res["threshold_2n_int"]


@pytest.mark.skipif(not FULL, reason="set PROP15590_FULL=1 for the slow p=7 blocks")
def test_p7_joint_deg46_kernel_dim_four():
    lb = lab7()
    _, bad6 = lb.check_deg6(samples=200)
    assert bad6 == [0, 0, 0, 0]
    j = lb.equivariant_system(deg6=True)
    assert not j["inconsistent"]
    assert j["kernel_dim"] == 4
