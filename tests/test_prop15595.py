"""Tests for Prop 15.595 — the delta-hierarchy; two open roots, not three."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general
from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15595 import (
    closes_with,
    hierarchy,
    leftover2_error_vs_signal,
    req_leftover1,
    req_leftover3,
    req_residual_i,
)

PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_scope_e1_closed_general_is_not_a_close():
    """The corrected global predicate and four-unit gate both remain open."""
    assert e1_closed_general() is False
    assert phi_F_ge_6_proved_general() is False  # leftover 1 still open
    from e1_main_chain_status import run_main_chain
    out = run_main_chain()
    assert out["L_status"] == "OPEN"
    assert out["docs"]["four_e1_units"]["closed"] is False


def test_hierarchy_strictly_ordered_all_primes():
    """leftover1 < leftover3 < residual_i: leftover 1 is the binding one."""
    for p in PRIMES:
        h = hierarchy(p)
        assert h["ordered"] is True
        assert req_leftover1(p) < req_leftover3(p, Fraction(1550, 100))
        assert req_leftover3(p, Fraction(1550, 100)) < req_residual_i(p)


def test_leftover1_threshold_tends_to_n_over_12():
    for p in (1009,):
        assert abs(float(req_leftover1(p) / (p * p + 1)) - 1 / 12) < 0.005


def test_n_over_12_closes_leftovers_1_and_3():
    c = Fraction(1, 12)
    res = closes_with(c, primes=(11, 13, 17, 23, 47))
    for p, r in res.items():
        assert r["leftover1"] is True
        assert r["leftover3"] is True
        assert r["residual_i"] is True


def test_leftover2_not_reachable_by_L2_delta_bound():
    """Error bar overtakes the signal at p=11 and diverges like ~p/11."""
    assert leftover2_error_vs_signal(5)["useless"] is False   # only small p
    for p in (11, 13, 17, 23):
        assert leftover2_error_vs_signal(p)["useless"] is True
    r11 = leftover2_error_vs_signal(11)["ratio"]
    r23 = leftover2_error_vs_signal(23)["ratio"]
    assert r23 > r11 > 1.0        # diverging, not merely large


def test_measured_delta_meets_threshold_only_from_p11():
    """Measured ||delta||^2/n = .9089/.2085/.01941 at p=5/7/11 vs req/n."""
    measured = {5: 0.9089, 7: 0.2085, 11: 0.01941}
    for p, m in measured.items():
        thr = float(req_leftover1(p) / (p * p + 1))
        assert (m <= thr) == (p == 11)
