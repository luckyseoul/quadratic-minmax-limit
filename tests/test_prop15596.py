"""Tests for Prop 15.596 — 15.406 Theorem C extended to p=11 (full ensemble)."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15596 import theorem_extension_p11


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_p11_full_ensemble_result():
    r = theorem_extension_p11()
    assert r["sampled"] is False
    assert r["nU"] + r["nUc"] == r["n_total_ensemble"] == 37457112
    assert r["closed"] is True
    assert r["ker_mixed"] == 0
    assert r["aff_mixed"] == 0
    assert r["solvable"] is True


def test_rank_deficient_not_full_span():
    """rank(B_U) < n/2 at p=11 -- Theorem C holds WITHOUT B_U spanning
    everything, unlike p=5,7. A future proof via 'full rank' would be
    false; do not attempt it."""
    r = theorem_extension_p11()
    assert r["rank_BU"] < r["n_over_2"]
    assert r["rank_equals_n_over_2"] is False
    assert r["ker_dim"] == r["n"] - r["rank_BU"]


def test_still_only_a_census_point():
    """One further prime does not close leftover 2 (fable.md acceptance)."""
    from e1_gmin_m4_prop15596 import main
    out = main()
    assert out["general_p_proved"] is False
    assert out["L_status"] == "OPEN"
    assert 11 in out["census_primes"] and len(out["census_primes"]) == 4
