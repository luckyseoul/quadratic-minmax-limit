"""Tests for Prop 15.132 — Max+-free dictionary; Aut; γ; dead envelopes."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_bud, m4f_flat  # noqa: E402
from e1_gmin_m4_prop15132 import (  # noqa: E402
    is_prime,
    main,
    moment_lp_max_ED4,
    pole_dmax_ED4,
    prove_aut_invariance,
    prove_dead_envelopes,
    prove_gamma0_mass,
    prove_gamma_parity,
    prove_residual_dictionary,
)


def test_residual_dictionary():
    primes = [p for p in range(5, 30) if is_prime(p)]
    r = prove_residual_dictionary(primes)
    assert r["proved"] is True
    for p in primes:
        assert m4f_bud(p) == m4f_flat(p) + room_hyp_over_24(p)
        if p >= 7:
            assert rho_min2(p) < room_hyp_over_24(p)
        if p == 5:
            assert rho_min2(p) > room_hyp_over_24(p)


def test_aut_and_gamma():
    assert prove_aut_invariance()["proved"] is True
    g = prove_gamma_parity()
    assert g["proved"] is True
    assert g["gamma_values"] == [-6, -4, -2, 0, 2, 4, 6]
    assert all(x % 2 == 0 for x in g["gamma_values"])
    # local Rayleigh at gamma=0 is 4p
    assert g["local_rayleigh_4p_minus_2gamma"]["7"]["0"] == 28


def test_gamma0_and_dead_envelopes():
    d = prove_gamma0_mass()
    assert d["proved"] is True
    assert d["p5"]["constant"] is True
    assert d["p5"]["zero_mass"] == 4350
    assert d["p7"]["n_unique_in_sample"] == 3

    primes = [p for p in range(5, 20) if is_prime(p)]
    e = prove_dead_envelopes(primes)
    assert e["proved"] is True
    for p in primes:
        assert e["by_p"][str(p)]["all_weak_vs_suf"] is True
        # moment LP exceeds suf
        lp = moment_lp_max_ED4(p)
        from e1_gmin_m4_prop15119 import ed4_flat_closed

        suf = float(ed4_flat_closed(p) + 24 * rho_min2(p))
        assert lp > suf
        # pole envelope with N=n exceeds suf
        assert pole_dmax_ED4(p, p * p + 1) > ed4_flat_closed(p) + 24 * rho_min2(p)


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["residual_dictionary_maxplus_free"] is True
    assert out["proved"]["aut_invariance_delta_Q_orbits"] is True
    assert out["proved"]["gamma_parity"] is True
    assert out["proved"]["dead_envelopes_p_ge_5"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15132.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.132"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
