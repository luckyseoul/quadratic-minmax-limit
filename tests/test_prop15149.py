"""Tests for Prop 15.149 — size-bias residual form."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15147 import d3, d4_suf  # noqa: E402
from e1_gmin_m4_prop15148 import C7  # noqa: E402
from e1_gmin_m4_prop15149 import (  # noqa: E402
    P5_DELTA2,
    P7_DELTA2,
    P_ind,
    e_flat,
    e_flat_closed,
    k_C7,
    k_act_from_delta2,
    k_flat,
    k_shift_flat,
    k_shift_flat_closed,
    k_suf,
    main,
    primes_ge,
    prove_independence_excess,
    prove_k_C7,
    prove_k_flat_shift,
    prove_open,
    prove_size_bias,
    prove_triple_program_and_census,
    size_bias_allowance,
)


def test_size_bias_forms():
    a = prove_size_bias(primes_ge(5, 31))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        assert k_suf(p) == 3 + Fraction(n - 3) * d4_suf(p) / d3(p)
        assert k_suf(p) - k_flat(p) == size_bias_allowance(p)
        assert size_bias_allowance(p) == (
            12 * rho_min2(p) / Fraction(n * (n - 2) * (n + 2))
        )
    # census residual
    assert k_act_from_delta2(5, P5_DELTA2) <= k_suf(5)
    assert k_act_from_delta2(7, P7_DELTA2) <= k_suf(7)
    assert k_act_from_delta2(5, P5_DELTA2) > k_flat(5)
    assert k_act_from_delta2(7, P7_DELTA2) > k_flat(7)


def test_k_flat_shift_and_independence():
    b = prove_k_flat_shift(primes_ge(5, 31))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert k_shift_flat(p) == k_shift_flat_closed(p)
        assert 0 < k_shift_flat(p) < 3
    c = prove_independence_excess(primes_ge(5, 31))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        assert d3(p) - Fraction(1, 8) == Fraction(3, 8 * (n - 1))
        assert e_flat(p) == e_flat_closed(p) > 0
        assert P_ind(p * p) > 0
        assert k_flat(p) > Fraction(n + 3, 2)  # > binomial size-bias


def test_k_C7_and_triple():
    d = prove_k_C7(primes_ge(5, 50))
    assert d["proved"] is True
    assert k_C7(7) == k_suf(7)
    for p in primes_ge(7, 30):
        assert k_C7(p) <= k_suf(p)
    e = prove_triple_program_and_census()
    assert e["proved_structure"] is True
    assert e["census_ok"] is True


def test_open_main():
    f = prove_open()
    assert f["residual_closed_general"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["size_bias_residual_form"] is True
    assert out["proved"]["independence_excesses"] is True
    assert out["proved"]["size_bias_bound_for_all_p_ge_7"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15149.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.149"
    assert data["status"]["gauss_sum_triple_program_ready"] is True
