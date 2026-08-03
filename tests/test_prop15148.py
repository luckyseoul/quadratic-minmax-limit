"""Tests for Prop 15.148 — relaxed-ULC residual calculus."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15147 import d4_suf, ulc_U  # noqa: E402
from e1_gmin_m4_prop15148 import (  # noqa: E402
    C7,
    C_act_from_delta2,
    C_flat,
    C_max,
    N_eval,
    P5_DELTA2,
    P7_DELTA2,
    kappa,
    main,
    primes_ge,
    prove_C_max_residual,
    prove_census_window,
    prove_linear_defect,
    prove_open,
    prove_uniform_criterion,
    relaxed_ULC_bound,
)


def test_linear_defect_and_C_max():
    a = prove_linear_defect(primes_ge(5, 31))
    assert a["proved"] is True
    b = prove_C_max_residual(primes_ge(5, 31))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert C_max(p) == C_flat(p) + kappa(p) * rho_min2(p)
        assert C_max(p) < 1
        assert relaxed_ULC_bound(p, C_max(p)) == d4_suf(p)
        # C_max closed form
        x = p * p
        Q = x**4 - 7 * x**3 + 71 * x**2 + 67 * x + 60
        assert C_max(p) == Fraction(Q, (x - 5) * (x - 2) * (x + 3) ** 2)


def test_uniform_C7_criterion():
    c = prove_uniform_criterion(primes_ge(5, 97))
    assert c["proved"] is True
    assert C_max(7) == C7 == Fraction(79923, 87373)
    assert N_eval(31) < 0
    assert N_eval(32) > 0
    for p in primes_ge(7, 50):
        assert C_max(p) >= C7
        assert relaxed_ULC_bound(p, C7) <= d4_suf(p)
    # equality at p=7
    assert relaxed_ULC_bound(7, C7) == d4_suf(7)


def test_census_window():
    d = prove_census_window()
    assert d["proved"] is True
    Ca5 = C_act_from_delta2(5, P5_DELTA2)
    Ca7 = C_act_from_delta2(7, P7_DELTA2)
    assert Ca5 < C_max(5)
    assert Ca7 < C_max(7)
    assert Ca7 <= C7
    # C=1 fails implication
    assert relaxed_ULC_bound(7, 1) > d4_suf(7)
    # window nonempty
    assert Ca7 <= C7
    # 9/10 under C7 but does not cover C_act(7)
    assert Fraction(9, 10) < Ca7
    assert Fraction(9, 10) < C7


def test_open_main():
    e = prove_open()
    assert e["residual_closed_general"] is False
    assert e["uniform_relaxed_ULC_proved"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["linear_defect_identity"] is True
    assert out["proved"]["uniform_C7_criterion_structure"] is True
    assert out["proved"]["uniform_relaxed_ULC_for_all_p_ge_7"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15148.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.148"
    assert data["closed_forms"]["C7"] == str(C7)
