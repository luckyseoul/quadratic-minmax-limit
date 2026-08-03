"""Tests for Prop 15.150 — srg triples, λ_e, π_e residual decomposition."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15123 import srg_degree  # noqa: E402
from e1_gmin_m4_prop15128 import KNOWN_N  # noqa: E402
from e1_gmin_m4_prop15150 import (  # noqa: E402
    LAMBDA_CERT,
    lambda_e,
    main,
    mixture_tot,
    n0_of,
    n1_of,
    n2_of,
    n3_of,
    pi_e,
    primes_ge,
    prove_lambda,
    prove_open,
    prove_pi,
    prove_residual_decomp,
    prove_triple_counts,
    solve_affine_lambda,
    triple_counts,
    triple_counts_from_srg,
)


def test_triple_counts():
    a = prove_triple_counts(primes_ge(3, 31))
    assert a["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        tc = triple_counts(p)
        assert tc == triple_counts_from_srg(p)
        assert sum(tc.values()) == comb(n_of(p), 3)
        assert tc[0] >= 0
        assert tc[3] == n3_of(p)
        assert tc[2] == n2_of(p)
        assert tc[1] == n1_of(p)
        assert tc[0] == n0_of(p)


def test_lambda_affine_formula():
    b = prove_lambda([5, 7])
    assert b["proved"] is True
    for p in (5, 7):
        N = KNOWN_N[p]
        A, B, lam = solve_affine_lambda(p, N)
        for e in range(4):
            assert lam[e] == lambda_e(p, e, N)
            assert lam[e] == LAMBDA_CERT[p][e]
        # λ_e = N(p+3-2e)/(8p)
        assert A == lambda_e(p, 0, N)
        assert B == lambda_e(p, 1, N) - lambda_e(p, 0, N)


def test_pi_e_maxplus_free():
    c = prove_pi(primes_ge(5, 31))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        pi = pi_e(p)
        assert sum(pi.values()) == 1
        assert all(pi[e] > 0 for e in range(4))
        # N-free: formula only uses srg params
        tc = triple_counts(p)
        tot = mixture_tot(p)
        assert tot == sum(tc[e] * (p + 3 - 2 * e) for e in range(4))
        assert tot == (p + 3) * comb(n_of(p), 3) - (n_of(p) - 2) * n_of(p) * srg_degree(
            p
        )


def test_residual_and_main():
    d = prove_residual_decomp()
    assert d["proved_structure"] is True
    assert d["census_ok"] is True
    e = prove_open()
    assert e["residual_closed_general"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["srg_triple_counts"] is True
    assert out["proved"]["lambda_e_affine_formula"] is True
    assert out["proved"]["pi_e_maxplus_free"] is True
    assert out["proved"]["m_e_bound_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15150.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.150"
    assert data["closed_forms"]["lambda_e"] == "N(p+3-2e)/(8p)"
