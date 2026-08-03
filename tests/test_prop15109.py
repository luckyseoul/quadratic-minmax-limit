"""Tests for Prop 15.109 — Φ–m4 identity; Aut δ; PF+rank obstruction."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15109 import (  # noqa: E402
    certify_aut_delta_p5,
    certify_phi_m4_identity,
    is_prime,
    kappa_B2_zero_diag,
    main,
    prove_kappa_B2_identity,
    prove_pf_rank_obstruction,
)
import numpy as np


def test_kappa_B2_identity():
    r = prove_kappa_B2_identity([8, 10, 12], trials=25)
    assert r["proved"] is True
    assert r["max_abs_err"] < 1e-9


def test_kappa_B2_unit_formula_spot():
    rng = np.random.default_rng(1)
    n = 10
    B = rng.normal(size=(n, n))
    B = 0.5 * (B + B.T)
    np.fill_diagonal(B, 0.0)
    B /= np.linalg.norm(B, "fro")
    # direct
    from itertools import combinations

    direct = 0.0
    for i, j, k, l in combinations(range(n), 4):
        t = B[i, j] * B[k, l] + B[i, k] * B[j, l] + B[i, l] * B[j, k]
        direct += t * t
    assert abs(direct - kappa_B2_zero_diag(B)) < 1e-10


def test_pf_rank_obstruction_all_primes():
    primes = [p for p in range(5, 50) if is_prime(p)]
    r = prove_pf_rank_obstruction(primes)
    assert r["proved"] is True
    for p in primes:
        assert r["by_p"][str(p)]["strict_lambda2_lt_d_over_N"] is True
        assert r["by_p"][str(p)]["remaining_after_m1_at_ceiling"] > 0


def test_phi_m4_identity_p3_p5():
    for p in (3, 5):
        c = certify_phi_m4_identity(p, ntrials=12)
        if c.get("skipped"):
            continue
        assert c["identity_holds"] is True
        assert c["le_16"] is True
        assert c["le_5_over_4"] is True


def test_aut_delta_p5_equality():
    c = certify_aut_delta_p5()
    if c.get("skipped"):
        return
    assert c["nullity_Aut_E4p"] == 1
    assert c["rho_min2_match"] is True
    assert c["equality_delta2_budget"] is True
    assert c["sum_rho2_eq_T"] is True


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["phi_m4_identity"] is True
    assert out["proved"]["pf_rank_obstruction"] is True
    assert out["proved"]["sum_rho2_le_T_for_all_p"] is False
    assert out["proved"]["sixteen_N_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15109.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.109"
    assert data["L_status"] == "OPEN"
