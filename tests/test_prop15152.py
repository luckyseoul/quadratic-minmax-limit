"""Tests for Prop 15.152 — free-param t-structure; multi-orbit; residual≡R4."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15149 import k_act_from_delta2, k_suf  # noqa: E402
from e1_gmin_m4_prop15150 import mixture_tot  # noqa: E402
from e1_gmin_m4_prop15151 import E_te, regular_set_identities  # noqa: E402
from e1_gmin_m4_prop15152 import (  # noqa: E402
    E_t3,
    MULTI_TYPE_WEIGHTS_P7,
    P5_DELTA2,
    P7_DELTA2,
    T_BY_K_P5,
    cov_k_Ck3_from_R4,
    main,
    primes_ge,
    prove_E_t3,
    prove_cov_reduction_p5,
    prove_free_param,
    prove_multitype_p7,
    prove_open,
    prove_residual_is_R4,
    prove_t3_closed_p5,
    t3_closed_p5,
    t_from_t3,
)


def test_free_param_and_t3_p5():
    a = prove_free_param()
    assert a["proved"] is True
    b = prove_t3_closed_p5()
    assert b["proved_census"] is True
    for k, t in T_BY_K_P5.items():
        if k < 3:
            continue
        rec = t_from_t3(k, 5, t[3])
        assert all(rec[e] == t[e] for e in range(4))
        assert regular_set_identities(k, 5, t)["ok"]
    assert t3_closed_p5(3) == 0
    assert t3_closed_p5(4) == 4
    assert t3_closed_p5(5) == 10
    assert t3_closed_p5(10) == 130


def test_p7_multitype_corrects_constancy():
    c = prove_multitype_p7()
    assert c["proved_census"] is True
    assert c["corrects_15_151_B"] is True
    assert c["n_multi_weights"] == 10
    assert set(MULTI_TYPE_WEIGHTS_P7) == set(range(16, 35, 2))
    assert max(MULTI_TYPE_WEIGHTS_P7.values()) == 8
    # free param yields nonneg integer t-vectors on samples
    for k, hist in c["t3_hist_sample"].items():
        for t3_s, _ in hist.items():
            t3 = int(t3_s)
            rec = t_from_t3(int(k), 7, t3)
            t_int = tuple(int(rec[e]) for e in range(4))
            assert regular_set_identities(int(k), 7, t_int)["ok"]
            assert all(x >= 0 for x in t_int)


def test_residual_R4_and_cov_reduction():
    d = prove_residual_is_R4()
    assert d["proved"] is True
    for p, d2 in ((5, P5_DELTA2), (7, P7_DELTA2)):
        row = d["by_p"][str(p)]
        assert row["match_k_act"] is True
        assert row["Cov_le_budget"] is True
        assert Fraction(row["E_mu"]) == k_act_from_delta2(p, d2)
        assert Fraction(row["E_mu"]) <= k_suf(p)
    e = prove_cov_reduction_p5()
    assert e["proved"] is True
    assert e["recon_ok"] is True
    assert e["gamma"] == {"0": -1, "1": 3, "2": -3, "3": 1}


def test_E_t3_open_and_main():
    f = prove_E_t3(primes_ge(5, 31))
    assert f["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert E_t3(p) == E_te(p, 3)
        assert E_t3(p) == Fraction(
            __import__("e1_gmin_m4_prop15150", fromlist=["n3_of"]).n3_of(p)
            * (p - 3),
            8 * p,
        )
    g = prove_open()
    assert g["residual_closed_general"] is False
    assert g["t_e_of_k_closed_form_all_p"] is False
    assert g["residual_sum_equals_R4_channel"] is True
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["free_param_t_reduction"] is True
    assert out["proved"]["t3_closed_p5"] is True
    assert out["proved"]["p7_multitype_nonconstancy"] is True
    assert out["proved"]["residual_sum_is_R4_channel"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15152.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.152"
    assert data["status"]["pure_te_of_k_blocked_p7"] is True
