"""Tests for Prop 15.142 — uniform affine law; partners; p=11 sample."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15142 import (  # noqa: E402
    AFFINE_AP_CENSUS,
    FOURTHS_PARTNER,
    main,
    prove_affine_all_S,
    prove_fourths_partners,
    prove_k_AP_split,
    prove_open,
    prove_p11_sample,
    prove_size12_fibre,
)


def test_affine_all_S():
    # p=11 full C(11,6)=462 is ok; keep all three for prop claim
    a = prove_affine_all_S([5, 7, 11])
    assert a["proved"] is True
    for p in (5, 7, 11):
        row = a["by_p"][str(p)]
        assert row["all_S_evec"] is True
        assert row["n_evec"] == AFFINE_AP_CENSUS[p]["n_evec"]
        assert row["n_k_AP"] == AFFINE_AP_CENSUS[p]["n_k_AP"]


def test_k_AP_and_fourths():
    b = prove_k_AP_split()
    assert b["proved"] is True
    assert AFFINE_AP_CENSUS[5]["n_non_AP"] == 0
    assert AFFINE_AP_CENSUS[7]["n_non_AP"] == 14
    c = prove_fourths_partners()
    assert c["proved"] is True
    assert c["by_p"]["5"]["hit_hs_base"] is True
    assert c["by_p"]["7"]["hit_hs_base"] is False
    assert FOURTHS_PARTNER[7]["any_hit"] is False


def test_size12_and_p11():
    d = prove_size12_fibre()
    assert d["proved"] is True
    assert d["explicit_T_evec"] is True
    assert d["k"] == 12
    e = prove_p11_sample()
    assert e["proved"] is True
    assert e["live_sizes"]["AP_interval"] == 330
    assert e["live_sizes"]["nonAP_012346"] == 660
    assert e["live_sizes"]["nonAP_012457"] == 132
    assert e["live_sizes"]["ystar"] == 3630


def test_open_and_main():
    f = prove_open()
    assert f["residual_closed_general"] is False
    assert f["residual_closed_p5_maxplus_free"] is True
    assert f["residual_closed_p7_maxplus_free"] is True
    assert f["residual_closed_p11"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["affine_all_S_p5_7_11"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15142.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.142"
