"""Tests for Prop 15.106 — rest-average-8 + kurtosis residual form."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of, wick_hi  # noqa: E402
from e1_gmin_m4_prop15105 import dim_Z, tr_Phi_Z  # noqa: E402
from e1_gmin_m4_prop15106 import (  # noqa: E402
    certify_rest_and_kurtosis,
    main,
    prove_general_status,
    prove_kurtosis_equivalence,
    prove_rest_average_8,
)


def test_rest_avg_equals_8_algebra():
    for p in [5, 7, 11, 13, 17, 19, 23, 31]:
        n = n_of(p)
        d = d_of(p)
        m = dim_Z(p)
        T = tr_Phi_Z(p)
        rest_dim = m - d
        assert rest_dim > 0
        assert Fraction(T - 16 * d, rest_dim) == 8


def test_rest_average_module():
    r = prove_rest_average_8([3, 5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    for p in [5, 7, 11]:
        assert r["by_p"][str(p)]["equals_8"] is True


def test_kurtosis_threshold():
    for p in [3, 5, 7, 11, 13]:
        n = n_of(p)
        assert Fraction(wick_hi(p), 4 * n * n) == 3 + Fraction(12, n)


def test_kurtosis_module():
    r = prove_kurtosis_equivalence([3, 5, 7, 11, 13, 17])
    assert r["proved"] is True


def test_general_status_open():
    g = prove_general_status()
    assert g["rest_average_8_proved"] is True
    assert g["kurtosis_equivalence_proved"] is True
    assert g["sixteen_N_proved_for_all_p"] is False
    assert g["kappa_bound_proved_for_all_p"] is False


def test_cert_if_available():
    for p in (3, 5):
        c = certify_rest_and_kurtosis(p)
        if c.get("skipped"):
            continue
        assert c["kurtosis_le_thr"] is True
        assert c["lambda_max_le_16"] is True
        assert c["kappa2_le_96n"] is True


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["rest_average_8"] is True
    assert out["proved"]["kurtosis_equivalence"] is True
    assert out["proved"]["sixteen_N_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15106.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.106"
    assert data["L_status"] == "OPEN"
