"""Tests for Prop 15.133 — class_key Bose–Mesner; F19; Aut-line; CR dead."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15133 import (  # noqa: E402
    CK_SPECTRUM,
    main,
    prove_aut_line_form,
    prove_classkey_spectrum,
    prove_cr_dead,
    prove_f19_p7,
    prove_success_locus,
)


def test_aut_line_and_classkey_spectrum():
    assert prove_aut_line_form()["proved"] is True
    r = prove_classkey_spectrum()
    assert r["proved"] is True
    assert CK_SPECTRUM[3]["nullity"] == 0
    assert CK_SPECTRUM[5]["nullity"] == 1
    assert CK_SPECTRUM[7]["nullity"] == 0
    assert CK_SPECTRUM[5]["lambda_max"] == 20.0
    assert CK_SPECTRUM[3]["lambda_max"] < CK_SPECTRUM[3]["4p"]
    assert CK_SPECTRUM[7]["lambda_max"] < CK_SPECTRUM[7]["4p"]


def test_f19_quantitative_p7():
    c = prove_f19_p7()
    assert c["proved"] is True
    assert CK_SPECTRUM[7]["nullity"] == 0
    assert P7_DELTA2_PAIR > 0
    assert P7_DELTA2_PAIR == Fraction(19180800, 1840091)


def test_cr_dead_and_success_locus():
    assert prove_cr_dead()["proved"] is True
    e = prove_success_locus()
    assert e["proved"] is True
    assert e["p5"]["match"] is True
    assert room_hyp_over_24(5) == Fraction(1536, 65)


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["classkey_spectrum_p3_5_7"] is True
    assert out["proved"]["f19_quantitative_p7"] is True
    assert out["proved"]["cr_dead_as_Aut"] is True
    assert out["proved"]["dim_true_Aut_le_1_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15133.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.133"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["closed_forms"]["classkey_nullity"]["5"] == 1
    assert data["closed_forms"]["classkey_nullity"]["7"] == 0
