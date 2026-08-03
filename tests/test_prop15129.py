"""Tests for Prop 15.129 — Jensen inequality; Hoffman geometry; r̄ formula."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed  # noqa: E402
from e1_gmin_m4_prop15127 import W_hoffman_int  # noqa: E402
from e1_gmin_m4_prop15128 import ED4_from_W, W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15129 import (  # noqa: E402
    certify_hoffman_geometry_p5,
    is_prime,
    main,
    prove_average_replication,
    prove_census_residual,
    prove_dictionary,
    prove_jensen,
    r_bar,
    sum_m4_sq_from_ED4,
)


def test_jensen_statement():
    a = prove_jensen()
    assert a["proved"] is True
    assert "room_hyp/24" in a["theorem"]


def test_dictionary_census():
    b = prove_dictionary()
    assert b["proved"] is True
    for p in (3, 5, 7):
        ED4 = ED4_from_W(W_CENSUS[p], n_of(p))
        sm = sum_m4_sq_from_ED4(p, ED4)
        row = b["by_p"][str(p)]
        assert row["match"] is True
        assert row["ED4_le_bud"] is True
        assert Fraction(row["ED4"]) == ED4


def test_replication_and_p7_nondesign():
    primes = [p for p in range(3, 40) if is_prime(p)]
    d = prove_average_replication(primes)
    assert d["proved"] is True
    assert d["p7_hoffman_not_1_design"] is True
    # p=5: r_bar = 3
    assert r_bar(5) == Fraction(3)
    assert r_bar(5).denominator == 1
    # p=7: r_bar = 11*8/50 = 44/25
    assert r_bar(7) == Fraction(11 * 8, 50)
    assert r_bar(7) == Fraction(44, 25)
    assert r_bar(7).denominator != 1
    # p=3: r_bar = 2
    assert r_bar(3) == Fraction(2)
    # consistency with formula
    for p in (3, 5, 7, 11):
        assert r_bar(p) == Fraction(W_hoffman_int(p) * (p + 1), n_of(p))


def test_geometry_census_main_open():
    c = certify_hoffman_geometry_p5()
    if not c.get("skipped"):
        assert c["match"] is True
        assert c["disjoint_pairs"] == 30
        assert c["disjoint_pairs_eq_W12_half"] is True

    e = prove_census_residual()
    assert e["proved_at_census_primes"] is True
    assert e["proved_for_all_p"] is False

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["jensen_coherent_mass"] is True
    assert out["proved"]["p7_hoffman_not_1_design"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["weight_enumerator_closed_general"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15129.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.129"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["proved"]["p7_hoffman_not_1_design"] is True
