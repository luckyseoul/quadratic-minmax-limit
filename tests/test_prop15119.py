"""Tests for Prop 15.119 — residual budgets; weight enum; halfspace pin."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import (  # noqa: E402
    certify_weight_enum,
    crude_ed4_envelope,
    ed4_bud_closed,
    ed4_flat_closed,
    ew2_bud,
    ew2_flat,
    is_prime,
    m4f_bud,
    m4f_bud_closed,
    m4f_flat,
    m4f_flat_closed,
    m4f_of_ed4,
    main,
    prove_ed4_budgets,
    prove_ew2_channel,
    prove_m4f_dictionary,
    prove_we_structure,
)


def test_ed4_budget_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_ed4_budgets(primes)
    assert r["proved"] is True
    # p=5 flat and bud
    assert ed4_flat_closed(5) == Fraction(43472, 5)
    assert ed4_bud_closed(5) == Fraction(120400, 13)


def test_ew2_channel_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_ew2_channel(primes)
    assert r["proved"] is True
    # gap = room_hyp/4 at p=5
    gap = ew2_bud(5) - ew2_flat(5)
    assert gap == 6 * room_hyp_over_24(5)
    assert gap == Fraction(24 * room_hyp_over_24(5), 4)


def test_m4f_dictionary_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_m4f_dictionary(primes)
    assert r["proved"] is True
    assert m4f_flat_closed(5) == Fraction(598, 5)
    assert m4f_bud_closed(5) == Fraction(1862, 13)
    assert m4f_bud(5) == m4f_flat(5) + room_hyp_over_24(5)
    # consistency with ED4 dictionary at p=5
    assert m4f_of_ed4(5, ed4_bud_closed(5)) == m4f_bud_closed(5)


def test_we_structure_and_envelope_weak():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_we_structure(primes)
    assert r["proved"] is True
    assert crude_ed4_envelope(5) > ed4_bud_closed(5)
    assert crude_ed4_envelope(7) > ed4_bud_closed(7)


def test_cert_weight_enum_and_main_open():
    for p in (3, 5):
        c = certify_weight_enum(p)
        if c.get("skipped"):
            continue
        assert c["hs_is_evec"] is True
        assert c["parity_ok"] is True
        assert c["dmax_ok"] is True
        assert c["spectrum_matches_expected"] is True
        assert c["ED4_le_bud"] is True
        assert c["saturates_bud"] is True
        assert c["m4f_hs_eq_bud"] is True
        assert c["EW2_le_bud"] is True
        assert c["delta2_le_hyp"] is True
        if p == 5:
            assert Fraction(c["orth_times_N"]) == 147456

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["ed4_budgets"] is True
    assert out["proved"]["ew2_channel"] is True
    assert out["proved"]["m4f_dictionary"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15119.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.119"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
