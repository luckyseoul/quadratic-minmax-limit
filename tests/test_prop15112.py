"""Tests for Prop 15.112 — design moments; conf ‖κ‖²; ED4 residual dictionary."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15112 import (  # noqa: E402
    certify_design_moments,
    conf_kappa2,
    ed4_bud,
    ed4_flat,
    ed4_suf,
    is_prime,
    main,
    prove_conf_kappa2,
    prove_ed4_dictionary,
)


def test_conf_kappa2_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_conf_kappa2(primes)
    assert r["proved"] is True
    # spot: p=5 direct was 40950
    assert conf_kappa2(5) == Fraction(40950)
    assert conf_kappa2(3) == Fraction(450)


def test_ed4_dictionary():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_ed4_dictionary(primes)
    assert r["proved"] is True
    for p in primes:
        if p < 7:
            continue
        assert ed4_suf(p) < ed4_bud(p)
        assert rho_min2(p) < delta2_budget_hyp(p)
        assert ed4_suf(p) - ed4_flat(p) == 24 * rho_min2(p)


def test_ed4_suf_spot():
    # p=5: equality case, actual ED4 = ED4_bud = ED4_suf? 
    # at p=5 ρ_min² > budget so suf > bud
    assert ed4_suf(5) > ed4_bud(5)
    assert ed4_flat(5) + 24 * delta2_budget_hyp(5) == ed4_bud(5)


def test_design_cert_p3_p5():
    for p in (3, 5):
        c = certify_design_moments(p)
        if c.get("skipped"):
            continue
        assert c["frame_EyyT_eq_2Pplus"] is True
        assert c["E_D2_eq_2n"] is True
        assert c["delta2_le_rho_min2"] is True
        assert c["ED4_le_suf"] is True


def test_design_cert_p7_if_available():
    c = certify_design_moments(7)
    if c.get("skipped"):
        return
    assert c["frame_EyyT_eq_2Pplus"] is True
    assert c["E_D2_eq_2n"] is True
    assert c["delta2_le_rho_min2"] is True
    assert c["ED4_le_suf"] is True
    assert c["all_halfspace_sums"] is True


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["conf_kappa2"] is True
    assert out["proved"]["ed4_dictionary"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15112.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.112"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
