"""Tests for Prop 15.156 — κ₄ residual dictionary."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15145 import ed4_suf_closed  # noqa: E402
from e1_gmin_m4_prop15153 import R4_from_W, mu4_from_R4  # noqa: E402
from e1_gmin_m4_prop15154 import eta_of_mu4, eta_suf  # noqa: E402
from e1_gmin_m4_prop15156 import (  # noqa: E402
    Es4_from_W,
    ed4_4des,
    kappa4_flat,
    kappa4_from_Es4,
    kappa4_from_eta,
    kappa4_suf,
    main,
    primes_ge,
    prove_bridge_eta,
    prove_census,
    prove_crude_and_lp,
    prove_design_orientation,
    prove_kappa4_dictionary,
    prove_open,
)


def test_kappa4_closed_forms():
    a = prove_kappa4_dictionary(primes_ge(5, 47))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        n = n_of(p)
        x = p * p
        assert kappa4_flat(p) == ed4_flat_closed(p) - 12 * n * n
        assert kappa4_suf(p) == ed4_suf_closed(p) - 12 * n * n
        assert kappa4_flat(p) == Fraction(16 * n * (x + 3), x - 5)
        assert kappa4_suf(p) == Fraction(
            4 * n * (9 * x * x + 22 * x - 15), x * (x - 5)
        )
        assert 0 < kappa4_flat(p) < kappa4_suf(p)


def test_bridge_and_census():
    b = prove_bridge_eta([5, 7])
    assert b["proved"] is True
    e = prove_census()
    assert e["proved_census"] is True
    for p in (5, 7):
        Es4 = Es4_from_W(p)
        k4 = kappa4_from_Es4(p, Es4)
        eta = eta_of_mu4(p, mu4_from_R4(p, R4_from_W(p)))
        assert k4 == kappa4_from_eta(p, eta)
        assert k4 <= kappa4_suf(p)
        n = n_of(p)
        n4 = n * (n - 1) * (n - 2) * (n - 3)
        assert kappa4_suf(p) == n4 * eta_suf(p) - 16 * n


def test_design_lp_crude():
    c = prove_design_orientation(primes_ge(5, 31))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17):
        assert ed4_4des(p) < ed4_flat_closed(p) < ed4_suf_closed(p)
        n, d = n_of(p), d_of(p)
        assert ed4_4des(p) == Fraction(3) * n**4 / Fraction(d * (d + 2))
    d = prove_crude_and_lp([5, 7, 11, 13, 17])
    assert d["proved_crude_dead"] is True
    assert d["proved_lp_dead"] is True
    for p in (5, 7, 11, 13, 17):
        assert d["lp"][str(p)]["ok"] is True
        assert d["lp"][str(p)]["ratio"] > 1.0


def test_open_and_main():
    g = prove_open()
    assert g["residual_closed_general"] is False
    assert g["kappa4_le_kappa4_suf_all_p"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["kappa4_dictionary_closed"] is True
    assert out["proved"]["bridge_kappa4_eta"] is True
    assert out["proved"]["design_4des_is_LB"] is True
    assert out["proved"]["crude_dead"] is True
    assert out["proved"]["moment_lp_dead"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15156.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.156"
    assert "kappa4_suf" in data["closed_forms"]
