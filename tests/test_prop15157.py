"""Tests for Prop 15.157 — Gegenbauer design-defect residual."""
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
from e1_gmin_m4_prop15156 import Es4_from_W, ed4_4des  # noqa: E402
from e1_gmin_m4_prop15157 import (  # noqa: E402
    Es4_from_mu_G4,
    a0_of,
    a2_of,
    a4_of,
    dim_harm,
    main,
    mu_G4_flat,
    mu_G4_from_Es4,
    mu_G4_suf,
    primes_ge,
    prove_bound_survey,
    prove_census,
    prove_defect_residual,
    prove_gegenbauer_coeffs,
    prove_open,
)


def test_gegenbauer_coeffs():
    a = prove_gegenbauer_coeffs([3, 5, 7, 9, 13, 25, 61, 85, 145])
    assert a["proved"] is True
    for d in (13, 25, 61, 85, 145, 181):
        assert a0_of(d) + a2_of(d) + a4_of(d) == 1
        assert a0_of(d) == Fraction(3, d * (d + 2))
        assert a2_of(d) == Fraction(6 * (d - 1), d * (d + 4))
        assert a4_of(d) == Fraction((d - 1) * (d + 1), (d + 2) * (d + 4))


def test_mu_G4_suf_and_reconstruction():
    c = prove_defect_residual(primes_ge(5, 47))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        n, d = n_of(p), d_of(p)
        x = p * p
        ms = mu_G4_suf(p)
        num = 4 * (21 * x**3 + 19 * x**2 + 35 * x - 75) * (x + 9)
        den = x * (x - 5) * (x + 1) ** 3 * (x + 3) * (x - 1)
        assert ms == Fraction(num, den)
        assert Es4_from_mu_G4(p, ms) == ed4_suf_closed(p)
        assert Es4_from_mu_G4(p, mu_G4_flat(p)) == ed4_flat_closed(p)
        assert Fraction(n**4) * a0_of(d) == ed4_4des(p)
        assert ms > mu_G4_flat(p) > 0


def test_census_defect():
    d = prove_census()
    assert d["proved_census"] is True
    for p in (5, 7):
        Es4 = Es4_from_W(p)
        mu = mu_G4_from_Es4(p, Es4)
        assert mu > 0  # not a 3-design
        assert mu <= mu_G4_suf(p)
        assert Es4_from_mu_G4(p, mu) == Es4
        assert mu_G4_flat(p) < mu


def test_bounds_open_main():
    e = prove_bound_survey()
    assert e["proved_survey"] is True
    assert e["details"]["bound_1_over_h4_false_p5"] is True
    assert e["details"]["bound_d_over_h4_too_weak_p7"] is True
    # h4 formula sanity
    assert dim_harm(13, 4) == 1729
    g = prove_open()
    assert g["residual_closed_general"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["gegenbauer_t4_coeffs"] is True
    assert out["proved"]["defect_residual_dictionary"] is True
    assert out["proved"]["census_mu_G4"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15157.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.157"
    assert "mu_G4_suf" in data["closed_forms"]
