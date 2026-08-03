"""Tests for Prop 15.146 — type-free R₄/μ₄ residual channel."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15124 import ed4_from_exact3  # noqa: E402
from e1_gmin_m4_prop15146 import (  # noqa: E402
    P5_DELTA2,
    P7_DELTA2,
    Pf_energy_ub,
    R4_bud,
    R4_flat,
    R4_suf,
    ed4_from_exact3_closed,
    main,
    mu4_flat,
    mu4_suf,
    primes_ge,
    prove_R4_dictionary,
    prove_census,
    prove_central_moments,
    prove_ed4_exact3,
    prove_open,
    prove_spectral_mass,
    w_star,
)


def test_ed4_exact3_closed():
    a = prove_ed4_exact3(primes_ge(5, 31))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert ed4_from_exact3_closed(p) == ed4_from_exact3(p)
        n = p * p + 1
        assert ed4_from_exact3_closed(p) == (
            -Fraction(n**4) + 28 * Fraction(n * n) - 40 * Fraction(n)
        )


def test_R4_dictionary_slacks():
    b = prove_R4_dictionary(primes_ge(5, 31))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17):
        assert R4_suf(p) - R4_flat(p) == Fraction(3, 2) * rho_min2(p)
        assert R4_bud(p) - R4_flat(p) == Fraction(3, 2) * room_hyp_over_24(p)
        # δ² = (2/3)(R4 - R4_flat) reconstruction for certified masses
        if p == 5:
            R4 = R4_flat(p) + Fraction(3, 2) * P5_DELTA2
            assert Fraction(2, 3) * (R4 - R4_flat(p)) == P5_DELTA2
            assert R4 <= R4_suf(p)
        if p == 7:
            R4 = R4_flat(p) + Fraction(3, 2) * P7_DELTA2
            assert Fraction(2, 3) * (R4 - R4_flat(p)) == P7_DELTA2
            assert R4 <= R4_suf(p)


def test_central_and_spectral():
    c = prove_central_moments(primes_ge(5, 31))
    assert c["proved"] is True
    for p in (5, 7, 11):
        assert mu4_suf(p) - mu4_flat(p) == Fraction(3, 2) * rho_min2(p)
    d = prove_spectral_mass(primes_ge(5, 31))
    assert d["proved"] is True
    assert d["dead_for_residual"] is True
    for p in (5, 7, 11, 13):
        assert w_star(p) == Fraction((p * p - 3) * (p * p - 4), p * p * (p * p - 1))
        assert Pf_energy_ub(p) > rho_min2(p)
        assert Pf_energy_ub(p) > room_hyp_over_24(p)


def test_census_open_main():
    e = prove_census()
    assert e["proved"] is True
    f = prove_open()
    assert f["residual_closed_general"] is False
    assert f["type_free_R4_le_R4_suf_all_p"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["R4_residual_dictionary"] is True
    assert out["proved"]["spectral_mass_w_star_too_weak"] is True
    assert out["proved"]["type_free_R4_le_R4_suf_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15146.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.146"
    assert data["status"]["jensen_spectral_mass_viable"] is False
