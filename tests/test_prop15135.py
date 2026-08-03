"""Tests for Prop 15.135 — coherent mass on G; halfspace char sums; residual OPEN."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15110 import e4_formula  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402
from e1_gmin_m4_prop15135 import (  # noqa: E402
    f_hs_char,
    main,
    prove_Ghs_incomplete,
    prove_halfspace_char,
    prove_moments_dont_pin,
    prove_spectral_coherent_mass,
)
from minmax_quadratic import halfspace_boolean_vector  # noqa: E402


def test_spectral_form():
    a = prove_spectral_coherent_mass()
    assert a["proved"] is True
    assert G_SPECTRUM[5]["nullity"] == 2
    assert G_SPECTRUM[7]["nullity"] == 7
    assert G_SPECTRUM[3]["nullity"] == 0


def test_halfspace_char_and_e4():
    b = prove_halfspace_char([5, 7, 11])
    assert b["proved"] is True
    for p in (5, 7, 11):
        hs = halfspace_boolean_vector(p)
        n = hs.shape[0]
        # full match for p=5,7
        if p <= 7:
            for S in combinations(range(n), 4):
                direct = int(hs[S[0]] * hs[S[1]] * hs[S[2]] * hs[S[3]])
                assert f_hs_char(p, S) == direct
        assert b["by_p"][str(p)]["e4_match"] is True
        assert b["by_p"][str(p)]["sum_f_hs"] == int(e4_formula(p))


def test_Ghs_dead_and_moments():
    c = prove_Ghs_incomplete()
    assert c["proved"] is True
    assert c["by_p"]["5"]["O_hs_with_antipodes"] == 60
    assert c["by_p"]["5"]["N_Max"] == 260
    assert c["by_p"]["7"]["O_hs_with_antipodes"] == 168
    assert c["dead_for_residual_UB"] is True
    assert prove_moments_dont_pin()["proved"] is True


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["spectral_coherent_mass_on_G"] is True
    assert out["proved"]["halfspace_char_formula"] is True
    assert out["proved"]["Ghs_incomplete_dead_UB"] is True
    assert out["proved"]["gauss_m4_on_orbits_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15135.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.135"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["closed_forms"]["nu_G"]["5"] == 2
    assert data["closed_forms"]["nu_G"]["7"] == 7
    # room still the hyp form
    assert room_hyp_over_24(5) == Fraction(1536, 65)
