"""Tests for Prop 15.219 — R-matrix decomposition; residual-(i) hinge still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15219 import (  # noqa: E402
    hinge_status_219,
    residual_i_closed_via_219,
    theorem_R_structure,
    theorem_Tr_RC_eq_pn_m4plus,
    theorem_cross_outside_formulas,
    theorem_formula_implies_two_over_n,
    theorem_inside_eq_4kappa_over_p,
)


def test_R_structure_proved():
    assert theorem_R_structure()["proved"] is True


def test_Tr_RC_proved():
    assert theorem_Tr_RC_eq_pn_m4plus()["proved"] is True


def test_inside_proved():
    assert theorem_inside_eq_4kappa_over_p()["proved"] is True


def test_cross_outside_algebra():
    t = theorem_cross_outside_formulas()
    assert t["proved"] is True
    # Explicit Fraction check at p=5, m=3/65, κ=1
    p, n = 5, 26
    m = Fraction(3, 65)
    k = 1
    inside = 4 * Fraction(k, p)
    cross = 8 * p * m - 8 * Fraction(k, p)
    outside = p * (n - 8) * m + 4 * Fraction(k, p)
    assert inside + cross + outside == p * n * m
    assert inside == Fraction(4, 5)


def test_formula_implies_two_over_n():
    t = theorem_formula_implies_two_over_n()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert is_prime(p)
        # naive majorant |4κ−φ|≤2p ⇒ |μ|≤2/n
        assert Fraction(2 * p, p * n_of(p)) == two_over_n(p)
        assert two_over_n(p) <= mu4_threshold(p)


def test_residual_i_still_open():
    assert residual_i_closed_via_219() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_219()
    assert h["R_structure"] is True
    assert h["Tr_RC"] is True
    assert h["inside_4kappa_over_p"] is True
    assert h["cross_outside"] is True
    assert h["formula_implies_two_over_n"] is True
    assert h["residual_i_closed_via_219"] is False
    assert "OPEN" in h["open_hinge"] or "1/(2p)" in h["open_hinge"]


def test_p5_census_formula_corroboration():
    """Optional census: at p=5, μ=(4κ−φ)/(pn) on |κ|=1 (not a general proof)."""
    import numpy as np
    from itertools import combinations

    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(np.int8)
    n = C.shape[0]
    Mp = boolean_evecs(C.astype(float), float(p), 0).astype(np.int8)
    Mm = boolean_evecs(C.astype(float), -float(p), 0).astype(np.int8)
    Np, Nm = len(Mp), len(Mm)
    max_abs = Fraction(0)
    n_k1 = 0
    for S in combinations(range(n), 4):
        a, b, c, d = S
        kap = int(
            C[a, b] * C[c, d]
            + C[a, c] * C[b, d]
            + C[a, d] * C[b, c]
        )
        if abs(kap) != 1:
            continue
        n_k1 += 1
        phi = sum(
            int(C[r, a] * C[r, b] * C[r, c] * C[r, d])
            for r in range(n)
            if r not in S
        )
        m4p = Fraction(int((Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]).sum()), Np)
        m4m = Fraction(int((Mm[:, a] * Mm[:, b] * Mm[:, c] * Mm[:, d]).sum()), Nm)
        mu = (m4p + m4m) / 2
        pred = Fraction(4 * kap - phi, p * n)
        assert mu == pred
        assert m4p == m4m
        max_abs = max(max_abs, abs(mu))
    assert n_k1 == 11700
    assert max_abs == Fraction(3, 65)
    assert max_abs <= mu4_threshold(5)
    assert max_abs <= two_over_n(5)
