"""Tests for Prop 15.239 — Per φ / Per κ Max+-free eigenforms; residual open."""
from __future__ import annotations

import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15238 import (  # noqa: E402
    mu_part_coeffs,
    p4_mupart_plus_2phi,
)
from e1_gmin_m4_prop15239 import (  # noqa: E402
    hinge_status_239,
    residual_i_closed_via_239,
    theorem_Per_kappa,
    theorem_Per_mupart_unconditional,
    theorem_Per_phi,
    theorem_delta_Per_on_kappa1,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_theorems_proved():
    assert theorem_Per_phi()["proved"] is True
    assert theorem_Per_kappa()["proved"] is True
    assert theorem_Per_mupart_unconditional()["proved"] is True
    assert theorem_delta_Per_on_kappa1()["proved_structure"] is True
    assert theorem_delta_Per_on_kappa1()["bound_proved_general"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_239() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_239()
    assert h["Per_phi"] is True
    assert h["Per_kappa"] is True
    assert h["Per_mupart_unconditional"] is True
    assert h["residual_i_closed_via_239"] is False


def _kap_phi(C, fours):
    a, b, c, d = fours.T
    kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    n = C.shape[0]
    ph = np.zeros(len(fours), dtype=np.int64)
    for i, S in enumerate(fours):
        s = set(S.tolist())
        t = 0
        for r in range(n):
            if r in s:
                continue
            pr = 1
            for j in S:
                pr *= int(C[r, j])
            t += pr
        ph[i] = t
    return kap, ph


def _per_apply(C, fours, S_row, vec, perms):
    acc = np.zeros(len(fours), dtype=np.int64)
    for sig in perms:
        pr = np.ones(len(fours), dtype=np.int64)
        for i in range(4):
            pr *= C[S_row[i], fours[:, sig[i]]]
        acc += pr
    return int(acc @ vec)


def test_live_Per_phi_kappa_p3_all_fours():
    """Drive permanent: Per φ=(2n+1)φ and Per κ=(n−1)²κ−6φ on all 4-sets p=3."""
    p = 3
    C = paley_conference_prime_power(p).astype(np.int32)
    n = C.shape[0]
    fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap, ph = _kap_phi(C, fours)
    perms = list(permutations(range(4)))
    for j in range(len(fours)):
        S = fours[j]
        Pph = _per_apply(C, fours, S, ph, perms)
        Pk = _per_apply(C, fours, S, kap, perms)
        assert Pph == (2 * n + 1) * int(ph[j])
        assert Pk == (n - 1) ** 2 * int(kap[j]) - 6 * int(ph[j])


def test_live_Per_eigenforms_p5_all_types():
    """One representative per (κ,φ) type at p=5."""
    p = 5
    C = paley_conference_prime_power(p).astype(np.int32)
    n = C.shape[0]
    fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap, ph = _kap_phi(C, fours)
    perms = list(permutations(range(4)))
    seen = {}
    for j in range(len(fours)):
        key = (int(kap[j]), int(ph[j]))
        if key not in seen:
            seen[key] = j
    assert len(seen) >= 8  # several κ,φ classes
    for (k, phi), j in seen.items():
        S = fours[j]
        Pph = _per_apply(C, fours, S, ph, perms)
        Pk = _per_apply(C, fours, S, kap, perms)
        assert Pph == (2 * n + 1) * phi
        assert Pk == (n - 1) ** 2 * k - 6 * phi
        # μ_part identity
        a, b = mu_part_coeffs(p)
        Pmp = a * Pk + b * Pph
        assert Pmp == p4_mupart_plus_2phi(p, k, phi)
