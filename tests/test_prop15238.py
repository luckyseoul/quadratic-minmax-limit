"""Tests for Prop 15.238 — Per μ_part identity; eigenforms certified small p; residual open."""
from __future__ import annotations

import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15238 import (  # noqa: E402
    hinge_status_238,
    mu_part_coeffs,
    p4_mupart_plus_2phi,
    per_kappa_form,
    per_mupart_via_eigen,
    per_phi_form,
    residual_i_closed_via_238,
    theorem_delta_Per_eigenvalue,
    theorem_per_eigen_certified_small_p,
    theorem_per_mupart_identity,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_per_mupart_identity_proved():
    assert theorem_per_mupart_identity()["proved"] is True


def test_delta_structure_no_bound():
    t = theorem_delta_Per_eigenvalue()
    assert t["proved_structure"] is True
    assert t["bound_proved_general"] is False


def test_per_eigen_not_general():
    assert theorem_per_eigen_certified_small_p()["proved"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_238() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_238()
    assert h["per_mupart_identity"] is True
    assert h["per_eigen_general"] is False
    assert h["residual_i_closed_via_238"] is False


def test_phi_coeff_identity_many_primes():
    """Drive mu_part_coeffs: −6a+b(2n+1)=p⁴b+2 for primes p≥5."""
    for p in range(5, 120):
        if not is_prime(p):
            continue
        a, b = mu_part_coeffs(p)
        n = p * p + 1
        assert -6 * a + b * (2 * n + 1) == (p**4) * b + 2


def test_scalar_FE_via_eigenforms():
    """per_mupart_via_eigen equals p4_mupart_plus_2phi on sample labels."""
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        for kappa in (-1, 1):
            for phi in (-2 * (p - 2), -2, 2, 2 * (p - 2)):
                assert per_mupart_via_eigen(p, kappa, phi) == p4_mupart_plus_2phi(
                    p, kappa, phi
                )


def _kappa_phi(C, fours):
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


def test_paley_p5_per_eigenforms_live():
    """Live Per@κ, Per@φ on each |κ|=1 type at Paley p=5 — drives permanent."""
    p = 5
    C = paley_conference_prime_power(p).astype(np.int32)
    n = C.shape[0]
    fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap, ph = _kappa_phi(C, fours)
    perms = list(permutations(range(4)))
    k1 = np.where(np.abs(kap) == 1)[0]
    seen = {}
    for j in k1:
        key = (int(kap[j]), int(ph[j]))
        if key not in seen:
            seen[key] = int(j)
    assert len(seen) == 4
    for (k, phi), j in seen.items():
        S = fours[j]
        Pk = Pph = 0
        for sig in perms:
            pr = np.ones(len(fours), dtype=np.int64)
            for i in range(4):
                pr *= C[S[i], fours[:, sig[i]]]
            Pk += int(pr @ kap)
            Pph += int(pr @ ph)
        assert Pk == per_kappa_form(p, k, phi)
        assert Pph == per_phi_form(p, phi)
        # identity A in live arithmetic
        a, b = mu_part_coeffs(p)
        Pmp = a * Pk + b * Pph
        assert Pmp == p4_mupart_plus_2phi(p, k, phi)


def test_paley_p3_per_eigenforms_live():
    """Live Per eigenforms at p=3 (full |κ|=1 type set)."""
    p = 3
    C = paley_conference_prime_power(p).astype(np.int32)
    n = C.shape[0]
    fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap, ph = _kappa_phi(C, fours)
    perms = list(permutations(range(4)))
    k1 = np.where(np.abs(kap) == 1)[0]
    seen = {}
    for j in k1:
        key = (int(kap[j]), int(ph[j]))
        if key not in seen:
            seen[key] = int(j)
    for (k, phi), j in seen.items():
        S = fours[j]
        Pk = Pph = 0
        for sig in perms:
            pr = np.ones(len(fours), dtype=np.int64)
            for i in range(4):
                pr *= C[S[i], fours[:, sig[i]]]
            Pk += int(pr @ kap)
            Pph += int(pr @ ph)
        assert Pk == per_kappa_form(p, k, phi)
        assert Pph == per_phi_form(p, phi)
