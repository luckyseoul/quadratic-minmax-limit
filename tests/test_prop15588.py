"""Prop 15.588 — profile classification + floor-as-class-function.

Load-bearing: every check must FAIL if the stated identity is wrong.
This prop must NOT close any leftover.
"""
from __future__ import annotations

import numpy as np
import pytest

import e1_gmin_m4_prop15588 as M

PRIMES = (5, 7)


# ------------------------------------------------------------------ Part I


def test_A_flat_marginals():
    r = M.theorem_A_flat_marginals(PRIMES)
    assert r["proved"], r
    for p in PRIMES:
        rec = r["per_prime"][p]
        assert rec["square"] != rec["nonsquare"]


def test_B_profile_bijection():
    r = M.theorem_B_profile_bijection(PRIMES)
    assert r["proved"], r
    assert r["per_prime"][5]["N"] == 260
    assert r["per_prime"][7]["N"] == 11452


def test_C_degree_bound():
    r = M.theorem_C_degree_bound(PRIMES)
    assert r["proved"], r
    assert r["per_prime"][7]["max_degree_per_k"][4] == 2
    assert r["per_prime"][7]["max_degree_per_k"][3] == 1


def test_D_k0_k2_empty():
    r = M.theorem_D_k2_empty(PRIMES)
    assert r["proved"], r
    for p in PRIMES:
        ks = r["per_prime"][p]["k_values"]
        assert 0 not in ks and 2 not in ks


def test_E_strata_counts_match_closed_forms():
    r = M.theorem_E_strata_counts(PRIMES)
    assert r["proved"], r
    assert r["per_prime"][5]["counts_eps_plus"] == {1: 30, 3: 100}
    assert r["per_prime"][7]["counts_eps_plus"] == {1: 140, 3: 1176, 4: 4410}
    # the repo's unclassified "full" family at p=7 is exactly k=4 = 90q
    assert r["per_prime"][7]["n_full_is_k_ge_4"] == 90 * 49


def test_F_translation_gauge():
    r = M.theorem_F_translation_gauge((7,))
    assert r["proved"], r
    assert r["per_prime"][7]["pure_reps"] == 90
    assert r["per_prime"][7]["one_rep_per_translation_class"] is True


# ----------------------------------------------------------------- Part II


def test_G_class_function_identity():
    r = M.theorem_G_class_function(PRIMES, sample=12)
    assert r["proved"], r
    for p in PRIMES:
        assert r["per_prime"][p]["max_abs_error"] < 1e-6


def test_H_eigenvalues_recovered_from_characters():
    h = M.corollary_H_eigenvalues_from_characters(5)
    assert h["proved"], h
    assert h["group_order"] == 15600
    assert h["n_classes"] == 30
    got = sorted(h["lambdas_from_characters"])
    want = sorted([80 / 13, 144 / 13, 176 / 13])
    assert np.allclose(got, want, atol=1e-9)


# ---------------------------------------------------------- fail-when-wrong


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_wrong_direction_class_breaks_reconstruction(p):
    """Using the FLAT class as the profile coordinates must destroy the
    reconstruction identity p*y = sum sigma_j(t_j) - (m-1) eps."""
    Y = M.maxplus(p)
    sq, nsq = M.directions(p)
    q = p * p
    T = [t for t, f in nsq]  # deliberately the flat class
    m = len(T)
    recon = np.zeros((len(Y), q), dtype=np.int64)
    for j in range(m):
        LS = np.stack([Y[:, 1:][:, T[j] == s].sum(axis=1) for s in range(p)], 1)
        recon += LS[:, T[j]]
    assert not (recon - (m - 1) * Y[:, [0]] == p * Y[:, 1:]).all()


@pytest.mark.parametrize("p", PRIMES)
@pytest.mark.parametrize("delta", [-2, 2])
def test_FWW_perturbing_the_offset_breaks_reconstruction(p, delta):
    """(m-1) is the only offset that works."""
    Y = M.maxplus(p)
    sq, _ = M.directions(p)
    q = p * p
    T = [t for t, f in sq]
    m = len(T)
    recon = np.zeros((len(Y), q), dtype=np.int64)
    P = M.profiles_of(p, Y)
    for j in range(m):
        recon += P[:, j, :][:, T[j]]
    assert not (recon - (m - 1 + delta) * Y[:, [0]] == p * Y[:, 1:]).all()


@pytest.mark.parametrize("bad_k", [1, 2, 3])
def test_FWW_degree_bound_is_tight_at_p7(bad_k):
    """k=4 profiles really do reach degree 2, so a bound of k-3 or lower is
    false: some k-active vector must violate deg <= bad_k - 2."""
    p = 7
    r = M.theorem_C_degree_bound((p,))
    worst = r["per_prime"][p]["max_degree_per_k"]
    assert worst[4] > max(0, bad_k - 2)


@pytest.mark.parametrize("p", PRIMES)
def test_FWW_n1d_closed_form_perturbation(p):
    """m*C(p,m) is the 1D count; neighbours are wrong."""
    r = M.theorem_E_strata_counts((p,))
    got = r["per_prime"][p]["counts_eps_plus"][1]
    assert got == M.n_1d_closed(p)
    assert got != M.n_1d_closed(p) + 1
    m = (p + 1) // 2
    from math import comb

    assert got != (m + 1) * comb(p, m)
    # note C(p, m-1) == C(p, m) identically here (p - m = m - 1), so perturb
    # the binomial by two instead of one
    assert got != m * comb(p, m - 2)
    assert got != (m - 1) * comb(p, m)


@pytest.mark.parametrize("p", PRIMES)
@pytest.mark.parametrize("bad", [-1, 1, "zero"])
def test_FWW_T_offset_breaks_class_function(p, bad):
    """T(g) = (1/N) sum <y,gy>^2 - 2n.  Any other constant breaks
    tr(Phi pi(g)) = T(g) on a non-identity element."""
    C = M.paley_conference_prime_power(p)
    n = C.shape[0]
    Y = M.maxplus(p).astype(np.float64)
    Bs = M.z_basis(p)
    Phi = M.phi_matrix(p)
    gens = M.aut_generators(p)
    g = gens[0]
    Mg = M.signed_matrix(g[0], g[1], n)
    ip = np.einsum("ai,ai->a", Y, (Mg @ Y.T).T)
    off = 0.0 if bad == "zero" else 2.0 * n + bad
    T_bad = float((ip ** 2).mean() - off)
    Pi = np.einsum("tij,uij->tu", Bs,
                   np.einsum("ik,ukl,jl->uij", Mg, Bs, Mg, optimize=True),
                   optimize=True)
    assert abs(T_bad - float(np.trace(Phi @ Pi))) > 1e-6


def test_FWW_character_recovery_needs_the_right_weights():
    """lambda_k = <T,chi_k>/<chi_k,chi_k> with CLASS-SIZE weights; dropping the
    weights must break the recovery."""
    p = 5
    h = M.corollary_H_eigenvalues_from_characters(p)
    assert h["max_recovery_error"] < 1e-9
    # unweighted recovery on the same data must not reproduce the spectrum
    import e1_gmin_m4_prop15588 as MM

    C = MM.paley_conference_prime_power(p)
    n = C.shape[0]
    Y = MM.maxplus(p).astype(np.float64)
    Bs = MM.z_basis(p)
    Phi = MM.phi_matrix(p)
    lam, V = np.linalg.eigh(Phi)
    gens = MM.aut_generators(p)
    rng = np.random.default_rng(1)
    idg = (np.arange(n, dtype=np.int64), np.ones(n, dtype=np.int64))
    Ts, chis0 = [], []
    groups = [np.where(np.abs(lam - v) < 1e-7)[0] for v in sorted(set(np.round(lam, 7)))]
    Ps = [V[:, g] @ V[:, g].T for g in groups]
    for _ in range(12):
        g = idg
        for _ in range(int(rng.integers(1, 6))):
            g = MM._compose(gens[int(rng.integers(0, len(gens)))], g)
        Mg = MM.signed_matrix(g[0], g[1], n)
        ip = np.einsum("ai,ai->a", Y, (Mg @ Y.T).T)
        Ts.append(float((ip ** 2).mean() - 2 * n))
        Pi = np.einsum("tij,uij->tu", Bs,
                       np.einsum("ik,ukl,jl->uij", Mg, Bs, Mg, optimize=True),
                       optimize=True)
        chis0.append(float(np.trace(Ps[0] @ Pi)))
    num = sum(t * c for t, c in zip(Ts, chis0))
    den = sum(c * c for c in chis0)
    assert abs(num / den - min(lam)) > 1e-6


# -------------------------------------------------------------------- flags


def test_no_leftover_flag_is_flipped():
    assert M.leftover_flags_unchanged()
    assert M.phi_F_ge_6_proved_general_via_15588() is False


def test_floor_stays_open_in_this_prop():
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    assert phi_F_ge_6_proved_general() is False
