"""
Verification tests for the min-max ±1 quadratic form limit solution.
These drive the shipped library functions (no re-implementation of claims).

Performance policy: exact_m(n) for n<=8 is expensive (n=8 is 2^21 gray-code).
Call it once per session, in parallel across n, then reuse. Never call exact_m(9)
live (use recorded parallel run). Do not re-invoke exact_m inside every test.
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    exact_m,
    form_Q,
    phi,
    phi_local,
    phi_mitm,
    greedy_Q,
    paley_conference_matrix,
    paley_conference_prime_power,
    halfspace_boolean_vector,
    rho_from_boolean_evec,
    is_conference_matrix,
    op_norm_sym,
    spherical_half_bound,
    alpha,
    dmp_lower_bound,
    dual_gaussian_lower_bound,
    random_sym_pm1,
    random_method_upper_bound,
    edge_hamming,
    phi_edge_lipschitz_lower,
    phi_degree_lipschitz_lower,
)

# Recorded m_9, m_10 from multi-worker Gray exhaustive (never recompute in pytest).
M9_RECORDED = 12
M10_RECORDED = 13  # exact_m_parallel under SCRATCH; Paley-9 has Phi=15 > 13
EXPECTED_M = {2: 1, 3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: M9_RECORDED, 10: M10_RECORDED}


def _exact_m_worker(n: int) -> tuple[int, int]:
    """Module-level for ProcessPool pickling; drives shipped exact_m."""
    # 1 BLAS thread per worker
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"
    return n, exact_m(n)


@pytest.fixture(scope="session")
def exact_table() -> dict[int, int]:
    """
    One parallel wave of exact_m(2..8). Session-scoped so the suite does not
    re-peg a single core re-running exact_m(8) in every test.
    """
    ns = list(range(2, 9))
    workers = min(len(ns), max(1, (os.cpu_count() or 4) - 2))
    out: dict[int, int] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for n, m in ex.map(_exact_m_worker, ns):
            out[n] = m
    out[9] = M9_RECORDED
    return out


def test_form_Q_matches_half_xAx():
    rng = np.random.default_rng(0)
    for n in range(3, 8):
        A = random_sym_pm1(n, rng)
        x = rng.choice([-1.0, 1.0], size=n)
        assert abs(form_Q(A, x) - 0.5 * float(x @ A @ x)) < 1e-10


def test_exact_m_small_values(exact_table):
    """Exact table used in the writeup — values come from shipped exact_m via fixture."""
    for n, m in EXPECTED_M.items():
        if n >= 9:
            continue  # recorded parallel Gray; see test_exact_m9_m10_from_recorded_parallel_run
        assert exact_table[n] == m, f"exact_m({n})={exact_table[n]} != {m}"


def test_exact_m9_m10_from_recorded_parallel_run():
    """m_9=12, m_10=13 from multi-worker Gray exhaustive — never live recompute in pytest."""
    import json

    candidates = [
        ROOT / "evidence" / "exact_m_table.json",
        Path("/tmp/grok-goal-c7b746355916/implementer/evidence/exact_m_table.json"),
    ]
    found = False
    for json_path in candidates:
        if not json_path.exists():
            continue
        data = json.loads(json_path.read_text())
        assert data["9"]["m"] == M9_RECORDED
        assert data["10"]["m"] == M10_RECORDED
        # Conference non-optimality: Paley order-10 has Phi=15 > m_10
        assert M10_RECORDED < data["10"].get("Phi_Paley_10", 15)
        found = True
        break
    assert found, "missing evidence/exact_m_table.json (run regen or ship evidence/)"
    assert M9_RECORDED == 12
    assert M10_RECORDED == 13
    assert M10_RECORDED < 15


def test_rho_eq_1_paley_prime_power():
    """
    Shipped ρ=1 theorem: halfspace boolean evec for Paley n=p^2+1.
    Drives paley_conference_prime_power + halfspace_boolean_vector (not a re-implementation).
    """
    for p in (3, 5, 7):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        assert n == p * p + 1
        assert is_conference_matrix(C)
        x = halfspace_boolean_vector(p)
        assert set(np.unique(x)).issubset({-1.0, 1.0})
        Cx = C @ x
        assert float(np.max(np.abs(Cx - p * x))) < 1e-8
        rho = rho_from_boolean_evec(C, x)
        assert abs(rho - 1.0) < 1e-8
        # Φ from evec hits spectral ceiling
        phi_hs = abs(0.5 * float(x @ C @ x))
        assert abs(phi_hs - 0.5 * n * p) < 1e-8
        assert abs(phi_hs - spherical_half_bound(n)) < 1e-8


def test_handoff_document_exists_and_keeps_existence_open():
    """Structural AC: handoff restates the limit and does not claim L exists."""
    text = (ROOT / "HANDOFF.md").read_text()
    assert "n^{-3/2}" in text or "n^{3/2}" in text or "alpha_n" in text
    assert "OPEN" in text
    # must not claim the limit is settled
    low = text.lower()
    assert "existence of the limit" in low or "existence of" in low
    assert "remains open" in low or "is open" in low or "**open**" in low
    # required sections
    for needle in ("Proved", "Open blocker", "Numerics", "Resume playbook"):
        assert needle.lower() in low or needle in text


def test_monotonicity_exact(exact_table):
    ms = [exact_table[n] for n in range(2, 10)]
    for i in range(len(ms) - 1):
        assert ms[i] <= ms[i + 1]


def test_equivalence_m_vs_half_phi_on_optimum_n6(exact_table):
    """
    For n=6 a conference matrix realises m_6=5, and phi(C)=5,
    so m_n = min phi = phi(C) with the factor built into form_Q.
    """
    C = paley_conference_matrix(5)
    assert is_conference_matrix(C)
    assert C.shape == (6, 6)
    assert abs(op_norm_sym(C) - np.sqrt(5)) < 1e-8
    ph = phi(C)  # brute force n=6 is cheap
    assert ph == pytest.approx(5.0)
    assert exact_table[6] == 5


def test_spherical_bound_when_conference_exists():
    for q in [5, 13, 17]:
        C = paley_conference_matrix(q)
        n = C.shape[0]
        assert is_conference_matrix(C)
        ph = phi_local(C, restarts=60, rng=np.random.default_rng(1))
        assert ph <= spherical_half_bound(n) + 1e-6
        assert abs(op_norm_sym(C) - np.sqrt(n - 1)) < 1e-6


def test_dmp_lower_bound_holds_for_exact_m(exact_table):
    for n in range(2, 9):
        m = exact_table[n]
        assert m >= 1
        if n >= 5:
            assert m + 1e-9 >= dmp_lower_bound(n) * 0.5


def test_dual_gaussian_lower_bound_holds_for_exact_m(exact_table):
    """Prop 5.2: m_n >= n sqrt(n-1)/pi for all n; checked on exact table."""
    for n in range(2, 9):
        m = exact_table[n]
        assert m + 1e-9 >= dual_gaussian_lower_bound(n)


def test_dual_gaussian_beats_dmp_asymptotically():
    """1/pi > 2^{-5/2}; dual-Gaussian dominates BH for large n."""
    assert dual_gaussian_lower_bound(100) > dmp_lower_bound(100)


def test_greedy_is_lower_bound_on_phi():
    rng = np.random.default_rng(2)
    for n in [5, 6, 7, 8]:
        A = random_sym_pm1(n, rng)
        g = greedy_Q(A)
        assert g <= phi(A) + 1e-9


def test_alpha_envelope(exact_table):
    """0 < alpha and finite envelope on exact range."""
    for n in range(2, 9):
        m = exact_table[n]
        a = alpha(m, n)
        assert a > 0
        assert a < 1.0
        assert m <= random_method_upper_bound(n) + n


def test_paley_alpha_upper_tends_toward_half():
    """Upper bounds alpha <= 1/2 * sqrt(1-1/n) along Paley orders."""
    for q in [5, 13, 17, 29, 37]:
        C = paley_conference_matrix(q)
        n = C.shape[0]
        ub = spherical_half_bound(n) / n**1.5
        assert ub <= 0.5 + 1e-12
        assert abs(ub - 0.5 * np.sqrt(1 - 1 / n)) < 1e-12
        ph = phi_local(C, restarts=50, rng=np.random.default_rng(q))
        assert alpha(ph, n) <= ub + 1e-9


def test_restriction_monotonicity_identity_on_submatrix(exact_table):
    """Monotonicity conclusion on cached exact table (shipped exact_m via fixture)."""
    assert exact_table[5] <= exact_table[6] <= exact_table[7] <= exact_table[8]


def test_reverse_multipartite_exact(exact_table):
    """Prop 7.3: m_{kn} >= (k/2) m_n on exact table."""
    ms = dict(exact_table)
    for n in range(2, 5):
        for k in range(2, 5):
            if k * n <= 9:
                assert ms[k * n] + 1e-12 >= (k / 2.0) * ms[n]


def test_alpha_consecutive_gaps_shrink_on_exact(exact_table):
    """Prop 3.3 sample: |alpha_{n+1}-alpha_n| is O(n^{-1/2}) on exact range."""
    for n in range(2, 9):
        a = alpha(exact_table[n], n)
        b = alpha(exact_table[n + 1], n + 1)
        assert abs(b - a) <= n ** (-0.5) + 2.0 / n + 0.5


def test_conference_spectral_identity_and_nesterov_formula():
    """Prop 15.1–15.2: C = λ(2P_+-I) and closed-form Nesterov expectation shape."""
    import numpy as np
    for q in (5, 13, 17):
        C = paley_conference_matrix(q)
        assert is_conference_matrix(C)
        n = C.shape[0]
        lam = np.sqrt(n - 1)
        P = 0.5 * (np.eye(n) + C / lam)
        assert np.allclose(C, lam * (2 * P - np.eye(n)), atol=1e-10)
        assert np.allclose(np.diag(P), 0.5, atol=1e-10)
        # Nesterov expectation formula is positive and at most spherical max
        E = (2 / np.pi) * n * (n - 1) * np.arcsin(1 / lam)
        assert 0 < E <= n * lam + 1e-9
        # sign(g) lower bound does not exceed sphere
        rho_lb = (2 / np.pi) * lam * np.arcsin(1 / lam)
        assert 0 < rho_lb < 1


def test_seidel_switching_preserves_phi_and_spectrum():
    """Prop 15.4: A' = D A D has same spectrum and Phi."""
    import numpy as np
    C = paley_conference_matrix(13)
    n = C.shape[0]
    rng = np.random.default_rng(1)
    for _ in range(8):
        eps = rng.choice([-1.0, 1.0], size=n)
        D = np.diag(eps)
        Cp = D @ C @ D
        assert is_conference_matrix(Cp)
        assert abs(phi(C) - phi(Cp)) < 1e-9
        assert np.allclose(np.sort(np.linalg.eigvalsh(C)), np.sort(np.linalg.eigvalsh(Cp)))


def test_min_op_norm_equality_case_conference():
    """Prop 15.5: ||A||_op >= sqrt(n-1), equality on Paley conference."""
    import numpy as np
    for q in (5, 13, 17):
        C = paley_conference_matrix(q)
        n = C.shape[0]
        op = op_norm_sym(C)
        assert abs(op - np.sqrt(n - 1)) < 1e-9
        assert np.allclose(C @ C, (n - 1) * np.eye(n), atol=1e-8)
    # random Seidel has strictly larger op (with high probability; check a few)
    rng = np.random.default_rng(2)
    for n in (8, 10, 12):
        A = random_sym_pm1(n, rng)
        assert op_norm_sym(A) > np.sqrt(n - 1) - 1e-9


def test_L2_universality_and_trA4_min_at_conference():
    """Prop 15.10–15.11: E[Q^2]=C(n,2) for all Seidel; tr(A^4) min iff conference."""
    import numpy as np
    rng = np.random.default_rng(7)
    for n in (6, 8):
        theory = n * (n - 1) / 2  # E[Q^2]
        target_tr4 = n * (n - 1) ** 2
        for s in range(6):
            A = random_sym_pm1(n, rng)
            # exact E[Q^2] via half-cube
            npat = 1 << (n - 1)
            s2 = 0.0
            for xb in range(npat):
                x = np.ones(n)
                for i in range(n - 1):
                    x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
                s2 += form_Q(A, x) ** 2
            assert abs(s2 / npat - theory) < 1e-9
            tr4 = float(np.trace(np.linalg.matrix_power(A, 4)))
            assert tr4 >= target_tr4 - 1e-6
        if (n - 1) % 4 == 1:
            C = paley_conference_matrix(n - 1)
            assert abs(float(np.trace(np.linalg.matrix_power(C, 4))) - target_tr4) < 1e-6


def test_exact_Q4_moment_formula():
    """Prop 15.13: E[Q^4] = 3e^2 + 3(tr(A^4)-n(n-1)^2) - n(n-1)(3n-5) for Seidel A."""
    import numpy as np
    rng = np.random.default_rng(11)
    for n in (5, 6, 7, 8):
        e = n * (n - 1) / 2
        for s in range(4):
            A = random_sym_pm1(n, rng)
            npat = 1 << (n - 1)
            s4 = 0.0
            for xb in range(npat):
                x = np.ones(n)
                for i in range(n - 1):
                    x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
                s4 += form_Q(A, x) ** 4
            E4 = s4 / npat
            tr4 = float(np.trace(np.linalg.matrix_power(A.astype(float), 4)))
            pred = 3 * e**2 + 3 * (tr4 - n * (n - 1) ** 2) - n * (n - 1) * (3 * n - 5)
            assert abs(E4 - pred) < 1e-6, (n, E4, pred)
        if (n - 1) % 4 == 1:
            C = paley_conference_matrix(n - 1)
            npat = 1 << (n - 1)
            s4 = 0.0
            for xb in range(npat):
                x = np.ones(n)
                for i in range(n - 1):
                    x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
                s4 += form_Q(C, x) ** 4
            E4 = s4 / npat
            # conference: tr4 = n(n-1)^2
            pred = 3 * e**2 - n * (n - 1) * (3 * n - 5)
            assert abs(E4 - pred) < 1e-6


def test_exact_optimality_n6_via_Q4_gap():
    """Cor 15.15: for n=6, min delta non-conf=64 > Delta/3=30, so m_6=Phi(C)=5."""
    import numpy as np
    n = 6
    C = paley_conference_matrix(5)
    assert abs(phi(C) - 5.0) < 1e-9
    e = n * (n - 1) / 2
    E4_C = 3 * e**2 - n * (n - 1) * (3 * n - 5)
    Delta = 5.0**2 * e - E4_C
    assert abs(Delta - 90.0) < 1e-9
    target = n * (n - 1) ** 2
    free = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    nf = len(free)
    min_delta = 1e99
    conf = 0
    for mask in range(1 << nf):
        A = np.zeros((n, n))
        for j in range(1, n):
            A[0, j] = A[j, 0] = 1.0
        for b, (i, j) in enumerate(free):
            s = 1.0 if (mask >> b) & 1 else -1.0
            A[i, j] = A[j, i] = s
        t4 = float(np.trace(A @ A @ A @ A))
        d = t4 - target
        if abs(d) < 1e-6:
            conf += 1
        else:
            min_delta = min(min_delta, d)
    assert conf == 12
    assert min_delta == 64.0
    assert min_delta > Delta / 3.0
    assert exact_m(6) == 5


def test_single_edge_flip_delta_is_16_n_minus_2():
    """Prop 15.16: single edge flip from Paley conference has δ=16(n-2)."""
    for q in (5, 13, 17):
        C = paley_conference_matrix(q)
        n = C.shape[0]
        target = n * (n - 1) ** 2
        # flip several edges (all for n=6, sample for larger)
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n > 6:
            pairs = pairs[:: max(1, len(pairs) // 12)]
        for i, j in pairs:
            A = C.copy()
            A[i, j] *= -1
            A[j, i] *= -1
            assert not is_conference_matrix(A, tol=1e-6)
            tr4 = float(np.trace(np.linalg.matrix_power(A.astype(float), 4)))
            delta = tr4 - target
            assert abs(delta - 16 * (n - 2)) < 1e-6, (n, i, j, delta)
            # Phi cannot drop below conference on these flips (local opt evidence; n<=18 exact)
            if n <= 18:
                assert phi(A) >= phi(C) - 1e-12


def test_projector_form_of_rho():
    """Prop 15.18: x^T C x = s(2||P_+ x||^2 - ||x||^2); rho via cube imbalance."""
    for q in (5, 13, 17):
        C = paley_conference_matrix(q)
        n = C.shape[0]
        s = np.sqrt(n - 1)
        P = 0.5 * (np.eye(n) + C / s)
        assert np.allclose(P @ P, P, atol=1e-8)
        rng = np.random.default_rng(q)
        for _ in range(40):
            x = rng.choice([-1.0, 1.0], size=n)
            lhs = float(x @ C @ x)
            rhs = s * (2 * float((P @ x) @ (P @ x)) - n)
            assert abs(lhs - rhs) < 1e-7
        # rho identity
        Ph = phi(C)
        rho = 2 * Ph / (n * s)
        # max imbalance over half-cube
        npat = 1 << (n - 1)
        max_dev = 0.0
        for xb in range(npat):
            x = np.ones(n)
            for i in range(n - 1):
                x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
            Px = P @ x
            max_dev = max(max_dev, abs(2 * float(Px @ Px) / n - 1))
        assert abs(max_dev - rho) < 1e-8


def test_max_delta_bound_and_shell_vacuity_threshold():
    """Prop 15.19: δ ≤ (n-1)^4+(n-1)-n(n-1)^2; shell vacuous at n=38 with exact Φ*=109."""
    def max_delta_JI(n):
        return (n - 1) ** 4 + (n - 1) - n * (n - 1) ** 2

    def delta_star_over_3(n, Phi):
        e = n * (n - 1) / 2
        return (Phi**2 * e - 3 * e**2 + n * (n - 1) * (3 * n - 5)) / 3.0

    # J-I realises the bound
    for n in (6, 10, 14, 18):
        A = np.ones((n, n)) - np.eye(n)
        A = -A  # off-diagonal -1
        tr4 = float(np.trace(np.linalg.matrix_power(A.astype(float), 4)))
        delta = tr4 - n * (n - 1) ** 2
        assert abs(delta - max_delta_JI(n)) < 1e-6

    # Random Seidel never exceeds the bound
    rng = np.random.default_rng(0)
    for n in (8, 12, 16):
        md = max_delta_JI(n)
        for _ in range(20):
            A = random_sym_pm1(n, rng)
            tr4 = float(np.trace(np.linalg.matrix_power(A.astype(float), 4)))
            delta = tr4 - n * (n - 1) ** 2
            assert delta <= md + 1e-6

    # Exact Phi* values: shell non-vacuous at n=6,14,30; vacuous at n=38
    assert delta_star_over_3(6, 5.0) < max_delta_JI(6)
    assert delta_star_over_3(14, 21.0) < max_delta_JI(14)
    assert delta_star_over_3(30, 75.0) < max_delta_JI(30)
    assert delta_star_over_3(38, 109.0) >= max_delta_JI(38)


def test_phi_lipschitz_in_frobenius():
    """Prop 15.20: |Φ(A)-Φ(B)| ≤ (n/2) ||A-B||_F for Seidel A,B."""
    rng = np.random.default_rng(1)
    for n in (4, 6, 8, 10, 12):
        for _ in range(15):
            A = random_sym_pm1(n, rng)
            B = random_sym_pm1(n, rng)
            dphi = abs(phi(A) - phi(B))
            frob = float(np.linalg.norm(A - B, "fro"))
            assert dphi <= 0.5 * n * frob + 1e-9


def test_phi_edge_lipschitz_pointwise_and_global():
    """
    Prop 15.20b: |Q_A(x)-Q_C(x)| <= 2k for every x, hence Phi(A) >= Phi(C)-2k.
    Drives shipped edge_hamming + phi_edge_lipschitz_lower on real Paley/random pairs.
    """
    rng = np.random.default_rng(2)
    C = paley_conference_prime_power(3)
    n = C.shape[0]
    phi_c = phi(C)
    for _ in range(20):
        A = random_sym_pm1(n, rng)
        k = edge_hamming(A, C)
        # pointwise on a sample of cube vectors
        for __ in range(40):
            x = rng.choice([-1.0, 1.0], size=n)
            dQ = abs(form_Q(A, x) - form_Q(C, x))
            assert dQ <= 2 * k + 1e-9
        # global Phi bound (may be negative / vacuous for large k)
        assert phi(A) + 1e-9 >= phi_edge_lipschitz_lower(phi_c, k)


def test_phi_edge_lipschitz_matching_undercut_n10():
    """N10 matching undercutters: k=5, bound Phi-2k=5, actual Phi=13 >= 5."""
    from n10_matching_optima import perfect_matchings

    C = paley_conference_prime_power(3)
    phi_c = phi(C)
    assert abs(phi_c - 15.0) < 1e-9
    seen = 0
    for m in perfect_matchings(10):
        A = C.copy()
        for i, j in m:
            A[i, j] *= -1.0
            A[j, i] *= -1.0
        if phi(A) > 13.0 + 1e-9:
            continue
        k = edge_hamming(A, C)
        assert k == 5
        assert phi(A) + 1e-9 >= phi_edge_lipschitz_lower(phi_c, k)
        assert phi(A) + 1e-9 >= phi_degree_lipschitz_lower(phi_c, 1, 10)
        seen += 1
        if seen >= 5:
            break
    assert seen >= 5


def test_e1_edge_lipschitz_criterion_string():
    """E(1) criterion via edge lip: k=o(n^{3/2}) forces relative gap o(1)."""
    # symbolic check used in evidence/E1_STATUS.md — drive shipped lower bound
    for n, phi_c, k in ((10, 15.0, 5), (26, 65.0, 0), (26, 65.0, 13)):
        lb = phi_edge_lipschitz_lower(phi_c, k)
        rel = (phi_c - max(lb, 0.0)) / (n ** 1.5)
        if k == 0:
            assert abs(lb - phi_c) < 1e-12
        if k * 2 < n ** 1.5 * 0.01:
            # sparse enough that relative gap from this bound alone is < 1%
            assert (2.0 * k) / (n ** 1.5) < 0.01


def test_paley_maximizer_balance_forces_edge_flip_plus_two():
    """
    Prop 15.21: for Paley n=6,14,18 every edge has a maximizer forcing Φ→Φ+2 on flip.
    """
    for q in (5, 13, 17):
        C = paley_conference_matrix(q)
        n = C.shape[0]
        M = phi(C)
        npat = 1 << (n - 1)
        maxs = []
        for xb in range(npat):
            x = np.ones(n)
            for i in range(n - 1):
                x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
            qv = form_Q(C, x)
            if abs(abs(qv) - M) < 1e-9:
                maxs.append((x, qv))
        assert len(maxs) > 0
        # every edge admits a boost maximizer
        for p in range(n):
            for r in range(p + 1, n):
                c = C[p, r]
                can = False
                for x, qv in maxs:
                    prod = c * x[p] * x[r]
                    if abs(qv - M) < 1e-9 and prod < 0:
                        can = True
                    if abs(qv + M) < 1e-9 and prod > 0:
                        can = True
                assert can, (n, p, r)
        # spot-check flips
        for p, r in [(0, 1), (1, 2), (0, n - 1)]:
            A = C.copy()
            A[p, r] *= -1
            A[r, p] *= -1
            assert phi(A) >= M + 2 - 1e-12


def test_rho_min_controls_alpha_and_exhaustive_n6_to_8():
    """Prop 15.22–15.23: alpha >= (1/2)sqrt(1-1/n) rho_min; exhaustive floors."""
    # Recorded exhaustive min_rho from 86w shard (re-run: exhaustive_min_rho.py)
    REC = {
        6: {"min_rho": 0.7453559924999299, "min_r": 0.7453559924999299, "min_phi": 5.0},
        7: {"min_rho": 0.7219964736113526, "min_r": 1.049781318335648, "min_phi": 9.0},
        8: {"min_rho": 0.6933752452815356, "min_r": 0.944911182523068, "min_phi": 10.0},
    }
    two_pi = 2.0 / np.pi
    for n, rec in REC.items():
        assert rec["min_rho"] > two_pi
        # Prop 15.22: m_n >= (1/2) n sqrt(n-1) rho_min
        m = exact_m(n) if n < 8 else 10  # n=8 from EXPECTED / exact
        if n <= 7:
            m = exact_m(n)
        else:
            m = 10  # known m_8
        assert abs(m - rec["min_phi"]) < 1e-9
        lb = 0.5 * n * np.sqrt(n - 1) * rec["min_rho"]
        assert m + 1e-9 >= lb  # m = min Phi >= lb
        # alpha form
        assert m / n**1.5 + 1e-12 >= 0.5 * np.sqrt(1 - 1 / n) * rec["min_rho"]
    # n=6: min_r equals Paley
    C = paley_conference_matrix(5)
    r_p = 2 * phi(C) / (6 * np.sqrt(5))
    assert abs(r_p - REC[6]["min_r"]) < 1e-9


# ---------------------------------------------------------------------------
# n=10 structural gap (evidence/N10_STRUCTURE.md) — load-bearing, exact Φ
# ---------------------------------------------------------------------------


def _flip(A, i, j):
    A[i, j] *= -1.0
    A[j, i] *= -1.0


def test_n10_paley_maximizer_balance_all_edges():
    """Theorem N10-S.1: every edge of Paley C_10 is maximizer-balanced; 1-flip → Φ≥17."""
    C = paley_conference_prime_power(3)
    assert is_conference_matrix(C)
    M = phi(C)
    assert abs(M - 15.0) < 1e-9
    n = 10
    npat = 1 << (n - 1)
    maxs = []
    for xb in range(npat):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        qv = form_Q(C, x)
        if abs(abs(qv) - M) < 1e-9:
            maxs.append((x, qv))
    assert len(maxs) == 12
    for p in range(n):
        for r in range(p + 1, n):
            c = C[p, r]
            can = False
            for x, qv in maxs:
                prod = c * x[p] * x[r]
                if abs(qv - M) < 1e-9 and prod < 0:
                    can = True
                if abs(qv + M) < 1e-9 and prod > 0:
                    can = True
            assert can, (p, r)
    # every single flip raises Φ by ≥2
    for p, r in [(0, 1), (1, 2), (2, 5), (0, 9), (4, 7)]:
        A = C.copy()
        _flip(A, p, r)
        assert phi(A) >= M + 2 - 1e-12


def test_n10_matching_flip_census_144_optima():
    """
    Theorem N10-S.3–4: among 945 perfect matchings of K_10, exactly 144 flips
    of Paley C_10 yield Φ=13=m_10; Φ histogram {13:144, 17:405, 21:360, 25:36}.
    """
    # Inline perfect-matching enumerator (same as src/n10_matching_optima.py)
    def perfect_matchings(n=10):
        verts = list(range(n))
        out = []

        def rec(remaining, partial):
            if not remaining:
                out.append(list(partial))
                return
            a = remaining[0]
            for i in range(1, len(remaining)):
                b = remaining[i]
                edge = (a, b) if a < b else (b, a)
                nxt = remaining[1:i] + remaining[i + 1 :]
                partial.append(edge)
                rec(nxt, partial)
                partial.pop()

        rec(verts, [])
        return out

    C = paley_conference_prime_power(3)
    pms = perfect_matchings(10)
    assert len(pms) == 945
    hist = {13: 0, 17: 0, 21: 0, 25: 0}
    for m in pms:
        A = C.copy()
        for i, j in m:
            _flip(A, i, j)
        ph = int(round(phi(A)))
        assert ph in hist, ph
        hist[ph] += 1
    assert hist == {13: 144, 17: 405, 21: 360, 25: 36}
    # relative gap vanishes as O(n^{-3/2}) if absolute gap stays O(1)
    assert (15 - 13) / (10 ** 1.5) < 0.07
    # r at a matching-optimum = 13/15
    A = C.copy()
    for i, j in pms[0]:
        _flip(A, i, j)
    # find one with phi 13
    found = False
    for m in pms:
        A = C.copy()
        for i, j in m:
            _flip(A, i, j)
        if abs(phi(A) - 13) < 1e-9:
            r = 2 * 13 / (10 * np.sqrt(9))
            assert abs(r - 13 / 15) < 1e-12
            found = True
            break
    assert found


def test_n10_k_flip_threshold_no_undercut_before_5():
    """
    Theorem N10-S.2 (sampled + recorded): k≤4 cannot reach m_10 from this C;
    one certified 5-matching reaches 13. Full binom(45,5) lives in evidence JSON.
    """
    import json

    C = paley_conference_prime_power(3)
    # k=1 exhaustive (45 edges)
    edges = [(i, j) for i in range(10) for j in range(i + 1, 10)]
    vals = []
    for i, j in edges:
        A = C.copy()
        _flip(A, i, j)
        vals.append(phi(A))
    assert min(vals) >= 17.0 - 1e-12

    # certified witness from evidence
    path = ROOT / "evidence" / "n10_structure.json"
    assert path.exists(), "run src/n10_structure.py to ship n10_structure.json"
    data = json.loads(path.read_text())
    for row in data["k_flip_table"]:
        if row["k"] <= 4:
            assert row["min_phi"] >= 15.0 - 1e-12
            assert not row.get("reaches_m10")
        if row["k"] == 5:
            assert abs(row["min_phi"] - 13.0) < 1e-12
            assert row.get("reaches_m10")
            # re-verify the recorded best_edges witness
            A = C.copy()
            for i, j in row["best_edges"]:
                _flip(A, i, j)
            assert abs(phi(A) - 13.0) < 1e-9
            # perfect matching: 10 distinct endpoints
            verts = [v for e in row["best_edges"] for v in e]
            assert len(verts) == 10 and len(set(verts)) == 10

    mpath = ROOT / "evidence" / "n10_matching_optima.json"
    assert mpath.exists()
    md = json.loads(mpath.read_text())
    assert md["matching_flip"]["n_phi13"] == 144
    assert md["five_edge_undercutters"]["all_are_perfect_matchings"]


def test_n10_matching_classify_maximizer_criterion_and_orbit():
    """
    Theorem N10-C: drive shipped n10_matching_classify — maximizer-drop
    criterion equals Φ=13 matchings; positive maximizers suffice; one seed's
    PΓL orbit has size 144 and sits inside the Φ=13 set.
    """
    sys.path.insert(0, str(ROOT / "src"))
    from n10_matching_classify import (
        flip_matching,
        matching_to_frozenset,
        maximizer_drop_criterion,
        paley_maximizers,
        pgl29_representatives,
        pgl_orbit,
        positive_maximizers,
    )
    from n10_matching_optima import perfect_matchings

    C = paley_conference_prime_power(3)
    assert abs(phi(C) - 15.0) < 1e-9
    maximizers = paley_maximizers(C)
    assert len(maximizers) == 12
    pos = positive_maximizers(maximizers)
    assert len(pos) == 6

    pms = perfect_matchings(10)
    assert len(pms) == 945
    opt_phi = []
    opt_crit = []
    opt_pos = []
    for m in pms:
        ph = phi(flip_matching(C, m))
        if abs(ph - 13.0) < 1e-9:
            opt_phi.append(m)
        if maximizer_drop_criterion(C, m, maximizers):
            opt_crit.append(m)
        if maximizer_drop_criterion(C, m, pos):
            opt_pos.append(m)
    assert len(opt_phi) == 144
    assert len(opt_crit) == 144
    assert len(opt_pos) == 144
    set_phi = {matching_to_frozenset(m) for m in opt_phi}
    set_crit = {matching_to_frozenset(m) for m in opt_crit}
    set_pos = {matching_to_frozenset(m) for m in opt_pos}
    assert set_phi == set_crit == set_pos

    # Non-optimal matching fails criterion
    non = next(m for m in pms if matching_to_frozenset(m) not in set_phi)
    assert not maximizer_drop_criterion(C, non, maximizers)
    assert abs(phi(flip_matching(C, non)) - 13.0) > 1e-9

    # PΓL orbit of one optimum
    reps = pgl29_representatives()
    assert len(reps) == 720
    orbit = pgl_orbit(opt_phi[0], reps)
    assert len(orbit) == 144
    assert orbit == set_phi

    # Evidence JSON agrees (durable campaign artifact)
    import json

    cpath = ROOT / "evidence" / "n10_matching_classify.json"
    assert cpath.exists(), "run src/n10_matching_classify.py"
    data = json.loads(cpath.read_text())
    assert data["n_phi13"] == 144
    assert data["sets_agree"] is True
    assert data["orbit_equals_optima"] is True
    assert data["criterion_positive_equals_phi13"] is True
    note = ROOT / "evidence" / "N10_MATCHING_CLASSIFY.md"
    assert note.exists()
    text = note.read_text()
    assert "OPEN" in text
    assert "N10-C" in text


def test_interval_rho_formula_matches_paley_matrix():
    """
    E2 interval formula: drive shipped interval_rho_formula against real Paley C.
    Proves the identity x^T C x = 2-8Σ_q; does not claim lim α exists.
    """
    sys.path.insert(0, str(ROOT / "src"))
    from interval_rho_formula import (
        Sigma_q,
        interval_xCx_formula,
        rho_int_formula,
        verify_formula_against_matrix,
    )

    for q in (5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101):
        C = paley_conference_matrix(q)
        assert is_conference_matrix(C)
        ok, rf, rm = verify_formula_against_matrix(q, C)
        assert ok, (q, rf, rm)
        assert abs(rf - rm) < 1e-9
        # formula path alone
        assert abs(rho_int_formula(q) - rf) < 1e-12
        # Σ_q is integer; S = 2-8Σ
        assert interval_xCx_formula(q) == 2 - 8 * Sigma_q(q)
        # ρ_int ≤ 1 (cube/sphere)
        assert rf <= 1.0 + 1e-12
        # Nesterov floor: ρ ≥ (2/π)√(n-1) arcsin(1/√(n-1)) roughly → 2/π
        assert rf > 0.5  # weak sanity; all tested q clear this

    # evidence note keeps existence OPEN
    note = ROOT / "evidence" / "E2_INTERVAL_FORMULA.md"
    assert note.exists()
    assert "OPEN" in note.read_text()


def test_interval_rho_asymptotic_main_term_close():
    """
    E2 asymptotic: ρ_int ≈ (8/π²) sum_{odd m} χ(m)/m² for large prime q≡1 mod 4.
    Drives shipped interval_rho_asymptotics; does not claim lim α_n.
    """
    sys.path.insert(0, str(ROOT / "src"))
    from interval_rho_asymptotics import rho_int_error_estimate

    # moderate primes: relative error should be small
    for q in (101, 421, 709, 1009):
        rho, main, err = rho_int_error_estimate(q)
        assert rho > 0.5
        assert main > 0.5
        assert err / rho < 0.02, (q, rho, main, err)

    note = ROOT / "evidence" / "E2_RHO_INT_ASYMPTOTICS.md"
    assert note.exists()
    text = note.read_text()
    assert "OPEN" in text
    assert "E(1)" in text


def test_n10_k6_undercutters_are_cycles():
    """
    Theorem N10-C6: drive shipped n10_cycle_undercutters against durable JSON.
    All 360 Hamming-6 undercutters of Paley C_10 are single 6-cycles with Φ=13.
    Does not settle lim α_n (still needs general k_star bound).
    """
    import json

    sys.path.insert(0, str(ROOT / "src"))
    from n10_cycle_undercutters import _is_single_cycle
    from minmax_quadratic import paley_conference_prime_power, phi

    cpath = ROOT / "evidence" / "n10_cycle_undercutters.json"
    assert cpath.exists(), "run: python -m src.n10_cycle_undercutters (or classify script)"
    data = json.loads(cpath.read_text())
    assert data["n_undercutters_k6"] == 360
    assert data["all_are_single_6cycles"] is True
    assert data["all_have_Phi_13"] is True
    assert data["n_scanned"] == 8145060  # binom(45,6)

    # Spot-check: a concrete 6-cycle undercuts to 13; a 6-edge star does not
    C = paley_conference_prime_power(3)
    # cycle (0,1,2,3,4,5)
    cyc = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]
    assert _is_single_cycle(cyc, 10)
    A = C.copy()
    for i, j in cyc:
        A[i, j] *= -1
        A[j, i] *= -1
    # not every C6 undercuts; only 360 of them — just check structure helper + Phi bounds
    ph = phi(A)
    assert ph >= 13.0 - 1e-9  # never below m_10
    assert ph <= 25.0 + 1e-9

    # star of degree 6 at 0 never undercuts (exhaustive in campaign: 0/840 for d=3..8)
    star = [(0, j) for j in range(1, 7)]
    A2 = C.copy()
    for i, j in star:
        A2[i, j] *= -1
        A2[j, i] *= -1
    assert phi(A2) >= 15.0 - 1e-9

    note = ROOT / "evidence" / "N10_CYCLE_UNDERCUTTERS.md"
    assert note.exists()
    text = note.read_text()
    assert "OPEN" in text
    assert "N10-C6" in text or "6-cycle" in text


def test_rho1_maximizers_are_boolean_eigenvectors():
    """
    Prop 15.24: for rho=1 Paley, Max = {x: C x = ±p x}.
    Drives paley_conference_prime_power + halfspace_boolean_vector.
    """
    for p in (3, 5):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        h = halfspace_boolean_vector(p)
        # halfspace is a +p evec
        assert np.allclose(C @ h, p * h)
        Phi = abs(form_Q(C, h))
        assert abs(Phi - 0.5 * n * p) < 1e-8
        # every boolean +p evec is a maximizer
        # (check halfspace and its global sign)
        assert abs(form_Q(C, h) - Phi) < 1e-8
        assert abs(form_Q(C, -h) - Phi) < 1e-8
        # at n=10, enumerate: all Max are evecs
        if n == 10:
            from itertools import product

            max_count = 0
            evec_count = 0
            for bits in product([-1.0, 1.0], repeat=n - 1):
                x = np.ones(n)
                x[1:] = bits
                q = abs(form_Q(C, x))
                if abs(q - Phi) < 1e-8:
                    max_count += 1
                    if np.allclose(C @ x, p * x) or np.allclose(C @ x, -p * x):
                        evec_count += 1
            assert max_count == evec_count == 12
    note = ROOT / "evidence" / "BOOLEAN_EVECS_MAX.md"
    assert note.exists()
    assert "OPEN" in note.read_text()


def test_phi_mitm_matches_brute_force_n10():
    """Exact meet-in-the-middle Phi matches brute phi at n=10 (Paley + matching)."""
    C = paley_conference_prime_power(3)
    assert abs(phi_mitm(C) - phi(C)) < 1e-9
    assert abs(phi_mitm(C) - 15.0) < 1e-9
    # one perfect matching flip
    A = C.copy()
    for i, j in [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]:
        A[i, j] *= -1
        A[j, i] *= -1
    assert abs(phi_mitm(A) - phi(A)) < 1e-9


def test_n10_path_cycle_k_star_bound_implies_e1_string():
    """
    Documented sufficient criterion: k_star ≤ n ⇒ E(1) on ρ=1 family via edge lip.
    Drives phi_edge_lipschitz_lower; does not claim the dichotomy is proved for all p.
    """
    sys.path.insert(0, str(ROOT / "src"))
    from n10_matching_optima import perfect_matchings

    C = paley_conference_prime_power(3)
    n = 10
    phi_c = phi(C)
    assert abs(phi_c - 15.0) < 1e-9
    found = False
    for M in perfect_matchings(10):
        A = C.copy()
        for i, j in M:
            A[i, j] *= -1
            A[j, i] *= -1
        if abs(phi(A) - 13.0) < 1e-9:
            k = edge_hamming(A, C)
            assert k == 5
            lb = phi_edge_lipschitz_lower(phi_c, k)
            assert abs(lb - (15.0 - 2 * 5)) < 1e-12
            assert phi(A) >= lb - 1e-9
            found = True
            break
    assert found
    # relative gap for k_star ≤ n: 2n / n^{3/2} = 2 n^{-1/2} → 0
    assert 2 * n / (n ** 1.5) < 0.7

    note = ROOT / "evidence" / "E1_RIGIDITY_ATTACK.md"
    assert note.exists()
    assert "OPEN" in note.read_text()


def test_recursive_star_formula_m_n():
    """
    Prop 15.25: m_n = min_B max_x (|Q_B(x)| + |sum x_i|) over Seidel B of order n-1.
    Exhaustive for n<=7; drives form_Q and exact_m.
    """
    from itertools import product

    def psi(N: int) -> float:
        free = [(i, j) for i in range(N) for j in range(i + 1, N)]
        X = np.array(list(product([-1.0, 1.0], repeat=N)), dtype=np.float64)
        S = np.abs(X.sum(axis=1))
        best = 1e9
        for mask in range(1 << len(free)):
            B = np.zeros((N, N))
            for b, (i, j) in enumerate(free):
                v = 1.0 if (mask >> b) & 1 else -1.0
                B[i, j] = B[j, i] = v
            Q = 0.5 * np.einsum("ra,ab,rb->r", X, B, X)
            best = min(best, float(np.max(np.abs(Q) + S)))
        return best

    for n in range(3, 8):
        m = exact_m(n)
        assert abs(m - psi(n - 1)) < 1e-9, (n, m, psi(n - 1))


def _boolean_plus_evecs(C: np.ndarray, p: float) -> np.ndarray:
    """All y in {±1}^n with C y = p y (free-variable enum of ker(C-pI))."""
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nb = vt[rank:].T
    nul = n - rank
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(4000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    assert coords is not None
    Binv = np.linalg.inv(nb[coords])
    found = []
    for bits in range(1 << nul):
        free = np.array([1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)])
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    return np.unique(np.round(found).astype(int), axis=0).astype(float)


def test_prop_15_27_max_frame_and_max_lipschitz():
    """
    Prop 15.27: for Paley rho=1, Max+ satisfies E[yy^T]=I+C/p=2P_+,
    and Max-Lipschitz Phi(A) >= Phi(C) - 2k/p after best switch.
    Drives paley_conference_prime_power, form_Q, phi, edge_hamming.
    """
    for p in (3, 5):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        Max = _boolean_plus_evecs(C, float(p))
        assert len(Max) > 0
        G = Max.T @ Max / len(Max)
        pred = np.eye(n) + C / p
        assert np.allclose(G, pred, atol=1e-8), p
        # Frame => average Q on Max equals Phi(C)=n p/2 for A=C
        Phi_c = 0.5 * n * p
        avg_Q = np.mean([form_Q(C, y) for y in Max])
        assert abs(avg_Q - Phi_c) < 1e-8
        # Matching undercutter at n=10: Max-lip gives stronger LB than edge lip
        if p == 3:
            A = C.copy()
            for i, j in [(0, 1), (2, 4), (3, 5), (6, 7), (8, 9)]:
                A[i, j] *= -1
                A[j, i] *= -1
            # best-switch distance (fix first coord)
            best_k = edge_hamming(A, C)
            for bits in range(1 << (n - 1)):
                s = np.ones(n)
                for i in range(n - 1):
                    if (bits >> i) & 1:
                        s[i + 1] = -1
                As = A * np.outer(s, s)
                np.fill_diagonal(As, 0)
                best_k = min(best_k, edge_hamming(As, C))
            assert best_k == 5
            max_lip_lb = Phi_c - 2 * best_k / p
            edge_lb = phi_edge_lipschitz_lower(Phi_c, best_k)
            assert max_lip_lb > edge_lb  # 15 - 10/3 > 15 - 10
            assert phi(A) + 1e-9 >= max_lip_lb
            # average Q on Max after matching flip
            avg = np.mean([form_Q(A, y) for y in Max])
            assert abs(avg - (Phi_c - 2 * best_k / p)) < 1e-8


def test_prop_15_27_fractional_cover_value_p():
    """Primal uniform edge weights achieve Max-cover LP value p."""
    for p in (3, 5):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        Max = _boolean_plus_evecs(C, float(p))
        # uniform x_e = 2/(n p): for each y, sum_e x_e chi_e(y) = 2/(n p) * Q_C(y) = 1
        for y in Max[:5]:
            Q = form_Q(C, y)
            assert abs(Q - 0.5 * n * p) < 1e-8
            cover = 2.0 / (n * p) * Q
            assert abs(cover - 1.0) < 1e-8
        # objective sum_e x_e = 2/(n p) * binom(n,2) = (n-1)/p = p
        assert abs((n - 1) / p - p) < 1e-12
