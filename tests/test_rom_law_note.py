"""Checks for evidence/NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md.

The first-moment lemma and the vector product inequality are proved in the
note.  These tests exercise their computational steps, the small exact
values, and the stored receipts.  Nothing here bounds any m_n beyond the
shipped table, and nothing here touches the open limit.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import betainc

from minmax_quadratic import phi, spherical_half_bound

ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / "evidence" / "rom_law_20261111"


def _paley(q: int) -> np.ndarray:
    path = ROOT / "scripts" / "rom_law" / "paley.py"
    spec = importlib.util.spec_from_file_location("rom_law_paley", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.paley(q)


def _phi_fast(A: np.ndarray) -> float:
    """Exact max |Q_A| over the projective cube, vectorised (n <= 20)."""
    n = A.shape[0]
    codes = np.arange(1 << (n - 1), dtype=np.int64)
    X = np.ones((codes.size, n))
    for i in range(n - 1):
        X[:, i + 1] = 1.0 - 2.0 * ((codes >> i) & 1)
    return float(np.max(np.abs(0.5 * np.einsum("ia,ab,ib->i", X, A, X))))


def _lift(C: np.ndarray, d: int) -> np.ndarray:
    B = C + d * np.eye(C.shape[0])
    return np.block([[C, B], [B, -C]])


def _spectral_ceiling(n: int) -> float:
    return 0.5 * n * math.sqrt(n - 1) if n % 4 == 2 else 0.5 * n**1.5


def test_first_moment_threshold_is_sqrt15_over_4():
    eps_star = (1.0 - math.sqrt(15.0) / 4.0) / 2.0
    assert abs(64.0 * eps_star * (1.0 - eps_star) - 1.0) < 1e-12
    assert abs(math.sqrt(15.0) / 8.0 - 0.48412291827592707) < 1e-15

    eps = (1.0 - (math.sqrt(15.0) / 4.0 + 0.01)) / 2.0
    base = 64.0 * eps * (1.0 - eps)
    assert 0.0 < base < 1.0

    def bound(n: int) -> float:
        return math.sqrt(n) / 4.0 * base ** (n / 4.0)

    assert bound(8000) < bound(4000) < 1e-6


def test_beta_small_ball_bound():
    for a in (2.0, 5.0, 12.5, 40.0):
        for eps in np.linspace(0.001, 0.499, 60):
            exact = float(betainc(a, a, eps))
            assert exact <= math.sqrt(a) * (4.0 * eps * (1.0 - eps)) ** a + 1e-15


def test_trace_identity_behind_the_product_inequality():
    rng = np.random.default_rng(20261111)
    for n in (5, 9, 12):
        U = np.triu(rng.choice([-1.0, 1.0], size=(n, n)), 1)
        A = U + U.T
        d = rng.uniform(0.5, 6.0, size=n)
        dp = rng.uniform(0.5, 6.0, size=n)
        delta = rng.normal(size=n)
        Delta = np.diag(delta)
        lhs = np.trace((np.diag(d) - A) @ Delta @ (np.diag(dp) + A) @ Delta)
        rhs = float(np.sum((1.0 + d * dp) * delta**2) - np.sum(delta) ** 2)
        assert abs(lhs - rhs) < 1e-9 * (1.0 + abs(rhs))


def test_scalar_minimum_under_the_harmonic_constraint():
    rng = np.random.default_rng(7)
    for n in range(2, 13):
        floor = n * math.sqrt(n - 1)
        for concentration in (0.05, 0.3, 1.0, 5.0):
            for _ in range(200):
                t = rng.dirichlet(np.full(n, concentration))
                t = np.clip(t, 1e-12, 1.0 - 1e-12)
                t *= rng.uniform(0.5, 1.0) / t.sum()
                assert np.sum(np.sqrt(1.0 / t - 1.0)) >= floor - 1e-9
        flat = np.full(n, 1.0 / n)
        assert abs(np.sum(np.sqrt(1.0 / flat - 1.0)) - floor) < 1e-9


def test_product_bound_is_attained_by_explicit_dual_certificates():
    for q in (5, 13):
        C = _paley(q)
        n = q + 1
        assert np.allclose(C @ C, (n - 1) * np.eye(n))
        root = math.sqrt(n - 1)
        assert np.linalg.eigvalsh(root * np.eye(n) - C).min() > -1e-9
        assert np.linalg.eigvalsh(root * np.eye(n) + C).min() > -1e-9
        # weak duality: V+ <= n*root/2 and V- <= n*root/2
        assert abs((n * root / 2.0) ** 2 - n * n * (n - 1) / 4.0) < 1e-9
        assert abs(n * root / 2.0 - spherical_half_bound(n)) < 1e-9
    for n in (7, 8):
        A = -(np.ones((n, n)) - np.eye(n))
        assert np.linalg.eigvalsh(np.eye(n) - A).min() > -1e-9
        assert np.linalg.eigvalsh((n - 1) * np.eye(n) + A).min() > -1e-9
        assert abs((n / 2.0) * (n * (n - 1) / 2.0) - n * n * (n - 1) / 4.0) < 1e-9


def test_boolean_product_falls_below_the_vector_hyperbola_at_order_ten():
    table = json.loads((ROOT / "evidence" / "exact_m_table.json").read_text())
    m10 = table["10"]["m"]
    assert m10 == 13
    assert m10 * m10 < 10 * 10 * 9 / 4


def test_small_paley_values_and_flat_lifts():
    C6, C14, C18 = _paley(5), _paley(13), _paley(17)
    assert phi(C6) == 5.0
    assert phi(C14) == 21.0
    assert _phi_fast(C14) == 21.0
    assert _phi_fast(C18) == 33.0
    for d in (-1, 1):
        K = _lift(C6, d)
        assert np.allclose(K, K.T) and np.all(np.diag(K) == 0)
        assert phi(K) == 18.0
    H = np.block([[C14 + np.eye(14), C14 - np.eye(14)],
                  [C14 - np.eye(14), -C14 - np.eye(14)]])
    assert np.allclose(H @ H, 28 * np.eye(28))
    assert abs(np.trace(H)) < 1e-12
    for d, expected in ((1, 20.0), (2, 20.0), (3, 19.0), (4, 17.0)):
        assert _phi_fast(C14[d:, d:]) == expected
    assert _phi_fast(C18[1:, 1:]) == 32.0
    assert _phi_fast(C18[2:, 2:]) == 32.0


def test_comparison_table_receipt():
    rows = json.loads((RECEIPTS / "comparison_table.json").read_text())
    exact = json.loads((ROOT / "evidence" / "exact_m_table.json").read_text())
    assert [row["n"] for row in rows] == [6, 8, 10, 12, 14, 16, 18, 20, 24, 26, 28, 30, 38]
    for row in rows:
        n = row["n"]
        assert abs(row["S"] - _spectral_ceiling(n)) < 1e-3
        assert abs(row["haar_times_S"] - row["haar_mean"] * _spectral_ceiling(n)) < 0.02
        assert row["best_seidel"] % 2 == (n * (n - 1) // 2) % 2
        assert row["largest_admissible_below"] % 2 == (n * (n - 1) // 2) % 2
        assert row["largest_admissible_below"] <= row["haar_times_S"]
        assert row["steps_below"] in (0, 1)
        assert -0.05 <= row["gap"] <= 2.7
        assert row["best_seidel"] < row["S"]
        if str(n) in exact:
            assert row["best_seidel"] == exact[str(n)]["m"]
    assert sum(row["steps_below"] == 0 for row in rows) == 9


def test_exact_evaluation_receipt():
    data = json.loads((RECEIPTS / "exact_evaluations.json").read_text())
    lift36 = data["lifts"]["order36 K(C18,C18-I)"]
    assert lift36["maxQ"] == 108 and lift36["minQ"] == -108
    assert lift36["Phi"] == 0.5 * 36**1.5
    assert data["restrictions"]["C38 minus 2"] == 108
    for key in ("C14 minus 1,2,3,4", "C18 minus 1,2,3,4", "C30 minus 1..6"):
        values = data["restrictions"][key]
        assert values == sorted(values, reverse=True)
    haar = [json.loads(line) for line in (RECEIPTS / "rom_results.jsonl").read_text().splitlines()]
    cap = math.sqrt(15.0) / 4.0
    for row in haar:
        assert 0.0 < row["one_sided_mean"] <= row["two_sided_mean"] < 1.0
        if row["n"] >= 10:
            assert row["one_sided_mean"] < cap
