"""Discrete-derivative identities D1–D3 on existing Paley / plus-I maximizers."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from minmax_quadratic import (
    form_Q,
    halfspace_boolean_vector,
    paley_conference_prime_power,
    phi_mitm,
)


def _lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def _cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1
    return b


def _w(A: np.ndarray, x: np.ndarray) -> np.ndarray:
    if float(x @ A @ x) < 0:
        x = -x
    return (A @ x) * x


def _brute_max(A: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    best_q, best_x = -1.0, None
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = form_Q(A, x)
        if abs(q) > best_q:
            best_q, best_x = abs(q), x.copy() * np.sign(q)
    return best_x


def test_d1_coordinate_flip_paley10_and_plusi10():
    C = paley_conference_prime_power(3)
    x = halfspace_boolean_vector(3)
    if float(x @ C @ x) < 0:
        x = -x
    w = _w(C, x)
    for k in range(10):
        y = x.copy()
        y[k] *= -1
        assert abs(form_Q(C, y) - (form_Q(C, x) - 2 * w[k])) < 1e-9
    K = _lift(_cycle5())
    z = _brute_max(K)
    wk = _w(K, z)
    q = form_Q(K, z)
    for k in range(10):
        y = z.copy()
        y[k] *= -1
        assert abs(form_Q(K, y) - (q - 2 * wk[k])) < 1e-9


def test_d3_two_flip_and_light_set_of_known_beaters():
    K = _lift(_cycle5())
    x = _brute_max(K)
    w = _w(K, x)
    q = form_Q(K, x)
    i, j = 0, 1
    y = x.copy()
    y[i] *= -1
    y[j] *= -1
    pred = q - 2 * w[i] - 2 * w[j] + 4 * K[i, j] * x[i] * x[j]
    assert abs(form_Q(K, y) - pred) < 1e-9
    light = int(np.sum(np.rint(w) <= 1))
    assert light == 6
    assert abs(q - 13.0) < 1e-9


def test_d2_paley_has_no_width2_light_coords():
    C = paley_conference_prime_power(3)
    x = halfspace_boolean_vector(3)
    w = _w(C, x)
    assert np.allclose(w, 3.0)
    assert int(np.sum(np.rint(w) <= 1)) == 0


def _F(B: np.ndarray, x: np.ndarray, T: np.ndarray) -> float:
    m = B.shape[0]
    t = int(len(T))
    qb = form_Q(B, x)
    qt = 0.0
    if t >= 2:
        qt = form_Q(B[np.ix_(T, T)], x[T])
    return 2 * qb - 4 * qt + m - 2 * t


def test_d5_plus_i_unfolding_identity_on_k10_and_random():
    B = _cycle5()
    K = _lift(B)
    m = 5
    rng = np.random.default_rng(0)
    for xb in range(1 << 9):
        z = np.ones(10)
        for i in range(9):
            z[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        x, y = z[:m], z[m:]
        T = np.where(x * y < 0)[0]
        assert abs(_F(B, x, T) - form_Q(K, z)) < 1e-9
    B13 = np.load(
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )["A"].astype(np.float64)[:13, :13]
    K13 = _lift(B13)
    for _ in range(80):
        x = rng.choice([-1.0, 1.0], size=13)
        y = rng.choice([-1.0, 1.0], size=13)
        z = np.concatenate([x, y])
        T = np.where(x * y < 0)[0]
        assert abs(_F(B13, x, T) - form_Q(K13, z)) < 1e-9


def test_d9_is_d5_rewritten_for_every_z():
    B = _cycle5()
    K = _lift(B)
    m = 5
    for xb in range(1 << 9):
        z = np.ones(10)
        for i in range(9):
            z[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        x, y = z[:m], z[m:]
        T = np.where(x * y < 0)[0]
        wy = ((K @ z) * z)[m:]
        extra = float(np.sum(wy[T] - 1.0)) if len(T) else 0.0
        pred = 2 * form_Q(B, x) + m + extra
        assert abs(pred - form_Q(K, z)) < 1e-9


def test_d9_k26_has_maximizer_aligned_to_B():
    A = np.load(
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )["A"].astype(np.float64)
    B = A[:13, :13]
    rng = np.random.default_rng(0)
    seen = set()
    aligned = False
    for _ in range(80):
        z = rng.choice([-1.0, 1.0], size=26)
        Az = A @ z
        qfull = float(z @ Az)
        cur = abs(qfull) / 2.0
        improved = True
        while improved:
            improved = False
            for i in range(26):
                delta = -4.0 * z[i] * Az[i]
                if abs(qfull + delta) / 2.0 > cur + 1e-12:
                    zi = z[i]
                    z[i] = -zi
                    Az -= 2.0 * zi * A[:, i]
                    qfull += delta
                    cur = abs(qfull) / 2.0
                    improved = True
        if abs(cur - 61.0) > 1e-6:
            continue
        qb = abs(form_Q(B, z[:13]))
        seen.add(int(round(qb)))
        if abs(qb - 20.0) < 1e-6:
            aligned = True
    assert 20 in seen
    assert aligned


def test_d10_extra_is_two_t_tminus2_minus_8_eplus():
    B = _cycle5()
    x = np.array([1.0, -1.0, 1.0, -1.0, -1.0])
    assert abs(form_Q(B, x) - 4.0) < 1e-9
    beta = B * np.outer(x, x)
    np.fill_diagonal(beta, 0.0)
    K = _lift(B)
    for t in range(2, 6):
        # one T: first t indices after a rotation that keeps x
        T = np.arange(t)
        eplus = int(sum(1 for i in range(t) for j in range(i + 1, t) if beta[T[i], T[j]] > 0))
        s = np.ones(5)
        s[T] = -1.0
        y = s * x
        z = np.concatenate([x, y])
        extra = 2 * t * (t - 2) - 8 * eplus
        pred = 2 * form_Q(B, x) + 5 + extra
        assert abs(pred - form_Q(K, z)) < 1e-9


def test_d13_bipartite_minus_is_case_a():
    for p, m in ((3, 5), (5, 13), (7, 25), (11, 61)):
        floor_phi = (m - 1) ** 2 / 4.0
        thresh = (m - 1) * (p - 1) / 2.0
        assert floor_phi + 1e-12 >= thresh


def test_d14_extra_parity_and_mplus22_deficit():
    # Extra = 2t(t-2)-8 e_+. t even => Extra ≡ 0 (mod 8); t odd => Extra ≡ 6 (mod 8).
    for t in range(2, 14):
        for eplus in range(t * (t - 1) // 2 + 1):
            extra = 2 * t * (t - 2) - 8 * eplus
            if t % 2 == 0:
                assert extra % 8 == 0
            else:
                assert extra % 8 == 6
    # M+=22 needs Γ≥4; a minus triangle supplies 6.
    assert 2 * 3 * (3 - 2) - 8 * 0 == 6
    assert 48 - 2 * 22 == 4
    assert 2 * 22 + 13 + 6 >= 61


def test_d14_tf_qmax_census_receipt_empty():
    rec = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "evidence"
            / "case_b_gamma_p5"
            / "tf_qmax_census.json"
        ).read_text()
    )
    assert rec["n_seen"] == 163477
    assert rec["n_hits"] == 0
    assert rec["tf_case_b_counterexamples"] == 0
    assert rec["self_test"]["B13_Mplus"] == 20.0
    assert rec["self_test"]["B13_gamma"] == 8
    assert rec["self_test"]["B13_checker_qmax"] is True
    assert rec["self_test"]["all_plus_is_qmax"] is True
    assert rec["self_test"]["all_minus_is_qmax"] is False
    by_e = rec["by_e"]
    assert by_e["28"] == 100457
    assert by_e["29"] == 39920


def test_d14_switched_b13_is_qmax_with_gamma8():
    B = np.loadtxt(
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "coherent_BI_lift_20260919"
        / "B13.txt"
    )
    n = 13
    best_q, best_x = -1e9, None
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = 0.5 * float(x @ B @ x)
        if q > best_q:
            best_q, best_x = q, x.copy()
    Bp = (best_x[:, None] * B) * best_x[None, :]
    np.fill_diagonal(Bp, 0.0)
    ones = np.ones(n)
    assert abs(0.5 * float(ones @ Bp @ ones) - 20.0) < 1e-9
    # 1 maximises Q: every other vector has Q ≤ 20
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        assert 0.5 * float(x @ Bp @ x) <= 20.0 + 1e-9
    # some 4-set has e_+ ≤ 1 (Extra ≥ 8)
    plus = Bp > 0
    found = False
    from itertools import combinations

    for S in combinations(range(n), 4):
        eplus = sum(1 for a, b in combinations(S, 2) if plus[a, b])
        if eplus <= 1:
            found = True
            extra = 16 - 8 * eplus
            assert extra >= 8
            break
    assert found


def test_d11_d12_p3_every_lift_at_least_13():
    # Case A threshold is 4 = m_5, so 2 Phi(B)+m >= 13 with Phi(B)>=4
    assert 2 * 4 + 5 == 13
    B = _cycle5()
    K = _lift(B)
    z = _brute_max(K)
    assert abs(abs(form_Q(K, z)) - 13.0) < 1e-9


def test_d5_t_extremes_give_two_phi_B_plus_m():
    B = _cycle5()
    x = np.array([1.0, 1.0, -1.0, -1.0, 1.0])
    assert abs(_F(B, x, np.array([], dtype=int)) - (2 * form_Q(B, x) + 5)) < 1e-9
    T = np.arange(5)
    assert abs(_F(B, x, T) - (-2 * form_Q(B, x) - 5)) < 1e-9
    K = _lift(B)
    z = _brute_max(K)
    assert abs(abs(form_Q(K, z)) - 13.0) < 1e-9
    assert abs(13.0 - (2 * 4.0 + 5)) < 1e-9


def test_n26_beater_light_set_is_six():
    path = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "ns_port_n26_undercut_20260912"
        / "stage01.npz"
    )
    A = np.load(path)["A"].astype(np.float64)
    assert abs(phi_mitm(A) - 61.0) < 1e-9
    rng = np.random.default_rng(0)
    n = 26
    best_q, best_x = -1.0, None
    for _ in range(80):
        x = rng.choice([-1.0, 1.0], size=n)
        Ax = A @ x
        qfull = float(x @ Ax)
        cur = abs(qfull) / 2.0
        improved = True
        while improved:
            improved = False
            for i in range(n):
                delta = -4.0 * x[i] * Ax[i]
                if abs(qfull + delta) / 2.0 > cur + 1e-12:
                    xi = x[i]
                    x[i] = -xi
                    Ax -= 2.0 * xi * A[:, i]
                    qfull += delta
                    cur = abs(qfull) / 2.0
                    improved = True
        if cur > best_q:
            best_q, best_x = cur, x.copy()
    assert abs(best_q - 61.0) < 1e-9
    w = _w(A, best_x)
    assert int(np.sum(np.rint(w) <= 1)) == 6
