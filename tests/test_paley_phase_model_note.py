"""Checks for evidence/NOTE_2026-11-13_PALEY_PHASE_MODEL.md."""
import json
import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts", "phase_model"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "parity_ceiling"))

from phase_lib import (admissible, exact_Q, rho_star, rho_star_all_t, split_model,  # noqa: E402
                       square_state, to_elements)
from paley_lib import paley, two_squares  # noqa: E402
from ceiling import levels_by_class  # noqa: E402

EVID = os.path.join(ROOT, "evidence", "phase_model_20261113")


@pytest.mark.parametrize("q", [13, 29, 61, 101, 197])
def test_lemmas_1_and_2(q):
    g, chi, pw, h, J = split_model(q)
    n1, M = q - 1, (q - 1) // 2
    S, T = np.meshgrid(np.arange(n1), np.arange(n1), indexing="ij")
    assert (chi[(pw[S] - pw[T]) % q] == ((-1) ** S) * h[(T - S) % n1]).all()
    k = np.arange(n1)
    assert (h[(-k) % n1] == ((-1) ** k) * h).all()
    off = (k != 0) & (k != M)
    assert np.allclose(np.abs(J[off]) ** 2, q)
    assert np.allclose(J[[0, M]], -1)
    assert np.allclose(J[(M - k) % n1], J)
    a, b = two_squares(q)
    assert abs(abs(J[M // 2].imag) - a) < 1e-8 and abs(abs(J[M // 2].real) - b) < 1e-8


@pytest.mark.parametrize("q", [29, 101, 181])
def test_proposition_3_against_dense_matrix(q):
    g, chi, pw, h, J = split_model(q)
    n1, M = q - 1, (q - 1) // 2
    k = np.arange(n1)
    C = np.asarray(paley(q), dtype=np.int64)
    rng = np.random.default_rng(q)
    for _ in range(25):
        z = rng.choice([-1, 1], size=n1)
        x0, xi = int(rng.choice([-1, 1])), int(rng.choice([-1, 1]))
        zh = np.fft.fft(z) / n1
        pred = (n1 * np.sum(J * zh * zh[(M - k) % n1])).real / 2 + xi * x0 + xi * int(z.sum()) \
            + x0 * int((((-1) ** np.arange(n1)) * z).sum())
        x = to_elements(q, pw, z, x0)
        full = np.concatenate(([xi], x))
        dense = int(full @ C @ full) // 2
        assert exact_Q(q, chi, xi, x) == dense
        assert abs(pred - dense) < 1e-6


@pytest.mark.parametrize("q", [61, 181, 1009, 2521])
def test_theorem_R_exact_integers(q):
    g, chi, pw, h, J = split_model(q)
    n1, M = q - 1, (q - 1) // 2
    phi = np.angle(J)
    rng = np.random.default_rng(7 * q)
    done = with_gcd = 0
    while done < 30:
        k = int(rng.integers(1, n1))
        gg = math.gcd(k, M)
        if (M // gg) % 2:
            continue
        delta = int(rng.integers(0, M // gg))
        z = square_state(q, k, delta)
        assert int(z.sum()) == 0 and int(z[0::2].sum()) == int(z[1::2].sum())
        Q = exact_Q(q, chi, 1, to_elements(q, pw, z, 1))
        assert abs(Q - (1 + 0.5 * n1 * math.sqrt(q) * rho_star(q, phi, k, delta))) < 1e-5
        done += 1
        with_gcd += gg > 1
    assert with_gcd > 0


def test_weights_sum_to_one():
    for M1 in (6, 30, 252, 5004):
        m = np.arange(1, M1, 2)
        assert abs(np.sum(4 / (M1 ** 2 * np.sin(np.pi * m / M1) ** 2)) - 1) < 1e-10


@pytest.mark.parametrize("q", [1009, 5009])
def test_corollary_R2(q):
    g, chi, pw, h, J = split_model(q)
    n1, M = q - 1, (q - 1) // 2
    phi = np.angle(J)
    checked = 0
    for k in admissible(q)[::2]:
        k = int(k)
        for r in (1, 2, 3):
            eps = np.array([(-1) ** ((m - 1) // 2) for m in range(1, 2 * r, 2)])
            s = np.array([math.sin(phi[(m * k) % n1]) for m in range(1, 2 * r, 2)])
            for sigma in (1, -1):
                eta = 1 - float((sigma * eps * s).min())
                if eta > 1:
                    continue
                tau = M // 2 if (M // 2) % 2 else M // 2 + 1
                delta = ((sigma * tau - k) // 2) % M
                assert (2 * delta + k - sigma * tau) % (2 * M) == 0
                bound = 1 - eta - 4 / (math.pi ** 2 * r) - 2 * math.pi ** 2 * r ** 2 / M ** 2
                assert rho_star(q, phi, k, delta) >= bound - 1e-9
                checked += 1
    assert checked > 100


def test_scan_records():
    recs = json.load(open(os.path.join(EVID, "scan.json")))
    assert [r["q"] for r in recs] == [101, 229, 509, 1009, 2017, 5009, 10009, 20021, 50021, 100049,
                                      200009, 500009, 1000033]
    for r in recs:
        q = r["q"]
        S = 0.5 * (q + 1) * math.sqrt(q)
        assert r["Q_square"] <= r["Q_polished"] < S
        assert abs(r["rho_square"] - r["Q_square"] / S) < 1e-12
        assert r["weyl_max_times_sqrt_q"] < 8
        for b in r["boxes"]:
            if b["haar_prediction"] > 200:
                assert 0.8 < b["count"] / b["haar_prediction"] < 1.25
        if q >= 5009:
            assert r["rho_square"] > 0.94
        if q >= 50021:
            assert r["rho_square"] > 0.955


@pytest.mark.parametrize("q", [101, 229, 509, 1009, 2017, 5009])
def test_scan_states_recomputed_exactly(q):
    rec = {r["q"]: r for r in json.load(open(os.path.join(EVID, "scan.json")))}[q]
    g, chi, pw, h, J = split_model(q)
    assert g == rec["primitive_root"]
    x = to_elements(q, pw, square_state(q, rec["k"], rec["delta"]), 1)
    assert exact_Q(q, chi, 1, x) == rec["Q_square"]


@pytest.mark.parametrize("q", [101, 229, 509])
def test_scan_is_the_square_wave_maximum(q):
    rec = {r["q"]: r for r in json.load(open(os.path.join(EVID, "scan.json")))}[q]
    g, chi, pw, h, J = split_model(q)
    phi = np.angle(J)
    best = -9.0
    for k in admissible(q):
        F = rho_star_all_t(q, phi, int(k))
        best = max(best, float(F[1::2].max()))
    assert round(1 + 0.5 * (q - 1) * math.sqrt(q) * best) == rec["Q_square"]


def test_known_exact_value_not_exceeded():
    rec = {r["q"]: r for r in json.load(open(os.path.join(EVID, "scan.json")))}[101]
    assert rec["Q_polished"] <= 501            # Phi(C102) = 501 (NOTE_2026-11-11)


def _graph6(s):
    data = [ord(c) - 63 for c in s.strip()]
    n, data = data[0], data[1:]
    assert n <= 62
    bits = [(d >> k) & 1 for d in data for k in range(5, -1, -1)]
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for j in range(1, n):
        for i in range(j):
            A[i, j] = A[j, i] = bits[idx]
            idx += 1
    return A


def test_mathon_order_46():
    top = levels_by_class(45, depth=1, use_row_identity=True)[0]
    assert top["Q"] == 147 and top["multisets"] == [{3: 10, 7: 33, 11: 3}]
    data = json.load(open(os.path.join(EVID, "mathon46.json")))
    assert len(data["graphs"]) == 6 and data["distinct_K4_profiles"] == 1
    for rec in data["graphs"]:
        A = _graph6(rec["graph6"])
        assert A.shape == (45, 45) and set(A.sum(1).tolist()) == {22}
        C = np.zeros((46, 46), dtype=np.int64)
        C[0, 1:] = 1
        C[1:, 0] = 1
        C[1:, 1:] = np.ones((45, 45), dtype=np.int64) - np.eye(45, dtype=np.int64) - 2 * A
        assert (C @ C == 45 * np.eye(46, dtype=np.int64)).all()
        x = np.array([1 if c == "+" else -1 for c in rec["witness"]], dtype=np.int64)
        assert int(x @ C @ x) // 2 == 147
        fields = sorted((x * (C @ x)).tolist())
        assert fields == [3] * 10 + [7] * 33 + [11] * 3
