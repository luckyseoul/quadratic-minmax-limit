"""Checks for evidence/NOTE_2026-11-12_QUARTIC_RESONANCE.md.

Theorems Q and L are proved in the note; these tests exercise the identities
in exact integer arithmetic and re-evaluate the stored witnesses. Nothing
here bounds any m_n from below or touches the open limit.
"""
from __future__ import annotations

import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts" / "parity_ceiling"
RECEIPTS = ROOT / "evidence" / "quartic_resonance_20261112"

sys.path.insert(0, str(SCRIPTS))
import ceiling as pc  # noqa: E402
import paley_lib as pl  # noqa: E402
import quartic as qr  # noqa: E402


def _primes(lo: int, hi: int):
    return [q for q in range(lo, hi, 4) if pl.is_prime(q)]


def test_jacobi_sum_lemma():
    rng = np.random.default_rng(7)
    for m in (13, 29, 61, 101):
        g = pl.prim_root(m)
        chi = pl.chi_table(m)
        psi = np.zeros(m, dtype=complex)
        v = 1
        for k in range(m - 1):
            psi[v] = 1j ** (k % 4)
            v = v * g % m
        J = sum(chi[(1 - z) % m] * psi[z] for z in range(2, m))
        a, b = pl.two_squares(m)
        assert abs(abs(J) ** 2 - m) < 1e-9
        assert round(abs(J.real)) == b and round(abs(J.imag)) == a
        Q = qr.core(m)
        for _ in range(5):
            alpha = complex(rng.normal(), rng.normal())
            g_alpha = (alpha * psi).real
            beta = (alpha * J).conjugate()
            assert np.allclose(Q @ g_alpha, (beta * psi).real, atol=1e-8)
            assert abs(g_alpha.sum()) < 1e-8 and abs((g_alpha * chi).sum()) < 1e-8


def test_theorem_q_square_wave_on_paley():
    for q in _primes(5, 300):
        a, b = pl.two_squares(q)
        val, x = qr.quartic_conference_state(q)
        assert val == (q - 1) * a // 2 + 1
        h = sorted(set(pl.fields(pl.paley(q), x).tolist()))
        assert set(h) <= {1, a + b + 2, a + b - 2, a - b} or set(h) <= {1, a - b + 2, a - b - 2, a + b}


def test_theorem_l_square_wave_pair_on_both_lifts():
    for m in _primes(5, 200):
        a, b = pl.two_squares(m)
        for d in (1, -1):
            val, x = qr.quartic_lift_state(m, d)
            assert val == (m - 1) * (a + b) + 1
            K = qr.lift(qr.core(m), d)
            assert int(x @ K @ x) == 2 * val
    for p, m in ((3, 5), (5, 13), (9, 41), (11, 61), (15, 113)):
        assert 2 * m - 1 == p * p and pl.is_prime(m)
        assert sum(pl.two_squares(m)) == p
        assert qr.quartic_lift_state(m, 1)[0] == m * p - (p - 1)


def test_lift_fields_lie_in_one_class_mod_4():
    rng = np.random.default_rng(11)
    for m in (13, 29):
        for d in (1, -1):
            K = qr.lift(qr.core(m), d)
            for _ in range(300):
                x = rng.choice([-1, 1], size=2 * m)
                assert len(set(((x * (K @ x)) % 4).tolist())) == 1


def test_gf25_lift_state_equals_the_exact_order_50_value():
    P, NR = 5, 2  # GF(25) = GF(5)[w]/(w^2 - 2)
    els = [(c0, c1) for c0 in range(P) for c1 in range(P)]

    def mul(x, y):
        return ((x[0] * y[0] + NR * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)

    def order(e):
        k, h = 1, e
        while h != (1, 0):
            h = mul(h, e)
            k += 1
        return k

    gen = next(e for e in els if e != (0, 0) and order(e) == 24)
    ind, h = {}, (1, 0)
    for k in range(24):
        ind[h] = k
        h = mul(h, gen)
    chi = {e: (0 if e == (0, 0) else (1 if ind[e] % 2 == 0 else -1)) for e in els}
    Q = np.array([[chi[((x[0] - y[0]) % P, (x[1] - y[1]) % P)] for y in els] for x in els], dtype=np.int64)
    assert (Q @ Q == 25 * np.eye(25, dtype=np.int64) - 1).all()

    def wave(rot, s0):
        return np.array([s0 if e == (0, 0) else (1 if ((ind[e] - rot) % 4) <= 1 else -1) for e in els], dtype=np.int64)

    for d in (1, -1):
        K = qr.lift(Q, d)
        best = max(
            abs(int(x @ K @ x)) // 2
            for ru in range(4)
            for rv in range(4)
            for v0 in (1, -1)
            for x in [np.concatenate([wave(ru, 1), wave(rv, v0)])]
        )
        assert best == 24 * 7 + 1 == 169


def test_order_198_exact_value():
    q = 197
    val, x = qr.quartic_conference_state(q)
    assert val == 1373 == int(x @ pl.paley(q) @ x) // 2
    assert pc.ceiling(q, use_row_identity=False) == 1373
    # the hand argument of Section 4
    assert max(t1 for t1 in range(-198, 199, 1) if t1 % 2 and (693 - 13 * t1) // 2 >= abs(t1)
               and ((693 - 13 * t1) // 2 - t1) % 2 == 0) == 43
    assert 1287 + 2 * 43 == 1373


def test_windows_receipt():
    rows = json.loads((RECEIPTS / "windows.json").read_text())
    assert [r["q"] for r in rows] == _primes(5, 1000) and len(rows) == 80
    for r in rows:
        q = r["q"]
        assert r["L"] <= r["U"] <= r["U_integral"] <= r["closed_form"] <= r["spectral"] + 1e-9
        a, b = pl.two_squares(q)
        assert r["L"] >= max((q + 1) * b // 2, (q - 1) * a // 2 + 1)
        if r["witness_x"] and q < 400:  # the larger witnesses are re-evaluated by the generating script
            x = np.array([1 if ch == "+" else -1 for ch in r["witness_x"]], dtype=np.int64)
            assert int(x @ pl.paley(q) @ x) // 2 == r["L"]
    assert [r["q"] for r in rows if r["exact"]] == [5, 13, 17, 29, 37, 41, 53, 101, 197]
    assert sum(r["L"] / r["spectral"] >= 0.97 for r in rows) == 28
    by_q = {r["q"]: r for r in rows}
    assert (by_q[173]["L"], by_q[173]["U"]) == (1131, 1133)
    assert (by_q[677]["L"], by_q[677]["U"]) == (8789, 8797)


def test_quartic_receipt():
    d = json.loads((RECEIPTS / "quartic.json").read_text())
    assert len(d["conference"]) == 80 and len(d["lift"]) == 59
    assert [c["p"] for c in d["crystalline"]] == [3, 5, 9, 11, 15, 19, 25, 29, 35, 39, 45, 49, 51, 59]
    for c in d["crystalline"]:
        assert c["lift_state"] == c["flat"] - (c["p"] - 1) and c["flat"] == c["m"] * c["p"]
    s = json.loads((RECEIPTS / "lift_search.json").read_text())
    assert s["m61_order122"]["square_wave_state"] == 661 > s["m61_order122"]["full_cube_annealing"]["best"]
    assert math.isclose(1373 / 198 ** 1.5, 0.4928, abs_tol=5e-5)
    assert cmath.isclose(cmath.exp(1j * math.pi / 4) ** 2, 1j, abs_tol=1e-12)
