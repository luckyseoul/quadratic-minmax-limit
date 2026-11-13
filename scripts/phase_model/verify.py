"""Exact checks of the phase model identities.  Writes JSON to stdout."""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase_lib import admissible, exact_Q, rho_star, split_model, square_state, to_elements  # noqa: E402
from paley_lib import paley, two_squares  # noqa: E402


def check(q, rng, trials=40):
    g, chi, pw, h, J = split_model(q)
    n1, M = q - 1, (q - 1) // 2
    k = np.arange(n1)
    rec = {"q": q}
    S, T = np.meshgrid(np.arange(n1), np.arange(n1), indexing="ij") if q <= 2600 else (None, None)
    if S is not None:
        rec["twisted_circulant"] = bool((chi[(pw[S] - pw[T]) % q] == ((-1) ** S) * h[(T - S) % n1]).all())
    rec["h_parity_symmetry"] = bool((h[(-k) % n1] == ((-1) ** k) * h).all())
    off = (k != 0) & (k != M)
    rec["abs_J_squared_is_q"] = bool(np.allclose(np.abs(J[off]) ** 2, q))
    rec["J_0_and_J_M_are_minus_1"] = bool(np.allclose(J[[0, M]], -1))
    rec["J_M_minus_k_equals_J_k"] = bool(np.allclose(J[(M - k) % n1], J))
    a, b = two_squares(q)
    rec["quartic_frequency"] = {"a": a, "b": b, "abs_Im_J_is_a": bool(abs(abs(J[M // 2].imag) - a) < 1e-6),
                                "abs_Re_J_is_b": bool(abs(abs(J[M // 2].real) - b) < 1e-6)}
    # Fourier identity against exact integers
    worst = 0.0
    dense = np.asarray(paley(q), dtype=np.int64) if q <= 2600 else None
    for _ in range(10):
        z = rng.choice([-1, 1], size=n1)
        x0, xi = int(rng.choice([-1, 1])), int(rng.choice([-1, 1]))
        zh = np.fft.fft(z) / n1
        two_star = (n1 * np.sum(J * zh * zh[(M - k) % n1])).real
        pred = two_star / 2 + xi * x0 + xi * int(z.sum()) + x0 * int((((-1) ** np.arange(n1)) * z).sum())
        x = to_elements(q, pw, z, x0)
        Q = exact_Q(q, chi, xi, x)
        if dense is not None:
            full = np.concatenate(([xi], x))
            assert Q == int(full @ dense @ full) // 2
        worst = max(worst, abs(pred - Q))
    rec["fourier_identity_max_error"] = worst
    # closed formula on random admissible and non-admissible frequencies
    phi = np.angle(J)
    worst, count, with_gcd = 0.0, 0, 0
    while count < trials:
        kk = int(rng.integers(1, n1))
        gg = math.gcd(kk, M)
        if (M // gg) % 2:
            continue
        delta = int(rng.integers(0, M // gg))
        z = square_state(q, kk, delta)
        assert int(z.sum()) == 0 and int(z[0::2].sum() - z[1::2].sum()) == 0
        Q = exact_Q(q, chi, 1, to_elements(q, pw, z, 1))
        worst = max(worst, abs(Q - (1 + 0.5 * n1 * math.sqrt(q) * rho_star(q, phi, kk, delta))))
        count += 1
        with_gcd += gg > 1
    rec["square_formula"] = {"trials": count, "with_gcd_above_1": int(with_gcd), "max_abs_error": worst}
    return rec


if __name__ == "__main__":
    rng = np.random.default_rng(20261113)
    qs = [int(a) for a in sys.argv[1:]] or [13, 29, 61, 101, 181, 197, 1009, 2521, 10009]
    print(json.dumps([check(q, rng) for q in qs], indent=1))
