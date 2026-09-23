"""Restriction modulus for Paley gaps: gamma_q cannot rise faster than M(n, N)."""
import numpy as np


def majorant(n: int, N: int) -> float:
    sn = np.sqrt(1.0 - 1.0 / n)
    sN = np.sqrt(1.0 - 1.0 / N)
    return 0.5 * (sN - sn * (n / N) ** 1.5)


def gamma(order: int, m: float) -> float:
    return 0.5 * np.sqrt(1.0 - 1.0 / order) - m / order ** 1.5


def _odd_primes(limit: int) -> list[int]:
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    sieve[4::2] = False
    for p in range(3, int(limit**0.5) + 1, 2):
        if sieve[p]:
            sieve[p * p :: 2 * p] = False
    return [int(p) for p in range(3, limit + 1, 2) if sieve[p]]


def test_majorant_dominates_every_admissible_pair():
    rng = np.random.default_rng(0)
    pairs = [(10, 26), (10, 50), (26, 50), (50, 122), (100, 110), (100, 1000), (1000, 1010)]
    for n, N in pairs:
        cap_n = 0.5 * n * np.sqrt(n - 1)
        cap_N = 0.5 * N * np.sqrt(N - 1)
        allowance = majorant(n, N)
        for _ in range(30):
            m_n = float(rng.uniform(0.0, cap_n))
            m_N = float(rng.uniform(m_n, cap_N))
            rise = gamma(N, m_N) - gamma(n, m_n)
            assert rise <= allowance + 1e-9


def test_consecutive_paley_rises_tend_to_zero():
    primes = _odd_primes(400)
    rises = [
        majorant(p * p + 1, q * q + 1) for p, q in zip(primes, primes[1:])
    ]
    assert rises[0] > 0.2
    assert rises[-1] < 0.05
    assert max(rises[-15:]) < max(rises[:15])


def test_fixed_seed_does_not_control_the_whole_tail():
    # Restriction's majorant tends to 1/2, which is no information about gamma_q.
    far = [majorant(10, q * q + 1) for q in (50, 100, 400)]
    assert far[0] < far[1] < far[2]
    assert far[-1] > 0.45
