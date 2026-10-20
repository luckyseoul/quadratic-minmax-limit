"""Integer sequences that obey the proved scalar axioms and do not converge.

The axioms are: m_n is a positive integer, m is nondecreasing, and
m_{n+1} - m_n <= n - 2 for every large n.  CORE section 4 gives
liminf alpha_n >= 1/pi and the conference bound gives limsup alpha_n <= 1/2.
Those constraints do not force the two ends to agree.  For any
0 < lambda <= Lambda <= 1/2 the sequence below has those exact ends.

The sequence is not a sequence of Seidel minima.  It is an obstruction to
reading convergence off the scalar axioms alone.
"""
from __future__ import annotations

import math

import numpy as np

def frequency(lambda_: float, lambda_upper: float) -> float:
    """Keep (3/2) lambda + g' positive.  Slow the wave when the gap is wide."""
    b = 0.5 * (lambda_upper - lambda_)
    if b == 0.0:
        return 1.0
    return min(1.0, 0.75 * lambda_ / b)


def ends(lambda_: float, lambda_upper: float) -> tuple[float, float, float]:
    """Return (a, b, lower derivative margin) for liminf lambda_ and limsup lambda_upper."""
    if not (0.0 < lambda_ <= lambda_upper <= 0.5):
        raise ValueError("need 0 < lambda <= Lambda <= 1/2")
    a = 0.5 * (lambda_upper + lambda_)
    b = 0.5 * (lambda_upper - lambda_)
    margin = 1.5 * lambda_ - b * frequency(lambda_, lambda_upper)
    if margin <= 0.0:
        raise ValueError("derivative margin is not positive")
    return a, b, margin


def high_slope(lambda_: float, lambda_upper: float) -> float:
    _a, b, _margin = ends(lambda_, lambda_upper)
    return 1.5 * lambda_upper + b * frequency(lambda_, lambda_upper)


def first_n(lambda_: float, lambda_upper: float) -> int:
    """First n at which the continuous increment is at least 1 and the step ceiling holds."""
    margin = ends(lambda_, lambda_upper)[2]
    ceiling = high_slope(lambda_, lambda_upper)
    n = 2
    while True:
        if margin * math.sqrt(n) >= 1.0 and ceiling * math.sqrt(n + 1) + 1.0 <= n - 2:
            return n
        n += 1
        if n > 10**7:
            raise RuntimeError("no starting index")


def values(lambda_: float, lambda_upper: float, n0: int, n1: int) -> np.ndarray:
    """m_n = floor(h(n) + 1/2) on n0 <= n < n1, h(n) = n^{3/2} g(ln n)."""
    a, b, _margin = ends(lambda_, lambda_upper)
    omega = frequency(lambda_, lambda_upper)
    n = np.arange(n0, n1, dtype=np.float64)
    h = np.power(n, 1.5) * (a + b * np.sin(omega * np.log(n)))
    return np.floor(h + 0.5)


def axioms_hold(lambda_: float, lambda_upper: float, n0: int, n1: int) -> dict:
    m = values(lambda_, lambda_upper, n0, n1)
    n = np.arange(n0, n1, dtype=np.int64)
    diffs = np.diff(m)
    steps_ok = bool(np.all(diffs <= (n[:-1] - 2)))
    mono_ok = bool(np.all(diffs >= 1))
    alpha = m / np.power(n.astype(np.float64), 1.5)
    return {
        "monotone": mono_ok,
        "step": steps_ok,
        "alpha_min": float(alpha.min()),
        "alpha_max": float(alpha.max()),
        "n0": n0,
        "n1": n1,
    }


def full_gap_pair() -> tuple[float, float]:
    """The CORE Gaussian floor and the conference ceiling."""
    return 1.0 / math.pi, 0.5
