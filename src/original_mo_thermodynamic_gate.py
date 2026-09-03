"""Exact arithmetic for the critical-temperature thermodynamic audit.

This module does not claim that the MathOverflow limit is closed.  It records
the finite-temperature entropy sandwich, the common-raw-temperature block
gap, and the exact second-moment expression used in the graphon no-go.
"""

from __future__ import annotations

from math import comb, cosh, log, sqrt, tanh


def softmax_entropy_error(c: float) -> float:
    """Return the uniform error in alpha_n <= pressure/c + error."""

    if c <= 0:
        raise ValueError("c must be positive")
    return log(2.0) / c


def common_beta_equal_split_gap(n: int, c: float) -> float:
    """Per-spin random-completion gap for N=2n at beta=c/sqrt(2n).

    If p_n(beta)=min_A log E exp(beta Q_A), the exact common-beta
    inequalities give

        f_n(c/sqrt(2)) <= f_{2n}(c)
          <= f_n(c/sqrt(2)) + common_beta_equal_split_gap(n,c),

    where f_n(c)=p_n(c/sqrt(n))/n.
    """

    if n < 1:
        raise ValueError("n must be positive")
    if c < 0:
        raise ValueError("c must be nonnegative")
    return 0.5 * n * log(cosh(c / sqrt(2.0 * n)))


def common_beta_equal_split_gap_limit(c: float) -> float:
    """Limit of ``common_beta_equal_split_gap(n,c)`` as n tends to infinity."""

    if c < 0:
        raise ValueError("c must be nonnegative")
    return c * c / 8.0


def random_second_moment_weight(n: int, overlap_sum: int, c: float) -> float:
    """Exact contribution for a pair with overlap sum ``overlap_sum``.

    The returned value is the pair's contribution to
    E_J[Z_J^2]/E_J[Z_J]^2 for Z_J=E_x exp(c Q_J(x)/sqrt(n)).
    """

    if n < 1:
        raise ValueError("n must be positive")
    if abs(overlap_sum) > n or (n - overlap_sum) % 2:
        raise ValueError("overlap_sum is not a sum of n signs")
    u = tanh(c / sqrt(float(n))) ** 2
    same = (n * n + overlap_sum * overlap_sum - 2 * n) // 4
    opposite = (n * n - overlap_sum * overlap_sum) // 4
    return (1.0 + u) ** same * (1.0 - u) ** opposite


def random_second_moment_ratio(n: int, c: float) -> float:
    """Exact binomial sum for E_J[Z_J^2]/E_J[Z_J]^2."""

    if n < 1:
        raise ValueError("n must be positive")
    total = 0.0
    for negatives in range(n + 1):
        overlap_sum = n - 2 * negatives
        total += comb(n, negatives) * random_second_moment_weight(
            n, overlap_sum, c
        )
    return total / (2.0**n)


def random_second_moment_uniform_bound(c: float) -> float:
    """Uniform bound valid for every n when 0 <= c < 1."""

    if not 0 <= c < 1:
        raise ValueError("the second-moment bound requires 0 <= c < 1")
    return 1.0 / sqrt(1.0 - c * c)


def graphon_pressure_gap(c: float) -> float:
    """Random high-temperature rate minus the conference upper rate."""

    if c < 0:
        raise ValueError("c must be nonnegative")
    return c * c / 4.0 - 0.5 * log(cosh(c))


def random_annealed_pressure(n: int, c: float) -> float:
    """n^{-1} log E_J E_x exp(c Q_J(x)/sqrt(n))."""

    if n < 1:
        raise ValueError("n must be positive")
    edges = n * (n - 1) // 2
    return edges * log(cosh(c / sqrt(float(n)))) / n
