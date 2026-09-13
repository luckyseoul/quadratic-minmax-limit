"""Exact polynomial corroboration of the gauge temperature normalization.

One fixed 3+3 signing is used. Laurent coefficients are integers: no
floating exponential, logarithm, cumulant truncation, or search over source
signings enters this check. The all-orders proof is Section 4.1 of
NOTE_2026-09-02_UPSTREAM_RELATIVE_GAUGE_BRIDGE.md.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import socket
import time


LEFT = ((0, 1, -1), (1, 0, 1), (-1, 1, 0))
RIGHT = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
CROSS = ((1, -1, 1), (-1, -1, 1), (1, 1, -1))


def projective(order: int):
    for tail in product((-1, 1), repeat=order - 1):
        yield (1,) + tail


def quadratic(matrix, spin):
    return sum(
        matrix[i][j] * spin[i] * spin[j]
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
    )


def rectangular(matrix, left, right):
    return sum(
        matrix[i][j] * left[i] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def even_polynomial(order, energy):
    # Z(log u) = sum_e coefficients[e] u^e / 2^order.
    coefficients = Counter()
    maximum = 0
    for spin in projective(order):
        value = energy(spin)
        maximum = max(maximum, abs(value))
        coefficients[value] += 1
        coefficients[-value] += 1
    assert sum(coefficients.values()) == 2 ** order
    assert max(coefficients) == maximum
    assert coefficients[maximum] >= 1
    assert all(coefficients[e] == coefficients[-e] for e in coefficients)
    return coefficients, maximum


def multiply(left, right):
    result = Counter()
    for a, count_a in left.items():
        for b, count_b in right.items():
            result[a + b] += count_a * count_b
    return result


def scaled(coefficients, factor):
    return Counter({e: factor * count for e, count in coefficients.items()})


def serializable(coefficients):
    return [[e, count] for e, count in sorted(coefficients.items())]


def main():
    started = time.monotonic()
    n, k = len(LEFT), len(RIGHT)
    order = n + k
    pa, _ = even_polynomial(n, lambda x: quadratic(LEFT, x))
    pb, _ = even_polynomial(k, lambda y: quadratic(RIGHT, y))
    pc, _ = even_polynomial(
        order, lambda z: rectangular(CROSS, z[:n], z[n:])
    )
    product_numerator = multiply(multiply(pa, pb), pc)
    gauge_numerator = Counter()
    positive_tau_numerator = Counter()
    checks = 0
    for alpha in projective(n):
        for beta_g in projective(k):
            for tau in (-1, 1):
                def energy(spin):
                    x, y = spin[:n], spin[n:]
                    ax = tuple(a * v for a, v in zip(alpha, x))
                    by = tuple(b * v for b, v in zip(beta_g, y))
                    return (
                        quadratic(LEFT, ax) + tau * quadratic(RIGHT, by)
                        + rectangular(CROSS, x, y)
                    )

                coefficients, maximum = even_polynomial(order, energy)
                gauge_numerator.update(coefficients)
                if tau == 1:
                    positive_tau_numerator.update(coefficients)
                # A distinct exact evaluation checks the norm sandwich at
                # beta=log(2), including its projective factor 2^order.
                z_at_two = sum(
                    (Fraction(2) ** e) * count
                    for e, count in coefficients.items()
                ) / (2 ** order)
                assert Fraction(2) ** maximum / (2 ** order) <= z_at_two
                assert z_at_two <= Fraction(2) ** maximum
                checks += 1

    assert checks == 2 ** (order - 1)
    # Average Z_g has denominator 2^(2*order-1); the product of the
    # three normalized source partitions has denominator 2^(2*order).
    assert scaled(gauge_numerator, 2) == product_numerator
    # Dropping tau leaves half as many gauges; this fixture then fails
    # the proposed normalization, so the check detects that omission.
    omitted_tau_is_different = (
        scaled(positive_tau_numerator, 4) != product_numerator
    )
    assert omitted_tau_is_different
    fixture = {"left": LEFT, "right": RIGHT, "cross": CROSS}
    fixture_bytes = json.dumps(
        fixture, sort_keys=True, separators=(",", ":")
    ).encode()
    return {
        "classification": "finite_exact_partition_identity_only",
        "host": socket.gethostname(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fixture": fixture,
        "fixture_sha256": hashlib.sha256(fixture_bytes).hexdigest(),
        "source_signings_checked": 1,
        "gauges_checked": checks,
        "arithmetic": "integer Laurent coefficients and exact rational evaluation",
        "full_gauge_average_equals_product": True,
        "omitting_tau_changes_average": omitted_tau_is_different,
        "norm_sandwich_checks_at_exp_beta_two": checks,
        "product_numerator": serializable(product_numerator),
        "temperature_model_inequalities_11_checked": False,
        "original_convergence_proved": False,
        "passed": True,
        "elapsed_seconds": time.monotonic() - started,
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
