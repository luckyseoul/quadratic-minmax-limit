"""Exact corroboration for the uniform-gauge tripling obstruction.

One fixed 4+4 source signing tests the randomized correction, including
conditional cancellation of its quadratic self-interaction. Rational
inequalities check the strict asymptotic constant gap. No source-signing
census or numerical approximation to the all-orders theorem is performed.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import socket
import time


FREE = ((0, 1, -1, 1), (1, 0, 1, -1),
        (-1, 1, 0, 1), (1, -1, 1, 0))
PINNED = ((0, -1, 1, 1), (-1, 0, 1, 1),
          (1, 1, 0, 1), (1, 1, 1, 0))
CROSS = ((1, 1, -1, 1), (1, -1, 1, -1),
         (-1, 1, 1, 1), (1, 1, 1, 1))
SEED_SET = (0, 1)


def qform(matrix, spin):
    return sum(matrix[i][j] * spin[i] * spin[j]
               for i in range(len(matrix))
               for j in range(i + 1, len(matrix)))


def a(order):
    if order == 0:
        return Fraction(0)
    return Fraction(order * comb(order - 1, (order - 1) // 2),
                    2 ** (order - 1))


def scalar_checks():
    # x^4(1-x)^4 = (1+x^2)q(x)-4, and integral_0^1 q=22/7.
    # Positivity of the integrand then gives pi<22/7 analytically.
    quotient = [4, 0, -4, 0, 5, -4, 1]
    reconstructed = [0] * 9
    for degree, coefficient in enumerate(quotient):
        reconstructed[degree] += coefficient
        reconstructed[degree + 2] += coefficient
    reconstructed[0] -= 4
    assert reconstructed == [0, 0, 0, 0, 1, -4, 6, -4, 1]
    assert sum(Fraction(v, d + 1) for d, v in enumerate(quotient)) == Fraction(22, 7)

    root_pi_upper = Fraction(17729, 10000)
    root_five_lower = Fraction(2236, 1000)
    root_three_upper = Fraction(17321, 10000)
    root_two_lower = Fraction(14142, 10000)
    assert root_pi_upper ** 2 > Fraction(22, 7)
    assert root_five_lower ** 2 < 5
    assert root_three_upper ** 2 > 3
    assert root_two_lower ** 2 < 2
    floor_lower = (2 + root_five_lower) / (2 * root_pi_upper)
    target_upper = (3 * root_three_upper - 2 * root_two_lower) / 2
    gap_lower = floor_lower - target_upper
    assert gap_lower > Fraction(1, 100)
    return {
        "pi_upper_integral_polynomial_identity": True,
        "floor_constant_lower": str(floor_lower),
        "target_difference_constant_upper": str(target_upper),
        "gap_lower": str(gap_lower),
        "gap_exceeds_one_over_100": True,
    }


def strategy_check():
    k, m = len(FREE), len(PINNED)
    seed = SEED_SET
    remaining = tuple(i for i in range(k) if i not in seed)
    pinned_spin = (1,) * m
    pinned_energy = qform(PINNED, pinned_spin)
    pinned_norm = max(
        abs(qform(PINNED, (1,) + tail))
        for tail in product((-1, 1), repeat=m - 1)
    )
    assert pinned_energy == pinned_norm == 4
    sigma = 1
    expected_total = Fraction(0)
    contexts = 0
    outcome_checks = 0
    nonzero_fixed_tau_contexts = 0
    for row_signs in product((-1, 1), repeat=m):
        h = tuple(sum(CROSS[i][j] * row_signs[j] for j in range(m))
                  for i in range(k))
        magnitude = tuple(abs(value) for value in h)
        sign_h = tuple(1 if value >= 0 else -1 for value in h)
        for seed_signs in product((-1, 1), repeat=len(seed)):
            xi_seed = dict(zip(seed, seed_signs))
            field = {
                i: sum(FREE[i][j] * xi_seed[j] for j in seed)
                for i in remaining
            }
            context_total = Fraction(0)
            context_tt = Fraction(0)
            fixed_tau_tt = {-1: Fraction(0), 1: Fraction(0)}
            context_weight = Fraction(0)
            for t in (-1, 1):
                for other_signs in product((-1, 1), repeat=len(remaining)):
                    xi = dict(xi_seed)
                    xi.update(zip(remaining, other_signs))
                    local = {i: xi[i] * magnitude[i] + t * field[i]
                             for i in remaining}
                    choices = [
                        (-1, 1) if local[i] == 0
                        else ((1,) if local[i] > 0 else (-1,))
                        for i in remaining
                    ]
                    ties = sum(local[i] == 0 for i in remaining)
                    weight = Fraction(1, 2 ** (1 + len(remaining) + ties))
                    for updated in product(*choices):
                        spin = dict(xi_seed)
                        spin.update(zip(remaining, updated))
                        x = tuple(spin[i] for i in range(k))
                        value = (pinned_norm + t * qform(FREE, x)
                                 + sum(xi[i] * magnitude[i] * x[i]
                                       for i in range(k)))
                        tt = t * sum(
                            FREE[i][j] * x[i] * x[j]
                            for i in remaining for j in remaining if i < j
                        )
                        # Embed the strategy in the actual gauge diamond.
                        alpha = tuple(xi[i] * sign_h[i] for i in range(k))
                        beta_g = tuple(row_signs[j] * pinned_spin[j]
                                       for j in range(m))
                        cross_value = sum(
                            x[i] * alpha[i] * CROSS[i][j]
                            * beta_g[j] * pinned_spin[j]
                            for i in range(k) for j in range(m)
                        )
                        assert cross_value == sum(
                            xi[i] * magnitude[i] * x[i] for i in range(k)
                        )
                        tau = t * sigma
                        internal = qform(FREE, x) + tau * pinned_energy
                        diamond_at_pinned_state = max(
                            abs(internal + cross_value),
                            abs(internal - cross_value)
                        )
                        assert diamond_at_pinned_state >= value
                        context_total += weight * value
                        context_tt += weight * tt
                        fixed_tau_tt[t] += 2 * weight * tt
                        context_weight += weight
                        outcome_checks += 1
            assert context_weight == 1
            assert context_tt == 0
            if any(value != 0 for value in fixed_tau_tt.values()):
                nonzero_fixed_tau_contexts += 1
            prediction = (pinned_norm + sum(magnitude[i] for i in seed)
                          + sum(max(magnitude[i], abs(field[i]))
                                for i in remaining))
            assert context_total == prediction
            expected_total += context_total
            contexts += 1
    expected_total /= contexts
    formula = pinned_norm + len(seed) * a(m) + len(remaining) * a(m + len(seed))
    assert expected_total == formula == Fraction(43, 4)
    assert nonzero_fixed_tau_contexts > 0
    return {
        "pinned_norm": pinned_norm,
        "conditional_contexts": contexts,
        "strategy_outcomes_checked": outcome_checks,
        "conditional_tt_cancellations": contexts,
        "nonzero_fixed_tau_contexts": nonzero_fixed_tau_contexts,
        "mean_strategy_score": str(expected_total),
        "universal_formula_on_fixture": str(formula),
        "all_strategy_outcomes_bounded_by_actual_diamond": True,
    }


def main():
    started = time.monotonic()
    fixture = {"free": FREE, "pinned": PINNED, "cross": CROSS,
               "seed_set": SEED_SET}
    record = {
        "classification": "finite_exact_strategy_and_scalar_corroboration",
        "host": socket.gethostname(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fixture": fixture,
        "fixture_sha256": hashlib.sha256(json.dumps(
            fixture, sort_keys=True, separators=(",", ":")
        ).encode()).hexdigest(),
        "source_signings_checked": 1,
        "strategy": strategy_check(),
        "scalar": scalar_checks(),
        "asymptotic_proof_scope": "analytic proof note, not finite enumeration",
        "original_convergence_proved": False,
        "passed": True,
        "elapsed_seconds": time.monotonic() - started,
    }
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
