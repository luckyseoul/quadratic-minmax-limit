#!/usr/bin/env python3
"""Independent exact checks for the joint H=5p endpoint.

The new flat branch uses the actual common P in 0..9, not any catalog of
mean-2p equality forms. These bounded symbolic/arithmetic checks supplement
the proof; old graph, slice, cube, and Boolean catalogs are not recomputed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import socket
import sys
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def const(value):
    return {(0, 0): F(value)} if value else {}


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            out[exponent] = out.get(exponent, F(0))+coefficient
    return {key: value for key, value in out.items() if value}


def scale(polynomial, coefficient):
    return {key: coefficient*value for key, value in polynomial.items() if coefficient*value}


def mul(left, right):
    out = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            exponent = (i+k, j+ell)
            out[exponent] = out.get(exponent, F(0))+a*b
    return {key: value for key, value in out.items() if value}


def show(polynomial):
    return {f'p^{i} P^{j}': str(value) for (i, j), value in sorted(polynomial.items())}


def check_flat_symbolic():
    p, parallel = {(1, 0): F(1)}, {(0, 1): F(1)}
    p1 = add(p, const(1))
    m = scale(p1, F(1, 2))
    edge = scale(p, 5)
    hT = add(mul(p1, parallel), scale(p, -5))
    hard = mul(m, parallel)
    opposite = add(edge, scale(hard, -1))
    next_q = add(const(9), scale(parallel, -1))
    forbidden_q = add(const(8), scale(parallel, -1))
    opposite_mean = lambda q: add(mul(p1, q), hT, scale(p, -3))
    surplus = add(opposite, scale(mul(m, next_q), -1))
    identities = {
        'hard_edges_from_global_T': add(scale(hard, 2), scale(add(edge, hT), -1)),
        'P_plus_Q_eight_scaled_mass': add(opposite_mean(forbidden_q), const(-8)),
        'P_plus_Q_nine_scaled_mass': add(opposite_mean(next_q), scale(p, -1), const(-9)),
        'surplus_is_m_minus_five': add(surplus, scale(m, -1), const(5)),
        'forced_rows_is_five': add(m, scale(surplus, -1), const(-5)),
    }
    require(all(not value for value in identities.values()), 'symbolic coefficient identity failed')
    # p=29+s. Positivity here certifies the exact P upper bound and nonempty
    # isolated chart throughout p>=29 without a list of prime samples.
    positivity = {'isolated_vertices': [552, 48, 1],
                  'five_p_minus_nine_m': ['10', '1/2'],
                  'ten_m_minus_five_p': [5]}
    require(29**2+1-10*29 == 552 and 2*29-10 == 48, 'isolated-chart shift failed')
    require(F(5*29)-9*F(30, 2) == 10 and 5-F(9, 2) == F(1, 2), 'P=9 budget shift failed')
    require(10*F(30, 2)-5*29 == 5, 'P=10 strict bound failed')
    return {'polynomial_variables': ['p', 'P'], 'all_coefficient_residuals_zero': True,
            'identity_residuals': {name: show(value) for name, value in identities.items()},
            'hard_sign_times_T': show(hT), 'opposite_total': show(opposite),
            'surplus': show(surplus), 'shifted_positive_polynomials_at_p29': positivity,
            'common_actual_P_range': list(range(10)), 'forced_mass': 'p+9',
            'forced_row_count': 5, 'equality_catalog_used': False,
            'P9_caution': 'Q=-1 is not an actual row; P=9 and Q>=0 already imply P+Q>=9',
            'classification': 'exact symbolic row identities, conditional on common-row normalization and local mass exclusions'}


def branch(p, name, u, low_quotient, parallel, forbidden_q, forbidden_mass, mass_offset, forced_count):
    m, q, edge = (p+1)//2, (p-1)//2, 5*p
    quotient_sum = p-u
    carry = quotient_sum-m*low_quotient
    mean = 2*u+(p+1)*low_quotient
    hT = (p+1)*parallel-3*p-mean
    hard = m*parallel+carry
    opposite = edge-hard
    next_q = forbidden_q+1
    low_mass = (p+1)*forbidden_q+hT-3*p
    next_mass = (p+1)*next_q+hT-3*p
    surplus = opposite-m*next_q
    require(carry >= 0 and 2*hard == edge+hT, 'common hard-row identity failed')
    require((low_mass, next_mass, m-surplus) == (forbidden_mass, p+mass_offset, forced_count), 'carried opposite ledger failed')
    require(0 < low_mass < p-3 and 0 <= surplus < m, 'mass or pigeonhole gate failed')
    candidates = [v for v in range(edge//m+1) if (v-parallel) % q == 0]
    require(candidates == [parallel], 'carried offset not rigid in edge budget')
    return {'branch': name, 'u': u, 'low_quotient': low_quotient,
            'quotient_sum': quotient_sum, 'quotient_carry': carry,
            'low_quotient_rows_at_least_if_all_at_least_low': m-carry,
            'hard_low_mean': mean, 'parallel_offset': parallel, 'hard_edges': hard,
            'opposite_edges': opposite, 'h_times_T': hT,
            'forbidden_Q': forbidden_q, 'forbidden_mass': low_mass,
            'forced_Q': next_q, 'forced_mass': next_mass,
            'surplus': surplus, 'forced_rows_at_least': m-surplus}


def check_p1_carries():
    rows = []
    for p in (29, 37, 41, 53):
        q = (p-1)//2
        specs = [('old_complement_literal', 0, 1, 5, 2, 6, 7, 4),
                 ('carried_complement_triple', q-3, 1, 2, 6, 14, 15, 8),
                 ('carried_XNOR_omitted_pair', q-2, 1, 3, 5, 12, 13, 7),
                 ('carried_XNOR_all_equal_triple', q-2, 1, 5, 3, 12, 13, 7),
                 ('carried_offset_four', q-1, 1, 4, 4, 10, 11, 6),
                 ('carried_offset_six', q-1, 1, 6, 2, 10, 11, 6),
                 ('quotient_zero_XNOR', q, 0, 4, 3, 8, 9, 5)]
        values = [branch(p, *spec) for spec in specs]
        require([v['forced_rows_at_least'] for v in values] == [4, 8, 7, 7, 6, 6, 5], 'p1 carried row counts changed')
        rows.append({'p': p, 't': q, 'k': 5*p-1, 'H_edge_count': 5*p, 'branches': values})
    return {'rows': rows, 'classification': 'independent representative-prime p1 carried-ledger arithmetic; prior local theorems are dependencies'}


def check_p3_carries():
    rows = []
    for p in (31, 43, 47, 59):
        q, m = (p-1)//2, (p+1)//2
        sharp = [branch(p, f'sharp_offset_{parallel}', q-2, 1, parallel,
                        8-parallel, 12, 13, 7) for parallel in (2, 3, 4, 5)]
        zeros = [branch(p, 'quotient_zero_XNOR', q, 0, 4, 3, 8, 9, 5),
                 branch(p, 'quotient_zero_complement_literal', q, 0, 3, 4, 8, 9, 5)]
        require(all(row['quotient_carry'] == 2 and row['low_quotient_rows_at_least_if_all_at_least_low'] == m-2 for row in sharp), 'p3 two-unit carry counted as one high row')
        # u=q-1 has at least m-1 quotient-one rows if no quotient-zero row;
        # each endpoint baseline leaves the already excluded local mass p-1.
        low_mean = 2*(q-1)+(p+1)
        quotient_sum = p-(q-1)
        require((low_mean-(p-1), 2*m-quotient_sum) == (p-1, m-1), 'p3 mass p-1 local branch arithmetic failed')
        rows.append({'p': p, 'sharp_branches': sharp, 'zero_quotient_branches': zeros,
                     'p_minus_one_branch': {'u': q-1, 'low_mean': low_mean, 'baseline_mean': p-1,
                                           'lift_mass': p-1, 'low_rows_at_least': m-1}})
    return {'rows': rows, 'classification': 'independent representative-prime p3 carry and zero-quotient arithmetic; no graph or local-equality census'}


def distribution(p, b):
    m = (p+1)//2
    return {r: F(comb(b, r)*comb(p-b, m-r), comb(p, m))
            for r in range(max(0, m-(p-b)), min(b, m)+1)}


def phase_zero_floor(p, b):
    """Independent tiny rational quadratic-majorant LP over parity data."""
    law = distribution(p, b)
    if len(law) < 3:
        return 2*p*sum(prob*(r % 2) for r, prob in law.items())
    best = None
    for nodes in combinations(law, 3):
        coefficients = [F(0)]*3
        for t in nodes:
            left, right = [v for v in nodes if v != t]
            weight = F(t % 2, (t-left)*(t-right))
            coefficients = [a+weight*b for a, b in zip(coefficients, (left*right, -left-right, 1))]
        values = {r: coefficients[0]+coefficients[1]*r+coefficients[2]*r*r for r in law}
        if all(values[r] >= r % 2 for r in law):
            value = 2*p*sum(law[r]*values[r] for r in law)
            best = value if best is None else min(best, value)
    require(best is not None, 'phase-zero LP has no feasible vertex')
    return best


def check_arm_floors():
    floor_rows = []
    for p in (29, 31):
        floors = {b: phase_zero_floor(p, b) for b in range(0, p, 2)}
        survivors = [b for b, floor in floors.items() if floor <= p+9]
        require(survivors == [0, 2, p-1], 'unexpected mass p+9 boundary')
        require([b for b, floor in floors.items() if floor <= 8] == [0], 'unexpected mass8 boundary')
        require(all(0 < p+9-floors[b] < p-3 for b in survivors if b), 'nonzero-boundary excess not subsharp')
        m = (p+1)//2
        omitted_bit_baseline = [(m-bit) % 2 for bit in (0, 1)]
        omitted_mean = F(p-m, p)*omitted_bit_baseline[0]+F(m, p)*omitted_bit_baseline[1]
        require(2*p*omitted_mean == floors[p-1], 'pointwise omitted-bit parity baseline mean mismatch')
        floor_rows.append({'p': p, 'phase_zero_LP_floors': {str(b): str(floor) for b, floor in floors.items()},
                           'mass_p_plus_nine_surviving_boundaries': survivors,
                           'omitted_bit_pointwise_parity_baseline_values': omitted_bit_baseline})
    flat_rows = []
    for p in (29, 31, 37, 41, 43, 53):
        m = (p+1)//2
        for parallel in range(10):
            threshold = 9-parallel
            opposite = 5*p-m*parallel
            surplus = opposite-m*threshold
            require((surplus, m-surplus) == (m-5, 5), 'arbitrary-P flat pigeonhole failed')
            require(all((p+1)*(parallel+Q)-8*p <= 8 for Q in range(max(0, threshold))), 'lower Q escaped mass8-or-negative split')
            require((p+1)*(parallel+threshold)-8*p == p+9, 'flat threshold mass failed')
            flat_rows.append({'p': p, 'actual_P': parallel, 'least_allowed_Q': threshold,
                              'forbidden_Q_is_legal': 8-parallel >= 0,
                              'forced_mass': p+9, 'forced_rows_at_least': m-surplus})
    return {'independent_floor_rows': floor_rows, 'flat_rows': flat_rows,
            'P_values_are_actual_counts_not_equality_offsets': True,
            'analytic_dependencies_not_reproved': ['phase-zero pointwise parity baseline subtraction', 'sharp integral lift floor', 'local p+9 theorem'],
            'classification': 'ARM64 rational-floor and integer-ledger portability, not an all-prime proof'}


MODES = {'flat-symbolic': check_flat_symbolic, 'p1-carries': check_p1_carries,
         'p3-carries': check_p3_carries, 'arm-floors': check_arm_floors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=MODES)
    args = parser.parse_args()
    result = MODES[args.mode]()
    print(json.dumps({'mode': args.mode, 'host': socket.gethostname(), 'architecture': platform.machine(),
                      'python': sys.version.split()[0], 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'passed': True, 'classification': 'independent exact joint H=5p bridge check', 'result': result}, sort_keys=True))


if __name__ == '__main__':
    main()
