#!/usr/bin/env python3
"""Independent bounded checks for the third p=1 mod4 post-band layer.

The fresh premise is the complement-triple punctured gap: the difference
from (r-2)^2 need not be nonnegative at r=0. No candidate-module imports,
Johnson-slice enumeration, previous p23 replay, or GPU computation is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import socket
import sys
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path


PRIME = 1_000_003


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def hypergeometric(n, k, active):
    return {r: F(comb(active, r)*comb(n-active, k-r), comb(n, k))
            for r in range(max(0, k-(n-active)), min(active, k)+1)}


def polynomial_product(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def check_quadrature():
    # With h=p-29, (p-5)(p-7)-16(p-4)=h^2+30h+128.
    numerator = polynomial_product([24, 1], [22, 1])
    numerator[0] -= 16*25
    numerator[1] -= 16
    require(numerator == [128, 30, 1], 'shifted strict-gap polynomial changed')
    rows = []
    for p in (29, 31, 33, 37, 41, 53, 101):
        m, n = (p+1)//2, p-3
        law = hypergeometric(p, m, 3)
        moments = [sum(prob*r**j for r, prob in law.items()) for j in range(3)]
        weights = [F(3*(p-3), 4*p), F(3, p), F(p-3, 4*p)]
        quadrature = [sum(w*r**j for r, w in zip((1, 2, 3), weights)) for j in range(3)]
        baseline = 2*p*sum(prob*(r-2)**2 for r, prob in law.items())
        require(moments == quadrature and all(w > 0 for w in weights), 'positive quadrature failed')
        require(baseline == 2*p-6, 'baseline scaled mean mismatch')
        lower = F((n-2)*(n-4), 4*n*(n-1))
        odd_floor = F(n-4, 4*(n-1))
        require(lower == F(n-2, n)*odd_floor, 'neighboring-slice averaging factor mismatch')
        section_average_ratios = []
        for k, selected in ((n//2+1, True), (n//2-1, False)):
            occurrence = k if selected else n-k
            section_size = comb(n-1, k-1) if selected else comb(n-1, k)
            ratio = F(occurrence*comb(n, k), n*section_size)
            require(ratio == 1, 'uniform section average is not the whole-slice average')
            require(k >= 3 and n-k >= 3, 'three-swap support argument lacks capacity')
            section_average_ratios.append(str(ratio))
        weighted_lower = (p-3)*lower
        require(weighted_lower == F((p-5)*(p-7), 4*(p-4)) > 4, 'punctured gap not strict')
        rows.append({'odd_order': p, 'prime_hypothesis_used': False,
                     'moments': list(map(str, moments)), 'quadrature_weights': list(map(str, weights)),
                     'baseline_scaled_mean': str(baseline), 'outside_order': n,
                     'neighboring_slice_nonzero_mean_lower_bound': str(lower),
                     'section_average_ratios': section_average_ratios,
                     'single_nonzero_contact_contribution_lower_bound': str(weighted_lower),
                     'strict_gap_above_four': str(weighted_lower-4)})
    return {'shifted_gap_numerator_coefficients': numerator, 'rows': rows,
            'analytic_dependency': 'odd-middle-slice integral quadratic floor and at most two identically zero sections'}


def monomials(n):
    return [0] + [1 << i for i in range(n)] + [(1 << i) | (1 << j) for i, j in combinations(range(n), 2)]


def add_pivot(vector, pivots):
    vector = [v % PRIME for v in vector]
    for column, pivot in sorted(pivots.items()):
        factor = vector[column]
        if factor:
            vector = [(a-factor*b) % PRIME for a, b in zip(vector, pivot)]
    leading = next((i for i, value in enumerate(vector) if value), None)
    if leading is None:
        return False
    inverse = pow(vector[leading], -1, PRIME)
    pivots[leading] = [v*inverse % PRIME for v in vector]
    return True


def odd_contact_kernels(n, k):
    columns = monomials(n)
    index = {mask: j for j, mask in enumerate(columns)}
    vectors = []
    base = [0]*len(columns)
    base[0] = -k
    for i in range(n):
        base[index[1 << i]] = 1
    vectors.append(base)
    for i in range(n):
        row = [0]*len(columns)
        row[index[1 << i]] = 1-k
        for j in range(n):
            if i != j:
                row[index[(1 << i) | (1 << j)]] = 1
        vectors.append(row)
    # On the three-bit odd half, z_i-z_j*z_k vanishes, z_i=2*x_i-1.
    for i in range(3):
        j, k_other = [v for v in range(3) if v != i]
        row = [0]*len(columns)
        row[0] = -2
        for v in range(3):
            row[index[1 << v]] = 2
        row[index[(1 << j) | (1 << k_other)]] = -4
        vectors.append(row)
    return vectors


def check_contact_kernel():
    n, k = 29, 15
    columns = monomials(n)
    vectors = odd_contact_kernels(n, k)
    kernel_pivots = {}
    for vector in vectors:
        add_pivot(vector, kernel_pivots)
    require(len(kernel_pivots) == 33, 'explicit kernels not independent')
    for mask in range(8):
        if mask.bit_count() % 2 == 1:
            row = [int(mask & m == m) for m in columns]
            require(all(sum(a*b for a, b in zip(row, vector)) == 0 for vector in vectors[-3:]), 'odd-contact character kernel mismatch')
    rng = random.Random(15772029)
    pivots, selected, examined = {}, [], 0
    while len(pivots) < 403 and examined < 2500:
        examined += 1
        r = rng.choice((1, 3))
        chosen = rng.sample(range(3), r) + rng.sample(range(3, n), k-r)
        mask = sum(1 << v for v in chosen)
        row = [int(mask & m == m) for m in columns]
        require(all(sum(a*b for a, b in zip(row, vector)) == 0 for vector in vectors), 'contact sample not annihilated by kernel')
        if add_pivot(row, pivots):
            selected.append(mask)
    require(len(pivots) == 403, 'contact rank lower certificate incomplete')
    return {'p': n, 'slice_weight': k, 'original_small_side_contact_layers': [1, 3],
            'raw_coefficients': len(columns), 'prime': PRIME,
            'slice_ideal_kernel_dimension': 30, 'odd_contact_extra_kernel_dimension': 3,
            'explicit_independent_kernel_vectors': len(kernel_pivots),
            'contact_rank': len(pivots), 'sampled_rows_examined': examined,
            'independent_contact_masks': selected,
            'classification': 'finite exact coefficient-rank certificate supplementing the general slice-kernel proof; not a slice census'}


def mobius_coefficients(values):
    output = list(values)
    for bit in range(3):
        for mask in range(8):
            if mask & (1 << bit):
                output[mask] -= output[mask ^ (1 << bit)]
    return output


def walsh_coefficients(values):
    return [F(sum(value * (-1)**(mask.bit_count()-(mask & vertex).bit_count())
                  for vertex, value in enumerate(values)), 8) for mask in range(8)]


def check_cube_equality():
    pairs = [3, 5, 6]
    rows = []
    for pair_values in product((0, 1), repeat=3):
        if sum(pair_values) > 1:
            continue
        difference = [0]*8
        difference[0] = -sum(pair_values)
        for mask, value in zip(pairs, pair_values):
            difference[mask] = value
        values = [(mask.bit_count()-2)**2+2*difference[mask] for mask in range(8)]
        delta = 4*sum(pair_values)
        require(mobius_coefficients(values)[7] == 0, 'equality truth table not quadratic')
        require(all(v >= 0 and v % 2 == mask.bit_count() % 2 for mask, v in enumerate(values)), 'nonnegativity or parity failed')
        signed = walsh_coefficients([3+2*v for v in values])
        require(signed[7] == 0 and all(v.denominator == 1 for v in signed), 'signed quadratic integrality failed')
        offset = signed[0]+sum(signed[1 << i] for i in range(3))
        if delta == 0:
            require(offset == 2, 'baseline offset changed')
        else:
            chosen_pair = pairs[pair_values.index(1)]
            missing_bit = 7 ^ chosen_pair
            expected = [F(0)]*8
            expected[0], expected[chosen_pair], expected[missing_bit] = F(5), F(1), F(-1)
            require(signed == expected and offset == 4, 'delta4 signed target is not 5+z_i*z_j-z_k')
        for p in (29, 33, 37, 53):
            law = hypergeometric(p, (p+1)//2, 3)
            mean = sum(law[r]*F(sum(values[mask] for mask in range(8) if mask.bit_count() == r), comb(3, r)) for r in law)
            require(2*p*mean == 2*p-6+delta, 'equality table mean mismatch')
        rows.append({'delta': delta, 'pair_values': list(pair_values), 'difference_truth_table': difference,
                     'A_truth_table': values, 'signed_target_Walsh_coefficients_mask_order': list(map(str, signed)),
                     'signed_offset': int(offset)})
    require(len(rows) == 4 and [row['delta'] for row in rows].count(4) == 3, 'equality form count changed')
    return {'forms': rows, 'delta_two_excluded_by': 'delta=4*sum(nonnegative integer pair values)',
            'negative_difference_at_empty_set_is_retained': True,
            'classification': 'complete residual three-bit equality check after the separately proved globalization'}


def check_branch_ledger():
    rows = []
    for p in (29, 37, 41, 53):
        m, q, edge = (p+1)//2, (p-1)//2, 5*p-2
        specs = [('literal', p+1, 5, m-2, 2, 6, 7, 5),
                 ('XNOR_zero_quotient', p-1, 4, m-1, 3, 8, 9, 6),
                 ('triple_delta0_two_unit_carry', 2*p-6, 2, 2, 6, 14, 15, 9),
                 ('XNOR_sharp_P3_one_unit_carry', 2*p-4, 3, 1, 5, 12, 13, 8),
                 ('XNOR_sharp_P5_one_unit_carry', 2*p-4, 5, 1, 3, 12, 13, 8),
                 ('literal_sharp_P4', 2*p-2, 4, 0, 4, 10, 11, 7),
                 ('literal_sharp_P6', 2*p-2, 6, 0, 2, 10, 11, 7),
                 ('triple_delta4_P4', 2*p-2, 4, 0, 4, 10, 11, 7)]
        prime_rows = []
        for name, low, parallel, carry, forbidden_q, forbidden_mass, mass_offset, forced_count in specs:
            hT = (p+1)*parallel-3*p-low
            hard, opposite = m*parallel+carry, edge-(m*parallel+carry)
            require(2*hard == edge+hT, 'hard edge identity failed')
            allowed_parallel = [v for v in range(edge//m+1) if (v-parallel) % q == 0]
            require(allowed_parallel == [parallel], 'parallel coefficient congruence not rigid')
            forced_q = forbidden_q+1
            actual_forbidden = (p+1)*forbidden_q+hT-3*p
            actual_forced = (p+1)*forced_q+hT-3*p
            surplus = opposite-m*forced_q
            require((actual_forbidden, actual_forced, surplus) == (forbidden_mass, p+mass_offset, m-forced_count), 'opposite row ledger mismatch')
            require(0 < forbidden_mass < p-3 and 0 <= surplus < m, 'floor or pigeonhole gate failed')
            floor_rows = {b: (0 if b == 0 else p+1 if b == 2 else p-1 if b == p-1 else 2*p-6 if b == 4 else 2*p) for b in range(0, p, 2)}
            nonzero = [(b, floor, actual_forced-floor) for b, floor in floor_rows.items() if b and floor <= actual_forced]
            require([v[0] for v in nonzero] == [2, p-1] and all(0 < v[2] < p-3 for v in nonzero), 'opposite nonzero-boundary floor gate failed')
            prime_rows.append({'branch': name, 'hard_low_mean': low, 'parallel_offset': parallel, 'quotient_carry': carry,
                               'hard_edges': hard, 'opposite_edges': opposite, 'h_times_T': hT,
                               'forbidden_Q': forbidden_q, 'forbidden_mass': actual_forbidden,
                               'forced_Q': forced_q, 'forced_mass': actual_forced, 'surplus': surplus,
                               'forced_rows_at_least': m-surplus, 'nonzero_boundary_floor_and_lift_rows': nonzero})
        rows.append({'p': p, 't': q-1, 'original_k': edge-1, 'H_edge_count': edge,
                     'guaranteed_isolated_vertices': p*p+1-2*edge, 'branches': prime_rows})
    return {'rows': rows, 'analytic_dependencies_not_reproved': ['isolated-chart common row and coefficient congruence', 'phase-zero exact floor table', 'sharp integral lift floor', 'local p+7,p+9,p+11,p+13,p+15 exclusions'],
            'classification': 'representative-prime exact integer/Fraction portability check, not an all-prime theorem'}


MODES = {'quadrature': check_quadrature, 'contact-kernel': check_contact_kernel,
         'cube-equality': check_cube_equality, 'branch-ledger': check_branch_ledger}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=MODES)
    args = parser.parse_args()
    result = MODES[args.mode]()
    print(json.dumps({'mode': args.mode, 'host': socket.gethostname(), 'architecture': platform.machine(),
                      'python': sys.version.split()[0], 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'passed': True, 'classification': 'independent bounded check of the new punctured-gap and third-layer bridge',
                      'result': result}, sort_keys=True))


if __name__ == '__main__':
    main()
