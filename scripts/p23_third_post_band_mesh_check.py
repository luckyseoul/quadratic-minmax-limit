#!/usr/bin/env python3
"""Independent bounded bridge checks for the p23,t11 endpoint.

Standard library only; no imports from the candidate proposition. These
arithmetic/certificate checks supplement the analytic proof, not replace it.
The contact-matrix certificate samples rows until its target rank is reached;
it does not enumerate the Johnson slice or repeat the five-set certificate.
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
from itertools import combinations
from math import comb
from pathlib import Path


N, K, PRIME = 23, 12, 1_000_003
CONTACTS = {
    4: ([0, 2, 4], [F(2, 23), F(18, 23), F(3, 23)]),
    6: ([0, 2, 4, 6], [F(1, 92), F(39, 92), F(51, 92), F(1, 92)]),
    8: ([0, 2, 4, 6, 8], [F(3, 115), F(3, 115), F(93, 115), F(13, 115), F(3, 115)]),
    10: ([0, 2, 4, 6, 8, 10], [F(3, 368), F(3, 368), F(3, 8), F(109, 184), F(3, 368), F(3, 368)]),
    12: ([2, 4, 6, 8, 10, 12], [F(3, 184), F(3, 184), F(77, 92), F(9, 92), F(3, 184), F(3, 184)]),
    14: ([4, 6, 8, 10, 12], [F(3, 230), F(83, 230), F(3, 5), F(3, 230), F(3, 230)]),
    16: ([6, 8, 10, 12], [F(1, 23), F(18, 23), F(3, 23), F(1, 23)]),
    18: ([8, 10, 12], [F(15, 46), F(15, 23), F(1, 46)]),
    20: ([10, 12], [F(18, 23), F(5, 23)]),
}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def distribution(b):
    return {r: F(comb(b, r) * comb(N-b, K-r), comb(N, K))
            for r in range(max(0, K-(N-b)), min(b, K)+1)}


def interpolate_parity(nodes, phase):
    """Lagrange interpolation in the coefficient basis 1,r,r^2."""
    output = [F(0)] * 3
    for t in nodes:
        others = [u for u in nodes if u != t]
        denominator = (t-others[0]) * (t-others[1])
        scale = F((t+phase) % 2, denominator)
        terms = [others[0]*others[1], -sum(others), 1]
        output = [a + scale*v for a, v in zip(output, terms)]
    return output


def exact_floor(b, phase):
    """Tiny independent rational LP: min 46 E[q], q(r)>=parity(r)."""
    law = distribution(b)
    if len(law) < 3:
        return 2*N*sum(prob*((r+phase) % 2) for r, prob in law.items())
    candidates = []
    for nodes in combinations(law, 3):
        coefficients = interpolate_parity(nodes, phase)
        values = {r: sum(c*r**j for j, c in enumerate(coefficients)) for r in law}
        if all(values[r] >= (r+phase) % 2 for r in law):
            candidates.append(2*N*sum(law[r]*values[r] for r in law))
    require(bool(candidates), 'parity LP has no vertex')
    return min(candidates)


def check_contacts():
    rows = []
    for b, (nodes, weights) in CONTACTS.items():
        law = distribution(b)
        moments = [sum(prob*r**j for r, prob in law.items()) for j in range(3)]
        reproduced = [sum(w*r**j for r, w in zip(nodes, weights)) for j in range(3)]
        require(all(w > 0 for w in weights), 'nonpositive contact weight')
        require(all(r in law and r % 2 == 0 for r in nodes), 'invalid contact node')
        require(moments == reproduced and moments[0] == 1, 'quadrature mismatch')
        floor = exact_floor(b, 1)
        require(floor == 46, 'phase-one floor is not 46')
        rows.append({'b': b, 'nodes': nodes, 'weights': list(map(str, weights)),
                     'moments': list(map(str, moments)), 'independent_LP_floor': str(floor)})
    return {'rows': rows, 'positive_degree_two_quadratures': len(rows)}


def monomials(d):
    return [0] + [1 << i for i in range(d)] + [(1 << i) | (1 << j) for i, j in combinations(range(d), 2)]


def add_pivot(vector, pivots, modulus=PRIME):
    vector = [x % modulus for x in vector]
    for column, pivot in sorted(pivots.items()):
        factor = vector[column]
        if factor:
            vector = [(a-factor*b) % modulus for a, b in zip(vector, pivot)]
    leading = next((i for i, value in enumerate(vector) if value), None)
    if leading is None:
        return False
    inverse = pow(vector[leading], -1, modulus)
    pivots[leading] = [value*inverse % modulus for value in vector]
    return True


def kernel_vectors(d):
    """24 slice-ideal vectors and three active-side even-contact kernels."""
    columns = monomials(N)
    index = {mask: j for j, mask in enumerate(columns)}
    vectors = []
    base = [0]*len(columns)
    base[0] = -K
    for i in range(N):
        base[index[1 << i]] = 1
    vectors.append(base)
    for i in range(N):
        row = [0]*len(columns)
        row[index[1 << i]] = 1-K
        for j in range(N):
            if i != j:
                row[index[(1 << i) | (1 << j)]] = 1
        vectors.append(row)
    if d == 4:
        for pair in ((0, 1), (0, 2), (0, 3)):
            complement = tuple(i for i in range(4) if i not in pair)
            row = [0]*len(columns)
            for sign, (i, j) in ((1, pair), (-1, complement)):
                row[index[1 << i]] -= 2*sign
                row[index[1 << j]] -= 2*sign
                row[index[(1 << i) | (1 << j)]] += 4*sign
            vectors.append(row)
    elif d == 3:
        for i in range(3):
            j, k = [v for v in range(3) if v != i]
            row = [0]*len(columns)
            row[index[1 << i]] = 2
            row[index[1 << j]] = row[index[1 << k]] = -2
            row[index[(1 << j) | (1 << k)]] = 4
            vectors.append(row)
    else:
        raise ValueError('small side must have size 3 or 4')
    return vectors


def small_kernel(d):
    columns = monomials(N)
    vectors = kernel_vectors(d)
    kernel_pivots = {}
    for vector in vectors:
        add_pivot(vector, kernel_pivots)
    require(len(kernel_pivots) == 27, 'kernel witnesses not independent')
    active_kernels = vectors[24:]
    for mask in range(1 << d):
        if mask.bit_count() % 2 == 0:
            row = [int(mask & m == m) for m in columns]
            require(all(sum(a*b for a, b in zip(row, v)) == 0 for v in active_kernels), 'active kernel witness does not vanish')
    rng = random.Random(1577100+d)
    pivots, selected, examined = {}, [], 0
    ranks_needed = len(columns)-len(vectors)
    while len(pivots) < ranks_needed and examined < 1500:
        examined += 1
        r = rng.choice(list(range(0, d+1, 2)))
        chosen = rng.sample(range(d), r) + rng.sample(range(d, N), K-r)
        mask = sum(1 << v for v in chosen)
        row = [int(mask & m == m) for m in columns]
        require(all(sum(a*b for a, b in zip(row, v)) == 0 for v in vectors), 'full contact annihilation failed')
        if add_pivot(row, pivots):
            selected.append(mask)
    require(len(pivots) == 250, 'contact rows did not reach rank 250')
    return {'active_side': d, 'raw_coefficients': len(columns), 'prime': PRIME,
            'explicit_independent_kernel_vectors': 27, 'slice_ideal_dimension': 24,
            'extra_contact_kernel_dimension': 3, 'sampled_rows_examined': examined,
            'contact_rank': len(pivots), 'independent_contact_masks': selected,
            'certificate_classification': 'finite exact matrix-rank certificate; not a slice census'}


def check_small_kernel():
    return {'boundaries': {'4': small_kernel(4), '20': small_kernel(3)}}


def check_cubes():
    rows = []
    for b in range(6, 20, 2):
        d = min(b, N-b)
        capacities = []
        covered = 0
        for r in range(d+1):
            available_unselected = N-d-(K-r)
            available_selected = K-r
            require(available_unselected >= r and available_selected >= d-r, 'swap pairing unavailable')
            covered += comb(d, r)*comb(N-d, K-r)
            capacities.append([r, available_unselected, r, available_selected, d-r, K-d])
        require(covered == comb(N, K), 'covering layer count mismatch')
        chars = monomials(d)
        differences = {a ^ c for a in chars for c in chars}
        vertices = [mask for mask in range(1 << d) if mask.bit_count() % 2 == 0]
        sums = {delta: sum((-1)**((mask & delta).bit_count()) for mask in vertices) for delta in differences}
        require(sums[0] == 2**(d-1) and all(value == 0 for delta, value in sums.items() if delta), 'Walsh Gram not diagonal')
        contact_original = sorted({r if b <= 11 else K-r for r in range(0, d+1, 2)})
        require(contact_original == CONTACTS[b][0], 'cube parity/contact mismatch')
        rows.append({'b': b, 'dimension': d, 'capacity_rows_r_unselected_need_selected_need_filler': capacities,
                     'covered_slice_configurations_via_Vandermonde': covered,
                     'contact_layers': contact_original, 'quadratic_characters': len(chars),
                     'Walsh_Gram_diagonal': 2**(d-1), 'Walsh_Gram_off_diagonal': 0})
    return {'rows': rows, 'construction': 'pair each small-side coordinate with an opposite-membership outside coordinate; keep the other 12-d selected coordinates fixed'}


def check_ledger():
    floors = {b: exact_floor(b, 0) for b in range(0, 23, 2)}
    require([b for b, floor in floors.items() if floor <= 32] == [0, 2, 22], 'unexpected phase-zero mass32 boundary')
    require(floors[0] == 0 and floors[2] == floors[22] == 24, 'phase-zero exceptional floors changed')
    require(0 < 32-24 < 20, 'boundary baseline excess no longer subsharp')
    offsets = []
    for offset in range(4, 9):
        candidates = [v for v in range(115//12+1) if (v-offset) % 11 == 0]
        require(candidates == [offset], 'offset rigidity changed')
        q0, q1 = 8-offset, 9-offset
        mass0, mass1 = 24*(offset+q0)-184, 24*(offset+q1)-184
        surplus = 115-12*offset-12*q1
        require((mass0, mass1, surplus, 12-surplus) == (8, 32, 7, 5), 'offset mass ledger changed')
        offsets.append({'offset_P': offset, 'forbidden_Q': q0, 'forced_Q': q1,
                        'forbidden_mass': mass0, 'forced_mass': mass1, 'surplus': surplus, 'forced_rows_at_least': 12-surplus})
    target = F(32, 92)
    junta = F(6*22*21*32, 23*23*24)
    require(F(23-7, 4) > 3 and F(32, 48) < F(3, 4), 'p+9 height numeric gate failed')
    require(junta < 7 and F(6, 23) < target < F(11, 23), 'p+9 Boolean numeric gate failed')
    return {'independent_phase_zero_LP_floors': {str(b): str(floor) for b, floor in floors.items()},
            'mass32_floor_surviving_boundaries': [0, 2, 22], 'nonzero_boundary_lift_excess': 8,
            'sharp_lift_floor_dependency': 20, 'offset_ledgers': offsets,
            'p_plus_nine_target_density': str(target), 'p_plus_nine_junta_bound': str(junta),
            'p_plus_nine_height_lower': '4', 'p_plus_nine_paired_mean_upper': '2/3',
            'analytic_dependencies_not_reproved': ['sharp lift floor', 'p+9 local exclusion', 'common-row normalization']}


MODES = {'contacts': check_contacts, 'small-kernel': check_small_kernel,
         'cubes': check_cubes, 'ledger': check_ledger}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=MODES)
    args = parser.parse_args()
    result = MODES[args.mode]()
    print(json.dumps({'mode': args.mode, 'host': socket.gethostname(), 'architecture': platform.machine(),
                      'python': sys.version.split()[0], 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'passed': True, 'classification': 'independent bounded bridge cross-check', 'result': result}, sort_keys=True))


if __name__ == '__main__':
    main()
