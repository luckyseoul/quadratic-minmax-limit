#!/usr/bin/env python3
"""Exact N=6 balanced-profile catalog; not an all-orders theorem.

Enumerate precisely 1024 first-row-positive signings and all 64 spin states
using direct edge loops. The pressure numerator is cosh(v)*P(X,Y), with
X=cosh(2u), Y=cosh(2v). Integer coefficientwise comparison after
X=1+p+q, Y=1+q gives sufficient dominance on u>=v>=0.
The undominated set is a sufficient candidate set, not automatically the
exact envelope winner set. All mathematical arithmetic is exact.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import socket
import time


N = 6
EDGES = tuple(itertools.combinations(range(N), 2))
FREE_EDGES = tuple((i, j) for i, j in EDGES if i != 0)
STATES = tuple(itertools.product((-1, 1), repeat=N))


def clean(poly):
    return {key: value for key, value in poly.items() if value}


def add(left, right, scale=1):
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, 0) + scale * value
    return clean(result)


def multiply(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, l), b in right.items():
            key = (i + k, j + l)
            result[key] = result.get(key, 0) + a * b
    return clean(result)


def scale(poly, factor):
    return clean({key: factor * value for key, value in poly.items()})


def evaluate(poly, x, y):
    return sum(value * x ** i * y ** j for (i, j), value in poly.items())


def terms(poly):
    return [[i, j, value] for (i, j), value in sorted(poly.items())]


def powers(base, degree):
    result = [{(0, 0): 1}]
    for _ in range(degree):
        result.append(multiply(result[-1], base))
    return result


def basis():
    one = {(0, 0): 1}
    x = {(1, 0): 1}
    y = {(0, 1): 1}
    ts = [one, x]
    for _ in range(2, 4):
        ts.append(add(scale(multiply(x, ts[-1]), 2), ts[-2], -1))
    rs = [one, {(0, 1): 2, (0, 0): -1}]
    for _ in range(2, 5):
        rs.append(add(scale(multiply(y, rs[-1]), 2), rs[-2], -1))
    return ts, rs


def cone_transform(poly):
    xp = powers({(0, 0): 1, (1, 0): 1, (0, 1): 1}, 3)
    yp = powers({(0, 0): 1, (0, 1): 1}, 4)
    result = {}
    for (i, j), value in poly.items():
        result = add(result, multiply(xp[i], yp[j]), value)
    return result


def matrix_from_mask(mask):
    matrix = [[0] * N for _ in range(N)]
    for j in range(1, N):
        matrix[0][j] = matrix[j][0] = 1
    for bit, (i, j) in enumerate(FREE_EDGES):
        matrix[i][j] = matrix[j][i] = -1 if mask & (1 << bit) else 1
    return matrix


def direct_joint(matrix):
    signed = Counter()
    for state in STATES:
        internal = 0
        cross = 0
        for i in range(N):
            for j in range(i + 1, N):
                contribution = matrix[i][j] * state[i] * state[j]
                if (i < 3) == (j < 3):
                    internal += contribution
                else:
                    cross += contribution
        signed[(internal, cross)] += 1
    absolute = Counter()
    for (internal, cross), count in signed.items():
        absolute[(abs(internal), abs(cross))] += count
    return signed, tuple(sorted((i, j, count) for (i, j), count in absolute.items()))


def polynomial(signature, ts, rs):
    result = {}
    for internal, cross, count in signature:
        result = add(result, multiply(ts[internal // 2], rs[(cross - 1) // 2]), count)
    return result


def cosh_power(base, power):
    return (base ** power + base ** (-power)) / 2


def nonnegative_nonzero(poly):
    return bool(poly) and all(value >= 0 for value in poly.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proof', required=True, type=Path)
    parser.add_argument('--expected-proof-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    proof_sha = hashlib.sha256(args.proof.read_bytes()).hexdigest()
    if proof_sha != args.expected_proof_sha256:
        raise AssertionError('proof hash mismatch')
    if args.output.exists():
        raise FileExistsError(args.output)
    started = time.monotonic()
    checks = Counter()

    def require(label, condition):
        checks[label] += 1
        if not condition:
            raise AssertionError(label)

    require('fixed_scope', len(EDGES) == 15 and len(FREE_EDGES) == 10 and len(STATES) == 64)
    ts, rs = basis()
    catalog = {}
    for mask in range(1 << len(FREE_EDGES)):
        matrix = matrix_from_mask(mask)
        signed, signature = direct_joint(matrix)
        require('state_count', sum(signed.values()) == 64)
        require('block_flip_symmetry', all(signed[(i, j)] == signed[(i, -j)] for i, j in signed))
        require('internal_cross_support', all(i in (0, 2, 4, 6) and j in (1, 3, 5, 7, 9)
                                              for i, j, count in signature))
        if signature not in catalog:
            catalog[signature] = {'masks': [], 'signed': signed}
        catalog[signature]['masks'].append(mask)

    ordered = sorted(catalog, key=lambda signature: catalog[signature]['masks'][0])
    polys = []
    cones = []
    rows = []
    for class_id, signature in enumerate(ordered):
        entry = catalog[signature]
        poly = polynomial(signature, ts, rs)
        cone = cone_transform(poly)
        require('zero_temperature_numerator', evaluate(poly, 1, 1) == 64)
        for p, q in ((0, 0), (1, 0), (0, 1), (2, 3)):
            require('cone_substitution', evaluate(cone, p, q) == evaluate(poly, 1 + p + q, 1 + q))
        for a, b in ((1, 1), (2, 1), (2, 2), (3, 2)):
            eu, ev = Fraction(a), Fraction(b)
            x, y = cosh_power(eu, 2), cosh_power(ev, 2)
            direct = sum(count * (eu ** i * ev ** j + eu ** (-i) * ev ** (-j)) / 2
                         for (i, j), count in entry['signed'].items())
            factored = sum(count * cosh_power(eu, i) * cosh_power(ev, j)
                           for i, j, count in signature)
            algebraic = cosh_power(ev, 1) * evaluate(poly, x, y)
            require('rational_direct_factorization', direct == factored)
            require('rational_chebyshev_formula', direct == algebraic)
        polys.append(poly)
        cones.append(cone)
        rows.append(dict(class_id=class_id, representative_mask=entry['masks'][0],
                         mask_count=len(entry['masks']), masks=entry['masks'],
                         representative_matrix=matrix_from_mask(entry['masks'][0]),
                         absolute_joint_signature=[list(row) for row in signature],
                         xy_polynomial_terms=terms(poly), pq_polynomial_terms=terms(cone)))

    domination = []
    for target, cone in enumerate(cones):
        lowers = []
        for source, other in enumerate(cones):
            if target == source:
                continue
            difference = add(cone, other, -1)
            if nonnegative_nonzero(difference):
                lowers.append(source)
        domination.append(lowers)
    candidates = [i for i, lowers in enumerate(domination) if not lowers]
    require('candidate_set_nonempty', bool(candidates))
    certificates = []
    for target in range(len(rows)):
        if target in candidates:
            continue
        sources = [source for source in candidates if source in domination[target]]
        require('minimal_candidate_covers_elimination', bool(sources))
        source = sources[0]
        difference = add(cones[target], cones[source], -1)
        certificates.append(dict(target_class=target, source_class=source,
                                 pq_difference_terms=terms(difference)))
    pair_differences = []
    for source, target in itertools.combinations(candidates, 2):
        difference = add(cones[target], cones[source], -1)
        pair_differences.append(dict(source_class=source, target_class=target,
                                    target_minus_source_pq_terms=terms(difference)))
    require('all_masks_covered_once', sorted(mask for row in rows for mask in row['masks']) == list(range(1024)))
    require('distinct_xy_polynomials', len({tuple(sorted(poly.items())) for poly in polys}) == len(rows))
    result = dict(status='PASS', classification='exact_fixed_order_candidate_reduction_not_all_orders',
                  n=N, blocks=[3, 3], signings=1024, states_per_signing=64,
                  convention=dict(free_edges=[list(edge) for edge in FREE_EDGES],
                                  bit_one_sign=-1, first_row_sign=1,
                                  polynomial='64 E cosh(u I+v C)=cosh(v) P(X,Y)',
                                  variables='X=cosh(2u), Y=cosh(2v), p=X-Y, q=Y-1',
                                  cone='u>=v>=0, equivalently p>=0 and q>=0',
                                  coefficient_terms='[power_of_first_variable,power_of_second_variable,integer_coefficient]'),
                  joint_signature_class_count=len(rows), catalog=rows,
                  coefficientwise_undominated_classes=candidates,
                  dominance_certificates=certificates,
                  candidate_pair_differences=pair_differences,
                  checks=dict(checks), check_count=sum(checks.values()),
                  inputs=dict(proof_sha256=proof_sha,
                              script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  environment=dict(hostname=socket.gethostname(), python=platform.python_version(),
                                   workers=1, arithmetic='stdlib exact integers and Fraction'),
                  elapsed_seconds=time.monotonic() - started)
    canonical = json.dumps(result, sort_keys=True, separators=(',', ':')).encode()
    result['canonical_result_sha256_before_this_field'] = hashlib.sha256(canonical).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: result[key] for key in ('status', 'joint_signature_class_count',
                                                'coefficientwise_undominated_classes', 'check_count',
                                                'elapsed_seconds', 'inputs', 'environment')}, sort_keys=True))


if __name__ == '__main__':
    main()
