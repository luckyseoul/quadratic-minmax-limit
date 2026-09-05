#!/usr/bin/env python3
"""Fixed N=4 regression of weighted-pressure identities; not theorem evidence.

Only 64 signings, 16 spin states, c=(0.5,2), t=(0.125,0.5,0.875).
No repository imports, optimization package, larger-order census, or writes.
"""

import argparse
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
import socket


N = 4
EDGES = tuple(itertools.combinations(range(N), 2))
INTERNAL = tuple(i for i, (a, b) in enumerate(EDGES) if a // 2 == b // 2)
CROSS = tuple(i for i in range(len(EDGES)) if i not in INTERNAL)
STATES = tuple(itertools.product((-1, 1), repeat=N))
CHARACTERS = tuple(tuple(x[i] * x[j] for i, j in EDGES) for x in STATES)
SIGNINGS = tuple(itertools.product((-1, 1), repeat=len(EDGES)))
CS = (0.5, 2.0)
TS = (0.125, 0.5, 0.875)
FD_STEP = 1e-6
MIN_TOL = 1e-11
IDENTITY_TOL = 2e-11


def profile(c, t):
    ui = c * math.sqrt((2.0 - t) / N)
    uc = c * math.sqrt(t / N)
    dui = -c / (2.0 * math.sqrt(N) * math.sqrt(2.0 - t))
    duc = c / (2.0 * math.sqrt(N) * math.sqrt(t))
    return ui, uc, dui, duc


def evaluate(a, c, t):
    ui, uc, dui, duc = profile(c, t)
    pairs = []
    for chi in CHARACTERS:
        pairs.append((sum(a[e] * chi[e] for e in INTERNAL),
                      sum(a[e] * chi[e] for e in CROSS)))
    energy = tuple(ui * i + uc * j for i, j in pairs)
    denom = math.fsum(math.cosh(h) for h in energy)
    z = denom / len(STATES)
    gamma = tuple(math.fsum(chi[e] * math.sinh(h)
                           for chi, h in zip(CHARACTERS, energy)) / denom
                  for e in range(len(EDGES)))
    r = tuple(a[e] * gamma[e] for e in range(len(EDGES)))
    counts = {}
    for i, j in pairs:
        key = (abs(i), abs(j))
        counts[key] = counts.get(key, 0) + 1
    signature = tuple(sorted((i, j, count) for (i, j), count in counts.items()))
    factor_z = math.fsum(math.cosh(ui * i) * math.cosh(uc * j)
                         for i, j in pairs) / len(STATES)
    factor_ri = math.fsum(i * math.sinh(ui * i) * math.cosh(uc * j)
                          for i, j in pairs) / denom
    factor_rc = math.fsum(j * math.sinh(uc * j) * math.cosh(ui * i)
                          for i, j in pairs) / denom
    ri = math.fsum(r[e] for e in INTERNAL)
    rc = math.fsum(r[e] for e in CROSS)
    return dict(f=math.log(z), z=z, factor_z=factor_z, r=r,
                ri=ri, rc=rc, factor_ri=factor_ri, factor_rc=factor_rc,
                derivative=dui * ri + duc * rc, signature=signature)


def minimizers(table):
    best = min(row['f'] for row in table.values())
    active = tuple(a for a in SIGNINGS if table[a]['f'] <= best + MIN_TOL)
    classes = frozenset(table[a]['signature'] for a in active)
    return best, active, classes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proof', required=True, type=Path)
    parser.add_argument('--expected-proof-sha256', required=True)
    args = parser.parse_args()
    proof_sha = hashlib.sha256(args.proof.read_bytes()).hexdigest()
    if proof_sha != args.expected_proof_sha256:
        raise AssertionError('proof input hash mismatch')
    worst = {}
    check_count = 0

    def check_close(name, lhs, rhs, tol=IDENTITY_TOL):
        nonlocal check_count
        error = abs(lhs - rhs)
        worst[name] = max(worst.get(name, 0.0), error)
        check_count += 1
        if error > tol * max(1.0, abs(lhs), abs(rhs)):
            raise AssertionError((name, lhs, rhs, error))

    def require(name, condition):
        nonlocal check_count
        check_count += 1
        if not condition:
            raise AssertionError(name)

    rows = []
    for c in CS:
        for t in TS:
            ui, uc, dui, duc = profile(c, t)
            u = tuple(ui if e in INTERNAL else uc for e in range(len(EDGES)))
            du = tuple(dui if e in INTERNAL else duc for e in range(len(EDGES)))
            table = {a: evaluate(a, c, t) for a in SIGNINGS}
            best, active, classes = minimizers(table)
            for a, row in table.items():
                check_close('group_cosh_factorization', row['z'], row['factor_z'])
                check_close('internal_group_derivative', row['ri'], row['factor_ri'])
                check_close('cross_group_derivative', row['rc'], row['factor_rc'])
                require('internal_group_nonnegative', row['ri'] >= -IDENTITY_TOL)
                require('cross_group_nonnegative', row['rc'] >= -IDENTITY_TOL)
                for e in range(len(EDGES)):
                    flipped = list(a)
                    flipped[e] *= -1
                    lhs = math.exp(table[tuple(flipped)]['f'] - row['f'])
                    rhs = math.cosh(2.0 * u[e]) - row['r'][e] * math.sinh(2.0 * u[e])
                    check_close('pure_edge_flip_ratio', lhs, rhs)
            for a in active:
                row = table[a]
                r = row['r']
                for e in range(len(EDGES)):
                    require('edge_local_optimality', r[e] <= math.tanh(u[e]) + IDENTITY_TOL)
                for group, ug in ((INTERNAL, ui), (CROSS, uc)):
                    require('group_l1_bound', math.fsum(abs(r[e]) for e in group)
                            <= 2.0 * len(group) * math.tanh(ug) + IDENTITY_TOL)
                h = tuple(math.tanh(v) - v / math.cosh(v) ** 2 for v in u)
                d = tuple(u[e] * (1.0 - r[e] ** 2) + h[e] - r[e]
                          for e in range(len(EDGES)))
                for e in range(len(EDGES)):
                    require('h_nonnegative_cubic_bound',
                            -IDENTITY_TOL <= h[e] <= (2.0 / 3.0) * u[e] ** 3 + IDENTITY_TOL)
                    require('defect_nonnegative', d[e] >= -IDENTITY_TOL)
                    check_close('scalar_defect_identity',
                                r[e], u[e] * (1.0 - r[e] ** 2) + h[e] - d[e])
                constant = math.fsum(u[e] * du[e] for e in range(len(EDGES)))
                square = -math.fsum(u[e] * du[e] * r[e] ** 2 for e in range(len(EDGES)))
                smooth = math.fsum(du[e] * h[e] for e in range(len(EDGES)))
                defect = math.fsum(du[e] * d[e] for e in range(len(EDGES)))
                check_close('variance_derivative_constant', constant, c * c / 4.0)
                check_close('frozen_derivative_defect_decomposition',
                            row['derivative'], constant + square + smooth - defect)
                require('uniform_signed_square_bound', abs(square) <= c ** 3 * math.sqrt(N) / 2.0 + IDENTITY_TOL)
                require('h_derivative_bound', abs(smooth) <= c ** 4 / 6.0 + IDENTITY_TOL)
            left = {a: evaluate(a, c, t - FD_STEP) for a in SIGNINGS}
            right = {a: evaluate(a, c, t + FD_STEP) for a in SIGNINGS}
            fl, _, cl = minimizers(left)
            fr, _, cr = minimizers(right)
            fd = dict(status='tied_or_changed_energy_class_skip',
                      central_classes=len(classes), left_classes=len(cl), right_classes=len(cr))
            if len(classes) == 1 and cl == classes and cr == classes:
                numerical = (fr - fl) / (2.0 * FD_STEP)
                analytic = table[active[0]]['derivative']
                check_close('isolated_class_minimum_finite_difference', numerical, analytic, tol=2e-7)
                fd = dict(status='one_joint_energy_class_checked', step=FD_STEP,
                          numerical=numerical, analytic=analytic, absolute_error=abs(numerical - analytic))
            rows.append(dict(c=c, t=t, minimum_pressure=best,
                             minimizing_signings=len(active), minimizing_joint_energy_classes=len(classes),
                             selected_signing=active[0], finite_difference=fd))
    result = dict(status='PASS', classification='finite_formula_regression_not_a_convergence_certificate',
                  n=N, signings_per_profile=len(SIGNINGS), states_per_signing=len(STATES),
                  prescribed_profiles=len(CS) * len(TS), checks=check_count,
                  worst_absolute_identity_errors=worst, profiles=rows,
                  inputs=dict(proof_sha256=proof_sha,
                              script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  environment=dict(hostname=socket.gethostname(), python=platform.python_version(),
                                   implementation=platform.python_implementation(), workers=1))
    canonical = json.dumps(result, sort_keys=True, separators=(',', ':')).encode()
    result['canonical_result_sha256_before_this_field'] = hashlib.sha256(canonical).hexdigest()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
