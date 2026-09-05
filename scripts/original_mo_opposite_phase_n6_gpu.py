#!/usr/bin/env python3
"""Fixed n=6 CUDA probe of actual opposite-phase Gibbs bilinear moments.

Enumerates all 1024 switching-normalized signings and 32 antipodal spin
representatives. Float64 pressure comparisons are numerical, not certified
optimizer classifications. T is evaluated directly as E[(x^T A y)^2],
independently of the CPU covariance/trace implementation. No larger census.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import socket
import time

import cupy as cp

N = 6
CS = (0.5, 1.0, 2.0, 4.0, 8.0)
TOL = 2e-11


def main():
    start = time.monotonic()
    edges = list(itertools.combinations(range(N), 2))
    free = [e for e, (i, j) in enumerate(edges) if i != 0]
    signs = []
    for mask in range(1 << len(free)):
        a = [1] * len(edges)
        for bit, e in enumerate(free):
            a[e] = 1 if not ((mask >> bit) & 1) else -1
        signs.append(a)
    states = [(1,) + tail for tail in itertools.product((-1, 1), repeat=N - 1)]
    chars = [[x[i] * x[j] for x in states] for i, j in edges]
    pairchars = [[x[i] * y[j] + x[j] * y[i]
                  for x in states for y in states] for i, j in edges]
    a_gpu = cp.asarray(signs, dtype=cp.float64)
    q = a_gpu @ cp.asarray(chars, dtype=cp.float64)
    b = (a_gpu @ cp.asarray(pairchars, dtype=cp.float64)).reshape(-1, len(states), len(states))
    b2 = b * b
    # Both integer expressions are exactly representable in float64.
    assert float(cp.max(cp.abs(q - cp.rint(q))).get()) == 0.0
    assert float(cp.max(cp.abs(b - cp.rint(b))).get()) == 0.0
    q_cpu = cp.asnumpy(q).astype(int)
    norm = cp.max(cp.abs(q), axis=1)
    rows = []
    for c in CS:
        beta = c / math.sqrt(N)
        hp = beta * q
        hm = -hp
        maxp = cp.max(hp, axis=1, keepdims=True)
        maxm = cp.max(hm, axis=1, keepdims=True)
        ep = cp.exp(hp - maxp)
        em = cp.exp(hm - maxm)
        sump = cp.sum(ep, axis=1, keepdims=True)
        summ = cp.sum(em, axis=1, keepdims=True)
        wp, wm = ep / sump, em / summ
        logzp = maxp[:, 0] + cp.log(sump[:, 0] / len(states))
        logzm = maxm[:, 0] + cp.log(summ[:, 0] / len(states))
        aval = (logzp + logzm) / 2
        aprime = cp.sum((wp - wm) * q, axis=1) / 2
        moment = cp.sum(wp[:, :, None] * b2 * wm[:, None, :], axis=(1, 2))
        values = cp.asnumpy(cp.stack((aval, aprime, moment, norm), axis=1))
        best = float(cp.min(aval).get())
        active = [i for i, v in enumerate(values) if v[0] <= best + TOL]
        others = [float(v[0]) for v in values if v[0] > best + TOL]
        types = {}
        for i in active:
            spectrum = sorted(q_cpu[i].tolist())
            reverse = sorted((-q_cpu[i]).tolist())
            key = tuple(min(spectrum, reverse))
            # Do not silently identify energy-spectrum types with isomorphism classes.
            if key not in types:
                types[key] = dict(representative_mask=i, candidate_masks=[],
                                  energy_spectrum=list(key),
                                  a=float(values[i, 0]), aprime=float(values[i, 1]),
                                  T=float(values[i, 2]), Phi=int(values[i, 3]),
                                  T_over_n2=float(values[i, 2]) / N**2,
                                  virial_gap_over_n=(beta * float(values[i, 2]) / 2 - float(values[i, 1])) / N)
            types[key]['candidate_masks'].append(i)
        active_t = [float(values[i, 2]) for i in active]
        rows.append(dict(c=c, beta=beta, minimum_a=best,
                         numerical_candidate_count=len(active),
                         next_pressure_gap=min(others) - best if others else None,
                         candidate_T_range=[min(active_t), max(active_t)],
                         candidate_spectral_types=list(types.values()),
                         maximum_T_over_n2=float(cp.max(moment).get()) / N**2))
    cp.cuda.Stream.null.synchronize()
    result = dict(classification='bounded_numerical_actual_host_probe_not_convergence_proof',
                  status='completed', n=N, switching_normalized_signings=len(signs),
                  antipodal_spin_representatives=len(states), c_values=list(CS),
                  pressure_tolerance=TOL, method='direct_pair_bilinear_squared_float64_cuda',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  hostname=socket.gethostname(), python=platform.python_version(),
                  cupy_version=cp.__version__,
                  cuda_device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                  elapsed_seconds=time.monotonic() - start, profiles=rows)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
