#!/usr/bin/env python3
"""Bounded new finite-step cross-block experiment, not an optimizer census.

One fixed order-six conference host (stored gauge mask 220), with fresh
Gaussian-sign and independent cross blocks. Exact spin sums in float64
compare individual, quenched, and sample-annealed pressures. Sample minima
are not exhaustive, and sample log means are not the true annealed value.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import socket
import time

import cupy as cp
import numpy as np

N = 6
CS = (0.5, 1.0, 2.0, 4.0, 8.0)
TS = (0.125, 0.25, 0.5, 1.0)


def logmeanexp(v, axis=-1):
    peak = cp.max(v, axis=axis, keepdims=True)
    return cp.squeeze(peak, axis=axis) + cp.log(cp.mean(cp.exp(v - peak), axis=axis))


def logcosh(v):
    a = cp.abs(v)
    return a + cp.log1p(cp.exp(-2 * a)) - math.log(2)


def host_matrix():
    a = [[0 if i == j else 1 for j in range(N)] for i in range(N)]
    for bit, (i, j) in enumerate(itertools.combinations(range(1, N), 2)):
        a[i][j] = a[j][i] = -1 if (220 >> bit) & 1 else 1
    return a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', type=int, default=8192)
    parser.add_argument('--batch', type=int, default=512)
    parser.add_argument('--seed', type=int, default=20260905)
    parser.add_argument('--out', type=Path, default=Path('result.json'))
    args = parser.parse_args()
    if args.samples < 1 or args.batch < 1:
        raise ValueError('positive sample and batch counts required')
    start = time.monotonic()
    a_list = host_matrix()
    a = cp.asarray(a_list, dtype=cp.float64)
    assert float(cp.max(cp.abs(a @ a - (N - 1) * cp.eye(N))).get()) == 0
    states = cp.asarray(list(itertools.product((-1, 1), repeat=N)), dtype=cp.float64)
    q = cp.einsum('xi,ij,xj->x', states, a, states) / 2
    assert float(cp.max(cp.abs(q - cp.rint(q))).get()) == 0
    # Row-major B vectorization; each column is the corresponding x_i y_j.
    chars = cp.einsum('xi,yj->xyij', states, states).reshape(64 * 64, N * N).T
    difference = (q[:, None] - q[None, :]).reshape(-1)
    rng_g = cp.random.RandomState(args.seed)
    rng_i = cp.random.RandomState(args.seed + 1)
    normal = rng_g.standard_normal((args.samples, N, N))
    # Cov(vec Z)=I+(A tensor A)/(N-1); alpha=0 for this conference host.
    z = (normal + a[None, :, :] @ normal @ a[None, :, :] / (N - 1)) / math.sqrt(2)
    gaussian = cp.where(z >= 0, 1, -1).astype(cp.float64)
    independent = (2 * rng_i.randint(0, 2, size=(args.samples, N, N)) - 1).astype(cp.float64)
    zero_gaussian_coordinates = int(cp.sum(z == 0).get())
    baseline = cp.stack((a + cp.eye(N), a - cp.eye(N)))
    candidates = {'gaussian': gaussian, 'independent': independent, 'coherent': baseline}
    profiles = []
    for c in CS:
        beta = c / math.sqrt(N)
        endpoint = float((logmeanexp(beta * q) + logmeanexp(-beta * q)).get())
        for t in TS:
            eta = beta * math.sqrt(1 - t / 2)
            gamma = beta * math.sqrt(t / 2)
            internal = eta * difference
            logweights = logcosh(internal)
            weights = cp.exp(logweights - cp.max(logweights))
            weights /= cp.sum(weights)
            # Both phases' actual covariance matrices have I +/- tau A form.
            ep = cp.exp(eta * q - cp.max(eta * q)); ep /= cp.sum(ep)
            em = cp.exp(-eta * q - cp.max(-eta * q)); em /= cp.sum(em)
            u = (states.T * ep) @ states
            v = (states.T * em) @ states
            tau = float(cp.trace(a @ u).get()) / (N * (N - 1))
            assert float(cp.max(cp.abs(u - cp.eye(N) - tau * a)).get()) < 3e-12
            assert float(cp.max(cp.abs(v - cp.eye(N) + tau * a)).get()) < 3e-12
            theory_mean_q = N * N - (2 / math.pi) * math.asin(1 / (N - 1)) * (tau * N * (N - 1)) ** 2
            profiles.append(dict(c=c, t=t, beta=beta, eta=eta, gamma=gamma,
                                 endpoint=endpoint, theory_gaussian_mean_q=theory_mean_q,
                                 _internal=internal, _weights=weights))
    replay_cases = []
    results = []
    for label, blocks in candidates.items():
        count = int(blocks.shape[0])
        f_all = cp.empty((count, len(profiles)), dtype=cp.float64)
        q_all = cp.empty_like(f_all)
        for first in range(0, count, args.batch):
            cross = blocks[first:first + args.batch].reshape(-1, N * N) @ chars
            assert float(cp.max(cp.abs(cross - cp.rint(cross))).get()) == 0
            for k, profile in enumerate(profiles):
                f_all[first:first + args.batch, k] = logmeanexp(
                    logcosh(profile['_internal'][None, :] + profile['gamma'] * cross), axis=1)
                q_all[first:first + args.batch, k] = (cross * cross) @ profile['_weights']
        f_cpu = cp.asnumpy(f_all)
        q_cpu = cp.asnumpy(q_all)
        blocks_cpu = cp.asnumpy(blocks).astype(np.int8)
        for k, profile in enumerate(profiles):
            best = int(np.argmin(f_cpu[:, k]))
            peak = float(np.max(f_cpu[:, k]))
            sample_annealed = peak + math.log(float(np.mean(np.exp(f_cpu[:, k] - peak))))
            row = {key: val for key, val in profile.items() if not key.startswith('_')}
            row.update(distribution=label, sample_count=count, minimum_F=float(f_cpu[best, k]),
                       maximum_F=peak, mean_F=float(np.mean(f_cpu[:, k])),
                       sample_log_mean_exp_F=sample_annealed,
                       minimum_gap=float(f_cpu[best, k]) - profile['endpoint'],
                       best_index=best, mean_q=float(np.mean(q_cpu[:, k])),
                       standard_error_mean_q=(float(np.std(q_cpu[:, k], ddof=1)) / math.sqrt(count)
                                              if count > 1 and label != "coherent" else None))
            results.append(row)
            for index in sorted(set((0, min(1, count - 1), best))):
                replay_cases.append(dict(case_id=f'{label}-c{profile["c"]}-t{profile["t"]}-i{index}',
                                         n=N, c=profile['c'], t=profile['t'],
                                         A=a_list, B=blocks_cpu[index].tolist(),
                                         F=float(f_cpu[index, k]), qbar=float(q_cpu[index, k]),
                                         endpoint=profile['endpoint'], candidate_label=label))
        print(json.dumps(dict(progress=label, blocks=count, profiles=len(profiles))), flush=True)
    archive = args.out.with_suffix('.npz')
    np.savez_compressed(archive, A=np.asarray(a_list, dtype=np.int8),
                        gaussian=cp.asnumpy(gaussian).astype(np.int8),
                        independent=cp.asnumpy(independent).astype(np.int8),
                        coherent=cp.asnumpy(baseline).astype(np.int8))
    cp.cuda.Stream.null.synchronize()
    result = dict(classification='bounded_new_finite_step_sample_not_optimizer_or_annealed_certificate',
                  status='completed', n=N, source_host_gauge_mask=220,
                  source_host_scope='fixed_conference_example_reused_not_new_minimizer_classification',
                  samples_per_random_distribution=args.samples, batch=args.batch, seed=args.seed,
                  spin_pairs=64 * 64, c_values=list(CS), t_values=list(TS),
                  zero_gaussian_coordinates=zero_gaussian_coordinates,
                  distribution_note='canonical Gaussian law has alpha=0; independent comparator; coherent A +/- I',
                  statistical_scope='sample minima are upper bounds only; sample log means can miss rare annealed events',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sample_archive_name=archive.name,
                  sample_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                  hostname=socket.gethostname(), python=platform.python_version(),
                  cupy_version=cp.__version__, numpy_version=np.__version__,
                  cuda_device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                  elapsed_seconds=time.monotonic() - start, profiles=results, replay_cases=replay_cases)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(status='completed', result=str(args.out),
                          result_sha256=hashlib.sha256(args.out.read_bytes()).hexdigest(),
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
