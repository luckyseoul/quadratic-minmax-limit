#!/usr/bin/env python3
"""N=6 optimized balanced-path probe and exact GPU joint-count catalog.

The finite grid does not certify an all-temperature envelope. Exact joint
counts are exported for an independent integer-polynomial dominance proof.
No host-norm classification, random cross blocks, or larger-order search.
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
EDGES = tuple(itertools.combinations(range(N), 2))
FREE = tuple(k for k, (i, j) in enumerate(EDGES) if i != 0)
INTERNAL = tuple(k for k, (i, j) in enumerate(EDGES) if i // 3 == j // 3)
CROSS = tuple(k for k in range(len(EDGES)) if k not in INTERNAL)
CS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
TS = tuple(k / 64.0 for k in range(65))
ATOL = 2e-11


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('result.json'))
    args = parser.parse_args()
    started = time.time()
    masks = np.arange(1024, dtype=np.int64)
    sign_rows = np.ones((1024, 15), dtype=np.int8)
    for bit, edge in enumerate(FREE):
        sign_rows[:, edge] = 1 - 2 * ((masks >> bit) & 1)
    states = np.asarray(list(itertools.product((-1, 1), repeat=N)), dtype=np.int8)
    characters = np.asarray([[int(x[i] * x[j]) for i, j in EDGES] for x in states], dtype=np.int8)
    a = cp.asarray(sign_rows, dtype=cp.float64)
    chi = cp.asarray(characters.T, dtype=cp.float64)
    internal = a[:, INTERNAL] @ chi[INTERNAL, :]
    cross = a[:, CROSS] @ chi[CROSS, :]
    for values in (internal, cross):
        assert float(cp.max(cp.abs(values - cp.rint(values))).get()) == 0.0
    abs_internal = cp.abs(internal).astype(cp.int64)
    abs_cross = cp.abs(cross).astype(cp.int64)
    assert bool(cp.all(abs_internal % 2 == 0).get())
    assert bool(cp.all(abs_cross % 2 == 1).get())
    bins = (abs_internal // 2) * 5 + (abs_cross - 1) // 2
    all_indices = cp.arange(1024, dtype=cp.int64)[:, None] * 20 + bins
    flat_counts = cp.bincount(all_indices.ravel(), minlength=1024 * 20)
    counts = cp.asnumpy(flat_counts.reshape(1024, 20)).astype(np.int64)
    assert np.all(counts.sum(axis=1) == 64)
    signatures, first_indices, inverse, multiplicities = np.unique(
        counts, axis=0, return_index=True, return_inverse=True, return_counts=True)
    catalog = [dict(signature_id=int(sid), representative_mask=int(first_indices[sid]),
                    mask_count=int(multiplicities[sid]), counts=row.tolist())
               for sid, row in enumerate(signatures)]
    rows = []
    values_archive = []
    derivatives_archive = []
    endpoint_summary = []
    for c in CS:
        c_rows = []
        for t in TS:
            ui = c * math.sqrt((2.0 - t) / N)
            uc = c * math.sqrt(t / N)
            energy = ui * internal + uc * cross
            absolute = cp.abs(energy)
            peak = cp.max(absolute, axis=1, keepdims=True)
            scaled = cp.exp(absolute - peak) / 2.0
            even = scaled * (1.0 + cp.exp(-2.0 * absolute))
            odd = scaled * cp.sign(energy) * (1.0 - cp.exp(-2.0 * absolute))
            denominator = cp.sum(even, axis=1)
            pressure = peak[:, 0] + cp.log(denominator / 64.0)
            ri = cp.sum(internal * odd, axis=1) / denominator
            dui = -c / (2.0 * math.sqrt(N) * math.sqrt(2.0 - t))
            if t == 0.0:
                # Right derivative of cosh(c sqrt(t/N) C) at zero.
                cross_second = cp.sum(cross * cross * even, axis=1) / denominator
                derivative = dui * ri + c * c * cross_second / (2.0 * N)
            else:
                rc = cp.sum(cross * odd, axis=1) / denominator
                duc = c / (2.0 * math.sqrt(N) * math.sqrt(t))
                derivative = dui * ri + duc * rc
            host_f = cp.asnumpy(pressure)
            host_d = cp.asnumpy(derivative)
            best = float(np.min(host_f))
            active = np.flatnonzero(host_f <= best + ATOL)
            exact_class_ids = np.unique(inverse[active]).tolist()
            selected = int(active[0])
            row = dict(c=c, t=t, internal_strength=ui, cross_strength=uc,
                       minimum_pressure=best, selected_mask=selected,
                       selected_signature_id=int(inverse[selected]),
                       active_mask_count=int(active.size),
                       active_signature_ids=exact_class_ids,
                       active_branch_derivative_min=float(np.min(host_d[active])),
                       active_branch_derivative_max=float(np.max(host_d[active])),
                       selected_branch_derivative=float(host_d[selected]),
                       derivative_scope=('right branch derivatives at t=0' if t == 0.0
                                         else 'ordinary frozen-branch derivative'))
            rows.append(row)
            c_rows.append(row)
            values_archive.append(host_f)
            derivatives_archive.append(host_d)
        # Every order-three signing has the same absolute-energy law.
        beta3 = c / math.sqrt(3.0)
        p3 = math.log((math.cosh(3.0 * beta3) + 3.0 * math.cosh(beta3)) / 4.0)
        endpoint_summary.append(dict(c=c, f0=c_rows[0]['minimum_pressure'],
                                     f1=c_rows[-1]['minimum_pressure'],
                                     endpoint_difference=c_rows[-1]['minimum_pressure'] - c_rows[0]['minimum_pressure'],
                                     twice_order_three_symmetric_pressure=2.0 * p3))
        print(json.dumps(dict(progress_c=c, profiles=len(c_rows),
                              signature_count=len(catalog))), flush=True)
    archive = args.output.with_suffix('.npz')
    np.savez_compressed(archive, masks=masks, sign_rows=sign_rows,
                        joint_counts=counts, signatures=signatures,
                        mask_to_signature_id=inverse,
                        pressures=np.asarray(values_archive),
                        derivatives=np.asarray(derivatives_archive),
                        c_values=np.asarray(CS), t_values=np.asarray(TS))
    cp.cuda.Stream.null.synchronize()
    result = dict(schema='original-mo-optimized-profile-gpu-v1', status='PASS',
                  classification='exact integer joint counts plus float64 finite grid; not an all-temperature or all-orders certificate',
                  n=N, block_sizes=[3, 3], gauge_rule='first row positive; free edges lexicographic, bit 1 means negative',
                  signing_count=1024, full_spin_state_count=64,
                  edges=[list(e) for e in EDGES], free_edge_indices=list(FREE),
                  internal_edge_indices=list(INTERNAL), cross_edge_indices=list(CROSS),
                  count_bin_order=[[i, j] for i in (0, 2, 4, 6) for j in (1, 3, 5, 7, 9)],
                  mask_to_signature_id=inverse.tolist(), signature_catalog=catalog,
                  c_values=list(CS), t_values=list(TS), profiles=rows,
                  endpoint_summary=endpoint_summary, minimizing_tolerance=ATOL,
                  sample_archive_name=archive.name,
                  sample_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                  execution=dict(host=socket.gethostname(), pid=__import__('os').getpid(),
                                 python=platform.python_version(), numpy=np.__version__, cupy=cp.__version__,
                                 device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                 started_unix_seconds=started, elapsed_seconds=time.time() - started))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(status='PASS', output=str(args.output),
                          result_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                          elapsed_seconds=result['execution']['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
