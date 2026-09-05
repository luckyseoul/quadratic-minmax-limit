#!/usr/bin/env python3
"""Compare frozen CPU/GPU N=6 joint catalogs without rerunning enumeration.

This is an exact integer comparison of existing stored outputs. It does
not replay floating-point optimizer selection or certify any all-order claim.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import socket


def read_hashed(path, expected):
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise AssertionError(('input_hash_mismatch', str(path), actual, expected))
    return json.loads(raw), actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu-result', required=True, type=Path)
    parser.add_argument('--expected-cpu-sha256', required=True)
    parser.add_argument('--gpu-result', required=True, type=Path)
    parser.add_argument('--expected-gpu-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cpu, cpu_sha = read_hashed(args.cpu_result, args.expected_cpu_sha256)
    gpu, gpu_sha = read_hashed(args.gpu_result, args.expected_gpu_sha256)
    checks = 0

    def require(label, condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(label)

    bins = [[i, j] for i in (0, 2, 4, 6) for j in (1, 3, 5, 7, 9)]
    require('cpu_pass', cpu['status'] == 'PASS')
    require('fixed_cpu_scope', cpu['n'] == 6 and cpu['signings'] == 1024
            and cpu['states_per_signing'] == 64 and cpu['blocks'] == [3, 3])
    require('gpu_bin_order', gpu['count_bin_order'] == bins)
    require('gpu_mask_count', len(gpu['mask_to_signature_id']) == 1024)
    cpu_by_mask = {}
    for row in cpu['catalog']:
        counts_by_bin = {(i, j): count for i, j, count in row['absolute_joint_signature']}
        require('cpu_unique_signature_bins', len(counts_by_bin) == len(row['absolute_joint_signature']))
        counts = [counts_by_bin.get(tuple(key), 0) for key in bins]
        require('cpu_class_state_count', sum(counts) == 64)
        require('cpu_class_mask_count', len(row['masks']) == row['mask_count'])
        for mask in row['masks']:
            require('cpu_mask_not_repeated', mask not in cpu_by_mask)
            cpu_by_mask[mask] = counts
    require('cpu_complete_mask_coverage', sorted(cpu_by_mask) == list(range(1024)))
    gpu_by_id = {}
    for row in gpu['signature_catalog']:
        signature_id = row['signature_id']
        require('gpu_signature_id_unique', signature_id not in gpu_by_id)
        require('gpu_count_vector_shape', len(row['counts']) == len(bins))
        require('gpu_integer_count_entries', all(type(value) is int and value >= 0 for value in row['counts']))
        require('gpu_class_state_count', sum(row['counts']) == 64)
        gpu_by_id[signature_id] = row
    require('same_catalog_size', len(gpu_by_id) == len(cpu['catalog']))
    gpu_multiplicity = {signature_id: 0 for signature_id in gpu_by_id}
    for mask, signature_id in enumerate(gpu['mask_to_signature_id']):
        require('gpu_mask_references_known_class', signature_id in gpu_by_id)
        require('direct_cpu_equals_gpu_histogram', cpu_by_mask[mask] == gpu_by_id[signature_id]['counts'])
        gpu_multiplicity[signature_id] += 1
    for signature_id, row in gpu_by_id.items():
        require('gpu_class_mask_count', gpu_multiplicity[signature_id] == row['mask_count'])
        require('gpu_representative_mask_matches',
                gpu['mask_to_signature_id'][row['representative_mask']] == signature_id)
    canonical_catalog = json.dumps([cpu_by_mask[mask] for mask in range(1024)], separators=(',', ':')).encode()
    result = dict(status='PASS', classification='exact_stored_catalog_comparison_no_enumeration_or_float_replay',
                  signings=1024, states_per_signing=64, histogram_bins_per_signing=len(bins),
                  exact_histogram_entries_compared=1024 * len(bins),
                  signature_classes=len(gpu_by_id), checks=checks,
                  canonical_mask_ordered_histograms_sha256=hashlib.sha256(canonical_catalog).hexdigest(),
                  inputs=dict(cpu_result_sha256=cpu_sha, gpu_result_sha256=gpu_sha,
                              script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  environment=dict(hostname=socket.gethostname(), python=platform.python_version(), workers=1))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
