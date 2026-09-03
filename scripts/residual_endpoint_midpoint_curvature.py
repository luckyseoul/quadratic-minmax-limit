#!/usr/bin/env python3
"""Endpoint--midpoint curvature seam for residual (ii).

For the normalized off-diagonal coefficient matrix W_L put

    E_L(a) = sum_{t != a} W_L(a,t),
    A_L(a) = sum_{s<t, s+t=2a} W_L(s,t).

The parallel-edge diagonal cancels from

    H_L = eta_L (E_L - 2 A_L) = Radon(kappa)_L,
    kappa = d_z - 2 g.

Thus p*kappa(x)=sum_L H_L(Lx), kappa mod 2 is the graph boundary, and
sum H_L(a)^2=p*sum kappa(x)^2 >= p*|boundary|.  ``coeff-check`` verifies
the cancellation on every one-edge line-label basis vector.  ``sample``
uses a CUDA/ROCm array backend to sample compact all-prime atom labels and
measure the exact integer curvature statistics.  Sampling is discovery,
not exhaustive evidence and not residual closure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np

from residual_midpoint_seam_gpu import balanced, directions_and_signs, is_prime


def coeff_check(p: int) -> dict[str, object]:
    """Exhaust every line-label pair and both signs coefficientwise."""
    if not is_prime(p) or p < 5:
        raise ValueError("p must be an odd prime")
    inv2 = pow(2, -1, p)
    errors = 0
    cases = 0
    for eta in (-1, 1):
        for tau in (-1, 1):
            for s in range(p):
                for t in range(s, p):
                    endpoint = np.zeros(p, dtype=np.int64)
                    anti = np.zeros(p, dtype=np.int64)
                    direct = np.zeros(p, dtype=np.int64)
                    if s != t:
                        w = eta * tau
                        endpoint[s] += w
                        endpoint[t] += w
                        midpoint = ((s + t) * inv2) % p
                        anti[midpoint] += w
                        direct[s] += tau
                        direct[t] += tau
                        direct[midpoint] -= 2 * tau
                    # For s=t the endpoint, midpoint, and direct curvature
                    # all vanish: this is precisely diagonal cancellation.
                    if not np.array_equal(eta * (endpoint - 2 * anti), direct):
                        errors += 1
                    cases += 1
    return {
        "command": "coeff-check",
        "classification": "exhaustive coefficient-basis validation of an algebraic identity",
        "host": platform.node(),
        "architecture": platform.machine(),
        "p": p,
        "eta_tau_line_label_cases": cases,
        "nonzero_errors": errors,
        "diagonal_cancellation_validated": errors == 0,
        "formula": "eta_L*(E_L-2*A_L)=Radon(d_z-2*g)_L",
    }


def distinct_pairs(rng: np.random.Generator, batch: int, p: int) -> np.ndarray:
    out = rng.integers(0, p, size=(batch, 2), dtype=np.int32)
    bad = out[:, 0] == out[:, 1]
    while np.any(bad):
        out[bad, 1] = rng.integers(0, p, size=int(np.count_nonzero(bad)), dtype=np.int32)
        bad = out[:, 0] == out[:, 1]
    return out


def distinct_triples(rng: np.random.Generator, batch: int, count: int, p: int) -> np.ndarray:
    if count == 0:
        return np.empty((batch, 0, 3), dtype=np.int32)
    out = rng.integers(0, p, size=(batch, count, 3), dtype=np.int32)
    bad = (out[:, :, 0] == out[:, :, 1]) | (out[:, :, 0] == out[:, :, 2]) | (out[:, :, 1] == out[:, :, 2])
    while np.any(bad):
        out[bad] = rng.integers(0, p, size=(int(np.count_nonzero(bad)), 3), dtype=np.int32)
        bad = (out[:, :, 0] == out[:, :, 1]) | (out[:, :, 0] == out[:, :, 2]) | (out[:, :, 1] == out[:, :, 2])
    return out


def add_edge_curvature(cp, array, direction: int, u, v, coefficient: int, p: int) -> None:
    batch_index = cp.arange(array.shape[0], dtype=cp.int32)
    midpoint = ((u + v) * pow(2, -1, p)) % p
    cp.add.at(array, (batch_index, direction, u), coefficient)
    cp.add.at(array, (batch_index, direction, v), coefficient)
    cp.add.at(array, (batch_index, direction, midpoint), -2 * coefficient)


def add_triangles(cp, array, direction: int, triples: np.ndarray, kind: str, p: int) -> None:
    for atom in range(triples.shape[1]):
        values = cp.asarray(triples[:, atom, :])
        a, b, c = values[:, 0], values[:, 1], values[:, 2]
        add_edge_curvature(cp, array, direction, a, b, 1, p)
        add_edge_curvature(cp, array, direction, a, c, 1 if kind == "plus" else -1, p)
        add_edge_curvature(cp, array, direction, b, c, 1 if kind == "plus" else -1, p)


def branch_parameters(p: int, branch: str, endpoint: str) -> dict[str, object]:
    m = (p + 1) // 2
    if branch == "B":
        if p % 4 != 1:
            raise ValueError("branch B requires p=1 mod 4")
        r = (p - 1) // 4
        lower, upper = 2 * r * r - 5 * r, 4 * r * r - 6 * r - 3
        t = lower if endpoint == "lower" else upper
        return {
            "r": r, "t": t, "hard_star_sign": 1,
            "hard_minus": balanced(t, m),
            "opposite_pairs": True,
            "opposite_plus": [r - 2] * m,
            "opposite_minus": [q - r for q in balanced(6 * r + t, m, r)],
        }
    if branch == "C":
        if p % 4 != 3:
            raise ValueError("branch C requires p=3 mod 4")
        r = (p - 3) // 4
        lower, upper = 2 * r * r - 4 * r - 2, 4 * r * r - 2 * r - 5
        t = lower if endpoint == "lower" else upper
        return {
            "r": r, "t": t, "hard_star_sign": -1,
            "hard_minus": balanced(t + 1, m),
            "opposite_pairs": False,
            "opposite_plus": [r - 1] * m,
            "opposite_minus": [q - r - 2 for q in balanced(10 * r + 6 + t, m, r + 2)],
        }
    raise ValueError("branch must be B or C")


def candidate_record(config: dict[str, object], index: int, K: np.ndarray, energy: int,
                     boundary_size: int, labels: np.ndarray, eta: np.ndarray, p: int) -> dict[str, object]:
    picked = {
        "hard_centres": config["hard_centres"][index].tolist(),
        "hard_minus_triples": [value[index].tolist() for value in config["hard_minus"]],
        "opposite_pairs": None if config["opposite_pairs"] is None else config["opposite_pairs"][index].tolist(),
        "opposite_plus_triples": [value[index].tolist() for value in config["opposite_plus"]],
        "opposite_minus_triples": [value[index].tolist() for value in config["opposite_minus"]],
    }
    encoded = json.dumps(picked, sort_keys=True, separators=(",", ":")).encode()
    numerator = np.sum((eta[:, None] * K)[np.arange(p + 1)[:, None], labels], axis=0)
    remainder = numerator % p
    divisible = not np.any(remainder)
    record: dict[str, object] = {
        "energy": energy,
        "boundary_size_from_hard_fibre_parity": boundary_size,
        "energy_minus_p_boundary": energy - p * boundary_size,
        "configuration_sha256": hashlib.sha256(encoded).hexdigest(),
        "curvature_rows_sha256": hashlib.sha256(K.astype("<i4").tobytes()).hexdigest(),
        "nonzero_radon_inversion_remainders_mod_p": int(np.count_nonzero(remainder)),
        "radon_image_divisibility": divisible,
        "configuration": picked,
    }
    if divisible:
        kappa = numerator // p
        boundary = (kappa % 2) != 0
        record.update({
            "kappa_sha256": hashlib.sha256(kappa.astype("<i8").tobytes()).hexdigest(),
            "kappa_total": int(kappa.sum()),
            "kappa_energy": int(np.dot(kappa, kappa)),
            "parseval_exact": energy == p * int(np.dot(kappa, kappa)),
            "kappa_odd_support": int(np.count_nonzero(boundary)),
            "kappa_parity_matches_hard_fibres": int(np.count_nonzero(boundary)) == boundary_size,
        })
    return record


def sample(args: argparse.Namespace) -> dict[str, object]:
    import cupy as cp

    p = args.p
    if not is_prime(p) or p < 29:
        raise ValueError("sample requires a held-out prime p>=29")
    directions, eta, nonsquare = directions_and_signs(p)
    hard = np.flatnonzero(eta == 1)
    opposite = np.flatnonzero(eta == -1)
    params = branch_parameters(p, args.branch, args.endpoint)
    hard_minus_counts = params["hard_minus"]
    opposite_plus_counts = params["opposite_plus"]
    opposite_minus_counts = params["opposite_minus"]
    rng = np.random.default_rng(args.seed)
    points = np.array([(x, y) for x in range(p) for y in range(p)], dtype=np.int32)
    labels = (directions[:, 0, None] * points[:, 0] + directions[:, 1, None] * points[:, 1]) % p

    extrema: dict[str, tuple[int, dict[str, object]] | None] = {
        "minimum_energy": None,
        "maximum_energy": None,
        "minimum_energy_minus_p_boundary": None,
        "maximum_energy_minus_p_boundary": None,
    }
    started = time.time()
    completed = 0
    for offset in range(0, args.samples, args.batch):
        batch = min(args.batch, args.samples - offset)
        K = cp.zeros((batch, p + 1, p), dtype=cp.int32)
        batch_index = cp.arange(batch, dtype=cp.int32)
        hard_centres = rng.integers(0, p, size=(batch, len(hard)), dtype=np.int32)
        hard_minus: list[np.ndarray] = []
        for slot, direction in enumerate(hard):
            sign = int(params["hard_star_sign"])
            K[:, int(direction), :] -= sign
            centre = cp.asarray(hard_centres[:, slot])
            cp.add.at(K, (batch_index, int(direction), centre), sign * p)
            triples = distinct_triples(rng, batch, int(hard_minus_counts[slot]), p)
            hard_minus.append(triples)
            add_triangles(cp, K, int(direction), triples, "minus", p)

        opposite_pairs = None
        if bool(params["opposite_pairs"]):
            opposite_pairs = np.empty((batch, len(opposite), 2), dtype=np.int32)
        opposite_plus: list[np.ndarray] = []
        opposite_minus: list[np.ndarray] = []
        for slot, direction in enumerate(opposite):
            if opposite_pairs is not None:
                pair = distinct_pairs(rng, batch, p)
                opposite_pairs[:, slot, :] = pair
                a, b = cp.asarray(pair[:, 0]), cp.asarray(pair[:, 1])
                K[:, int(direction), :] += 2
                cp.add.at(K, (batch_index, int(direction), a), -p)
                cp.add.at(K, (batch_index, int(direction), b), -p)
                add_edge_curvature(cp, K, int(direction), a, b, 1, p)
            plus = distinct_triples(rng, batch, int(opposite_plus_counts[slot]), p)
            minus = distinct_triples(rng, batch, int(opposite_minus_counts[slot]), p)
            opposite_plus.append(plus); opposite_minus.append(minus)
            add_triangles(cp, K, int(direction), plus, "plus", p)
            add_triangles(cp, K, int(direction), minus, "minus", p)

        energy = cp.sum(K.astype(cp.int64) ** 2, axis=(1, 2))
        hard_labels = cp.asarray(labels[hard])
        centres_gpu = cp.asarray(hard_centres)
        matches = cp.sum(hard_labels[None, :, :] == centres_gpu[:, :, None], axis=1)
        # Keep this tiny reduction on the host.  CuPy 14's CUB reduction
        # currently pulls CUDA-13 fp4/fp6/fp8 headers that NVRTC cannot parse
        # for the V100 (sm_70); the large coefficient construction and energy
        # reduction remain on-device.
        matches_host = cp.asnumpy(matches)
        boundary_host = np.sum(
            ((len(hard) - matches_host) & 1) != 0, axis=1, dtype=np.int64
        )
        energy_host = cp.asnumpy(energy)
        defect_host = energy_host - p * boundary_host
        config = {
            "hard_centres": hard_centres,
            "hard_minus": hard_minus,
            "opposite_pairs": opposite_pairs,
            "opposite_plus": opposite_plus,
            "opposite_minus": opposite_minus,
        }
        K_host_cache: dict[int, np.ndarray] = {}
        choices = {
            "minimum_energy": int(np.argmin(energy_host)),
            "maximum_energy": int(np.argmax(energy_host)),
            "minimum_energy_minus_p_boundary": int(np.argmin(defect_host)),
            "maximum_energy_minus_p_boundary": int(np.argmax(defect_host)),
        }
        values = {
            "minimum_energy": int(energy_host[choices["minimum_energy"]]),
            "maximum_energy": int(energy_host[choices["maximum_energy"]]),
            "minimum_energy_minus_p_boundary": int(defect_host[choices["minimum_energy_minus_p_boundary"]]),
            "maximum_energy_minus_p_boundary": int(defect_host[choices["maximum_energy_minus_p_boundary"]]),
        }
        for key, local_index in choices.items():
            old = extrema[key]
            better = old is None or (values[key] < old[0] if key.startswith("minimum") else values[key] > old[0])
            if not better:
                continue
            if local_index not in K_host_cache:
                K_host_cache[local_index] = cp.asnumpy(K[local_index])
            record = candidate_record(
                config, local_index, K_host_cache[local_index], int(energy_host[local_index]),
                int(boundary_host[local_index]), labels, eta, p,
            )
            record["global_sample_index"] = offset + local_index
            extrema[key] = (values[key], record)
        completed += batch

    cp.cuda.runtime.deviceSynchronize()
    device = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())
    device_name = device["name"].decode() if isinstance(device["name"], bytes) else str(device["name"])
    return {
        "command": "sample",
        "classification": "non-exhaustive randomized compact-atom curvature search; not residual closure",
        "host": platform.node(), "architecture": platform.machine(),
        "gpu_backend": "cupy", "gpu_device": device_name,
        "p": p, "branch": args.branch, "endpoint": args.endpoint,
        "r": params["r"], "t": params["t"], "nonsquare": nonsquare,
        "seed": args.seed, "samples_requested": args.samples,
        "samples_completed": completed, "batch": args.batch,
        "elapsed_seconds": time.time() - started,
        "hard_directions": hard.tolist(), "opposite_directions": opposite.tolist(),
        "hard_minus_counts": hard_minus_counts,
        "opposite_plus_counts": opposite_plus_counts,
        "opposite_minus_counts": opposite_minus_counts,
        "extrema": {key: value[1] for key, value in extrema.items() if value is not None},
        "identity_under_test": "H=eta*(E-2*A)=Radon(kappa), sum(H^2)=p*sum(kappa^2)>=p*|boundary|",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("coeff-check")
    check.add_argument("--p", type=int, default=29)
    check.add_argument("--output", type=Path, required=True)
    probe = sub.add_parser("sample")
    probe.add_argument("--p", type=int, required=True)
    probe.add_argument("--branch", choices=("B", "C"), required=True)
    probe.add_argument("--endpoint", choices=("lower", "upper"), default="lower")
    probe.add_argument("--seed", type=int, required=True)
    probe.add_argument("--samples", type=int, default=32768)
    probe.add_argument("--batch", type=int, default=512)
    probe.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = coeff_check(args.p) if args.command == "coeff-check" else sample(args)
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
