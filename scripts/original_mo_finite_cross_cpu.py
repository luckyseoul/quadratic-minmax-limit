#!/usr/bin/env python3
"""Independent bounded CPU regression for finite paired cross-block pressure.

Default scope: n=2, all 16 integral cross blocks, all 16 spin pairs, and
the four equally likely canonical Gaussian-sign blocks [[r,s],[s,r]].
Optional --gpu-result replays only the selected n=6 blocks in replay_cases;
it does not search, classify, or resample any order-six signing.

Exhaustive state sums use float64 stable log formulas, not interval arithmetic.
The output is finite numerical evidence, never an asymptotic certificate.
"""

import os

for _thread_setting in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
):
    os.environ[_thread_setting] = "1"

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import socket
import sys
import time

import numpy as np


CS = (0.5, 1.0, 2.0, 4.0, 8.0)
TS = (0.125, 0.25, 0.5, 1.0)
ATOL = 2e-11
LOG2 = math.log(2.0)


def logcosh(values):
    absolute = np.abs(np.asarray(values, dtype=np.float64))
    return absolute + np.log1p(np.exp(-2.0 * absolute)) - LOG2


def logsumexp(values):
    values = np.asarray(values, dtype=np.float64)
    largest = float(np.max(values))
    return largest + math.log(float(np.sum(np.exp(values - largest))))


def logmeanexp(values):
    values = np.asarray(values, dtype=np.float64)
    return logsumexp(values) - math.log(values.size)


def all_states(n):
    return np.asarray(list(itertools.product((-1.0, 1.0), repeat=n)))


def energies(A, states):
    # Direct unordered-edge sum, independent of a covariance-trace formula.
    answer = np.zeros(states.shape[0], dtype=np.float64)
    for i in range(A.shape[0]):
        for j in range(i + 1, A.shape[0]):
            answer += A[i, j] * states[:, i] * states[:, j]
    return answer


def endpoint(energy, beta):
    return logmeanexp(beta * energy) + logmeanexp(-beta * energy)


def branch_from_cross(energy, cross, beta, t):
    eta = beta * math.sqrt(1.0 - t / 2.0)
    gamma = beta * math.sqrt(t / 2.0)
    internal = eta * (energy[:, None] - energy[None, :])
    return logmeanexp(logcosh(internal + gamma * cross))


def close_record(checks, name, actual, expected):
    difference = abs(float(actual) - float(expected))
    checks.append({"name": name, "absolute_error": difference})
    if not math.isfinite(difference) or difference > ATOL:
        raise AssertionError(f"{name}: {actual!r} != {expected!r}; error={difference}")


def require_order(left, right, name):
    if not math.isfinite(left) or not math.isfinite(right) or left > right + ATOL:
        raise AssertionError(f"{name}: {left!r} exceeds {right!r}")


def cost_summary(costs, matrices):
    values = np.asarray(costs, dtype=np.float64)
    best = float(np.min(values))
    minimizers = [
        matrix.astype(int).tolist()
        for value, matrix in zip(values, matrices)
        if abs(float(value) - best) <= ATOL
    ]
    return {
        "count": len(matrices), "minimum_F": best,
        "mean_F": float(np.mean(values)), "log_mean_exp_F": logmeanexp(values),
        "maximum_F": float(np.max(values)), "minimizers_within_atol": minimizers,
    }


def exhaustive_n2():
    A = np.asarray([[0.0, 1.0], [1.0, 0.0]])
    states = all_states(2)
    energy = energies(A, states)
    blocks = [np.asarray(bits).reshape(2, 2) for bits in
              itertools.product((-1.0, 1.0), repeat=4)]
    gaussian_blocks = [np.asarray([[r, s], [s, r]]) for r, s in
                       itertools.product((-1.0, 1.0), repeat=2)]
    all_cross = [states @ B @ states.T for B in blocks]
    gaussian_cross = np.asarray([states @ B @ states.T for B in gaussian_blocks])
    # Column-vectorized Gaussian covariance: z11=z22 and z21=z12.
    sigma = np.eye(4) + np.kron(A, A)
    checks = []
    profiles = []
    for c in CS:
        beta = c / math.sqrt(2.0)
        initial = endpoint(energy, beta)
        close_record(checks, f"n2_endpoint_c{c}", initial, 2.0 * float(logcosh(beta)))
        for t in TS:
            label = f"n2_c{c}_t{t}"
            eta = beta * math.sqrt(1.0 - t / 2.0)
            gamma = beta * math.sqrt(t / 2.0)
            costs = [branch_from_cross(energy, cross, beta, t) for cross in all_cross]
            gaussian_costs = [branch_from_cross(energy, cross, beta, t)
                              for cross in gaussian_cross]
            all_summary = cost_summary(costs, blocks)
            gaussian_summary = cost_summary(gaussian_costs, gaussian_blocks)
            # Two independent closed formulas for the annealed partition sums.
            iid_closed = 2.0 * float(logcosh(eta)) + 4.0 * float(logcosh(gamma))
            gaussian_closed = logsumexp([
                float(logcosh(2.0 * eta)) - LOG2,
                float(logcosh(4.0 * gamma)) - 2.0 * LOG2,
                -2.0 * LOG2,
            ])
            close_record(checks, label + "_iid_annealed",
                         all_summary["log_mean_exp_F"], iid_closed)
            for block_index, value in enumerate(gaussian_costs):
                close_record(checks, label + f"_canonical_F_{block_index}",
                             value, gaussian_closed)
            for group_name, summary in (("all", all_summary), ("gaussian", gaussian_summary)):
                require_order(summary["minimum_F"], summary["mean_F"], group_name + " min/mean")
                require_order(summary["mean_F"], summary["log_mean_exp_F"], group_name + " Jensen")
            require_order(all_summary["minimum_F"], gaussian_summary["minimum_F"],
                          "all cross-block minimum versus Gaussian support")

            internal = eta * (energy[:, None] - energy[None, :])
            log_weights = logcosh(internal) - logsumexp(logcosh(internal))
            conditional = []
            log_mgfs = np.empty((4, 4))
            variances = np.empty((4, 4))
            for xi, x in enumerate(states):
                for yi, y in enumerate(states):
                    bilinear = gaussian_cross[:, xi, yi]
                    actual_log_mgf = logmeanexp(gamma * bilinear)
                    mean = float(np.mean(bilinear))
                    variance = float(np.mean((bilinear - mean) ** 2))
                    first = x[0] * y[0] + x[1] * y[1]
                    second = x[0] * y[1] + x[1] * y[0]
                    exact_log_mgf = float(logcosh(gamma * first) + logcosh(gamma * second))
                    vector = np.kron(y, x)
                    close_record(checks, label + f"_pair{xi}_{yi}_mean", mean, 0.0)
                    close_record(checks, label + f"_pair{xi}_{yi}_mgf",
                                 actual_log_mgf, exact_log_mgf)
                    close_record(checks, label + f"_pair{xi}_{yi}_variance",
                                 variance, float(vector @ sigma @ vector))
                    log_mgfs[xi, yi] = actual_log_mgf
                    variances[xi, yi] = variance
                    conditional.append({
                        "x": x.astype(int).tolist(), "y": y.astype(int).tolist(),
                        "internal_energy_difference": float(energy[xi] - energy[yi]),
                        "actual_internal_mixture_weight": float(math.exp(log_weights[xi, yi])),
                        "log_mgf": actual_log_mgf, "variance_bilinear": variance,
                        "half_gamma_squared_variance": 0.5 * gamma * gamma * variance,
                        "log_mgf_minus_half_variance": actual_log_mgf - 0.5 * gamma * gamma * variance,
                    })
            internal_endpoint = endpoint(energy, eta)
            annealed_from_conditionals = internal_endpoint + logsumexp(log_weights + log_mgfs)
            close_record(checks, label + "_conditional_annealed",
                         annealed_from_conditionals, gaussian_summary["log_mean_exp_F"])
            averaged_variance = float(np.sum(np.exp(log_weights) * variances))
            close_record(checks, label + "_actual_phase_variance",
                         averaged_variance, 4.0 / math.cosh(eta) ** 2)
            profiles.append({
                "n": 2, "c": c, "t": t, "beta": beta, "eta": eta, "gamma": gamma,
                "endpoint_2a_beta": initial, "internal_endpoint_2a_eta": internal_endpoint,
                "all_16_cross_blocks": all_summary,
                "canonical_gaussian_4_support": gaussian_summary,
                "all_16_F_values_in_block_order": costs,
                "canonical_gaussian_F_values_in_support_order": gaussian_costs,
                "minimum_all_gap_from_endpoint": all_summary["minimum_F"] - initial,
                "minimum_gaussian_gap_from_endpoint": gaussian_summary["minimum_F"] - initial,
                "gaussian_annealed_gap_from_endpoint": gaussian_summary["log_mean_exp_F"] - initial,
                "canonical_gaussian_conditional": conditional,
                "actual_phase_average_conditional_variance": averaged_variance,
                "actual_phase_mean_conditional_log_mgf": float(np.sum(np.exp(log_weights) * log_mgfs)),
                "gaussian_annealed_from_conditionals": annealed_from_conditionals,
                "gaussian_quadratic_mgf_proxy": internal_endpoint +
                    logsumexp(log_weights + 0.5 * gamma * gamma * variances),
            })
    return {
        "classification": "exhaustive n2 float64 finite regression, not an asymptotic certificate",
        "A": A.astype(int).tolist(), "spin_states": states.astype(int).tolist(),
        "spin_pair_count": 16, "all_cross_block_count": 16,
        "all_cross_blocks_in_order": [B.astype(int).tolist() for B in blocks],
        "canonical_gaussian_support_in_order": [B.astype(int).tolist() for B in gaussian_blocks],
        "canonical_gaussian_support_probabilities": [0.25] * 4,
        "canonical_gaussian_covariance_column_vectorized": sigma.tolist(),
        "profiles": profiles, "check_count": len(checks),
        "max_absolute_formula_error": max(item["absolute_error"] for item in checks),
        "formula_checks": checks,
    }


def validate_complete_signing(A, n, label):
    if A.shape != (n, n) or not np.array_equal(A, A.T):
        raise ValueError(label + " must be a symmetric n-by-n matrix")
    if not np.array_equal(np.diag(A), np.zeros(n)):
        raise ValueError(label + " must have zero diagonal")
    entries = A[~np.eye(n, dtype=bool)]
    if not np.all(np.isin(entries, (-1.0, 1.0))):
        raise ValueError(label + " must have off-diagonal signs")


def replay_gpu_winners(path):
    payload_bytes = path.read_bytes()
    payload = json.loads(payload_bytes)
    cases = payload.get("replay_cases")
    if not isinstance(cases, list) or not 1 <= len(cases) <= 256:
        raise ValueError("GPU JSON must supply 1..256 explicitly selected replay_cases")
    states = all_states(6)  # Full 64 states, no antipodal reduction.
    output = []
    checks = []
    for index, case in enumerate(cases):
        if int(case.get("n", 6)) != 6:
            raise ValueError("Only n=6 GPU winner replay is authorized")
        A = np.asarray(case["A"], dtype=np.float64)
        B = np.asarray(case["B"], dtype=np.float64)
        validate_complete_signing(A, 6, "replay A")
        if B.shape != (6, 6) or not np.all(np.isin(B, (-1.0, 1.0))):
            raise ValueError("replay B must have 36 freely chosen sign entries")
        c, t = float(case["c"]), float(case["t"])
        if not math.isfinite(c) or c <= 0.0 or not math.isfinite(t) or not 0.0 <= t <= 1.0:
            raise ValueError("Invalid finite replay temperature/profile")
        beta = c / math.sqrt(6.0)
        energy = energies(A, states)
        cross = states @ B @ states.T
        value = branch_from_cross(energy, cross, beta, t)
        initial = endpoint(energy, beta)
        case_id = str(case.get("case_id", index))
        close_record(checks, "n6_" + case_id + "_F", value, float(case["F"]))
        close_record(checks, "n6_" + case_id + "_endpoint_factorization",
                     branch_from_cross(energy, cross, beta, 0.0), initial)
        endpoint_key = "endpoint_2a_beta" if "endpoint_2a_beta" in case else "endpoint"
        if endpoint_key in case:
            close_record(checks, "n6_" + case_id + "_GPU_endpoint",
                         initial, float(case[endpoint_key]))
        output.append({
            "case_id": case_id, "candidate_label": case.get("candidate_label"),
            "n": 6, "c": c, "t": t, "A": A.astype(int).tolist(),
            "B": B.astype(int).tolist(), "F_cpu": value, "F_gpu": float(case["F"]),
            "absolute_F_difference": abs(value - float(case["F"])),
            "endpoint_2a_beta_cpu": initial, "gap_from_endpoint_cpu": value - initial,
            "spin_state_count": 64, "spin_pair_count": 4096,
        })
    return {
        "classification": "independent full-state replay of selected n6 blocks only",
        "input_json_sha256": hashlib.sha256(payload_bytes).hexdigest(),
        "input_json_name": path.name, "case_count": len(cases), "cases": output,
        "check_count": len(checks), "formula_checks": checks,
        "max_absolute_error": max(item["absolute_error"] for item in checks),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu-result", type=Path, help="JSON with selected n=6 replay_cases")
    parser.add_argument("--skip-n2", action="store_true", help="replay GPU winners without repeating n=2")
    parser.add_argument("--output", type=Path, help="write result JSON here, otherwise stdout")
    args = parser.parse_args()
    if args.skip_n2 and not args.gpu_result:
        parser.error("--skip-n2 requires --gpu-result")
    started = time.time()
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = {
        "schema_version": 1, "source_sha256": source_hash,
        "execution": {
            "host": socket.gethostname(), "pid": os.getpid(),
            "python": platform.python_version(), "numpy": np.__version__,
            "cpu_workers": 1, "blas_threads": 1,
            "started_unix_seconds": started,
            "argv": sys.argv,
        },
        "absolute_tolerance": ATOL, "n2": None if args.skip_n2 else exhaustive_n2(),
        "n6_gpu_winner_replay": replay_gpu_winners(args.gpu_result) if args.gpu_result else None,
    }
    result["execution"]["elapsed_seconds"] = time.time() - started
    result["status"] = "PASS"
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
        print(json.dumps({"status": "PASS", "output": str(args.output),
                          "result_sha256": hashlib.sha256(encoded.encode()).hexdigest(),
                          "n2_checks": (result["n2"] or {}).get("check_count", 0),
                          "n6_replay_cases": (result["n6_gpu_winner_replay"] or {}).get("case_count", 0),
                          "elapsed_seconds": result["execution"]["elapsed_seconds"]}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
