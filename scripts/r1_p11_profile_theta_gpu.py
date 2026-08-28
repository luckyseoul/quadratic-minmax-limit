#!/usr/bin/env python3
"""Exact p=11 ordinary dual-theta coefficients from profile glue.

The ten-dimensional finite-field Fourier filter is first quotiented by plane
translations and nonzero scalars.  The remaining weighted six-tuples are
split into two triples.  For each additive character, a sparse weighted
matrix contraction replaces a loop over all tuple-polynomial convolutions.
Several 31-bit primes give exact integer coefficients by CRT; the unrestricted
six-profile count is a rigorous reconstruction bound.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

from r1_p11_profile_dual_orbits import P


DEFAULT_MODULI = (
    1999999871,
    1999999321,
    1999999013,
    1999998947,
    1999997957,
)
KNOWN_COEFFICIENTS = {0: 1, 11: 244, 20: 16104, 24: 14762, 27: 442860}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primitive_eleventh_root(modulus: int) -> int:
    if (modulus - 1) % P:
        raise ValueError(f"modulus {modulus} is not 1 modulo 11")
    for candidate in range(2, 1000):
        root = pow(candidate, (modulus - 1) // P, modulus)
        if root != 1 and pow(root, P, modulus) == 1:
            return root
    raise ArithmeticError(f"failed to find an eleventh root modulo {modulus}")


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def unpack_triples(keys: np.ndarray) -> np.ndarray:
    return np.stack(
        [((keys >> (10 * index)) & 1023).astype(np.uint16) for index in range(3)],
        axis=1,
    )


def tuple_topology(keys: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    prefix_keys = keys & ((1 << 30) - 1)
    suffix_keys = keys >> 30
    unique_prefix, prefix_index = np.unique(prefix_keys, return_inverse=True)
    unique_suffix, suffix_index = np.unique(suffix_keys, return_inverse=True)
    return (
        unpack_triples(unique_prefix),
        unpack_triples(unique_suffix),
        prefix_index.astype(np.uint32),
        suffix_index.astype(np.uint32),
    )


def convolve_power(coefficients: list[int], power: int, max_k: int) -> list[int]:
    result = [1] + [0] * max_k
    for _ in range(power):
        result = [
            sum(result[index] * coefficients[degree - index] for index in range(degree + 1))
            for degree in range(max_k + 1)
        ]
    return result


def unrestricted_profile_products(counts: np.ndarray, max_k: int) -> list[list[int]]:
    products = []
    for residue_sum in range(P):
        coefficients = [
            sum(int(value) for value in counts[0, residue_sum, degree])
            for degree in range(max_k + 1)
        ]
        products.append(convolve_power(coefficients, 6, max_k))
    return products


def balanced_energy(common_sum: int) -> int:
    quotient, residue = divmod(common_sum, P)
    return P * quotient * quotient + 2 * quotient * residue + residue


def theta_bound_and_terms(
    unrestricted: list[list[int]], max_e: int, max_k: int
) -> tuple[
    list[int],
    list[int],
    list[int],
    list[list[tuple[int, int, int]]],
]:
    common_sum_bound = math.isqrt(P * max_e)
    bounds = [0] * (max_e + 1)
    second_moment_bounds = [0] * (max_e + 1)
    fourth_moment_bounds = [0] * (max_e + 1)
    terms: list[list[tuple[int, int, int]]] = [[] for _ in range(max_e + 1)]
    for exponent in range(max_e + 1):
        for common_sum in range(-common_sum_bound, common_sum_bound + 1):
            minimum = 12 * balanced_energy(common_sum) - common_sum * common_sum
            difference = exponent - minimum
            if difference < 0 or difference % 4:
                continue
            excess = difference // 4
            if excess > max_k:
                continue
            residue = common_sum % P
            terms[exponent].append((common_sum, residue, excess))
            count = unrestricted[residue][excess]
            bounds[exponent] += count
            second_moment_bounds[exponent] += common_sum**2 * count
            fourth_moment_bounds[exponent] += common_sum**4 * count
    return bounds, second_moment_bounds, fourth_moment_bounds, terms


CUDA_SOURCE = r"""
extern "C" {

__device__ __forceinline__ unsigned long long mul_mod(
    unsigned long long left,
    unsigned long long right,
    unsigned long long modulus
) {
    return (left * right) % modulus;
}

__global__ void build_triples(
    unsigned long long count,
    const unsigned short *ids,
    const unsigned long long *single,
    int stride,
    unsigned long long modulus,
    unsigned long long *output
) {
    unsigned long long row = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= count) return;
    unsigned long long current[PROFILE_STRIDE_CAPACITY];
    unsigned long long following[PROFILE_STRIDE_CAPACITY];
    unsigned int first = ids[3 * row];
    for (int degree = 0; degree < stride; ++degree)
        current[degree] = single[(unsigned long long)first * stride + degree];
    for (int factor = 1; factor < 3; ++factor) {
        unsigned int identifier = ids[3 * row + factor];
        for (int degree = 0; degree < stride; ++degree) {
            unsigned long long total = 0ull;
            for (int left_degree = 0; left_degree <= degree; ++left_degree) {
                unsigned long long term = mul_mod(
                    current[left_degree],
                    single[(unsigned long long)identifier * stride + degree - left_degree],
                    modulus
                );
                total += term;
                if (total >= modulus) total %= modulus;
            }
            following[degree] = total % modulus;
        }
        for (int degree = 0; degree < stride; ++degree)
            current[degree] = following[degree];
    }
    for (int degree = 0; degree < stride; ++degree)
        output[row * stride + degree] = current[degree];
}

__global__ void sparse_weighted_multiply(
    unsigned long long edge_count,
    const unsigned int *prefix_index,
    const unsigned int *suffix_index,
    const unsigned long long *weights,
    const unsigned long long *suffix_polynomials,
    int stride,
    unsigned long long modulus,
    unsigned long long *output
) {
    unsigned long long edge = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (edge >= edge_count) return;
    unsigned long long prefix = prefix_index[edge];
    unsigned long long suffix = suffix_index[edge];
    unsigned long long weight = weights[edge] % modulus;
    for (int degree = 0; degree < stride; ++degree) {
        unsigned long long contribution = mul_mod(
            weight, suffix_polynomials[suffix * stride + degree], modulus
        );
        atomicAdd(&output[prefix * stride + degree], contribution);
    }
}

__global__ void reduce_entries(
    unsigned long long count,
    unsigned long long modulus,
    unsigned long long *entries
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index < count) entries[index] %= modulus;
}

__global__ void bilinear_matrix(
    const unsigned long long *prefix_polynomials,
    const unsigned long long *weighted_suffix,
    unsigned long long prefix_count,
    int stride,
    unsigned long long modulus,
    unsigned long long *matrix
) {
    extern __shared__ unsigned long long partial[];
    int entry = blockIdx.x;
    int left_degree = entry / stride;
    int right_degree = entry % stride;
    unsigned long long total = 0ull;
    for (unsigned long long row = threadIdx.x; row < prefix_count; row += blockDim.x) {
        unsigned long long term = mul_mod(
            prefix_polynomials[row * stride + left_degree],
            weighted_suffix[row * stride + right_degree],
            modulus
        );
        total += term;
        if (total >= modulus) total %= modulus;
    }
    partial[threadIdx.x] = total % modulus;
    __syncthreads();
    for (int offset = blockDim.x / 2; offset; offset >>= 1) {
        if (threadIdx.x < offset) {
            partial[threadIdx.x] += partial[threadIdx.x + offset];
            if (partial[threadIdx.x] >= modulus) partial[threadIdx.x] %= modulus;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) matrix[entry] = partial[0] % modulus;
}

}
"""


class GPUContractor:
    def __init__(
        self,
        prefix_triples: np.ndarray,
        suffix_triples: np.ndarray,
        prefix_index: np.ndarray,
        suffix_index: np.ndarray,
        weights: np.ndarray,
        max_k: int,
        device: int,
    ) -> None:
        import cupy as cp

        cp.cuda.Device(device).use()
        self.cp = cp
        self.max_k = max_k
        self.stride = max_k + 1
        self.prefix_count = len(prefix_triples)
        self.suffix_count = len(suffix_triples)
        self.edge_count = len(weights)
        self.prefix_triples = cp.asarray(prefix_triples)
        self.suffix_triples = cp.asarray(suffix_triples)
        self.prefix_index = cp.asarray(prefix_index)
        self.suffix_index = cp.asarray(suffix_index)
        self.weights = cp.asarray(weights)
        if self.stride < 1:
            raise ValueError("profile polynomial stride must be positive")
        source = CUDA_SOURCE.replace("PROFILE_STRIDE_CAPACITY", str(self.stride))
        module = cp.RawModule(code=source, options=("--std=c++11",))
        self.build = module.get_function("build_triples")
        self.sparse = module.get_function("sparse_weighted_multiply")
        self.reduce = module.get_function("reduce_entries")
        self.bilinear = module.get_function("bilinear_matrix")
        self.block = 256

    def triple_polynomials(self, triples: object, count: int, single: object, modulus: int) -> object:
        cp = self.cp
        output = cp.empty((count, self.stride), dtype=cp.uint64)
        self.build(
            ((count + self.block - 1) // self.block,),
            (self.block,),
            (
                np.uint64(count), triples, single, np.int32(self.stride),
                np.uint64(modulus), output,
            ),
        )
        return output

    def contract(self, single: np.ndarray, modulus: int) -> np.ndarray:
        cp = self.cp
        d_single = cp.asarray(single, dtype=cp.uint64)
        prefix = self.triple_polynomials(
            self.prefix_triples, self.prefix_count, d_single, modulus
        )
        suffix = self.triple_polynomials(
            self.suffix_triples, self.suffix_count, d_single, modulus
        )
        weighted = cp.zeros((self.prefix_count, self.stride), dtype=cp.uint64)
        self.sparse(
            ((self.edge_count + self.block - 1) // self.block,),
            (self.block,),
            (
                np.uint64(self.edge_count), self.prefix_index, self.suffix_index,
                self.weights, suffix, np.int32(self.stride), np.uint64(modulus),
                weighted,
            ),
        )
        entry_count = self.prefix_count * self.stride
        self.reduce(
            ((entry_count + self.block - 1) // self.block,),
            (self.block,),
            (np.uint64(entry_count), np.uint64(modulus), weighted),
        )
        matrix = cp.empty((self.stride, self.stride), dtype=cp.uint64)
        self.bilinear(
            (self.stride * self.stride,),
            (self.block,),
            (
                prefix, weighted, np.uint64(self.prefix_count), np.int32(self.stride),
                np.uint64(modulus), matrix,
            ),
            shared_mem=self.block * 8,
        )
        cp.cuda.runtime.deviceSynchronize()
        host = cp.asnumpy(matrix)
        return np.asarray(
            [
                sum(int(host[left, degree - left]) for left in range(degree + 1))
                % modulus
                for degree in range(self.stride)
            ],
            dtype=np.uint64,
        )


def character_profile(
    counts: np.ndarray, residue_sum: int, root: int, character: int, modulus: int
) -> np.ndarray:
    powers = np.asarray(
        [pow(root, character * phase, modulus) for phase in range(P)], dtype=np.uint64
    )
    values = counts[:, residue_sum].astype(np.uint64) * powers
    return values.sum(axis=2, dtype=np.uint64) % np.uint64(modulus)


def crt(residues: list[int], moduli: list[int]) -> int:
    value = 0
    product = 1
    for residue, modulus in zip(residues, moduli, strict=True):
        correction = (residue - value) % modulus
        correction = correction * pow(product % modulus, -1, modulus) % modulus
        value += product * correction
        product *= modulus
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tuples", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--max-e", type=int, default=87)
    parser.add_argument("--moduli", default=",".join(map(str, DEFAULT_MODULI)))
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    moduli = [int(value) for value in args.moduli.split(",") if value]
    if not moduli or any(not is_prime(modulus) for modulus in moduli):
        raise ValueError("all moduli must be prime")

    with np.load(args.tuples) as archive:
        keys = np.asarray(archive["keys"], dtype=np.uint64)
        weights = np.asarray(archive["weights"], dtype=np.uint64)
    with np.load(args.profiles) as archive:
        counts = np.asarray(archive["counts"], dtype=np.uint64)
    max_k = counts.shape[2] - 1
    if counts.shape[:2] != (604, P) or counts.shape[3] != P:
        raise ArithmeticError(f"unexpected profile table shape {counts.shape}")

    prefix_triples, suffix_triples, prefix_index, suffix_index = tuple_topology(keys)
    contractor = GPUContractor(
        prefix_triples, suffix_triples, prefix_index, suffix_index, weights,
        max_k, args.device,
    )
    unrestricted = unrestricted_profile_products(counts, max_k)
    bounds, second_bounds, fourth_bounds, terms = theta_bound_and_terms(
        unrestricted, args.max_e, max_k
    )
    residues_by_exponent: list[list[int]] = [[] for _ in range(args.max_e + 1)]
    second_residues_by_exponent: list[list[int]] = [
        [] for _ in range(args.max_e + 1)
    ]
    fourth_residues_by_exponent: list[list[int]] = [
        [] for _ in range(args.max_e + 1)
    ]
    modulus_rows = []

    for modulus in moduli:
        modulus_started = time.monotonic()
        root = primitive_eleventh_root(modulus)
        fixed_sum = np.zeros((P, max_k + 1), dtype=np.uint64)
        for residue_sum in range(P):
            nonzero = np.zeros(max_k + 1, dtype=np.uint64)
            for character in range(1, P):
                single = character_profile(
                    counts, residue_sum, root, character, modulus
                )
                contracted = contractor.contract(single, modulus)
                nonzero = (nonzero + contracted) % np.uint64(modulus)
            inverse_filter = pow(pow(P, 10, modulus), -1, modulus)
            for excess in range(max_k + 1):
                numerator = unrestricted[residue_sum][excess] % modulus
                numerator = (numerator + int(nonzero[excess])) % modulus
                fixed_sum[residue_sum, excess] = numerator * inverse_filter % modulus

        theta_residues = []
        second_residues = []
        fourth_residues = []
        for exponent in range(args.max_e + 1):
            value = sum(
                int(fixed_sum[residue, excess])
                for _common_sum, residue, excess in terms[exponent]
            ) % modulus
            second = sum(
                common_sum**2 * int(fixed_sum[residue, excess])
                for common_sum, residue, excess in terms[exponent]
            ) % modulus
            fourth = sum(
                common_sum**4 * int(fixed_sum[residue, excess])
                for common_sum, residue, excess in terms[exponent]
            ) % modulus
            theta_residues.append(value)
            second_residues.append(second)
            fourth_residues.append(fourth)
            residues_by_exponent[exponent].append(value)
            second_residues_by_exponent[exponent].append(second)
            fourth_residues_by_exponent[exponent].append(fourth)
        failed_known = {
            str(exponent): {"expected": expected, "residue": theta_residues[exponent]}
            for exponent, expected in KNOWN_COEFFICIENTS.items()
            if exponent <= args.max_e and theta_residues[exponent] != expected % modulus
        }
        if failed_known:
            raise ArithmeticError(f"known theta calibration failed modulo {modulus}: {failed_known}")
        modulus_rows.append(
            {
                "modulus": modulus,
                "primitive_eleventh_root": root,
                "known_coefficients_match": True,
                "elapsed_seconds": time.monotonic() - modulus_started,
            }
        )
        print(json.dumps(modulus_rows[-1]), flush=True)

    modulus_product = math.prod(moduli)
    maximum_bound = max(bounds)
    maximum_second_bound = max(second_bounds)
    maximum_fourth_bound = max(fourth_bounds)
    maximum_reconstruction_bound = max(
        maximum_bound, maximum_second_bound, maximum_fourth_bound
    )
    if modulus_product <= maximum_reconstruction_bound:
        raise ArithmeticError(
            f"CRT product {modulus_product} does not exceed bound "
            f"{maximum_reconstruction_bound}"
        )
    coefficients = [crt(residues, moduli) for residues in residues_by_exponent]
    second_moments = [
        crt(residues, moduli) for residues in second_residues_by_exponent
    ]
    fourth_moments = [
        crt(residues, moduli) for residues in fourth_residues_by_exponent
    ]
    if any(coefficient > bound for coefficient, bound in zip(coefficients, bounds)):
        bad = next(
            index for index, (coefficient, bound) in enumerate(zip(coefficients, bounds))
            if coefficient > bound
        )
        raise ArithmeticError(
            f"CRT coefficient at {bad} exceeds unrestricted bound: {coefficients[bad]} > {bounds[bad]}"
        )
    failed_exact = {
        str(exponent): {"expected": expected, "actual": coefficients[exponent]}
        for exponent, expected in KNOWN_COEFFICIENTS.items()
        if exponent <= args.max_e and coefficients[exponent] != expected
    }
    if failed_exact:
        raise ArithmeticError(f"exact known theta calibration failed: {failed_exact}")

    for name, values, value_bounds in (
        ("second", second_moments, second_bounds),
        ("fourth", fourth_moments, fourth_bounds),
    ):
        if any(value > bound for value, bound in zip(values, value_bounds)):
            bad = next(
                index
                for index, (value, bound) in enumerate(zip(values, value_bounds))
                if value > bound
            )
            raise ArithmeticError(
                f"CRT {name} moment at {bad} exceeds unrestricted bound: "
                f"{values[bad]} > {value_bounds[bad]}"
            )

    n = P * P + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    second_moment_tight_frame_checks = []
    raw_trace = []
    harmonic_trace = []
    for exponent, (count, second, fourth) in enumerate(
        zip(coefficients, second_moments, fourth_moments)
    ):
        expected_second = 2 * P * exponent * count
        second_ok = n * second == expected_second
        if not second_ok:
            raise ArithmeticError(
                f"common-sum second moment fails the tight-frame identity at "
                f"{exponent}: {n * second} != {expected_second}"
            )
        if count:
            second_moment_tight_frame_checks.append(exponent)
        radius_sq = Fraction(exponent, 2 * P)
        tau = Fraction(n * (count * exponent * exponent - fourth), 4 * P * P * (P * P - 1))
        rho = Fraction(2 * count, d * (d + 2)) * radius_sq**2
        raw_trace.append(tau)
        harmonic_trace.append(tau - zdim * rho)

    known_raw_trace = {
        11: Fraction(0),
        # Directly summing ||proj_Z(xx^T)||^2 over the classified shells
        # gives these unscaled lattice-shell traces.  The older Prop. 15.665
        # audit accidentally combined H(x/2) at e=20,24 with the radial
        # correction for H(x); those two harmonic rows differ by 16.
        20: Fraction(89792, 11),
        24: Fraction(7076),
        27: Fraction(538752),
    }
    failed_trace = {
        str(exponent): {"expected": str(expected), "actual": str(raw_trace[exponent])}
        for exponent, expected in known_raw_trace.items()
        if exponent <= args.max_e and raw_trace[exponent] != expected
    }
    if failed_trace:
        raise ArithmeticError(f"known raw-trace calibration failed: {failed_trace}")

    device_properties = contractor.cp.cuda.runtime.getDeviceProperties(args.device)
    gpu_name = device_properties["name"]
    if isinstance(gpu_name, bytes):
        gpu_name = gpu_name.decode()
    report = {
        "experiment": "r1_p11_profile_theta_gpu",
        "status": "complete_exact_theta_prefix",
        "p": P,
        "gpu_name": str(gpu_name),
        "max_exponent": args.max_e,
        "max_excess_parameter_k": max_k,
        "weighted_tuple_edges": int(len(keys)),
        "prefix_triple_types": int(len(prefix_triples)),
        "suffix_triple_types": int(len(suffix_triples)),
        "moduli": modulus_rows,
        "crt_modulus_product": modulus_product,
        "maximum_unrestricted_bound": maximum_bound,
        "maximum_unrestricted_second_moment_bound": maximum_second_bound,
        "maximum_unrestricted_fourth_moment_bound": maximum_fourth_bound,
        "crt_product_exceeds_every_bound": True,
        "known_coefficients": {str(key): value for key, value in KNOWN_COEFFICIENTS.items()},
        "known_coefficients_match": True,
        "theta_coefficients": coefficients,
        "common_sum_second_moments": second_moments,
        "common_sum_fourth_moments": fourth_moments,
        "second_moment_tight_frame_identity_checked_at_nonempty_exponents": (
            second_moment_tight_frame_checks
        ),
        "raw_quartic_trace_coefficients": [str(value) for value in raw_trace],
        "harmonic_trace_coefficients": [str(value) for value in harmonic_trace],
        "known_raw_trace_coefficients": {
            str(exponent): str(value) for exponent, value in known_raw_trace.items()
        },
        "known_raw_trace_coefficients_match": True,
        "unrestricted_bounds": bounds,
        "unrestricted_second_moment_bounds": second_bounds,
        "unrestricted_fourth_moment_bounds": fourth_bounds,
        "elapsed_seconds": time.monotonic() - started,
        "inputs": {
            "tuples": str(args.tuples),
            "tuples_sha256": sha256(args.tuples),
            "profiles": str(args.profiles),
            "profiles_sha256": sha256(args.profiles),
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
