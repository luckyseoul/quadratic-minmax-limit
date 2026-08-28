#!/usr/bin/env python3
"""Exact p=11 broad-channel shell masses from marked profile contractions.

Let ``O`` be the projected square-circle tensor operator.  For a shell raw
operator ``R_e`` this computes enough exact marked profile data to recover

    tr(R_e Pi_kernel), tr(R_e Pi_low), tr(R_e Pi_high).

The key reductions are as follows.  For square circles S put
``z_S=(x.w_S)^2`` and let ``M`` be point-circle incidence.  Evaluation on
projected tensors is the orthogonal projection ``y=Pi_ker(M^T) z``.  On
``ker(M^T)`` the circle Gram matrix is

    G = p^3(p+1) I + 2 p^2 A_2,

where ``A_2`` is adjacency by two-point circle intersection.  Therefore
the low broad mass is already determined by ``y^T A_2 y``.  The required
raw statistics reduce to:

* the sum of squared profile excesses ``sum_j k_j^2``;
* ``U4(a)=sum_c (sum_s eta(s-c)a_s)^4``.

Both are additive one-profile marks and are contracted over the complete
10-dimensional glue code with the same five-prime CRT as the ordinary
counter.  No floating-point value enters the reconstruction.
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
from r1_p11_profile_theta_gpu import (
    CUDA_SOURCE,
    DEFAULT_MODULI,
    GPUContractor,
    balanced_energy,
    crt,
    is_prime,
    primitive_eleventh_root,
    theta_bound_and_terms,
    tuple_topology,
    unrestricted_profile_products,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def character_profile_modular(
    tables: np.ndarray,
    residue_sum: int,
    root: int,
    character: int,
    modulus: int,
) -> np.ndarray:
    """Character transform with a reduction before every uint64 product.

    The ordinary count entries are small enough for the older vectorized
    multiply-and-sum implementation.  U4 entries reach 10^13 at k=30, so a
    31-bit character multiplier would overflow uint64 before the final
    reduction.  Reducing each table entry and each phase contribution first
    keeps every product below 4*10^18 and every addition below 4*10^9.
    """
    result = np.zeros(tables.shape[:1] + tables.shape[2:3], dtype=np.uint64)
    modulus_u = np.uint64(modulus)
    for phase in range(P):
        phase_power = np.uint64(pow(root, character * phase, modulus))
        values = tables[:, residue_sum, :, phase] % modulus_u
        contribution = values * phase_power % modulus_u
        result = (result + contribution) % modulus_u
    return result


MARKED_CUDA_SOURCE = r"""
extern "C" {

__device__ __forceinline__ unsigned long long marked_mul_mod(
    unsigned long long left,
    unsigned long long right,
    unsigned long long modulus
) {
    return (left * right) % modulus;
}

__device__ __forceinline__ void marked_add(
    unsigned long long *total,
    unsigned long long value,
    unsigned long long modulus
) {
    *total += value;
    if (*total >= modulus) *total %= modulus;
}

__global__ void build_triples_two_marked(
    unsigned long long count,
    const unsigned short *ids,
    const unsigned long long *single,
    const unsigned long long *mark1,
    const unsigned long long *mark2,
    int stride,
    unsigned long long modulus,
    unsigned long long *ordinary_output,
    unsigned long long *mark1_output,
    unsigned long long *mark2_output
) {
    unsigned long long row = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= count) return;
    unsigned long long ordinary[PROFILE_STRIDE_CAPACITY];
    unsigned long long first_mark[PROFILE_STRIDE_CAPACITY];
    unsigned long long second_mark[PROFILE_STRIDE_CAPACITY];
    unsigned long long next_ordinary[PROFILE_STRIDE_CAPACITY];
    unsigned long long next_first[PROFILE_STRIDE_CAPACITY];
    unsigned long long next_second[PROFILE_STRIDE_CAPACITY];
    unsigned int first = ids[3 * row];
    for (int degree = 0; degree < stride; ++degree) {
        unsigned long long offset = (unsigned long long)first * stride + degree;
        ordinary[degree] = single[offset];
        first_mark[degree] = mark1[offset];
        second_mark[degree] = mark2[offset];
    }
    for (int factor = 1; factor < 3; ++factor) {
        unsigned int identifier = ids[3 * row + factor];
        for (int degree = 0; degree < stride; ++degree) {
            unsigned long long total = 0ull;
            unsigned long long total_first = 0ull;
            unsigned long long total_second = 0ull;
            for (int left_degree = 0; left_degree <= degree; ++left_degree) {
                int right_degree = degree - left_degree;
                unsigned long long offset =
                    (unsigned long long)identifier * stride + right_degree;
                marked_add(&total, marked_mul_mod(
                    ordinary[left_degree], single[offset], modulus
                ), modulus);
                marked_add(&total_first, marked_mul_mod(
                    first_mark[left_degree], single[offset], modulus
                ), modulus);
                marked_add(&total_first, marked_mul_mod(
                    ordinary[left_degree], mark1[offset], modulus
                ), modulus);
                marked_add(&total_second, marked_mul_mod(
                    second_mark[left_degree], single[offset], modulus
                ), modulus);
                marked_add(&total_second, marked_mul_mod(
                    ordinary[left_degree], mark2[offset], modulus
                ), modulus);
            }
            next_ordinary[degree] = total % modulus;
            next_first[degree] = total_first % modulus;
            next_second[degree] = total_second % modulus;
        }
        for (int degree = 0; degree < stride; ++degree) {
            ordinary[degree] = next_ordinary[degree];
            first_mark[degree] = next_first[degree];
            second_mark[degree] = next_second[degree];
        }
    }
    for (int degree = 0; degree < stride; ++degree) {
        unsigned long long offset = row * stride + degree;
        ordinary_output[offset] = ordinary[degree];
        mark1_output[offset] = first_mark[degree];
        mark2_output[offset] = second_mark[degree];
    }
}

}
"""


class TwoMarkedGPUContractor(GPUContractor):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        source = MARKED_CUDA_SOURCE.replace(
            "PROFILE_STRIDE_CAPACITY", str(self.stride)
        )
        module = self.cp.RawModule(code=source, options=("--std=c++11",))
        self.build_two_marked = module.get_function("build_triples_two_marked")

    def triple_data(
        self,
        triples: object,
        count: int,
        single: object,
        mark1: object,
        mark2: object,
        modulus: int,
    ) -> tuple[object, object, object]:
        cp = self.cp
        outputs = tuple(
            cp.empty((count, self.stride), dtype=cp.uint64) for _ in range(3)
        )
        self.build_two_marked(
            ((count + self.block - 1) // self.block,),
            (self.block,),
            (
                np.uint64(count),
                triples,
                single,
                mark1,
                mark2,
                np.int32(self.stride),
                np.uint64(modulus),
                *outputs,
            ),
        )
        return outputs

    def _weighted_suffix(self, suffix: object, modulus: int) -> object:
        cp = self.cp
        output = cp.zeros((self.prefix_count, self.stride), dtype=cp.uint64)
        self.sparse(
            ((self.edge_count + self.block - 1) // self.block,),
            (self.block,),
            (
                np.uint64(self.edge_count),
                self.prefix_index,
                self.suffix_index,
                self.weights,
                suffix,
                np.int32(self.stride),
                np.uint64(modulus),
                output,
            ),
        )
        entry_count = self.prefix_count * self.stride
        self.reduce(
            ((entry_count + self.block - 1) // self.block,),
            (self.block,),
            (np.uint64(entry_count), np.uint64(modulus), output),
        )
        return output

    def _bilinear_polynomial(
        self, prefix: object, weighted_suffix: object, modulus: int
    ) -> np.ndarray:
        cp = self.cp
        matrix = cp.empty((self.stride, self.stride), dtype=cp.uint64)
        self.bilinear(
            (self.stride * self.stride,),
            (self.block,),
            (
                prefix,
                weighted_suffix,
                np.uint64(self.prefix_count),
                np.int32(self.stride),
                np.uint64(modulus),
                matrix,
            ),
            shared_mem=self.block * 8,
        )
        host = cp.asnumpy(matrix)
        return np.asarray(
            [
                sum(int(host[left, degree - left]) for left in range(degree + 1))
                % modulus
                for degree in range(self.stride)
            ],
            dtype=np.uint64,
        )

    def contract_two_marked(
        self,
        single: np.ndarray,
        mark1: np.ndarray,
        mark2: np.ndarray,
        modulus: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        cp = self.cp
        device_arrays = tuple(
            cp.asarray(value, dtype=cp.uint64) for value in (single, mark1, mark2)
        )
        prefix = self.triple_data(
            self.prefix_triples,
            self.prefix_count,
            *device_arrays,
            modulus,
        )
        suffix = self.triple_data(
            self.suffix_triples,
            self.suffix_count,
            *device_arrays,
            modulus,
        )
        weighted = tuple(self._weighted_suffix(value, modulus) for value in suffix)
        ordinary = self._bilinear_polynomial(prefix[0], weighted[0], modulus)
        first = (
            self._bilinear_polynomial(prefix[1], weighted[0], modulus)
            + self._bilinear_polynomial(prefix[0], weighted[1], modulus)
        ) % np.uint64(modulus)
        second = (
            self._bilinear_polynomial(prefix[2], weighted[0], modulus)
            + self._bilinear_polynomial(prefix[0], weighted[2], modulus)
        ) % np.uint64(modulus)
        cp.cuda.runtime.deviceSynchronize()
        return ordinary, first, second


def convolve(left: list[int], right: list[int], max_k: int) -> list[int]:
    return [
        sum(left[index] * right[degree - index] for index in range(degree + 1))
        for degree in range(max_k + 1)
    ]


def power_polynomial(values: list[int], exponent: int, max_k: int) -> list[int]:
    result = [1] + [0] * max_k
    for _ in range(exponent):
        result = convolve(result, values, max_k)
    return result


def unrestricted_marked_products(
    counts: np.ndarray, u4: np.ndarray, max_k: int
) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    ordinary_products = []
    square_excess_products = []
    u4_products = []
    degrees = list(range(max_k + 1))
    for residue_sum in range(P):
        ordinary = [
            sum(int(value) for value in counts[0, residue_sum, degree])
            for degree in degrees
        ]
        # Removing the phase condition makes both tables independent of the
        # rich residue polynomial.  Assert that rather than assume it.
        weighted_u4 = [
            sum(int(value) for value in u4[0, residue_sum, degree])
            for degree in degrees
        ]
        ordinary_power5 = power_polynomial(ordinary, 5, max_k)
        ordinary_products.append(convolve(ordinary_power5, ordinary, max_k))
        square_excess_products.append(
            [
                6 * value
                for value in convolve(
                    ordinary_power5,
                    [degree * degree * ordinary[degree] for degree in degrees],
                    max_k,
                )
            ]
        )
        u4_products.append(
            [
                6 * value
                for value in convolve(ordinary_power5, weighted_u4, max_k)
            ]
        )
    return ordinary_products, square_excess_products, u4_products


def projection_norm_sum(
    p: int,
    exponent: int,
    shell_count: int,
    coordinate_fourth_sum: Fraction,
) -> Fraction:
    n = p * p + 1
    radius = Fraction(exponent, 2 * p)
    incidence_perp = Fraction(p * p - 1, 2)
    incidence_constant = Fraction(p * (p + 1) ** 2, 2)
    u_norm = p**4 * (
        (n - 4) * shell_count * radius**2 + 4 * coordinate_fourth_sum
    )
    u_sum_squared = (
        shell_count * (p * p * (p * p - 1) * radius) ** 2
    )
    return (
        u_norm / incidence_perp
        + (Fraction(1, incidence_constant) - Fraction(1, incidence_perp))
        * u_sum_squared
        / n
    )


def broad_masses(
    exponent: int,
    count: int,
    common_fourth: int,
    profile_legendre_fourth: int,
    relation_two_raw: int,
) -> dict[str, Fraction]:
    p = P
    n = p * p + 1
    d = n // 2
    circle_count = p * n // 2
    relation_two_valency = p * (p * p - 1) // 4
    relation_two_point_eigenvalue = Fraction((p - 1) ** 2, 4)
    circle_low_eigenvalue = p**3 * (p - 1)
    circle_high_eigenvalue = p**3 * (p + 1)
    radius = Fraction(exponent, 2 * p)
    coordinate_fourth = Fraction(n * common_fourth, 16 * p**4)
    z_column_norm = projection_norm_sum(
        p, exponent, count, coordinate_fourth
    )
    z_sum = p * p * (p - 1) * radius
    z_square = Fraction(n, p + 1) * profile_legendre_fourth
    z_adjacency = d * relation_two_raw
    projected_square = z_square - z_column_norm
    projected_adjacency = (
        z_adjacency
        - relation_two_point_eigenvalue * z_column_norm
        - Fraction(
            relation_two_valency - relation_two_point_eigenvalue,
            circle_count,
        )
        * count
        * z_sum**2
    )
    low = -projected_adjacency / (p * circle_low_eigenvalue)
    high = (
        projected_square - circle_low_eigenvalue * low
    ) / circle_high_eigenvalue
    raw_trace = Fraction(
        n * (count * exponent * exponent - common_fourth),
        4 * p * p * (p * p - 1),
    )
    kernel = raw_trace - low - high
    return {
        "kernel": kernel,
        "low": low,
        "high": high,
        "raw_trace": raw_trace,
        "z_square": z_square,
        "z_adjacency": z_adjacency,
        "z_column_projection_norm": z_column_norm,
        "projected_circle_evaluation_norm": projected_square,
        "projected_intersection_two_adjacency": projected_adjacency,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tuples", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--ordinary-moments", type=Path, required=True)
    parser.add_argument("--max-e", type=int, default=120)
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
        u4 = np.asarray(archive["legendre_fourth"], dtype=np.uint64)
    if counts.shape != u4.shape or counts.shape[:2] != (1007, P):
        raise ArithmeticError(f"unexpected marked profile shape {counts.shape}")
    max_k = counts.shape[2] - 1
    unphased_count = counts.sum(axis=3, dtype=np.uint64)
    unphased_u4 = u4.sum(axis=3, dtype=np.uint64)
    if not np.all(unphased_count == unphased_count[0:1]):
        raise ArithmeticError("unphased count depends on rich profile type")
    if not np.all(unphased_u4 == unphased_u4[0:1]):
        raise ArithmeticError("unphased U4 depends on rich profile type")

    prefix_triples, suffix_triples, prefix_index, suffix_index = tuple_topology(keys)
    contractor = TwoMarkedGPUContractor(
        prefix_triples,
        suffix_triples,
        prefix_index,
        suffix_index,
        weights,
        max_k,
        args.device,
    )
    unrestricted, unrestricted_s2, unrestricted_u4 = unrestricted_marked_products(
        counts, u4, max_k
    )
    bounds, _second_bounds, _fourth_bounds, terms = theta_bound_and_terms(
        unrestricted, args.max_e, max_k
    )

    residue_count = [[0] * len(moduli) for _ in range(args.max_e + 1)]
    residue_s2 = [[0] * len(moduli) for _ in range(args.max_e + 1)]
    residue_profile_u4 = [[0] * len(moduli) for _ in range(args.max_e + 1)]
    # The low-channel relation-two statistic depends on the individual
    # common-sum term, so retain fixed arrays modulo each prime until the
    # exponent assembly below.
    fixed_by_modulus = []
    modulus_rows = []
    degree_square = np.arange(max_k + 1, dtype=np.uint64) ** 2
    inverse_denominators = []

    for modulus_index, modulus in enumerate(moduli):
        modulus_started = time.monotonic()
        root = primitive_eleventh_root(modulus)
        fixed_count = np.zeros((P, max_k + 1), dtype=np.uint64)
        fixed_s2 = np.zeros_like(fixed_count)
        fixed_u4 = np.zeros_like(fixed_count)
        inverse_filter = pow(pow(P, 10, modulus), -1, modulus)
        inverse_denominators.append(inverse_filter)
        for residue_sum in range(P):
            nonzero_count = np.zeros(max_k + 1, dtype=np.uint64)
            nonzero_s2 = np.zeros_like(nonzero_count)
            nonzero_u4 = np.zeros_like(nonzero_count)
            for character in range(1, P):
                single = character_profile_modular(
                    counts, residue_sum, root, character, modulus
                )
                mark_s2 = single * degree_square[None, :] % np.uint64(modulus)
                mark_u4 = character_profile_modular(
                    u4, residue_sum, root, character, modulus
                )
                contracted = contractor.contract_two_marked(
                    single, mark_s2, mark_u4, modulus
                )
                nonzero_count = (nonzero_count + contracted[0]) % np.uint64(modulus)
                nonzero_s2 = (nonzero_s2 + contracted[1]) % np.uint64(modulus)
                nonzero_u4 = (nonzero_u4 + contracted[2]) % np.uint64(modulus)
            for excess in range(max_k + 1):
                fixed_count[residue_sum, excess] = (
                    (unrestricted[residue_sum][excess] + int(nonzero_count[excess]))
                    % modulus
                    * inverse_filter
                    % modulus
                )
                fixed_s2[residue_sum, excess] = (
                    (unrestricted_s2[residue_sum][excess] + int(nonzero_s2[excess]))
                    % modulus
                    * inverse_filter
                    % modulus
                )
                fixed_u4[residue_sum, excess] = (
                    (unrestricted_u4[residue_sum][excess] + int(nonzero_u4[excess]))
                    % modulus
                    * inverse_filter
                    % modulus
                )
        fixed_by_modulus.append((fixed_count, fixed_s2, fixed_u4))
        for exponent in range(args.max_e + 1):
            residue_count[exponent][modulus_index] = sum(
                int(fixed_count[residue, excess])
                for _common_sum, residue, excess in terms[exponent]
            ) % modulus
            residue_s2[exponent][modulus_index] = sum(
                int(fixed_s2[residue, excess])
                for _common_sum, residue, excess in terms[exponent]
            ) % modulus
            residue_profile_u4[exponent][modulus_index] = sum(
                int(fixed_u4[residue, excess])
                for _common_sum, residue, excess in terms[exponent]
            ) % modulus
        modulus_rows.append(
            {
                "modulus": modulus,
                "primitive_eleventh_root": root,
                "elapsed_seconds": time.monotonic() - modulus_started,
            }
        )
        print(json.dumps(modulus_rows[-1]), flush=True)

    modulus_product = math.prod(moduli)
    count_values = [crt(values, moduli) for values in residue_count]
    s2_values = [crt(values, moduli) for values in residue_s2]
    profile_u4_values = [crt(values, moduli) for values in residue_profile_u4]
    s2_bounds = [6 * max_k * max_k * value for value in bounds]
    maximum_profile_u4 = 11 * (math.isqrt(P * ((P - 1) + 2 * max_k)) ** 4)
    profile_u4_bounds = [6 * maximum_profile_u4 * value for value in bounds]
    if modulus_product <= max(max(s2_bounds), max(profile_u4_bounds)):
        raise ArithmeticError("CRT product does not dominate marked-profile bounds")
    for name, values, value_bounds in (
        ("count", count_values, bounds),
        ("squared-excess", s2_values, s2_bounds),
        ("profile-U4", profile_u4_values, profile_u4_bounds),
    ):
        if any(value > bound for value, bound in zip(values, value_bounds)):
            bad = next(
                index
                for index, (value, bound) in enumerate(zip(values, value_bounds))
                if value > bound
            )
            raise ArithmeticError(
                f"CRT {name} value at {bad} exceeds bound: {values[bad]} > {value_bounds[bad]}"
            )

    # Reconstruct the relation-two profile statistic term by term.  Its
    # residues are assembled after CRT because the common sum enters the
    # integer polynomial multiplying each fixed profile count.
    fixed_exact = []
    for residue in range(P):
        exact_rows = []
        for excess in range(max_k + 1):
            exact_rows.append(
                tuple(
                    crt(
                        [int(fixed[part][residue, excess]) for fixed in fixed_by_modulus],
                        moduli,
                    )
                    for part in range(3)
                )
            )
        fixed_exact.append(exact_rows)

    relation_two_raw = []
    for exponent in range(args.max_e + 1):
        total = 0
        for common_sum, residue, excess in terms[exponent]:
            count, sum_k_squared, _weighted_u4 = fixed_exact[residue][excess]
            base = balanced_energy(common_sum)
            constant = P * base - common_sum * common_sum
            total_u2 = 6 * constant + 2 * P * excess
            sum_u2_squared = (
                count * (6 * constant * constant + 4 * P * constant * excess)
                + 4 * P * P * sum_k_squared
            )
            total += count * total_u2 * total_u2 - sum_u2_squared
        relation_two_raw.append(total)

    ordinary_report = json.loads(args.ordinary_moments.read_text())
    ordinary_counts = [int(value) for value in ordinary_report["theta_coefficients"]]
    common_fourth = [
        int(value) for value in ordinary_report["common_sum_fourth_moments"]
    ]
    if ordinary_counts[: args.max_e + 1] != count_values:
        mismatch = next(
            index
            for index, (left, right) in enumerate(
                zip(ordinary_counts, count_values)
            )
            if left != right
        )
        raise ArithmeticError(
            f"rich counter fails ordinary coefficient audit at {mismatch}"
        )

    mass_rows = []
    for exponent in range(args.max_e + 1):
        masses = broad_masses(
            exponent,
            count_values[exponent],
            common_fourth[exponent],
            profile_u4_values[exponent],
            relation_two_raw[exponent],
        )
        if any(masses[channel] < 0 for channel in ("kernel", "low", "high")):
            raise ArithmeticError(
                f"negative reconstructed broad mass at exponent {exponent}: {masses}"
            )
        mass_rows.append(
            {
                "exponent": exponent,
                "shell_count": count_values[exponent],
                "sum_profile_excess_squares": s2_values[exponent],
                "profile_legendre_fourth": profile_u4_values[exponent],
                "relation_two_raw": relation_two_raw[exponent],
                **{key: str(value) for key, value in masses.items()},
            }
        )

    # Independent calibration against the four classified shell operators.
    import sys

    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "src"))
    from e1_gmin_m4_prop15665 import broad_channel_dimensions, p11_early_shell_audit

    dimensions = broad_channel_dimensions(P)
    early = {row["scaled_norm"]: row for row in p11_early_shell_audit()}
    for exponent, row in early.items():
        for channel in ("kernel", "low", "high"):
            label = f"circle-{channel}"
            expected = Fraction(row["raw_channel_eigenvalues"][label]) * dimensions[label]
            actual = Fraction(mass_rows[exponent][channel])
            if actual != expected:
                raise ArithmeticError(
                    f"classified-shell {channel} mass mismatch at {exponent}: "
                    f"{actual} != {expected}"
                )

    report = {
        "experiment": "r1_p11_channel_moments_gpu",
        "status": "complete_exact_broad_channel_mass_prefix",
        "p": P,
        "max_exponent": args.max_e,
        "max_excess_parameter_k": max_k,
        "weighted_tuple_edges": int(len(keys)),
        "prefix_triple_types": int(len(prefix_triples)),
        "suffix_triple_types": int(len(suffix_triples)),
        "moduli": modulus_rows,
        "crt_modulus_product": modulus_product,
        "maximum_count_bound": max(bounds),
        "maximum_squared_excess_bound": max(s2_bounds),
        "maximum_profile_legendre_fourth_bound": max(profile_u4_bounds),
        "crt_product_exceeds_every_bound": True,
        "ordinary_coefficients_match_prior_exact_profile_report": True,
        "classified_shell_channel_masses_match": True,
        "broad_channel_dimensions": dimensions,
        "rows": mass_rows,
        "elapsed_seconds": time.monotonic() - started,
        "inputs": {
            "tuples": str(args.tuples),
            "tuples_sha256": sha256(args.tuples),
            "profiles": str(args.profiles),
            "profiles_sha256": sha256(args.profiles),
            "ordinary_moments": str(args.ordinary_moments),
            "ordinary_moments_sha256": sha256(args.ordinary_moments),
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "rows"}, indent=2))
    print("NONZERO_BROAD_MASSES=")
    for row in mass_rows:
        if row["shell_count"]:
            print(
                row["exponent"], row["kernel"], row["low"], row["high"],
                flush=True,
            )


if __name__ == "__main__":
    main()
