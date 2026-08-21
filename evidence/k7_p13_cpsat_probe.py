#!/usr/bin/env python3
"""Coupled CP-SAT probe for the p=13, k=7 Max+ profile stratum.

The independent depressed-quintic type lists are far too large for a
Cartesian product.  This model keeps one finite-domain profile per square
direction and imposes all constraints at once:

* the degree 3, 2, and 1 homogeneous coefficient kernels;
* the constant reconstruction congruence and total normalized energy 21;
* the 169 simultaneous Boolean ridge-reconstruction equations.

Every solution is a translation-gauged epsilon=+1 Max+ representative with
nonzero degree-five kernel scalar.  The script is a finite diagnostic, not a
general-p QVAR proof.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k5_p29_coefficient_sieve import (  # noqa: E402
    homogeneous_matrix,
    kernel_modp,
    square_directions,
)
from k7_quintic_profile_probe import relevant_quintic_types  # noqa: E402

P = 13
K = 7


def lifted_profile_rows(p: int, cutoff: int) -> dict[int, list[tuple[int, ...]]]:
    """Return allowed (c,d,e,f,energy,h(0),...,h(p-1)) rows by leading term."""
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    rows: dict[int, set[tuple[int, ...]]] = defaultdict(set)
    for energy, leading, cubic, quadratic, linear, constant in (
        relevant_quintic_types(p, cutoff)
    ):
        polynomial = (
            leading * s**5
            + cubic * s**3
            + quadratic * s**2
            + linear * s
            + constant
        ) % p
        centered = np.where(polynomial <= midpoint, polynomial, polynomial - p)
        replacements = int(np.sum(centered) // p)
        endpoints = np.where(polynomial == midpoint)[0]
        for endpoint_set in itertools.combinations(endpoints, replacements):
            profile = centered.copy()
            if endpoint_set:
                profile[list(endpoint_set)] -= p
            if int(np.sum(profile)) != 0:
                raise RuntimeError("nonzero profile sum")
            if int(np.sum(profile * profile)) != 2 * p * energy:
                raise RuntimeError("profile energy mismatch")
            rows[int(leading)].add(
                (
                    int(cubic),
                    int(quadratic),
                    int(linear),
                    int(constant),
                    int(energy),
                    *map(int, profile),
                )
            )
    return {leading: sorted(values) for leading, values in rows.items()}


def add_multiple_of_p(
    model: cp_model.CpModel,
    expression: cp_model.LinearExpr,
    p: int,
    upper: int,
    name: str,
) -> None:
    quotient = model.new_int_var(0, upper // p + 1, name)
    model.add(expression == p * quotient)


def build_model(top_scalar: int) -> tuple[
    cp_model.CpModel,
    list[cp_model.IntVar],
    dict,
]:
    p = P
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    cutoff = total - (K - 1)
    square = square_directions(p)
    if len(square) != K:
        raise RuntimeError(f"expected {K} square directions, got {len(square)}")
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    top_kernel = kernel_modp(homogeneous_matrix(forms, 5, p), p)
    if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
        raise RuntimeError("unexpected degree-five kernel")
    leading = top_scalar * top_kernel[0] % p
    rows_by_leading = lifted_profile_rows(p, cutoff)

    model = cp_model.CpModel()
    coefficient_vars: dict[str, list[cp_model.IntVar]] = {
        key: [] for key in ("cubic", "quadratic", "linear", "constant")
    }
    energies: list[cp_model.IntVar] = []
    profile_values: list[list[cp_model.IntVar]] = []
    domain_sizes = []

    for direction in range(K):
        c = model.new_int_var(0, p - 1, f"c_{direction}")
        d = model.new_int_var(0, p - 1, f"d_{direction}")
        e = model.new_int_var(0, p - 1, f"e_{direction}")
        f = model.new_int_var(0, p - 1, f"f_{direction}")
        energy = model.new_int_var(1, cutoff, f"energy_{direction}")
        values = [
            model.new_int_var(-midpoint - 1, midpoint, f"h_{direction}_{s}")
            for s in range(p)
        ]
        allowed = rows_by_leading[int(leading[direction])]
        model.add_allowed_assignments([c, d, e, f, energy, *values], allowed)
        coefficient_vars["cubic"].append(c)
        coefficient_vars["quadratic"].append(d)
        coefficient_vars["linear"].append(e)
        coefficient_vars["constant"].append(f)
        energies.append(energy)
        profile_values.append(values)
        domain_sizes.append(len(allowed))

    model.add(sum(energies) == total)
    for degree, key in ((3, "cubic"), (2, "quadratic"), (1, "linear")):
        matrix = homogeneous_matrix(forms, degree, p)
        for row_index, matrix_row in enumerate(matrix):
            expression = sum(
                int(matrix_row[j]) * coefficient_vars[key][j] for j in range(K)
            )
            add_multiple_of_p(
                model,
                expression,
                p,
                int(np.sum(matrix_row)) * (p - 1),
                f"q_{key}_{row_index}",
            )

    required_rho_constant_sum = (-(K + 1) * pow(2, p - 2, p)) % p
    required_depressed_constant_sum = (
        required_rho_constant_sum - K * midpoint
    ) % p
    constant_sum = sum(coefficient_vars["constant"])
    constant_quotient = model.new_int_var(0, K, "q_constant")
    model.add(
        constant_sum
        == required_depressed_constant_sum + p * constant_quotient
    )

    negative = [model.new_bool_var(f"negative_{x}") for x in range(p * p)]
    for x in range(p * p):
        ridge_sum = sum(
            profile_values[j][int(coordinates[j][x])] for j in range(K)
        )
        model.add(ridge_sum + p * negative[x] == midpoint)

    metadata = {
        "p": p,
        "k": K,
        "top_scalar": top_scalar,
        "leading": list(map(int, leading)),
        "normalized_total_T": total,
        "relevant_cutoff": cutoff,
        "profile_domain_sizes": domain_sizes,
        "required_depressed_constant_sum_mod_p": required_depressed_constant_sum,
    }
    return model, negative, metadata


class SolutionProbe(cp_model.CpSolverSolutionCallback):
    def __init__(
        self,
        negative: list[cp_model.IntVar],
        kernel_real: np.ndarray,
        kernel_imag: np.ndarray,
        max_solutions: int,
        stop_at_limit: bool,
    ) -> None:
        super().__init__()
        self.negative = negative
        self.kernel_real = kernel_real.astype(np.int64)
        self.kernel_imag = kernel_imag.astype(np.int64)
        self.max_solutions = max_solutions
        self.stop_at_limit = stop_at_limit
        self.count = 0
        self.abs_z_sq_histogram: Counter[int] = Counter()
        self.min_abs_z_sq: int | None = None
        self.best_real: int | None = None
        self.best_imag: int | None = None
        self.best_negative_indices: list[int] | None = None

    def on_solution_callback(self) -> None:
        bits = np.fromiter(
            (self.value(value) for value in self.negative),
            dtype=np.int64,
            count=len(self.negative),
        )
        real = int(bits @ self.kernel_real @ bits)
        imag = int(bits @ self.kernel_imag @ bits)
        abs_sq = real * real + imag * imag
        self.count += 1
        self.abs_z_sq_histogram[abs_sq] += 1
        if self.min_abs_z_sq is None or abs_sq < self.min_abs_z_sq:
            self.min_abs_z_sq = abs_sq
            self.best_real = real
            self.best_imag = imag
            self.best_negative_indices = list(map(int, np.flatnonzero(bits)))
        if self.count <= 10 or self.count % 1000 == 0:
            print(
                f"solution={self.count} |Zpsi|^2={abs_sq} min={self.min_abs_z_sq}",
                flush=True,
            )
        if self.stop_at_limit and self.count >= self.max_solutions:
            self.stop_search()


def add_linf_quartic_objective(
    model: cp_model.CpModel,
    negative: list[cp_model.IntVar],
    kernel_real: np.ndarray,
    kernel_imag: np.ndarray,
) -> cp_model.IntVar:
    """Minimize max(|Re Zpsi|, |Im Zpsi|) using exact pair products."""
    real_terms = []
    imag_terms = []
    q = len(negative)
    for i in range(q):
        for j in range(i + 1, q):
            real_coefficient = 2 * int(kernel_real[i, j])
            imag_coefficient = 2 * int(kernel_imag[i, j])
            if not real_coefficient and not imag_coefficient:
                continue
            pair = model.new_bool_var(f"negative_pair_{i}_{j}")
            model.add(pair <= negative[i])
            model.add(pair <= negative[j])
            model.add(pair >= negative[i] + negative[j] - 1)
            if real_coefficient:
                real_terms.append(real_coefficient * pair)
            if imag_coefficient:
                imag_terms.append(imag_coefficient * pair)
    bound = q * (q - 1)
    real = model.new_int_var(-bound, bound, "Zpsi_real")
    imag = model.new_int_var(-bound, bound, "Zpsi_imag")
    model.add(real == sum(real_terms))
    model.add(imag == sum(imag_terms))
    abs_real = model.new_int_var(0, bound, "abs_Zpsi_real")
    abs_imag = model.new_int_var(0, bound, "abs_Zpsi_imag")
    model.add_abs_equality(abs_real, real)
    model.add_abs_equality(abs_imag, imag)
    linf = model.new_int_var(0, bound, "Zpsi_linf")
    model.add_max_equality(linf, [abs_real, abs_imag])
    model.minimize(linf)
    return linf


def scan(
    top_scalar: int,
    time_limit: float,
    max_solutions: int,
    workers: int,
    seed: int,
    optimize_z_linf: bool,
) -> dict:
    start = time.time()
    model, negative, metadata = build_model(top_scalar)
    build_seconds = time.time() - start
    kernel_real, kernel_imag = quartic_kernel(P)
    objective = None
    if optimize_z_linf:
        objective = add_linf_quartic_objective(
            model, negative, kernel_real, kernel_imag
        )
    callback = SolutionProbe(
        negative,
        kernel_real,
        kernel_imag,
        max_solutions=max_solutions,
        stop_at_limit=not optimize_z_linf,
    )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    enumerate_all = max_solutions > 1 and not optimize_z_linf
    solver.parameters.enumerate_all_solutions = enumerate_all
    solver.parameters.num_search_workers = 1 if enumerate_all else workers
    solver.parameters.cp_model_presolve = True
    solver.parameters.random_seed = seed
    solver.parameters.randomize_search = True
    status = solver.solve(model, callback)
    threshold_num = 3 * P * P * (P * P - 1)
    threshold_den = 16
    return {
        **metadata,
        "algorithm": "CP-SAT coupled coefficient/profile/Boolean reconstruction",
        "status": solver.status_name(status),
        "build_seconds": build_seconds,
        "solve_seconds": solver.wall_time,
        "solutions_seen": callback.count,
        "workers": 1 if enumerate_all else workers,
        "random_seed": seed,
        "optimized_Zpsi_linf": optimize_z_linf,
        "best_Zpsi_linf": (
            int(round(solver.objective_value))
            if optimize_z_linf and status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
            else None
        ),
        "best_bound_Zpsi_linf": (
            int(np.floor(solver.best_objective_bound))
            if optimize_z_linf and status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
            else None
        ),
        "stopped_at_solution_limit": callback.count >= max_solutions,
        "min_abs_Zpsi_sq_seen": callback.min_abs_z_sq,
        "best_Zpsi_seen": {
            "real": callback.best_real,
            "imag": callback.best_imag,
            "negative_indices": callback.best_negative_indices,
        },
        "abs_Zpsi_sq_histogram_seen": {
            str(value): count
            for value, count in sorted(callback.abs_z_sq_histogram.items())
        },
        "QVAR_threshold": f"{threshold_num}/{threshold_den}",
        "pointwise_QVAR_counterexample_seen": bool(
            callback.min_abs_z_sq is not None
            and threshold_den * callback.min_abs_z_sq < threshold_num
        ),
        "complete_enumeration": bool(
            enumerate_all
            and status == cp_model.OPTIMAL
            and callback.count < max_solutions
        ),
        "finite_diagnostic_only": True,
    }


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-scalar", type=int, default=1)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--max-solutions", type=int, default=10_000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--optimize-z-linf", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 1 <= args.top_scalar < P:
        raise ValueError("top scalar must lie in 1..12")
    report = scan(
        args.top_scalar,
        args.time_limit,
        args.max_solutions,
        args.workers,
        args.seed,
        args.optimize_z_linf,
    )
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    main()
