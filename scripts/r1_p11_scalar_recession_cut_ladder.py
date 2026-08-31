#!/usr/bin/env python3
"""Find the exact scalar shell cuts needed to remove p=11 LP recession rays.

Only the homogeneous coefficient functionals matter to a recession cone, so
the numerical values assigned to newly fixed scalar coefficients do not affect
this calculation.  We use zero right-hand sides, add the earliest scalar
coefficient moved by each exact PPL ray, and rebuild until both target senses
have no improving recession direction or a pure-harmonic ray remains.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
import time
from fractions import Fraction
from pathlib import Path

from r1_p11_scalar_coupled_exact_lp import (
    Constraint,
    ExactModel,
    build_scalar_trace_budget_model,
    scalar_affine_reduction,
)
from r1_p11_scalar_coupled_lp import CHANNELS, load_rows
from r1_p11_scalar_coupled_ppl import (
    canonical_rows,
    solve_unbounded_ray,
    verify_unbounded_ray,
)
from r1_p11_trace_coupled_exact_lp import (
    parse_solution,
    verify_certificate as verify_qsopt_certificate,
    write_lp,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_unbounded_ray_qsopt(
    model: ExactModel,
    objective: tuple[Fraction, ...],
    esolver: Path,
) -> tuple[str, tuple[Fraction, ...] | None, dict[str, object]]:
    """Solve the normalized recession feasibility LP with QSopt_ex."""
    recession = ExactModel(
        target_base=Fraction(),
        target=(Fraction(),) * len(objective),
        constraints=tuple(
            Constraint(
                f"ray_{constraint.name}",
                constraint.coefficients,
                constraint.sense,
                Fraction(),
            )
            for constraint in model.constraints
        )
        + (Constraint("ray_objective_normalization", objective, "=", Fraction(-1)),),
        fixed_checks=(),
    )
    with tempfile.TemporaryDirectory(prefix="r1-qsopt-ray-", dir="/tmp") as raw_dir:
        directory = Path(raw_dir)
        lp_path = directory / "recession.lp"
        solution_path = directory / "recession.sol"
        written_objective = write_lp(lp_path, recession, "minimum")
        if any(written_objective):
            raise ArithmeticError("recession feasibility LP acquired an objective")
        result = subprocess.run(
            [str(esolver), "-L", "-O", str(solution_path), str(lp_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        log = result.stdout + result.stderr
        if result.returncode != 0 or not solution_path.is_file():
            raise RuntimeError(f"QSopt_ex recession solve failed:\n{log}")
        solution_text = solution_path.read_text()
        metadata: dict[str, object] = {
            "backend": "QSopt_ex",
            "lp_sha256": sha256(lp_path),
            "solver_log_tail": log.splitlines()[-12:],
        }
        if "status = INFEASIBLE" in solution_text:
            if "Problem Is Infeasible" not in log:
                raise ArithmeticError("QSopt_ex solution/log infeasibility mismatch")
            return "unfeasible", None, metadata
        if "status OPTIMAL" not in solution_text or "Problem Solved Exactly" not in log:
            raise ArithmeticError(f"unexpected QSopt_ex recession status:\n{solution_text}\n{log}")
        metadata["certificate"] = verify_qsopt_certificate(
            recession, written_objective, solution_path
        )
        _value, variables_by_name, _pi = parse_solution(solution_path)
        ray = tuple(
            variables_by_name.get(f"y{index}", Fraction())
            for index in range(1, len(objective) + 1)
        )
        return "optimized", ray, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scalar-qrows", type=Path, required=True)
    parser.add_argument("--affine-prefix", type=Path, required=True)
    parser.add_argument("--scalar-half-target-rows", type=Path, required=True)
    parser.add_argument("--scalar-half-target-first", type=Fraction, required=True)
    parser.add_argument("--p", type=int, default=11)
    parser.add_argument("--case", default="circle-kernel-principal")
    parser.add_argument("--initial-fixed-exponent", type=int, action="append", default=[])
    parser.add_argument(
        "--initial-fixed-pivots-through",
        type=int,
        help="also fix every scalar identity-pivot coefficient through this exponent",
    )
    parser.add_argument("--max-rounds", type=int, default=20)
    parser.add_argument(
        "--ray-backend", choices=("ppl", "qsopt"), default="qsopt"
    )
    parser.add_argument("--esolver", type=Path, default=Path("/usr/bin/esolver"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    fixed_through = 2 * (args.p + 3)
    if any(exponent <= fixed_through for exponent in args.initial_fixed_exponent):
        parser.error(
            f"initial fixed exponents must exceed the proved prefix {fixed_through}"
        )
    if args.max_rounds < 1:
        parser.error("--max-rounds must be positive")

    scalar_rows = load_rows(args.scalar_qrows)
    scalar_half_target_rows = load_rows(args.scalar_half_target_rows)
    harmonic_paths = {
        channel: Path(f"{args.affine_prefix}{channel}_qrows_exact_20260827.txt")
        for channel in CHANNELS
    }
    harmonic_rows = {
        channel: load_rows(path) for channel, path in harmonic_paths.items()
    }
    scalar_width = len(scalar_rows[0])
    identity_pivot_exponents: list[int] = []
    for column in range(scalar_width):
        expected = [Fraction()] * scalar_width
        expected[column] = Fraction(1)
        matches = [
            exponent for exponent, row in enumerate(scalar_rows) if row == expected
        ]
        if len(matches) != 1:
            raise ArithmeticError(
                f"scalar coordinate {column} has {len(matches)} identity-pivot rows"
            )
        identity_pivot_exponents.append(matches[0])
    initial_fixed_exponents = set(args.initial_fixed_exponent)
    if args.initial_fixed_pivots_through is not None:
        initial_fixed_exponents.update(
            exponent
            for exponent in identity_pivot_exponents
            if fixed_through < exponent <= args.initial_fixed_pivots_through
        )
    fixed_counts = {
        exponent: Fraction() for exponent in sorted(initial_fixed_exponents)
    }
    report: dict[str, object] = {
        "experiment": "r1_p11_scalar_recession_cut_ladder",
        "status": "running",
        "p": args.p,
        "case": args.case,
        "method": (
            f"{args.ray_backend} normalized recession rays; earliest moved scalar "
            "coefficient added after each round"
        ),
        "ray_backend": args.ray_backend,
        "scalar_qrows": str(args.scalar_qrows),
        "scalar_qrows_sha256": sha256(args.scalar_qrows),
        "scalar_half_target_rows": str(args.scalar_half_target_rows),
        "scalar_half_target_rows_sha256": sha256(args.scalar_half_target_rows),
        "scalar_half_target_first": str(args.scalar_half_target_first),
        "harmonic_qrows": {
            channel: {"path": str(path), "sha256": sha256(path)}
            for channel, path in harmonic_paths.items()
        },
        "initial_fixed_exponents": sorted(fixed_counts),
        "scalar_identity_pivot_exponents": identity_pivot_exponents,
        "initial_fixed_pivots_through": args.initial_fixed_pivots_through,
        "rounds": [],
    }
    started_all = time.monotonic()

    for round_index in range(args.max_rounds):
        model, representatives = build_scalar_trace_budget_model(
            scalar_rows,
            harmonic_rows,
            args.p,
            args.case,
            parity_fourth_moment=True,
            scalar_half_target_rows=scalar_half_target_rows,
            scalar_half_target_first=args.scalar_half_target_first,
            scalar_fixed_counts=fixed_counts,
        )
        _base, scalar_matrix, _pivots = scalar_affine_reduction(
            scalar_rows, args.p, fixed_through, fixed_counts
        )
        scalar_dimension = len(scalar_matrix[0])
        if scalar_dimension > len(model.target):
            raise ArithmeticError("scalar dimension exceeds the full model dimension")

        round_row: dict[str, object] = {
            "round": round_index,
            "fixed_exponents": sorted(fixed_counts),
            "variable_count": len(model.target),
            "scalar_variable_count": scalar_dimension,
            "canonical_constraint_count": len(canonical_rows(model)),
            "symmetry_representatives": representatives,
            "senses": {},
        }
        new_exponents: set[int] = set()
        pure_harmonic = False
        all_closed = True

        for sense in ("minimum", "maximum"):
            objective = (
                model.target
                if sense == "minimum"
                else tuple(-value for value in model.target)
            )
            started = time.monotonic()
            if args.ray_backend == "qsopt":
                ray_status, ray, solver_metadata = solve_unbounded_ray_qsopt(
                    model, objective, args.esolver
                )
            else:
                ray_status, ray = solve_unbounded_ray(model, objective)
                solver_metadata = {"backend": "PPL"}
            elapsed = time.monotonic() - started
            if ray_status == "unfeasible":
                round_row["senses"][sense] = {
                    "status": "exact_no_improving_recession_ray",
                    "solver_status": ray_status,
                    "seconds": elapsed,
                    "solver_metadata": solver_metadata,
                }
                continue
            all_closed = False
            if ray_status != "optimized" or ray is None:
                raise ArithmeticError(f"unexpected {sense} ray status: {ray_status}")

            certificate = verify_unbounded_ray(model, objective, ray)
            scalar_ray = ray[-scalar_dimension:] if scalar_dimension else ()
            coefficient_direction = [
                sum(
                    (coefficient * value for coefficient, value in zip(row, scalar_ray)),
                    Fraction(),
                )
                for row in scalar_matrix
            ]
            support = [
                exponent
                for exponent, value in enumerate(coefficient_direction)
                if value
            ]
            if not support:
                pure_harmonic = True
            else:
                new_exponents.add(support[0])
            round_row["senses"][sense] = {
                "status": "exact_unbounded_ray",
                "solver_status": ray_status,
                "seconds": elapsed,
                "nonzero_full_variables": certificate["nonzero_ray_variables"],
                "nonzero_scalar_variables": sum(value != 0 for value in scalar_ray),
                "first_scalar_support_exponents": support[:20],
                "selected_cut_exponent": support[0] if support else None,
                "objective_derivative": certificate["objective_derivative"],
                "recession_constraints_verified": certificate[
                    "recession_constraints_verified"
                ],
                "solver_metadata": solver_metadata,
                "ray": certificate["ray"],
            }

        report["rounds"].append(round_row)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(
            f"round={round_index} fixed={sorted(fixed_counts)} "
            f"new={sorted(new_exponents)} closed={all_closed} "
            f"pure_harmonic={pure_harmonic}",
            flush=True,
        )
        if all_closed:
            report["status"] = "exact_recession_cone_closed_for_both_senses"
            break
        if pure_harmonic:
            report["status"] = "pure_harmonic_recession_ray_survives"
            break
        genuinely_new = new_exponents - fixed_counts.keys()
        if not genuinely_new:
            raise ArithmeticError("ray support produced no new scalar cut")
        fixed_counts.update({exponent: Fraction() for exponent in genuinely_new})
    else:
        report["status"] = "round_limit_reached"

    report["final_fixed_exponents"] = sorted(fixed_counts)
    report["elapsed_seconds"] = time.monotonic() - started_all
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "status": report["status"],
        "final_fixed_exponents": report["final_fixed_exponents"],
        "round_count": len(report["rounds"]),
        "elapsed_seconds": report["elapsed_seconds"],
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
