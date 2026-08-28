from fractions import Fraction as F
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from r1_p11_trace_coupled_exact_lp import (
    build_shellwise_conserved_model,
    component_cases,
    Constraint,
    ExactModel,
    maximum_weighted_variance_exact,
    verify_certificate,
    write_lp,
)
from r1_p11_scalar_coupled_exact_lp import scalar_affine_reduction


def toy_model() -> ExactModel:
    return ExactModel(
        target_base=F(7, 5),
        target=(F(1, 3), F(2)),
        constraints=(
            Constraint("c1", (F(1), F(1)), ">=", F(1, 2)),
            Constraint("c2", (F(1), F(0)), "<=", F(2)),
            Constraint("c3", (F(0), F(1)), "<=", F(3)),
        ),
        fixed_checks=(),
    )


def test_exact_primal_dual_verifier_for_both_senses(tmp_path):
    model = toy_model()
    minimum = tmp_path / "minimum.sol"
    minimum.write_text(
        """status = OPTIMAL
status OPTIMAL
\tValue = -7/3
VARS:
y1 = 2
y2 = -3/2
REDUCED COST:
PI:
c1 = 2
c2 = -5/3
SLACK:
"""
    )
    minimum_check = verify_certificate(model, model.target, minimum)
    assert minimum_check["solver_objective"] == "-7/3"
    assert minimum_check["primal_constraints_verified"] == 3
    assert minimum_check["dual_stationarity_equations_verified"] == 2

    maximum = tmp_path / "maximum.sol"
    maximum.write_text(
        """status = OPTIMAL
status OPTIMAL
\tValue = -20/3
VARS:
y1 = 2
y2 = 3
REDUCED COST:
PI:
c2 = -1/3
c3 = -2
SLACK:
"""
    )
    maximum_objective = tuple(-value for value in model.target)
    maximum_check = verify_certificate(model, maximum_objective, maximum)
    assert maximum_check["solver_objective"] == "-20/3"


def test_grouped_exact_variance_enumeration():
    intervals = [(F(0), F(2), 3)] * 3 + [(F(-1), F(1), 2)]
    result = maximum_weighted_variance_exact(intervals, F(5))
    assert result["normalized_mean"] == F(5, 11)
    assert result["normalized_variance_max"] == F(398, 33)
    assert result["interval_symmetry_groups"] == 2
    assert len(result["maximizer"]) == 4


def test_shellwise_model_uses_exact_symmetry_quotient_and_mass_identity():
    length = 30
    base = [F()] * length
    matrix = [(F(),)] * (length - 1) + [(F(1),)]
    reduction = (base, matrix, F(), (F(1),))
    reductions = {
        "circle-kernel": reduction,
        "circle-low": reduction,
        "circle-high": reduction,
    }
    counts = [0] * length
    counts[-1] = 1
    trace = [F()] * length
    cases = component_cases(11)
    model, representatives = build_shellwise_conserved_model(
        reductions, cases, 0, counts, trace, 11
    )

    assert len(representatives) == 5
    assert [row["multiplicity"] for row in representatives] == [1, 9, 1, 2, 2]
    assert model.target == (F(1), F(), F(), F(), F())
    mass = next(row for row in model.constraints if row.name == "s29_mass")
    assert mass.coefficients == (F(122), F(1098), F(61), F(244), F(244))
    assert mass.rhs == 0
    assert len([row for row in model.constraints if row.name.endswith("_positive")]) == 5
    assert len(model.fixed_checks) == 29


def test_scalar_affine_reduction_is_basis_independent():
    rows = [
        [F(1), F(1), F(0)],
        [F(0), F(0), F(0)],
        [F(0), F(0), F(0)],
        [F(1), F(-1), F(1)],
        [F(2), F(3), F(4)],
    ]
    base, matrix, pivots = scalar_affine_reduction(rows, 11, 3)

    assert pivots == (0, 1)
    assert base[:4] == [F(1), F(0), F(0), F(0)]
    assert matrix[:4] == [(F(0),)] * 4
    assert base[4] == F(5, 2)
    assert matrix[4] == (F(9, 2),)


def test_scalar_affine_reduction_accepts_additional_exact_shell_count():
    rows = [
        [F(1), F(0), F(0)],
        [F(0), F(1), F(0)],
        [F(0), F(0), F(0)],
        [F(0), F(0), F(1)],
        [F(1), F(1), F(1)],
    ]
    base, matrix, pivots = scalar_affine_reduction(
        rows, 11, 2, {3: F(7)}
    )

    assert pivots == (0, 1, 2)
    assert base == [F(1), F(0), F(0), F(7), F(8)]
    assert matrix == [()] * len(rows)


def test_qsopt_lp_writer_supports_zero_objective(tmp_path):
    model = ExactModel(
        target_base=F(),
        target=(F(),),
        constraints=(Constraint("unit", (F(1),), "=", F()),),
        fixed_checks=(),
    )
    path = tmp_path / "feasibility.lp"
    write_lp(path, model, "minimum")

    assert "obj: 0 y1" in path.read_text()
