"""Exact algebra for the conditional all-size threshold-valley reduction.

This is not a residual-(ii) closure. The geometric hypotheses and proofs are
in NOTE_2026-09-04_THRESHOLD_VALLEY_ACTIVE_GEOMETRY.md. No eigenshell,
signing, or prime census is performed here.
"""
from fractions import Fraction


def _odd_parameter(p: int) -> None:
    if type(p) is not int or p < 5 or p % 2 != 1:
        raise ValueError("need an odd integer p >= 5; Paley applicability is separate")


def valley_parameters(p: int, shell_floor: int) -> tuple[Fraction, Fraction]:
    """Return lambda, gamma for the r=3,4,5 conditional shell floors."""
    _odd_parameter(p)
    if type(shell_floor) is not int or shell_floor not in (3, 4, 5):
        raise ValueError("shell floor must be 3, 4, or 5")
    parameter = Fraction(shell_floor - 2, p + shell_floor - 2)
    margin = 2 * (p - 2) * parameter
    return parameter, margin


def interpolated_slack(
    p: int, shell_floor: int, conference_defect: int, active_slack: int
) -> Fraction:
    """M minus the signed interpolated score for one admissible row.

Here M=Phi(C)-4, delta=Phi(C)-eps Q_C(x), and
Delta=M-eps Q_A(x). Hypotheses: Delta>=0; either delta=0 with
Delta>=2r-4, or delta>=2p. Scores use the unordered-edge normalization.
"""
    parameter, _ = valley_parameters(p, shell_floor)
    if type(conference_defect) is not int or type(active_slack) is not int:
        raise ValueError("row defects and slacks must be integers")
    if active_slack < 0:
        raise ValueError("A row exceeds the asserted norm")
    if conference_defect == 0:
        if active_slack < 2 * shell_floor - 4:
            raise ValueError("conference-shell row violates the shell floor")
    elif conference_defect < 2 * p:
        raise ValueError("row violates the off-shell conference gap")
    return (1 - parameter) * active_slack + parameter * conference_defect - 4 * parameter


def restoration_error_second_moment(
    p: int, shell_floor: int, h: int, restored_edges: int
) -> Fraction:
    """Exact E[q_(A_D-A_lambda)^2], by distinct Walsh edge characters."""
    parameter, _ = valley_parameters(p, shell_floor)
    if type(h) is not int or type(restored_edges) is not int:
        raise ValueError("edge counts must be integers")
    if not 0 <= restored_edges <= h:
        raise ValueError("need 0 <= restored_edges <= h")
    return 4 * (parameter * parameter * h + (1 - 2 * parameter) * restored_edges)


def uniform_triangle_rounding_blocked(p: int, shell_floor: int, h: int) -> bool:
    """Whether Parseval rules out error norm <= gamma for every nonempty D.

True blocks only the uniform triangle-inequality certificate. It does NOT
say the restored signing has norm greater than M: row-dependent slack may
permit it.
"""
    _, margin = valley_parameters(p, shell_floor)
    if type(h) is not int or h < 1:
        raise ValueError("need positive disagreement size")
    return restoration_error_second_moment(p, shell_floor, h, 1) > margin * margin


def first_shell_cut_lower_bound(
    p: int, distance: int, contains_exceptional: bool, boolean_anchor: bool
) -> int:
    """Signed C-cut lower bound between same-phase first-shell states.

For distance one at a good exceptional coordinate this bound is negative;
the exact exceptional local field -p separately excludes a zero cut.
"""
    _odd_parameter(p)
    if type(distance) is not int or not 1 <= distance <= p * p + 1:
        raise ValueError("invalid nonzero Hamming distance")
    if type(contains_exceptional) is not bool or type(boolean_anchor) is not bool:
        raise ValueError("anchor case switches must be Boolean")
    if not contains_exceptional:
        return distance * (p - distance - 1)
    shift = -2 if boolean_anchor else 2
    return (distance + shift) * (p - distance + 1)


def restoration_row_slack(active_slack: int, signed_restored_sum: int) -> int:
    """M-eps Q_(A_D)(x); all signed rows must have nonnegative output."""
    if type(active_slack) is not int or active_slack < 0:
        raise ValueError("active slack must be a nonnegative integer")
    if type(signed_restored_sum) is not int:
        raise ValueError("signed restored-edge sum must be an integer")
    return active_slack - 2 * signed_restored_sum


SCOPE = {
    "conditional_fractional_valley": True,
    "integral_restoration_theorem": False,
    "all_size_minimal_witness_exclusion": False,
    "residual_ii_closed": False,
    "e1_closed_general": False,
    "L": False,
}
