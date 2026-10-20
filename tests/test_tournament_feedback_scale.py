"""A regular tournament can force a super-cube feedback allowance."""
from __future__ import annotations

from minmax_quadratic import directed_triangle_count, nearly_regular_tournament


def _feedback_lower_bound(order: int) -> float:
    tournament = nearly_regular_tournament(order)
    return directed_triangle_count(tournament) / float(order - 2)


def test_regular_feedback_allowance_exceeds_the_cube_scale():
    crossed = False
    for order in range(5, 41, 2):
        allowance = 4.0 * _feedback_lower_bound(order)
        if allowance > order ** 1.5:
            crossed = True
            later = order + 4
            assert 4.0 * _feedback_lower_bound(later) > later ** 1.5
            break
    assert crossed
