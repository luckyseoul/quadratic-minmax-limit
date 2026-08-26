import itertools

from scripts.p7_no_infinity_unsaturated_cpsat import direction_target_options
from scripts.p7_no_infinity_unsaturated_partition_retry import (
    initial_intervals,
    interval_leaves,
)
from scripts.p7_unsaturated_slack_catalog import (
    exact_slack_catalog_values,
    exact_target_catalog_rows,
)


def test_unsaturated_catalog_counts_are_complete():
    expected = {
        (4, 1, 14): 36,
        (0, 0, 8): 1764,
        (2, 1, 14): 1764,
        (2, 0, 16): 1764,
        (4, 0, 16): 2233,
    }
    for key, count in expected.items():
        slacks = exact_slack_catalog_values(*key)
        targets = exact_target_catalog_rows(*key)
        assert len(slacks) == len(set(slacks)) == count
        assert len(targets) == len(set(targets)) == count


def test_negative_infinity_odd_fibre_catalog_counts_are_complete():
    expected = {
        (1, 1, 6): 1,
        (5, 1, 6): 1,
        (1, 1, 14): 1764,
        (5, 1, 14): 1764,
        (3, 1, 14): 36,
    }
    for key, count in expected.items():
        slacks = exact_slack_catalog_values(*key)
        targets = exact_target_catalog_rows(*key)
        assert len(slacks) == len(set(slacks)) == count
        assert len(targets) == len(set(targets)) == count


def test_three_1764_catalogs_are_universal_excess_translates():
    points = tuple(itertools.combinations(range(7), 4))
    excess = set(exact_slack_catalog_values(0, 0, 8))
    for phase, scaled_mean in ((0, 16), (1, 14)):
        parity = tuple(
            (sum(vertex in point for vertex in (0, 1)) + phase) & 1
            for point in points
        )
        translated = {
            tuple(row[index] - parity[index] for index in range(35))
            for row in exact_slack_catalog_values(2, phase, scaled_mean)
        }
        assert translated == excess


def test_b4_catalog_splits_as_1764_plus_448_plus_21():
    points = tuple(itertools.combinations(range(7), 4))
    odd_set = {0, 1, 2, 3}
    minimum = tuple(
        (len(set(point) & odd_set) - 2) ** 2 for point in points
    )
    regular = shallow = deep = 0
    for row in exact_slack_catalog_values(4, 0, 16):
        difference = tuple(row[index] - minimum[index] for index in range(35))
        negative = tuple(value for value in difference if value < 0)
        if not negative:
            assert all(value % 2 == 0 for value in difference)
            regular += 1
        elif negative == (-2,):
            shallow += 1
        elif negative == (-4,):
            deep += 1
        else:
            raise AssertionError(f"unexpected b=4 deficit pattern {negative}")
    assert (regular, shallow, deep) == (1764, 448, 21)


def test_fixed_elevation_target_tables_have_the_expected_scope():
    cases = {
        (2, 1, frozenset((0, 1)), 24): (1765, {6, 14}),
        (0, 0, frozenset(), 24): (1765, {0, 8}),
        (2, 0, frozenset((0, 1)), 24): (1765, {8, 16}),
        (4, 0, frozenset((0, 1, 2, 3)), 24): (2234, {8, 16}),
        (4, 1, frozenset((0, 1, 2, 3)), 32): (36, {14}),
    }
    for (b, phase, odd_set, floor_sum), (count, means) in cases.items():
        rows = direction_target_options(b, phase, set(odd_set), floor_sum)
        assert len(rows) == count
        assert {row[0] for row in rows} == means


def test_adaptive_catalog_partition_is_disjoint_and_resumable():
    roots = initial_intervals(10, 4)
    assert roots == ((0, 4), (4, 8), (8, 10))
    latest = {
        (0, 4): {"result": {"solver_status": "UNKNOWN"}},
        (0, 2): {"result": {"solver_status": "INFEASIBLE"}},
        (2, 4): {"result": {"solver_status": "UNKNOWN"}},
        (2, 3): {"result": {"solver_status": "INFEASIBLE"}},
        (3, 4): {"result": {"solver_status": "INFEASIBLE"}},
        (4, 8): {"result": {"solver_status": "INFEASIBLE"}},
        (8, 10): {"result": {"solver_status": "INFEASIBLE"}},
    }
    leaves = interval_leaves(roots, latest, 1)
    assert leaves == [
        (0, 2, "INFEASIBLE"),
        (2, 3, "INFEASIBLE"),
        (3, 4, "INFEASIBLE"),
        (4, 8, "INFEASIBLE"),
        (8, 10, "INFEASIBLE"),
    ]
    assert sum(stop - start for start, stop, _status in leaves) == 10
