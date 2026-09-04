"""CPU-side gates for the independent p=31 OpenCL centre search."""

from __future__ import annotations

import numpy as np

from scripts.residual_branch_c_center_transverse_gpu import build_exact_tables
from scripts.residual_branch_c_center_transverse_opencl import (
    V100_REFERENCE_GRAPH_SHA256,
    _lexicographic_choice,
    _reference_cpu_replays,
)


def test_reference_graphs_replay_before_opencl_search() -> None:
    result = _reference_cpu_replays(build_exact_tables())
    baseline = result["published_5068_194"]
    winner = result["v100_4604_178"]
    assert (
        baseline["total_positive_l1_excess"],
        baseline["maximum_l1_excess"],
    ) == (5068, 194)
    assert (
        winner["total_positive_l1_excess"],
        winner["maximum_l1_excess"],
        winner["graph_sha256"],
    ) == (4604, 178, V100_REFERENCE_GRAPH_SHA256)


def test_host_selector_obeys_objective_validity_and_stable_ties() -> None:
    totals = np.asarray([[8, 7, 7], [3, 4, 2]], dtype=np.int32)
    maxima = np.asarray([[2, 5, 4], [9, 1, 8]], dtype=np.int32)
    valid = np.asarray([[True, True, True], [True, True, False]])
    assert _lexicographic_choice(totals, maxima, valid, "total").tolist() == [2, 0]
    assert _lexicographic_choice(totals, maxima, valid, "max").tolist() == [0, 1]
