import sys
from fractions import Fraction
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from r1_p11_trace_reconstruct import solve_affine_prefix  # noqa: E402


def test_trace_affine_prefix_exact_rref_and_pivot_exponents():
    matrix = [
        (Fraction(1), Fraction(0)),
        (Fraction(2), Fraction(0)),
        (Fraction(0), Fraction(3)),
    ]
    coordinates, pivots, full_rank = solve_affine_prefix(
        matrix, [Fraction(2), Fraction(4), Fraction(9)], fixed_through=2
    )
    assert coordinates == [Fraction(2), Fraction(3)]
    assert pivots == [0, 2]
    assert full_rank == 2


def test_trace_affine_prefix_rejects_inconsistent_dependent_row():
    matrix = [(Fraction(1),), (Fraction(2),)]
    with pytest.raises(ArithmeticError, match="inconsistent"):
        solve_affine_prefix(
            matrix, [Fraction(2), Fraction(5)], fixed_through=1
        )


def test_profile_cuda_local_capacity_tracks_requested_stride():
    source = (SCRIPTS / "r1_p11_profile_theta_gpu.py").read_text()
    assert "unsigned long long current[PROFILE_STRIDE_CAPACITY]" in source
    assert 'CUDA_SOURCE.replace("PROFILE_STRIDE_CAPACITY", str(self.stride))' in source
    assert "unsigned long long current[22]" not in source
