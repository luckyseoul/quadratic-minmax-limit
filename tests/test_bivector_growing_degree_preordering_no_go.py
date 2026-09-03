"""Checks for the growing-degree squared-bivector preordering obstruction."""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "evidence" / "NOTE_2026-09-02_BIVECTOR_GROWING_DEGREE_PREORDERING_NO_GO.md"


def kappa() -> float:
    return math.pi**2 / (4.0 * (math.sqrt(2.0) - 1.0) ** 2)


def criterion(n: int, half_degree: int) -> float:
    d = half_degree
    return (2 * d + 1) ** d * ((1.0 + kappa() / (n - 1)) ** d - 1.0)


def first_certified_order(half_degree: int) -> int:
    d = half_degree
    c_d = (2 * d + 1) ** d
    denominator = (1.0 + 1.0 / c_d) ** (1.0 / d) - 1.0
    # Strict inequality is required.  Advancing from the real threshold by
    # two makes this insensitive to the threshold being integral.
    return math.floor(kappa() / denominator) + 2


def test_exact_constant_and_fixed_degree_thresholds() -> None:
    assert math.isclose(kappa(), 14.38106750045589, rel_tol=0.0, abs_tol=1e-12)
    assert first_certified_order(1) == 45
    assert first_certified_order(2) == 728
    assert first_certified_order(3) == 14814
    for d in range(1, 6):
        n0 = first_certified_order(d)
        assert criterion(n0, d) < 1.0
        assert criterion(n0 - 2, d) >= 1.0


def test_growing_degree_examples_are_monotone_in_order() -> None:
    for d in (1, 2, 3, 4, 5):
        n0 = first_certified_order(d)
        assert criterion(2 * n0, d) < criterion(n0, d) < 1.0


def test_note_records_scope_and_quantitative_gate() -> None:
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "degree-`2D` preordering",
        "(2D+1)^D",
        "Bonami--Beckner",
        "Khintchine",
        "log n",
        "different lifted encoding",
        "not an orientation construction",
        "raw polynomial degree",
        "For the signed outgoing-half system",
    )
    for token in required:
        assert token in text
