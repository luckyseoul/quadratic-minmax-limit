"""Fekete's lemma does not apply to the certified Seidel minima."""
from __future__ import annotations

import json
from pathlib import Path

from minmax_quadratic import root_scale, two_thirds_scale


def _heights():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    return {
        int(order): int(row["m"])
        for order, row in payload.items()
        if order.isdigit()
    }


def test_root_scale_is_neither_superadditive_nor_subadditive():
    height = _heights()
    assert root_scale(height[6], 6) < root_scale(height[3], 3) + root_scale(height[3], 3)
    assert root_scale(height[10], 10) > root_scale(height[5], 5) + root_scale(height[5], 5)


def test_two_thirds_scale_is_not_subadditive_on_certified_splits():
    height = _heights()
    assert two_thirds_scale(height[4]) > two_thirds_scale(height[2]) + two_thirds_scale(height[2])
    assert two_thirds_scale(height[8]) > two_thirds_scale(height[3]) + two_thirds_scale(height[5])
