"""A step ceiling compatible with the certified table does not beat 1/2."""
from __future__ import annotations

import json
from pathlib import Path

from minmax_quadratic import step_equilibrium_level


def test_certified_steps_force_the_equilibrium_above_one_half():
    payload = json.loads(
        (Path(__file__).resolve().parents[1] / "evidence" / "exact_m_table.json").read_text()
    )
    worst = None
    for order in range(2, 15):
        step = payload[str(order + 1)]["m"] - payload[str(order)]["m"]
        level = step_equilibrium_level(step, order)
        if worst is None or level > worst:
            worst = level
    assert worst > 0.5
