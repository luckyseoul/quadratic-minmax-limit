#!/usr/bin/env python3
"""Prop. 15.691 -- fractional-moment no-go for the c=2 free-energy target."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def adversarial_upper_rate(c: float) -> float:
    if c <= 2 * math.sqrt(math.log(2)):
        raise ValueError("optimized interior theta requires c>2*sqrt(log 2)")
    return -(c / 2 - math.sqrt(math.log(2))) ** 2


def sufficient_lower_rate(c: float) -> float:
    return c / 2 - c * c / 4


def fractional_moment_barrier() -> float:
    return math.log(2) / (math.sqrt(math.log(2)) - 0.5)


def theorem_record() -> dict[str, object]:
    upper_two = adversarial_upper_rate(2.0)
    target_two = sufficient_lower_rate(2.0)
    upper_three = adversarial_upper_rate(3.0)
    target_three = sufficient_lower_rate(3.0)
    barrier = fractional_moment_barrier()
    if not upper_two < target_two or not upper_three > target_three:
        raise ArithmeticError("fractional-moment comparison changed")
    return {
        "proposition": "15.691",
        "adversarial_upper_rate": "-(c/2-sqrt(log 2))^2",
        "sufficient_lower_rate": "c/2-c^2/4",
        "c2": {
            "upper_rate": upper_two,
            "required_rate": target_two,
            "old_target_false": True,
        },
        "barrier_c": barrier,
        "all_c_below_barrier_excluded": True,
        "c3": {
            "upper_rate": upper_three,
            "required_rate": target_three,
            "target": "log P_a(tanh(3/sqrt(n))) >= -3n/4-o(n)",
            "excluded_by_fractional_moment": False,
        },
        "closes_limit": False,
        "L_status": "OPEN",
        "proved": True,
    }


def main() -> None:
    target = ROOT / "evidence" / "e1_gmin_m4_prop15691.json"
    target.write_text(json.dumps(theorem_record(), indent=2) + "\n")
    print(f"Prop. 15.691: c=2 free-energy target false; wrote {target}")


if __name__ == "__main__":
    main()
