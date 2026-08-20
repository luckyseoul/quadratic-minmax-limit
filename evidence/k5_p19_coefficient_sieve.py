#!/usr/bin/env python3
"""Exact coefficient, endpoint, and QVAR sieve for the p=19, k=5 stratum."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_p23_coefficient_sieve import scan_prime  # noqa: E402


def main() -> dict:
    report = scan_prime(19)
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
