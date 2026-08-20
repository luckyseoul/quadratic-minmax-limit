#!/usr/bin/env python3
"""Exact coefficient, endpoint, and QVAR sieve for the p=29, k=6 stratum."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k6_coefficient_sieve_fast import scan_prime_fast  # noqa: E402


def main() -> dict:
    report = scan_prime_fast(29, minimum=14)
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
