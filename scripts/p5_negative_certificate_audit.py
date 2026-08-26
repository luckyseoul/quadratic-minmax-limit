#!/usr/bin/env python3
"""Audit the 33-class p=5 negative two-point certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_negative_profile_cpsat import count_candidates  # noqa: E402
from p5_negative_symmetry_certificate import CASES, validate_case_cover  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(certificate_path: Path) -> dict:
    certificate = json.loads(certificate_path.read_text())
    assert validate_case_cover()
    assert len(count_candidates()) == certificate["profile_count"] == 24
    assert len(CASES) == certificate["placement_orbit_count"] == 33
    assert certificate["coverage_validated"] is True
    assert certificate["status"] == "complete_all_infeasible"
    assert certificate["completed_count"] == 33
    assert certificate["unknown_count"] == certificate["feasible_count"] == 0
    assert len(certificate["rows"]) == 33
    by_index = {row["case_index"]: row for row in certificate["rows"]}
    assert set(by_index) == set(range(33))
    for index, record in enumerate(CASES):
        pp, np, x, y, orbit, pe, ne = record
        row = by_index[index]
        assert (
            row["positive_profile"],
            row["negative_profile"],
            row["positive_parallel_baseline"],
            row["negative_parallel_baseline"],
            row["placement_orbit"],
            row["positive_exception"],
            row["negative_exception"],
        ) == (pp, np, x, y, orbit, pe, ne)
        assert row["solver_status"] == "INFEASIBLE"
        assert row["feasible"] is False
        assert row["constraints"] == {
            "baseline_K": False,
            "pair_coefficients": False,
        }
    return {
        "experiment": "p5_negative_certificate_audit",
        "proved": True,
        "arithmetic_profile_count": 24,
        "placement_orbit_count": 33,
        "infeasible_count": 33,
        "unknown_count": 0,
        "feasible_count": 0,
        "certificate_sha256": sha256(certificate_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.certificate)
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
