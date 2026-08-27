#!/usr/bin/env python3
"""Create a compact hash summary of an independent p=7 size-six replay."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import socket
from pathlib import Path


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value) -> str:
    encoded = json.dumps(value, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--floor", type=Path, required=True)
    parser.add_argument("--orbits", type=Path, required=True)
    parser.add_argument("--ordinary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    floor = json.loads(args.floor.read_text())
    orbits = json.loads(args.orbits.read_text())
    ordinary = json.loads(args.ordinary.read_text())
    out = {
        "experiment": "p7_size6_independent_summary",
        "status": "complete_compact_hash_summary",
        "host": socket.gethostname(),
        "floor": {
            "path": str(args.floor),
            "file_sha256": file_hash(args.floor),
            "backend": floor["backend"],
            "device": floor["device"],
            "checked_boundaries": floor["checked_boundaries"],
            "floor_surviving_boundaries": floor["floor_surviving_boundaries"],
            "floor_rejected_boundaries": floor["floor_rejected_boundaries"],
            "direction_odd_fibre_histogram": floor["direction_odd_fibre_histogram"],
            "survivor_sha256": floor["survivor_sha256"],
        },
        "orbits": {
            "path": str(args.orbits),
            "file_sha256": file_hash(args.orbits),
            "candidate_boundaries": orbits["candidate_boundaries"],
            "orbit_count": orbits["orbit_count"],
            "orbit_size_sum": orbits["orbit_size_sum"],
            "ordered_orbits_canonical_sha256": canonical_hash(orbits["orbits"]),
            "profile_histogram_canonical_sha256": canonical_hash(orbits["profile_histogram"]),
        },
        "ordinary": {
            "path": str(args.ordinary),
            "file_sha256": file_hash(args.ordinary),
            "ordinary_orbits_in_source": ordinary["ordinary_orbits_in_source"],
            "deep_deficit_orbits_in_source": ordinary["deep_deficit_orbits_in_source"],
            "processed_ordinary_orbits": ordinary["processed_ordinary_orbits"],
            "elevation_cases": ordinary["elevation_cases"],
            "modular_infeasible_cases": ordinary["modular_infeasible_cases"],
            "surviving_cases": ordinary["surviving_cases"],
            "floor_pair_counts": ordinary["floor_pair_counts"],
            "catalog_pattern_counts": ordinary["catalog_pattern_counts"],
        },
    }
    temporary = args.output.with_name(f".{args.output.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(out, indent=2) + "\n")
    os.replace(temporary, args.output)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
