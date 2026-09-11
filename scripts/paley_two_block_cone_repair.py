"""Finite cone-ranked repair of a q=13 Paley two-block lift.

Every trial cross block is accepted only after ``paley_two_block_score`` has
enumerated its entire 2^27 projective cube.  The cone is a ranking device for
integral one- and two-entry flips; it supplies neither a fractional signing
nor an asymptotic assertion.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

navier_root = Path(os.environ.get(
    "NAVIER_TECHNIQUES",
    str(Path(__file__).resolve().parents[1] / "navier-stokes-techniques"),
))
sys.path.insert(0, str(navier_root))
from stress_cone.stress_cone import NotInConeError, general_cone_representation


N = 14


def score(exe: Path, block: np.ndarray, directory: Path) -> dict:
    path = directory / "block.txt"
    np.savetxt(path, block, fmt="%d")
    return json.loads(subprocess.check_output([str(exe), str(path)], text=True))


def spins(row: dict) -> tuple[np.ndarray, np.ndarray]:
    x = np.ones(N, dtype=np.int8)
    x[1:] = [1 if (row["xmask"] >> i) & 1 else -1 for i in range(N - 1)]
    y = np.array([1 if (row["ymask"] >> i) & 1 else -1 for i in range(N)], dtype=np.int8)
    return x, y


def active_constraints(row: dict) -> tuple[np.ndarray, np.ndarray]:
    directions = []
    defect = []
    for state in row["active"]:
        x, y = spins(state)
        directions.append((x, y, int(state["energy"])))
        defect.append(int(state["value"]) - (int(row["phi"]) - 2))
    return np.asarray(defect, dtype=float), directions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scorer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=731)
    parser.add_argument("--stages", type=int, default=8)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("refusing to overwrite evidence")
    rng = np.random.default_rng(args.seed)
    block = rng.choice(np.array([-1, 1], dtype=np.int8), size=(N, N))
    records: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="paley-cross-repair.") as raw:
        directory = Path(raw)
        current = score(args.scorer, block, directory)
        for stage in range(args.stages):
            # Each active state gives an exact flip derivative.  The cone
            # representation asks whether those simultaneous linearized
            # deficits can be covered; it only ranks integral proposals.
            defect, states = active_constraints(current)
            directions = np.empty((N * N, len(states)), dtype=float)
            for column, (x, y, energy) in enumerate(states):
                delta = -2 * block.reshape(-1) * np.outer(x, y).reshape(-1)
                directions[:, column] = -np.sign(energy) * delta
            try:
                rep = general_cone_representation(directions, defect)
                ranked = np.argsort(-rep.coefficients, kind="stable")[:16]
                cone = "fractional ranking only"
            except NotInConeError:
                ranked = np.argsort(-directions.sum(axis=1), kind="stable")[:16]
                cone = "fallback descent ranking"
            candidates = [(int(e),) for e in ranked]
            candidates += [(int(a), int(b)) for k, a in enumerate(ranked) for b in ranked[k + 1 :]]
            best, selected = current, None
            for edges in candidates:
                trial = block.copy()
                for edge in edges:
                    trial.reshape(-1)[edge] *= -1
                result = score(args.scorer, trial, directory)
                if result["phi"] < best["phi"]:
                    best, selected = result, edges
            records.append({"stage": stage, "before": current, "after": best,
                            "selected_edges": selected, "cone": cone,
                            "evaluated_candidates": len(candidates)})
            if selected is None:
                break
            for edge in selected:
                block.reshape(-1)[edge] *= -1
            current = best
    args.output.write_text(json.dumps({
        "classification": "finite exact two-block cone-repair diagnostic; no all-orders claim",
        "seed": args.seed, "stages": records, "final": current,
        "final_block": block.astype(int).tolist(),
    }, indent=2) + "\n")
    print(json.dumps({"initial_phi": records[0]["before"]["phi"], "final_phi": current["phi"], "stages": len(records)}))


if __name__ == "__main__":
    main()
