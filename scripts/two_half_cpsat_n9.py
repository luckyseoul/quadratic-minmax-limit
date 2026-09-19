#!/usr/bin/env python3
"""CP-SAT min_R B(A,R) for an exact order-9 minimizer (m_9=12).

Uses the existing complete finite minimax in original_mo_two_half_geometry.
ILS supplies A and a hint orientation; CP-SAT must return OPTIMAL.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from original_mo_two_half_geometry import (  # noqa: E402
    analyze_orientation,
    bits_from_skew,
    solve_cpsat,
    skew_from_bits,
)
from two_half_extend_n9_n10 import ils_orientation, ils_signing  # noqa: E402

OUT = ROOT / "evidence" / "two_half_n9_n10_20260919"
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    phi, a = ils_signing(9, 12, seed=0)
    print(f"A Phi={phi}", flush=True)
    assert phi == 12
    np.save(OUT / "A9_exact.npy", a)
    print("short ILS orientation (80 steps, 1 worker) ...", flush=True)
    b, r = ils_orientation(a, seed=0, steps=80)
    hint = bits_from_skew(r)
    print(f"ILS hint B={b} bits={hint}", flush=True)

    print("CP-SAT n=9 workers=86 ...", flush=True)
    solved = solve_cpsat(a, 12, hint, workers=86, max_time=600.0)
    target = 2 * math.sqrt(2) * 12
    rec = {
        "n": 9,
        "M": 12,
        "target_2sqrt2_M": target,
        "cpsat": solved,
        "zero_error_passes": solved["objective_B"] is not None
        and solved["objective_B"] <= target + 1e-9,
        "normalized_excess": None
        if solved["objective_B"] is None
        else (solved["objective_B"] - target) / (9 ** 1.5),
        "A": a.tolist(),
        "seconds": round(time.time() - t0, 3),
    }
    (OUT / "n9_cpsat.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in rec if k != "A"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
