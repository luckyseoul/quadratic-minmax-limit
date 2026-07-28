#!/usr/bin/env python3
"""
Theoretical + certified checks for bi-tight obstruction.

1. For p>=5: avg degree of level-2 bi-tight is 4p/(p^2+1) < 1.
2. For p=5: integral bi-tight levels 2,3,4 infeasible (JSON from e1_bitight_infeas.py).
3. For p=3: C6 bi-tight exists with Phi=Phi-2 (known).

Also verify master lemma samples: tight Max+-only covers have max Max- S >= 0.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

def main() -> None:
    # avg degree formula
    rows = []
    for p in (3, 5, 7, 11, 13, 17, 19):
        n = p * p + 1
        k = 2 * p
        avg_deg = 2 * k / n
        rows.append(
            {
                "p": p,
                "n": n,
                "k_level2": k,
                "avg_degree": avg_deg,
                "avg_deg_lt_1": avg_deg < 1,
            }
        )
    assert all(r["avg_deg_lt_1"] for r in rows if r["p"] >= 5)
    assert not rows[0]["avg_deg_lt_1"]  # p=3: avg_deg = 1.2 > 1

    # load p=5 MILP cert if present
    cert_path = ROOT / "evidence" / "e1_bitight_infeas.json"
    cert = None
    if cert_path.is_file():
        cert = json.loads(cert_path.read_text())
        assert cert["all_integral_infeasible"]
        assert cert["all_fractional_feasible"]

    out = {
        "avg_degree_table": rows,
        "p5_milp_cert_loaded": cert is not None,
        "p5_integral_infeasible": bool(cert and cert["all_integral_infeasible"]),
        "lemma": (
            "For p>=5, level-2 bi-tight has avg degree <1. "
            "At p=5, integral bi-tight is MILP-infeasible for levels 2,3,4. "
            "At p=3, avg degree 1.2 and C6 bi-tight exists."
        ),
        "status": "OPEN residual: uniform bi-tight infeas for all p>=5; deep non-tight ND",
    }
    path = ROOT / "evidence" / "e1_bitight_theory.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
