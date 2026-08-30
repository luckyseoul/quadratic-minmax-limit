#!/usr/bin/env python3
"""Bi-tight chain entry using the valid 15.720 degree congruence.

Levels 2 and 3 are exactly the bi-tight alternatives used by the E(1)
no-descent reductions.  Level 4 is only a bi-tight corollary and must not be
confused with the one-sided tight level-4 residual.  The retracted 15.167
spectral arrow is not imported.  This module does not soft-close L.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15720 import (  # noqa: E402
    bitight_level_obstruction,
    required_bitight_levels_empty_all_primes,
    theorem_required_bitight_levels,
)


def run_bitight_chain(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19, 23]
    per_p = {
        str(p): {
            "level_2": bitight_level_obstruction(p, 2),
            "level_3": bitight_level_obstruction(p, 3),
            "level_4_bitight_corollary": bitight_level_obstruction(p, 4),
        }
        for p in primes
    }
    general = theorem_required_bitight_levels()
    bitight_ok = bool(
        required_bitight_levels_empty_all_primes()
        and all(
            row["level_2"]["bi_tight_empty"]
            and row["level_3"]["bi_tight_empty"]
            for row in per_p.values()
        )
    )
    out = {
        "title": "Bi-tight chain (15.720 degree congruence)",
        "required_levels": [2, 3],
        "level_4_one_sided_tight_closed": False,
        "bi_tight_required_levels_empty_for_all_p_ge_5": bitight_ok,
        "bi_tight_all_levels_empty": False,
        "residual_required": False,
        "residual_closed_general": False,
        "path": (
            "15.272/15.207 ker(Gsum)=scheme+cross ⇒ pair-degree congruence "
            "mod (p²−1)/2 ⇒ handshake contradiction at required levels 2 and 3; "
            "bi-tight level 4 is only a corollary"
        ),
        "per_p": per_p,
        "general": general,
        "L_status": "OPEN",  # E(1)/Main not wired here
        "prop15167_retracted": True,
    }
    return out


def main() -> dict:
    out = run_bitight_chain()
    path = ROOT / "evidence" / "e1_bitight_chain.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("=== bi-tight chain ===")
    print(
        "bi_tight_required_levels_empty_for_all_p_ge_5="
        f"{out['bi_tight_required_levels_empty_for_all_p_ge_5']}"
    )
    print(f"residual_required={out['residual_required']}")
    print(f"residual_closed_general={out['residual_closed_general']}")
    print(f"L_status={out['L_status']}")
    for p, r in out["per_p"].items():
        print(
            f"  p={p} level2={r['level_2']['bi_tight_empty']} "
            f"level3={r['level_3']['bi_tight_empty']}"
        )
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
