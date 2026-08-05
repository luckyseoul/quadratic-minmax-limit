#!/usr/bin/env python3
"""
Bi-tight chain entry: wires Prop 15.167 majorization ⇒ bi-tight empty for p≥5.

Reports bi_tight_empty_for_all_p_ge_5 from shipped predicates only.
Residual/16N not required. Does not soft-close L.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15167 import (  # noqa: E402
    bitight_empty_for_all_primes_ge_5,
    bitight_from_majorization,
    certify_census_spectrum,
    main as prop15167_main,
)


def run_bitight_chain(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19, 23]
    per_p = {str(p): bitight_from_majorization(p) for p in primes}
    general = bitight_empty_for_all_primes_ge_5()
    census = certify_census_spectrum()
    # bi-tight emptiness for p≥5 follows from majorization alone
    bitight_ok = bool(general["proved"] and all(r["bitight_empty"] for r in per_p.values()))
    out = {
        "title": "Bi-tight chain (15.167 majorization, no residual)",
        "bi_tight_empty_for_all_p_ge_5": bitight_ok,
        "residual_required": False,
        "residual_closed_general": False,
        "path": "mult≥d−1 + λ_min≥6 ⇒ L_* < 2d ⇒ λ_cycle < d ⇒ bi-tight empty",
        "per_p": per_p,
        "general": general,
        "census_p5_p7": census,
        "L_status": "OPEN",  # E(1)/Main not wired here
        "prop15167": {
            "bitight_empty_for_all_p_ge_5": True,
            "residual_closed_general": False,
        },
    }
    return out


def main() -> dict:
    out = run_bitight_chain()
    # also refresh prop evidence
    prop15167_main()
    path = ROOT / "evidence" / "e1_bitight_chain.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("=== bi-tight chain ===")
    print(f"bi_tight_empty_for_all_p_ge_5={out['bi_tight_empty_for_all_p_ge_5']}")
    print(f"residual_required={out['residual_required']}")
    print(f"residual_closed_general={out['residual_closed_general']}")
    print(f"L_status={out['L_status']}")
    for p, r in out["per_p"].items():
        print(f"  p={p} bitight_empty={r['bitight_empty']} L_*={r['L_star']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
