#!/usr/bin/env python3
"""W2 for t=-2: candidate p≡1 or 7 (mod 8). Discriminate at 31,47."""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from w2_split_involution import eval_mat  # noqa: E402


def run_one(p):
    r = eval_mat(p, 1, 0, p - 2, p - 1)
    r["name"] = f"p={p} t=-2"
    r["p_mod_8"] = p % 8
    return r


def main():
    primes = (7, 17, 23, 31, 41, 47)
    print("n", len(primes), flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(run_one, p) for p in primes]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(
                {
                    "p": r["p"],
                    "mod8": r["p_mod_8"],
                    "inU": r["inU_y"],
                    "W2": r["W2"],
                    "eigen": r["eigen_minus"],
                },
                flush=True,
            )
    dest = ROOT / "evidence" / "w2_t_m2.json"
    dest.write_text(json.dumps(rows, indent=2))
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
