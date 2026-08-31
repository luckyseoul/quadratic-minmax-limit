#!/usr/bin/env python3
"""W2 t=±i (i^2=-1) at p≡1 (mod 4). Independent primes, ProcessPool."""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from w2_split_involution import ab, eval_mat  # noqa: E402


def jobs():
    out = []
    for p in (5, 13, 17, 29, 37, 41):
        a, b = ab(p)
        i = (b * pow(a, p - 2, p)) % p
        out.append((f"p={p} t=i", p, 1, 0, i, p - 1))
        out.append((f"p={p} t=-i", p, 1, 0, (p - i) % p, p - 1))
        out.append((f"p={p} t=-2", p, 1, 0, p - 2, p - 1))
    return out


def run_one(job):
    name, p, A, B, C, D = job
    r = eval_mat(p, A, B, C, D)
    r["name"] = name
    return r


def main():
    js = jobs()
    print(f"n={len(js)}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=18) as ex:
        futs = [ex.submit(run_one, j) for j in js]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(
                {
                    "name": r["name"],
                    "eigen": r["eigen_minus"],
                    "inU": r["inU_y"],
                    "W2": r["W2"],
                    "ABCD": r["ABCD"],
                },
                flush=True,
            )
    dest = ROOT / "evidence" / "w2_t_i.json"
    dest.write_text(json.dumps(rows, indent=2))
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
