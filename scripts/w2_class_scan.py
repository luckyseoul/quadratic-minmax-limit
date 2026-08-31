#!/usr/bin/env python3
"""W2: all det-normalized split involutions (tr=0, det=-1). Independent mats."""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from w2_split_involution import eval_mat  # noqa: E402


def mats(p):
    """One representative per {±} of [[α,β],[γ,-α]] with α²+βγ=1."""
    out = []
    seen = set()
    for a in range(p):
        for b in range(p):
            need = (1 - (a * a) % p) % p
            if b == 0:
                if need != 0:
                    continue
                cs = range(p)
            else:
                cs = [((need) * pow(b, p - 2, p)) % p]
            for c in cs:
                d = (-a) % p
                key = tuple(min((a, b, c, d), ((-a) % p, (-b) % p, (-c) % p, (-d) % p)))
                if key in seen:
                    continue
                seen.add(key)
                out.append(key)
    return out


def run_one(args):
    p, A, B, C, D = args
    r = eval_mat(p, A, B, C, D)
    return r


def scan(p, workers=86):
    ms = mats(p)
    print(f"p={p} n_class={len(ms)} expect~{p*(p+1)//2}", flush=True)
    jobs = [(p, A, B, C, D) for (A, B, C, D) in ms]
    rows = []
    n_inU = 0
    n_w2 = 0
    hits = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        done = 0
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            done += 1
            if r["inU_y"]:
                n_inU += 1
            if r["W2"]:
                n_w2 += 1
                hits.append(r["ABCD"])
                print(f"  HIT {r['ABCD']} wt={r['wt']}", flush=True)
            if done % 50 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} inU={n_inU} W2={n_w2}", flush=True)
    return {
        "p": p,
        "n_class": len(ms),
        "n_eigen": sum(1 for r in rows if r["eigen_minus"]),
        "n_inU": n_inU,
        "n_W2": n_w2,
        "hits": hits[:20],
    }


def main():
    out = {}
    for p in (17, 31):
        out[str(p)] = scan(p)
        print("SUMMARY", out[str(p)], flush=True)
    dest = ROOT / "evidence" / "w2_class_scan.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
