#!/usr/bin/env python3
"""Piecewise linear stay box on each (2/p)_8 class. ProcessPool."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from w1_octic import extra_cands, oct8, sfree
from w1_residual_ab import eps_from_pref, qr_prefix, residuals

ROOT = Path(__file__).resolve().parents[1]


def worker(row):
    p, a, b = row
    m = (p - 1) // 2
    pref = qr_prefix(p)
    table = [eps_from_pref(p, d, pref) for d in range(p)]
    i = (b * pow(a, p - 2, p)) % p
    o = oct8(p)
    ag = a if a % 4 == 1 else -a
    hits = []
    for u in range(-4, 5):
        for v in range(-4, 5):
            for w in range(-4, 5):
                for k in range(-8, 9):
                    d = (u * a + v * b + w * i + k) % p
                    if d >= m + 1 and table[d] == 1:
                        hits.append((u, v, w, k))
    return {
        "p": p,
        "a": a,
        "b": b,
        "ag": ag,
        "oct8": o,
        "sfree": sfree(a),
        "hits": hits,
        "n_hits": len(hits),
    }


def main():
    rows_in = residuals(20000, strict=True)
    print(f"n={len(rows_in)}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=86) as ex:
        futs = [ex.submit(worker, r) for r in rows_in]
        for f in as_completed(futs):
            rows.append(f.result())
    rows.sort(key=lambda r: r["p"])
    plus = [r for r in rows if r["oct8"] == 1]
    minus = [r for r in rows if r["oct8"] == -1]
    print(f"oct+ {len(plus)} oct- {len(minus)}", flush=True)

    def inter(group):
        sets = [set(r["hits"]) for r in group]
        return set.intersection(*sets) if sets else set()

    ip = inter(plus)
    im = inter(minus)
    print(f"box ∩ oct+ n={len(ip)} sample={sorted(ip)[:12]}", flush=True)
    print(f"box ∩ oct- n={len(im)} sample={sorted(im)[:12]}", flush=True)

    # also split by ag%8
    by = {}
    for r in rows:
        by.setdefault((r["oct8"], r["ag"] % 8), []).append(r)
    print("---- oct8, ag%8 box intersections ----", flush=True)
    splits = {}
    for key, g in sorted(by.items()):
        s = inter(g)
        splits[str(key)] = {"n_primes": len(g), "n_forms": len(s), "sample": sorted(s)[:8], "ps": [x["p"] for x in g[:6]]}
        print(f"  {key} n_p={len(g):2d} n_forms={len(s):4d} ps={splits[str(key)]['ps']}", flush=True)

    dest = ROOT / "evidence" / "w1_octic_box.json"
    dest.write_text(
        json.dumps(
            {
                "n_plus": len(plus),
                "n_minus": len(minus),
                "n_inter_plus": len(ip),
                "n_inter_minus": len(im),
                "inter_plus": sorted(ip)[:40],
                "inter_minus": sorted(im)[:40],
                "splits": splits,
            },
            indent=2,
        )
    )
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
