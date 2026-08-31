#!/usr/bin/env python3
"""W1 residual: class function of (2/p)_8 / Gaussian π=a+8ci. ProcessPool."""
from __future__ import annotations

import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from w1_residual_ab import (
    chi,
    eps_from_pref,
    named,
    qr_prefix,
    residuals,
)

ROOT = Path(__file__).resolve().parents[1]


def oct8(p: int) -> int:
    r = pow(2, (p - 1) // 8, p)
    if r == 1:
        return 1
    if r == p - 1:
        return -1
    raise ValueError(f"2^((p-1)/8) not ±1 at p={p}: {r}")


def sfree(n: int) -> int:
    n = abs(n)
    r = 1
    x = n
    p = 2
    while p * p <= x:
        c = 0
        while x % p == 0:
            x //= p
            c += 1
        if c % 2:
            r *= p
        p = p + 1 if p == 2 else p + 2
    if x > 1:
        r *= x
    return r


def extra_cands(p, a, b, i):
    ag = a if a % 4 == 1 else -a
    c = b // 8
    sf = sfree(a)
    half = pow(2, p - 2, p)
    o = oct8(p)
    out = {
        "c": c % p,
        "-c": (-c) % p,
        "ag": ag % p,
        "-ag": (-ag) % p,
        "sf": sf % p,
        "-sf": (-sf) % p,
        "ag+8c": (ag + b) % p,
        "ag-8c": (ag - b) % p,
        "-ag-8c": (-ag - b) % p,
        "o*(-ag)": (-o * ag) % p,
        "o*(-b)": (-o * b) % p,
        "o*(-sf)": (-o * sf) % p,
        "o*a/2": (o * a * half) % p,
        "ag+o": (ag + o) % p,
        "-ag+o": (-ag + o) % p,
        "c*o": (c * o) % p,
        "-c*o": (-c * o) % p,
    }
    return out


def worker(row):
    p, a, b = row
    m = (p - 1) // 2
    pref = qr_prefix(p)
    o = oct8(p)
    i = (b * pow(a, p - 2, p)) % p
    cands = named(p, a, b)
    cands.update(extra_cands(p, a, b, i))
    stay = {}
    for k, v in cands.items():
        v %= p
        if v == 0:
            stay[k] = None
            continue
        up = m + 1 <= v <= p - 1
        stay[k] = eps_from_pref(p, v, pref) if up else None
    ag = a if a % 4 == 1 else -a
    return {
        "p": p,
        "a": a,
        "b": b,
        "c": b // 8,
        "ag": ag,
        "ag_mod_8": ag % 8,
        "a_mod_8": a % 8,
        "c_mod_2": (b // 8) % 2,
        "sfree_a": sfree(a),
        "oct8": o,
        "leg2_ag": 1 if ((ag * ag - 1) // 8) % 2 == 0 else -1,
        "stay": stay,
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
    keys = sorted(set().union(*(r["stay"].keys() for r in rows)))

    def split(pred):
        return [r for r in rows if pred(r)]

    plus = split(lambda r: r["oct8"] == 1)
    minus = split(lambda r: r["oct8"] == -1)
    print(f"oct8=+ {len(plus)} first={[r['p'] for r in plus[:8]]}", flush=True)
    print(f"oct8=- {len(minus)} first={[r['p'] for r in minus[:8]]}", flush=True)

    print("---- stay ε=1 ALWAYS on oct=+ / oct=- (upper only) ----", flush=True)
    always_plus = []
    always_minus = []
    for k in keys:
        def tally(group):
            ones = [r["p"] for r in group if r["stay"].get(k) == 1]
            zeros = [r["p"] for r in group if r["stay"].get(k) == 0]
            none = [r["p"] for r in group if r["stay"].get(k) is None]
            return ones, zeros, none

        op, zp, np_ = tally(plus)
        om, zm, nm = tally(minus)
        tagp = "ALWAYS" if op and not zp and not np_ else (
            "NEVER" if zp and not op and not np_ else ("LOWER" if np_ == plus else "MIXED")
        )
        tagm = "ALWAYS" if om and not zm and not nm else (
            "NEVER" if zm and not om and not nm else ("LOWER" if nm == minus else "MIXED")
        )
        if tagp == "ALWAYS":
            always_plus.append(k)
        if tagm == "ALWAYS":
            always_minus.append(k)
        if tagp == "ALWAYS" or tagm == "ALWAYS":
            print(f"  {k:16s} oct+:{tagp:7s} n1={len(op):2d} n0={len(zp):2d}  oct-:{tagm:7s} n1={len(om):2d} n0={len(zm):2d}", flush=True)

    print("ALWAYS oct+", always_plus, flush=True)
    print("ALWAYS oct-", always_minus, flush=True)

    # pairwise class function: (d+, d-) named
    pairs = []
    for kp in keys:
        for km in keys:
            ok = True
            for r in plus:
                if r["stay"].get(kp) != 1:
                    ok = False
                    break
            if not ok:
                continue
            for r in minus:
                if r["stay"].get(km) != 1:
                    ok = False
                    break
            if ok:
                pairs.append((kp, km))
    print(f"PAIR always n={len(pairs)} sample={pairs[:20]}", flush=True)

    # a mod 8 after ag, vs oct8 and ε(-ag)
    print("---- first 16: p a ag ag%8 oct8 stay(-ag) stay(-a) stay(a/2) ----", flush=True)
    for r in rows[:16]:
        print(
            r["p"],
            r["a"],
            r["ag"],
            r["ag_mod_8"],
            r["oct8"],
            r["stay"].get("-ag"),
            r["stay"].get("-a"),
            r["stay"].get("a/2"),
            r["stay"].get("o*(-ag)"),
            flush=True,
        )

    dest = ROOT / "evidence" / "w1_octic.json"
    dest.write_text(
        json.dumps(
            {
                "n": len(rows),
                "n_oct_plus": len(plus),
                "n_oct_minus": len(minus),
                "always_plus": always_plus,
                "always_minus": always_minus,
                "pairs": pairs[:50],
                "n_pairs": len(pairs),
                "rows": [
                    {
                        "p": r["p"],
                        "a": r["a"],
                        "ag": r["ag"],
                        "c": r["c"],
                        "oct8": r["oct8"],
                        "leg2_ag": r["leg2_ag"],
                        "stay": {k: v for k, v in r["stay"].items() if v is not None},
                    }
                    for r in rows
                ],
            },
            indent=2,
        )
    )
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
