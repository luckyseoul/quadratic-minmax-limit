#!/usr/bin/env python3
"""W1 residual p=a^2+64c^2: named d from (a,b,i). ProcessPool over primes."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    if n % 3 == 0:
        return n == 3
    d = 5
    while d * d <= n:
        if n % d == 0 or n % (d + 2) == 0:
            return False
        d += 6
    return True


def residuals(pmax: int, strict: bool = True):
    """p=a^2+64c^2. strict: residual W1, p≡1 or 49 (mod 120)."""
    out = []
    c = 1
    while 64 * c * c + 1 < pmax:
        a = 1
        while True:
            p = a * a + 64 * c * c
            if p > pmax:
                break
            if is_prime(p) and (not strict or p % 120 in (1, 49)):
                out.append((p, a, 8 * c))
            a += 2
        c += 1
    out.sort()
    return out


def chi(p, x):
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def qr_prefix(p):
    """pref[n] = #QR in {1,...,n}."""
    pref = [0] * p
    acc = 0
    for n in range(1, p):
        if chi(p, n) == 1:
            acc += 1
        pref[n] = acc
    return pref


def qr_interval(pref, lo, hi, p):
    """#QR in [lo, hi] ∩ {1,...,p-1}, 0-based inclusive, wrapping not allowed."""
    lo = max(lo, 1)
    hi = min(hi, p - 1)
    if hi < lo:
        return 0
    return pref[hi] - pref[lo - 1]


def eps_from_pref(p, d, pref):
    """ε(d)=|QR ∩ (S Δ (S+d))| mod 2, S={0,...,m}."""
    d %= p
    if d == 0:
        return 0
    m = (p - 1) // 2
    nqr = 0
    if d <= m:
        nqr += qr_interval(pref, 1, d - 1, p)
        nqr += qr_interval(pref, m + 1, d + m, p)
    else:
        k = p - d
        nqr += qr_interval(pref, m - k + 1, m, p)
        nqr += qr_interval(pref, p - k, p - 1, p)
    return nqr % 2


def eps_d(p, d):
    return eps_from_pref(p, d, qr_prefix(p))


def ab_from(p, a, b):
    # unique odd a>0, even b>0 already; flip a sign so a≡1 (mod 4) when wanted
    a1 = a if a % 4 == 1 else p - a  # not in F_p; a is integer <sqrt(p)
    a1 = a if a % 4 == 1 else -a
    return a, b, a1


def named(p, a, b):
    inv = lambda t: pow(t % p, p - 2, p)
    i = (b * inv(a)) % p  # i^2 ≡ -1
    half = inv(2)
    cands = {
        "a": a % p,
        "-a": (-a) % p,
        "b": b % p,
        "-b": (-b) % p,
        "a+b": (a + b) % p,
        "a-b": (a - b) % p,
        "b-a": (b - a) % p,
        "(a+b)/2": ((a + b) * half) % p,
        "(a-b)/2": ((a - b) * half) % p,
        "a/2": (a * half) % p,
        "b/2": (b * half) % p,
        "b/4": (b * inv(4)) % p,
        "b/8": (b * inv(8)) % p,
        "2a": (2 * a) % p,
        "2b": (2 * b) % p,
        "4a": (4 * a) % p,
        "i": i,
        "-i": (-i) % p,
        "2i": (2 * i) % p,
        "i+1": (i + 1) % p,
        "i-1": (i - 1) % p,
        "ia": (i * a) % p,
        "ib": (i * b) % p,
        "a+i": (a + i) % p,
        "b+i": (b + i) % p,
        "a*b": (a * b) % p,
        "(p-1)/16": ((p - 1) // 16) % p if p % 16 == 1 else None,
        "(p-1)/24": ((p - 1) // 24) % p,
        "(p-1)/40": ((p - 1) // 40) % p if p % 40 == 1 else None,
        "(p-1)/48": ((p - 1) // 48) % p if p % 48 == 1 else None,
        "(p-1)/80": ((p - 1) // 80) % p if p % 80 == 1 else None,
        "(p-1)/120": ((p - 1) // 120) % p if p % 120 == 1 else None,
        "eighth": ((p - 1) // 8) % p,
        "-eighth": (-((p - 1) // 8)) % p,
    }
    # Gauss a ≡ 1 mod 4 representative as integer in {±a}
    ag = a if a % 4 == 1 else -a
    cands["ag"] = ag % p
    cands["-ag"] = (-ag) % p
    cands["ag+b"] = (ag + b) % p
    cands["ag-b"] = (ag - b) % p
    return {k: v for k, v in cands.items() if v is not None}


def worker(row):
    p, a, b = row
    m = (p - 1) // 2
    pref = qr_prefix(p)
    table = [eps_from_pref(p, d, pref) for d in range(p)]
    cands = named(p, a, b)
    eps = {k: table[v % p] for k, v in cands.items()}
    upper = {k: bool(m + 1 <= (cands[k] % p) <= p - 1) for k in cands}
    stay_eps = {k: (eps[k] if upper[k] else None) for k in cands}
    i = cands["i"]
    hits = [d for d in range(m + 1, p) if table[d] == 1]
    forms = []
    for u in range(-4, 5):
        for v in range(-4, 5):
            for w in range(-4, 5):
                for k in range(-8, 9):
                    d = (u * a + v * b + w * i + k) % p
                    if d >= m + 1 and table[d] == 1:
                        forms.append((u, v, w, k, d))
    return {
        "p": p,
        "a": a,
        "b": b,
        "p_mod_120": p % 120,
        "i": i,
        "i2": (i * i) % p,
        "eps": eps,
        "upper": upper,
        "stay_eps": stay_eps,
        "cands": cands,
        "n_hits": len(hits),
        "forms": forms,
    }


def main():
    rows = residuals(20000, strict=True)
    print(f"n_residual={len(rows)} first={rows[:12]}", flush=True)
    out_rows = []
    with ProcessPoolExecutor(max_workers=86) as ex:
        futs = [ex.submit(worker, r) for r in rows]
        for f in as_completed(futs):
            out_rows.append(f.result())
    out_rows.sort(key=lambda r: r["p"])
    # which named keys have ε=1 on EVERY residual prime
    keys = sorted(set().union(*(r["eps"].keys() for r in out_rows)))
    print("primes", [r["p"] for r in out_rows], flush=True)
    print("---- named STAY ε (upper half only) ----", flush=True)
    always = []
    never = []
    mixed = []
    for k in keys:
        ones = [r["p"] for r in out_rows if r["stay_eps"].get(k) == 1]
        zeros = [r["p"] for r in out_rows if r["stay_eps"].get(k) == 0]
        none = [r["p"] for r in out_rows if r["stay_eps"].get(k) is None]
        if none:
            tag = "LOWER"
        elif not zeros:
            tag = "ALWAYS"
            always.append(k)
        elif not ones:
            tag = "NEVER"
            never.append(k)
        else:
            tag = "MIXED"
            mixed.append((k, ones[:6], zeros[:6]))
        print(
            f"  {k:16s} {tag:7s} n1={len(ones):2d} n0={len(zeros):2d} lower={len(none):2d}",
            flush=True,
        )
    print("ALWAYS", always, flush=True)
    print("NEVER", never, flush=True)
    for r in out_rows[:6]:
        print(
            f"brute p={r['p']} a={r['a']} b={r['b']} hits={r['n_hits']} forms={len(r['forms'])}",
            flush=True,
        )
    # intersection of (u,v,w,k) that hit ε=1 at every residual prime
    keys_per = [set((u, v, w, k) for u, v, w, k, _d in r["forms"]) for r in out_rows]
    common = set.intersection(*keys_per) if keys_per else set()
    global_ok = [{"form": f"{u}*a+{v}*b+{w}*i+{k}", "uvwk": [u, v, w, k]} for (u, v, w, k) in sorted(common)]
    print("GLOBAL_LINEAR_ALWAYS n=", len(global_ok), global_ok[:30], flush=True)
    dest = ROOT / "evidence" / "w1_residual_ab.json"
    dest.write_text(
        json.dumps(
            {
                "always": always,
                "never": never,
                "mixed": mixed,
                "global_linear": global_ok,
                "rows": [
                    {
                        "p": r["p"],
                        "a": r["a"],
                        "b": r["b"],
                        "eps": r["eps"],
                    }
                    for r in out_rows
                ],
            },
            indent=2,
        )
    )
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
