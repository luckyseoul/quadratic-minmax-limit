#!/usr/bin/env python3
"""W1 p≡1 mod 8: test named d = -g (least primitive root) and friends."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def least_primroot(p):
    fac = []
    m = p - 1
    d = 2
    mm = m
    while d * d <= mm:
        if mm % d == 0:
            fac.append(d)
            while mm % d == 0:
                mm //= d
        d += 1
    if mm > 1:
        fac.append(mm)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    return None


def qr_sym_odd(p, d):
    m = (p - 1) // 2
    S = set(range(m + 1))
    Sd = {(x + d) % p for x in S}
    delta = S.symmetric_difference(Sd)
    nqr = sum(1 for x in delta if x != 0 and pow(x, (p - 1) // 2, p) == 1)
    return nqr % 2


def worker(p):
    g = least_primroot(p)
    m = (p - 1) // 2
    cands = {
        "-g": (p - g) % p,
        "-2g": (p - 2 * g) % p,
        "g": g % p,
        "-5": (p - 5) % p,
        "-3": (p - 3) % p,
        "i": next(x for x in range(p) if (x * x) % p == p - 1),
    }
    cands["-i"] = (p - cands["i"]) % p
    upper = set(range(m + 1, p))
    rec = {"p": p, "g": g, "in_upper": {}, "eps": {}}
    for name, d in cands.items():
        rec["in_upper"][name] = d in upper
        rec["eps"][name] = qr_sym_odd(p, d) if d else None
    # least d in upper with eps=1
    rec["least_upper_hit"] = next(d for d in range(m + 1, p) if qr_sym_odd(p, d) == 1)
    rec["n_hits"] = sum(1 for d in range(m + 1, p) if qr_sym_odd(p, d) == 1)
    return rec


def main():
    primes = [p for p in range(17, 200) if is_prime(p) and p % 8 == 1]
    with ProcessPoolExecutor(max_workers=min(len(primes), 40)) as ex:
        futs = [ex.submit(worker, p) for p in primes]
        for f in futs:
            r = f.result()
            print(
                f"p={r['p']} g={r['g']} least_hit={r['least_upper_hit']} "
                f"n_hits={r['n_hits']} -g_up={r['in_upper']['-g']} -g_eps={r['eps']['-g']} "
                f"-3={r['eps']['-3']} -5={r['eps']['-5']} -i={r['eps']['-i']}",
                flush=True,
            )


if __name__ == "__main__":
    main()
