"""Arithmetic ceiling for Phi(C), C any symmetric conference matrix of order n = q+1.

For x in {+-1}^n put h_i = x_i (C x)_i.  The note proves:
  (A) every h_i is odd and all h_i are congruent mod 4;
  (B) sum h_i^2 = q n;       (C) sum h_i = 2 Q_C(x);
  (R) row identity: with D = diag(x) C diag(x),  D h = q 1.  Hence for each i there is a
      set M_i of (q - h_i)/2 indices j != i with  sum_{j in M_i} h_j = (sum h - h_i - q)/2.

closed_form(q)  : max over the two classes of min_{kappa} n (q + kappa^2 - 4) / (4 kappa), rounded
                  down to an odd integer.  Uses (A)-(C) and (h-c)(h-c-4) >= 0 only.
ceiling(q)      : largest Q for which an integer field multiset satisfying (A), (B), (C), (R) exists.
Both are upper bounds for Phi(C) for EVERY symmetric conference matrix of order q+1.

    python3 scripts/parity_ceiling/ceiling.py            # table
    python3 scripts/parity_ceiling/ceiling.py --json out.json
"""
from __future__ import annotations

import json
import math
import sys


def closed_form(q: int) -> int:
    n = q + 1
    best = 0.0
    for r in (1, 3):
        vals = [n * (q + k * k - 4) / (4 * k) for k in range(r, 4 * int(math.isqrt(q)) + 8, 4)]
        best = max(best, min(vals))
    u = math.floor(best + 1e-9)
    return u if u % 2 == 1 else u - 1


def multisets(n: int, T1: int, T2: int) -> list[dict[int, int]]:
    """All multisets {value: count} of nonzero integers, at most n of them, with sum T1 and sum of squares T2."""
    out: list[dict[int, int]] = []
    if T2 < 0 or (T2 - T1) % 2:
        return out
    vals: list[int] = []
    for a in range(math.isqrt(T2), 0, -1):
        vals += [a, -a]

    def rec(k: int, rem: int, s1: int, s2: int, cur: dict[int, int]) -> None:
        if s2 == 0:
            if s1 == 0:
                out.append(dict(cur))
            return
        if k == len(vals) or abs(s1) > s2:
            return
        v = vals[k]
        for c in range(min(rem, s2 // (v * v)), -1, -1):
            if c:
                cur[v] = c
            rec(k + 1, rem - c, s1 - c * v, s2 - c * v * v, cur)
            if c:
                del cur[v]

    rec(0, n, T1, T2, {})
    return out


def row_identity_ok(levels: dict[int, int], q: int) -> bool:
    """Necessary condition (R) for a field multiset {field value: count}."""
    n = sum(levels.values())
    S = sum(v * c for v, c in levels.items())
    items = sorted(levels.items())
    for a, _ in items:
        if (q - a) % 2 or (S - a - q) % 2:
            return False
        need_k, need_s = (q - a) // 2, (S - a - q) // 2
        if not 0 <= need_k <= n - 1:
            return False
        reach: dict[int, set[int]] = {0: {0}}
        for v, c in items:
            cap = c - (1 if v == a else 0)
            new: dict[int, set[int]] = {}
            for k, sums in reach.items():
                for m in range(cap + 1):
                    if k + m > need_k:
                        break
                    tgt = new.setdefault(k + m, set())
                    add = m * v
                    for s in sums:
                        tgt.add(s + add)
            reach = new
        if need_s not in reach.get(need_k, ()):
            return False
    return True


def levels_by_class(q: int, depth: int = 1, use_row_identity: bool = True):
    """For each residue class r in {1,3}: the top `depth` admissible values of Q with their field multisets."""
    n = q + 1
    res = []
    for r in (1, 3):
        c = r
        while (c + 4) ** 2 <= q:
            c += 4
        assert (q - c * c) % 4 == 0
        mp = (q - c * c) // 4
        found = 0
        T1 = (n * mp) // (2 * c + 4) + 2
        while found < depth and T1 > -n:
            num = n * mp - 2 * c * T1  # = 4 T2
            if num % 4 == 0:
                T2 = num // 4
                if T2 >= abs(T1) and (T2 - T1) % 2 == 0:
                    good = []
                    for m in multisets(n, T1, T2):
                        lev = {c + 4 * t: cnt for t, cnt in m.items()}
                        rest = n - sum(m.values())
                        if rest:
                            lev[c] = rest
                        if not use_row_identity or row_identity_ok(lev, q):
                            good.append(dict(sorted(lev.items())))
                    if good:
                        res.append({"Q": (n * c + 4 * T1) // 2, "class": r, "T1": T1, "T2": T2, "multisets": good})
                        found += 1
            T1 -= 1
    res.sort(key=lambda z: -z["Q"])
    return res


def ceiling(q: int, use_row_identity: bool = True) -> int:
    return levels_by_class(q, 1, use_row_identity)[0]["Q"]


TABLE_Q = [5, 9, 13, 17, 25, 29, 37, 41, 49, 53, 61, 73, 89, 97, 101, 109, 113, 173, 229, 293]

if __name__ == "__main__":
    rows = []
    for q in TABLE_Q:
        n = q + 1
        top = levels_by_class(q, 1, True)[0]
        rows.append(
            {
                "q": q,
                "n": n,
                "spectral": 0.5 * n * math.sqrt(q),
                "closed_form": closed_form(q),
                "integral_ceiling": ceiling(q, False),
                "ceiling": top["Q"],
                "class": top["class"],
                "example_fields": {str(k): v for k, v in top["multisets"][0].items()},
                "n_multisets": len(top["multisets"]),
            }
        )
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump(rows, f, indent=1)
    print("  q    n  spectral  closed  integral  ceiling  example field multiset at the ceiling")
    for r in rows:
        print(
            f"{r['q']:4d} {r['n']:4d}  {r['spectral']:8.2f}  {r['closed_form']:6d}  {r['integral_ceiling']:8d}  {r['ceiling']:7d}  {r['example_fields']}"
        )
