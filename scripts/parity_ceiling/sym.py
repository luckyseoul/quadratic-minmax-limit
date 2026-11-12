"""Best Boolean states of the Paley conference matrix C_{q+1} fixed by an automorphism of prime order l.

PSL(2,q) acts on R^{q+1} by signed permutations that preserve x^T C x:
    g v_P = lambda_P(g) v_{gP},  v_inf = (1,0), v_j = (j,1),   (g.x)(gP) = chi(lambda_P(g)) x(P).
Every x whose stabiliser in PSL(2,q) (mod global sign) is nontrivial is fixed, up to sign, by an element of
prime order l, and all subgroups of order l are conjugate.  So for each prime l dividing (q-1)/2 or (q+1)/2
the script builds the orbit quotient R (m x m, m = number of orbits) and maximises it:
  exact Gray-code enumeration when m <= EXACT_LIMIT, annealing otherwise (lower bound only).
The unipotent case l = q has two orbits and value q; it is omitted.

    gcc -O2 -o qopt scripts/parity_ceiling/qopt.c -lm      # build once, from the repo root
    python3 scripts/parity_ceiling/sym.py 61 73 [--json out.json]
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

from paley_lib import chi_table, is_prime, paley, prim_root

EXACT_LIMIT = 34
QOPT = os.environ.get("QOPT", "./qopt")


def signed_perm(g, q, chi):
    (a, b), (c, d) = g
    n = q + 1
    pi = np.zeros(n, dtype=np.int64)
    s = np.zeros(n, dtype=np.int64)

    def img(u, v):
        if v % q != 0:
            return 1 + (u * pow(v, q - 2, q)) % q, v % q
        return 0, u % q

    P, lam = img(a, c)
    pi[0], s[0] = P, chi[lam]
    for j in range(q):
        P, lam = img(a * j + b, c * j + d)
        pi[1 + j], s[1 + j] = P, chi[lam]
    return pi, s


def _mul(A, B, q):
    return (
        ((A[0][0] * B[0][0] + A[0][1] * B[1][0]) % q, (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % q),
        ((A[1][0] * B[0][0] + A[1][1] * B[1][0]) % q, (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % q),
    )


def _pow(A, e, q):
    R = ((1, 0), (0, 1))
    while e:
        if e & 1:
            R = _mul(R, A, q)
        A = _mul(A, A, q)
        e >>= 1
    return R


def element_of_order(q, l, chi):
    g0 = prim_root(q)
    if (q - 1) // 2 % l == 0:  # split torus: j -> t^2 j
        t = pow(g0, (q - 1) // (2 * l), q)
        return ((t, 0), (0, pow(t, q - 2, q)))
    assert (q + 1) // 2 % l == 0
    for tau in range(q):  # companion matrix of x^2 - tau x + 1, irreducible iff tau^2 - 4 is a nonsquare
        if chi[(tau * tau - 4) % q] == -1:
            h = _pow(((0, q - 1), (1, tau)), (q + 1) // l, q)
            if h not in (((1, 0), (0, 1)), ((q - 1, 0), (0, q - 1))):
                return h
    raise ValueError((q, l))


def quotient(C, pi, s):
    """Orbit vectors U (rows) of the Boolean solutions of g.x = eps x, and R = U C U^T; None if there are none."""
    n = len(pi)
    for eps in (1, -1):
        seen = np.zeros(n, bool)
        U = []
        ok = True
        for P in range(n):
            if seen[P]:
                continue
            u = np.zeros(n, dtype=np.int64)
            cur, val = P, 1
            while True:
                seen[cur] = True
                u[cur] = val
                val = val * eps * s[cur]
                cur = pi[cur]
                if cur == P:
                    break
            if val != 1:
                ok = False
                break
            U.append(u)
        if ok:
            U = np.array(U)
            return eps, U, U @ C @ U.T
    return None


def best_invariant(q, l, C, chi):
    pi, s = signed_perm(element_of_order(q, l, chi), q, chi)
    assert (C[np.ix_(pi, pi)] * np.outer(s, s) == C).all(), "not an automorphism"
    res = quotient(C, pi, s)
    if res is None:
        return None
    _, U, R = res
    m = len(U)
    const = int(np.trace(R)) // 2
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"{m}\n")
        for row in R:
            f.write(" ".join(map(str, row.tolist())) + "\n")
        fn = f.name
    try:
        if m <= EXACT_LIMIT:
            out = subprocess.run([QOPT, "exact", fn], capture_output=True, text=True, check=True).stdout.split()
            d = dict(t.split("=") for t in out if "=" in t)
            return {"l": l, "orbits": m, "mode": "exact", "max": int(d["max"]) + const, "min": int(d["min"]) + const}
        out = subprocess.run([QOPT, "sa", fn, "4000", "400", "7"], capture_output=True, text=True, check=True).stdout
        mx = int(out.split("max>=")[1].split()[0])
        mn = int(out.split("min<=")[1].split()[0])
        return {"l": l, "orbits": m, "mode": "anneal", "max": mx + const, "min": mn + const}
    finally:
        os.unlink(fn)


def symmetric_best(q):
    C, chi = paley(q), chi_table(q)
    ls = sorted({p for p in range(2, q) if is_prime(p) and ((q - 1) // 2 % p == 0 or (q + 1) // 2 % p == 0)})
    rows = [r for r in (best_invariant(q, l, C, chi) for l in ls) if r is not None]
    best = max(max(r["max"], -r["min"]) for r in rows)
    return {"q": q, "n": q + 1, "best": best, "all_exact": all(r["mode"] == "exact" for r in rows), "by_order": rows}


if __name__ == "__main__":
    qs = [int(a) for a in sys.argv[1:] if a.isdigit()]
    out = []
    for q in qs:
        r = symmetric_best(q)
        out.append(r)
        S = 0.5 * (q + 1) * math.sqrt(q)
        for row in r["by_order"]:
            print(f"q={q} l={row['l']:3d} orbits={row['orbits']:3d} {row['mode']:6s} max={row['max']} min={row['min']}", flush=True)
        print(f"q={q}: best state with a nontrivial stabiliser = {r['best']} ({r['best'] / S:.4f})"
              f"{'  [all quotients exact]' if r['all_exact'] else ''}\n", flush=True)
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump(out, f, indent=1)
