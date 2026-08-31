#!/usr/bin/env python3
"""W2: p-dependent split involutions (tr=0, det=-1). Independent (p,t)."""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import krylov_g, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def eval_mat(p, A, B, C, D):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    pi = _mobius_perm(p, A, B, C, D)
    inv = np.empty_like(pi)
    inv[pi] = np.arange(len(pi))
    y = np.zeros_like(z)
    for j in range(q + 1):
        src = int(inv[j])
        if j == 0:
            sw = chi(C) if C else 1
            if sw == 0:
                sw = 1
        else:
            lin = add(mul(C, j - 1), D)
            sw = chi(lin)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    Cmat = paley_conference_prime_power(p)
    yy = y.astype(np.float64)
    em = bool(np.max(np.abs(Cmat @ yy + p * yy)) < 1e-6)
    yb = ((1 - y) // 2).astype(np.uint8)
    inU_y = bool(int(yb[0]) == 1 and int(yb[1]) == 0)
    w2 = None
    wt = int(((bits ^ yb) & 1).sum())
    if em and inU_y:
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        gamma, _, _, _ = named_gamma(p)
        _, facs = _g_factors(p)
        d = (bits ^ yb) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        c = krylov_g(wfn, gamma, mul, gen, q, (q - 1) // 2)
        if c is not None:
            cl = list(map(int, c))
            w2 = all(_poly_gcd(cl, f) == [1] for f in facs)
    return {
        "p": p,
        "ABCD": [A, B, C, D],
        "eigen_minus": em,
        "inU_y": inU_y,
        "W2": w2,
        "wt": wt,
    }


def ab(p):
    if p % 4 != 1:
        return None, None
    for a in range(1, p, 2):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = int(round(b2**0.5))
        if b * b == b2 and b % 2 == 0 and b > 0:
            return a, b
    return None, None


def named_jobs(p):
    jobs = []
    # family π(x)=x/(t x - 1) = (1,0,t,-1)
    jobs.append(("t=1 x/(x-1)", p, 1, 0, 1, p - 1))
    m = (p - 1) // 2
    jobs.append(("t=m pole-2", p, 1, 0, m, p - 1))
    jobs.append(("t=2", p, 1, 0, 2, p - 1))
    jobs.append(("t=-2", p, 1, 0, p - 2, p - 1))
    jobs.append(("t=3", p, 1, 0, 3, p - 1))
    # smallest nsq
    nsq = next(t for t in range(2, p) if pow(t, (p - 1) // 2, p) == p - 1)
    jobs.append((f"t=nsq{nsq}", p, 1, 0, nsq, p - 1))
    if p % 4 == 1:
        a, b = ab(p)
        i = (b * pow(a, p - 2, p)) % p
        jobs.append(("t=i", p, 1, 0, i, p - 1))
        jobs.append(("t=-i", p, 1, 0, (p - i) % p, p - 1))
        jobs.append(("t=a", p, 1, 0, a % p, p - 1))
        jobs.append(("t=b", p, 1, 0, b % p, p - 1))
        jobs.append(("t=1+i", p, 1, 0, (1 + i) % p, p - 1))
        jobs.append(("t=2i", p, 1, 0, (2 * i) % p, p - 1))
        # Cayley-like (i, -1; 1, -i) det = -i²+1=2, not -1
        # (i,1; 1, -i): det = -i²-1=0.  (i,1; -1, -i): det=-i²+1=2
        # Want α²+βγ=1, α=i: -1 + βγ=1 ⇒ βγ=2.
        jobs.append(("alpha=i beta=1 gamma=2", p, i, 1, 2, (p - i) % p))
        jobs.append(("alpha=i beta=2 gamma=1", p, i, 2, 1, (p - i) % p))
        jobs.append(("alpha=0 beta=1 gamma=1", p, 0, 1, 1, 0))  # π=1/x  det=-1
        jobs.append(("alpha=0 beta=-1 gamma=1", p, 0, p - 1, 1, 0))  # π=-1/x
    return jobs


def run_one(job):
    name, p, A, B, C, D = job
    r = eval_mat(p, A, B, C, D)
    r["name"] = name
    return r


def main():
    primes = (5, 7, 11, 13, 17, 19, 23)
    jobs = []
    for p in primes:
        jobs.extend(named_jobs(p))
    print(f"n_jobs={len(jobs)}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=86) as ex:
        futs = [ex.submit(run_one, j) for j in jobs]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            print(
                {
                    "p": r["p"],
                    "name": r["name"],
                    "eigen": r["eigen_minus"],
                    "inU": r["inU_y"],
                    "W2": r["W2"],
                },
                flush=True,
            )
    rows.sort(key=lambda r: (r["p"], r["name"]))
    dest = ROOT / "evidence" / "w2_split_involution.json"
    dest.write_text(json.dumps(rows, indent=2))
    # summary: names that hit W2 at every p where inU
    print("---- W2 hits ----", flush=True)
    for r in rows:
        if r["W2"]:
            print(r["p"], r["name"], r["ABCD"], flush=True)


if __name__ == "__main__":
    main()
