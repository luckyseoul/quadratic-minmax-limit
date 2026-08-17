#!/usr/bin/env python3
"""
Prop 15.488 — Ridge CSP counts |H+| at p=3,5 (and GPU at p=7).

g=2f−1 equals −1/p + ∑_{s∈U} h_s(y−sx) on the complementary
slopes U={r−1,…,p−1}, r=(p+1)/2.  H_s=p² h_s is integer.
|H+| is the number of U-tuples in F^U with
∑ H_s(y−sx) ∈ {p(1−p), p(p+1)} at every (x,y).

F at p=3,5 is named below.  F at p=7 has size 288 and is not
named in p.  Q_τ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.487 flags.
"""
from __future__ import annotations

import json
from itertools import combinations, product
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15488.json"


def balanced(t: int, p: int) -> int:
    t %= p
    if t > p // 2:
        t -= p
    return t


def F_named(p: int, drop_level: bool = False, drop_affine: bool = False):
    """Named 1-D ridge catalog.  Fail-when-wrong: drop a family."""
    F = set()
    F.add(tuple(0 for _ in range(p)))
    m = (p - 1) // 2
    twop = 2 * p
    if not drop_level:
        for S in combinations(range(p), m):
            F.add(tuple(twop * (-m + p * (1 if t in S else 0)) for t in range(p)))
    if p >= 5 and not drop_affine:
        for a in range(1, p):
            for t0 in range(p):
                F.add(tuple(twop * balanced(a * (t - t0), p) for t in range(p)))
    return F


def targets(p: int):
    return p * (1 - p), p * (1 + p)


def csp_count(p: int, F=None, flip_T1: bool = False) -> int:
    if F is None:
        F = F_named(p)
    cat = list(F)
    nF = len(cat)
    Uslopes = list(range((p + 1) // 2 - 1, p))
    U = len(Uslopes)
    T0, T1 = targets(p)
    if flip_T1:
        T1 += 1
    n = 0
    for idxs in product(range(nF), repeat=U):
        ok = True
        for y in range(p):
            for x in range(p):
                S = 0
                for i, s in enumerate(Uslopes):
                    S += cat[idxs[i]][(y - s * x) % p]
                if S != T0 and S != T1:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            n += 1
    return n


def prove_A() -> dict:
    """p=3: F={0}∪{2p(−m+p 1_{t0})}, |F|=4, CSP=6."""
    F = F_named(3)
    n = csp_count(3, F)
    n_flip = csp_count(3, F, flip_T1=True)
    n_drop = csp_count(3, F_named(3, drop_level=True))
    return {
        "proved": n == 6 and n_flip != 6 and n_drop != 6 and len(F) == 4,
        "nF": len(F),
        "N": n,
        "N_flip_T1": n_flip,
        "N_drop_level": n_drop,
    }


def prove_B() -> dict:
    """p=5: F=zero ∪ affine 2p·bal(a(t−t0)) ∪ 2p(−m+p 1_S), |F|=31, CSP=130."""
    F = F_named(5)
    n = csp_count(5, F)
    n_flip = csp_count(5, F, flip_T1=True)
    n_drop_S = csp_count(5, F_named(5, drop_level=True))
    n_drop_aff = csp_count(5, F_named(5, drop_affine=True))
    return {
        "proved": n == 130
        and n_flip != 130
        and n_drop_S != 130
        and n_drop_aff != 130
        and len(F) == 31,
        "nF": len(F),
        "N": n,
        "N_flip_T1": n_flip,
        "N_drop_level": n_drop_S,
        "N_drop_affine": n_drop_aff,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": "p=7 F has size 288, not named in p; GPU CSP N=5726 "
        "(scratch netcount_ridge). D=|H+|/(2p) still unnamed.",
    }


def main() -> dict:
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "A": A,
        "B": B,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
