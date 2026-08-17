#!/usr/bin/env python3
"""
Prop 15.489 — Leftover p=7 ridges are 2p(χ*μ) for small mean-zero μ.

F(p) = {0}
      ∪ {2p(−m+p 1_S): |S|=m}
      ∪ {2p bal(a(t−t0)): a≠0}          (p≥5)
      ∪ {2p (χ_p * μ): μ:F_p→{1−m..m−1}, ∑μ=0, ∑μ²∈2{1..m}}.

At p=3,5 the CSP count is 6,130.  At p=7, |F|=638 and GPU CSP
N=5726; dropping the χ-family leaves |F|=78 and N=1316≠5726.
Q_τ still unnamed (D=|H+|/(2p) is this CSP, not a short Jacobi
closed form).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.488 flags.
"""
from __future__ import annotations

import json
from itertools import combinations, product
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15488 import csp_count, targets

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15489.json"


def balanced(t: int, p: int) -> int:
    t %= p
    if t > p // 2:
        t -= p
    return t


def chi_p(t: int, p: int) -> int:
    t %= p
    if t == 0:
        return 0
    return 1 if pow(t, (p - 1) // 2, p) == 1 else -1


def conv_chi(mu, p: int) -> tuple:
    twop = 2 * p
    return tuple(
        twop
        * sum(mu[c] * chi_p(t - c, p) for c in range(p))
        for t in range(p)
    )


def F_named(
    p: int,
    drop_level: bool = False,
    drop_affine: bool = False,
    drop_chi: bool = False,
):
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
    if p >= 5 and not drop_chi:
        lo, hi = 1 - m, m - 1
        allowed_ss = {2 * k for k in range(1, m + 1)}
        for mu in product(range(lo, hi + 1), repeat=p):
            if sum(mu) != 0:
                continue
            if sum(x * x for x in mu) not in allowed_ss:
                continue
            if not any(mu):
                continue
            F.add(conv_chi(mu, p))
    return F


def prove_A() -> dict:
    """p=3,5 CSP with the enlarged F still hits 6,130. Fail flip T1."""
    n3 = csp_count(3, F_named(3))
    n5 = csp_count(5, F_named(5))
    n5f = csp_count(5, F_named(5), flip_T1=True)
    return {
        "proved": n3 == 6 and n5 == 130 and n5f != 130,
        "N3": n3,
        "N5": n5,
        "nF5": len(F_named(5)),
        "N5_flip": n5f,
    }


def prove_B() -> dict:
    """p=7: χ*μ family is necessary. |F|=638; drop χ ⇒ 78."""
    n_full = len(F_named(7))
    n_drop = len(F_named(7, drop_chi=True))
    return {
        "proved": n_full == 638 and n_drop == 78,
        "nF": n_full,
        "nF_drop_chi": n_drop,
        "gpu_N_full": 5726,
        "gpu_N_drop_chi": 1316,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": "F is named as a generating set; D=|H+|/(2p) is the "
        "ridge CSP, not a short Jacobi value. Floor OPEN.",
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
