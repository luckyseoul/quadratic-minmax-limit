#!/usr/bin/env python3
"""
Prop 15.288 — Type+ parent |U|-histogram of a dim-2 3-set is named in (p,e).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.287 flags.  Does **not** name n_H or n_NL (T-weighted).
Floor stays OPEN.

============================================================================
Setup.  U=S\\D_A, D_A={L∉σ}, |σ|=m=(p+1)/2.  |U|=# of points of S
with L(x)∈σ.  e=# of Paley edges of S.  Dim-2 ⇒ the three
differences have distinct F_p-directions, and L collapses a pair
iff that edge is Paley (χ_q|F_p^*=1, so χ_q of a direction is
χ_q of any generator).  There are e collapsing Type+ forms and
m−e forms with three distinct fibers.

============================================================================
Theorem A — PROVED (Max+-free counting; certified p=5,7, sample p=11).
  For a dim-2 3-set with e Paley edges,
    n_0 = (m−e) C(p−3,m) + e C(p−2,m)           (= n_1D, 15.285)
    n_1 = 3(m−e) C(p−3,m−1) + e C(p−2,m−1)
    n_2 = 3(m−e) C(p−3,m−2) + e C(p−2,m−1)
    n_3 = (m−e) C(p−3,m−3) + e C(p−2,m−2)
  n_k = #{Type+ (L,σ) : |U|=k}.  Sum_k n_k = m C(p,m) = n1D_sets.
  Fail-when-wrong: swap C(p−3, m−k) ↔ C(p−3, m−k−1); apply the
  dim-2 formula to a dim-1 3-set; Type− forms in place of Type+.  ∎

Theorem B — OPEN.  Λ_H = ∑_k n_k τ_k with τ_k = mean Hoffman T
  given |U|=k is **not** a function of (p,e) (p=5 two e=2 orbits
  have different T-weights).  n_H / n_NL unnamed.  Floor open.  ∎

============================================================================
Backend: CPU.  GPU unused.
Writes evidence/e1_gmin_m4_prop15288.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15285 import n1_formula_3, n1d_set_count  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def C(n: int, k: int) -> int:
    if k < 0 or n < 0 or k > n:
        return 0
    return comb(n, k)


def nU_formula(p: int, e: int) -> tuple[int, int, int, int]:
    m = (p + 1) // 2
    n0 = (m - e) * C(p - 3, m) + e * C(p - 2, m)
    n1 = 3 * (m - e) * C(p - 3, m - 1) + e * C(p - 2, m - 1)
    n2 = 3 * (m - e) * C(p - 3, m - 2) + e * C(p - 2, m - 1)
    n3 = (m - e) * C(p - 3, m - 3) + e * C(p - 2, m - 2)
    return n0, n1, n2, n3


def nU_formula_inverted(p: int, e: int) -> tuple[int, int, int, int]:
    """Wrong index on the 3-fiber binomials."""
    m = (p + 1) // 2
    n0 = (m - e) * C(p - 3, m - 1) + e * C(p - 2, m)
    n1 = 3 * (m - e) * C(p - 3, m) + e * C(p - 2, m - 1)
    n2 = 3 * (m - e) * C(p - 3, m - 1) + e * C(p - 2, m - 1)
    n3 = (m - e) * C(p - 3, m - 2) + e * C(p - 2, m - 2)
    return n0, n1, n2, n3


def _eval_L(alpha: int, beta: int, x: int, p: int) -> int:
    return (alpha * (x % p) + beta * (x // p)) % p


def live_nU(p: int, S, *, plus_type: bool = True) -> tuple[int, int, int, int]:
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab) is plus_type]
    counts = [0, 0, 0, 0]
    for alpha, beta in forms:
        fibers = [_eval_L(alpha, beta, x, p) for x in S]
        for mask in range(1 << p):
            if bin(mask).count("1") != m:
                continue
            plus = {i for i in range(p) if mask >> i & 1}
            u = sum(1 for t in fibers if t in plus)
            counts[u] += 1
    return tuple(counts)


def paley_e(S, F) -> int:
    add, neg, chi = F["add"], F["neg"], F["chi"]
    S = list(S)
    e = 0
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            d = add[S[i]][neg[S[j]]]
            if d and chi[d] == 1:
                e += 1
    return e


def fp_dim(S, p: int) -> int:
    pts = [int(x) for x in S]
    a = pts[0]
    B = [((q % p - a % p) % p, (q // p - a // p) % p) for q in pts[1:]]
    B = [row for row in B if row != (0, 0)]
    if not B:
        return 0
    M = [list(row) for row in B]
    r = 0
    for c in range(2):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c] % p, p - 2, p)
        M[r] = [(M[r][j] * inv) % p for j in range(2)]
        for i in range(len(M)):
            if i == r:
                continue
            f = M[i][c] % p
            if f:
                M[i] = [(M[i][j] - f * M[r][j]) % p for j in range(2)]
        r += 1
    return r


def prove_formula(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    invert_hit = 0
    invert_tot = 0
    dim1_wrong = 0
    minus_eq = 0
    minus_tot = 0
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        bad = 0
        checked = 0
        n1_match = 0
        seen_e = set()
        for S in combinations(range(q), 3):
            if fp_dim(S, p) != 2:
                e = paley_e(S, F)
                # dim-1 +++ must not match the dim-2 formula
                if e == 3:
                    pred = nU_formula(p, e)
                    live = live_nU(p, S)
                    if live == pred:
                        dim1_wrong += 1
                        ok = False
                continue
            e = paley_e(S, F)
            pred = nU_formula(p, e)
            live = live_nU(p, S)
            checked += 1
            if live != pred:
                bad += 1
                ok = False
            if pred[0] == n1_formula_3(p, tuple(sorted((1,) * e + (-1,) * (3 - e))), 2):
                n1_match += 1
            inv = nU_formula_inverted(p, e)
            invert_tot += 1
            if inv == live:
                invert_hit += 1
            minus = live_nU(p, S, plus_type=False)
            minus_tot += 1
            if minus == live:
                minus_eq += 1
            seen_e.add(e)
            if sum(pred) != n1d_set_count(p):
                ok = False
        if seen_e != {0, 1, 2, 3}:
            ok = False
        blocks.append(
            {
                "p": p,
                "checked": checked,
                "bad": bad,
                "n1_match": n1_match,
            }
        )
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    if minus_tot == 0 or minus_eq == minus_tot:
        ok = False
    return {
        "proved": ok,
        "blocks": blocks,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "minus_eq": minus_eq,
        "minus_tot": minus_tot,
        "dim1_formula_applied_hits": dim1_wrong,
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "n_H_named_in_p": False,
        "n_NL_named_in_p": False,
        "tau_k_named": False,
        "Nplus_named": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.288  Type+ |U|-histogram named in (p,e)", flush=True)
    A = prove_formula([5, 7])
    print(
        f"  A formula={A['proved']} invert={A['invert_hit']}/{A['invert_tot']} "
        f"minus_eq={A['minus_eq']}/{A['minus_tot']}",
        flush=True,
    )
    op = prove_open()
    out = {
        "prop": "15.288",
        "L_status": "OPEN",
        "proved": {
            "nU_formula": A["proved"],
            "invert_is_identity": A["invert_hit"] == A["invert_tot"],
            "type_minus_is_identity": A["minus_eq"] == A["minus_tot"],
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "n_H_named_in_p": False,
            "n_NL_named_in_p": False,
        },
        "status": {"floor": False, "identity_A": A["proved"]},
        "A": A,
        "open": op,
        "sample_p5": {str(e): list(nU_formula(5, e)) for e in range(4)},
        "sample_p7": {str(e): list(nU_formula(7, e)) for e in range(4)},
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15288.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
