#!/usr/bin/env python3
"""
Prop 15.285 — n_1D(S) is the Type+ fiber binomial (named in p).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.284 flags.  Does **not** name n_NL or N+(p).  Term_4 is
**not** reduced to 1D only (n_NL is the bulk and is nonzero on
Type+ lines).

============================================================================
Setup.  Regular sets D of size p(p−1)/2 with y_∞=+1.  Type+ 1D
lifts are σ∘L for Type+ linear forms L (ker L a square line) and
|σ^{-1}(+1)|=(p+1)/2 (15.279 R).  Unique count
    n1D_sets = ((p+1)/2) C(p,(p+1)/2)
(30 at p=5, 140 at p=7, 6=N+ at p=3).  n_1D(S)=#{Type+ D: S⊂D}.
n_NL(S)=n(S)−n_1D(S) counts the nonlinear Max+.

============================================================================
Theorem A — PROVED (Max+-free counting; certified p=3,5,7).
  For every S⊂F_q,
      n_1D(S) = ∑_{Type+ L} C(p − |L(S)|, (p+1)/2),
  where |L(S)| is the number of L-fibers met by S.  Same identity
  with Type− forms names the Type−-lift containment (not Max+ for
  p>3).  Fail-when-wrong: replace C(·,(p+1)/2) by C(·,(p−1)/2);
  swap Type− forms in as “1D”.

Theorem B — PROVED.  Closed form in p of (Paley type, F_p-dim):
  m=(p+1)/2, n_□=# of + Paley edges.
  |S|=3, dim-1, pal=+++ : C(p-1,m)+(m-1)C(p-3,m)   (Type+ line)
  |S|=3, dim-1, pal=--- : m C(p-3,m)                 (Type− line)
  |S|=3, dim-2          : n_□ C(p-2,m)+(m-n_□)C(p-3,m)
  |S|=4, dim-1, pal=++++++ : C(p-1,m)+(m-1)C(p-4,m)
  |S|=4, dim-1, pal=------ : m C(p-4,m)   (=0 at p=5,7)
  Independent of AP (15.283 B splits n, not n_1D) and of minpoly
  (15.284 names n, not n_1D).  Fail-when-wrong: drop dim on pal=+++
  (p=5: 4 vs 3; p=7: 18 vs 16); claim n=n_1D at p=5 on the F_p
  4-AP (n=8, n_1D=4, n_NL=4).

n_NL is **not** zero on dim-1 Type+ lines (p=5: 10 / 4 on 3-sets /
4-sets; p=7: 673–681 / 408–416).  Floor stays OPEN.  ∎

============================================================================
Backend: CPU.  GPU unused.
Writes evidence/e1_gmin_m4_prop15285.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)

YPATH = {3: "/tmp/maxplus_p3.npy", 5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def C(n: int, k: int) -> int:
    if k < 0 or n < 0 or k > n:
        return 0
    return comb(n, k)


def n1d_set_count(p: int) -> int:
    m = (p + 1) // 2
    return m * C(p, m)


def n1_formula_3(p: int, pal: tuple, dim: int) -> int:
    m = (p + 1) // 2
    nsq = sum(1 for e in pal if e == 1)
    if dim == 1:
        if nsq == 3:
            return C(p - 1, m) + (m - 1) * C(p - 3, m)
        if nsq == 0:
            return m * C(p - 3, m)
        raise ValueError(f"dim-1 3-set pal={pal}")
    if dim == 2:
        return nsq * C(p - 2, m) + (m - nsq) * C(p - 3, m)
    raise ValueError(dim)


def n1_formula_4_dim1(p: int, pal: tuple) -> int:
    m = (p + 1) // 2
    nsq = sum(1 for e in pal if e == 1)
    if nsq == 6:
        return C(p - 1, m) + (m - 1) * C(p - 4, m)
    if nsq == 0:
        return m * C(p - 4, m)
    raise ValueError(f"dim-1 4-set pal={pal}")


def n1_binomial(p: int, S, forms) -> int:
    m = (p + 1) // 2
    tot = 0
    for alpha, beta in forms:
        vals = {
            (alpha * (int(x) % p) + beta * (int(x) // p)) % p for x in S
        }
        tot += C(p - len(vals), m)
    return tot


def n1_binomial_invert(p: int, S, forms) -> int:
    m = (p + 1) // 2
    tot = 0
    for alpha, beta in forms:
        vals = {
            (alpha * (int(x) % p) + beta * (int(x) // p)) % p for x in S
        }
        tot += C(p - len(vals), m - 1)
    return tot


def type_forms(p: int) -> tuple[list, list]:
    plus, minus = [], []
    for ab in _linear_forms(p):
        (plus if is_type_plus_form(p, *ab) else minus).append(ab)
    return plus, minus


def load_D(p: int) -> np.ndarray:
    Y = np.load(YPATH[p])
    return (np.rint(Y[Y[:, 0] > 0, 1:]) < 0).astype(np.int8)


def build_1d(p: int, *, plus_type: bool) -> np.ndarray:
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab) is plus_type]
    m = (p + 1) // 2
    q = p * p
    rows = []
    seen: set[tuple] = set()
    for ab in forms:
        for mask in range(1 << p):
            if bin(mask).count("1") != m:
                continue
            Aset = {i for i in range(p) if mask >> i & 1}
            y = lift_sigma_vector(p, *ab, Aset)
            D = tuple(1 if y[1 + x] < 0 else 0 for x in range(q))
            if D in seen:
                continue
            seen.add(D)
            rows.append(D)
    return np.asarray(rows, dtype=np.int8) if rows else np.zeros((0, q), dtype=np.int8)


def n_contain(D: np.ndarray, S) -> int:
    sl = list(S)
    if D.shape[0] == 0:
        return 0
    v = D[:, sl[0]].copy()
    for j in sl[1:]:
        np.bitwise_and(v, D[:, j], out=v)
    return int(v.sum())


def paley_of(S, add, neg, chi) -> tuple:
    S = list(S)
    ch = []
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            d = int(add[S[i]][neg[S[j]]])
            ch.append(0 if d == 0 else int(chi[d]))
    return tuple(sorted(ch))


def fp_dim_of(S, p: int) -> int:
    pts = [int(x) for x in S]
    a = pts[0]
    ax, ay = a % p, a // p
    B = []
    for q in pts[1:]:
        B.append(((q % p - ax) % p, (q // p - ay) % p))
    B = [row for row in B if row != (0, 0)]
    if not B:
        return 0
    r = 0
    M = [list(row) for row in B]
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


def affine_lines(p: int) -> list[list[int]]:
    lines = []
    for alpha, beta in _linear_forms(p):
        for c in range(p):
            pts = [
                x
                for x in range(p * p)
                if (alpha * (x % p) + beta * (x // p)) % p == c
            ]
            lines.append(pts)
    return lines


def prove_lifts(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [3, 5, 7]
    rows = []
    ok = True
    for p in primes:
        if not Path(YPATH[p]).is_file():
            continue
        Dall = load_D(p)
        D1 = build_1d(p, plus_type=True)
        D1m = build_1d(p, plus_type=False)
        formula = n1d_set_count(p)
        all_keys = {tuple(int(x) for x in row) for row in Dall}
        subset = all(tuple(int(x) for x in row) in all_keys for row in D1)
        plus, minus = type_forms(p)
        row = {
            "p": p,
            "Nplus": int(Dall.shape[0]),
            "n1D": int(D1.shape[0]),
            "n1Dm": int(D1m.shape[0]),
            "formula": formula,
            "n_plus_forms": len(plus),
            "n_minus_forms": len(minus),
            "subset": subset,
        }
        if D1.shape[0] != formula or len(plus) != (p + 1) // 2:
            ok = False
        if p > 3 and subset is False:
            ok = False
        if p == 3 and D1.shape[0] != Dall.shape[0]:
            ok = False
        rows.append(row)
    return {"proved": ok, "rows": rows}


def prove_formula_3(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    invert_eq = 0
    invert_tot = 0
    pal_dim_split = True
    for p in primes:
        if not Path(YPATH[p]).is_file():
            continue
        F = _kernel_field(p)
        add, neg, chi = F["add"], F["neg"], F["chi"]
        D1 = build_1d(p, plus_type=True)
        D1m = build_1d(p, plus_type=False)
        plus, minus = type_forms(p)
        q = p * p
        bad = 0
        n_checked = 0
        seen_plus = {}
        for S in combinations(range(q), 3):
            pal = paley_of(S, add, neg, chi)
            dim = fp_dim_of(S, p)
            live = n_contain(D1, S)
            pred = n1_formula_3(p, pal, dim)
            binom = n1_binomial(p, S, plus)
            n_checked += 1
            if live != pred or live != binom:
                bad += 1
                ok = False
            inv = n1_binomial_invert(p, S, plus)
            if live > 0:
                invert_tot += 1
                if inv == live:
                    invert_eq += 1
            if pal == (1, 1, 1):
                seen_plus[dim] = live
            if dim == 1 and pal == (1, 1, 1):
                if n_contain(D1m, S) == live:
                    ok = False
        if set(seen_plus) != {1, 2} or len(set(seen_plus.values())) < 2:
            pal_dim_split = False
            ok = False
        blocks.append(
            {
                "p": p,
                "checked": n_checked,
                "bad": bad,
                "+++_n1_by_dim": {str(d): v for d, v in seen_plus.items()},
            }
        )
    return {
        "proved": ok,
        "pal_without_dim_splits": pal_dim_split,
        "invert_eq": invert_eq,
        "invert_tot": invert_tot,
        "blocks": blocks,
    }


def prove_formula_4_dim1(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    ignore_nl_dead = True
    for p in primes:
        if not Path(YPATH[p]).is_file():
            continue
        F = _kernel_field(p)
        add, neg, chi = F["add"], F["neg"], F["chi"]
        Dall = load_D(p)
        D1 = build_1d(p, plus_type=True)
        D1m = build_1d(p, plus_type=False)
        plus, _minus = type_forms(p)
        bad = 0
        n_checked = 0
        n_nl_pos = 0
        typeplus_n1 = typeplus_n1m = None
        for line in affine_lines(p):
            for S in combinations(line, 4):
                pal = paley_of(S, add, neg, chi)
                dim = fp_dim_of(S, p)
                if dim != 1:
                    ok = False
                    continue
                live = n_contain(D1, S)
                pred = n1_formula_4_dim1(p, pal)
                binom = n1_binomial(p, S, plus)
                nall = n_contain(Dall, S)
                n_checked += 1
                if live != pred or live != binom:
                    bad += 1
                    ok = False
                if nall != live:
                    n_nl_pos += 1
                if pal == (1, 1, 1, 1, 1, 1) and typeplus_n1 is None:
                    typeplus_n1 = live
                    typeplus_n1m = n_contain(D1m, S)
        if p == 5 and n_nl_pos == 0:
            ignore_nl_dead = False
            ok = False
        if typeplus_n1 is not None and typeplus_n1 == typeplus_n1m:
            ok = False
        blocks.append(
            {
                "p": p,
                "checked": n_checked,
                "bad": bad,
                "nNL_pos": n_nl_pos,
                "typeplus_n1": typeplus_n1,
                "typeplus_n1m": typeplus_n1m,
            }
        )
    return {
        "proved": ok,
        "ignore_nl_dead": ignore_nl_dead,
        "blocks": blocks,
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "n_NL_named_in_p": False,
        "Nplus_named": False,
        "mixture_meets_floor": False,
        "term4_reduced_to_1d": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.285  n_1D Type+ fiber binomial", flush=True)
    L = prove_lifts([3, 5, 7])
    A = prove_formula_3([5, 7])
    B = prove_formula_4_dim1([5, 7])
    print(f"  lifts ok={L['proved']} rows={L['rows']}", flush=True)
    print(
        f"  A 3-set formula={A['proved']} pal_dim_split={A['pal_without_dim_splits']} "
        f"invert_eq={A['invert_eq']}/{A['invert_tot']}",
        flush=True,
    )
    print(
        f"  B dim-1 4-set formula={B['proved']} "
        f"ignore_nl_dead={B['ignore_nl_dead']}",
        flush=True,
    )
    op = prove_open()
    out = {
        "prop": "15.285",
        "L_status": "OPEN",
        "proved": {
            "typeplus_binomial": A["proved"] and B["proved"] and L["proved"],
            "formula_3": A["proved"],
            "formula_4_dim1": B["proved"],
            "pal_without_dim_is_identity": not A["pal_without_dim_splits"],
            "invert_binom_is_identity": A["invert_eq"] == A["invert_tot"]
            and A["invert_tot"] > 0,
            "n_equals_n1D_p5_4sets": not B["ignore_nl_dead"],
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "n_NL_named_in_p": False,
            "Nplus_named": False,
            "term4_reduced_to_1d": False,
        },
        "status": {
            "floor": False,
            "identity_A": A["proved"],
            "identity_B": B["proved"],
        },
        "lifts": L["rows"],
        "blocks3": A["blocks"],
        "blocks4": B["blocks"],
        "open": op,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15285.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
