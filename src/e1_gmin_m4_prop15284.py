#!/usr/bin/env python3
"""
Prop 15.284 — dim-2 3-set Aut_∞ type is (Paley, minpoly of ρ=v/u).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.283 flags.  Does **not** name the integer n(S) in p, so
the floor stays OPEN.

============================================================================
Setup.  Aut_∞ key of a 3-set is (u,v) after trans/×□/Frob (15.283).
ρ=v u^{-1}∈F_q.  dim_{F_p}span{u,v}=2 ⇔ ρ∉F_p.  Galois
    Tr(ρ)=(2 ρ_0 + ia ρ_1) mod p,   N(ρ)=ρ^{p+1}
are the coefficients of the minpoly X²−Tr X+N (shipped field
_kernel_field / _make_field).

============================================================================
Theorem A — PROVED at p=5 (6 dim-2 orbits) and p=7 (10 dim-2 orbits).
  n(S) is a function of (Paley triple, Tr(ρ), N(ρ)).
  pal+Tr alone SPLITS at p=7 (two pairs).
  pal+N alone SPLITS at p=5 and p=7.
  Digit-sum of the index is not the Galois trace (coincided at p=5;
  not used).
  Fail-when-wrong: pal-only; pal+Tr; pal+N; invert minpoly (swap Tr,N).
  Values n(pal,Tr,N) are a census table, not a formula in p.  ∎

============================================================================
Backend: CPU, C(q,3) is small.  GPU unused.
Writes evidence/e1_gmin_m4_prop15284.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _load_D(p: int) -> np.ndarray:
    Y = np.load(YPATH[p])
    return (np.rint(Y[Y[:, 0] > 0, 1:]) < 0).astype(np.int8)


def _pow_row(x: int, e: int, mul: list) -> int:
    if x == 0:
        return 0 if e > 0 else 1
    acc, b, ee = 1, int(x), int(e)
    while ee:
        if ee & 1:
            acc = mul[acc][b]
        b = mul[b][b]
        ee >>= 1
    return acc


def _fp_dim(pts, p: int) -> int:
    B = [[int(x) % p, int(x) // p] for x in pts if int(x)]
    if not B:
        return 0
    r = 0
    for c in range(2):
        piv = next((i for i in range(r, len(B)) if B[i][c] % p), None)
        if piv is None:
            continue
        B[r], B[piv] = B[piv], B[r]
        inv = pow(B[r][c] % p, p - 2, p)
        B[r] = [(B[r][j] * inv) % p for j in range(2)]
        for i in range(len(B)):
            if i == r:
                continue
            f = B[i][c] % p
            if f:
                B[i] = [(B[i][j] - f * B[r][j]) % p for j in range(2)]
        r += 1
    return r


def _paley3(S, add, neg, chi) -> tuple:
    a, b, c = (int(x) for x in S)

    def ch(x, y):
        d = int(add[x][neg[y]])
        return 0 if d == 0 else int(chi[d])

    return tuple(sorted((ch(a, b), ch(a, c), ch(b, c))))


def _aut_key3(S, add, mul, neg, frob, squares) -> tuple:
    pts = [int(x) for x in S]
    best = None
    for src in pts:
        t = int(neg[src])
        raw = [int(add[x][t]) for x in pts if x != src]
        for s in squares:
            scaled = tuple(sorted(int(mul[u][s]) for u in raw))
            if best is None or scaled < best:
                best = scaled
            fr = tuple(sorted(int(frob[int(mul[u][s])]) for u in raw))
            if fr < best:
                best = fr
    return best


def _orbit_table(p: int) -> dict:
    F = _kernel_field(p)
    q = p * p
    add, mul, neg, chi = F["add"], F["mul"], F["neg"], F["chi"]
    inv, trt = F["inv"], F["trt"]
    frob = [_pow_row(x, p, mul) for x in range(q)]
    squares = F["squares"]
    D = _load_D(p)
    acc: dict = defaultdict(list)
    for S in combinations(range(q), 3):
        nS = int(D[:, list(S)].min(axis=1).sum())
        key = _aut_key3(S, add, mul, neg, frob, squares)
        pal = _paley3(S, add, neg, chi)
        dim = _fp_dim(key, p)
        Tr = Nr = None
        if key and dim == 2:
            u, v = int(key[0]), int(key[1])
            rho = mul[v][inv[u]]
            Tr = int(trt[rho])
            Nr = int(mul[rho][frob[rho]])
        acc[key].append((nS, pal, dim, Tr, Nr))
    rows = []
    for key, vs in acc.items():
        ns = {t[0] for t in vs}
        rows.append(
            {
                "key": key,
                "n": next(iter(ns)),
                "span": max(ns) - min(ns),
                "pal": vs[0][1],
                "dim": vs[0][2],
                "Tr": vs[0][3],
                "N": vs[0][4],
                "|o|": len(vs),
            }
        )
    return {"p": p, "Nplus": int(D.shape[0]), "rows": rows}


def _split(rows, proj) -> tuple[int, int]:
    g: dict = defaultdict(set)
    for r in rows:
        g[proj(r)].add(r["n"])
    return len(g), sum(1 for s in g.values() if len(s) > 1)


def prove_minpoly_type(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    pal_tr_splits = False
    pal_n_splits = False
    for p in primes:
        tab = _orbit_table(p)
        d2 = [r for r in tab["rows"] if r["dim"] == 2]
        assert d2 and all(r["span"] == 0 for r in tab["rows"])
        n_ptn, s_ptn = _split(d2, lambda r: (r["pal"], r["Tr"], r["N"]))
        n_pt, s_pt = _split(d2, lambda r: (r["pal"], r["Tr"]))
        n_pn, s_pn = _split(d2, lambda r: (r["pal"], r["N"]))
        n_p, s_p = _split(d2, lambda r: r["pal"])
        if s_ptn:
            ok = False
        if s_pt:
            pal_tr_splits = True
        if s_pn:
            pal_n_splits = True
        blocks.append(
            {
                "p": p,
                "n_dim2": len(d2),
                "pal_Tr_N_classes": n_ptn,
                "pal_Tr_N_split": s_ptn,
                "pal_Tr_split": s_pt,
                "pal_N_split": s_pn,
                "pal_split": s_p,
                "orbits": [
                    {
                        "n": r["n"],
                        "pal": list(r["pal"]),
                        "Tr": r["Tr"],
                        "N": r["N"],
                        "key": list(r["key"]),
                    }
                    for r in sorted(d2, key=lambda r: (-r["n"], r["key"]))
                ],
            }
        )
    return {
        "proved": ok,
        "pal_Tr_splits": pal_tr_splits,
        "pal_N_splits": pal_n_splits,
        "blocks": blocks,
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "n_values_named_in_p": False,
        "Nplus_named": False,
        "mixture_meets_floor": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.284  dim-2 3-set type = (Paley, minpoly ρ)", flush=True)
    A = prove_minpoly_type([5, 7])
    print(
        f"  A pal+Tr+N determines n={A['proved']} "
        f"pal+Tr splits={A['pal_Tr_splits']} pal+N splits={A['pal_N_splits']}",
        flush=True,
    )
    op = prove_open()
    out = {
        "prop": "15.284",
        "L_status": "OPEN",
        "proved": {
            "minpoly_determines_n": A["proved"],
            "pal_Tr_is_identity": not A["pal_Tr_splits"],
            "pal_N_is_identity": not A["pal_N_splits"],
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "n_values_named_in_p": False,
        },
        "status": {"floor": False, "identity_A": A["proved"]},
        "blocks": A["blocks"],
        "open": op,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15284.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
