#!/usr/bin/env python3
"""
Prop 15.286 — n_H = 1D⊙Hoffman-nc containment; n_NL = n_H at p=5.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.285 flags.  Does **not** name n_NL in p, nor N+(p).
Does **not** close the floor.

============================================================================
Setup.  n(S)=n_1D(S)+n_NL(S) (15.285).  Type+ 1D lifts y_A.  Hoffman
norm-circle T of C'=D_{y_A} C D_{y_A} (15.138–15.141) gives
    D = {finite x : (y_A ⊙ z_T)_x = −1},
the 1D⊙nc family.  n_H(S)=#{such D ∋ S}.  Labeled count
    Λ_H(S)=∑_{Type+ A} #{Hoffman T : S ⊂ D_A Δ T}.
Overcount κ=Λ_H/n_H when n_H>0.

============================================================================
Theorem A — PROVED at p=5 (all Aut_∞ 3-orbits).
  n_NL(S)=n_H(S).  The 100 1D⊙nc sets are exactly the nonlinear
  Max+ (y_∞=+1).  Fail-when-wrong: Type− parents as “1D⊙nc”
  (n_H^− ≠ n_NL on +++ dim-2); claim n_H=n_NL at p=7
  (n_rest>0 on all 10 dim-2 orbits).  ∎

Theorem B — PROVED at p=5 and p=7 on every dim-2 3-orbit.
  κ=Λ_H/n_H is Aut_∞-constant: 3 at p=5, 1 at p=7.
  Fail-when-wrong: κ=1 at p=5; κ=3 at p=7.  The interpolant
  8−p is not a theorem (p=11 would be negative).  ∎

Theorem C — OPEN.  n_H is not a stock binomial / linear form in
  (p, #□-edges, Tr, N).  n_NL is not named in p.  1D⊙nc does not
  exhaust NL at p=7 (1764 vs 5586).  Floor / phi_F_ge_6 stay False.  ∎

============================================================================
Backend: CPU ProcessPool over Type+ parents.  GPU unused.
Writes evidence/e1_gmin_m4_prop15286.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import norm_circle_z_on  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15284 import _aut_key3, _fp_dim, _paley3  # noqa: E402
from e1_gmin_m4_prop15285 import build_1d, n1d_set_count  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _pow(x, e, mul):
    if x == 0:
        return 0 if e > 0 else 1
    acc, b, ee = 1, int(x), int(e)
    while ee:
        if ee & 1:
            acc = mul[acc][b]
        b = mul[b][b]
        ee >>= 1
    return acc


def _plus_sets(p: int):
    m = (p + 1) // 2
    return [
        {i for i in range(p) if mask >> i & 1}
        for mask in range(1 << p)
        if bin(mask).count("1") == m
    ]


def typeplus_y(p: int, *, plus_type: bool = True):
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab) is plus_type]
    seen = {}
    for ab in forms:
        for Aset in _plus_sets(p):
            y = np.asarray(lift_sigma_vector(p, *ab, Aset), dtype=np.float64)
            if y[0] < 0:
                y = -y
            key = tuple(int(np.sign(v)) for v in y)
            if key not in seen:
                seen[key] = y
    return list(seen.values())


def _parent_circles(payload):
    y, C, p, q, mul_t, add_t, neg_t, frob_tab = payload

    def mul(u, v):
        return mul_t[int(u)][int(v)]

    def add(u, v):
        return add_t[int(u)][int(v)]

    def neg(u):
        return neg_t[int(u)]

    def frob(u):
        return frob_tab[int(u)]

    Cp = np.diag(y) @ C @ np.diag(y)
    ncs = norm_circle_z_on(Cp, p, mul, add, neg, frob, q)
    Ds = []
    for item in ncs:
        yy = y * item["z"]
        if yy[0] < 0:
            yy = -yy
        Ds.append(tuple(1 if yy[1 + u] < 0 else 0 for u in range(q)))
    return Ds


def build_1d_odot_nc(p: int, *, plus_type: bool = True) -> np.ndarray:
    C = paley_conference_prime_power(p).astype(np.float64)
    F = _kernel_field(p)
    q = p * p
    mul_t, add_t, neg_t = F["mul"], F["add"], F["neg"]
    frob_tab = [_pow(x, p, mul_t) for x in range(q)]
    parents = typeplus_y(p, plus_type=plus_type)
    payloads = [(y, C, p, q, mul_t, add_t, neg_t, frob_tab) for y in parents]
    uniq: set[tuple] = set()
    W = min(default_workers(), max(1, len(payloads)))
    with ProcessPoolExecutor(max_workers=W) as ex:
        for part in ex.map(_parent_circles, payloads):
            uniq.update(part)
    if not uniq:
        return np.zeros((0, q), dtype=np.int8)
    return np.asarray(sorted(uniq), dtype=np.int8)


def _load_D(p: int) -> np.ndarray:
    Y = np.load(YPATH[p])
    return (np.rint(Y[Y[:, 0] > 0, 1:]) < 0).astype(np.int8)


def _n_contain(D: np.ndarray, S) -> int:
    sl = list(S)
    if D.shape[0] == 0:
        return 0
    v = D[:, sl[0]].copy()
    for j in sl[1:]:
        np.bitwise_and(v, D[:, j], out=v)
    return int(v.sum())


def _aut3_counts(p: int, Dall, D1, DH):
    F = _kernel_field(p)
    q = p * p
    add, mul, neg, chi = F["add"], F["mul"], F["neg"], F["chi"]
    frob = [_pow(x, p, mul) for x in range(q)]
    squares = [r for r in range(1, q) if chi[r] == 1]
    acc: dict = defaultdict(list)
    for S in combinations(range(q), 3):
        key = _aut_key3(S, add, mul, neg, frob, squares)
        acc[key].append(
            (
                _paley3(S, add, neg, chi),
                _fp_dim(key, p),
                _n_contain(Dall, S),
                _n_contain(D1, S),
                _n_contain(DH, S),
            )
        )
    rows = []
    for key, vs in acc.items():
        pal, dim, n, n1, nH = vs[0]
        rows.append(
            {
                "key": key,
                "pal": pal,
                "dim": dim,
                "n": n,
                "n1": n1,
                "nNL": n - n1,
                "nH": nH,
                "|o|": len(vs),
            }
        )
    return rows


def prove_p5_nl_is_1dnc() -> dict:
    p = 5
    Dall = _load_D(p)
    D1 = build_1d(p, plus_type=True)
    DH = build_1d_odot_nc(p, plus_type=True)
    DHm = build_1d_odot_nc(p, plus_type=False)
    rows = _aut3_counts(p, Dall, D1, DH)
    rows_m = _aut3_counts(p, Dall, D1, DHm)
    eq = all(r["nNL"] == r["nH"] for r in rows)
    minus_eq = all(a["nNL"] == b["nH"] for a, b in zip(rows, rows_m))
    # +++ dim-2 must differ under Type− swap
    plus_d2 = next(r for r in rows if r["dim"] == 2 and r["pal"] == (1, 1, 1))
    minus_d2 = next(r for r in rows_m if r["dim"] == 2 and r["pal"] == (1, 1, 1))
    swap_fails = plus_d2["nH"] != minus_d2["nH"]
    return {
        "proved": bool(eq and swap_fails and not minus_eq),
        "nH_sets": int(DH.shape[0]),
        "Nplus": int(Dall.shape[0]),
        "n1D": int(D1.shape[0]),
        "all_orbits_nH_eq_nNL": eq,
        "type_minus_also_eq": minus_eq,
        "type_minus_swap_changes_nH": swap_fails,
        "n_orbits": len(rows),
    }


def _all_circles(p, F):
    q = p * p
    mul, add, neg = F["mul"], F["add"], F["neg"]
    frob = [_pow(x, p, mul) for x in range(q)]
    out = []
    for t in range(q):
        for c in range(1, p):
            pts = []
            for u in range(q):
                d = add[u][neg[t]]
                if mul[d][frob[d]] == c:
                    pts.append(u)
            if len(pts) == p + 1:
                out.append(frozenset(pts))
    return out


def _label_one(payload):
    D, y, C, circles, reps, p = payload
    q = p * p
    Cp = np.diag(y) @ C @ np.diag(y)
    hoff = []
    for pts in circles:
        Sidx = [1 + u for u in pts]
        good = all(
            Cp[Sidx[i], Sidx[j]] > 0
            for i in range(len(Sidx))
            for j in range(i + 1, len(Sidx))
        )
        if not good:
            continue
        z = np.ones(q + 1)
        for i in Sidx:
            z[i] = -1.0
        if np.allclose(Cp @ z, p * z, atol=1e-5):
            hoff.append(pts)
    Dset = {u for u in range(q) if D[u]}
    out = []
    for key, S in reps:
        U = frozenset(x for x in S if x not in Dset)
        V = frozenset(x for x in S if x in Dset)
        n_h = sum(1 for pts in hoff if U <= pts and pts.isdisjoint(V))
        out.append((key, n_h))
    return out


def prove_overcount(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    expect = {5: 3.0, 7: 1.0}
    for p in primes:
        F = _kernel_field(p)
        C = paley_conference_prime_power(p).astype(np.float64)
        circles = _all_circles(p, F)
        parents = typeplus_y(p, plus_type=True)
        q = p * p
        add, mul, neg, chi = F["add"], F["mul"], F["neg"], F["chi"]
        frob = [_pow(x, p, mul) for x in range(q)]
        squares = [r for r in range(1, q) if chi[r] == 1]
        reps_map = {}
        for S in combinations(range(q), 3):
            key = _aut_key3(S, add, mul, neg, frob, squares)
            if _fp_dim(key, p) != 2 or key in reps_map:
                continue
            reps_map[key] = S
        reps = list(reps_map.items())
        parent_D = []
        for y in parents:
            D = tuple(1 if y[1 + u] < 0 else 0 for u in range(q))
            parent_D.append((D, y))
        payloads = [(D, y, C, circles, reps, p) for D, y in parent_D]
        labeled = defaultdict(int)
        W = min(default_workers(), max(1, len(payloads)))
        with ProcessPoolExecutor(max_workers=W) as ex:
            for part in ex.map(_label_one, payloads):
                for key, n_h in part:
                    labeled[key] += n_h
        Dall = _load_D(p)
        D1 = build_1d(p, plus_type=True)
        DH = build_1d_odot_nc(p, plus_type=True)
        rows = _aut3_counts(p, Dall, D1, DH)
        kappas = []
        for r in rows:
            if r["dim"] != 2 or r["nH"] == 0:
                continue
            kap = labeled[r["key"]] / r["nH"]
            kappas.append(kap)
        const = bool(kappas) and max(kappas) - min(kappas) < 1e-9
        kap0 = float(kappas[0]) if kappas else None
        want = expect[p]
        if not (const and kap0 is not None and abs(kap0 - want) < 1e-9):
            ok = False
        p7_eq = None
        if p == 7:
            d2 = [r for r in rows if r["dim"] == 2]
            p7_eq = all(r["nH"] == r["nNL"] for r in d2)
            if p7_eq:
                ok = False
        blocks.append(
            {
                "p": p,
                "kappa": kap0,
                "constant": const,
                "n_d2": sum(1 for r in rows if r["dim"] == 2),
                "p7_nH_eq_nNL": p7_eq,
                "nH_sets": int(DH.shape[0]),
            }
        )
    return {"proved": ok, "blocks": blocks}


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "n_NL_named_in_p": False,
        "n_H_named_in_p": False,
        "Nplus_named": False,
        "kappa_named_in_p": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.286  n_H = 1D⊙Hoffman-nc; n_NL=n_H at p=5", flush=True)
    A = prove_p5_nl_is_1dnc()
    print(
        f"  A p=5 n_NL=n_H={A['proved']} sets={A['nH_sets']} "
        f"swap_changes={A['type_minus_swap_changes_nH']}",
        flush=True,
    )
    B = prove_overcount([5, 7])
    print(f"  B overcount={B['proved']} blocks={B['blocks']}", flush=True)
    op = prove_open()
    out = {
        "prop": "15.286",
        "L_status": "OPEN",
        "proved": {
            "p5_nNL_eq_nH": A["proved"],
            "overcount_constant": B["proved"],
            "p7_nH_eq_nNL": False,
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "n_NL_named_in_p": False,
            "n_H_named_in_p": False,
            "Nplus_named": False,
        },
        "status": {"floor": False, "identity_A": A["proved"], "identity_B": B["proved"]},
        "A": A,
        "B": B,
        "open": op,
        "n1d_set_count_5": n1d_set_count(5),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15286.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
