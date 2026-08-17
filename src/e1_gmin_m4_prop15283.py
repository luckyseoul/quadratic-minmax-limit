#!/usr/bin/env python3
"""
Prop 15.283 — n(S)=#{regular D∋S} is Aut_∞-invariant for |S|≤4.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.282 flags.  Does **not** name n(S) as a Paley-type formula
(false) and does **not** close the Aut_∞ mixture / floor.

============================================================================
Setup.  Aut_∞=⟨trans, ×□, Frob⟩ on F_q, q=p².  Regular sets D of
size p(p−1)/2 with y_∞=+1.  n(S)=#{such D: S⊂D}.
2-set n is already χ-class (15.279 V / 2-design).

============================================================================
Theorem A — PROVED at p=5 (full C(q,k)) and p=7 (full C(q,3) and
full C(q,4)).
  n(S) is constant on Aut_∞-orbits for |S|=2,3,4.
  Paley edge-type of a 3-set or 4-set does **not** determine n(S)
  (span ≥1 at p=5, ≥16 at p=7 on 3-sets; larger on 4-sets).
  Fail-when-wrong: drop Frob from the key (orbits split and n
  still constant on the finer key — not a disproof); claim Paley
  type determines n (false).  Unsigned inv is not Aut (15.278).
  Floor / N+(p) / mixture still OPEN: Aut_∞-constancy reduces
  Term_4 to a sum over Aut_∞ 4-orbits (34 at p=5, 114 at p=7)
  whose values are not yet a Max+-free formula in p.  ∎

Theorem B — PROVED at p=5,7.  If S lies on an F_p-line (AG(2,p)
  dim-1), then n(S) is a function of (Paley type, 3-AP, 4-AP).
  Split=0 on all dim-1 Aut_∞ orbits.  The same (Paley, AP) key
  does **not** determine n on dim-2 4-sets (span 10–55 at p=7).
  Fail-when-wrong: apply B to a spanning (dim-2) 4-set.  ∎

============================================================================
Backend: CPU ProcessPool over 4-sets.  GPU unused.
Writes evidence/e1_gmin_m4_prop15283.json
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

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _load_D(p: int) -> np.ndarray:
    Y = np.load(YPATH[p])
    return (np.rint(Y[Y[:, 0] > 0, 1:]) < 0).astype(np.int8)


def _pow_mul(x: int, e: int, mul: np.ndarray) -> int:
    if x == 0:
        return 0 if e > 0 else 1
    acc, b, ee = 1, int(x), int(e)
    while ee:
        if ee & 1:
            acc = int(mul[acc, b])
        b = int(mul[b, b])
        ee >>= 1
    return acc


def _paley_type(S, add, neg, chi) -> tuple:
    S = list(S)
    chs = []
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            d = int(add[S[i], neg[S[j]]])
            chs.append(0 if d == 0 else int(chi[d]))
    return tuple(sorted(chs))


def _aut_key(S, add, mul, neg, frob, squares) -> tuple:
    pts = [int(x) for x in S]
    best = None
    for src in pts:
        t = int(neg[src])
        raw = [int(add[x, t]) for x in pts if x != src]
        for s in squares:
            scaled = tuple(sorted(int(mul[u, s]) for u in raw))
            if best is None or scaled < best:
                best = scaled
            fr = tuple(sorted(int(frob[int(mul[u, s])]) for u in raw))
            if fr < best:
                best = fr
    return best


def _shard(payload):
    D, sets, add, mul, neg, frob, squares, chi = payload
    out = []
    for S in sets:
        nS = int(D[:, list(S)].min(axis=1).sum())
        out.append(
            (
                _aut_key(S, add, mul, neg, frob, squares),
                _paley_type(S, add, neg, chi),
                nS,
            )
        )
    return out


def _census(p: int, k: int, *, sample: int | None = None, seed: int = 0) -> dict:
    F = _kernel_field(p)
    q = p * p
    add = np.array(F["add"], dtype=np.int32)
    mul = np.array(F["mul"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    chi = np.array(F["chi"], dtype=np.int8)
    frob = np.array([_pow_mul(x, p, mul) for x in range(q)], dtype=np.int32)
    squares = [r for r in range(1, q) if chi[r] == 1]
    D = _load_D(p)
    if sample is None:
        sets = list(combinations(range(q), k))
    else:
        rng = np.random.default_rng(seed)
        seen = set()
        sets = []
        while len(sets) < sample:
            S = tuple(sorted(rng.choice(q, k, replace=False).tolist()))
            if S not in seen:
                seen.add(S)
                sets.append(S)
    W = min(default_workers(), max(2, len(sets) // 32))
    shards = np.array_split(np.arange(len(sets)), W)
    with ProcessPoolExecutor(max_workers=W) as ex:
        parts = list(
            ex.map(
                _shard,
                [
                    (
                        D,
                        [sets[int(i)] for i in sh],
                        add,
                        mul,
                        neg,
                        frob,
                        squares,
                        chi,
                    )
                    for sh in shards
                    if len(sh)
                ],
            )
        )
    aut: dict = defaultdict(list)
    pal: dict = defaultdict(list)
    for part in parts:
        for key, pt, nS in part:
            aut[key].append(nS)
            pal[pt].append(nS)
    aut_span = sum(1 for vs in aut.values() if max(vs) != min(vs))
    pal_span = {
        str(t): int(max(vs) - min(vs)) for t, vs in pal.items() if max(vs) != min(vs)
    }
    return {
        "p": p,
        "k": k,
        "Nplus": int(D.shape[0]),
        "n_sets": len(sets),
        "aut_orbits": len(aut),
        "aut_spanning": aut_span,
        "aut_constant": aut_span == 0,
        "paley_types": len(pal),
        "paley_spanning_types": len(pal_span),
        "paley_type_constant": len(pal_span) == 0,
        "paley_spans": pal_span,
        "workers": W,
    }


def prove_aut_invariance(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = []
    ok = True
    paley_dead = True
    for p in primes:
        if p == 5:
            r3 = _census(5, 3)
            r4 = _census(5, 4)
        else:
            r3 = _census(7, 3)
            r4 = _census(7, 4)
        rows.append({"3": r3, "4": r4})
        if not (r3["aut_constant"] and r4["aut_constant"]):
            ok = False
        if r3["paley_type_constant"] or r4["paley_type_constant"]:
            paley_dead = False
    return {
        "proved": ok,
        "paley_type_determines_n": (not paley_dead),
        "rows": rows,
    }


def _collinear_fp(pts: list[int], p: int) -> bool:
    if len(pts) < 3:
        return True
    a = pts[0]
    ax, ay = a % p, a // p
    found = False
    bx = by = 0
    for q in pts[1:]:
        dx, dy = (q % p - ax) % p, (q // p - ay) % p
        if dx or dy:
            bx, by = dx, dy
            found = True
            break
    if not found:
        return True
    for q in pts[1:]:
        dx, dy = (q % p - ax) % p, (q // p - ay) % p
        if (bx * dy - by * dx) % p != 0:
            return False
    return True


def _ap_flags(pts: list[int], p: int) -> tuple[bool, bool]:
    """3-AP / 4-AP after identifying the F_p-line with F_p."""
    a = pts[0]
    ax, ay = a % p, a // p
    direction = None
    for q in pts[1:]:
        dx, dy = (q % p - ax) % p, (q // p - ay) % p
        if dx or dy:
            direction = (dx, dy)
            break
    if direction is None:
        return False, False
    dx, dy = direction
    ts = []
    for q in pts:
        qx, qy = (q % p - ax) % p, (q // p - ay) % p
        if dx:
            t = (qx * pow(dx, p - 2, p)) % p
        else:
            t = (qy * pow(dy, p - 2, p)) % p
        ts.append(t)
    S = set(ts)
    ap3 = ap4 = False
    for t in S:
        for d in range(1, p):
            if (t + d) % p in S and (t + 2 * d) % p in S:
                ap3 = True
                if (t + 3 * d) % p in S:
                    ap4 = True
    return ap3, ap4


def prove_dim1_ap(primes: list[int] | None = None) -> dict:
    """B. dim-1 n is a function of (Paley, AP); dim-2 is not."""
    if primes is None:
        primes = [5, 7]
    rows = []
    ok = True
    dim2_not = True
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        add = np.array(F["add"], dtype=np.int32)
        neg = np.array(F["neg"], dtype=np.int32)
        chi = np.array(F["chi"], dtype=np.int8)
        D = _load_D(p)
        g1: dict = defaultdict(set)
        g2: dict = defaultdict(set)
        for S in combinations(range(q), 4):
            nS = int(D[:, list(S)].min(axis=1).sum())
            pal = _paley_type(S, add, neg, chi)
            if _collinear_fp(list(S), p):
                ap3, ap4 = _ap_flags(list(S), p)
                g1[(pal, ap3, ap4)].add(nS)
            else:
                g2[pal].add(nS)
        split1 = sum(1 for s in g1.values() if len(s) > 1)
        split2 = sum(1 for s in g2.values() if len(s) > 1)
        if split1:
            ok = False
        if split2 == 0:
            dim2_not = False
        rows.append(
            {
                "p": p,
                "dim1_classes": len(g1),
                "dim1_split": split1,
                "dim2_pal_split": split2,
            }
        )
    return {
        "proved": ok,
        "dim2_paley_does_not_determine": dim2_not,
        "rows": rows,
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "Nplus_named": False,
        "mixture_meets_floor": False,
        "n_S_named_formula": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.283  n(S) Aut_∞-invariant for |S|≤4", flush=True)
    A = prove_aut_invariance([5, 7])
    B = prove_dim1_ap([5, 7])
    print(
        f"  A aut-constant={A['proved']} "
        f"paley_determines_n={A['paley_type_determines_n']}",
        flush=True,
    )
    print(
        f"  B dim1-AP named={B['proved']} "
        f"dim2_paley_splits={B['dim2_paley_does_not_determine']}",
        flush=True,
    )
    for block in A["rows"]:
        for k, r in block.items():
            print(
                f"    p={r['p']} |S|={k} orbits={r['aut_orbits']} "
                f"span={r['aut_spanning']} paley_span_types={r['paley_spanning_types']}",
                flush=True,
            )
    op = prove_open()
    out = {
        "prop": "15.283",
        "L_status": "OPEN",
        "proved": {
            "aut_invariance_k_le_4": A["proved"],
            "paley_type_determines_n": A["paley_type_determines_n"],
            "dim1_ap_names_n": B["proved"],
            "dim2_paley_does_not_determine": B["dim2_paley_does_not_determine"],
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "Nplus_named": False,
            "mixture_meets_floor": False,
        },
        "status": {
            "floor": False,
            "identity_A": A["proved"],
            "identity_B": B["proved"],
        },
        "census": A["rows"],
        "dim1": B["rows"],
        "open": op,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15283.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
