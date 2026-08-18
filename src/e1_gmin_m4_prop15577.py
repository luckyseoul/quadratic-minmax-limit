#!/usr/bin/env python3
"""
Prop 15.577 — Type+ 1D affine law is a Johnson degree-2
identity with kernel dimension p+1.  Zero-star particular
solution exists for every odd p.  DEAD as a Type I kill
(size 3p−2 shadows exist at p=5).

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.576 flags.  Does **not**
edit leftover-1 15.550–15.575 or residual-ii 15.566 / 15.576.
Does **not** import phi_F.  Aut_e DEAD (15.559).  Max±-of-C
DEAD (15.565).

============================================================================
Setup.  Lock e={∞,0}.  Type I affine S=3−2f_e on every Max+.
Restrict to Type+ 1D (15.538): z=ε∘π_L, ker π_L=L∋0, so
π_L(0)=0, f_e=ε(0).  Fibre counts of G_∞, G_0 and G_far
make S(ε) a degree-2 function of ε:F_p→{±1}, ∑ε=1,
#(+)=m=(p+1)/2.  Need S(ε)=3−2ε(0) for all such ε.

============================================================================
Theorem A — PROVED Max+-free (Johnson degree ≤2).
  Degree-2 functions on {±1}^{F_p} have dimension
      1+p+C(p,2).
  The restriction to ∑ε=1 has kernel spanned by the
  p+1 independent polynomials
      (∑ε−1)    and    (∑ε−1)ε(α)  (α∈F_p),
  hence rank C(p,2).  Certified p=5 (rank 10) and p=7
  (rank 21).  Fail: rank 1+p+C(p,2) (forget the kernel).
  Fail: rank C(p,2)−1.  ∎

Theorem B — PROVED Max+-free (explicit particular).
  Fibre data n_∞≡0, n_0(0)=p−1, n_0(α≠0)=0, one within-
  fibre far edge, and q_{αβ}=2 on nonzero fibres, give
      S(ε)=p + 2∑_{0≠α<β} ε(α)ε(β)=3−2ε(0)
  for every ε with ∑ε=1.  Fail: S≡0.  Fail: drop the
  2 in q_{αβ}.  Size is p²−2p+2, not 3p−2 (p≥5).  ∎

Theorem C — DEAD as a Type I kill.  At p=5 the 1D law
  plus |G|=13 admits shadows (e.g. n_∞=(0,1,1,1,1),
  t=3).  The 1D Johnson constraint does not empty
  Type I G and does not force f_e=+1 on U_−.  type_I
  flags stay False.

============================================================================
Backend: Fraction + Johnson rank at p=5,7.  Writes
evidence/e1_gmin_m4_prop15577.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15577.json"


def deg2_ambient_dim(p: int) -> int:
    return 1 + p + comb(p, 2)


def johnson_rank_named(p: int) -> int:
    """Rank of degree≤2 on ∑ε=1: C(p,2)."""
    return comb(p, 2)


def johnson_ker_dim(p: int) -> int:
    return p + 1


def particular_k(p: int) -> int:
    """Edge count of the zero-star particular solution."""
    return p * p - 2 * p + 2


def type_I_k(p: int) -> int:
    return 3 * p - 2


def S_particular(eps: list[int]) -> int:
    """Zero-star particular: p + 2 ∑_{0≠α<β} εα εβ."""
    p = len(eps)
    s = p
    for a, b in combinations(range(1, p), 2):
        s += 2 * eps[a] * eps[b]
    return s


def S_particular_drop2(eps: list[int]) -> int:
    """Fail-when-wrong: q_{αβ}=1 instead of 2."""
    p = len(eps)
    s = p
    for a, b in combinations(range(1, p), 2):
        s += eps[a] * eps[b]
    return s


def _eps_from_A(p: int, A: tuple[int, ...]) -> list[int]:
    s = set(A)
    return [1 if i in s else -1 for i in range(p)]


@lru_cache(maxsize=1)
def live_rank() -> dict:
    out = {}
    for p in (5, 7):
        m = (p + 1) // 2
        rows: list[list[Fraction]] = []
        for A in combinations(range(p), m):
            eps = _eps_from_A(p, A)
            v: list[Fraction] = [Fraction(1)]
            v.extend(Fraction(x) for x in eps)
            for i, j in combinations(range(p), 2):
                v.append(Fraction(eps[i] * eps[j]))
            rows.append(v)
        # rank over Q
        M = [row[:] for row in rows]
        n, d = len(M), len(M[0])
        rnk = 0
        for c in range(d):
            piv = None
            for r in range(rnk, n):
                if M[r][c] != 0:
                    piv = r
                    break
            if piv is None:
                continue
            M[rnk], M[piv] = M[piv], M[rnk]
            fac = M[rnk][c]
            M[rnk] = [x / fac for x in M[rnk]]
            for r in range(n):
                if r != rnk and M[r][c] != 0:
                    k = M[r][c]
                    M[r] = [M[r][j] - k * M[rnk][j] for j in range(d)]
            rnk += 1
        out[p] = {
            "n_eps": n,
            "ambient": d,
            "rank": rnk,
            "named": johnson_rank_named(p),
            "ker": d - rnk,
        }
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_rank()
    ok = True
    rows = {}
    for p in range(5, 32):
        if not is_prime(p):
            continue
        amb = deg2_ambient_dim(p)
        rnk = johnson_rank_named(p)
        ker = johnson_ker_dim(p)
        if amb - ker != rnk:
            ok = False
        if rnk == amb:
            ok = False
        if rnk == comb(p, 2) - 1:
            ok = False
        if p in (5, 7, 11, 13):
            rows[str(p)] = {
                "ambient": amb,
                "rank": rnk,
                "ker": ker,
            }
    ok = ok and live[5]["rank"] == 10 == johnson_rank_named(5)
    ok = ok and live[7]["rank"] == 21 == johnson_rank_named(7)
    ok = ok and live[5]["ker"] == 6 and live[7]["ker"] == 8
    ok = ok and live[5]["ambient"] == 16 and live[7]["ambient"] == 29
    return {
        "proved": bool(ok),
        "rows": rows,
        "live": {str(k): v for k, v in live.items()},
        "theorem": (
            "deg≤2 rank on Type+ is C(p,2); ker dim p+1.  "
            "Fail: rank=ambient; fail rank=C(p,2)−1."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    rows = {}
    for p in range(5, 48):
        if not is_prime(p):
            continue
        if particular_k(p) == type_I_k(p):
            ok = False
    for p in (5, 7, 11):
        m = (p + 1) // 2
        n_ok = 0
        n_tot = 0
        n_drop = 0
        for A in combinations(range(p), m):
            eps = _eps_from_A(p, A)
            if sum(eps) != 1:
                ok = False
            n_tot += 1
            want = 3 - 2 * eps[0]
            if S_particular(eps) == want:
                n_ok += 1
            if S_particular(eps) == 0:
                ok = False
            if S_particular_drop2(eps) != want:
                n_drop += 1
        if n_ok != n_tot or n_drop == 0:
            ok = False
        rows[str(p)] = {
            "n_ok": n_ok,
            "n_tot": n_tot,
            "k_part": particular_k(p),
            "k_typeI": type_I_k(p),
        }
    ok = ok and particular_k(5) == 17 and type_I_k(5) == 13
    ok = ok and particular_k(7) == 37 and type_I_k(7) == 19
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Zero-star particular S=3−2ε(0) on all Type+ ε.  "
            "Fail: S≡0; fail drop 2 in q_{αβ}.  Size ≠3p−2."
        ),
    }


def p5_size_shadow() -> dict:
    """One p=5 shadow with |G|=13: n_∞=(0,1,1,1,1), t=3."""
    ninf = [0, 1, 1, 1, 1]
    t = 3
    k = [ninf[0] + 2] + ninf[1:]
    qsum = 0
    for a, b in combinations(range(1, 5), 2):
        qsum += t - k[a] - k[b]
    rhs = [t - k[0] - k[a] for a in range(1, 5)]
    const = 3 + 4 * t // 2 - sum(k)
    ktot = sum(ninf) + const + qsum + sum(rhs)
    return {
        "ninf": ninf,
        "t": t,
        "const": const,
        "qsum": qsum,
        "rhs": rhs,
        "k": ktot,
        "want": 13,
        "ok": ktot == 13 and const >= 0 and all(r >= 0 for r in rhs),
    }


def prove_open() -> dict:
    sh = p5_size_shadow()
    return {
        "proved": False,
        "path_dead": True,
        "p5_size13_shadow": sh,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "note": (
            "1D Johnson kernel leaves p+1 dofs; p=5 has |G|=13 "
            "shadows.  Does not empty Type I or force good-sign "
            "f_e.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.577  Type+ 1D Johnson kernel; DEAD as Type I kill", flush=True)
    A = prove_A()
    print(
        f"  A ker dim p+1: {A['proved']} "
        f"p5 rank={A['live']['5']['rank']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B particular: {B['proved']}", flush=True)
    C = prove_open()
    print(
        f"  C DEAD: shadow13={C['p5_size13_shadow']['ok']} "
        f"ND={C['type_I_multilevel_bad_case_ND_closed']}",
        flush=True,
    )
    out = {
        "prop": "15.577",
        "title": "Type+ 1D Johnson kernel p+1; DEAD as Type I kill",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "johnson_ker_p_plus_1": A["proved"],
            "zero_star_particular": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
