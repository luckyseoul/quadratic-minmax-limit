#!/usr/bin/env python3
"""
Prop 15.289 — Hoffman T on a Type+ parent vanishes unless U is a
Paley clique; n(ω) is named for every ω⊂S (dim-2 3-set).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.288 flags.  Does **not** name τ_ω, n_H, n_NL, or N+(p).
Floor stays OPEN.

============================================================================
Setup.  U = {x∈S : L(x)∈σ} = S \\ D_A on a Type+ parent y=σ∘L
(15.288).  Hoffman T of C'=D_y C D_y (15.138–15.141, 15.286).
S ⊂ D_A Δ T ⇔ U ⊂ T and T ∩ (S\\U)=∅.  ω ranges over subsets of S.

============================================================================
Theorem A — PROVED (Max+-free: C' clique + Type+ sign).
  If U ⊂ T for a Hoffman circle T of Type+ C_y, then y|_U ≡ +1, so
      C'_ab = y_a C_ab y_b = χ_q(b−a)    (a,b∈U),
  and T a C'-clique forces χ_q(b−a)=+1.  Hence U is a Paley clique.
  Corollary: τ_ω = 0 unless ω is a Paley clique.  In particular
      τ_3=0 unless e=3,   τ_2(ab)=0 unless ab is Paley,
      and e=0 ⇒ τ_2=τ_3=0.
  Fail-when-wrong: a non-clique U with n_H>0; claim τ_2=0 on every
  Paley pair (e≥1).  ∎

Theorem B — PROVED (Max+-free counting; certified p=5,7).
  m=(p+1)/2, e=# Paley edges.  Non-collapsing Type+ L: m−e of them.
  Collapsing L: one per Paley edge (two points share a fiber).
    n(ω) = (m−e) C(p−3, m−|ω|) + extra(ω)
  extra(∅)=e C(p−2,m);
  extra({x})=C(p−2,m−1) if the edge opposite x is Paley, else 0;
  extra({a,b})=C(p−2,m−1) if ab is Paley, else 0;
  extra(S)=e C(p−2,m−2).
  ∑_{|ω|=k} n(ω) = n_k of 15.288.  Fail-when-wrong: drop extra;
  invert the 3-fiber binomial.  ∎

Theorem C — OPEN.  On Paley-clique ω, τ_ω is not a function of
  (p,e) (15.288 B) and is not named in p.  Λ_H=∑_ω n(ω)τ_ω still
  needs those T-weights.  n_H / n_NL unnamed.  Floor open.  ∎

============================================================================
Backend: CPU ProcessPool over Type+ parents.  GPU unused.
Writes evidence/e1_gmin_m4_prop15289.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from math import comb
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
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15284 import _aut_key3, _fp_dim, _paley3  # noqa: E402
from e1_gmin_m4_prop15285 import n1d_set_count  # noqa: E402
from e1_gmin_m4_prop15286 import _all_circles, _pow, typeplus_y  # noqa: E402
from e1_gmin_m4_prop15288 import nU_formula  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from workers import default_workers  # noqa: E402


def C(n: int, k: int) -> int:
    if k < 0 or n < 0 or k > n:
        return 0
    return comb(n, k)


def is_paley_clique(omega, pal_pairs) -> bool:
    pts = list(omega)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if frozenset((pts[i], pts[j])) not in pal_pairs:
                return False
    return True


def extra_collapsing(p: int, e: int, k: int, omega_is_special: bool) -> int:
    """Collapsing extra in Theorem B.  omega_is_special: opp-edge Paley
    (k=1) or the pair itself Paley (k=2)."""
    m = (p + 1) // 2
    if k == 0:
        return e * C(p - 2, m)
    if k == 3:
        return e * C(p - 2, m - 2)
    if k in (1, 2):
        return C(p - 2, m - 1) if omega_is_special else 0
    return 0


def n_omega_formula(
    p: int, e: int, k: int, *, omega_is_special: bool
) -> int:
    m = (p + 1) // 2
    return (m - e) * C(p - 3, m - k) + extra_collapsing(
        p, e, k, omega_is_special
    )


def n_omega_formula_inverted(
    p: int, e: int, k: int, *, omega_is_special: bool
) -> int:
    m = (p + 1) // 2
    return (m - e) * C(p - 3, m - k - 1) + extra_collapsing(
        p, e, k, omega_is_special
    )


def _pal_pairs(S, add, neg, chi):
    S = [int(x) for x in S]
    out = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            d = add[S[i]][neg[S[j]]]
            if d and chi[d] == 1:
                out.add(frozenset((S[i], S[j])))
    return out


def _opp_edge_paley(S, x, pal_pairs) -> bool:
    others = [int(y) for y in S if int(y) != int(x)]
    if len(others) != 2:
        return False
    return frozenset(others) in pal_pairs


def _dim2_reps(p: int):
    F = _kernel_field(p)
    q = p * p
    add, mul, neg, chi = F["add"], F["mul"], F["neg"], F["chi"]
    frob = [_pow(x, p, mul) for x in range(q)]
    squares = [r for r in range(1, q) if chi[r] == 1]
    reps = {}
    for S in combinations(range(q), 3):
        key = _aut_key3(S, add, mul, neg, frob, squares)
        if _fp_dim(key, p) != 2 or key in reps:
            continue
        pal = _paley3(S, add, neg, chi)
        pairs = _pal_pairs(S, add, neg, chi)
        e = len(pairs)
        reps[key] = {
            "S": S,
            "pal": pal,
            "pairs": pairs,
            "e": e,
        }
    return reps


def live_n_omega(p: int, S, pal_pairs):
    from e1_gmin_m4_prop15279 import (  # local to keep worker import light
        _linear_forms,
        is_type_plus_form,
    )

    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    counts: dict[frozenset, int] = defaultdict(int)
    S = [int(x) for x in S]
    for alpha, beta in forms:
        fibers = [((alpha * (x % p) + beta * (x // p)) % p) for x in S]
        for mask in range(1 << p):
            if bin(mask).count("1") != m:
                continue
            plus = {i for i in range(p) if mask >> i & 1}
            omega = frozenset(S[i] for i, t in enumerate(fibers) if t in plus)
            counts[omega] += 1
    return counts


def prove_nomega(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    invert_hit = invert_tot = 0
    drop_extra_hit = drop_extra_tot = 0
    for p in primes:
        reps = _dim2_reps(p)
        bad = 0
        checked = 0
        for rec in reps.values():
            S, pairs, e = rec["S"], rec["pairs"], rec["e"]
            live = live_n_omega(p, S, pairs)
            pred_sum = [0, 0, 0, 0]
            live_sum = [0, 0, 0, 0]
            for r in range(4):
                for omega in combinations(S, r):
                    om = frozenset(omega)
                    k = r
                    if k == 1:
                        special = _opp_edge_paley(S, next(iter(om)), pairs)
                    elif k == 2:
                        special = om in pairs
                    else:
                        special = False
                    pred = n_omega_formula(p, e, k, omega_is_special=special)
                    got = int(live.get(om, 0))
                    checked += 1
                    if pred != got:
                        bad += 1
                        ok = False
                    pred_sum[k] += pred
                    live_sum[k] += got
                    inv = n_omega_formula_inverted(
                        p, e, k, omega_is_special=special
                    )
                    invert_tot += 1
                    if inv == got:
                        invert_hit += 1
                    drop = (p + 1) // 2 - e
                    drop *= C(p - 3, ((p + 1) // 2) - k)
                    drop_extra_tot += 1
                    if drop == got and extra_collapsing(p, e, k, special):
                        drop_extra_hit += 1
            if tuple(pred_sum) != nU_formula(p, e):
                ok = False
            if tuple(live_sum) != nU_formula(p, e):
                ok = False
            if sum(pred_sum) != n1d_set_count(p):
                ok = False
        if not reps:
            ok = False
        blocks.append(
            {
                "p": p,
                "n_orbits": len(reps),
                "checked": checked,
                "bad": bad,
            }
        )
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "blocks": blocks,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "drop_extra_hit": drop_extra_hit,
        "drop_extra_tot": drop_extra_tot,
    }


def _vanish_one(payload):
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
    for key, S, pairs in reps:
        U = frozenset(x for x in S if x not in Dset)
        V = frozenset(x for x in S if x in Dset)
        n_h = sum(1 for pts in hoff if U <= pts and pts.isdisjoint(V))
        clique = is_paley_clique(U, pairs)
        out.append((key, tuple(sorted(U)), n_h, clique, len(U)))
    return out


def prove_vanish(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    blocks = []
    ok = True
    for p in primes:
        F = _kernel_field(p)
        C = paley_conference_prime_power(p).astype(np.float64)
        circles = _all_circles(p, F)
        parents = typeplus_y(p, plus_type=True)
        q = p * p
        recs = _dim2_reps(p)
        reps = [(key, rec["S"], rec["pairs"]) for key, rec in recs.items()]
        parent_D = []
        for y in parents:
            D = tuple(1 if y[1 + u] < 0 else 0 for u in range(q))
            parent_D.append((D, y))
        payloads = [(D, y, C, circles, reps, p) for D, y in parent_D]
        nonclique_pos = 0
        nonclique_tot = 0
        paley_pair_pos = 0
        paley_pair_tot = 0
        e3_u3_pos = 0
        e3_tot = 0
        e0_u23 = 0
        W = min(default_workers(), max(1, len(payloads)))
        with ProcessPoolExecutor(max_workers=W) as ex:
            for part in ex.map(_vanish_one, payloads):
                for key, _U, n_h, clique, uk in part:
                    e = recs[key]["e"]
                    if not clique:
                        nonclique_tot += 1
                        if n_h:
                            nonclique_pos += 1
                    if uk == 2 and clique:
                        paley_pair_tot += 1
                        if n_h:
                            paley_pair_pos += 1
                    if e == 3 and uk == 3:
                        e3_tot += 1
                        if n_h:
                            e3_u3_pos += 1
                    if e == 0 and uk >= 2 and n_h:
                        e0_u23 += 1
        if nonclique_pos:
            ok = False
        if paley_pair_tot == 0 or paley_pair_pos == 0:
            ok = False
        if e3_tot == 0 or e3_u3_pos == 0:
            ok = False
        if e0_u23:
            ok = False
        blocks.append(
            {
                "p": p,
                "n_parents": len(parents),
                "n_orbits": len(recs),
                "nonclique_pos": nonclique_pos,
                "nonclique_tot": nonclique_tot,
                "paley_pair_pos": paley_pair_pos,
                "paley_pair_tot": paley_pair_tot,
                "e3_u3_pos": e3_u3_pos,
                "e3_tot": e3_tot,
                "e0_u23_pos": e0_u23,
                "workers": W,
            }
        )
    return {"proved": ok, "blocks": blocks}


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "tau_omega_named": False,
        "n_H_named_in_p": False,
        "n_NL_named_in_p": False,
        "Nplus_named": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_ge_budget": bool(rhat_ge_budget_proved_general()),
    }


def main() -> dict:
    print("Prop 15.289  Paley-clique vanishing + n(ω)", flush=True)
    B = prove_nomega([5, 7])
    print(
        f"  B n(ω)={B['proved']} invert={B['invert_hit']}/{B['invert_tot']} "
        f"drop_extra={B['drop_extra_hit']}/{B['drop_extra_tot']}",
        flush=True,
    )
    A = prove_vanish([5, 7])
    print(f"  A vanish={A['proved']} blocks={A['blocks']}", flush=True)
    op = prove_open()
    out = {
        "prop": "15.289",
        "L_status": "OPEN",
        "proved": {
            "paley_clique_vanishing": A["proved"],
            "n_omega_formula": B["proved"],
            "invert_is_identity": B["invert_hit"] == B["invert_tot"],
            "phi_F_ge_6_proved_general": op["phi_F_ge_6"],
            "inequality_proved": False,
            "tau_omega_named": False,
            "n_H_named_in_p": False,
            "n_NL_named_in_p": False,
        },
        "status": {
            "floor": False,
            "identity_A": A["proved"],
            "identity_B": B["proved"],
        },
        "A": A,
        "B": B,
        "open": op,
        "backend": "ProcessPool CPU; GPU unused",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15289.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
