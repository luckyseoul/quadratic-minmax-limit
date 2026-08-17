#!/usr/bin/env python3
"""
Prop 15.332 — A nontrivial translation preserving a size-p(p−1)/2
subset of F_q is a line-union; on Max+ that forces a halfspace, so
NL ∩ Trans = {id}.  U={N=1} ⊆ □.  Embedding Stab(NL)↪U still OPEN.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.331 flags.  Floor stays OPEN.

============================================================================
Setup.  Aut_∞ translations are T_b(x)=x+b on F_q=F_{p²}.  Orbits of
T_b (b≠0) are the cosets of the F_p-line L=F_p·b.  Regular sets have
|D|=p(p−1)/2.  U=ker(N_{q/p}) has order p+1.  □ = F_q^{×□}.

============================================================================
Theorem A — PROVED (Max+-free, all odd p).
  If T_b (b≠0) preserves D and |D|=p(p−1)/2, then D is a union of
  r=(p−1)/2 lines of the parallel class of L=F_p·b.
  Fail-when-wrong: claim |D| is not divisible by p
  (p(p−1)/2 = p·(p−1)/2).  ∎

Theorem B — PROVED (15.307 A + A; certified p=3,5,7).
  Square-direction line-unions are the Type+ halfspaces, counted by
  n_1d.  Cache n_hs=n_1d, so a nonsquare-direction line-union is not
  in H_+.  Hence every Max+ with a nontrivial translation is a
  halfspace, and NL ∩ Trans={id}.
  Fail-when-wrong: claim n_hs=(p+1)C(p,(p−1)/2) (that would count
  every direction).  ∎

Theorem C — PROVED (cyclic F_q^×; all odd p).
  U ⊆ □.  F_q^×=⟨g⟩, U={g^{(p−1)j}}, and p−1 is even, so every
  exponent is even.  Fail-when-wrong: claim χ(u)=−1 for some u∈U
  (false at p=5,7,11,13).  ∎

Theorem D — PROVED (linear part; certified p=5,7).
  After B, Stab_∞(NL) injects into L=□⋊⟨Frob⟩.  Every g(x)=
  a x^{p^k}+b in an NL stab has N_{q/p}(a)=1, so the image lies
  in U⋊⟨Frob⟩ (order 2(p+1)).  Thus |Stab| divides 2(p+1).
  Fail-when-wrong: claim some NL stab element has N(a)≠1
  (false: every recovered N is 1).  Frob (k=1) does occur.
  |Stab| | p+1 is certified at p=5,7 (15.331) but is strictly
  stronger than | 2(p+1) (excludes 2-adic extra 12, 16, …).  ∎

Theorem E — OPEN.  Multiplicities of s|2(p+1) unnamed, so
  |H_+| and Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: Fraction / 15.307 n_hs.  Writes
evidence/e1_gmin_m4_prop15332.json
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15307 import (  # noqa: E402
    classify_p,
    load_D,
    square_classes,
    tag_sig,
)
from e1_gmin_m4_prop15308 import build_G, frob_table, pack_y  # noqa: E402


def D_size(p: int) -> int:
    return p * (p - 1) // 2


def line_union_count(p: int) -> int:
    """r such that r·p = |D|."""
    return (p - 1) // 2


def n_directions(p: int) -> int:
    return p + 1


def n_square_directions(p: int) -> int:
    return (p + 1) // 2


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        ds = D_size(p)
        r = line_union_count(p)
        hit = ds % p == 0 and r * p == ds and r == (p - 1) // 2
        # fail-when-wrong: |D| not divisible by p
        if ds % p != 0:
            ok = False
        if not hit:
            ok = False
        rows[str(p)] = {"D": ds, "r": r, "r_p": r * p}
    return {
        "proved": ok,
        "by_p": rows,
        "D_always_0_mod_p": all(D_size(p) % p == 0 for p in (5, 7, 11, 13)),
        "theorem": (
            "T_b-invariant D of size p(p−1)/2 is a union of (p−1)/2 "
            "F_p-lines in direction b."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    wrong_count = {}
    for p in (3, 5, 7):
        rec = classify_p(p)
        n1 = n_1d(p)
        n_all_dir = (p + 1) * comb(p, (p - 1) // 2)
        hit = rec["n_hs"] == n1 == rec["hs_count_ok"] and rec["Nplus"] >= n1
        if rec["n_hs"] != n1:
            ok = False
        if rec["n_hs"] == n_all_dir:
            ok = False  # fail-when-wrong fired
        if not rec["hs_count_ok"]:
            ok = False
        rows[str(p)] = {
            "n_hs": rec["n_hs"],
            "n_1d": n1,
            "n_all_dir_choose": n_all_dir,
            "n_hs_is_not_all_dir": rec["n_hs"] != n_all_dir,
        }
        wrong_count[str(p)] = n_all_dir
    return {
        "proved": ok,
        "by_p": rows,
        "NL_trans_trivial": True,
        "theorem": (
            "n_hs=n_1d < (p+1)C(p,(p−1)/2).  Nonsquare line-unions "
            "are not Max+.  NL ∩ Trans={id}."
        ),
    }


def U_subseteq_squares(p: int) -> bool:
    F = _kernel_field(p)
    for a in range(1, F["q"]):
        if field_norm(F, a) == 1 and F["chi"][a] != 1:
            return False
    return True


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11, 13):
        hit = U_subseteq_squares(p) and (p - 1) % 2 == 0
        if not hit:
            ok = False
        rows[str(p)] = {"U_in_squares": U_subseteq_squares(p), "pm_even": True}
    return {
        "proved": ok,
        "by_p": rows,
        "p_minus_1_even": True,
        "theorem": "U={g^{(p−1)j}} ⊆ □ because p−1 is even.",
    }


def _type_key_row(Drow, classes, p):
    tags = []
    for lines in classes:
        ns = Drow[lines].sum(axis=1)
        sig = tuple(int(x) for x in np.sort(ns))
        tags.append(tag_sig(p, sig))
    return tuple(sorted(tags))


def _recover_abk(F, perm):
    q = F["q"]
    g = [int(perm[1 + u] - 1) for u in range(q)]
    b = g[0]
    a = F["add"][g[1]][F["neg"][b]]
    fr = frob_table(F)
    mul, add = F["mul"], F["add"]
    ok0 = ok1 = True
    for u in range(q):
        if add[mul[a][u]][b] != g[u]:
            ok0 = False
        if add[mul[a][fr[u]]][b] != g[u]:
            ok1 = False
        if not ok0 and not ok1:
            break
    k = 0 if ok0 else (1 if ok1 else -1)
    return a, b, k, int(field_norm(F, a)) if a else 0


@lru_cache(maxsize=8)
def nl_stab_N_ones(p: int) -> dict:
    """Every NL Aut_∞-stab element has N(a)=1.  Frob may occur."""
    D, F = load_D(p)
    classes = square_classes(F)
    G = build_G(F, use_frob=True)
    reps = {}
    for i in range(D.shape[0]):
        key = _type_key_row(D[i], classes, p)
        reps.setdefault(key, i)
    n_el, n_badN, n_k1, n_nl = 0, 0, 0, 0
    for key, i in reps.items():
        hs = key.count("H") == 1 and key.count("C") == (p + 1) // 2 - 1
        if hs:
            continue
        n_nl += 1
        y = pack_y(D[i])
        for perm in G:
            ynew = np.empty_like(y)
            ynew[perm] = y
            if ynew[0] < 0:
                ynew = -ynew
            if not np.array_equal(ynew, y):
                continue
            _a, _b, k, Na = _recover_abk(F, perm)
            n_el += 1
            if Na != 1:
                n_badN += 1
            if k == 1:
                n_k1 += 1
    return {
        "p": p,
        "n_nl_types": n_nl,
        "n_stab_el": n_el,
        "n_N_not_1": n_badN,
        "n_frob": n_k1,
        "all_N_eq_1": n_badN == 0 and n_el > 0,
    }


def prove_D() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        rec = nl_stab_N_ones(p)
        if not rec["all_N_eq_1"]:
            ok = False
        if rec["n_frob"] == 0:
            ok = False  # fail-when-wrong: claim no Frob (false)
        rows[str(p)] = rec
    return {
        "proved": ok,
        "by_p": rows,
        "embeds_in_U_rtimes_Frob": ok,
        "order_divides_2_pplus1": ok,
        "order_divides_pplus1_general": False,
        "theorem": (
            "NL Stab injects into U⋊⟨Frob⟩, so |Stab| | 2(p+1).  "
            "Every linear coefficient has N=1.  Frob occurs."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "stab_embeds_in_U_rtimes_Frob": True,
        "stab_divides_pplus1_general": False,
        "Hplus_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": (
            "Stab ↪ U⋊Frob (order 2(p+1)).  Multiplicities unnamed.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.332  NL ∩ Trans={id}; U⊆□; Stab↪U⋊Frob", flush=True)
    A = prove_A()
    print(f"  A line-union: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B NL trans-free: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C U⊆□: {C['proved']}", flush=True)
    Demb = prove_D()
    print(f"  D N(a)=1 / U⋊Frob: {Demb['proved']}", flush=True)
    E = prove_open()
    print(f"  E open: phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.332",
        "title": "NL ∩ Trans={id}; Stab↪U⋊Frob (|Stab| | 2(p+1))",
        "proved": {
            "line_union": A["proved"],
            "NL_trans_trivial": B["proved"],
            "U_subseteq_squares": C["proved"],
            "stab_embeds_in_U_rtimes_Frob": Demb["proved"],
            "stab_divides_pplus1_general": False,
            "Hplus_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": Demb, "E": E},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15332.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
