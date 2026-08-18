#!/usr/bin/env python3
"""
Prop 15.500 — Strictly central Max+ and F_p^×-stabilizers of
NL T-orbits, certified at the live caches only.  D = {finite
x : z_x = −1}.  D is strictly central when D = −D.

    |NL ∩ {D=−D}| = 0 (p=5), 42 (p=7).

These two values match p(p−1) at p=7≡3 (mod 4) and 0 at
p=5≡1.  That matching is **not** a p-law: p=3≡3 has NL empty
so 0≠6.  Fail: include Type+ 1D (6≠0 at p=5, 54≠42 at p=7);
apply p(p−1) at p=5; extend the p≡3 formula to p=3.

Each of the 42 is the unique strictly central representative
of its T-orbit (odd char + T free).  Those T-orbits are the
−1-stab ones at p=7.  2p = 14 scale-orbits of size 3 is read
off 42, not an independent count.

    |1D ∩ {D=−D}| = (p²−1)/4

at p=3,5,7 (2, 6, 12).  Fail: n_1d (6≠30, 12≠140).

On NL T-orbits, F_p^×-stab ⊆ {±1} is certified at p=5,7 only:
  p=5, ord=4: no k-subset is scale-invariant (k≢0,1 mod 4).
  p=7, ord=6: same arithmetic obstruction.
  p=7, ord=3: d|(p−1)/2, arithmetic empty fails; 11440
  invariant k-sets, 8 lie in H_+, all Type+ 1D, 0 NL.
  Fail: claim those 8 are NL.
The congruence obstruction is **not** general (p=11 d=5,
p=13 d=3,6 all have k≡0 mod d).  No p-uniform reason that
⟨a⟩-invariant H+ is forced 1D.

n_T = (# central NL T-orbits) + nS_free(p−1) at these two
primes.  nS_free unnamed (1 and 12).  D unnamed.  phi_F not
imported.  Does **not** flip a leftover flag.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.499 flags.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15294 import _fp_dirs, _is_1d
from e1_gmin_m4_prop15298 import YPATH
from e1_gmin_m4_prop15499 import _pack

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15500.json"


@lru_cache(maxsize=4)
def all_designs(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    z = np.sign(Yp[:, 1 : 1 + q].astype(np.float64))
    D = (z < 0).astype(np.uint8)
    dirs = _fp_dirs(p)
    one = np.array([_is_1d(z[i], dirs, F["add"], p) for i in range(len(z))])
    return F, D, one


def central_mask(F, D: np.ndarray) -> np.ndarray:
    neg = np.array(F["neg"], dtype=np.int32)
    return np.all(D == D[:, neg], axis=1)


def hyp_nl_central(p: int) -> int:
    return p * (p - 1) if (p % 4 == 3) else 0


def hyp_1d_central(p: int) -> int:
    return (p * p - 1) // 4


def mul_ord_fp(a: int, p: int) -> int:
    acc, o = a % p, 1
    while acc != 1:
        acc = (acc * a) % p
        o += 1
        if o > p:
            return o
    return o


def scale_invariant_census(p: int, a: int) -> dict:
    """k-subsets of F_q invariant under x↦a x, vs H_+ / NL."""
    F, D, one = all_designs(p)
    q = F["q"]
    mul = F["mul"]
    k = p * (p - 1) // 2
    o = mul_ord_fp(a, p)
    masks_all = set(int(m) for m in _pack(D))
    masks_nl = set(int(m) for m in _pack(D[~one]))
    masks_1d = set(int(m) for m in _pack(D[one]))
    seen: set[int] = set()
    full: list[tuple[int, ...]] = []
    has0 = False
    for x in range(q):
        if x in seen:
            continue
        oset: list[int] = []
        y = x
        for _ in range(o + 2):
            if y in seen:
                break
            oset.append(y)
            seen.add(y)
            y = mul[a][y]
        if len(oset) == 1 and oset[0] == 0:
            has0 = True
        elif len(oset) == o:
            full.append(tuple(oset))
    hits = {"with0": 0, "without0": 0, "with0_nl": 0, "without0_nl": 0, "without0_1d": 0}
    n_cand = 0

    def tally(m: int, key: str) -> None:
        if m in masks_all:
            hits[key] += 1
        if m in masks_nl:
            hits[key + "_nl"] += 1
        if key == "without0" and m in masks_1d:
            hits["without0_1d"] += 1

    if k % o == 0:
        need = k // o
        if 0 <= need <= len(full):
            for comb in combinations(range(len(full)), need):
                m = 0
                for i in comb:
                    for x in full[i]:
                        m |= 1 << x
                n_cand += 1
                tally(m, "without0")
    if has0 and (k - 1) % o == 0:
        need = (k - 1) // o
        if 0 <= need <= len(full):
            for comb in combinations(range(len(full)), need):
                m = 1 << 0
                for i in comb:
                    for x in full[i]:
                        m |= 1 << x
                n_cand += 1
                tally(m, "with0")
    return {"a": a, "ord": o, "n_full": len(full), "n_cand": n_cand, "hits": hits}


def prove_A() -> dict:
    """NL strictly central count. Fail: include 1D; p=5 formula; p=3."""
    ok = True
    rows = {}
    for p in (5, 7):
        F, D, one = all_designs(p)
        c = central_mask(F, D)
        n_nl = int((c & ~one).sum())
        n_all = int(c.sum())
        hyp = hyp_nl_central(p)
        if n_nl != hyp:
            ok = False
        rows[p] = {
            "n_nl_central": n_nl,
            "n_all_central": n_all,
            "hyp": hyp,
            "n_1d_central": int((c & one).sum()),
        }
    if rows[5]["n_nl_central"] == 5 * 4:
        ok = False
    if rows[7]["n_nl_central"] == rows[7]["n_all_central"]:
        ok = False
    if rows[7]["n_nl_central"] == 0:
        ok = False
    # p=3: NL empty, p≡3 formula would say 6
    F3, D3, one3 = all_designs(3)
    c3 = central_mask(F3, D3)
    n3_nl = int((c3 & ~one3).sum())
    rows[3] = {"n_nl_central": n3_nl, "nH": int(len(D3)), "n1d": int(one3.sum())}
    if n3_nl == hyp_nl_central(3):
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """1D strictly central = (p²−1)/4. Fail: n_1d."""
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F, D, one = all_designs(p)
        c = central_mask(F, D)
        got = int((c & one).sum())
        hyp = hyp_1d_central(p)
        if got != hyp:
            ok = False
        if got == n_1d(p):
            ok = False
        rows[p] = {"n_1d_central": got, "hyp": hyp, "n_1d": n_1d(p)}
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """NL T-orbit F_p^×-stab ⊆ {±1}. Fail: p=7 order-3 hits are NL."""
    ok = True
    rows = {}
    # p=5: a=2,3 have order 4; k=10 ≢ 0,1 (mod 4)
    r5 = scale_invariant_census(5, 2)
    if r5["n_cand"] != 0 or r5["hits"]["without0_nl"] != 0:
        ok = False
    rows[5] = r5
    # p=7 order 6 (a=3): arithmetic empty
    r76 = scale_invariant_census(7, 3)
    if r76["n_cand"] != 0:
        ok = False
    rows["7_ord6"] = r76
    # p=7 order 3 (a=2): 8 Max+, all 1D, 0 NL
    r73 = scale_invariant_census(7, 2)
    if r73["hits"]["without0_nl"] != 0 or r73["hits"]["with0_nl"] != 0:
        ok = False
    if r73["hits"]["without0"] != 8:
        ok = False
    if r73["hits"]["without0_1d"] != 8:
        ok = False
    if r73["hits"]["without0"] == r73["hits"]["without0_nl"]:
        ok = False
    rows["7_ord3"] = r73
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "nS_free_named": False,
        "note": (
            "NL-central count and stab⊆{±1} certified p=5,7. "
            "n_T = 1_{p≡3} p(p−1) + nS_free(p−1) still has nS_free "
            "unnamed. D unnamed."
        ),
    }


def main() -> dict:
    print("15.500 central Max+ / F_p^× stab", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
