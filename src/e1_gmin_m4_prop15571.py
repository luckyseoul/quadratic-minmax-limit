#!/usr/bin/env python3
"""
Prop 15.571 — A6 pair product is a Gauss expansion; the
amp(3−χ(2))/2 guess dies at p=23.

Nondegenerate s1+s2≠0 (always when p≡3 (mod 4)):
    k̂(g)=1+2(x+y+xy),  x=ω^{−gs1}, y=ω^{−gs2}
    E_g[ẑ²ẑ(−2g)] = 8τ³χ(−2)/(p−1) Gauss[(1+2(x+y+xy))²
        (1+2(x^{−2}+y^{−2}+(xy)^{−2}))].
Degenerate s2=−s1 (only p≡1): k̂=1+2(x+x^{−1}), separate
expansion.  Fail: treat p=5 as a 3-point plus.
amp(3−χ(2))/2 matches p=7,11,19 and dies at p=23 (9/5).
amp(p−3)/4 dies at p=19.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.570 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.570.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.568 A6: +1 on t+{s1,s2,s1+s2}, −1 off t.  χ(−1)=−1
⇒ −s ∉ □ if s∈□, so no degenerates for p≡3.  15.567 ẑ=2G_χ k̂.

============================================================================
Theorem A — PROVED (monomial Gauss; p=5,7,11,13).
  Nondeg pairs match the (1+2(x+y+xy)) expansion.
  Degenerate pairs match the (1+2(x+x^{−1})) expansion.
  Fail: apply the 3-point expansion on p=5 (s1+s2=0).  ∎

Theorem B — PROVED (family means; p=7,11,19,23).
  A6=amp(3−χ(2))/2 at p=7,11,19; ratio=9/5 at p=23.
  A6=amp(p−3)/4 at p=7,11; dies at p=19.
  Fail: either guess as a p-law.  ∎

Theorem C — OPEN.  A6 family mean is not a short Jacobi
  integer in p.  A_full / Q_τ unnamed.  Do not interpolate
  (p−5)/15.  phi_F not imported.

============================================================================
Backend: field Gauss O(p²) pairs + family means.  Writes
evidence/e1_gmin_m4_prop15571.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15567 import chip, gauss_table, prod_dbl_k, twochi_amp_named
from e1_gmin_m4_prop15568 import fam_A6_named, mean_prod, squares

EV = ROOT / "evidence" / "e1_gmin_m4_prop15571.json"


def gauss_eval(p: int, tau: complex, m: int) -> complex:
    if m % p == 0:
        return 0j
    return chip(p, m) * tau


def a6_nondeg_named(p: int, s1: int, s2: int, tau: complex) -> complex:
    basis = ((0, 0, 1), (1, 0, 2), (0, 1, 2), (1, 1, 2))
    left: dict[tuple[int, int], int] = {}
    for i1, j1, c1 in basis:
        for i2, j2, c2 in basis:
            key = (i1 + i2, j1 + j2)
            left[key] = left.get(key, 0) + c1 * c2
    acc: dict[tuple[int, int], int] = {}
    for (i, j), c0 in left.items():
        for di, dj, cf in ((0, 0, 1), (-2, 0, 2), (0, -2, 2), (-2, -2, 2)):
            key = (i + di, j + dj)
            acc[key] = acc.get(key, 0) + c0 * cf
    S = 0j
    for (i, j), coef in acc.items():
        S += coef * gauss_eval(p, tau, (i * s1 + j * s2) % p)
    return 8 * tau**3 * chip(p, (-2) % p) / (p - 1) * S


def a6_deg_named(p: int, s1: int, tau: complex) -> complex:
    left: dict[int, int] = {}
    for e1, c1 in ((0, 1), (1, 2), (-1, 2)):
        for e2, c2 in ((0, 1), (1, 2), (-1, 2)):
            left[e1 + e2] = left.get(e1 + e2, 0) + c1 * c2
    acc: dict[int, int] = {}
    for e, c0 in left.items():
        for f, cf in ((0, 1), (-2, 2), (2, 2)):
            acc[e + f] = acc.get(e + f, 0) + c0 * cf
    S = 0j
    for e, coef in acc.items():
        S += coef * gauss_eval(p, tau, (e * s1) % p)
    return 8 * tau**3 * chip(p, (-2) % p) / (p - 1) * S


def a6_guess_3mchi2(p: int) -> float:
    return twochi_amp_named(p) * (3 - chip(p, 2)) / 2.0


def a6_guess_pm3_over4(p: int) -> float:
    return twochi_amp_named(p) * (p - 3) / 4.0


def pair_k(p: int, plus: set[int]) -> tuple[int, ...]:
    k = [0] * p
    for i in range(p):
        if i == 0:
            continue
        k[i] = 1 if i in plus else -1
    return tuple(k)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        G, omega = gauss_table(p)
        tau = G[1]
        S = sorted(squares(p))
        nd_err = 0.0
        dg_err = 0.0
        n_nd = n_dg = 0
        mis_as_3pt = 0.0
        for s1, s2 in combinations(S, 2):
            if (s1 + s2) % p == 0:
                plus = {s1, s2}
                got = prod_dbl_k(pair_k(p, plus), p, G, omega)
                pred = a6_deg_named(p, s1, tau)
                dg_err = max(dg_err, abs(got - pred))
                mis_as_3pt = max(mis_as_3pt, abs(got - a6_nondeg_named(p, s1, s2, tau)))
                n_dg += 1
            else:
                plus = {s1, s2, (s1 + s2) % p}
                got = prod_dbl_k(pair_k(p, plus), p, G, omega)
                pred = a6_nondeg_named(p, s1, s2, tau)
                nd_err = max(nd_err, abs(got - pred))
                n_nd += 1
        match = nd_err < 1e-8 and (n_dg == 0 or dg_err < 1e-8)
        fail3 = (p % 4 == 1 and mis_as_3pt > 1.0) or p % 4 == 3
        if not (match and fail3):
            ok = False
        rows[str(p)] = {
            "n_nondeg": n_nd,
            "n_deg": n_dg,
            "nondeg_err": float(nd_err),
            "deg_err": float(dg_err),
            "mis_3pt_err": float(mis_as_3pt),
            "match": bool(match),
            "fail_3pt_on_deg": bool(fail3),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "A6 nondeg/deg pair products are the two Gauss expansions. "
            "Fail 3-point expansion on p=5 degenerates."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (7, 11, 19, 23):
        G, omega = gauss_table(p)
        live = mean_prod(fam_A6_named(p), p, G, omega)
        g1 = a6_guess_3mchi2(p)
        g2 = a6_guess_pm3_over4(p)
        amp = twochi_amp_named(p)
        ratio = live / amp if abs(amp) > 1e-9 else float("nan")
        if p in (7, 11, 19):
            hit1 = abs(live - g1) < 1e-6
        else:
            hit1 = abs(live - g1) > 1.0
        if p in (7, 11):
            hit2_small = abs(live - g2) < 1e-6
            die2 = True
        else:
            hit2_small = True
            die2 = abs(live - g2) > 1.0
        if not (hit1 and hit2_small and die2):
            ok = False
        rows[str(p)] = {
            "live": float(live),
            "amp": amp,
            "ratio": float(ratio),
            "guess_3mchi2": g1,
            "guess_pm3_4": g2,
            "g1_ok_or_dies_p23": bool(hit1),
            "g2_ok_small_or_dies": bool(hit2_small and die2),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "amp(3−χ(2))/2 dies at p=23 (ratio 9/5). "
            "amp(p−3)/4 dies at p=19. Neither is a p-law."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A6_named_in_p": False,
        "A_full_dbl_named_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "A6 pair products are named expansions. Family mean is not "
            "a short p-integer. Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.571  A6 expansions; two amp-guesses die", flush=True)
    A = prove_A()
    print(f"  A expansions: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B guesses die: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.571",
        "title": "A6 pair Gauss expansions; amp(3−χ(2))/2 dies at p=23",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "a6_pair_expansion": A["proved"],
            "a6_amp_guesses_die": B["proved"],
            "A6_named_in_p": False,
            "Q_tau_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "A6 nondeg = Gauss[(1+2(x+y+xy))^2 (1+2(x^{-2}+y^{-2}+(xy)^{-2}))]; "
            "amp(3-chi(2))/2 dies at p=23"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
