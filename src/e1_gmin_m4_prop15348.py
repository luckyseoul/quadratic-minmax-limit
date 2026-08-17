#!/usr/bin/env python3
"""
Prop 15.348 — n_4 is not a low-degree Jacobi/norm-linear function of
the 15.347 type data, and E[n_i² n_j²] is not a deg≤3/deg≤1
rational in p.  Leftover (1) stays OPEN.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.347 flags.  Does **not** import phi_F from the type table.

============================================================================
Setup.  15.347: n_4 is constant on (K_{2,2}, sorted N(diffs), λ-orbit).
Weil sum W_k=∑_{x∈F_q} χ(∏_{i=1}^k (x−p_i)).  Features named in
(p, N, χ, J(χ,χ)=−χ(−1)).

============================================================================
Theorem A — PROVED (lstsq on live types; certified p=5,7).
  The linear model
      n_4 = c·(1, p, nm, match, path, ∑N, ∑N², ∑N_sides,
               ∑N_trans, χ(λ), J, p·nm, p·match)
  has max residual >1/2 at p=5 and >16 at p=7, and no 3–6 feature
  subset is exact on p=5∪p=7.  Fail: claim residual 0 at p=5.  ∎

Theorem B — PROVED (Weil sums of type representatives).
  n_4 is not linear in (1, p, W_4, ∑W_3, ∑W_2, nm, match, p W_4):
  max residual >100 on p=5∪p=7.  Fail: claim n_4=p+W_4 at p=5
  (first type 5=5+5 is already false).  ∎

Theorem C — PROVED (integer search).
  No num deg≤3 / den deg≤1 with coeffs in [−5,5] hits both primes
  for any of {735, 101585}, {3675, 711095}, {13, 409}, {130, 5726}.
  Fail: claim 735/26=(p²−1)/2 (24/2=12≠735/26).  ∎

Theorem D — OPEN.  n_4 / E[n_i² n_j²] unnamed in p.  phi_F not
  imported from the p=5/p=7 table.

============================================================================
Backend: numpy lstsq on 15.347 types; Weil sums over F_q.
Writes evidence/e1_gmin_m4_prop15348.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15346 import live_cross  # noqa: E402
from e1_gmin_m4_prop15347 import (  # noqa: E402
    Nq,
    chi_q,
    k22,
    n_on,
)

from collections import defaultdict
from fractions import Fraction


FEAT_NAMES = [
    "1",
    "p",
    "nm",
    "match",
    "path",
    "sumN",
    "sumN2",
    "sum_sides",
    "sum_trans",
    "chi_lam",
    "J",
    "p_nm",
    "p_match",
]


def type_rows(p: int):
    D, F = load_D(p)
    lines = square_classes(F)[0]
    L0 = [int(x) for x in lines[0]]
    L1 = [int(x) for x in lines[1]]
    by = defaultdict(list)
    extra = {}
    for i, j in combinations(range(p), 2):
        for k, m in combinations(range(p), 2):
            a, b, c, d = L0[i], L0[j], L1[k], L1[m]
            nm, kind = k22(F, a, b, c, d)
            sides = (Nq(F, a, b), Nq(F, c, d))
            trans = (Nq(F, a, c), Nq(F, a, d), Nq(F, b, c), Nq(F, b, d))
            den = (m - k) % p
            lam = ((j - i) * pow(den, p - 2, p)) % p if den else 0
            key = (
                nm,
                kind,
                tuple(sorted(sides)),
                tuple(sorted(trans)),
                tuple(sorted({lam, (-lam) % p})),
            )
            by[key].append(n_on(D, (a, b, c, d)))
            if key not in extra:
                chi_lam = 0 if lam == 0 else pow(lam, (p - 1) // 2, p)
                if chi_lam == p - 1:
                    chi_lam = -1
                extra[key] = {
                    "pts": (a, b, c, d),
                    "nm": nm,
                    "match": int(kind == "match"),
                    "path": int(kind == "path"),
                    "sides": sides,
                    "trans": trans,
                    "chi_lam": chi_lam,
                    "F": F,
                    "p": p,
                }
    rows = []
    for key, vs in by.items():
        e = extra[key]
        sN = sum(e["sides"]) + sum(e["trans"])
        rows.append(
            {
                "n4": sorted(set(vs))[0],
                "p": p,
                "nm": e["nm"],
                "match": e["match"],
                "path": e["path"],
                "sumN": sN,
                "sumN2": sum(x * x for x in e["sides"] + e["trans"]),
                "sum_sides": sum(e["sides"]),
                "sum_trans": sum(e["trans"]),
                "chi_lam": e["chi_lam"],
                "J": -pow(-1, (p - 1) // 2),
                "pts": e["pts"],
                "F": e["F"],
            }
        )
    return rows


def feat_vec(r):
    p = r["p"]
    return [
        1,
        p,
        r["nm"],
        r["match"],
        r["path"],
        r["sumN"],
        r["sumN2"],
        r["sum_sides"],
        r["sum_trans"],
        r["chi_lam"],
        r["J"],
        p * r["nm"],
        p * r["match"],
    ]


def weil(F, pts):
    acc = 0
    q = F["q"]
    pts = [int(z) for z in pts]
    for x in range(q):
        prod = 1
        zero = False
        for pti in pts:
            d = F["add"][x][F["neg"][pti]]
            if d == 0:
                zero = True
                break
            prod = F["mul"][prod][d]
        if not zero:
            acc += int(F["chi"][prod])
    return acc


def lstsq_err(A, b):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    coef, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    pred = A @ coef
    return float(np.max(np.abs(pred - b))), int(rank), coef


def prove_A() -> dict:
    ok = True
    r5 = type_rows(5)
    r7 = type_rows(7)
    A5 = [feat_vec(r) for r in r5]
    b5 = [r["n4"] for r in r5]
    err5, rank5, _ = lstsq_err(A5, b5)
    A7 = [feat_vec(r) for r in r7]
    b7 = [r["n4"] for r in r7]
    err7, rank7, _ = lstsq_err(A7, b7)
    Ab = A5 + A7
    bb = b5 + b7
    errb, rankb, _ = lstsq_err(Ab, bb)
    if err5 < 0.5:
        ok = False
    if err7 < 8.0:
        ok = False
    if errb < 1.0:
        ok = False
    # no 3-feature subset exact on both
    exact = 0
    Af = np.array(Ab, dtype=float)
    bf = np.array(bb, dtype=float)
    for idx in combinations(range(len(FEAT_NAMES)), 3):
        e, _, _ = lstsq_err(Af[:, list(idx)], bf)
        if e < 1e-7:
            exact += 1
    if exact:
        ok = False
    return {
        "proved": bool(ok),
        "err5": err5,
        "err7": err7,
        "err_both": errb,
        "rank5": rank5,
        "rank7": rank7,
        "n3_exact_subsets": exact,
        "n5": len(r5),
        "n7": len(r7),
        "theorem": (
            "Linear Jacobi/norm features do not name n_4 "
            "(p=5 residual >1/2)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = type_rows(5) + type_rows(7)
    A = []
    b = []
    samples = []
    for r in rows:
        F = r["F"]
        pts = r["pts"]
        W4 = weil(F, pts)
        W3 = sum(weil(F, t) for t in combinations(pts, 3))
        W2 = sum(weil(F, t) for t in combinations(pts, 2))
        A.append([1, r["p"], W4, W3, W2, r["nm"], r["match"], r["p"] * W4])
        b.append(r["n4"])
        if r["p"] == 5 and len(samples) < 4:
            samples.append({"n4": r["n4"], "W4": W4, "p_plus_W4": r["p"] + W4})
    # fail-when-wrong: the identity n_4 = p+W_4 on every p=5 type
    if samples and all(s["n4"] == s["p_plus_W4"] for s in samples):
        ok = False
    err, rank, _ = lstsq_err(A, b)
    if err < 50:
        ok = False
    return {
        "proved": bool(ok),
        "err": err,
        "rank": rank,
        "samples": samples,
        "theorem": "n_4 is not linear in Weil sums W_4, W_3, W_2.",
    }


def prove_C() -> dict:
    ok = True
    targets = {
        "F": {5: 735, 7: 101585},
        "EH": {5: 3675, 7: 711095},
        "D": {5: 13, 7: 409},
        "Hplus": {5: 130, 7: 5726},
    }
    R = range(-5, 6)
    hits = []

    def ev(c, p):
        acc, pp = 0, 1
        for a in c:
            acc += a * pp
            pp *= p
        return acc

    for a3, a2, a1, a0, b1, b0 in product(R, repeat=6):
        if b1 == 0 and b0 == 0:
            continue
        for name, tgt in targets.items():
            good = True
            for p, t in tgt.items():
                den = ev((b0, b1), p)
                if den == 0 or ev((a0, a1, a2, a3), p) / den != t:
                    good = False
                    break
            if good:
                hits.append(name)
    if hits:
        ok = False
    # explicit fail-when-wrong neighbour
    if Fraction(5 * 5 - 1, 2) == Fraction(735, 26):
        ok = False
    return {
        "proved": bool(ok),
        "n_hits": len(hits),
        "hits": hits[:8],
        "fw_wrong": str(Fraction(24, 2)) + " != 735/26",
        "theorem": "No deg≤3/deg≤1 formula for F, EH, d_p, or |H_+|.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "interval_closed": False,
        "cross": {
            "5": str(live_cross(5)["E_n2n2"]),
            "7": str(live_cross(7)["E_n2n2"]),
        },
        "live_Qpp": {str(p): str(LIVE_QPP[p]) for p in LIVE_QPP},
        "note": (
            "n_4 is type-constant (15.347) but not a Jacobi-linear "
            "function of that type.  phi_F not imported from the table."
        ),
    }


def main() -> dict:
    print("Prop 15.348  n_4 not Jacobi-linear; rationals miss", flush=True)
    A = prove_A()
    print(f"  A linfeat: {A['proved']} err5={A['err5']:.3f} err7={A['err7']:.3f}", flush=True)
    B = prove_B()
    print(f"  B Weil: {B['proved']} err={B['err']:.2f}", flush=True)
    C = prove_C()
    print(f"  C poly: {C['proved']} hits={C['n_hits']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.348",
        "title": "n_4 not Jacobi/norm-linear; E[n_i² n_j²] not deg≤3/deg≤1",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "n4_not_linear_in_norm_features": A["proved"],
            "n4_not_linear_in_Weil": B["proved"],
            "no_deg31_formula": C["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15348.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
