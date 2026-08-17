#!/usr/bin/env python3
"""
Prop 15.401 — Mix-aware floor: λ is not constant on 15.307
D-types (except p=5 and the p=7 size-294 type at j=2).
Two size-1176 NL types have mean 8/3<6; the 294-type has
λ_2≡0.  1D at p=7 is not λ-constant (56 AP zeros).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.400 flags.

============================================================================
Setup.  15.307 D-types are square-dir signature multisets.
15.400: p=5 is {0}^100 ∪ {80/3}^30; p=7 NL does not vanish.
15.292: S_1D/q² ≥ 4p(p−1)/(p−2); n_1d = m C(p,m).

============================================================================
Theorem A — PROVED (cache fold; certified p=5,7).
  λ is D-type-constant at p=5 (RRR ≡0, CCH ≡80/3) and on the
  p=7 size-294 type at j=2 (≡0).  Every p=7 size-1176 type
  splits 392+392+392 (nunique=3).  Fail: claim nunique=1
  on a p=7 1176-type.  ∎

Theorem B — PROVED (A; p=7 j=2 type-means).
  Size-weighted means are 8/3, 40/3, 8, 8/3, 56/3, 0, 168/5
  and recover ensemble 3360/409.  Two 1176-types have mean
  8/3<6 (all 2352 rows below 6).  Size-294 has λ_2≡0.
  Fail: claim every NL type-mean is ≥6.  ∎

Theorem C — PROVED (1D split; p=7).
  The n_1d=140 type has mean 168/5 = 4p(p−1)/(p−2) but is
  not constant: 56 rows (AP) have λ=0 and three 28-blocks
  sit above.  Fail: claim every 1D lift has λ=S_1D.  ∎

Theorem D — OPEN.  Type-means other than 1D-average and the
  vanishing types are not a Gauss/Jacobi name in p.  A
  1D-only lower bound is <1 at p=7.  phi_F not imported.

============================================================================
Backend: 15.307 signatures + 15.304 G_U fold on p=5,7 caches.
Writes evidence/e1_gmin_m4_prop15401.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15302 import nontrivial_s  # noqa: E402
from e1_gmin_m4_prop15304 import gu_vec  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15400 import one_d_weight_lb, pointwise_floor, s1d_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15401_catalog.json"

P7_J2_MEANS = (
    Fraction(8, 3),
    Fraction(40, 3),
    Fraction(8),
    Fraction(8, 3),
    Fraction(56, 3),
    Fraction(0),
    Fraction(168, 5),
)


def nl_type_means_ge_6() -> bool:
    """True only if every NL D-type has mean λ≥6. It does not."""
    return False


@lru_cache(maxsize=8)
def _table(p: int, j: int) -> list[dict]:
    D, F = load_D(p)
    classes = square_classes(F)
    packs = []
    for lines in classes:
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1)
        packs.append([tuple(int(x) for x in row) for row in np.sort(ns, axis=1)])
    tags = [
        tuple(sorted(tag_sig(p, packs[di][i]) for di in range(len(classes))))
        for i in range(D.shape[0])
    ]
    N = load_N(p, F).astype(np.float64)
    s = nontrivial_s(F)
    perm = [F["mul"][F["inv"][s]][d] for d in range(F["q"])]
    dN = N - N[:, perm]
    xi = omega_set(F)[0]
    q2 = float(F["q"] ** 2)
    lam = (8.0 / (p + 1)) * np.abs(dN @ gu_vec(F, j, xi)) ** 2 / q2
    by: dict[tuple, list[float]] = defaultdict(list)
    for t, v in zip(tags, lam):
        by[t].append(float(v))
    rows = []
    for t, vs in sorted(by.items(), key=lambda kv: -len(kv[1])):
        a = np.array(vs)
        snapped = {round(x, 6) for x in a}
        rows.append(
            {
                "tag": str(t),
                "n": int(len(vs)),
                "mean": float(a.mean()),
                "min": float(a.min()),
                "max": float(a.max()),
                "nunique": len(snapped),
                "n0": int(np.sum(np.abs(a) < 1e-8)),
                "n_lt6": int(np.sum(a < 6 - 1e-9)),
            }
        )
    return rows


def _near(x: float, fr: Fraction) -> bool:
    return abs(x - float(fr)) < 1e-6


def prove_A() -> dict:
    ok = True
    t5 = _table(5, 2)
    t7 = _table(7, 2)
    # p=5 two types, each constant
    ok = ok and len(t5) == 2
    for r in t5:
        ok = ok and r["nunique"] == 1
        if r["n"] == 100:
            ok = ok and _near(r["mean"], Fraction(0))
        if r["n"] == 30:
            ok = ok and _near(r["mean"], Fraction(80, 3))
    # p=7: 1176-types split; 294 is constant 0
    n1176_split = 0
    n294_const0 = 0
    for r in t7:
        if r["n"] == 1176:
            ok = ok and r["nunique"] == 3
            if r["nunique"] == 1:
                ok = False
            n1176_split += 1
        if r["n"] == 294:
            ok = ok and r["nunique"] == 1 and _near(r["mean"], Fraction(0))
            n294_const0 += 1
    ok = ok and n1176_split == 4 and n294_const0 == 1
    return {
        "proved": bool(ok),
        "p5": t5,
        "p7_n1176_split": n1176_split,
        "p7_n294_const0": n294_const0,
        "theorem": (
            "λ constant on p=5 types and p=7 size-294. "
            "p=7 1176-types split 392³. Fail: nunique=1 there."
        ),
    }


def prove_B() -> dict:
    ok = True
    t7 = _table(7, 2)
    means = sorted((r["n"], r["mean"], r["n_lt6"]) for r in t7)
    # two 1176-types at 8/3, all rows <6
    n_83 = sum(1 for n, m, b in means if n == 1176 and _near(m, Fraction(8, 3)))
    n_403 = sum(1 for n, m, _b in means if n == 1176 and _near(m, Fraction(40, 3)))
    n_8 = sum(1 for n, m, _b in means if n == 1176 and _near(m, Fraction(8)))
    ok = ok and n_83 == 2 and n_403 == 1 and n_8 == 1
    rows_below = sum(r["n"] for r in t7 if r["n"] == 1176 and _near(r["mean"], Fraction(8, 3)))
    ok = ok and rows_below == 2352
    for r in t7:
        if r["n"] == 1176 and _near(r["mean"], Fraction(8, 3)):
            ok = ok and r["n_lt6"] == 1176
        if r["n"] == 294:
            ok = ok and r["n0"] == 294
        if r["n"] == 588:
            ok = ok and _near(r["mean"], Fraction(56, 3))
        if r["n"] == 140:
            ok = ok and _near(r["mean"], Fraction(168, 5))
    # weighted mix recovers 3360/409
    acc = sum(r["n"] * r["mean"] for r in t7)
    mix = acc / HPLUS[7]
    ok = ok and abs(mix - float(Fraction(3360, 409))) < 1e-6
    if n_83 == 0:
        ok = False
    # fail-when-wrong: every NL type-mean ≥6
    nl_ge6 = all(r["mean"] >= 6 - 1e-9 for r in t7 if r["n"] != n_1d(7))
    ok = ok and not nl_ge6
    if nl_ge6:
        ok = False
    return {
        "proved": bool(ok),
        "n_mean_8_3": n_83,
        "rows_at_8_3": rows_below,
        "mix": str(Fraction(3360, 409)),
        "nl_type_means_ge_6": nl_ge6,
        "means": [
            {"n": r["n"], "mean": r["mean"], "n_lt6": r["n_lt6"]} for r in t7
        ],
        "theorem": (
            "p=7 j=2 type-means 8/3,40/3,8,8/3,56/3,0,168/5. "
            "Fail: every NL type-mean ≥6."
        ),
    }


def prove_C() -> dict:
    ok = True
    t7 = _table(7, 2)
    hs = [r for r in t7 if r["n"] == n_1d(7)]
    ok = ok and len(hs) == 1
    r = hs[0]
    ok = ok and _near(r["mean"], s1d_named(7))
    ok = ok and r["n0"] == 56
    ok = ok and r["nunique"] == 4
    ok = ok and r["mean"] != 0
    if r["n0"] == 0 or r["nunique"] == 1:
        ok = False
    return {
        "proved": bool(ok),
        "n_1d": n_1d(7),
        "mean": r["mean"],
        "s1d": str(s1d_named(7)),
        "n0_AP": r["n0"],
        "nunique": r["nunique"],
        "theorem": (
            "p=7 1D mean is S_1D=168/5 but 56 AP rows have λ=0. "
            "Fail: every 1D lift has λ=S_1D."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "nl_type_means_ge_6": nl_type_means_ge_6(),
        "pointwise_floor": pointwise_floor(),
        "one_d_lb_p7": str(one_d_weight_lb(7)),
        "one_d_lb_p7_lt_6": one_d_weight_lb(7) < 6,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Mix-aware table: two 1176-types sit below 6; 294 vanishes. "
            "Type-means unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.401  mix-aware λ by 15.307 D-type; NL means can be <6", flush=True)
    A = prove_A()
    print(f"  A splits: {A['proved']} n1176_split={A['p7_n1176_split']}", flush=True)
    B = prove_B()
    print(f"  B means: {B['proved']} n_8/3={B['n_mean_8_3']} mix={B['mix']}", flush=True)
    C = prove_C()
    print(f"  C 1D split: {C['proved']} n0_AP={C['n0_AP']}", flush=True)
    D = prove_open()
    print(f"  D open: nl_ge6={D['nl_type_means_ge_6']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": {"p7_n1176_split": A["p7_n1176_split"]}, "B": B["means"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.401",
        "title": "Mix-aware: NL D-type means can be <6; 1D not λ-constant at p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "D_type_not_lambda_constant": A["proved"],
            "NL_type_means_can_be_lt_6": B["proved"],
            "one_d_not_constant_p7": C["proved"],
            "nl_type_means_ge_6": D["nl_type_means_ge_6"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15401.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
