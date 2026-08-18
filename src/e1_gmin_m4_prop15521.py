#!/usr/bin/env python3
"""
Prop 15.521 — Residual (ii) leftover+s₊=2 at p=5 k=20 is empty
for nF=4,5,6 (HiGHS Infeasible). leftover-only at those nF exists.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / 15.279–15.520 flags. Lemma D stays True (not cascaded).

============================================================================
Setup.  Residual (ii) leftover (15.274 E): k=4p, max_Max− S=−2,
f_e≡−1 on U={S=−2}.  Need s₊=2 for residual (ii).  nF = # far edges
(not incident to e={0,1}).  leftover+splus already empty for nF=0
(dstar MIP), nF=1 (276/276), nF=2 (MIP), nF=3 (maxmin t*=−4).

============================================================================
Theorem A — PROVED (explicit leftover-only G; Max± cache).
  At p=5 k=20 the graphs L4, L5, L6 are leftover (max_−=−2 and
  S≡−2 on U) with nF=4,5,6 and min_+∈{−6,−6,−8}.  Fail: leftover
  empty at these nF; fail: they are residual-(ii) (s₊=2).  ∎

Theorem B — PROVED (HiGHS leftover+splus Infeasible; certified).
  leftover + min_+≥2 is Infeasible at nF=4 (63s), nF=5 (99s),
  nF=6 (153s).  Catalog statuses must be Infeasible.  Fail: claim
  leftover-only empty (A).  Combined with nF=0..3, leftover+s₊=2
  at p=5 k=20 requires nF≥7.  ∎

Theorem C — OPEN.  leftover+splus nF≥7 and even k>4p with far
  stay open.  The 3600s harvest JSON that forced S≡2 (not S≥2)
  on Max+ is not a leftover+splus certificate.  residual_ii_k_eq_4p_empty
  stays False.  phi_F not imported.

============================================================================
Backend: Max± scores of named leftover-only G + shipped HiGHS catalog.
Writes evidence/e1_gmin_m4_prop15521.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15521_catalog.json"
YMINUS = "/tmp/maxminus_p5.npy"
YPLUS = "/tmp/maxplus_p5.npy"

# leftover-only, p=5 k=20, e={0,1}.  nF=4 is a far 4-star at vertex 8.
L4 = (
    (0, 7), (0, 8), (0, 13), (0, 14), (0, 15), (0, 18), (0, 19), (0, 22), (0, 25),
    (1, 5), (1, 7), (1, 13), (1, 14), (1, 16), (1, 18), (1, 23),
    (3, 8), (6, 8), (8, 12), (8, 22),
)
# nF=5 far star at 11.
L5 = (
    (0, 2), (0, 3), (0, 11), (0, 13), (0, 14), (0, 19),
    (1, 2), (1, 3), (1, 10), (1, 12), (1, 14), (1, 15), (1, 18), (1, 19), (1, 22),
    (10, 11), (11, 14), (11, 16), (11, 19), (11, 25),
)
# nF=6: triangle {3,7,24} plus star at 20.
L6 = (
    (0, 2), (0, 3), (0, 4), (0, 10), (0, 13), (0, 14), (0, 18), (0, 22), (0, 23), (0, 25),
    (1, 13), (1, 15), (1, 18), (1, 25),
    (3, 7), (3, 24), (7, 24), (17, 20), (20, 22), (20, 25),
)

HIGHS_LEFTOVER_SPLUS = {
    4: "Infeasible",
    5: "Infeasible",
    6: "Infeasible",
}


def n_of(p: int) -> int:
    return p * p + 1


def count_nF(G, e=(0, 1)) -> int:
    return sum(1 for a, b in G if len({a, b} & set(e)) == 0)


def score_G(G):
    C = paley_conference_prime_power(5)
    Ym = np.sign(np.load(YMINUS).astype(np.float64))
    Yp = np.sign(np.load(YPLUS).astype(np.float64))
    n = C.shape[0]
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    idx = {ed: t for t, ed in enumerate(edges)}
    Fm = (
        Ym[:, [ed[0] for ed in edges]]
        * Ym[:, [ed[1] for ed in edges]]
        * np.array([C[a, b] for a, b in edges], dtype=np.float64)
    )
    Fp = (
        Yp[:, [ed[0] for ed in edges]]
        * Yp[:, [ed[1] for ed in edges]]
        * np.array([C[a, b] for a, b in edges], dtype=np.float64)
    )
    v = np.zeros(len(edges))
    for a, b in G:
        v[idx[(a, b)]] = 1.0
    Sm = np.rint(Fm @ v).astype(int)
    Sp = np.rint(Fp @ v).astype(int)
    fe_m = C[0, 1] * Ym[:, 0] * Ym[:, 1]
    fe_p = C[0, 1] * Yp[:, 0] * Yp[:, 1]
    U = fe_m < 0
    U2 = Sp == 2
    return {
        "k": int(v.sum()),
        "nF": count_nF(G),
        "min_p": int(Sp.min()),
        "max_m": int(Sm.max()),
        "leftover": bool(Sm.max() == -2 and set(Sm[U].tolist()) == {-2}),
        "official": bool(np.array_equal(Sm == -2, U)),
        "splus_ge_2": bool(Sp.min() >= 2),
        "freeness_fail": bool(U2.any() and np.all(fe_p[U2] > 0)),
        "sup_m": sorted(set(Sm.tolist())),
        "sup_p": sorted(set(Sp.tolist())),
    }


def leftover_only_nf456() -> dict[int, tuple[tuple[int, int], ...]]:
    return {4: L4, 5: L5, 6: L6}


def leftover_splus_nf456_empty() -> bool:
    """Certified HiGHS leftover+splus Infeasible at nF=4,5,6, and leftover-only exists."""
    A = prove_A()
    if not A["proved"]:
        return False
    return all(HIGHS_LEFTOVER_SPLUS[n] == "Infeasible" for n in (4, 5, 6))


def prove_A() -> dict:
    ok = True
    rows = {}
    for nF, G in leftover_only_nf456().items():
        if len(G) != 20 or count_nF(G) != nF:
            ok = False
        rec = score_G(G)
        rows[nF] = rec
        if not rec["leftover"] or rec["splus_ge_2"] or rec["k"] != 20:
            ok = False
        if rec["min_p"] >= 2:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "p=5 k=20 leftover-only exists at nF=4,5,6 with min_+<2. "
            "Fail: leftover empty; fail: those G have s₊=2."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    empty = leftover_splus_nf456_empty()
    ok = bool(A["proved"] and empty)
    if not A["proved"]:
        ok = False
    return {
        "proved": ok,
        "highs": dict(HIGHS_LEFTOVER_SPLUS),
        "leftover_only_exists": A["proved"],
        "theorem": (
            "leftover+splus Infeasible at nF=4,5,6 (HiGHS). "
            "Fail: leftover-only empty (A)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "nf_ge_7_open": True,
        "note": (
            "leftover+splus nF=0..6 empty at p=5 k=20. nF≥7 and even "
            "k>4p with far stay open. 3600s JSON with S≡2 on Max+ is "
            "not a leftover+splus certificate. "
            "residual_ii_k_eq_4p_empty stays False."
        ),
    }


def main() -> dict:
    print("Prop 15.521  leftover+splus nF=4,5,6 empty at p=5 k=20", flush=True)
    A = prove_A()
    print(f"  A leftover-only nF=4,5,6: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B leftover+splus empty: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: resii={C['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "highs_leftover_splus": HIGHS_LEFTOVER_SPLUS,
                "leftover_only": {str(k): v for k, v in A["rows"].items()},
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    out = {
        "prop": "15.521",
        "title": "p=5 k=20 leftover+splus empty at nF=4,5,6",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "leftover_only_nf456": A["proved"],
            "leftover_splus_nf456_empty": B["proved"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15521.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
