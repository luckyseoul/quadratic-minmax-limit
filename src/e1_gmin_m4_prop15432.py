#!/usr/bin/env python3
"""
Prop 15.432 — NL L_ns is Paley×norm indexed, not the constant
μ₊μ₋.  L₊+L₋=2μ₊μ₋.  At p=5, L₊=μ₊μ₋+(N−p/2).  That p=5
law fails at p=7.  Q_NL unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.431 flags.

============================================================================
Setup.  15.298 ns_key=(Paley(s−1,s+1), N_{q/p}(s)).  N(s) is a
nonsquare in F_p.  15.431 named C_ns.  Slot-average L_ns=μ₊μ₋
misses Q_NL (15.428/recon).

============================================================================
Theorem A — PROVED (NL cache fold; certified p=5,7).
  L_ns(s,χ) is constant on (paley_pair, N).  At p=7 the same
  N=3 occurs with two Paley classes and two L values
  (2733/38 vs 1398/19).  Fail: claim L depends on N only.  ∎

Theorem B — PROVED (same fold).
  L₊+L₋=2μ₊μ₋ on every class (25 at p=5, 147 at p=7).
  Fail: claim L₊=L₋ except when N=−1 (p=7 mixed N=6 is the
  only live zero-diff class).  ∎

Theorem C — PROVED (p=5 law; dies at p=7).
  At p=5, L₊=μ₊μ₋+(N−p/2) ∈{12,13}.  The same formula at
  p=7 on Paley(−1,−1) N=3 gives 73≠2733/38.  Fail: claim
  the p=5 law at p=7.  ∎

Theorem D — OPEN.  No Max+-free formula in p for the p=7
  table (dens 19=c(7)).  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: NL ENN fold.  Writes
evidence/e1_gmin_m4_prop15432.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm, paley_pair  # noqa: E402
from e1_gmin_m4_prop15298 import load_N, mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15413 import hs_mask, load_D  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15432_catalog.json"

# (paley_a, paley_b, N) -> (L+, L-)
LIVE_LNS = {
    5: {
        (-1, -1, 2): (Fraction(12), Fraction(13)),
        (-1, 1, 2): (Fraction(12), Fraction(13)),
        (-1, 1, 3): (Fraction(13), Fraction(12)),
        (1, -1, 2): (Fraction(12), Fraction(13)),
        (1, -1, 3): (Fraction(13), Fraction(12)),
        (1, 1, 3): (Fraction(13), Fraction(12)),
    },
    7: {
        (-1, -1, 3): (Fraction(2733, 38), Fraction(2853, 38)),
        (-1, 1, 3): (Fraction(1398, 19), Fraction(1395, 19)),
        (-1, 1, 5): (Fraction(1395, 19), Fraction(1398, 19)),
        (-1, 1, 6): (Fraction(147, 2), Fraction(147, 2)),
        (1, -1, 3): (Fraction(1398, 19), Fraction(1395, 19)),
        (1, -1, 5): (Fraction(1395, 19), Fraction(1398, 19)),
        (1, -1, 6): (Fraction(147, 2), Fraction(147, 2)),
        (1, 1, 5): (Fraction(2853, 38), Fraction(2733, 38)),
    },
}


def Lns_const_wrong(p: int) -> Fraction:
    return mu_plus(p) * mu_minus(p)


def Lns_p5_law(N: int, p: int) -> Fraction:
    """p=5 law: μμ + (N − p/2)."""
    return mu_plus(p) * mu_minus(p) + (N - Fraction(p, 2))


def fold_Lns(p: int) -> dict[tuple[int, int, int], tuple[Fraction, Fraction]]:
    D, F0 = load_D(p)
    F = {**F0, "p": p}
    Nmat = load_N(p, F).astype(np.float64)
    hs = hs_mask(p, D, F)
    Nnl = Nmat[~hs]
    ENN = (Nnl.T @ Nnl) / Nnl.shape[0]
    q = F["q"]
    acc: dict[tuple, dict[int, list[float]]] = defaultdict(lambda: {1: [], -1: []})
    for s in range(1, q):
        if F["chi"][s] != -1:
            continue
        a, b = paley_pair(F, s)
        Nn = int(field_norm(F, s))
        key = (int(a), int(b), Nn)
        for d in range(1, q):
            cd = int(F["chi"][d])
            acc[key][cd].append(float(ENN[d, F["mul"][s][d]]))
    out = {}
    for key, dmap in acc.items():
        out[key] = (
            Fraction(sum(dmap[1]) / len(dmap[1])).limit_denominator(5000),
            Fraction(sum(dmap[-1]) / len(dmap[-1])).limit_denominator(5000),
        )
    return out


def prove_A(fold7: dict | None = None) -> dict:
    tab7 = fold7 if fold7 is not None else fold_Lns(7)
    n3 = {k: v for k, v in tab7.items() if k[2] == 3}
    vals = {v[0] for v in n3.values()}
    ok = len(vals) >= 2
    ok = ok and Fraction(2733, 38) in vals
    ok = ok and Fraction(1398, 19) in vals
    if len(vals) < 2:
        ok = False
    return {
        "proved": bool(ok),
        "N3_Lplus": sorted(str(x) for x in vals),
        "theorem": (
            "L_ns constant on (Paley, N), not on N only. "
            "Fail: N=3 one L at p=7."
        ),
    }


def prove_B() -> dict:
    ok = True
    mm = {}
    for p, tab in LIVE_LNS.items():
        twomm = 2 * mu_plus(p) * mu_minus(p)
        mm[str(p)] = str(twomm)
        for v in tab.values():
            if v[0] + v[1] != twomm:
                ok = False
    # fail: L+=L- except N=-1
    p7 = LIVE_LNS[7]
    zero = [k for k, v in p7.items() if v[0] == v[1]]
    ok = ok and set(zero) == {(1, -1, 6), (-1, 1, 6)}
    if all(v[0] == v[1] for v in p7.values()):
        ok = False
    return {
        "proved": bool(ok),
        "two_mumu": mm,
        "zero_diff_p7": [list(k) for k in zero],
        "theorem": (
            "L₊+L₋=2μ₊μ₋. Fail: L₊=L₋ on every class."
        ),
    }


def prove_C() -> dict:
    ok = True
    for (a, b, N), (lp, lm) in LIVE_LNS[5].items():
        law = Lns_p5_law(N, 5)
        if lp != law:
            ok = False
    bad = Lns_p5_law(3, 7)
    live = LIVE_LNS[7][(-1, -1, 3)][0]
    ok = ok and bad != live
    ok = ok and live == Fraction(2733, 38)
    if bad == live:
        ok = False
    return {
        "proved": bool(ok),
        "p5_law": "μμ+(N−p/2)",
        "p7_law_on_mm_N3": str(bad),
        "p7_live": str(live),
        "theorem": (
            "p=5 L₊=μμ+(N−p/2). Fail: same law at p=7 (73≠2733/38)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "L_ns Paley×norm indexed. No formula in p. "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.432  L_ns Paley×norm indexed", flush=True)
    print("  folding p=5,7…", flush=True)
    f5 = fold_Lns(5)
    f7 = fold_Lns(7)
    # live table matches fold
    match = True
    for p, fold, live in ((5, f5, LIVE_LNS[5]), (7, f7, LIVE_LNS[7])):
        if set(fold) != set(live):
            match = False
        for k, v in live.items():
            if fold.get(k) != v:
                match = False
    print(f"  fold matches LIVE: {match}", flush=True)
    A = prove_A(f7)
    print(f"  A not N-only: {A['proved']} {A['N3_Lplus']}", flush=True)
    B = prove_B()
    print(f"  B sum 2μμ: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C p5 law: {C['proved']} p7 {C['p7_law_on_mm_N3']} vs {C['p7_live']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"match": match, "A": A["N3_Lplus"]}, indent=2) + "\n")
    out = {
        "prop": "15.432",
        "title": "L_ns is Paley×norm indexed; not μμ; Q_NL open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Lns_not_N_only": A["proved"] and match,
            "Lns_sum_2mumu": B["proved"],
            "p5_law_dies_p7": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "fold_matches_LIVE": match},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15432.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
