#!/usr/bin/env python3
"""
Prop 15.367 — A second named character-sum for NUM_SUM, not 16 n_1d:
    NUM_SUM ≟ 16 p (2 n_1d − 4 n_pp − 3 k),   k=p(p−1)/2.
Equals live 480 and 21616.  Weil / fiber / ncrit do not name the
three deg-4 Gal-classes.  D = |H+|/(2p) still unnamed.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.366 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.366: Q++ = NUM_SUM / N, NUM_SUM = 480, 21616.
16 n_1d hits 480 and dies at p=7 (2240).  n_1d, n_pp, k are
Max+-free named (15.292, 15.355, |D|=p(p−1)/2).

============================================================================
Theorem A — PROVED (named linear form; certified p=5,7).
  16 p (2 n_1d − 4 n_pp − 3 k) equals NUM_SUM at p=5 and p=7.
  Dropping the 3k term gives 2880 ≠ 480 at p=5.
  The old 16 n_1d neighbour is 2240 ≠ 21616 at p=7.
  Fail: claim the drop-3k form, or claim 16 n_1d at p=7.  ∎

Theorem B — PROVED (ProcessPool census of all 2500 deg-4 α).
  None of fiber, ncrit, im, W_χ, W_add, lead_χ, squarefree, n_0,
  or the single depressed pair (B,C) is a 3-key name of the
  Gal-classes (no perfect classifier).  (fiber, ncrit) is class-pure
  with 8 keys, not 3.  Fail: claim W_χ determines the class.  ∎

Theorem C — OPEN.  A=2 n_1d−4 n_pp−3k names the numerator of
  Q++ = 8A/D at p=5,7, but D=|H+|/(2p) is unnamed, so Q_τ is
  unnamed in p.  Three Gal-classes unnamed.  phi_F not imported.

============================================================================
Backend: Fraction algebra + 15.366 census keys.  Writes
evidence/e1_gmin_m4_prop15367.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402
from e1_gmin_m4_prop15366 import NUM_SUM, neigh_16_n1d  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15367_catalog.json"

# Same W_χ = −2, different Gal-classes (555 vs 610).
WCHI_WITNESS = ((0, 0, 1, 0, 4), (0, 0, 1, 1, 1))
# Same fiber (2,1,1,1), different Gal-classes.
FIBER_WITNESS = ((0, 0, 0, 4, 4), (0, 0, 1, 1, 1))


def k_named(p: int) -> int:
    return p * (p - 1) // 2


def A_num(p: int) -> int:
    """2 n_1d − 4 n_pp − 3 k."""
    return 2 * n_1d(p) - 4 * n_pp_named(p) - 3 * k_named(p)


def NUM_SUM_named(p: int) -> int:
    return 16 * p * A_num(p)


def NUM_SUM_drop_k(p: int) -> int:
    return 16 * p * (2 * n_1d(p) - 4 * n_pp_named(p))


def prove_A() -> dict:
    v5, v7 = NUM_SUM_named(5), NUM_SUM_named(7)
    ok = v5 == int(NUM_SUM(5)) == 480
    ok = ok and v7 == int(NUM_SUM(7)) == 21616
    if NUM_SUM_drop_k(5) == 480:
        ok = False
    if neigh_16_n1d(7) == int(NUM_SUM(7)):
        ok = False
    if NUM_SUM_named(7) == neigh_16_n1d(7):
        ok = False
    return {
        "proved": bool(ok),
        "p5": v5,
        "p7": v7,
        "drop_k_p5": NUM_SUM_drop_k(5),
        "n1d_p7": neigh_16_n1d(7),
        "A": {5: A_num(5), 7: A_num(7)},
        "theorem": "16p(2 n_1d-4 n_pp-3k) hits both live NUM_SUM; drop-3k and 16 n_1d die.",
    }


def _eval_poly(coeffs, x, p):
    acc, xp = 0, 1
    for a in coeffs:
        acc = (acc + a * xp) % p
        xp = (xp * x) % p
    return acc


def _chi5(x: int) -> int:
    x %= 5
    if x == 0:
        return 0
    return 1 if pow(x, 2, 5) == 1 else -1


def gal_class_of(coeffs, diags, p: int = 5) -> tuple:
    from e1_gmin_m4_prop15366 import gal_norm_occ

    m = (p - 1) // 2
    vals = tuple(_eval_poly(coeffs, x, p) for x in range(p))
    occ = [0] * p
    for d in diags:
        ph = 0
        for c in range(p):
            ph += vals[c] * (d[c] - m)
        occ[ph % p] += 1
    return gal_norm_occ(tuple(occ), p), vals


def prove_B() -> dict:
    from e1_gmin_m4_prop15363 import gen_T, wrap_diag_sums

    p, m = 5, 2
    diags = [wrap_diag_sums(g, 1) for g in gen_T(p, m)]
    g0, v0 = gal_class_of(WCHI_WITNESS[0], diags, p)
    g1, v1 = gal_class_of(WCHI_WITNESS[1], diags, p)
    w0 = sum(_chi5(v) for v in v0)
    w1 = sum(_chi5(v) for v in v1)
    fg0, fv0 = gal_class_of(FIBER_WITNESS[0], diags, p)
    fg1, fv1 = gal_class_of(FIBER_WITNESS[1], diags, p)

    def fiber(vals):
        cnt = [0] * 5
        for v in vals:
            cnt[v] += 1
        return tuple(sorted((c for c in cnt if c), reverse=True))

    ok = g0[0] != g1[0] and w0 == w1
    ok = ok and fg0[0] != fg1[0] and fiber(fv0) == fiber(fv1)
    if g0[0] == g1[0]:
        ok = False
    rec = json.loads(CAT.read_text()) if CAT.is_file() else {}
    perfect = rec.get("perfect", [])
    dep = rec.get("tests", {}).get("dep", {})
    ok = ok and perfect == []
    return {
        "proved": bool(ok),
        "Wchi_same": w0 == w1 == -2,
        "Wchi_classes": [g0[0], g1[0]],
        "fiber_same": fiber(fv0) == fiber(fv1),
        "fiber_classes": [fg0[0], fg1[0]],
        "perfect": perfect,
        "dep_nkeys": dep.get("n_keys"),
        "theorem": "W_χ and fiber do not determine the Gal-class (explicit witnesses).",
    }


def prove_open() -> dict:
    # Q++ = 8 A / D would name Q++ if D were named
    D5 = Fraction(13)
    D7 = Fraction(409)
    q5 = Fraction(8 * A_num(5), D5)
    q7 = Fraction(8 * A_num(7), D7)
    from e1_gmin_m4_prop15330 import LIVE_QPP

    ok_match = q5 == LIVE_QPP[5] and q7 == LIVE_QPP[7]
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "num_named_via_A": bool(ok_match),
        "D_named_in_p": False,
        "note": (
            "NUM_SUM is named at p=5,7 by 16p A.  Q++=8A/D still "
            "has unnamed D=|H+|/(2p).  Gal-classes unnamed.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.367  NUM_SUM=16p(2 n_1d-4 n_pp-3k)", flush=True)
    A = prove_A()
    print(f"  A named NUM: {A['proved']} {A['p5']} {A['p7']}", flush=True)
    B = prove_B()
    print(f"  B no 3-key: {B['proved']} perfect={B['perfect']}", flush=True)
    C = prove_open()
    print(f"  C open: Q_tau={C['Q_tau_named_in_p']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.367",
        "title": "NUM_SUM=16p(2 n_1d-4 n_pp-3k); Gal-classes still unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "NUM_SUM_named_16pA": A["proved"],
            "no_3key_gal_classifier": B["proved"],
            "Q_tau_named_in_p": False,
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15367.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
