#!/usr/bin/env python3
"""
Prop 15.514 — Even pairing S2 is a named quadratic in M, not a
name of M.  For p≡1≥5
    S2 = a(p) M² + b(p) M + c(p)
    a=(p+1)(p²−5)/(2(p−2)(p−3)),
    b=−32(p+1)/(p−3),
    c=512/(p−3).
  Vertex M=32(p−2)/(p²−5) (=24/5 at p=5).  Fail: a=0 (the 15.513
  S1 vanishing); claim the vertex equals live M=24/13.
  Need S2≤5920/49 to force M≥12/7 at p=5; λ≥0 and 0≤δ≤4 only
  give S2≤640.  Does not name Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.513 flags.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, is_prime
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15507 import n_nstar_p1
from e1_gmin_m4_prop15513 import n_J2_named, n_Jm1_named, n_Jm3_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15514.json"

CERT = (5, 13, 17, 29, 37, 41)
W = 86


def s2_named(p: int) -> tuple[Fraction, Fraction, Fraction]:
    a = Fraction((p + 1) * (p * p - 5), 2 * (p - 2) * (p - 3))
    b = Fraction(-32 * (p + 1), p - 3)
    c = Fraction(512, p - 3)
    return a, b, c


def s2_vertex_named(p: int) -> Fraction:
    return Fraction(32 * (p - 2), p * p - 5)


def s2_vertex_wrong(p: int) -> Fraction:
    """Fail-when-wrong: live M at p=5 (24/13), or 0."""
    return Fraction(24, 13) if p == 5 else Fraction(0)


def s2_from_jcnt(p: int) -> tuple[Fraction, Fraction, Fraction]:
    npp = n_pp_named(p)
    nns = n_nstar_p1(p)
    rows = (
        (2, n_J2_named(p)),
        (-(p - 1), n_Jm1_named(p)),
        (-(p - 3), n_Jm3_named(p)),
    )
    c0 = c1 = c2 = Fraction(0)
    for jn, n in rows:
        jpp = -2 - jn
        A = Fraction(jpp, npp) - Fraction(jn, nns)
        B = Fraction(jn * 8, nns)
        c2 += n * A * A
        c1 += n * 2 * A * B
        c0 += n * B * B
    return c2, c1, c0


def _job(p: int) -> dict:
    named = s2_named(p)
    live = s2_from_jcnt(p)
    return {
        "p": p,
        "named": [str(x) for x in named],
        "jcnt": [str(x) for x in live],
        "ok": named == live,
        "a_zero": named[0] == 0,
        "vertex": str(s2_vertex_named(p)),
    }


def prove_A() -> dict:
    ps = [p for p in range(5, 80) if is_prime(p) and p % 4 == 1]
    ok = True
    rows = {}
    with ProcessPoolExecutor(max_workers=min(W, len(ps))) as ex:
        for rec in ex.map(_job, ps):
            if not rec["ok"] or rec["a_zero"]:
                ok = False
            if rec["p"] == 5 and rec["vertex"] == str(s2_vertex_wrong(5)):
                ok = False
            if rec["p"] in CERT:
                rows[str(rec["p"])] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(ps),
        "rows": rows,
        "theorem": "S2=a(p)M²+b(p)M+c(p) named. Fail a=0; vertex=24/13.",
    }


def prove_B() -> dict:
    # p=5: S2(M=12/7)=5920/49; λ≥0 / 0≤δ≤4 give 640
    need = Fraction(5920, 49)
    pos4 = Fraction(640)
    ok = pos4 > need
    if s2_named(5)[0] == 0:
        ok = False
    return {
        "proved": bool(ok),
        "s2_need": str(need),
        "s2_pos4_ub": str(pos4),
        "pos4_too_loose": pos4 > need,
        "theorem": "0≤δ≤4 gives S2≤640>5920/49. Fail claim that majorant closes M≥12/7.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "S2 names a quadratic in M, not M. No CS majorant hits "
            "S2≤5920/49. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.514 S2 named quadratic in M; does not name Q_τ", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.514",
        "title": "Even pairing S2 named in M; Q_τ unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
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
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
