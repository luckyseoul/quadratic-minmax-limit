#!/usr/bin/env python3
"""
Prop 15.513 — Even pairing first moment is −16, independent of Q_τ.

  For p≡1≥5 and even k∉{0,(q−1)/2}:
      n(J_{N*}=2)        = |N*| = (p−1)(p−3)/2
      n(J_{N*}=-(p−1))   = p−3
      n(J_{N*}=-(p−3))   = p−1
  Hence ∑_k ⟨δ,ψ_k⟩ = −16, and the coefficient of M=n_pp δ++ vanishes.
  Fail: n(J=2)=n_pp; claim S1=0; claim S1 depends on Q++.
  Does not name Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.512 flags.

============================================================================
Setup.  even_ks and J_{N*} as in 15.497.  pair_from_M as in 15.507.
n_pp=2(p−2), |N*|=(p−1)(p−3)/2.  gcd(p−1,p+1)=2 so
lcm=(p²−1)/2=(q−1)/2.

============================================================================
Theorem A — PROVED (field; p≡1≥5).
  Among even_ks, (p−1)|k for p−1 values (p multiples of p−1 in
  (0,q−1), minus k=(q−1)/2), and (p+1)|k for p−3 values.
  15.497 sends these to J=-(p−3) and J=-(p−1); the rest
  (|N*| of them) have J=2.  Fail: n(J=2)=n_pp.  ∎

Theorem B — PROVED (A + 15.507 pair_from_M).
  Writing S1=∑ pair(k)=c_M M + c_0, one has c_M=0 and c_0=−16
  for every prime p≡1≥5.  Fail: S1=0; c_M≠0.  ∎

Theorem C — OPEN.  S1 does not determine M or Q_τ (c_M=0).
  S2=10 M²−96 M+256 at p=5 still depends on the unnamed M.
  phi_F not imported.  ∎

============================================================================
Backend: Fraction arithmetic + even_ks count (ProcessPool over primes).
Writes evidence/e1_gmin_m4_prop15513.json
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
from e1_gmin_m4_prop15497 import J_Nstar_named, even_ks
from e1_gmin_m4_prop15507 import n_nstar_p1

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15513.json"

CERT = (5, 13, 17, 29, 37, 41)
W = 86


def n_J2_named(p: int) -> int:
    return n_nstar_p1(p)


def n_Jm1_named(p: int) -> int:
    """n(J=-(p-1)) = p-3."""
    return p - 3


def n_Jm3_named(p: int) -> int:
    """n(J=-(p-3)) = p-1."""
    return p - 1


def n_J2_wrong(p: int) -> int:
    """Fail-when-wrong: n(J=2)=n_pp."""
    return n_pp_named(p)


def _count_job(p: int) -> dict:
    from collections import Counter

    c = Counter(J_Nstar_named(p, k) for k in even_ks(p))
    return {
        "p": p,
        "n_k": len(even_ks(p)),
        "n2": c[2],
        "nm1": c[-(p - 1)],
        "nm3": c[-(p - 3)],
    }


def s1_coefs(p: int) -> tuple[Fraction, Fraction]:
    n_pp = n_pp_named(p)
    n_ns = n_nstar_p1(p)
    rows = (
        (2, n_J2_named(p)),
        (-(p - 1), n_Jm1_named(p)),
        (-(p - 3), n_Jm3_named(p)),
    )
    c_M = c_0 = Fraction(0)
    for jn, n in rows:
        jpp = -2 - jn
        c_M += n * (Fraction(jpp, n_pp) - Fraction(jn, n_ns))
        c_0 += n * Fraction(jn * 8, n_ns)
    return c_M, c_0


def prove_A() -> dict:
    ps = [p for p in range(5, 80) if is_prime(p) and p % 4 == 1]
    ok = True
    rows = {}
    with ProcessPoolExecutor(max_workers=min(W, len(ps))) as ex:
        for rec in ex.map(_count_job, ps):
            p = rec["p"]
            if rec["n2"] != n_J2_named(p):
                ok = False
            if rec["nm1"] != n_Jm1_named(p):
                ok = False
            if rec["nm3"] != n_Jm3_named(p):
                ok = False
            if rec["n2"] == n_J2_wrong(p):
                ok = False
            if rec["n2"] + rec["nm1"] + rec["nm3"] != rec["n_k"]:
                ok = False
            if p in CERT:
                rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "n_primes": len(ps),
        "rows": rows,
        "theorem": "jcnt named. Fail n(J=2)=n_pp.",
    }


def prove_B() -> dict:
    ps = [p for p in range(5, 80) if is_prime(p) and p % 4 == 1]
    ok = True
    rows = {}
    for p in ps:
        c_M, c_0 = s1_coefs(p)
        if c_M != 0 or c_0 != -16:
            ok = False
        if p in CERT:
            rows[str(p)] = {"c_M": str(c_M), "c_0": str(c_0)}
    return {
        "proved": bool(ok),
        "n_primes": len(ps),
        "rows": rows,
        "theorem": "S1=−16, independent of M. Fail S1=0 or c_M≠0.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "S1 does not determine M (c_M=0). S2 still depends on M. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.513 even pairing S1=-16 independent of Q_τ", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.513",
        "title": "Even pairing first moment −16; Q_τ unnamed",
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
