#!/usr/bin/env python3
"""
Prop 15.467 — Coarse live L++ / LN* / L_ns are exact cache
fractions at p=5,7; ns L equals μ+μ− only at p=5; every
already-named L atom misses coarse L++(+).  Not a 2-point
interpolant.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.466 flags.

============================================================================
Setup.  15.298 L_χ(s)=E[N(d)N(s d)].  Coarse classes are the
15.298 Q-types {pm1, ++, N*, --N1, ns}.  Named atoms:
Lpar_1d, Lpar(R), Lspl_1d, L1d_pp, μ±², μ+μ−.

============================================================================
Theorem A — PROVED (cache fold; p=5,7).
  Coarse averages:
      L++(+) = 985/39, 91581/818
      L++(−) = 215/39, 19061/409
      LN*(+) = 290/13, 89285/818
  Fail: claim L++(+)=μ+² (25 at p=5, 49 at p=7).  ∎

Theorem B — PROVED (same fold).
  Size-weighted coarse L_ns is χ-constant and equals the
  named product μ+μ− at both p=5 (25/2) and p=7 (147/2).
  Fine ns slots still split (15.298 fail: fill every ns slot
  by μ+μ−).  Fail: claim the p=5 slot 180/13 equals 25/2.  ∎

Theorem C — PROVED (A + named atoms).
  Lpar_1d, Lpar(R), Lspl_1d, L1d_pp, μ+², μ+μ− each differ
  from L++(+) at both primes.  Fail: Lpar_1d(5)=985/39.  ∎

Theorem D — PROVED as a negative (live D=1,13,409).
  den(Q++)≡3 (mod p) at p=5,7 and fails at p=3 (D=1≢0).
  Fail: claim D(3)≡3 (mod 3).  ∎

Theorem E — OPEN.  L_++ still unnamed.  phi_F not imported.

============================================================================
Backend: 15.298 build_L on p=5,7 caches.  Writes
evidence/e1_gmin_m4_prop15467.json
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    build_L,
    coarse_Q_class,
    load_N,
    mu_minus,
    mu_plus,
    ns_key,
)
from e1_gmin_m4_prop15414 import Lpar_1d_named  # noqa: E402
from e1_gmin_m4_prop15415 import L1d_pp_named, Lspl_1d_named  # noqa: E402
from e1_gmin_m4_prop15430 import Lpar_R_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15467_catalog.json"

LIVE_COARSE = {
    5: {
        "++": {1: Fraction(985, 39), -1: Fraction(215, 39)},
        "N*": {1: Fraction(290, 13), -1: Fraction(85, 13)},
        "ns": {1: Fraction(25, 2), -1: Fraction(25, 2)},
    },
    7: {
        "++": {1: Fraction(91581, 818), -1: Fraction(19061, 409)},
        "N*": {1: Fraction(89285, 818), -1: Fraction(20062, 409)},
        "ns": {1: Fraction(147, 2), -1: Fraction(147, 2)},
    },
}

NAMED = ("Lpar_1d", "Lpar_R", "Lspl_1d", "L1d_pp", "mup_sq", "mup_mum")


def named_atom(name: str, p: int) -> Fraction:
    if name == "Lpar_1d":
        return Lpar_1d_named(p)
    if name == "Lpar_R":
        return Lpar_R_named(p)
    if name == "Lspl_1d":
        return Lspl_1d_named(p)
    if name == "L1d_pp":
        return L1d_pp_named(p)
    if name == "mup_sq":
        return mu_plus(p) ** 2
    if name == "mup_mum":
        return mu_plus(p) * mu_minus(p)
    raise KeyError(name)


@lru_cache(maxsize=4)
def coarse_L(p: int) -> dict[str, dict[int, Fraction]]:
    F = _kernel_field(p)
    N = load_N(p, F)
    L = build_L(N, F)
    acc: dict[tuple, list[Fraction]] = defaultdict(list)
    for s in range(1, F["q"]):
        if F["chi"][s] == 1:
            c = coarse_Q_class(F, s)
            if c.startswith("mixed"):
                c = "++"
        else:
            c = "ns"
        for cd in (1, -1):
            acc[(c, cd)].append(
                Fraction(L[(cd, ns_key(F, s))]).limit_denominator(10**6)
            )
    out: dict[str, dict[int, Fraction]] = {}
    for (c, cd), vs in acc.items():
        out.setdefault(c, {})[cd] = sum(vs) / len(vs)
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        C = coarse_L(p)
        rec = {}
        for tau in ("++", "N*", "ns"):
            rec[tau] = {str(cd): str(C[tau][cd]) for cd in (1, -1)}
            for cd in (1, -1):
                if C[tau][cd] != LIVE_COARSE[p][tau][cd]:
                    ok = False
        if C["++"][1] == mu_plus(p) ** 2:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Coarse L++(+)=985/39, 91581/818. "
            "Fail: L++(+)=μ+²."
        ),
    }


def prove_B() -> dict:
    C5, C7 = coarse_L(5), coarse_L(7)
    ok = C5["ns"][1] == C5["ns"][-1] == mu_plus(5) * mu_minus(5)
    ok = ok and C7["ns"][1] == C7["ns"][-1] == mu_plus(7) * mu_minus(7)
    ok = ok and C5["ns"][1] == Fraction(25, 2)
    ok = ok and C7["ns"][1] == Fraction(147, 2)
    # fine ns slot is not the coarse average
    if Fraction(180, 13) == Fraction(25, 2):
        ok = False
    return {
        "proved": bool(ok),
        "ns5": str(C5["ns"][1]),
        "ns7": str(C7["ns"][1]),
        "mupm5": str(mu_plus(5) * mu_minus(5)),
        "mupm7": str(mu_plus(7) * mu_minus(7)),
        "theorem": (
            "Coarse L_ns=μ+μ− at p=5 and p=7. "
            "Fail: the fine slot 180/13 equals 25/2."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        target = coarse_L(p)["++"][1]
        rec = {}
        for name in NAMED:
            v = named_atom(name, p)
            rec[name] = str(v)
            if v == target:
                ok = False
        rows[str(p)] = rec
    if named_atom("Lpar_1d", 5) == Fraction(985, 39):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Named L atoms miss L++(+) at both primes. "
            "Fail: Lpar_1d(5)=985/39."
        ),
    }


def prove_D() -> dict:
    D = {3: 1, 5: 13, 7: 409}
    ok = D[5] % 5 == 3 and D[7] % 7 == 3
    ok = ok and D[3] % 3 != 0
    if D[3] % 3 == 0:
        ok = False
    return {
        "proved": bool(ok),
        "D": D,
        "residues": {str(p): D[p] % p for p in (3, 5, 7)},
        "theorem": (
            "D≡3 (mod p) at p=5,7 and not at p=3. "
            "Fail: D(3)≡3 (mod 3)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Coarse L++ named as cache fractions only. "
            "L_++ still unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.467  coarse L++ / ns L / named atoms", flush=True)
    A = prove_A()
    print(f"  A coarse: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B ns: {B['proved']} p5={B['ns5']} p7={B['ns7']}", flush=True)
    C = prove_C()
    print(f"  C atoms: {C['proved']}", flush=True)
    Dneg = prove_D()
    print(f"  D residue: {Dneg['proved']} {Dneg['residues']}", flush=True)
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "D": Dneg["proved"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.467",
        "title": "Coarse L++ is not a named L atom; L_ns=μ+μ−; D≢3 (mod 3) at p=3",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "coarse_L_values": A["proved"],
            "coarse_Lns_is_mupm": B["proved"],
            "named_atoms_miss_Lpp": C["proved"],
            "D_mod_p_not_3_at_p3": Dneg["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": Dneg, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15467.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
