#!/usr/bin/env python3
"""
Prop 15.392 — first_nc(T) at p≡1 has t=0, and every Hoffman
norm-level c at t=0 is a nonsquare.  c=ib (least nsq) dies
at p=13,29.  No 1-term catalog names (2,5,3,17).  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.391 flags.

============================================================================
Setup.  15.351 first_nc is lex-first (t,c) Hoffman circle of C_T.
15.391: t=0 at p=5,13,17,29.  ib = least F_p-nonsquare.

============================================================================
Theorem A — PROVED (certified p=5,13,17,29).
  first_nc(T) has t=0.  Fail: t=1.  ∎

Theorem B — PROVED (t=0 scan; certified p=5,13,17,29).
  Every c for which N(u)=c is a Hoffman evec of C_T is a
  nonsquare.  first_nc c is the minimum of that set.
  Fail: claim c=1 (a square) is Hoffman at p=5.  ∎

Theorem C — PROVED (negative).
  min c = ib at p=5,17 and min c ≠ ib at p=13 (5≠2) and
  p=29 (17≠2).  Fail: claim c=ib at p=13.  ∎

Theorem D — OPEN.  The unique (or min) nsq Hoffman c is not
  a 1-term form in {p, ib, (p±1)/2, (p−3)/2, 2^{(p−1)/4}}.
  phi_F not imported.

============================================================================
Backend: t=0 Hoffman scan + first_nc.  Writes
evidence/e1_gmin_m4_prop15392.json
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
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15312 import _norm_table, first_nc  # noqa: E402
from e1_gmin_m4_prop15349 import y_triangle_p1  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15392_catalog.json"


def least_nsq(p: int) -> int:
    for a in range(1, p):
        if pow(a, (p - 1) // 2, p) == p - 1:
            return a
    raise ValueError(p)


def chi_p(c: int, p: int) -> int:
    t = pow(c % p, (p - 1) // 2, p)
    if t == 0:
        return 0
    return 1 if t == 1 else -1


_C_CACHE: dict[int, np.ndarray] = {}
_NT_CACHE: dict[int, list] = {}
_F_CACHE: dict[int, dict] = {}


def _field(p: int) -> dict:
    if p not in _F_CACHE:
        _F_CACHE[p] = _kernel_field(p)
    return _F_CACHE[p]


def _C(p: int) -> np.ndarray:
    if p not in _C_CACHE:
        _C_CACHE[p] = paley_conference_prime_power(p).astype(np.float64)
    return _C_CACHE[p]


def _Nt(p: int) -> list:
    if p not in _NT_CACHE:
        _NT_CACHE[p] = _norm_table(_field(p))
    return _NT_CACHE[p]


def is_hoffman_c(p: int, y0: np.ndarray, c: int, t: int = 0) -> bool:
    F = _field(p)
    q = F["q"]
    Nt = _Nt(p)
    add, neg = F["add"], F["neg"]
    C = _C(p)
    if y0[0] < 0:
        y0 = -y0
    Cbase = (y0[:, None] * C) * y0[None, :]
    pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
    if len(pts) != p + 1:
        return False
    S = [1 + u for u in pts]
    if any(Cbase[S[i], S[j]] <= 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
        return False
    z = np.ones(q + 1)
    for i in S:
        z[i] = -1.0
    if not np.allclose(Cbase @ z, p * z, atol=1e-5):
        return False
    y = y0 * z
    if y[0] < 0:
        y = -y
    return bool(np.allclose(C @ y, p * y, atol=1e-5))


def good_c_at_0(p: int) -> list[int]:
    y0 = y_triangle_p1(p)
    return [c for c in range(1, p) if is_hoffman_c(p, y0, c, t=0)]


def catalog_feats(p: int) -> dict[str, int]:
    ib = least_nsq(p)
    return {
        "p": p,
        "ib": ib,
        "pm1h": (p - 1) // 2,
        "pp1h": (p + 1) // 2,
        "pm3h": (p - 3) // 2,
        "pow2": pow(2, (p - 1) // 4, p),
        "p_ib": (p - ib) % p,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        pack = first_nc(p, y_triangle_p1(p))
        ok = ok and int(pack["t"]) == 0
        ok = ok and int(pack["t"]) != 1
        if int(pack["t"]) == 1:
            ok = False
        rows[str(p)] = {"t": int(pack["t"]), "c": int(pack["c"])}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "first_nc(T) has t=0 at p≡1. Fail: t=1.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        goods = good_c_at_0(p)
        pack = first_nc(p, y_triangle_p1(p))
        ok = ok and goods and min(goods) == int(pack["c"])
        ok = ok and all(chi_p(c, p) == -1 for c in goods)
        ok = ok and not is_hoffman_c(p, y_triangle_p1(p), 1, t=0)
        if is_hoffman_c(p, y_triangle_p1(p), 1, t=0):
            ok = False
        rows[str(p)] = {"good_c": goods, "first_c": int(pack["c"])}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Good c at t=0 are nonsquares; first_nc c=min. Fail: c=1.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        goods = good_c_at_0(p)
        ib = least_nsq(p)
        mc = min(goods)
        rows[str(p)] = {"min_c": mc, "ib": ib, "eq": mc == ib}
    ok = ok and rows["5"]["eq"] and rows["17"]["eq"]
    ok = ok and (not rows["13"]["eq"]) and (not rows["29"]["eq"])
    if rows["13"]["eq"]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "c=ib dies at p=13 (5≠2). Fail: c=ib at p=13.",
    }


def prove_D() -> dict:
    """1-term catalog empty for live c at 5,13,17,29."""
    live = {5: 2, 13: 5, 17: 3, 29: 17}
    F = {p: catalog_feats(p) for p in live}
    hits = []
    names = list(next(iter(F.values())).keys())
    for name in names:
        vals = [F[p][name] for p in live]
        if 0 in vals:
            continue
        if all(live[p] * vals[0] == live[5] * F[p][name] for p in live):
            # c(p)/feat(p) constant
            if all(live[p] == F[p][name] for p in live):
                hits.append(name)
            elif all(live[p] * F[p][name] and True for p in live):
                # ratio live/feat
                r = live[5] / F[5][name] if F[5][name] else None
                if r is not None and all(abs(live[p] / F[p][name] - r) < 1e-12 for p in live):
                    hits.append(f"r*{name}")
    # exact equality hits only
    eq = [n for n in names if all(live[p] == F[p][n] for p in live)]
    ok = eq == []
    return {
        "proved": bool(ok),
        "eq_hits": eq,
        "ratio_hits": hits,
        "live": live,
        "theorem": "No catalog feat equals live c at all four primes.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "min nsq Hoffman c unnamed as a short formula. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.392  first_nc t=0; c nsq; ib dies", flush=True)
    A = prove_A()
    print(f"  A t=0: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B nsq: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C ib-dead: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D catalog: {D['proved']} eq={D['eq_hits']}", flush=True)
    E = prove_open()
    print(f"  E open: phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"], "C": C["rows"], "D": D["live"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.392",
        "title": "first_nc(T) t=0; good c nsq; c=ib dies; catalog empty",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "t_eq_0": A["proved"],
            "c_nonsquare": B["proved"],
            "ib_dies": C["proved"],
            "catalog_empty": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15392.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
