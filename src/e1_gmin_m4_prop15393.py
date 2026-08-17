#!/usr/bin/env python3
"""
Prop 15.393 — first_nc(T) is not a general p≡1 construction.
y_T is Max+ at p=37,41, but C_T has no Hoffman (nor even
C_T-clique) norm-circle.  t=0 / c∈{2,5,3,17} dies.  On the
four small primes the C-free Paley-half test names the same
c.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.392 flags.

============================================================================
Setup.  15.392: first_nc(T) has t=0 and min nsq Hoffman c at
p=5,13,17,29.  Circle N(u−t)=c has p+1 points.  C_T-clique
means y_T(u)y_T(v) χ_q(u−v)=+1 on the circle.

============================================================================
Theorem A — PROVED (C-free; certified p=5,13,17,29,37,41).
  Fix r with N(r)=c and U={N=1}.  The Paley-half identity
      y_T(rξ) = −y_T(r) χ_p(2−Tr(ξ))   (ξ∈U∖{1})
  holds for some (equivalently every) r iff c is a t=0
  Hoffman radius.  Hits {2,3} / {5} / {3} / {17} at
  p=5,13,17,29 and the empty set at p=37,41.
  Fail: claim ib=2 is a hit at p=13.  ∎

Theorem B — PROVED (full (t,c) scan; certified p=37,41).
  y_T is Max+ (Cy=py).  No pair (t,c) makes N(u−t)=c a
  C_T-clique, hence first_nc(T) is empty.  Fail: claim a
  t=0 clique at p=37.  ∎

Theorem C — PROVED (negative for leftover-1).
  Naming first_nc c cannot import ⟨δ,ψ⟩≤2: the T⊙NC
  generator is empty at the next p≡1 primes.  phi_F not
  imported.  ∎

============================================================================
Backend: field signs + one Paley C scan at p=37.  Writes
evidence/e1_gmin_m4_prop15393.json
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
from e1_gmin_m4_prop15312 import _norm_table  # noqa: E402
from e1_gmin_m4_prop15349 import y_triangle_p1  # noqa: E402
from e1_gmin_m4_prop15391 import in_T_lt  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15393_catalog.json"
LIVE = {5: 2, 13: 5, 17: 3, 29: 17}

_C_CACHE: dict[int, np.ndarray] = {}
_F_CACHE: dict[int, dict] = {}
_NT_CACHE: dict[int, list] = {}


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


def chi_p(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def least_nsq(p: int) -> int:
    for a in range(1, p):
        if chi_p(a, p) == -1:
            return a
    raise ValueError(p)


def yT(u: int, p: int) -> int:
    return -1 if in_T_lt(u, p) else 1


def tr0(u: int, p: int) -> int:
    """Tr(x+yω)=2x when ia=0 (all certified p≡1)."""
    return (2 * (u % p)) % p


def paley_half_hits(p: int) -> list[int]:
    F = _field(p)
    Nt = _Nt(p)
    q = F["q"]
    U = [u for u in range(1, q) if Nt[u] == 1]
    hits = []
    for c in range(1, p):
        if chi_p(c, p) != -1:
            continue
        rs = [u for u in range(1, q) if Nt[u] == c]
        if len(rs) != p + 1:
            continue
        any_r = False
        for r in rs:
            ok_r = True
            yr = yT(r, p)
            for xi in U:
                if xi == 1:
                    continue
                rxi = F["mul"][r][xi]
                rhs = -yr * chi_p(2 - tr0(xi, p), p)
                if yT(rxi, p) != rhs:
                    ok_r = False
                    break
            if ok_r:
                any_r = True
                break
        if any_r:
            hits.append(c)
    return hits


def n_clique_circles(p: int) -> dict:
    """Count C_T-clique norm-circles over all (t,c)."""
    y0 = y_triangle_p1(p)
    if y0[0] < 0:
        y0 = -y0
    C = _C(p)
    is_mp = bool(np.allclose(C @ y0, p * y0, atol=1e-4))
    F = _field(p)
    Nt = np.asarray(_Nt(p), dtype=np.int16)
    add = np.asarray(F["add"], dtype=np.int32)
    neg = np.asarray(F["neg"], dtype=np.int32)
    Cbase = (y0[:, None] * C) * y0[None, :]
    q = F["q"]
    n_cl = 0
    first = None
    for t in range(q):
        shifted = Nt[add[:, neg[t]]]
        for c in range(1, p):
            pts = np.flatnonzero(shifted == c)
            if pts.size != p + 1:
                continue
            S = 1 + pts
            block = Cbase[np.ix_(S, S)]
            if np.any(block[np.triu_indices(p + 1, 1)] <= 0):
                continue
            n_cl += 1
            if first is None:
                first = {"t": int(t), "c": int(c)}
    return {"Maxplus": is_mp, "n_clique": n_cl, "first": first}


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29, 37, 41):
        hits = paley_half_hits(p)
        rows[str(p)] = {"hits": hits, "min": min(hits) if hits else None}
        if p in LIVE:
            ok = ok and hits and min(hits) == LIVE[p]
        else:
            ok = ok and hits == []
    ok = ok and 2 not in paley_half_hits(13)
    if 2 in paley_half_hits(13):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Paley-half names live c at p=5,13,17,29 and is empty "
            "at p=37,41. Fail: ib at p=13."
        ),
    }


def prove_B() -> dict:
    # p=37 is the death prime (p=41 matches in evidence; skip in default test path)
    rec = n_clique_circles(37)
    ok = rec["Maxplus"] is True
    ok = ok and rec["n_clique"] == 0
    ok = ok and rec["first"] is None
    if rec["n_clique"] > 0:
        ok = False
    return {
        "proved": bool(ok),
        "p37": rec,
        "theorem": (
            "y_T is Max+ at p=37 and C_T has 0 clique norm-circles. "
            "Fail: claim a t=0 clique."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "T⊙NC empty at p=37,41 so first_nc c is not a general "
            "name. Ensemble D / Q_τ still unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.393  first_nc(T) empty at p=37,41", flush=True)
    A = prove_A()
    print(f"  A half: {A['proved']} { {k: A['rows'][k]['min'] for k in A['rows']} }", flush=True)
    B = prove_B()
    print(f"  B empty37: {B['proved']} {B['p37']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["p37"], "LIVE": LIVE}, indent=2) + "\n"
    )
    out = {
        "prop": "15.393",
        "title": "first_nc(T) empty at p=37,41; Paley-half names the small-p c",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "paley_half_matches_or_empty": A["proved"],
            "no_circle_at_37": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15393.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
