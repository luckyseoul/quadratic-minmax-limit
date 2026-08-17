#!/usr/bin/env python3
"""
Prop 15.484 — TR L++ is not a formula in p.  The 2-point interpolant
    L = μ+² − χ_p(−1) + (1−χ_p(−1))/12
hits 24 and 1337/12 and dies on constructed p=11 TR Max+
(Cy=py).  Raw T is not the TR 4-point.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.483 flags.

============================================================================
Setup.  15.483: live-cache TR L++ is the single value 24 (p=5) and
1337/12 (p=7).  15.477: AGL(2,p)·T ∩ Max+ = TR.

============================================================================
Theorem A — PROVED (algebra + cache; p=5,7).
  L_interp(p)=μ+²−χ_p(−1)+(1−χ_p(−1))/12 equals live TR L++ at
  p=5,7.  Fail: drop the 1/12 term (445/4 ≠ 1337/12 at p=7).  ∎

Theorem B — PROVED (AGL construction; Cy=py).
  AGL(2,11)·T contains Max+ TR rows with at least two distinct
  L++ values, neither equal to L_interp(11)=9089/12.
  Fail: claim L_interp(11) is TR L++.  ∎

Theorem C — PROVED (same construction + cache).
  TR L++ is not a function of p: one value at p=5 and p=7 (full
  live cache), ≥2 values among Cy=py AGL images at p=11.
  Fail: claim a single TR L++ at p=11.  ∎

Theorem D — PROVED (triangle mask).
  Raw T has L++=18 at p=5 ≠ 24.  Fail: claim raw T is the TR
  4-point.  ∎

Theorem E — OPEN.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.477 triangle + Paley conference Cy=py.  Writes
evidence/e1_gmin_m4_prop15484.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15476 import is_triple_rainbow, live_delta2_and_sigs  # noqa: E402
from e1_gmin_m4_prop15477 import triangle_mask  # noqa: E402
from e1_gmin_m4_prop15478 import I_pp  # noqa: E402
from e1_gmin_m4_prop15483 import fold  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15484_catalog.json"


def chi_m1(p: int) -> int:
    return 1 if p % 4 == 1 else -1


def L_tr_interp(p: int) -> Fraction:
    mu = Fraction(p * (p - 1), 4)
    ch = chi_m1(p)
    return mu * mu - ch + Fraction(1 - ch, 12)


def L_tr_drop12(p: int) -> Fraction:
    mu = Fraction(p * (p - 1), 4)
    return mu * mu - chi_m1(p)


def L_on_mask(F, mask: np.ndarray) -> Fraction:
    add = np.array(F["add"], dtype=np.int32)
    D = mask.astype(bool)
    N = np.array([(D & D[add[:, d]]).sum() for d in range(F["q"])], dtype=np.int64)
    s0 = next(r for r in range(1, F["q"]) if I_pp(F, r) > 0.5)
    plus = [d for d in range(1, F["q"]) if F["chi"][d] == 1]
    tot = int(sum(int(N[d]) * int(N[F["mul"][s0][d]]) for d in plus))
    return Fraction(tot, len(plus))


def live_cache_L_tr(p: int) -> Fraction:
    F = _kernel_field(p)
    _d2, sigs, _ = live_delta2_and_sigs(p)
    N = load_N(p, F)
    s0 = next(r for r in range(1, F["q"]) if I_pp(F, r) > 0.5)
    plus = [d for d in range(1, F["q"]) if F["chi"][d] == 1]
    sd = [F["mul"][s0][d] for d in plus]
    tr = [i for i, s in enumerate(sigs) if is_triple_rainbow(s, p)]
    tot = int((N[np.ix_(tr, plus)] * N[np.ix_(tr, sd)]).sum())
    return Fraction(tot, len(tr) * len(plus))


@lru_cache(None)
def agl_tr_maxplus_L(p: int, max_tries: int = 200, seed: int = 1, want: int = 4):
    """L++ of AGL(T) images with Cy=py (true Max+)."""
    from e1_gmin_m4_prop15473 import e_matrix
    from e1_gmin_m4_prop15474 import sigma_named
    from e1_gmin_m4_prop15476 import line_blocks
    from minmax_quadratic import paley_conference_prime_power

    rng = np.random.default_rng(seed)
    F = _kernel_field(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    pts = [(x, y) for x in range(p) for y in range(x)]
    EM = e_matrix(F)
    chi = np.array(F["chi"], dtype=np.int32)
    sigma = sigma_named(p)
    Omega = chi == sigma
    Omega[0] = False
    blocks = line_blocks(p, q)
    out = []
    for _t in range(max_tries):
        while True:
            A = rng.integers(0, p, size=(2, 2))
            if (int(A[0, 0]) * int(A[1, 1]) - int(A[0, 1]) * int(A[1, 0])) % p:
                break
        b = rng.integers(0, p, size=2)
        m = np.zeros(q, dtype=np.float64)
        for x, y in pts:
            nx = (int(A[0, 0]) * x + int(A[0, 1]) * y + int(b[0])) % p
            ny = (int(A[1, 0]) * x + int(A[1, 1]) * y + int(b[1])) % p
            m[nx + ny * p] = 1.0
        hat = m @ EM.T
        off = ~Omega
        off[0] = False
        if float(np.max(np.abs(hat[off])) ** 2) > 1e-6:
            continue
        types = [tuple(sorted(int(x) for x in (m @ B.T))) for B in blocks]
        sig = tuple(sorted((tt, n) for tt, n in Counter(types).items()))
        if not is_triple_rainbow(sig, p):
            continue
        y = np.ones(q + 1, dtype=np.float64)
        y[1:] = 1.0 - 2.0 * m
        if float(np.max(np.abs(C @ y - p * y))) > 1e-6:
            continue
        out.append(L_on_mask(F, m))
        if len(out) >= want:
            break
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = live_cache_L_tr(p)
        named = L_tr_interp(p)
        drop = L_tr_drop12(p)
        if live != named:
            ok = False
        if p == 7 and drop == named:
            ok = False
        rows[str(p)] = {"live": str(live), "named": str(named), "drop12": str(drop)}
    # also 15.483 fold mean as a float check
    if abs(fold(5)["L_tr"] - 24.0) > 1e-8:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "L_interp=μ+²−χ+(1−χ)/12 hits live TR L++ at p=5,7. "
            "Fail: drop 1/12 (445/4≠1337/12)."
        ),
    }


def prove_B() -> dict:
    Ls = agl_tr_maxplus_L(11, max_tries=200, seed=1, want=4)
    pred = L_tr_interp(11)
    ok = len(Ls) >= 2
    if any(L == pred for L in Ls):
        ok = False
    if pred == Fraction(9089, 12) and Fraction(9089, 12) in Ls:
        ok = False
    return {
        "proved": bool(ok),
        "n": len(Ls),
        "values": [str(L) for L in Ls],
        "pred": str(pred),
        "theorem": (
            "p=11 AGL Max+ TR L++ misses L_interp=9089/12. "
            "Fail: claim the interpolant at p=11."
        ),
    }


def prove_C() -> dict:
    Ls = agl_tr_maxplus_L(11, max_tries=200, seed=1, want=4)
    nuniq = len(set(Ls))
    ok = nuniq >= 2
    return {
        "proved": bool(ok),
        "nuniq_p11": nuniq,
        "values": [str(L) for L in Ls],
        "theorem": (
            "TR L++ is not a function of p (≥2 values at p=11). "
            "Fail: one TR L++ at p=11."
        ),
    }


def prove_D() -> dict:
    F = _kernel_field(5)
    raw = L_on_mask(F, triangle_mask(5))
    ok = raw == 18 and raw != 24
    return {
        "proved": bool(ok),
        "raw5": str(raw),
        "tr5": "24",
        "theorem": "Raw T L++=18≠24 at p=5. Fail: raw T is the TR 4-point.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "TR L++ is not a formula in p. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.484  TR L++ is not a formula in p", flush=True)
    A = prove_A()
    print(f"  A interp: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B p11 miss: {B['proved']} pred={B['pred']} got={B['values']}", flush=True)
    C = prove_C()
    print(f"  C not a fn of p: {C['proved']} nuniq={C['nuniq_p11']}", flush=True)
    D = prove_D()
    print(f"  D raw T: {D['proved']} {D['raw5']}", flush=True)
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["proved"], "B": B["proved"], "C": C["proved"], "D": D["proved"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.484",
        "title": "TR L++ is not a formula in p; 2-point interpolant dies at p=11",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "interp_hits_live_57": A["proved"],
            "interp_dies_p11": B["proved"],
            "TR_L_not_fn_of_p": C["proved"],
            "raw_T_not_TR_4pt": D["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15484.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
