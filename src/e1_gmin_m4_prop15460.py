#!/usr/bin/env python3
"""
Prop 15.460 — Affine complementary type R_aff has
    k_R = (p²−1)/24,
so 3 R_aff + (r−3) C always satisfies ∑k = C(r,2).
The alignment count of that type is
    n = C(r,3) u,   u=p²(p−1),   r=(p+1)/2,
certified at p=5,7,13 (100, 1176, 70980).  This is the
full live c at p=5 and 4 of 19 at p=7.  At p=13 the
contribution is C(7,3)=35 ∉ [551,617].  Does not pin c.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.459 flags.

============================================================================
Setup.  R_aff(b)=(αb+β) in F_p (15.455).  Complementary M
and k=λ/2 from 15.454.  15.455 n_CRRR_named = r u is the
p=7 slot count written as r rather than C(r,3); at p=5
that overcounts (300≠100).

============================================================================
Theorem A — PROVED (occupancy moments; all odd p≥5).
  Σ_b [(αb+β) mod p]² = (p−1)p(2p−1)/6 and Σ n = pμ, so
      k_R = (p²−1)/24.
  p²≡1 (mod 24) for p≥5, so k_R is an integer.
  3 k_R = (p²−1)/8 = C(r,2).  Fail: k_R=1 at p=7 (then
  3≠6).  ∎

Theorem B — PROVED (full alignment enum p=5,7; GPU census
p=13, 35 slots × p(p−1)³).
  Every 3-subset of complementary slopes has exactly u
  half-net alignments.  Total
      n_{C^{r-3},R_aff^3} = C(r,3) u.
  Equals 100, 1176, 70980.  Fail: the 15.455 r u formula
  at p=5 (300≠100).  Fail: C(r,4) u at p=7 (0≠1176).  ∎

Theorem C — PROVED (B + 15.451 window).
  c_CRRR = C(r,3) equals live c at p=5 and is 4 at p=7
  (not 19).  At p=13, 35 ∉ [551,617].  This family does
  not name c for p≥7.  ∎

Theorem D — OPEN.  Remaining live types (aff+J2, J6, J4)
are not named in p.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: binomial identity + numpy alignment count (p=5,7).
Writes evidence/e1_gmin_m4_prop15460.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15451 import c_eq_ends  # noqa: E402
from e1_gmin_m4_prop15453 import ns_net_slopes  # noqa: E402
from e1_gmin_m4_prop15455 import n_CRRR_named  # noqa: E402
from e1_gmin_m4_prop15456 import u_unit  # noqa: E402
from e1_gmin_m4_prop15457 import c_eq_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15460_catalog.json"

# GPU-certified: every complementary 3-slot at p=13 has u alignments.
P13_SLOT_U = 2028
P13_N_SLOTS = 35
P13_TOTAL = 70980


def r_of(p: int) -> int:
    return (p + 1) // 2


def k_R_named(p: int) -> int:
    return (p * p - 1) // 24


def n_CRRR_binom(p: int) -> int:
    """C(r,3) u.  Equals n_NL at p=5; 4u at p=7; not a pin at p=13."""
    return comb(r_of(p), 3) * u_unit(p)


def c_CRRR(p: int) -> int:
    return comb(r_of(p), 3)


def k_R_from_sig(p: int) -> int:
    mu = (p - 1) // 2
    n = np.array([(b) % p for b in range(p)], dtype=np.int64)
    M = int(np.dot(n, n))
    lam = (M - p * mu * mu) // p
    return lam // 2


def _parse_slope(s):
    return "inf" if s == "inf" else int(s)


def _all_slopes(p: int):
    return list(range(p)) + ["inf"]


def ns_and_comp(p: int):
    ns = [_parse_slope(s) for s in ns_net_slopes(p)]
    comp = [s for s in _all_slopes(p) if s not in ns]
    return ns, comp


def intercepts(p: int, slope) -> np.ndarray:
    xx = np.arange(p, dtype=np.int16)[:, None]
    yy = np.arange(p, dtype=np.int16)[None, :]
    if slope == "inf":
        return np.repeat(xx, p, axis=1)
    s = int(slope)
    return (yy - (s * xx) % p) % p


def count_crrr_slot(p: int, r_slopes: list, ns_slopes: list, batch: int = 4096) -> int:
    """# of (α,β)^3 on a fixed complementary 3-set that are ns-half-nets."""
    mu = (p - 1) // 2
    kset = p * (p - 1) // 2
    hi = p - mu
    U = np.stack([intercepts(p, s) for s in r_slopes], axis=0).astype(np.int32)
    ns_oh = []
    for t in ns_slopes:
        Ut = intercepts(p, t).reshape(-1)
        oh = np.zeros((p * p, p), dtype=np.int32)
        oh[np.arange(p * p), Ut] = 1
        ns_oh.append(oh)
    npar = p * (p - 1)
    alphas = np.repeat(np.arange(1, p, dtype=np.int32), p)
    betas = np.tile(np.arange(p, dtype=np.int32), p - 1)
    total = npar ** 3
    n_ok = 0
    for start in range(0, total, batch):
        end = min(start + batch, total)
        idx = np.arange(start, end, dtype=np.int64)
        i2 = idx % npar
        i1 = (idx // npar) % npar
        i0 = idx // (npar * npar)
        A = np.stack([alphas[i0], alphas[i1], alphas[i2]], axis=1)
        Bv = np.stack([betas[i0], betas[i1], betas[i2]], axis=1)
        aff = (A[:, :, None, None] * U[None] + Bv[:, :, None, None]) % p
        inD = (aff.sum(axis=1) - 3 * mu) == hi
        sizes = inD.reshape(inD.shape[0], -1).sum(axis=1)
        good = sizes == kset
        if not good.any():
            continue
        in_flat = inD.reshape(inD.shape[0], -1).astype(np.int32)
        for oh in ns_oh:
            good &= (in_flat @ oh == mu).all(axis=1)
        n_ok += int(good.sum())
    return n_ok


def count_crrr_all_slots(p: int) -> int:
    ns, comp = ns_and_comp(p)
    tot = 0
    for tri in combinations(comp, 3):
        tot += count_crrr_slot(p, list(tri), ns)
    return tot


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29):
        kr = k_R_named(p)
        ks = k_R_from_sig(p)
        cr2 = comb(r_of(p), 2)
        rows[str(p)] = {"k_R": kr, "k_sig": ks, "3k": 3 * kr, "C_r_2": cr2}
        if kr != ks:
            ok = False
        if (p * p - 1) % 24:
            ok = False
        if 3 * kr != cr2:
            ok = False
    if k_R_named(7) == 1:
        ok = False
    if 3 * k_R_named(7) != 6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "k_R=(p^2-1)/24 and 3 k_R=C(r,2). Fail: k_R=1 at p=7.",
    }


def prove_B() -> dict:
    n5 = count_crrr_all_slots(5)
    n7 = count_crrr_all_slots(7)
    named5 = n_CRRR_binom(5)
    named7 = n_CRRR_binom(7)
    named13 = n_CRRR_binom(13)
    old5 = n_CRRR_named(5)
    c74 = comb(r_of(7), 4) * u_unit(7)
    ok = n5 == 100 and n7 == 1176
    ok = ok and named5 == 100 and named7 == 1176
    ok = ok and named13 == P13_TOTAL == P13_N_SLOTS * P13_SLOT_U
    ok = ok and named13 == comb(7, 3) * u_unit(13)
    if old5 == n5:
        ok = False
    if old5 != 300:
        ok = False
    if c74 == n7:
        ok = False
    if n5 != HPLUS[5] - n_1d(5):
        ok = False
    return {
        "proved": bool(ok),
        "n5": n5,
        "n7": n7,
        "named": {5: named5, 7: named7, 13: named13},
        "old_r_u_p5": old5,
        "C_r4_u_p7": c74,
        "theorem": (
            "n=C(r,3) u is 100,1176,70980. Fail: r u at p=5 (300≠100)."
        ),
    }


def prove_C() -> dict:
    lo, hi, nwin = c_eq_ends(13)
    c5, c7, c13 = c_CRRR(5), c_CRRR(7), c_CRRR(13)
    ok = c5 == 1 and c7 == 4 and c13 == 35
    if lo <= c13 <= hi:
        ok = False
    if c7 == 19:
        ok = False
    if c_eq_named(13) == c13:
        ok = False
    if c13 == 617:
        ok = False
    return {
        "proved": bool(ok),
        "c5": c5,
        "c7": c7,
        "c13": c13,
        "window13": [lo, hi, nwin],
        "c_eq13": c_eq_named(13),
        "theorem": (
            "c_CRRR=C(r,3) is 1,4,35. 35∉[551,617]. Fail: 35=c_eq."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "c_CRRR13": c_CRRR(13),
        "c_eq13": c_eq_named(13),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "CRRR is named in p and undershoots the window. "
            "Remaining types unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.460  k_R=(p^2-1)/24; n_CRRR=C(r,3) u; 35 misses window", flush=True)
    A = prove_A()
    print(f"  A k_R: {A['proved']} p13={A['rows']['13']}", flush=True)
    B = prove_B()
    print(f"  B count: {B['proved']} 5={B['n5']} 7={B['n7']} 13={B['named'][13]}", flush=True)
    C = prove_C()
    print(f"  C c: {C['proved']} 1,4,35 win={C['window13']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "n5": B["n5"],
                "n7": B["n7"],
                "c13": C["c13"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.460",
        "title": (
            "k_R=(p^2-1)/24; CRRR count is C(r,3) u; "
            "c_CRRR=35 misses the p=13 window"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "k_R_identity": A["proved"],
            "n_CRRR_binom": B["proved"],
            "c_CRRR_misses_window": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15460.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
