#!/usr/bin/env python3
"""
Prop 15.511 — 1D summand of involution Fix.

  The number of Type+ 1D Max+ fixed by any g_b: u↦−u+b equals
      n_1D,inv = ((p+1)/2) · C((p−1)/2, ⌊(p−1)/4⌋)
  = (# Type+ forms) × (#{A⊂F_p : |A|=m, A=−A}), m=(p+1)/2.
  Live 1D flags: 2,6,12 at p=3,5,7.
  Fail: n_1d/p (2,6,20≠12 at p=7); drop the (p+1)/2 factor (1,2,3).
  Full Fix is 2,6,54 — NL leftover 0,0,42 unnamed.
  Does not name D or Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.510 flags.

============================================================================
Setup.  Type+ forms: n_type_plus_forms=(p+1)/2 (15.292).  A Type+
lift has y_x=σ(L(x)) with |plus|=m and ker L a square line.
g_0: u↦−u sends L(x)↦−L(x), so the lift is g_0-invariant iff
plus=−plus.  Conjugates g_b have the same 1D Fix (class-constant).

============================================================================
Theorem A — PROVED (field; odd p).
  #{A⊂F_p : |A|=m, A=−A} = C((p−1)/2, ⌊(p−1)/4⌋).
  Fail: C((p−1)/2, 1) (hits 1,2,3; dies at p=11: 5≠10).  ∎

Theorem B — PROVED (field; 15.292).
  n_1D,inv := n_type_plus_forms(p) · (A) = ((p+1)/2)·C((p−1)/2,⌊(p−1)/4⌋).
  Fail: n_1d/p; drop the form-count.  ∎

Theorem C — PROVED (cache; certified p=3,5,7).
  On H_+, the number of 1D designs fixed by g_0 (hence any g_b)
  equals n_1D,inv (2,6,12).  Fail: claim this equals full Fix
  (12≠54 at p=7).  ∎

Theorem D — OPEN.  NL leftover of Fix is 0,0,42.  Full Fix
  formula 15.509 C is still a 3-point cache identity.
  D / Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: field subset count + cache 1D flags (ProcessPool).
Writes evidence/e1_gmin_m4_prop15511.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import n_type_plus_forms
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15307 import load_D
from e1_gmin_m4_prop15308 import pack_y
from e1_gmin_m4_prop15360 import is_1d_z
from e1_gmin_m4_prop15508 import apply_map
from e1_gmin_m4_prop15509 import fix_named as fix_full_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15511.json"

CERT = (3, 5, 7)
FIELD_PS = (3, 5, 7, 11, 13, 17, 19)
W = 86


def m_half(p: int) -> int:
    return (p + 1) // 2


def n_even_m_subsets_named(p: int) -> int:
    return comb((p - 1) // 2, (p - 1) // 4)


def n_even_m_subsets_wrong(p: int) -> int:
    """Fail-when-wrong: C((p-1)/2, 1)=(p-1)/2 (hits p=3,5,7; dies at p=11)."""
    return (p - 1) // 2


def n_1d_inv_named(p: int) -> int:
    return n_type_plus_forms(p) * n_even_m_subsets_named(p)


def n_1d_inv_wrong_n1dp(p: int) -> int:
    return n_1d(p) // p


def n_1d_inv_wrong_drop_m(p: int) -> int:
    return n_even_m_subsets_named(p)


def n_even_m_subsets_enum(p: int) -> int:
    m = m_half(p)
    n = 0
    for A in combinations(range(p), m):
        s = set(A)
        if s == {(-x) % p for x in s}:
            n += 1
    return n


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in FIELD_PS:
        live = n_even_m_subsets_enum(p)
        named = n_even_m_subsets_named(p)
        wrong = n_even_m_subsets_wrong(p)
        if live != named:
            ok = False
        # C((p-1)/2,1) interpolates p=3,5,7 and dies at p=11 (5≠10).
        if p == 11 and named == wrong:
            ok = False
        rows[str(p)] = {"live": live, "named": named, "wrong": wrong}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "#{A:|A|=m, A=−A}=C((p−1)/2,⌊(p−1)/4⌋). Fail C(*,1).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in FIELD_PS:
        named = n_1d_inv_named(p)
        w1 = n_1d_inv_wrong_n1dp(p)
        w2 = n_1d_inv_wrong_drop_m(p)
        if named != n_type_plus_forms(p) * n_even_m_subsets_named(p):
            ok = False
        if p == 7 and named == w1:
            ok = False
        if named == w2:
            ok = False
        rows[str(p)] = {"named": named, "wrong_n1dp": w1, "wrong_drop_m": w2}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "n_1D,inv=((p+1)/2)·C((p−1)/2,⌊(p−1)/4⌋). Fail n_1d/p.",
    }


def _one_d_job(spec):
    p, lo, hi = spec
    D, F = load_D(p)
    Ys = np.stack([pack_y(D[i]) for i in range(D.shape[0])])
    g = apply_map(F, int(F["neg1"]), 0)
    n1 = 0
    nfix = 0
    for i in range(lo, hi):
        if not np.all(Ys[i, g] == Ys[i]):
            continue
        nfix += 1
        z = Ys[i, 1 : 1 + F["q"]]
        if is_1d_z(z, p):
            n1 += 1
    return n1, nfix


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        D, _ = load_D(p)
        H = int(D.shape[0])
        sw = max(1, (H + W - 1) // W)
        jobs = [(p, lo, min(H, lo + sw)) for lo in range(0, H, sw)]
        n1 = nfix = 0
        with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
            for a, b in ex.map(_one_d_job, jobs):
                n1 += a
                nfix += b
        named = n_1d_inv_named(p)
        full = fix_full_named(p)
        if n1 != named:
            ok = False
        if p == 7 and named == full:
            ok = False
        rows[str(p)] = {
            "n_1d_live": n1,
            "n_fix_live": nfix,
            "named": named,
            "full_fix": full,
            "nl_left": nfix - n1,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "1D Fix(g_0)=n_1D,inv. Fail claim = full Fix.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "nl_leftover": {str(p): fix_full_named(p) - n_1d_inv_named(p) for p in CERT},
        "note": (
            "1D involution Fix is named. NL leftover 0,0,42 unnamed. "
            "D / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.511 1D involution Fix = ((p+1)/2)·C((p-1)/2,⌊(p-1)/4⌋)", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "prop": "15.511",
        "title": "1D involution Fix named; NL leftover / D unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
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
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
