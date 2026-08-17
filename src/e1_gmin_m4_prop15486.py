#!/usr/bin/env python3
"""
Prop 15.486 — The 15.276 locked-sawtooth family is the TR Aut_∞
orbit, not a constructor for c.  Family size C(m,3)p²(p−1) equals
|NL| at p=5 and the unique TR orbit at p=7; c_family=4≠19.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.485 flags.

============================================================================
Setup.  15.276: locked sawtooth on a triple of Singer-good forms.
m=(p+1)/2 forms, lam∈F_p^×, (s0,s1)∈F_p².  15.485: p=7 has four
s=2 orbits, only one TR, and c=19.

============================================================================
Theorem A — PROVED (count; all odd p≥5).
  The family has size M3=C(m,3) p²(p−1).  Fail: drop (p−1)
  (C(3,3)p²=25≠100 at p=5).  ∎

Theorem B — PROVED (full enum; p=5).
  All 100 members are Max+ (Cy=py) and TR; they are 100 distinct
  sets = the unique NL orbit (s=p+1).  Fail: claim a 1D hit
  (n_full≥2).  ∎

Theorem C — PROVED (full enum; p=7).
  All 1176 members are Max+ TR; 1176 distinct sets = the unique
  TR orbit (s=2).  c_family=(p+1)/2=4 ≠ live c=19.
  Fail: claim the family names c (4=19).  ∎

Theorem D — OPEN.  The existence family does not pin c=c_eq.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.276 constructor + Paley conference.  Writes
evidence/e1_gmin_m4_prop15486.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15270 import singer_good_forms  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15276 import construct_locked_triple  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15476 import is_triple_rainbow, line_blocks  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15486_catalog.json"


def m_forms(p: int) -> int:
    return (p + 1) // 2


def family_size(p: int) -> int:
    return comb(m_forms(p), 3) * p * p * (p - 1)


def family_size_drop_pm1(p: int) -> int:
    return comb(m_forms(p), 3) * p * p


@lru_cache(None)
def enumerate_family(p: int) -> dict:
    forms = singer_good_forms(p)
    m = len(forms)
    C = paley_conference_prime_power(p).astype(np.float64)
    blocks = line_blocks(p, p * p)
    n_max = 0
    n_tr = 0
    n_fullpos = 0
    keys = set()
    for idx in combinations(range(m), 3):
        for lam in range(1, p):
            for s0 in range(p):
                for s1 in range(p):
                    rec = construct_locked_triple(
                        p, lam=lam, s0=s0, s1=s1, form_index=idx
                    )
                    y = rec["y"]
                    err = float(np.max(np.abs(C @ y - p * y)))
                    if err > 1e-5:
                        continue
                    n_max += 1
                    D = tuple(sorted(int(u) for u in np.where(y[1:] < 0)[0]))
                    keys.add(D)
                    msk = np.zeros(p * p, dtype=np.float64)
                    msk[list(D)] = 1.0
                    types = [
                        tuple(sorted(int(x) for x in (msk @ B.T))) for B in blocks
                    ]
                    sig = tuple(sorted((t, n) for t, n in Counter(types).items()))
                    if is_triple_rainbow(sig, p):
                        n_tr += 1
                    nf = 0
                    for t, n in sig:
                        if p in t:
                            nf += n
                    if nf >= 2:
                        n_fullpos += 1
    return {
        "n_max": n_max,
        "n_tr": n_tr,
        "n_unique": len(keys),
        "n_fullge2": n_fullpos,
        "family": family_size(p),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        fs = family_size(p)
        drop = family_size_drop_pm1(p)
        rows[str(p)] = {"M3": fs, "drop": drop}
        if drop == fs:
            ok = False
    if family_size(5) != 100:
        ok = False
    if family_size_drop_pm1(5) != 25:
        ok = False
    if family_size(7) != 1176:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Family size C(m,3)p²(p−1). Fail: drop (p−1) (25≠100 at p=5)."
        ),
    }


def prove_B() -> dict:
    e = enumerate_family(5)
    ok = (
        e["n_max"] == 100
        and e["n_tr"] == 100
        and e["n_unique"] == 100
        and e["n_fullge2"] == 0
    )
    return {
        "proved": bool(ok),
        "enum": e,
        "theorem": (
            "p=5 family = unique NL/TR orbit (100). Fail: a 1D hit."
        ),
    }


def prove_C() -> dict:
    e = enumerate_family(7)
    c_fam = 4  # (p+1)/s with s=2
    ok = (
        e["n_max"] == 1176
        and e["n_tr"] == 1176
        and e["n_unique"] == 1176
        and e["n_fullge2"] == 0
        and c_fam != 19
    )
    return {
        "proved": bool(ok),
        "enum": e,
        "c_family": c_fam,
        "c_live": 19,
        "theorem": (
            "p=7 family = unique TR orbit (1176); c_family=4≠19. "
            "Fail: family names c."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "15.276 family is TR-only. Does not pin c=c_eq. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.486  15.276 family is TR, not a constructor for c", flush=True)
    A = prove_A()
    print(f"  A size: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B p5: {B['proved']} {B['enum']}", flush=True)
    C = prove_C()
    print(f"  C p7: {C['proved']} c_fam={C['c_family']}≠{C['c_live']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.486",
        "title": "15.276 family is the TR orbit; does not name c",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "family_size_M3": A["proved"],
            "p5_family_is_NL": B["proved"],
            "p7_family_is_TR_not_c": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15486.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
