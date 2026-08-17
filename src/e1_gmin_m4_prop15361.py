#!/usr/bin/env python3
"""
Prop 15.361 — QR half-net count: n_AP is named; every live Max+ is
a half-net; n_1d + n_Hoff names |H_+| at p=5 and fails at p=3,7.
Geometric 3-feature catalog has 0 hits for N, n_NL, D.
D still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.360 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (binomial arithmetic).
  # of k-APs in F_p is p(p−1)/2.  Times r=(p+1)/2 QR directions:
      n_AP = p(p²−1)/4 = |G|/(4p)
  Equals n_1d at p=3,5 and is 84 < 140 at p=7.
  Fail: replace #APs by C(p,m) (then n_AP=n_1d at p=7).  ∎

Theorem B — PROVED (cache line-sums; p=5,7).
  Every y_∞=+1 Max+ is a QR half-net: at least r const-sum-1
  directions.  130/130 and 5726/5726.  Fail: claim every NL has
  p such directions.  ∎

Theorem C — PROVED as a negative (live cardinalities).
  n_1d + p²(p−1) equals 130 at p=5 and equals 24≠6 at p=3,
  434≠5726 at p=7.  The 1D+Hoffman mix is not |H_+|.
  Fail: claim 434=5726.  ∎

Theorem D — PROVED (ProcessPool 3-feature census).
  0 products f g / h, 0 products f g h, 0 sums f+g on the
  geometric catalog hit N∈{6,130,5726}, n_NL∈{0,100,5586},
  or D∈{1,13,409}.  Fail: claim N=n_1d+n_Hoff is such a hit
  (it misses p=7).  ∎

Theorem E — OPEN.  The half-net count N is not a Max+-free
  formula in p.  phi_F not imported.

============================================================================
Backend: Fraction algebra + cache line-sums.  Writes
evidence/e1_gmin_m4_prop15361.json
"""
from __future__ import annotations

import json
import sys
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15308 import AP_orbit, G_order  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15361_catalog.json"
N_LIVE = {3: 6, 5: 130, 7: 5726}


def n_AP_named(p: int) -> int:
    return p * (p * p - 1) // 4


def n_AP_from_qr_aps(p: int) -> int:
    r = (p + 1) // 2
    n_aps = p * (p - 1) // 2
    return r * n_aps


def n_hoff_named(p: int) -> int:
    return p * p * (p - 1)


def n_1d_plus_hoff(p: int) -> int:
    return n_1d(p) + n_hoff_named(p)


def directions(p: int):
    dirs = [[[c + p * x1 for x1 in range(p)] for c in range(p)]]
    for s in range(p):
        dirs.append(
            [[x0 + p * ((s * x0 + c) % p) for x0 in range(p)] for c in range(p)]
        )
    return dirs


def n_const1(z, p: int) -> int:
    n = 0
    for lines in directions(p):
        if all(int(z[ln].sum()) == 1 for ln in lines):
            n += 1
    return n


def halfnet_census(p: int, path: str) -> dict:
    Y = np.load(path)
    q = p * p
    Z = Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)
    r = (p + 1) // 2
    n_ok = 0
    mins = []
    for z in Z:
        c = n_const1(z, p)
        mins.append(c)
        if c >= r:
            n_ok += 1
    return {
        "n": int(Z.shape[0]),
        "n_halfnet": n_ok,
        "const1_min": int(min(mins)),
        "const1_max": int(max(mins)),
        "r": r,
    }


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (3, 5, 7, 11, 13, 17, 19):
        a = n_AP_named(p)
        b = n_AP_from_qr_aps(p)
        c = AP_orbit(p)
        rows[str(p)] = {"n_AP": a, "from_qr": b, "AP_orbit": c, "n_1d": n_1d(p)}
        if a != b or a != c:
            ok = False
    if n_AP_named(5) != n_1d(5):
        ok = False
    if n_AP_named(7) == n_1d(7):
        ok = False
    # fail-when-wrong: C(p,m) in place of #APs
    p = 7
    fake = ((p + 1) // 2) * comb(p, (p + 1) // 2)
    if n_AP_named(7) == fake:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "fake_Cpm_p7": fake,
        "theorem": "n_AP=p(p²−1)/4=|G|/(4p); not n_1d at p=7.",
    }


def prove_B() -> dict:
    recs = {}
    ok = True
    for p, path in ((5, "/tmp/maxplus_p5.npy"), (7, "/tmp/maxplus_p7.npy")):
        rec = halfnet_census(p, path)
        recs[str(p)] = rec
        if rec["n_halfnet"] != N_LIVE[p] or rec["n"] != N_LIVE[p]:
            ok = False
        if rec["const1_min"] < rec["r"]:
            ok = False
    if recs.get("5", {}).get("const1_min") == 5 and recs.get("7", {}).get("const1_min") == 7:
        ok = False
    return {
        "proved": bool(ok),
        "rows": recs,
        "theorem": "Every live y_∞=+1 Max+ is a QR half-net.",
    }


def prove_C() -> dict:
    rows = {str(p): n_1d_plus_hoff(p) for p in (3, 5, 7)}
    ok = rows["5"] == 130 and rows["3"] != 6 and rows["7"] != 5726
    if rows["7"] == 5726:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "live": N_LIVE,
        "theorem": "n_1d+n_Hoff hits p=5 and misses p=3,7.",
    }


def prove_D() -> dict:
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    nN = int(cat.get("n_hits_N", -1))
    nNL = int(cat.get("n_hits_NL", -1))
    nD = int(cat.get("n_hits_D", -1))
    ok = nN == 0 and nNL == 0 and nD == 0
    if n_1d_plus_hoff(7) == 5726:
        ok = False
    if nN != 0:
        ok = False
    return {
        "proved": bool(ok),
        "n_hits_N": nN,
        "n_hits_NL": nNL,
        "n_hits_D": nD,
        "theorem": "No 3-feature geometric product names N, n_NL, or D.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "D_named_in_p": False,
        "Q_NL_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": (
            "Half-nets certified.  n_AP named.  1D+Hoffman is not N.  "
            "D unnamed.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.361  half-net certified; n_AP named; 1D+Hoff dies", flush=True)
    A = prove_A()
    print(f"  A n_AP: {A['proved']} p7={A['rows']['7']}", flush=True)
    B = prove_B()
    print(f"  B halfnet: {B['proved']} { {k: v['n_halfnet'] for k,v in B['rows'].items()} }", flush=True)
    C = prove_C()
    print(f"  C 1d+Hoff: {C['proved']} {C['rows']}", flush=True)
    D = prove_D()
    print(f"  D catalog: {D['proved']} N={D['n_hits_N']} NL={D['n_hits_NL']} D={D['n_hits_D']}", flush=True)
    E = prove_open()
    print(f"  E open: D={E['D_named_in_p']} phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.361",
        "title": "QR half-net certified; n_AP named; 1D+Hoffman not N",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "n_AP_named_in_p": A["proved"],
            "halfnet_certified": B["proved"],
            "one_d_plus_hoff_not_N": C["proved"],
            "no_3feature_name": D["proved"],
            "D_named_in_p": False,
            "Q_NL_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15361.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
