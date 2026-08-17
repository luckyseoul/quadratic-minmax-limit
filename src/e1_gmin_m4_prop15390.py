#!/usr/bin/env python3
"""
Prop 15.390 — R-secant meet m(j) is the F_p-character of the
circle-line discriminant
    Δ = Tr(w s^p)² − 4 N(s)(N(w)−c),
w = u_j − t.  k=|T∩S∩ℓ| is not a function of (n,m).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.389 flags.

============================================================================
Setup.  15.389: n(j)=aj+b on R-dirs (p≡1).  first_nc(T) circle
S={N(u−t)=c}, |S|=p+1.  Line ℓ_j = u_j + F_p·s.

============================================================================
Theorem A — PROVED (circle-line quadratic; certified p=5,13,17).
  N(w+λs)=N(w)+λ Tr(w s^p)+λ² N(s).  The number of roots λ∈F_p
  is
      m = 1_{Δ=0} + 2·1_{χ_p(Δ)=+1}.
  Matches live |S∩ℓ| on every R-line (15+39+51).  Fail: flip χ
  (then 6 mismatches at each of p=5,13,17).  ∎

Theorem B — PROVED (negative; (n,m)↛k).
  The same (n,m) takes two k-values (10 split pairs across
  p=5,13,17), e.g. (3,2) has k∈{0,2} and (2,1) has k∈{0,1}.
  Fail: claim k is a function of (n,m).  ∎

Theorem C — OPEN.  A Max+-free rule for k (or a formula for c)
  still unnamed.  t=0 on first_nc(T) at p=5,13,17 is not a c-namer.
  phi_F not imported.

============================================================================
Backend: first_nc(T) + field disc.  Writes
evidence/e1_gmin_m4_prop15390.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15312 import _norm_table, first_nc  # noqa: E402
from e1_gmin_m4_prop15351 import y_T  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15390_catalog.json"


def chi_p(x: int, p: int) -> int:
    x %= p
    if x == 0:
        return 0
    t = pow(x, (p - 1) // 2, p)
    return 1 if t == 1 else -1


def frob_p(F: dict, z: int) -> int:
    if z == 0:
        return 0
    q, p = F["q"], F["p"]
    units = [0] * (q - 1)
    g, cur = F["g"], 1
    for j in range(q - 1):
        units[j] = cur
        cur = F["mul"][cur][g]
    return units[(F["dlog"][z] * p) % (q - 1)]


def disc_line(F: dict, Nt, u0: int, s: int, t: int, c: int) -> int:
    """Δ = Tr(w s^p)² − 4 N(s)(N(w)−c) in F_p."""
    p = F["p"]
    w = F["add"][u0][F["neg"][t]]
    Nw = int(Nt[w])
    Ns = int(Nt[s])
    sp = frob_p(F, s)
    B = int(F["trt"][F["mul"][w][sp]])
    return (B * B - 4 * Ns * ((Nw - c) % p)) % p


def m_from_disc(disc: int, p: int) -> int:
    if disc == 0:
        return 1
    return 2 if chi_p(disc, p) == 1 else 0


def m_from_disc_flip(disc: int, p: int) -> int:
    if disc == 0:
        return 1
    return 2 if chi_p(disc, p) == -1 else 0


def r_dir_rows(p: int):
    F = _kernel_field(p)
    Nt = _norm_table(F)
    y0 = y_T(p)
    pack = first_nc(p, y0)
    D = y0[1 : 1 + p * p] < 0
    S = {int(u) for u in pack["circle"]}
    t, c = int(pack["t"]), int(pack["c"])
    rows = []
    for lines in square_classes(F):
        ns = D[lines].sum(axis=1)
        if tag_sig(p, tuple(int(x) for x in sorted(ns))) != "R":
            continue
        s = int(F["add"][int(lines[0][1])][F["neg"][int(lines[0][0])]])
        for j, ln in enumerate(lines):
            ln = [int(x) for x in ln]
            u0 = ln[0]
            disc = disc_line(F, Nt, u0, s, t, c)
            m_live = sum(1 for x in ln if x in S)
            k_live = sum(1 for x in ln if x in S and D[x])
            rows.append(
                {
                    "n": int(ns[j]),
                    "m_live": m_live,
                    "k_live": k_live,
                    "disc": disc,
                    "m_pred": m_from_disc(disc, p),
                    "m_flip": m_from_disc_flip(disc, p),
                    "t": t,
                    "c": c,
                }
            )
    return pack, rows


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17):
        pack, recs = r_dir_rows(p)
        n_ok = sum(1 for r in recs if r["m_pred"] == r["m_live"])
        n_flip = sum(1 for r in recs if r["m_flip"] == r["m_live"])
        ok = ok and n_ok == len(recs) and len(recs) > 0
        ok = ok and n_flip < len(recs)
        if n_flip == len(recs):
            ok = False
        rows[str(p)] = {
            "t": int(pack["t"]),
            "c": int(pack["c"]),
            "n": len(recs),
            "n_ok": n_ok,
            "n_flip": n_flip,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "m=1_{Δ=0}+2·1_{χ_p(Δ)=+1}. Fail: flip χ."
        ),
    }


def prove_B() -> dict:
    ok = True
    km = defaultdict(set)
    for p in (5, 13, 17):
        _pack, recs = r_dir_rows(p)
        for r in recs:
            km[(r["n"], r["m_live"])].add(r["k_live"])
    splits = {f"{n},{m}": sorted(v) for (n, m), v in km.items() if len(v) > 1}
    ok = ok and len(splits) >= 1
    ok = ok and (3, 2) in km and km[(3, 2)] != {0} and km[(3, 2)] != {2}
    if len(splits) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "n_split": len(splits),
        "splits": splits,
        "theorem": "k is not a function of (n,m). Fail: unique k per (n,m).",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "m named from disc; k and c unnamed. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.390  m from circle-line disc; k not f(n,m)", flush=True)
    A = prove_A()
    print(f"  A disc: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B k-split: {B['proved']} n_split={B['n_split']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["rows"], "splits": B["splits"]}, indent=2) + "\n")
    out = {
        "prop": "15.390",
        "title": "m(j) from circle-line disc; k not a function of (n,m)",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "m_from_disc": A["proved"],
            "k_not_fn_nm": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15390.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
