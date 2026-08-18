#!/usr/bin/env python3
"""
Prop 15.517 — p=5 H_+ W-histogram is the 15.323 two-orbit law.

  W=∑_{u∈U} N(u) on H_+.  At p=5 the values are
      {EW−2p, EW−(p+1), EW+(p+1), EW+2p} = {20,24,36,40}
  with masses 15:50:50:15.  Hence Var(W)=Var_2orb=660/13.
  Fail: claim the same 4-point support at p=7 (13 values);
  claim Var=Var_2orb at p=7 (21504/409 < 280/3).
  Does not prove Var≤Var_2orb for p≡1.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.516 flags.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15290 import field_norm
from e1_gmin_m4_prop15298 import load_N
from e1_gmin_m4_prop15323 import EW_named, Var_2orb

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15517.json"

W = 86
CERT = (5, 7)


def support_named(p: int) -> tuple[int, int, int, int]:
    ew = EW_named(p)
    return (ew - 2 * p, ew - (p + 1), ew + (p + 1), ew + 2 * p)


def _shard(spec):
    p, lo, hi = spec
    F = _kernel_field(p)
    N = load_N(p, F)
    U = [d for d in range(1, F["q"]) if int(field_norm(F, d)) == 1]
    return [int(sum(int(N[i, u]) for u in U)) for i in range(lo, hi)]


def fold(p: int) -> dict:
    F = _kernel_field(p)
    H = int(load_N(p, F).shape[0])
    sw = max(1, (H + W - 1) // W)
    jobs = [(p, lo, min(H, lo + sw)) for lo in range(0, H, sw)]
    ws = []
    with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
        for part in ex.map(_shard, jobs):
            ws.extend(part)
    mean = Fraction(sum(ws), len(ws))
    var = Fraction(sum((w - mean) ** 2 for w in ws), len(ws))
    return {
        "n": len(ws),
        "mean": str(mean),
        "var": str(var),
        "hist": {str(k): v for k, v in sorted(Counter(ws).items())},
        "nuniq": len(set(ws)),
        "support": sorted(set(ws)),
    }


def prove_A() -> dict:
    rec = fold(5)
    want = {20: 15, 24: 50, 36: 50, 40: 15}
    live = {int(k): v for k, v in rec["hist"].items()}
    ok = live == want
    if rec["var"] != str(Var_2orb(5)):
        ok = False
    if rec["mean"] != str(EW_named(5)):
        ok = False
    if tuple(rec["support"]) != support_named(5):
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "want": want,
        "Var_2orb": str(Var_2orb(5)),
        "theorem": "p=5 H_+ W-hist is 2-orbit {20,24,36,40}=15:50:50:15.",
    }


def prove_B() -> dict:
    rec = fold(7)
    ok = rec["nuniq"] != 4
    if rec["var"] == str(Var_2orb(7)):
        ok = False
    if rec["var"] != "21504/409":
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "Var_2orb": str(Var_2orb(7)),
        "theorem": "p=7 has 13 W-values; Var<Var_2orb. Fail 4-point / Var=Var_2orb.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Two-orbit W-law is cache at p=5, false as support at p=7. "
            "Var≤Var_2orb unproved for p≡1. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.517 p=5 H_+ W is 2-orbit; p=7 is not", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.517",
        "title": "p=5 H_+ W-hist is 2-orbit; Var=Var_2orb; p=7 fails",
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
    print(
        f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
