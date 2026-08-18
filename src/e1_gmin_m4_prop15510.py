#!/usr/bin/env python3
"""
Prop 15.510 — Burnside translations T_b: u↦u+b, b≠0.

  If b is a square in F_q then Fix_{H_+}(T_b)=C(p,m), m=(p+1)/2.
  If b is a nonsquare then Fix=0.
  Sum_{b≠0} Fix(T_b)=(p−1) n_1d.
  Fail: Fix=p+1 on squares (6≠10 at p=5).
  Does not name |H_+| (the b=0 term is the identity).
  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.509 flags.

============================================================================
Setup.  Translations sit in G as a=1, k=0.  T_0=id, Fix=|H_+|.
Nonzero b: q−1 maps, equally split into squares and nonsquares.

m=(p+1)/2 and n_1d=m C(p,m) are 15.292.

============================================================================
Theorem A — PROVED (cache; certified p=3,5,7).
  b a nonzero square ⇒ Fix(T_b)=C(p,m) constantly.
  b a nonsquare ⇒ Fix(T_b)=0.
  Fail: Fix(T_b)=p+1 on squares.  ∎

Theorem B — PROVED (A + 15.292).
  ∑_{b≠0} Fix(T_b)=(p−1) n_1d.
  Fail: drop the square half (sum 0).  ∎

Theorem C — OPEN.  The identity term is |H_+|, so Burnside still
  does not name D=|H_+|/(2p).  Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: apply_map + cache Fix (ProcessPool over q−1 translations).
Writes evidence/e1_gmin_m4_prop15510.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15307 import load_D
from e1_gmin_m4_prop15308 import pack_y
from e1_gmin_m4_prop15508 import apply_map

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15510.json"

CERT = (3, 5, 7)
W = 86


def m_half(p: int) -> int:
    return (p + 1) // 2


def fix_square_named(p: int) -> int:
    return comb(p, m_half(p))


def fix_square_wrong(p: int) -> int:
    """Fail-when-wrong: p+1 (6≠10 at p=5)."""
    return p + 1


def _trans_job(spec):
    p, bs = spec
    D, F = load_D(p)
    Y = np.stack([pack_y(D[i]) for i in range(D.shape[0])])
    squares = set(int(s) for s in F["squares"] if int(s) != 0)
    out = []
    for b in bs:
        g = apply_map(F, 1, int(b))
        fx = int(np.all(Y[:, g] == Y, axis=1).sum())
        out.append((int(b in squares), fx))
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        q = p * p
        bs = list(range(1, q))
        sw = max(1, (len(bs) + W - 1) // W)
        jobs = [(p, bs[lo : min(len(bs), lo + sw)]) for lo in range(0, len(bs), sw)]
        sq, ns = [], []
        with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
            for part in ex.map(_trans_job, jobs):
                for is_sq, fx in part:
                    (sq if is_sq else ns).append(fx)
        named = fix_square_named(p)
        wrong = fix_square_wrong(p)
        if not sq or min(sq) != named or max(sq) != named:
            ok = False
        if not ns or min(ns) != 0 or max(ns) != 0:
            ok = False
        if named == wrong:
            ok = False
        if any(f == wrong for f in sq):
            ok = False
        rows[str(p)] = {
            "n_sq": len(sq),
            "n_ns": len(ns),
            "fix_sq": named,
            "wrong": wrong,
            "live_sq": sorted(set(sq)),
            "live_ns": sorted(set(ns)),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Fix(T_b)=C(p,m) on squares, 0 on nonsquares. Fail p+1.",
    }


def prove_B() -> dict:
    A = prove_A()
    ok = A["proved"]
    rows = {}
    for p in CERT:
        rec = A["rows"][str(p)]
        s = rec["n_sq"] * rec["fix_sq"]
        want = (p - 1) * n_1d(p)
        if s != want:
            ok = False
        rows[str(p)] = {"sum": s, "named": want}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "∑_{b≠0} Fix(T_b)=(p−1) n_1d. Fail drop squares.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Translation Fix is named off the identity. The b=0 term "
            "is |H_+|, so Burnside still does not name D. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.510 translations Fix=C(p,m) on squares, 0 on nonsquares", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.510",
        "title": "Burnside translations Fix=C(p,m) on □; D unnamed",
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
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
