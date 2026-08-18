#!/usr/bin/env python3
"""
Prop 15.515 — p=5 incidence of 26−7 Q_y^{++}.

  ∑_y (26 − 7 Q_y^{++}) = n_1d − k = 4p = 20.
  1D rows (n=30) have 26−7 Q=122/9>0; NL rows (n=100) have
  −58/15<0.  Ensemble 26−7 Q++=2/13>0 is the average, not a
  per-y count.  Fail: claim every summand is ≥0 (NL is negative);
  claim the sum equals n_1d (30≠20).
  Does not name Q_τ for p>5.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.514 flags.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15290 import omega_set, zhat_omega
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15300 import merge_cls
from e1_gmin_m4_prop15307 import load_D
from e1_gmin_m4_prop15360 import is_1d_z
from e1_gmin_m4_prop15367 import k_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15515.json"

P = 5
W = 86


def sum_named(p: int) -> int:
    return n_1d(p) - k_named(p)


def sum_named_4p(p: int) -> int:
    return 4 * p


def sum_wrong_n1d(p: int) -> int:
    """Fail-when-wrong: the sum equals n_1d."""
    return n_1d(p)


def _shard(spec):
    lo, hi = spec
    D, F = load_D(P)
    Omega = omega_set(F)
    om = {int(xi): i for i, xi in enumerate(Omega)}
    q2 = float(F["q"] ** 2)
    rpp = [int(r) for r in F["squares"] if merge_cls(F, int(r)) == "++"]
    out = []
    for i in range(lo, hi):
        z = np.where(D[i], -1.0, 1.0)
        u = np.abs(zhat_omega(z, F, Omega)) ** 2
        sm = 0.0
        for r in rpp:
            for xi in Omega:
                sm += u[om[int(xi)]] * u[om[F["mul"][r][int(xi)]]]
        qpp = sm / (len(Omega) * len(rpp) * q2)
        sl = Fraction(26 - 7 * qpp).limit_denominator(P**4)
        out.append((bool(is_1d_z(z, P)), sl))
    return out


def fold() -> dict:
    D, _ = load_D(P)
    H = int(D.shape[0])
    sw = max(1, (H + W - 1) // W)
    jobs = [(lo, min(H, lo + sw)) for lo in range(0, H, sw)]
    rows = []
    with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
        for part in ex.map(_shard, jobs):
            rows.extend(part)
    tot = sum((sl for _, sl in rows), Fraction(0))
    n1 = sum(1 for one, _ in rows if one)
    nneg = sum(1 for _, sl in rows if sl < 0)
    by1 = Counter(str(sl) for one, sl in rows if one)
    byn = Counter(str(sl) for one, sl in rows if not one)
    return {
        "n": len(rows),
        "sum": str(tot),
        "n_1d": n1,
        "n_neg": nneg,
        "by_1d": dict(by1),
        "by_NL": dict(byn),
    }


def prove_A() -> dict:
    rec = fold()
    tot = Fraction(rec["sum"])
    ok = tot == sum_named(P) == sum_named_4p(P)
    if tot == sum_wrong_n1d(P):
        ok = False
    if rec["n_neg"] == 0:
        ok = False
    if rec["n_1d"] != n_1d(P):
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "named": sum_named(P),
        "wrong_n1d": sum_wrong_n1d(P),
        "theorem": "∑(26−7 Q_y)=n_1d−k=4p. Fail all≥0; sum=n_1d.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "p=5 sum is named and positive, but NL summands are "
            "negative, so this is not a nonnegative count. Not a "
            "p-law for Q++≤26/7. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.515 p=5 ∑(26-7 Q_y)=n_1d-k; NL negative", flush=True)
    A, D = prove_A(), prove_open()
    out = {
        "prop": "15.515",
        "title": "p=5 26−7 Q++ incidence; not a nonnegative count",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": D,
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
    print(f"A={A['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
