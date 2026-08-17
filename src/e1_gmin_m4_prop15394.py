#!/usr/bin/env python3
"""
Prop 15.394 — Constructors that survive p=37: T is still
3 R + (n_sq−3) C, Q_T Fejer still equals 72832/5985, and
QR0⊙NC is nonempty (first (t,c)=(p, p−2)).  AP⊙NC and
T⊙NC are empty.  Ensemble Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.393 flags.

============================================================================
Setup.  15.387 T occupancy; 15.350 Q_T Fejer; 15.309 AP⊙NC
empty at p=13; 15.393 T⊙NC empty at p=37.

============================================================================
Theorem A — PROVED (signatures; certified p=5,13,17,29,37,41).
  y_T has exactly 3 R square-dirs and n_sq−3 C-dirs at p=37
  and p=41 (same law as 15.387).  Fail: all-R at p=37
  (then n_C=0, live n_C=16).  ∎

Theorem B — PROVED (15.350 formula at p=37).
  Qpp_from_L(37)=72832/5985.  Drop-16 is 4552/5985.
  Fail: drop the 16.  ∎

Theorem C — PROVED (circle scan; certified p=37).
  AP-hs has 0 C_T-clique norm-circles.  QR0-hs has first
  clique at (t,c)=(37,35)=(p,p−2) and that circle is Hoffman.
  T⊙NC remains empty.  Fail: claim AP nonempty, or first
  QR0 c=ib=2.  ∎

Theorem D — OPEN.  QR0⊙NC exists at p=37 but is not the
  ensemble.  D / Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: sig_y + Qpp_from_L + one Paley C scan.  Writes
evidence/e1_gmin_m4_prop15394.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15312 import _norm_table  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_drop16, Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15351 import sig_y, y_T  # noqa: E402
from e1_gmin_m4_prop15387 import n_sq_named  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15394_catalog.json"
QT37 = Fraction(72832, 5985)

_C_CACHE: dict[int, np.ndarray] = {}


def _C(p: int) -> np.ndarray:
    if p not in _C_CACHE:
        _C_CACHE[p] = paley_conference_prime_power(p).astype(np.float64)
    return _C_CACHE[p]


def least_nsq(p: int) -> int:
    for a in range(1, p):
        if pow(a, (p - 1) // 2, p) == p - 1:
            return a
    raise ValueError(p)


def n_cliques_of(p: int, y0: np.ndarray, stop_first: bool = False) -> dict:
    y0 = np.asarray(y0, dtype=np.float64)
    if y0[0] < 0:
        y0 = -y0
    F = _kernel_field(p)
    Nt = np.asarray(_norm_table(F), dtype=np.int16)
    add = np.asarray(F["add"], dtype=np.int32)
    neg = np.asarray(F["neg"], dtype=np.int32)
    C = _C(p)
    Cbase = (y0[:, None] * C) * y0[None, :]
    q = F["q"]
    is_mp = bool(np.allclose(C @ y0, p * y0, atol=1e-4))
    n_cl = 0
    first = None
    first_hoff = False
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
                z = np.ones(q + 1)
                z[S] = -1.0
                y = y0 * z
                if y[0] < 0:
                    y = -y
                first_hoff = bool(np.allclose(C @ y, p * y, atol=1e-4))
                if stop_first:
                    return {
                        "Maxplus": is_mp,
                        "n_clique": n_cl,
                        "first": first,
                        "first_hoffman": first_hoff,
                    }
    return {
        "Maxplus": is_mp,
        "n_clique": n_cl,
        "first": first,
        "first_hoffman": first_hoff,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29, 37, 41):
        tags = sig_y(p, y_T(p))
        nR, nC = tags.count("R"), tags.count("C")
        rows[str(p)] = {"nR": nR, "nC": nC, "n_sq": n_sq_named(p)}
        ok = ok and nR == 3 and nC == n_sq_named(p) - 3
        if p == 37 and nC == 0:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "T is 3R+(n_sq−3)C at p=37,41. Fail: all-R at p=37.",
    }


def prove_B() -> dict:
    q = Qpp_from_L(37)
    d = Qpp_drop16(37)
    ok = q == QT37
    ok = ok and d != QT37
    ok = ok and d == Fraction(4552, 5985)
    if d == QT37:
        ok = False
    return {
        "proved": bool(ok),
        "Qpp": str(q),
        "drop16": str(d),
        "theorem": "Q_T(37)=72832/5985. Fail: drop the 16.",
    }


def prove_C() -> dict:
    p = 37
    yAP = affine_halfspace(p, (0, 1), set(range((p + 1) // 2)))
    yQR = affine_halfspace(p, (0, 1), qr0(p))
    ap = n_cliques_of(p, yAP, stop_first=False)
    qr = n_cliques_of(p, yQR, stop_first=True)
    ok = ap["Maxplus"] and qr["Maxplus"]
    ok = ok and ap["n_clique"] == 0 and ap["first"] is None
    ok = ok and qr["first"] == {"t": p, "c": p - 2}
    ok = ok and qr["first_hoffman"] is True
    ok = ok and qr["first"]["c"] != least_nsq(p)
    if ap["n_clique"] > 0 or qr["first"]["c"] == 2:
        ok = False
    return {
        "proved": bool(ok),
        "AP": {k: ap[k] for k in ("Maxplus", "n_clique", "first")},
        "QR0": qr,
        "theorem": (
            "AP⊙NC empty at p=37; QR0⊙NC first (p,p−2) is Hoffman. "
            "Fail: AP nonempty or c=ib."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "QR0⊙NC survives p=37; T/AP NC do not. Ensemble D / Q_τ "
            "unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.394  p=37 constructors: T occupancy, Q_T, QR0⊙NC", flush=True)
    A = prove_A()
    print(f"  A 3R: {A['proved']} p37={A['rows']['37']}", flush=True)
    B = prove_B()
    print(f"  B Q_T: {B['proved']} {B['Qpp']}", flush=True)
    C = prove_C()
    print(f"  C NC: {C['proved']} AP={C['AP']} QR0={C['QR0']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": {"Qpp": B["Qpp"]}, "C": {"AP": C["AP"], "QR0": C["QR0"]}}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.394",
        "title": "p=37: T occupancy and Q_T live; QR0⊙NC nonempty; AP/T NC empty",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "T_3R_at_37": A["proved"],
            "Q_T_fejer_37": B["proved"],
            "QR0_NC_survives_AP_empty": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15394.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
