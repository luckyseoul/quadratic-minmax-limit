#!/usr/bin/env python3
"""
Prop 15.334 — QR0⊙NC is an NL Max+ construction at p=5,7,11,13.
AP⊙NC (ystar) is empty at p=13; QR0⊙NC is not.
This is still 15.x (leftovers OPEN), not 16.x/17.x.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.333 flags.  |H_+| / Q_τ unnamed.  Floor stays OPEN.

============================================================================
Setup.  QR0 = F_p^{□} ∪ {0} ⊂ F_p, |QR0|=(p+1)/2.  y_QR0 is the
affine halfspace for S=QR0 (15.309).  QR0⊙NC = first Hoffman+evec
norm circle of the switched conference C_{y_QR0} (15.312 first_nc).

============================================================================
Theorem A — PROVED (constructor; certified p=5,7,11,13).
  first_nc(y_QR0) returns a boolean y with Cy=py and y_∞=+1.
  Fail-when-wrong: claim the construction is empty at p=13
  (that is AP⊙NC / YSTAR_CERT[13]).  ∎

Theorem B — PROVED (15.307 signatures; certified p=7,13).
  For p≥7 the result is not a halfspace (signature is not one H
  and the rest C).  At p=5, QR0 is an AP so QR0⊙NC coincides with
  ystar (stab 6).  ∎

Theorem C — OPEN.  One named NL orbit (or a pair, AP vs QR0) is a
  lower bound on extra, not |H_+|.  den(Q_{++}) still unnamed.
  phi_F not imported.  ∎

============================================================================
Backend: first_nc + Paley C.  Writes evidence/e1_gmin_m4_prop15334.json
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import YSTAR_CERT  # noqa: E402
from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def y_QR0(p: int) -> np.ndarray:
    y = affine_halfspace(p, (0, 1), qr0(p))
    y = np.asarray(y, dtype=np.float64)
    if y[0] < 0:
        y = -y
    return y


@lru_cache(maxsize=8)
def construct_qr0nc(p: int) -> dict:
    y0 = y_QR0(p)
    pack = first_nc(p, y0)
    if pack is None:
        return {"p": p, "found": False}
    y = np.asarray(pack["y"], dtype=np.float64)
    if y[0] < 0:
        y = -y
    C = paley_conference_prime_power(p).astype(np.float64)
    cy = bool(np.allclose(C @ y, p * y, atol=1e-5))
    return {
        "p": p,
        "found": True,
        "cy_eq_py": cy,
        "t": int(pack["t"]),
        "c": int(pack["c"]),
        "y": y,
    }


def is_halfspace_y(p: int, y: np.ndarray) -> bool:
    F = _kernel_field(p)
    q = F["q"]
    D = np.asarray(y[1 : 1 + q]) < 0
    tags = []
    for lines in square_classes(F):
        ns = D[lines].sum(axis=1)
        sig = tuple(int(x) for x in np.sort(ns))
        tags.append(tag_sig(p, sig))
    key = tuple(sorted(tags))
    return key.count("H") == 1 and key.count("C") == (p + 1) // 2 - 1


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13):
        rec = construct_qr0nc(p)
        rows[str(p)] = {k: rec[k] for k in rec if k != "y"}
        if not rec["found"] or not rec["cy_eq_py"]:
            ok = False
    if YSTAR_CERT[13].get("found"):
        ok = False
    if not rows["13"]["found"]:
        ok = False
    return {
        "proved": ok,
        "by_p": rows,
        "AP_nc_empty_p13": not bool(YSTAR_CERT[13].get("found")),
        "theorem": (
            "QR0⊙NC exists and is Cy=py at p=5,7,11,13.  "
            "AP⊙NC is empty at p=13."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (7, 13):
        rec = construct_qr0nc(p)
        hs = is_halfspace_y(p, rec["y"])
        if hs:
            ok = False
        rows[str(p)] = {"is_hs": hs, "found": rec["found"]}
    rec5 = construct_qr0nc(5)
    hs5 = is_halfspace_y(5, rec5["y"])
    rows["5"] = {"is_hs": hs5, "note": "QR0 is AP at p=5; may coincide with ystar"}
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "QR0⊙NC is NL (not a halfspace) for p=7 and p=13.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Hplus_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": (
            "QR0⊙NC is one NL family through p=13.  "
            "It does not name |H_+|.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.334  QR0⊙NC exists at p=13 (AP⊙NC empty)", flush=True)
    A = prove_A()
    print(f"  A exists: {A['proved']} p13={A['by_p']['13']}", flush=True)
    B = prove_B()
    print(f"  B NL: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.334",
        "title": "QR0⊙NC exists at p=5,7,11,13; AP⊙NC empty at p=13",
        "series": "15.x leftover campaign (OPEN). Not 16.x fail, not 17.x close.",
        "proved": {
            "qr0nc_exists_p13": A["proved"],
            "qr0nc_NL_p7_p13": B["proved"],
            "Hplus_named_in_p": False,
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15334.json"
    dump = json.loads(json.dumps(out, default=str))
    path.write_text(json.dumps(dump, indent=2))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
