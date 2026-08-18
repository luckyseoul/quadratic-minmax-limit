#!/usr/bin/env python3
"""
Prop 15.527 — Free F_q^+ orbits of H_+ split as affine-R plus a
non-R leftover: n_free = n_R + n_X with
    (n_R, n_X) = (4, 0) at p=5 and (16, 98) at p=7.
p=5 recovers 15.505 (NL is all R).  Fail: n_R=p−1 at p=7 (6≠16);
fail n_X=0 at p=7.  n_R and n_X are not p-laws.  Q_τ unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.526 flags.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15526 import is_1d_row, load_Hplus

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15527.json"


def is_affine_R(y: np.ndarray, p: int) -> bool:
    """Some F_p-parallel occupancy is an affine function of the label."""
    q = p * p
    seen: set[tuple[int, int]] = set()
    dirs: list[tuple[int, int]] = []
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0:
                continue
            g = (1, (b * pow(a, p - 2, p)) % p) if a else (0, 1)
            if g in seen:
                continue
            seen.add(g)
            dirs.append(g)
    for da, db in dirs:
        ca, cb = (0, 1) if da else (1, 0)
        occ = [0] * p
        for x in range(q):
            u, v = x % p, x // p
            lab = (ca * u + cb * v) % p
            if y[1 + x] < 0:
                occ[lab] += 1
        B = occ[0]
        A = occ[1] - occ[0]
        if all(occ[c] == A * c + B for c in range(p)):
            return True
    return False


def orbit_min_key(y: np.ndarray, p: int) -> bytes:
    q = p * p
    bits = y[1:]
    mk = None
    for t0 in range(p):
        for t1 in range(p):
            z = np.empty_like(bits)
            for x in range(q):
                a, b = x % p, x // p
                src = ((a - t0) % p) + ((b - t1) % p) * p
                z[x] = bits[src]
            kb = z.tobytes()
            if mk is None or kb < mk:
                mk = kb
    return mk


def split_free_orbits(p: int) -> dict:
    H = load_Hplus(p)
    q = p * p
    reps: dict[bytes, np.ndarray] = {}
    for i in range(len(H)):
        y = H[i]
        if is_1d_row(y, p):
            continue
        mk = orbit_min_key(y, p)
        if mk not in reps:
            reps[mk] = y
    n_R = sum(1 for y in reps.values() if is_affine_R(y, p))
    n_X = len(reps) - n_R
    n_free = (HPLUS[p] - n_1d(p)) // q
    return {
        "n_free": n_free,
        "n_reps": len(reps),
        "n_R": n_R,
        "n_X": n_X,
    }


def n_R_wrong_pm1(p: int) -> int:
    return p - 1


def prove_A() -> dict:
    rows = {}
    ok = True
    want = {5: (4, 0), 7: (16, 98)}
    for p in (5, 7):
        rec = split_free_orbits(p)
        wr, wx = want[p]
        row_ok = (
            rec["n_reps"] == rec["n_free"]
            and rec["n_R"] == wr
            and rec["n_X"] == wx
            and rec["n_R"] + rec["n_X"] == rec["n_free"]
        )
        if not row_ok:
            ok = False
        rows[str(p)] = rec
        rows[str(p)]["ok"] = row_ok
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_free=n_R+n_X with (4,0) at p=5 and (16,98) at p=7. "
            "p=5 NL is all affine-R (15.505)."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    r7 = A["by_p"]["7"]["n_R"]
    x7 = A["by_p"]["7"]["n_X"]
    ok = (
        A["proved"]
        and r7 != n_R_wrong_pm1(7)
        and x7 != 0
        and A["by_p"]["5"]["n_R"] == n_R_wrong_pm1(5)
    )
    return {
        "proved": bool(ok),
        "p7_n_R": r7,
        "p7_n_X": x7,
        "p7_pm1": n_R_wrong_pm1(7),
        "theorem": (
            "Fail n_R=p−1 at p=7 (16≠6). Fail n_X=0 at p=7. "
            "p−1 equals n_R only at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_R_named_in_p": False,
        "n_X_named_in_p": False,
        "n_free_named_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "n_R=4,16 and n_X=0,98 are not p-laws. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.527 n_free = n_R + n_X (4+0 / 16+98)", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.527",
        "title": "n_free = n_R + n_X; p=5 all R, p=7 is 16+98",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
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
