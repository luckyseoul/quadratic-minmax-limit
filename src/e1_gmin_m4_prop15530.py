#!/usr/bin/env python3
"""
Prop 15.530 — n_R=(p−3)² at p=3,5,7.  Fail 2^{p−3} (1≠0 at p=3).
Fail n_R=p−1 (2≠0 at p=3; 6≠16 at p=7).  Not proved for p≥11.
Does not name Q_τ in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.529 flags.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15526 import is_1d_row
from e1_gmin_m4_prop15527 import is_affine_R, orbit_min_key

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15530.json"

YPATH = {
    3: "/tmp/maxplus_p3.npy",
    5: "/tmp/maxplus_p5.npy",
    7: "/tmp/maxplus_p7.npy",
}


def load_Hplus(p: int) -> np.ndarray:
    Y = np.sign(np.load(YPATH[p]).astype(np.float64))
    return Y[Y[:, 0] > 0]


def n_R_named(p: int) -> int:
    return (p - 3) ** 2


def n_R_wrong_pow2(p: int) -> int:
    """Fail-when-wrong interpolant 2^{p−3}."""
    return 2 ** (p - 3)


def n_R_wrong_pm1(p: int) -> int:
    return p - 1


def n_R_live(p: int) -> int:
    H = load_Hplus(p)
    reps: set[bytes] = set()
    for i in range(len(H)):
        y = H[i]
        if is_1d_row(y, p) or not is_affine_R(y, p):
            continue
        reps.add(orbit_min_key(y, p))
    return len(reps)


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (3, 5, 7):
        live = n_R_live(p)
        named = n_R_named(p)
        row_ok = live == named
        if p == 3:
            row_ok = row_ok and live == 0 and n_1d(p) == len(load_Hplus(p))
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "live": live,
            "named": named,
            "pow2": n_R_wrong_pow2(p),
            "pm1": n_R_wrong_pm1(p),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "n_R=(p−3)² at p=3,5,7 (H_+ affine-R free orbits).",
    }


def prove_B() -> dict:
    """Fail 2^{p−3} and p−1."""
    ok = (
        n_R_named(3) != n_R_wrong_pow2(3)
        and n_R_live(3) != n_R_wrong_pow2(3)
        and n_R_named(3) != n_R_wrong_pm1(3)
        and n_R_named(7) != n_R_wrong_pm1(7)
        and n_R_named(5) == n_R_wrong_pm1(5)
        and n_R_named(5) == n_R_wrong_pow2(5)
        and n_R_named(7) == n_R_wrong_pow2(7)
    )
    return {
        "proved": bool(ok),
        "p3_pow2": n_R_wrong_pow2(3),
        "p3_named": n_R_named(3),
        "p7_pm1": n_R_wrong_pm1(7),
        "p7_named": n_R_named(7),
        "theorem": (
            "Fail 2^{p−3} (p=3: 1≠0). Fail p−1 (p=3: 2≠0; p=7: 6≠16). "
            "Those interpolants agree with (p−3)² only on {5,7}."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_R_general": False,
        "n_X_named_in_p": False,
        "n_free_named_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "(p−3)² is certified at p=3,5,7, not a bijection for p≥11. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.530 n_R=(p-3)^2 at p=3,5,7; fail 2^{p-3} at p=3", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.530",
        "title": "n_R=(p−3)² at p=3,5,7; 2^{p−3} dies at p=3",
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
