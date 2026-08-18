#!/usr/bin/env python3
"""
Prop 15.526 — F_q^+ partitions H_+ into Type+ 1D (size n_1d)
and free NL orbits (size q).  Hence
    D = n_1d/(2p) + (p/2) n_free
with n_1d named (15.292) and n_free still unnamed (4 at p=5,
114 at p=7).  Fail: n_free=0; n_free=p−1 at p=7; |H_+|=n_1d.
Does not name ensemble Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.525 flags.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _linear_forms
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15522 import D_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15526.json"

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def load_Hplus(p: int) -> np.ndarray:
    """y_∞=+1 slice of the Max+ cache."""
    Y = np.sign(np.load(YPATH[p]).astype(np.float64))
    return Y[Y[:, 0] > 0]


def is_1d_row(y: np.ndarray, p: int) -> bool:
    """True iff y_x is constant on the fibers of a nonzero F_p-linear form."""
    q = p * p
    for alpha, beta in _linear_forms(p):
        val: dict[int, float] = {}
        ok = True
        for x in range(q):
            t = (alpha * (x % p) + beta * (x // p)) % p
            bit = float(y[1 + x])
            prev = val.get(t)
            if prev is None:
                val[t] = bit
            elif prev != bit:
                ok = False
                break
        if ok and len(val) > 1:
            return True
    return False


def translation_orbits(H: np.ndarray, p: int) -> dict:
    """Orbits of F_q^+ acting by y_x ↦ y_{x−t}."""
    q = p * p
    nH = len(H)

    def key(row: np.ndarray) -> bytes:
        return np.ascontiguousarray(row.astype(np.int8)).tobytes()

    table = {key(H[i]): i for i in range(nH)}
    T1 = [((x % p + 1) % p) + (x // p) * p for x in range(q)]
    Tw = [(x % p) + ((x // p + 1) % p) * p for x in range(q)]

    def push(y: np.ndarray, perm: list[int]) -> np.ndarray:
        z = np.empty_like(y)
        z[0] = 1.0
        for x in range(q):
            z[1 + perm[x]] = y[1 + x]
        return z

    seen = np.zeros(nH, dtype=np.int8)
    sizes: list[int] = []
    n1d_in = 0
    n1d_out = 0
    for i in range(nH):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = 1
        orb: list[int] = []
        while stack:
            j = stack.pop()
            orb.append(j)
            y = H[j]
            for perm in (T1, Tw):
                z = push(y, perm)
                k = table.get(key(z))
                if k is None or seen[k]:
                    continue
                seen[k] = 1
                stack.append(k)
        sizes.append(len(orb))
        ones = sum(1 for j in orb if is_1d_row(H[j], p))
        n1d_in += ones
        n1d_out += len(orb) - ones
    c = Counter(sizes)
    n_free = int(c.get(q, 0))
    n_per = int(c.get(p, 0))
    return {
        "nH": nH,
        "n_orbits": len(sizes),
        "sizes": {str(k): int(v) for k, v in sorted(c.items())},
        "n_free_orbits": n_free,
        "n_periodic_orbits": n_per,
        "n_1d_in_H": n1d_in,
        "n_non1d_in_H": n1d_out,
    }


def n_free_wrong_zero() -> int:
    return 0


def n_free_wrong_pm1(p: int) -> int:
    return p - 1


def D_drop_nfree(p: int) -> float:
    return n_1d(p) / (2 * p)


def D_from_nfree(p: int, n_free: int) -> float:
    return n_1d(p) / (2 * p) + (p / 2) * n_free


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        H = load_Hplus(p)
        rec = translation_orbits(H, p)
        want_1d = n_1d(p)
        q = p * p
        row_ok = (
            rec["nH"] == HPLUS[p]
            and rec["n_1d_in_H"] == want_1d
            and rec["n_non1d_in_H"] == HPLUS[p] - want_1d
            and rec["n_periodic_orbits"] * p == want_1d
            and rec["n_free_orbits"] * q == HPLUS[p] - want_1d
            and rec["nH"] != want_1d
        )
        if not row_ok:
            ok = False
        rows[str(p)] = rec
        rows[str(p)]["n_1d_named"] = want_1d
        rows[str(p)]["ok"] = row_ok
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "F_q^+ on H_+: Type+ 1D is the size-p orbits (count n_1d); "
            "NL is size-q free orbits. Fail |H_+|=n_1d."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    ok = A["proved"]
    rows = {}
    for p in (5, 7):
        nf = A["by_p"][str(p)]["n_free_orbits"]
        D = D_from_nfree(p, nf)
        drop = D_drop_nfree(p)
        live = D_named(p)
        # p−1 equals live n_free at p=5 (4=4) and dies at p=7 (6≠114)
        row_ok = (
            abs(D - live) < 1e-12
            and nf != n_free_wrong_zero()
            and abs(drop - live) > 0.5
        )
        if p == 7 and nf == n_free_wrong_pm1(p):
            row_ok = False
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_free": nf,
            "D_from_nfree": D,
            "D_live": live,
            "D_drop_nfree": drop,
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "D=n_1d/(2p)+(p/2)n_free. Fail n_free=0 (D=3 at p=5); "
            "fail n_free=p−1 at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "n_free_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "n_free is 4 at p=5 and 114 at p=7, not a named p-law. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.526 H_+ = n_1d + q n_free; D = n_1d/(2p)+(p/2)n_free", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.526",
        "title": "H_+ = n_1d + q n_free; D named in n_free",
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
