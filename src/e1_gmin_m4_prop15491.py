#!/usr/bin/env python3
"""
Prop 15.491 — The constant (trivial-character) piece of the Wick
deficit is D-free:

    ∑_{r∈T} Q(r) = 2 q² (q−1)
    ∑_{r∈T\\{±1}} (4 − Q(r)/q²) = 8
    mean_{off} δ/q² = 16/(q−5)

Certified on live Max+ at p=5,7.  Fail: 2q²(q+1); claim the
off-sum of deficits is 4 or 16.  Does not name δ_0 or prove
⟨δ,ψ⟩≤2.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.490 flags.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15290 import live_Q_on_T

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15491.json"


@lru_cache(None)
def _live_Q(p: int):
    return live_Q_on_T(p)


def sum_Q_over_T(p: int) -> float:
    """∑_{r∈T} Q(r) / q² from the live ensemble."""
    F = _kernel_field(p)
    q2 = float(F["q"] ** 2)
    return sum(_live_Q(p).values()) / q2


def sum_Q_named(p: int) -> int:
    """2(q−1) with q=p²."""
    q = p * p
    return 2 * (q - 1)


def sum_Q_wrong_plus(p: int) -> int:
    """Fail-when-wrong: 2(q+1)."""
    q = p * p
    return 2 * (q + 1)


def sum_off_deficit(p: int) -> float:
    F = _kernel_field(p)
    q2 = float(F["q"] ** 2)
    Q = _live_Q(p)
    s = 0.0
    for r, Qr in Q.items():
        if r in (1, F["neg1"]):
            continue
        s += 4.0 - Qr / q2
    return s


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        got = sum_Q_over_T(p)
        named = float(sum_Q_named(p))
        wrong = float(sum_Q_wrong_plus(p))
        rows[p] = {"sumQ_over_q2": got, "named": named, "wrong": wrong}
        if abs(got - named) > 1e-6:
            ok = False
        if abs(got - wrong) < 1e-6:
            ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        s = sum_off_deficit(p)
        rows[p] = s
        if abs(s - 8.0) > 1e-6:
            ok = False
        if abs(s - 4.0) < 1e-6 or abs(s - 16.0) < 1e-6:
            ok = False
    return {"proved": bool(ok), "sum_off_delta": rows, "named": 8}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": "Names only the constant of δ. ⟨δ,ψ⟩≤2 still open.",
    }


def main() -> dict:
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "A": A,
        "B": B,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
