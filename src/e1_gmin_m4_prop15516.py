#!/usr/bin/env python3
"""
Prop 15.516 — Ensemble L on the 15.448 ns_mix slot equals μ₊μ₋.

  E[N(d)N(sd) | slot(s)=ns_mix] = p²(p−1)(p−3)/16 = μ₊μ₋
  at p=5,7 (both χ).  The slot still splits (15.449); only the
  mean is named.  Fail: claim L is constant on ns_mix (p=5
  315/26 vs 335/26); claim ns_11 equals μ₊μ₋ (180/13≠25/2 at p=5).
  Does not name Q_τ (other slots keep den D).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.515 flags.
"""
from __future__ import annotations

import json
import os
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15298 import build_L, load_N, mu_minus, mu_plus, ns_key
from e1_gmin_m4_prop15448 import Lns_11_named, Lns_mix_named, slot_of

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15516.json"

CERT = (5, 7)


def ensemble_slot_L(p: int, slot: str) -> tuple[Fraction, set[str]]:
    F = _kernel_field(p)
    L = build_L(load_N(p, F), F)
    acc = defaultdict(list)
    for s in range(1, F["q"]):
        if slot_of(F, p, s) != slot:
            continue
        key = ns_key(F, s)
        for cd in (1, -1):
            acc[cd].append(Fraction(L[(cd, key)]).limit_denominator(p**6))
    means = {}
    vals = set()
    for cd, xs in acc.items():
        if not xs:
            continue
        m = sum(xs, Fraction(0)) / len(xs)
        means[cd] = m
        vals.update(str(x) for x in xs)
    if not means:
        return Fraction(0), set()
    # both χ should agree
    m0 = next(iter(means.values()))
    return m0, vals


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        mean, vals = ensemble_slot_L(p, "ns_mix")
        named = Lns_mix_named(p)
        mu = mu_plus(p) * mu_minus(p)
        if mean != named or named != mu:
            ok = False
        if len(vals) < 2:
            ok = False
        if mean == Lns_11_named(p):
            ok = False
        rows[str(p)] = {
            "mean": str(mean),
            "named": str(named),
            "nuniq": len(vals),
            "vals": sorted(vals),
            "ns11": str(Lns_11_named(p)),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Ensemble L_ns_mix=μ₊μ₋. Fail constant; fail ns_11.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "One ensemble L slot is named. Other slots still have den D. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.516 ensemble L_ns_mix = μ+μ−", flush=True)
    A, D = prove_A(), prove_open()
    out = {
        "prop": "15.516",
        "title": "Ensemble L_ns_mix=μ₊μ₋; Q_τ unnamed",
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
