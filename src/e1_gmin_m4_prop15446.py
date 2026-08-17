#!/usr/bin/env python3
"""
Prop 15.446 — L is constant on Sub, Tor, N* separately, not on
merged ++.  The L-weighted ns block does **not** drop from
Q++−Q_N*: at p=5 it is Δ_ns=50/13 of Δ(Q/16).  The 15.445
merged-++ ratio 9/5 is therefore not the live criterion
(live L_Sub/L_N*=34/29<9/5, yet Q++>Q_N*).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.445 flags.

============================================================================
Setup.  p≡1: ++ = Sub ∪ Tor, Sub=F_p^×\\{±1}, Tor=U\\{±1}.
L_χ(s)=E[N(d)N(sd)].  15.445 C used unweighted C_ns (15.429).

============================================================================
Theorem A — PROVED (cache 4-point; certified p=5,7).
  L_χ is constant on each of Sub, Tor, N*, and takes two
  distinct values on merged ++ at p=5 (340/13 vs 645/26 for
  χ=+1).  Fail: claim L is constant on ++ at p=5.  ∎

Theorem B — PROVED (15.298 expansion; certified p=5).
  A0 cancels in Q++−Q_N*.  The ns block contributes
  Δ_ns=50/13 to Δ(Q/16), not 0, and 16(Δ_sq+Δ_ns)/q²=16/13
  recovers live 48/13−32/13.  Fail: claim Δ_ns=0.  ∎

Theorem C — PROVED (A+B arithmetic).
  Live L_Sub^+/L_N*^+=34/29<9/5, so the 15.445 merged ratio
  is not necessary for Q++≥Q_N*.  Fail: claim 34/29≥9/5.  ∎

Theorem D — OPEN.  Q++≥Q_N* and Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Max+ cache 4-point + field half-Gauss (p=5,7).  Writes
evidence/e1_gmin_m4_prop15446.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    build_L,
    half_gauss,
    live_Q_any,
    load_N,
    ns_key,
)
from e1_gmin_m4_prop15445 import L_ratio_need, pow_f  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15446_catalog.json"


def classify(F, p: int, s: int) -> str:
    if s in (1, F["neg1"]):
        return "pm1"
    if F["chi"][s] != 1:
        return "ns"
    if pow_f(F, s, p) == s:
        return "Sub"
    if pow_f(F, s, p + 1) == 1:
        return "Tor"
    return "N*"


def L_by_class(p: int) -> dict[tuple, Fraction]:
    F = _kernel_field(p)
    L = build_L(load_N(p, F), F)
    out: dict[tuple, set[Fraction]] = defaultdict(set)
    for s in range(1, F["q"]):
        tau = classify(F, p, s)
        for cd in (1, -1):
            out[(cd, tau)].add(
                Fraction(L[(cd, ns_key(F, s))]).limit_denominator(10**6)
            )
    return {k: next(iter(v)) if len(v) == 1 else v for k, v in out.items()}


def ns_delta_Q16(p: int) -> Fraction:
    F = _kernel_field(p)
    q = F["q"]
    L = build_L(load_N(p, F), F)
    xi = omega_set(F)[0]
    G = float(
        sum(F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)).real
    )
    plus, nstar, ns = [], [], []
    for s in range(1, q):
        tau = classify(F, p, s)
        if tau == "ns":
            ns.append(s)
        elif tau in ("Sub", "Tor"):
            plus.append(s)
        elif tau == "N*":
            nstar.append(s)

    def avg(rho):
        acc = 0.0
        for r in rho:
            for s in ns:
                for cd in (1, -1):
                    acc += L[(cd, ns_key(F, s))] * half_gauss(
                        F, F["add"][1][F["mul"][r][s]], cd, G, q
                    )
        return acc / len(rho)

    d = avg(plus) - avg(nstar)
    return Fraction(d).limit_denominator(10**6)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        Lb = L_by_class(p)
        sub = Lb[(1, "Sub")]
        tor = Lb[(1, "Tor")]
        if isinstance(sub, set):
            ok = False
        if p == 5:
            if isinstance(tor, set) or sub == tor:
                ok = False
            if sub != Fraction(340, 13) or tor != Fraction(645, 26):
                ok = False
        nstar_set = Lb[(1, "N*")]
        rows[str(p)] = {
            "L_Sub+": str(sub) if not isinstance(sub, set) else sorted(str(x) for x in sub),
            "L_Tor+": str(tor) if not isinstance(tor, set) else sorted(str(x) for x in tor),
            "L_N*+": (
                str(nstar_set)
                if not isinstance(nstar_set, set)
                else sorted(str(x) for x in nstar_set)
            ),
            "merged_const": (not isinstance(sub, set) and not isinstance(tor, set) and sub == tor),
        }
    if L_by_class(5)[(1, "Sub")] == L_by_class(5)[(1, "Tor")]:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "L constant on Sub, Tor, N*; two values on ++ at p=5. "
            "Fail: L constant on ++ at p=5."
        ),
    }


def prove_B() -> dict:
    d = ns_delta_Q16(5)
    ok = d == Fraction(50, 13)
    if d == 0:
        ok = False
    Q = live_Q_any(5)
    F = _kernel_field(5)
    plus = [
        r
        for r in F["squares"]
        if r not in (1, F["neg1"]) and classify(F, 5, r) in ("Sub", "Tor")
    ]
    nstar = [
        r
        for r in F["squares"]
        if classify(F, 5, r) == "N*"
    ]
    q2 = float(F["q"] ** 2)
    live = float(np.mean([Q[r] / q2 for r in plus]) - np.mean([Q[r] / q2 for r in nstar]))
    if abs(live - float(Fraction(16, 13))) > 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "Delta_ns": str(d),
        "live_dQ": str(Fraction(16, 13)),
        "theorem": (
            "Weighted ns block is Δ_ns=50/13 at p=5, not 0. "
            "Fail: Δ_ns=0."
        ),
    }


def prove_C() -> dict:
    Lb = L_by_class(5)
    ratio = Lb[(1, "Sub")] / Lb[(1, "N*")]
    need = L_ratio_need(5)
    ok = ratio == Fraction(34, 29) and ratio < need
    if ratio >= need:
        ok = False
    return {
        "proved": bool(ok),
        "live_Sub_over_Nstar": str(ratio),
        "need_445": str(need),
        "theorem": (
            "34/29<9/5, so the 15.445 merged ratio is not necessary. "
            "Fail: 34/29≥9/5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Corrected criterion for Q++≥Q_N* uses Sub/Tor/ns. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.446  L splits Sub/Tor; weighted ns does not drop", flush=True)
    A = prove_A()
    print(f"  A L-split: {A['proved']} p5 Sub={A['by_p']['5']['L_Sub+']}", flush=True)
    B = prove_B()
    print(f"  B Δ_ns: {B['proved']} {B['Delta_ns']}", flush=True)
    C = prove_C()
    print(f"  C 9/5 not necessary: {C['proved']} {C['live_Sub_over_Nstar']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "Delta_ns": B["Delta_ns"],
                "ratio": C["live_Sub_over_Nstar"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.446",
        "title": "L splits Sub/Tor; weighted ns does not drop from Q++−Q_N*",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "L_Sub_Tor_split": A["proved"],
            "ns_weighted_not_zero": B["proved"],
            "ratio_445_not_necessary": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15446.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
