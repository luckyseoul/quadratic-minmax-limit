#!/usr/bin/env python3
"""
Prop 15.494 — The 3-linear A_r of 1̂_D on Ω vanishes iff χ(r+1)=−1.

    A_r = E[1̂_D(ξ) 1̂_D(rξ) 1̂_D(−(1+r)ξ)],   ξ∈Ω, r∈T\\{±1}.

On Max+, ẑ is supported on {0}∪Ω and 1̂_D(α)=−ẑ(α)/2 for α≠0.
χ_q(−1)=1, so −(1+r)ξ∈Ω ⇔ χ(r+1)=1.  Thus χ(r+1)=−1 ⇒
−(1+r)ξ∉Ω∪{0} ⇒ 1̂_D(−(1+r)ξ)=0 pointwise ⇒ A_r=0.

Fail-when-wrong: claim A_r=0 on χ(r+1)=+1 (p=5 Paley ++ sub
has A_r=125/26≠0); invert to χ(r−1)=−1 (p=5 (−1,1) N1 is
625/26≠0).  Does not name Q_τ (Q is constant on ++ at p=7
while A_r spreads on ++ sub).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.493 flags.
"""
from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15290 import omega_set, paley_pair, type_key

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15494.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def minus_one_plus_r(F, r: int) -> int:
    """−(1+r) in F_q."""
    return F["neg"][F["add"][1][r]]


def in_Omega(F, Omega_set: set[int], alpha: int) -> bool:
    return alpha != 0 and alpha in Omega_set


def support_vanishes(F, r: int, xi: int, Omega_set: set[int]) -> bool:
    """True iff 1̂_D(−(1+r)ξ) must vanish by supp ẑ ⊂ {0}∪Ω."""
    sxi = F["mul"][minus_one_plus_r(F, r)][xi]
    return not in_Omega(F, Omega_set, sxi)


def prove_A() -> dict:
    """Max+-free: χ(r+1)=−1 ⇔ −(1+r)Ω ∩ Ω = ∅."""
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        Omega = omega_set(F)
        Oset = set(Omega)
        xi = Omega[0]
        n_m1 = 0
        n_p1 = 0
        bad_m1 = 0
        bad_p1 = 0
        for r in F["squares"]:
            if r in (1, F["neg1"]):
                continue
            _a, ch_rp1 = paley_pair(F, r)
            vanishes = support_vanishes(F, r, xi, Oset)
            # also check every xi in Ω
            all_v = all(support_vanishes(F, r, eta, Oset) for eta in Omega)
            if ch_rp1 == -1:
                n_m1 += 1
                if (not vanishes) or (not all_v):
                    bad_m1 += 1
                    ok = False
            else:
                n_p1 += 1
                if vanishes or all_v:
                    bad_p1 += 1
                    ok = False
        rows[p] = {
            "n_chi_rp1_m1": n_m1,
            "n_chi_rp1_p1": n_p1,
            "bad_m1": bad_m1,
            "bad_p1": bad_p1,
        }
        if n_m1 == 0 or n_p1 == 0:
            ok = False
    return {"proved": bool(ok), "by_p": rows}


def _hats(F, D, alphas):
    mul = F["mul"]
    out = {}
    for a in alphas:
        e = np.array([_e_tr(F, mul[a][x]) for x in range(F["q"])])
        out[a] = (D.astype(np.float64) * e).sum(axis=1)
    return out


def live_A3_table(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    D = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0
    Omega = omega_set(F)
    xi = Omega[0]
    Toff = [r for r in F["squares"] if r not in (1, F["neg1"])]
    need = {xi}
    for r in Toff:
        need.add(F["mul"][r][xi])
        need.add(F["mul"][minus_one_plus_r(F, r)][xi])
    hats = _hats(F, D, need)
    by = defaultdict(list)
    for r in Toff:
        rxi = F["mul"][r][xi]
        sxi = F["mul"][minus_one_plus_r(F, r)][xi]
        Ar = np.mean(hats[xi] * hats[rxi] * hats[sxi]).real
        _a, ch = paley_pair(F, r)
        by[(ch, type_key(F, r))].append(float(Ar))
    packed = {}
    for (ch, tk), vs in by.items():
        a = np.asarray(vs, dtype=np.float64)
        packed[str(tk)] = {
            "chi_rp1": int(ch),
            "n": len(vs),
            "mean": float(a.mean()),
            "spread": float(a.max() - a.min()),
        }
    return packed


def prove_B() -> dict:
    """Certified p=5,7: live A_r=0 on χ(r+1)=−1; nonzero on a ++ type."""
    ok = True
    rows = {}
    for p in (5, 7):
        tab = live_A3_table(p)
        rows[p] = tab
        for rec in tab.values():
            if rec["chi_rp1"] == -1:
                if abs(rec["mean"]) > 1e-8 or rec["spread"] > 1e-8:
                    ok = False
        # fail-when-wrong: vanish on ++ (χ(r+1)=+1)
        plus = [rec for rec in tab.values() if rec["chi_rp1"] == 1]
        if not plus or all(abs(rec["mean"]) < 1e-8 for rec in plus):
            ok = False
    # explicit p=5 witnesses
    t5 = rows[5]
    if abs(t5["(1, 1, 'sub')"]["mean"] - 125 / 26) > 1e-6:
        ok = False
    if abs(t5["(-1, 1, 'N1')"]["mean"] - 625 / 26) > 1e-6:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_imported": False,
        "Q_tau_named": False,
        "note": (
            "A_r=0 on χ(r+1)=−1 is support of ẑ, not a name of Q_τ. "
            "At p=7 A_r spreads on ++ sub while Q++ is constant. "
            "⟨δ,ψ⟩≤2 stays OPEN."
        ),
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
