#!/usr/bin/env python3
"""
Prop 15.497 — J_{N*}(ψ_k) is named by divisibility, and on a
two-type (p≡1) ensemble the pairing is

    ⟨δ,ψ_k⟩ = −8 + 16/D · (A + (p−4) J_{N*}(k))

with A from 15.367 and D=|H+|/(2p).  Here Q_{++}/q²=8A/D
and Q_{N*}/q²=8(A−2(p−4))/D as in 15.490 (not the
un-normalized Q).  δ=4−Q/q².  Worst even pairing is
−8+16(A+2p−8)/D, ≤2 iff D≥8(A+2p−8)/5.  Live D=13≥64/5
at p=5 is a spot-check, not a name of D.  J_{N*} is
Max+-free.  D unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.496 flags.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec
from e1_gmin_m4_prop15290 import field_norm, in_subfield, live_Q_on_T
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15367 import A_num

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15497.json"

CERT_J = (5, 7, 11, 13, 17, 19)


def even_ks(p: int) -> list[int]:
    q1 = p * p - 1
    half = q1 // 2
    return [k for k in range(2, q1, 2) if k != half]


def J_Nstar_named(p: int, k: int) -> int:
    """Restriction of ψ_k to F_p^× and to U={N=1}.

    ∑_{F_p^×} ψ = p−1 if (p−1)|k else 0.
    ∑_{U} ψ     = p+1 if (p+1)|k else 0.
    ∑_{N1} ψ    = ∑_U ψ − 2  (even ⇒ ψ(±1)=1).
    T = F_p^× ⊔ N1 ⊔ N* and ∑_T ψ = 0, so
        J_{N*} = −∑_{F_p^×} ψ − ∑_{N1} ψ.
    Exclusive for even k ∉ {0,(q−1)/2}: lcm(p−1,p+1)=(q−1)/2.
    """
    if k % (p - 1) == 0:
        return -(p - 3)
    if k % (p + 1) == 0:
        return -(p - 1)
    return 2


def J_Nstar_swapped(p: int, k: int) -> int:
    """Fail-when-wrong: swap the (p−1)|(p+1) cases."""
    if k % (p - 1) == 0:
        return -(p - 1)
    if k % (p + 1) == 0:
        return -(p - 3)
    return 2


def nstar_set(F) -> list[int]:
    out = []
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        if in_subfield(F, r):
            continue
        if field_norm(F, r) == 1:
            continue
        out.append(int(r))
    return out


def J_Nstar_live(F, k: int, ns: list[int]) -> float:
    psi = _psi_vec(F, k)
    return float(sum(psi[r] for r in ns).real)


def pair_named(p: int, k: int, D: int) -> Fraction:
    """p≡1 two-type identity."""
    A = A_num(p)
    j = J_Nstar_named(p, k)
    return Fraction(-8) + Fraction(16 * (A + (p - 4) * j), D)


def pair_signflip(p: int, k: int, D: int) -> Fraction:
    """Fail-when-wrong: flip the sign of (p−4) J_{N*}."""
    A = A_num(p)
    j = J_Nstar_named(p, k)
    return Fraction(-8) + Fraction(16 * (A - (p - 4) * j), D)


def max_pair_named(p: int, D: int) -> Fraction:
    """Worst even pairing: J_{N*}=2."""
    A = A_num(p)
    return Fraction(-8) + Fraction(16 * (A + 2 * (p - 4)), D)


def D_need_named(p: int) -> Fraction:
    """D ≥ 8(A+2p−8)/5  ⇔  max pairing ≤ 2."""
    return Fraction(8 * (A_num(p) + 2 * p - 8), 5)


def D_live(p: int) -> int:
    return HPLUS[p] // (2 * p)


def live_pair(p: int, k: int) -> float:
    F = _kernel_field(p)
    Q = live_Q_on_T(p)
    q2 = float(F["q"] ** 2)
    psi = _psi_vec(F, k)
    s = 0.0
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        s += (4.0 - Q[r] / q2) * psi[r].real
    return s


def prove_A() -> dict:
    """J_{N*} named. Fail: swapped divisibility."""
    ok = True
    rows = {}
    for p in CERT_J:
        F = _kernel_field(p)
        ns = nstar_set(F)
        rec = {}
        swap_hit = False
        for k in even_ks(p):
            got = J_Nstar_live(F, k, ns)
            named = J_Nstar_named(p, k)
            rec[str(k)] = {"live": got, "named": named}
            if abs(got - named) > 1e-8:
                ok = False
            if abs(got - J_Nstar_swapped(p, k)) > 1e-6:
                swap_hit = True
        if not swap_hit:
            ok = False
        rows[p] = {
            "n_N*": len(ns),
            "values": sorted({int(round(v["live"])) for v in rec.values()}),
            "named_set": sorted({2, -(p - 3), -(p - 1)}),
        }
        if rows[p]["values"] != rows[p]["named_set"]:
            ok = False
    # explicit p=5 k=4: (p−1)|4 ⇒ −(p−3)=−2, not −4
    if J_Nstar_named(5, 4) != -2 or J_Nstar_swapped(5, 4) != -4:
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """Pairing identity at p=5 (two-type). Fail: flip (p−4) J sign."""
    ok = True
    p, D = 5, D_live(5)
    rows = {}
    flip_hit = False
    max_live = -1e9
    for k in even_ks(p):
        lv = live_pair(p, k)
        nm = float(pair_named(p, k, D))
        rows[str(k)] = {"live": lv, "named": nm}
        if abs(lv - nm) > 1e-8:
            ok = False
        if abs(lv - float(pair_signflip(p, k, D))) > 1e-6:
            flip_hit = True
        max_live = max(max_live, lv)
    if not flip_hit:
        ok = False
    mx = max_pair_named(p, D)
    if abs(max_live - float(mx)) > 1e-8:
        ok = False
    if mx != Fraction(24, 13):
        ok = False
    return {
        "proved": bool(ok),
        "max_live": max_live,
        "max_named": str(mx),
        "n_checked": len(rows),
    }


def prove_C() -> dict:
    """max pairing ≤2 ⇔ D ≥ 8(A+2p−8)/5. Live D at p=5. Not a D name."""
    ok = True
    p, D = 5, D_live(5)
    need = D_need_named(p)
    if need != Fraction(64, 5):
        ok = False
    if D < need:
        ok = False
    # fail: claim D_1d=3 suffices
    if 3 >= float(need):
        ok = False
    mx = max_pair_named(p, D)
    if mx > 2:
        ok = False
    if max_pair_named(p, 3) <= 2:
        ok = False
    return {
        "proved": bool(ok),
        "D_live": D,
        "D_need": str(need),
        "max_pair": str(mx),
        "D_1d_too_small": True,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "J_{N*} named; pairing reduced to D. D=|H+|/(2p) unnamed. "
            "p≡3 still has a --N1 term. phi_F not imported."
        ),
    }


def main() -> dict:
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "A": A,
        "B": B,
        "C": C,
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
