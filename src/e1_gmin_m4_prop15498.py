#!/usr/bin/env python3
"""
Prop 15.498 — p≡3 pairing term.  On N1, -- ∩ N1 is nonempty
iff p≡3, and

    J_{--}(ψ_k) = (p+1)/2   if (p+1)|k
                = −(p+1)/2  if k ≡ (p+1)/2 (mod p+1)
                = 0         otherwise.

On the live p=7 three-type ensemble, Q_{--}/q²=8(A−3p)/D and

    ⟨δ,ψ_k⟩ = −8 + 16/D (A + (p−4) J_{N*} + (3p/2) J_{--}).

Max+-free J_{--}.  A−3p is certified at p=7, not a general
name of A_{--}.  D unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.497 flags.
"""
from __future__ import annotations

import json
import os
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec
from e1_gmin_m4_prop15290 import field_norm, in_subfield, live_Q_on_T, paley_pair
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15497 import J_Nstar_named, even_ks

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15498.json"

CERT_J = (7,)


def J_mm_named(p: int, k: int) -> int:
    """Character sum of ψ_k on -- ∩ N1, p≡3.

    j=k mod (p+1) is ψ|U.  j=0: ψ trivial on U, J=|--N1|=(p+1)/2
    (15.355).  j=(p+1)/2: ψ|U is the quadratic character of U,
    equal to χ(u−1) on N1 when p≡3 (mix empty), hence
    J=−|--N1|.  Otherwise the even U-characters sum to 0 on --.
    p≡1: --N1 empty, J=0.
    """
    if p % 4 != 3:
        return 0
    m = p + 1
    if k % m == 0:
        return m // 2
    if k % m == m // 2:
        return -(m // 2)
    return 0


def J_mm_wrong_half(p: int, k: int) -> int:
    """Fail-when-wrong: (p−1)/2 in place of (p+1)/2."""
    if p % 4 != 3:
        return 0
    m = p + 1
    h = (p - 1) // 2
    if k % m == 0:
        return h
    if k % m == m // 2:
        return -h
    return 0


def mm_n1_set(F) -> list[int]:
    out = []
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        if field_norm(F, r) != 1 or in_subfield(F, r):
            continue
        a, b = paley_pair(F, r)
        if a == -1 and b == -1:
            out.append(int(r))
    return out


def J_mm_live(F, k: int, rs: list[int]) -> float:
    psi = _psi_vec(F, k)
    return float(sum(psi[r] for r in rs).real)


def D_live(p: int) -> int:
    return HPLUS[p] // (2 * p)


def pair_p3_named(p: int, k: int, D: int) -> Fraction:
    """Certified p=7: A_{--}=A−3p.  Q/q² as in 15.490."""
    A = A_num(p)
    M = A + (p - 4) * J_Nstar_named(p, k) + Fraction(3 * p, 2) * J_mm_named(p, k)
    return Fraction(-8) + Fraction(16, D) * M


def pair_p3_A2p(p: int, k: int, D: int) -> Fraction:
    """Fail-when-wrong: A_{--}=A−2p."""
    A = A_num(p)
    M = A + (p - 4) * J_Nstar_named(p, k) + p * J_mm_named(p, k)
    return Fraction(-8) + Fraction(16, D) * M


def live_pair(p: int, k: int, F=None, Q=None) -> float:
    if F is None:
        F = _kernel_field(p)
    if Q is None:
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
    ok = True
    rows = {}
    for p in CERT_J:
        F = _kernel_field(p)
        rs = mm_n1_set(F)
        wrong_hit = False
        rec_vals = set()
        for k in even_ks(p):
            got = J_mm_live(F, k, rs)
            named = J_mm_named(p, k)
            rec_vals.add(int(round(got)))
            if abs(got - named) > 1e-8:
                ok = False
            if abs(got - J_mm_wrong_half(p, k)) > 1e-6:
                wrong_hit = True
        if not wrong_hit:
            ok = False
        rows[p] = {
            "n": len(rs),
            "n_named": (p + 1) // 2,
            "values": sorted(rec_vals),
        }
        if len(rs) != (p + 1) // 2:
            ok = False
    # p≡1 empty
    if J_mm_named(5, 2) != 0:
        ok = False
    F5 = _kernel_field(5)
    if mm_n1_set(F5):
        ok = False
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """p=7 pairing identity. Fail: A−2p."""
    ok = True
    p, D = 7, D_live(7)
    F = _kernel_field(p)
    Q = live_Q_on_T(p)
    flip_hit = False
    max_live = -1e9
    for k in even_ks(p):
        lv = live_pair(p, k, F=F, Q=Q)
        nm = float(pair_p3_named(p, k, D))
        if abs(lv - nm) > 1e-8:
            ok = False
        if abs(lv - float(pair_p3_A2p(p, k, D))) > 1e-6:
            flip_hit = True
        max_live = max(max_live, lv)
    if not flip_hit:
        ok = False
    if abs(max_live - 200 / 409) > 1e-8:
        ok = False
    return {
        "proved": bool(ok),
        "max_live": max_live,
        "max_named": "200/409",
        "A_minus_3p": A_num(7) - 21,
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "J_{--} named for p≡3. A−3p only certified at p=7. "
            "D unnamed. phi_F not imported."
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
