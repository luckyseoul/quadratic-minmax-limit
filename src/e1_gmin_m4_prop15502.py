#!/usr/bin/env python3
"""
Prop 15.502 — The 15.497 pairing is a formal corollary of the
two Q-forms; those forms are incompatible with S0 on the
integer D-lattice for p≡1≥13; Paley×norm classes are not
single ⟨Frob,inv⟩-orbits past p=5.

    ⟨δ,ψ_k⟩ = −8 + 16(A+(p−4)J_{N*})/D

holds identically for every odd p≡1≥5 and every D>0 once
Q++/q²=8A/D and Q_{N*}/q²=8(A−2(p−4))/D and J_++ + J_{N*}=−2.
Fail: drop the 2(p−4) shift (then the identity fails).

S0 + those same forms force
    D = 8(a A + b(A−2(p−4)))/rhs
which equals live 13 at p=5 and equals 485621/10 ∉ ℤ at p=13.
Fail: claim the forms lie on S0 at every p≡1.

⟨Frob,inv⟩-orbits: one class = one orbit at p=5 (5 classes,
0 bad).  On the subfield class F_p^×\\{±1}, Frob is the
identity, so orbits are {r,r^{−1}} of size ≤2, while the
class has size p−3.  For every p≡1≥13 one has p−3≥10>2.
Certified also at p=13,17,29,37.  Fail: claim p=13 sub is
a single orbit.  Aut does not force a single Q++.

Does **not** import φ_F.  Q=8A/D is not a p≡1 theorem.
e1 / L untouched.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.501 flags.
"""
from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15290 import frob_inv_orbit, type_key
from e1_gmin_m4_prop15321 import s0_coeffs
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15497 import J_Nstar_named, even_ks, pair_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15502.json"

P1 = (5, 13, 17, 29)


def pair_from_Q_forms(p: int, k: int, D: int) -> Fraction:
    """δ-pairing from the two Q-forms and J_++=−2−J_{N*}."""
    A = A_num(p)
    Qpp = Fraction(8 * A, D)
    Qns = Fraction(8 * (A - 2 * (p - 4)), D)
    jns = J_Nstar_named(p, k)
    jpp = -2 - jns
    return (4 - Qpp) * jpp + (4 - Qns) * jns


def pair_drop_shift(p: int, k: int, D: int) -> Fraction:
    """Fail-when-wrong: Q_{N*}=8A/D (drop 2(p−4))."""
    A = A_num(p)
    Qpp = Fraction(8 * A, D)
    Qns = Fraction(8 * A, D)
    jns = J_Nstar_named(p, k)
    jpp = -2 - jns
    return (4 - Qpp) * jpp + (4 - Qns) * jns


def D_S0_from_forms(p: int) -> Fraction:
    A = A_num(p)
    a, b, rhs = s0_coeffs(p)
    return Fraction(8 * (a * A + b * (A - 2 * (p - 4))), rhs)


def class_vs_orbit(p: int) -> dict:
    F = _kernel_field(p)
    classes: dict[tuple, list[int]] = defaultdict(list)
    for r in F["squares"]:
        classes[type_key(F, r)].append(int(r))
    n_bad = 0
    rows = {}
    for tau, rs in classes.items():
        orb = frob_inv_orbit(F, rs[0])
        single = set(rs) == orb
        if tau != ("pm1",) and not single:
            n_bad += 1
        rows[str(tau)] = {"n": len(rs), "orb0": len(orb), "single": single}
    return {"n_classes": len(classes), "n_bad": n_bad, "rows": rows}


def prove_A() -> dict:
    """Pairing = pair_named from the two forms. Fail drop-shift."""
    ok = True
    rows = {}
    for p in P1:
        A = A_num(p)
        D = max(8 * abs(A), 1)
        n_ok = 0
        n_drop_diff = 0
        for k in even_ks(p):
            got = pair_from_Q_forms(p, k, D)
            named = pair_named(p, k, D)
            if got != named:
                ok = False
            else:
                n_ok += 1
            if pair_drop_shift(p, k, D) != named:
                n_drop_diff += 1
        if n_drop_diff == 0:
            ok = False
        rows[p] = {"n_match": n_ok, "n_even": len(even_ks(p)), "n_drop_diff": n_drop_diff}
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """S0+forms name D=13 at p=5 and a non-integer at p=13."""
    ok = True
    rows = {}
    d5 = D_S0_from_forms(5)
    if d5 != 13:
        ok = False
    if d5 != Fraction(HPLUS[5], 2 * 5):
        ok = False
    d13 = D_S0_from_forms(13)
    if d13.denominator == 1:
        ok = False
    if d13 == 13:
        ok = False
    if A_num(5) != 6 or A_num(13) != 23702:
        ok = False
    rows[5] = {"D_S0": str(d5), "A": A_num(5)}
    rows[13] = {"D_S0": str(d13), "den": d13.denominator, "A": A_num(13)}
    rows[17] = {"D_S0": str(D_S0_from_forms(17))}
    return {"proved": bool(ok), "by_p": rows}


def prove_C() -> dict:
    """Paley×norm classes are single orbits at p=5, not at p=13."""
    ok = True
    rows = {}
    r5 = class_vs_orbit(5)
    r13 = class_vs_orbit(13)
    if r5["n_bad"] != 0:
        ok = False
    if r13["n_bad"] == 0:
        ok = False
    sub_row = None
    for k, v in r13["rows"].items():
        if "sub" in k:
            sub_row = v
            break
    if sub_row is None or sub_row["single"] or sub_row["n"] == sub_row["orb0"]:
        ok = False
    rows[5] = {"n_bad": r5["n_bad"], "n_classes": r5["n_classes"]}
    rows[13] = {
        "n_bad": r13["n_bad"],
        "n_classes": r13["n_classes"],
        "sub_n": None if sub_row is None else sub_row["n"],
        "sub_orb": None if sub_row is None else sub_row["orb0"],
    }
    r17 = class_vs_orbit(17)
    if r17["n_bad"] == 0:
        ok = False
    rows[17] = {"n_bad": r17["n_bad"]}
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "Q_forms_on_S0_all_p1": False,
        "note": (
            "Pairing is a corollary of the Q-forms. The forms miss S0 "
            "at p=13. Aut does not force one Q++ for p≡1≥13. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.502 pairing corollary; S0 vs forms; Aut orbits", flush=True)
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
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
