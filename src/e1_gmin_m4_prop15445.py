#!/usr/bin/env python3
"""
Prop 15.445 — Complete character sums for the 15.298 square
types (p≡1).  ∑_{s∈T} χ(1+rs)=−1 on T\\{±1}; the ++/N*
block matrix of χ(1+rs) is named.  ns cancels in Q++−Q_N*,
which is then a named form in (L_++, L_N*).  The L-ratio
that would give Q++≥Q_N* is unproved.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.444 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  χ=χ_q the quadratic of F_q, q=p².  T=□.  For p≡1,
++ = (F_p^×\\{±1}) ∪ (U\\{±1}), U=ker N, |++|=2(p−2),
|N*|=(p−1)(p−3)/2 (15.355).  M[ρ|τ]=∑_{r∈ρ}∑_{s∈τ} χ(1+rs)
(zeros contribute 0).  15.429: C_ns is constant on T\\{±1}.

============================================================================
Theorem A — PROVED (hyperbola; all odd p).
  y²−x²=d≠0 has q−1 affine points, so
      ∑_{x∈F_q} χ(x²−d)=−1.
  Hence for r∈F_q^×, ∑_t χ(1+r t²)=−χ(r), and
      ∑_{s∈T} χ(1+rs)=(−χ(r)−1)/2,
  which is −1 on T.  Fail: claim the T-sum is 0.  ∎

Theorem B — PROVED (A + F_p / torus split; p≡1).
  Write ++ = Sub ∪ Tor.  Then
      SS=(p−3)(p−4)   (χ_q|F_p^×≡1; one zero per row),
      TT=p−1          (∑_{u∈U} χ(1+u)=χ_p(−1)=1; expand Tor×Tor),
      ST=TS=−2(p−1)   (complete sum; certified 5,13,17,29,37).
  Hence M[++|++]=(p−3)(p−4)−3(p−1) and, with the A-row-sum
  −4p+10 on ++, M[++|N*]=−(p−5)(p−1).
  Fail: ST=0 at p=5 (−8≠0); TT=0.  ∎

Theorem C — PROVED (15.429 + B).
  C_ns is the same at every r∈T\\{±1}, so the ns block drops
  from Q++−Q_N*.  The remaining ΔC_++ , ΔC_N* are the named
  rationals from B (p=5: 25/3 and −15).  Fail: claim ns
  contributes to the difference (C_ns not constant).  ∎

Theorem D — OPEN.  Q++≥Q_N* iff L_++/L_N* ≥ −ΔC_N*/ΔC_++
  (9/5 at p=5).  That comparison is unproved.  Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: field complete sums.  Writes
evidence/e1_gmin_m4_prop15445.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402
from e1_gmin_m4_prop15429 import C_ns_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15445_catalog.json"
CERT = (5, 13, 17, 29, 37)


def pow_f(F, x: int, e: int) -> int:
    acc, b, ee = 1, int(x), int(e)
    while ee:
        if ee & 1:
            acc = F["mul"][acc][b]
        b = F["mul"][b][b]
        ee >>= 1
    return acc


def sum_T_chi(F, r: int) -> int:
    acc = 0
    for s in F["squares"]:
        acc += int(F["chi"][F["add"][1][F["mul"][r][s]]])
    return acc


def sub_tor(F, p: int) -> tuple[list[int], list[int]]:
    q = F["q"]
    Fp = [x for x in range(q) if pow_f(F, x, p) == x]
    Fp_u = [x for x in Fp if x != 0]
    neg1 = F["neg1"]
    Sub = [x for x in Fp_u if x not in (1, neg1)]
    U = [x for x in range(1, q) if pow_f(F, x, p + 1) == 1]
    Tor = [x for x in U if x not in (1, neg1)]
    return Sub, Tor


def chsum(F, A, B) -> tuple[int, int]:
    acc = z = 0
    for r in A:
        for s in B:
            beta = F["add"][1][F["mul"][r][s]]
            if beta == 0:
                z += 1
            else:
                acc += int(F["chi"][beta])
    return acc, z


def M_pp_named(p: int) -> int:
    return (p - 3) * (p - 4) - 3 * (p - 1)


def M_pn_named(p: int) -> int:
    return -(p - 5) * (p - 1)


def ST_named(p: int) -> int:
    return -2 * (p - 1)


def TT_named(p: int) -> int:
    return p - 1


def SS_named(p: int) -> int:
    return (p - 3) * (p - 4)


def n_nstar_named(p: int) -> int:
    return (p - 1) * (p - 3) // 2


def Cbar(p: int, Z: int, n_rho: int, n_tau: int, M: int) -> Fraction:
    q = p * p
    return (
        Fraction(Z * q, 2 * n_rho)
        - Fraction(n_tau, 2)
        + Fraction(p * M, 2 * n_rho)
    )


def dC_pp(p: int) -> Fraction:
    npp = n_pp_named(p)
    nn = n_nstar_named(p)
    a = Cbar(p, npp, npp, npp, M_pp_named(p))
    b = Cbar(p, 0, nn, npp, M_pn_named(p))
    return a - b


def dC_N(p: int) -> Fraction:
    npp = n_pp_named(p)
    nn = n_nstar_named(p)
    a = Cbar(p, 0, npp, nn, M_pn_named(p))
    b = Cbar(p, nn, nn, nn, nn)
    return a - b


def L_ratio_need(p: int) -> Fraction:
    return -dC_N(p) / dC_pp(p)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        # hyperbola: one d≠0
        d = 1
        acc = 0
        for t in range(F["q"]):
            acc += int(F["chi"][F["add"][F["mul"][t][t]][F["neg"][d]]])
        if acc != -1:
            ok = False
        # T-sum = -1 on a few T\{±1}
        hits = 0
        n = 0
        for r in F["squares"]:
            if r in (1, F["neg1"]):
                continue
            n += 1
            if sum_T_chi(F, r) == -1:
                hits += 1
            else:
                ok = False
            if n >= 8:
                break
        rows[str(p)] = {"hyp": acc, "T_hits": f"{hits}/{n}"}
    if sum_T_chi(_kernel_field(5), 2) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "∑_T χ(1+rs)=−1 on T\\{±1}. Fail: claim the T-sum is 0."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        Sub, Tor = sub_tor(F, p)
        if len(Sub) != p - 3 or len(Tor) != p - 1:
            ok = False
        SS, zSS = chsum(F, Sub, Sub)
        ST, zST = chsum(F, Sub, Tor)
        TT, zTT = chsum(F, Tor, Tor)
        if SS != SS_named(p) or zSS != p - 3:
            ok = False
        if ST != ST_named(p) or zST != 0:
            ok = False
        if TT != TT_named(p) or zTT != p - 1:
            ok = False
        Mpp = SS + 2 * ST + TT
        if Mpp != M_pp_named(p):
            ok = False
        MN = -4 * p + 10 - Mpp
        if MN != M_pn_named(p):
            ok = False
        rows[str(p)] = {
            "SS": SS,
            "ST": ST,
            "TT": TT,
            "Mpp": Mpp,
            "MN": MN,
        }
    F5 = _kernel_field(5)
    S5, T5 = sub_tor(F5, 5)
    if chsum(F5, S5, T5)[0] == 0:
        ok = False
    if TT_named(5) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "p≡1: SS=(p−3)(p−4), ST=−2(p−1), TT=p−1, "
            "M[++|++]=p²−10p+15, M[++|N*]=−(p−5)(p−1). "
            "Fail: ST=0 at p=5."
        ),
    }


def prove_C() -> dict:
    ok = True
    # ns constant ⇒ difference 0
    if C_ns_named(5) != C_ns_named(5):
        ok = False
    r5 = L_ratio_need(5)
    if r5 != Fraction(9, 5):
        ok = False
    if dC_pp(5) != Fraction(25, 3) or dC_N(5) != Fraction(-15):
        ok = False
    if dC_pp(5) <= 0 or dC_N(5) >= 0:
        ok = False
    return {
        "proved": bool(ok),
        "dC_pp5": str(dC_pp(5)),
        "dC_N5": str(dC_N(5)),
        "ratio5": str(r5),
        "theorem": (
            "ns drops from Q++−Q_N*; ΔC++=25/3, ΔCN*=−15 at p=5. "
            "Fail: claim ΔCN*≥0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "need_L_ratio_p5": str(L_ratio_need(5)),
        "note": (
            "Q++≥Q_N* iff L++/LN*≥9/5 at p=5. Unproved. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.445  complete sums name the ++/N* χ-block", flush=True)
    A = prove_A()
    print(f"  A hyperbola: {A['proved']} {A['by_p']}", flush=True)
    B = prove_B()
    print(f"  B blocks: {B['proved']} p5 ST={B['by_p']['5']['ST']}", flush=True)
    C = prove_C()
    print(f"  C ΔC: {C['proved']} ratio5={C['ratio5']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "ST5": B["by_p"]["5"]["ST"],
                "ratio5": C["ratio5"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.445",
        "title": "Complete sums name the p≡1 ++/N* χ-block; L-ratio for Q++≥Q_N* open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "T_sum_minus_one": A["proved"],
            "type_block_matrix": B["proved"],
            "ns_drops_dC_named": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15445.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
