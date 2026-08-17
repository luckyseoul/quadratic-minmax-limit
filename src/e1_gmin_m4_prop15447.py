#!/usr/bin/env python3
"""
Prop 15.447 — The p≡1 Sub/Tor/N* χ-block is named by complete
sums, and Q_τ is the half-Gauss image of the split L through
that block.  Dropping ST misses live Q++.  L_τ still unnamed
in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.446 flags.

============================================================================
Setup.  p≡1.  Sub=F_p^×\\{±1}, Tor=U\\{±1}, N*=T\\{±1}\\(++).
M[ρ|τ]=∑_{r∈ρ}∑_{s∈τ} χ(1+rs).  15.445: hyperbola ∑_T χ(1+rs)=−1;
SS=(p−3)(p−4), ST=−2(p−1), TT=p−1.  15.298: Q=16(A0+∑ L S).
A0=k(k+p)=p²(q−1)/4.

============================================================================
Theorem A — PROVED (15.445 + row sums; p≡1).
  For r∈Tor, ∑_{T\\{±1}} χ=−1 and the Sub/Tor rows give
      ∑_{s∈N*} χ(1+rs)=0, hence TN*=NT*=0.
  For r∈Sub, ∑_{T\\{±1}} χ=−3 and SS, ST give
      SN*=NS*=−(p−5)(p−1).
  Fail: SN*=0 at p=13 (−96≠0).  ∎

Theorem B — PROVED (A + 15.445 −1-sum on N* rows).
  ∑_{r∈N*} (χ(r+1)+χ(r−1))=−2(p−1) (certified 5,13,17), hence
      NN*=|N*|=(p−1)(p−3)/2.
  Fail: NN*=0 at p=5.  ∎

Theorem C — PROVED (15.298 + A,B C̄; certified p=5).
  C̄_{χ,τ}(ρ) from (Z, M) recovers live Q++=48/13 and Q_N*=32/13
  when the split L is inserted.  Fail: set ST=0 (miss 48/13).  ∎

Theorem D — OPEN.  L_τ is not a formula in p, so Q_τ is named
  only as this GJ-linear image of L.  phi_F not imported.

============================================================================
Backend: field complete sums + p=5 Max+ L.  Writes
evidence/e1_gmin_m4_prop15447.json
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
from e1_gmin_m4_prop15429 import A0_named  # noqa: E402
from e1_gmin_m4_prop15445 import (  # noqa: E402
    M_pn_named,
    SS_named,
    ST_named,
    TT_named,
    chsum,
    n_nstar_named,
    n_pp_named,
    pow_f,
    sub_tor,
)

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15447_catalog.json"
CERT = (5, 13, 17)


def SN_named(p: int) -> int:
    return -(p - 5) * (p - 1)


def NN_named(p: int) -> int:
    return n_nstar_named(p)


def parts(F, p: int) -> dict[str, list[int]]:
    Sub, Tor = sub_tor(F, p)
    neg1 = F["neg1"]
    plus = set(Sub) | set(Tor) | {1, neg1}
    Nstar = [x for x in F["squares"] if x not in plus]
    return {"Sub": Sub, "Tor": Tor, "N*": Nstar}


def nstar_paley_sum(F, Nstar: list[int]) -> int:
    acc = 0
    for r in Nstar:
        acc += int(F["chi"][F["add"][r][1]]) + int(
            F["chi"][F["add"][r][F["neg"][1]]]
        )
    return acc


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        P = parts(F, p)
        ST, _ = chsum(F, P["Sub"], P["Tor"])
        TN, _ = chsum(F, P["Tor"], P["N*"])
        SN, _ = chsum(F, P["Sub"], P["N*"])
        if ST != ST_named(p):
            ok = False
        if TN != 0:
            ok = False
        if SN != SN_named(p):
            ok = False
        rows[str(p)] = {"ST": ST, "TN": TN, "SN": SN}
    if SN_named(13) == 0:
        ok = False
    F13 = _kernel_field(13)
    P13 = parts(F13, 13)
    if chsum(F13, P13["Sub"], P13["N*"])[0] == 0:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "TN*=0 and SN*=−(p−5)(p−1) for p≡1. Fail: SN*=0 at p=13."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        P = parts(F, p)
        pal = nstar_paley_sum(F, P["N*"])
        NN, z = chsum(F, P["N*"], P["N*"])
        want_pal = -2 * (p - 1)
        if pal != want_pal:
            ok = False
        if NN != NN_named(p) or z != NN_named(p):
            ok = False
        rows[str(p)] = {"paley_sum": pal, "NN": NN, "Z": z}
    F5 = _kernel_field(5)
    P5 = parts(F5, 5)
    if chsum(F5, P5["N*"], P5["N*"])[0] == 0:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "∑_{N*}(χ(r+1)+χ(r−1))=−2(p−1) ⇒ NN*=|N*|. Fail: NN*=0 at p=5."
        ),
    }


def Cbar_sq(p: int, chi: int, tau: str, rho: str, st: int | None = None) -> Fraction:
    sizes = {"Sub": p - 3, "Tor": p - 1, "N*": n_nstar_named(p)}
    st_val = ST_named(p) if st is None else st
    Mtab = {
        ("Sub", "Sub"): SS_named(p),
        ("Sub", "Tor"): st_val,
        ("Tor", "Sub"): st_val,
        ("Tor", "Tor"): TT_named(p),
        ("Sub", "N*"): SN_named(p),
        ("N*", "Sub"): SN_named(p),
        ("Tor", "N*"): 0,
        ("N*", "Tor"): 0,
        ("N*", "N*"): NN_named(p),
    }
    Ztab = {
        ("Sub", "Sub"): p - 3,
        ("Tor", "Tor"): p - 1,
        ("N*", "N*"): n_nstar_named(p),
    }
    if rho == "++":
        nS, nT = sizes["Sub"], sizes["Tor"]
        return (
            nS * Cbar_sq(p, chi, tau, "Sub", st=st)
            + nT * Cbar_sq(p, chi, tau, "Tor", st=st)
        ) / (nS + nT)
    n_rho = sizes[rho]
    n_tau = sizes[tau]
    m = Mtab.get((rho, tau), 0)
    z = Ztab.get((rho, tau), 0)
    q = p * p
    G = p
    # z terms of (q-1)/2 and (n_tau-z) terms of (-1 + chi*G*χ(β))/2
    # sum_χ(β) over nonzero pairs = M
    return (
        Fraction(z * q, 2 * n_rho)
        - Fraction(n_tau, 2)
        + Fraction(chi * G * m, 2 * n_rho)
    )


def recover_Q(p: int, drop_SN: bool = False, drop_ST: bool = False) -> dict[str, Fraction]:
    """Q/q² = 16(A0 + ∑ L C)/q². Named C̄ on Sub/Tor/N*; live C on pm1/ns."""
    F = _kernel_field(p)
    q = F["q"]
    Ltab = build_L(load_N(p, F), F)
    P = parts(F, p)
    xi = omega_set(F)[0]
    G = float(
        sum(F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)).real
    )

    def classify(s: int) -> str:
        if s in (1, F["neg1"]):
            return "pm1"
        if F["chi"][s] != 1:
            return "ns"
        if pow_f(F, s, p) == s:
            return "Sub"
        if pow_f(F, s, p + 1) == 1:
            return "Tor"
        return "N*"

    def L_of(tau: str, cd: int) -> float:
        for s in range(1, q):
            if classify(s) == tau:
                return float(Ltab[(cd, ns_key(F, s))])
        raise KeyError(tau)

    plus = P["Sub"] + P["Tor"]
    out: dict[str, Fraction] = {}
    for name, rho in (("++", plus), ("N*", P["N*"])):
        acc = float(A0_named(p))
        for tau in ("Sub", "Tor", "N*"):
            for cd in (1, -1):
                st = 0 if drop_ST else None
                if drop_SN:
                    cbar = float(_Cbar_drop_SN(p, cd, tau, name))
                elif name == "++":
                    cbar = float(Cbar_sq(p, cd, tau, "++", st=st))
                else:
                    cbar = float(Cbar_sq(p, cd, tau, "N*", st=st))
                acc += L_of(tau, cd) * cbar
        # pm1 + ns: live half-Gauss (L not constant on ns)
        for r in rho:
            for s in range(1, q):
                tau = classify(s)
                if tau not in ("pm1", "ns"):
                    continue
                for cd in (1, -1):
                    acc += Ltab[(cd, ns_key(F, s))] * half_gauss(
                        F, F["add"][1][F["mul"][r][s]], cd, G, q
                    ) / len(rho)
        out[name] = Fraction(16 * acc / (q * q)).limit_denominator(10**6)
    return out


def _Cbar_drop_SN(p: int, chi: int, tau: str, rho_name: str) -> Fraction:
    """C̄ with SN* forced to 0 (fail-when-wrong)."""
    sizes = {"Sub": p - 3, "Tor": p - 1, "N*": n_nstar_named(p), "++": n_pp_named(p)}
    Mtab = {
        ("Sub", "Sub"): SS_named(p),
        ("Sub", "Tor"): ST_named(p),
        ("Tor", "Sub"): ST_named(p),
        ("Tor", "Tor"): TT_named(p),
        ("Sub", "N*"): 0,  # dropped
        ("N*", "Sub"): 0,
        ("Tor", "N*"): 0,
        ("N*", "Tor"): 0,
        ("N*", "N*"): NN_named(p),
    }
    Ztab = {
        ("Sub", "Sub"): p - 3,
        ("Tor", "Tor"): p - 1,
        ("N*", "N*"): n_nstar_named(p),
    }

    def one(rho: str) -> Fraction:
        n_rho = sizes[rho]
        n_tau = sizes[tau]
        m = Mtab.get((rho, tau), 0)
        z = Ztab.get((rho, tau), 0)
        q = p * p
        return (
            Fraction(z * q, 2 * n_rho)
            - Fraction(n_tau, 2)
            + Fraction(chi * p * m, 2 * n_rho)
        )

    if rho_name == "++":
        nS, nT = sizes["Sub"], sizes["Tor"]
        return (nS * one("Sub") + nT * one("Tor")) / (nS + nT)
    return one("N*")


def prove_C() -> dict:
    ok = True
    rec = recover_Q(5, drop_ST=False)
    rec_bad = recover_Q(5, drop_ST=True)
    if rec["++"] != Fraction(48, 13):
        ok = False
    if rec["N*"] != Fraction(32, 13):
        ok = False
    if rec_bad["++"] == Fraction(48, 13):
        ok = False
    return {
        "proved": bool(ok),
        "Qpp": str(rec["++"]),
        "QN": str(rec["N*"]),
        "Qpp_dropST": str(rec_bad["++"]),
        "theorem": (
            "Named C̄ + split L recovers 48/13 and 32/13. "
            f"Fail: ST=0 gives {rec_bad['++']}≠48/13."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Q_τ is the GJ-linear image of split L. L_τ unnamed in p. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.447  Q_τ = half-Gauss of split L through named C̄", flush=True)
    A = prove_A()
    print(f"  A SN*/TN*: {A['proved']} p13 SN={A['by_p']['13']['SN']}", flush=True)
    B = prove_B()
    print(f"  B NN*: {B['proved']} p5 NN={B['by_p']['5']['NN']}", flush=True)
    C = prove_C()
    print(f"  C recover: {C['proved']} Q++={C['Qpp']} drop={C['Qpp_dropST']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "Qpp": C["Qpp"],
                "dropST": C["Qpp_dropST"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.447",
        "title": "Q_τ is the half-Gauss image of split L through the named Sub/Tor/N* χ-block",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "SN_TN_named": A["proved"],
            "NN_named": B["proved"],
            "Q_from_split_L": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15447.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
