#!/usr/bin/env python3
"""
Prop 15.449 — L_NL on the 15.448 χ(s±1) slots is a function
of (slot, N_{q/p}(s)).  The constant square/unmixed-ns values
have essential den c(7)=19, which is not a Gauss/Jacobi
integer.  ystar is all of NL at p=5 and is not L_NL at p=7.
L_NL / Q_τ stay unnamed in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.448 flags.

============================================================================
Setup.  15.448 names every Type+ L_1D slot.  Ensemble
L=w L_1D+(1−w)L_NL, w=n_1d/|H+|.  NL = H+ \\ Type+.  N(s)=s^{p+1}.
15.430 L_par(R) and 15.138 ystar are Max+-free NL pieces.
15.439: 19 is essential in Q_NL(7) and not GJ.

============================================================================
Theorem A — PROVED (NL cache; certified p=5,7).
  L_NL is constant on (15.448 slot, N_{q/p}(s)).  The 15.448
  slot alone splits: ns_mix takes {12,13} at p=5 and three
  values at p=7; sq_mm splits at p=7.  Fail: claim L_NL
  constant on ns_mix.  ∎

Theorem B — PROVED (NL cache + 15.430).
  On slots with a single N-class,
      p=5: Sub=sq_mm=q−1, pm1=sq_mix=q+1,
           ns_11=(q+1)/2, ns_mm=(q−1)/2,
      p=7: Sub=2113/19, pm1=4435/38, sq_11=673/6,
           sq_mix=6269/57, ns_11=2853/38, ns_mm=2733/38.
  Fail: L_Sub^{NL}(7)=L_par(R)(7)=133.  ∎

Theorem C — PROVED (reduced Fraction + 15.330 S7).
  gcd(2113,19)=1, so 19 is essential in L_Sub^{NL}(7).
  19∉named_gj_ints(7).  Same 19 in 4435/38, 2853/38, 2733/38,
  6269/57.  Fail: claim 19|2113; claim 19∈S7.  ∎

Theorem D — PROVED (15.138 orbit; certified p=5,7).
  The ystar Aut_∞ orbit is all of NL at p=5 (size 100) and
  has L_Sub=1321/12 ≠ 2113/19 at p=7.  Fail: L_ystar=L_NL
  at p=7.  ∎

Theorem E — OPEN.  L_NL is not a field GJ formula in p
  (essential 19).  Ensemble L_τ still has w=n_1d/|H+|.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: p=5,7 NL fold + ystar Aut orbit.  Writes
evidence/e1_gmin_m4_prop15449.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15308 import build_G  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS, named_gj_ints  # noqa: E402
from e1_gmin_m4_prop15430 import Lpar_R_named  # noqa: E402
from e1_gmin_m4_prop15445 import pow_f  # noqa: E402
from e1_gmin_m4_prop15448 import slot_of  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15449_catalog.json"


def L_from_N(N: np.ndarray, F) -> dict[int, Fraction]:
    ENN = (N.astype(np.float64).T @ N) / N.shape[0]
    out: dict[int, Fraction] = {}
    for s in range(1, F["q"]):
        acc = 0.0
        n = 0
        for d in range(1, F["q"]):
            if F["chi"][d] != 1:
                continue
            acc += ENN[d, F["mul"][s][d]]
            n += 1
        out[s] = Fraction(acc / n).limit_denominator(10**7)
    return out


def type_plus_mask(p: int, F, Yplus: np.ndarray) -> np.ndarray:
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    plus = list(combinations(range(p), m))
    lifts = {
        tuple(np.sign(lift_sigma_vector(p, a, b, set(A))).astype(np.int8))
        for a, b in forms
        for A in plus
    }
    return np.array([tuple(row.astype(np.int8)) in lifts for row in Yplus])


def nl_L(p: int):
    from e1_gmin_m4_prop15298 import YPATH

    F = _kernel_field(p)
    Y = np.load(YPATH[p])
    Yp = np.sign(Y[Y[:, 0] > 0].astype(np.float64))
    N = load_N(p, F)
    mask = type_plus_mask(p, F, Yp)
    if int(mask.sum()) != n_1d(p):
        raise RuntimeError(f"1D mask {int(mask.sum())} != {n_1d(p)}")
    return F, L_from_N(N[~mask], F), int((~mask).sum())


def bucket_slot_N(F, p: int, L: dict[int, Fraction]):
    by: dict[tuple, set[Fraction]] = defaultdict(set)
    by_slot: dict[str, set[Fraction]] = defaultdict(set)
    for s, v in L.items():
        sl = slot_of(F, p, s)
        nu = int(pow_f(F, s, p + 1))
        by[(sl, nu)].add(v)
        by_slot[sl].add(v)
    return by, by_slot


def orbit_ys(p: int, y0) -> np.ndarray:
    F = _kernel_field(p)
    G = build_G(F, use_frob=True)
    y0 = np.asarray(y0, dtype=np.float64)
    seen: dict[bytes, np.ndarray] = {}
    for perm in G:
        ynew = np.empty_like(y0)
        ynew[perm] = y0
        if ynew[0] < 0:
            ynew = -ynew
        seen[ynew.tobytes()] = ynew
    return np.stack(list(seen.values()))


def N_from_Y(Y: np.ndarray, F) -> np.ndarray:
    Ds = np.sign(Y)[:, 1 : 1 + F["q"]] < 0
    add = np.array(F["add"], dtype=np.int32)
    return np.stack(
        [(Ds & Ds[:, add[:, d]]).sum(axis=1) for d in range(F["q"])], axis=1
    )


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, L, nnl = nl_L(p)
        by, by_slot = bucket_slot_N(F, p, L)
        split_slot = {sl: sorted(vs) for sl, vs in by_slot.items() if len(vs) > 1}
        split_sn = {str(k): sorted(vs) for k, vs in by.items() if len(vs) > 1}
        if split_sn:
            ok = False
        if "ns_mix" not in split_slot:
            ok = False
        if p == 5 and split_slot.get("ns_mix") != [12, 13]:
            ok = False
        if p == 7 and "sq_mm" not in split_slot:
            ok = False
        rows[str(p)] = {
            "n_NL": nnl,
            "split_slots": {sl: [str(v) for v in vs] for sl, vs in split_slot.items()},
            "n_classes": len(by),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "L_NL constant on (15.448 slot, N). Fail: constant on ns_mix."
        ),
    }


def prove_B() -> dict:
    ok = True
    want5 = {
        "Sub": Fraction(24),
        "sq_mm": Fraction(24),
        "pm1": Fraction(26),
        "sq_mix": Fraction(26),
        "ns_11": Fraction(13),
        "ns_mm": Fraction(12),
    }
    want7 = {
        "Sub": Fraction(2113, 19),
        "pm1": Fraction(4435, 38),
        "sq_11": Fraction(673, 6),
        "sq_mix": Fraction(6269, 57),
        "ns_11": Fraction(2853, 38),
        "ns_mm": Fraction(2733, 38),
    }
    live = {}
    for p, want in ((5, want5), (7, want7)):
        F, L, _nnl = nl_L(p)
        _by, by_slot = bucket_slot_N(F, p, L)
        row = {}
        for sl, val in want.items():
            vs = by_slot.get(sl, set())
            row[sl] = [str(v) for v in sorted(vs)]
            if vs != {val}:
                ok = False
        live[str(p)] = row
        if p == 5:
            q = p * p
            if want["Sub"] != q - 1 or want["pm1"] != q + 1:
                ok = False
    if Lpar_R_named(7) == want7["Sub"]:
        ok = False
    if Lpar_R_named(7) != 133:
        ok = False
    return {
        "proved": bool(ok),
        "live": live,
        "Lpar_R_7": str(Lpar_R_named(7)),
        "theorem": (
            "Constant-slot L_NL as listed. Fail: L_Sub^{NL}(7)=L_par(R)=133."
        ),
    }


def prove_C() -> dict:
    ok = True
    if gcd(2113, 19) != 1:
        ok = False
    if 19 in named_gj_ints(7):
        ok = False
    if 2113 % 19 == 0:
        ok = False
    s7 = named_gj_ints(7)
    for bad in (19, 38, 57, 171):
        if bad in s7 or -bad in s7:
            ok = False
    return {
        "proved": bool(ok),
        "gcd": gcd(2113, 19),
        "in_S7": 19 in named_gj_ints(7),
        "theorem": (
            "19 essential in 2113/19 and not in named_gj_ints(7). "
            "Fail: 19|2113; fail: 19∈S7."
        ),
    }


def prove_D() -> dict:
    ok = True
    # p=5: ystar orbit = NL
    F5 = _kernel_field(5)
    Y5 = orbit_ys(5, construct_ystar(5)["y"])
    if len(Y5) != HPLUS[5] - n_1d(5):
        ok = False
    L5 = L_from_N(N_from_Y(Y5, F5), F5)
    _by5, sl5 = bucket_slot_N(F5, 5, L5)
    if sl5.get("Sub") != {24}:
        ok = False
    # p=7: ystar Sub ≠ NL Sub
    F7 = _kernel_field(7)
    Y7 = orbit_ys(7, construct_ystar(7)["y"])
    L7 = L_from_N(N_from_Y(Y7, F7), F7)
    _by7, sl7 = bucket_slot_N(F7, 7, L7)
    ysub = next(iter(sl7["Sub"])) if sl7.get("Sub") else None
    if ysub != Fraction(1321, 12):
        ok = False
    if ysub == Fraction(2113, 19):
        ok = False
    return {
        "proved": bool(ok),
        "ystar5": len(Y5),
        "ystar7": len(Y7),
        "L_Sub_ystar7": str(ysub),
        "L_Sub_NL7": "2113/19",
        "theorem": (
            "ystar=NL at p=5; L_Sub(ystar)=1321/12≠2113/19 at p=7. "
            "Fail: ystar names L_NL."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "L_NL on 15.448 slots is not a field GJ formula (essential 19). "
            "Ensemble L_τ still has w=n_1d/|H+|. Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.449  L_NL on 15.448 slots; 19 essential, not GJ", flush=True)
    A = prove_A()
    print(f"  A constancy: {A['proved']} splits5={A['by_p']['5']['split_slots']}", flush=True)
    B = prove_B()
    print(f"  B live: {B['proved']} LparR7={B['Lpar_R_7']}", flush=True)
    C = prove_C()
    print(f"  C den19: {C['proved']} gcd={C['gcd']} inS7={C['in_S7']}", flush=True)
    D = prove_D()
    print(
        f"  D ystar: {D['proved']} p7={D['L_Sub_ystar7']} vs {D['L_Sub_NL7']}",
        flush=True,
    )
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']} e1={E['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "D": D["proved"],
                "ystar7": D["L_Sub_ystar7"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.449",
        "title": (
            "L_NL is a function of (15.448 slot, N); essential den 19 "
            "is not GJ; ystar is not L_NL at p=7"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "const_on_slot_N": A["proved"],
            "live_constant_slots": B["proved"],
            "den19_not_GJ": C["proved"],
            "ystar_not_NL": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15449.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
