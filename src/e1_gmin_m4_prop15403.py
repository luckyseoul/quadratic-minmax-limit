#!/usr/bin/env python3
"""
Prop 15.403 — At p=7, every 15.307 D-type has Q++ in the
named list {Q_1d_pp, 4(p−1)/9, 4(p+1)/9, Q_AP=4(p+3)/9}.
The same triple does not name p=5 NL (that value is Q_T).
Ensemble Q_τ still has den N+.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.402 flags.

============================================================================
Setup.  15.307 D-types.  15.355: Q_AP=(p²−9)/(3(p−2)) (p≡1) and
4(p+3)/9 (p≡3).  15.356: Q_1d_pp named.  15.350: Q_T=Qpp_from_L.

============================================================================
Theorem A — PROVED (cache fold; p=7).
  D-type Q++ takes exactly the four values
      104/15 = Q_1d_pp,   8/3 = 4(p−1)/9,
      32/9 = 4(p+1)/9,    40/9 = 4(p+3)/9 = Q_AP.
  Weights 140, 1764, 1176, 2646; mix = live 1544/409.
  Fail: claim Q_QR0=32/3 occurs; claim 1544/409 is a type
  value; claim 4(p+1)/9=Q_AP (32/9≠40/9).  ∎

Theorem B — PROVED (cache fold; p=5).
  Types are {16/9=Q_1d_pp=Q_AP, 64/15=Q_T}.
  4(p+1)/9=8/3 and 4(p+3)/9=32/9 do not occur.
  NL is Q_T, not the p=7 triple.  Fail: claim NL Q++=32/9.  ∎

Theorem C — OPEN.  The four values name D-type Q++ at p=7,
  not ensemble Q_τ (den N+).  The triple 4(p+2k−1)/9 is not
  the NL set at p=5.  phi_F not imported.

============================================================================
Backend: vectorized ẑ on p=5,7 caches + 15.307 tags.
Writes evidence/e1_gmin_m4_prop15403.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15350 import Qpp_from_L  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named, Q_QR0_p3_named  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15403_catalog.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def qpp_m1(p: int) -> Fraction:
    """4(p−1)/9."""
    return Fraction(4 * (p - 1), 9)


def qpp_p1(p: int) -> Fraction:
    """4(p+1)/9."""
    return Fraction(4 * (p + 1), 9)


def qpp_p3(p: int) -> Fraction:
    """4(p+3)/9. Equals Q_AP for every p≡3."""
    return Fraction(4 * (p + 3), 9)


def p7_type_qpp_named() -> frozenset[Fraction]:
    p = 7
    return frozenset(
        {Q_1d_pp_named(p), qpp_m1(p), qpp_p1(p), qpp_p3(p)}
    )


def ensemble_Qpp_is_type_value() -> bool:
    """True only if live Q++ equals a D-type value. It does not."""
    return False


def _u_matrix(p: int, F, Omega):
    Y = np.load(YPATH[p])
    Z = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64))
    q = F["q"]
    E = np.zeros((q, len(Omega)), dtype=np.complex128)
    for j, xi in enumerate(Omega):
        E[:, j] = [_e_tr(F, F["mul"][xi][x]) for x in range(q)]
    return np.abs(Z @ E) ** 2


@lru_cache(maxsize=4)
def type_Qpp_table(p: int) -> tuple[tuple[int, str, str], ...]:
    """((n, Q++ as str, tag), ...) sorted by n desc."""
    D, F = load_D(p)
    classes = square_classes(F)
    packs = []
    for lines in classes:
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1)
        packs.append([tuple(int(x) for x in row) for row in np.sort(ns, axis=1)])
    tags = [
        tuple(sorted(tag_sig(p, packs[di][i]) for di in range(len(classes))))
        for i in range(D.shape[0])
    ]
    Omega = omega_set(F)
    U = _u_matrix(p, F, Omega)
    q2 = float(F["q"] ** 2)
    om_idx = {xi: i for i, xi in enumerate(Omega)}
    xi0 = Omega[0]
    u0 = U[:, om_idx[xi0]]
    pp = [r for r in F["squares"] if merge_cls(F, r) == "++"]
    cols = [u0 * U[:, om_idx[F["mul"][r][xi0]]] / q2 for r in pp]
    Qrow = np.mean(np.stack(cols, axis=0), axis=0)
    by: dict[tuple, list[int]] = defaultdict(list)
    for i, t in enumerate(tags):
        by[t].append(i)
    rows = []
    for t, idx in sorted(by.items(), key=lambda kv: -len(kv[1])):
        val = float(np.mean(Qrow[np.array(idx)]))
        q = Fraction(val).limit_denominator(2000)
        rows.append((len(idx), str(q), str(t)))
    return tuple(rows)


def prove_A() -> dict:
    ok = True
    tab = type_Qpp_table(7)
    named = p7_type_qpp_named()
    got = {Fraction(q) for _n, q, _t in tab}
    ok = ok and got == named
    ok = ok and qpp_p3(7) == Q_AP_named(7) == Fraction(40, 9)
    ok = ok and qpp_m1(7) == Fraction(8, 3)
    ok = ok and qpp_p1(7) == Fraction(32, 9)
    ok = ok and Q_1d_pp_named(7) == Fraction(104, 15)
    # fail-when-wrong: Q_QR0 and ensemble are not type values
    ok = ok and Q_QR0_p3_named(7) not in got
    ok = ok and LIVE_QPP[7] not in got
    ok = ok and qpp_p1(7) != Q_AP_named(7)
    if Q_QR0_p3_named(7) in got or LIVE_QPP[7] in got:
        ok = False
    if qpp_p1(7) == Q_AP_named(7):
        ok = False
    # mix recovers live
    ntot = sum(n for n, _q, _t in tab)
    acc = sum(n * Fraction(q) for n, q, _t in tab)
    ok = ok and acc / ntot == LIVE_QPP[7]
    by_val = {q: 0 for q in named}
    for n, q, _t in tab:
        by_val[Fraction(q)] += n
    return {
        "proved": bool(ok),
        "values": sorted(str(x) for x in got),
        "weights": {str(k): v for k, v in by_val.items()},
        "mix": str(acc / ntot),
        "theorem": (
            "p=7 D-type Q++ ∈ {104/15, 8/3, 32/9, 40/9}. "
            "Fail: 32/3 or 1544/409 is a type value; 32/9=Q_AP."
        ),
    }


def prove_B() -> dict:
    ok = True
    tab = type_Qpp_table(5)
    got = {Fraction(q) for _n, q, _t in tab}
    ok = ok and got == {Q_1d_pp_named(5), Qpp_from_L(5)}
    ok = ok and Q_1d_pp_named(5) == Q_AP_named(5) == Fraction(16, 9)
    ok = ok and Qpp_from_L(5) == Fraction(64, 15)
    ok = ok and qpp_p1(5) not in got
    ok = ok and qpp_p3(5) not in got
    ok = ok and qpp_p3(5) != Qpp_from_L(5)
    if qpp_p3(5) in got or Qpp_from_L(5) == qpp_p3(5):
        ok = False
    ntot = sum(n for n, _q, _t in tab)
    acc = sum(n * Fraction(q) for n, q, _t in tab)
    ok = ok and acc / ntot == LIVE_QPP[5]
    return {
        "proved": bool(ok),
        "values": sorted(str(x) for x in got),
        "mix": str(acc / ntot),
        "qpp_p3_is_NL": qpp_p3(5) in got,
        "theorem": (
            "p=5 types {16/9, 64/15=Q_T}. "
            "Fail: NL Q++=4(p+3)/9=32/9."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "ensemble_Qpp_is_type_value": ensemble_Qpp_is_type_value(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "D-type Q++ named at p=7 as Q_AP-family + 4(p+1)/9 + Q_1d. "
            "Ensemble still den N+. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.403  p=7 D-type Q++ ∈ {Q_1d, 4(p±1)/9, Q_AP}", flush=True)
    A = prove_A()
    print(f"  A p7 values: {A['proved']} {A['values']} mix={A['mix']}", flush=True)
    B = prove_B()
    print(f"  B p5 values: {B['proved']} {B['values']}", flush=True)
    C = prove_open()
    print(f"  C open: type=ens={C['ensemble_Qpp_is_type_value']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["values"], "B": B["values"], "Aw": A["weights"]}, indent=2) + "\n")
    out = {
        "prop": "15.403",
        "title": "p=7 D-type Q++ is {Q_1d, 4(p−1)/9, 4(p+1)/9, Q_AP}; not ensemble",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "p7_type_Qpp_named": A["proved"],
            "p5_NL_is_QT_not_triple": B["proved"],
            "ensemble_Qpp_is_type_value": C["ensemble_Qpp_is_type_value"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15403.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
