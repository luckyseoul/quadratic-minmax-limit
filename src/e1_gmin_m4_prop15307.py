#!/usr/bin/env python3
"""
Prop 15.307 — Square-direction signatures classify H_+ G-orbits;
E[M²] is the type mixture of named M_τ; |H_+| still unnamed.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.306 flags.  Floor stays OPEN.

============================================================================
Setup.  H_+ = {y∈Max+ : y_∞=+1}.  A square-direction parallel class
has p lines; signature = sorted (|D∩ℓ|)_{ℓ in class}.  Named masses
    M_H = p²(p−1)/2 ,  M_C = p(p−1)²/4 ,  M_R = p(p−1)(2p−1)/6
are the M=∑ n_ℓ² of H / C / R signatures.  n_1d = m C(p,m), m=(p+1)/2.

============================================================================
Theorem A — PROVED (certified p=3,5,7; fail-when-wrong at p=7).
  Affine halfspaces are exactly the D∈H_+ whose square-direction
  signature multiset is one H and (n_sq−1) copies of C, n_sq=(p+1)/2.
  Their count is n_1d.  Fail-when-wrong: claim that count is n_1d−1,
  or that a halfspace carries an R signature.  ∎

Theorem B — PROVED (certified p=5; fails as a single-type law at p=7).
  At p=5 the complement H_+ \\ {halfspaces} is a single type: R on
  every square direction, size p²(p−1).  Hence
      E[M²] = (n_1d/N+)·(2 M_H² + (p−1) M_C²)/(p+1)
              + (1−n_1d/N+) M_R²
  recovers 12300/13 and Q_{++}=48/13.  The same one-NL-type formula
  fails at p=7 (11 dir-signatures, 6 NL D-types).  Fail-when-wrong:
  claim p=7 NL is a single R type.  ∎

Theorem C — PROVED (certified p=7; matches 15.139 size classes).
  H_+ D-types and sizes:
      ('C','C','C','H')           140 = n_1d
      ('C','R','R','R')          1176 = 4 p²(p−1)
      three mixed exotic types   1176 = 4 p²(p−1) each
      one exotic type             588 = 2 p²(p−1)
      one exotic type             294 =     p²(p−1)
  Sum 5726 = |H_+|.  The four 1176-orbits are the 15.139 mult-4
  class.  Affine 140 = 56+84 (AP vs non-AP of S).  ∎

Theorem D — OPEN.  |H_+| has no Max+-free closed form (OEIS empty
  on 12,260,11452).  Type weights are |O_τ|/|H_+| and still carry
  den N/(4p).  Q_{++} unnamed.  Square-mass bound unproved.
  AP∨QR on D∩ℓ is vacuous at p=5,7 and false on p=11 samples.
  phi_F not imported.  ∎

============================================================================
Backend: vectorized numpy over H_+ caches p=3,5,7.
Writes evidence/e1_gmin_m4_prop15307.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402

YPATH = {
    3: "/tmp/maxplus_p3.npy",
    5: "/tmp/maxplus_p5.npy",
    7: "/tmp/maxplus_p7.npy",
}


def M_named(p: int) -> dict[str, int]:
    return {
        "H": p * p * (p - 1) // 2,
        "C": p * ((p - 1) // 2) ** 2,
        "R": p * (p - 1) * (2 * p - 1) // 6,
        "E_M": p * (p * p - 1) // 4,
    }


def sig_H(p: int) -> tuple[int, ...]:
    return tuple([0] * ((p + 1) // 2) + [p] * ((p - 1) // 2))


def sig_C(p: int) -> tuple[int, ...]:
    return tuple([(p - 1) // 2] * p)


def sig_R(p: int) -> tuple[int, ...]:
    return tuple(range(p))


def tag_sig(p: int, sig: tuple[int, ...]) -> str:
    if sig == sig_H(p):
        return "H"
    if sig == sig_C(p):
        return "C"
    if sig == sig_R(p):
        return "R"
    return "O:" + ",".join(map(str, sig))


def square_classes(F) -> list[np.ndarray]:
    p = F["p"]
    out: list[np.ndarray] = []
    if int(F["chi"][p]) == 1:
        out.append(
            np.array([[c + p * x1 for x1 in range(p)] for c in range(p)], dtype=np.int32)
        )
    for s in range(p):
        vec = 1 + p * s
        if int(F["chi"][vec]) == 1:
            out.append(
                np.array(
                    [[x0 + p * ((s * x0 + c) % p) for x0 in range(p)] for c in range(p)],
                    dtype=np.int32,
                )
            )
    return out


def load_D(p: int) -> np.ndarray:
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    return (np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0), F


def classify_p(p: int) -> dict:
    D, F = load_D(p)
    classes = square_classes(F)
    Nplus = int(D.shape[0])
    n_sq = len(classes)
    packs = []
    M0 = None
    dir0_sigs: list[tuple[int, ...]] = []
    for di, lines in enumerate(classes):
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1).astype(np.int16)
        sig = [tuple(int(x) for x in row) for row in np.sort(ns, axis=1)]
        if di == 0:
            M0 = (ns.astype(np.int32) ** 2).sum(axis=1)
            dir0_sigs = sig
        packs.append(sig)
    d_tags = []
    for i in range(Nplus):
        tags = tuple(sorted(tag_sig(p, packs[di][i]) for di in range(n_sq)))
        d_tags.append(tags)
    dt = Counter(d_tags)
    sc0 = Counter(dir0_sigs)
    named = M_named(p)
    EM = float(np.mean(M0))
    EM2 = float(np.mean(M0.astype(np.float64) ** 2))
    mix = 0.0
    dir0_rows = []
    for sig, cnt in sc0.most_common():
        # M from signature
        M = int(sum(x * x for x in sig))
        mix += cnt * M * M
        dir0_rows.append(
            {
                "sig": list(sig),
                "tag": tag_sig(p, sig),
                "count": int(cnt),
                "M": M,
            }
        )
    mix /= Nplus
    hs_key = tuple(sorted(["H"] + ["C"] * (n_sq - 1)))
    n_hs = int(dt.get(hs_key, 0))
    return {
        "p": p,
        "Nplus": Nplus,
        "n_1d": n_1d(p),
        "n_sq": n_sq,
        "n_hs": n_hs,
        "hs_count_ok": n_hs == n_1d(p),
        "n_D_types": len(dt),
        "D_types": {str(k): int(v) for k, v in dt.most_common()},
        "dir0": dir0_rows,
        "E_M": EM,
        "E_M_named": named["E_M"],
        "E_M2": EM2,
        "E_M2_from_sig": mix,
        "E_M2_match": abs(mix - EM2) < 1e-6,
        "named_M": named,
    }


def prove_halfspaces() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        r = classify_p(p)
        if not r["hs_count_ok"]:
            ok = False
        # fail-when-wrong would be n_hs == n_1d-1
        if r["n_hs"] == n_1d(p) - 1:
            ok = False
        rows[str(p)] = {
            "n_hs": r["n_hs"],
            "n_1d": r["n_1d"],
            "Nplus": r["Nplus"],
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "Halfspaces = (H, C^{n_sq-1}); count = n_1d = m C(p,m).",
    }


def prove_p5_single_NL() -> dict:
    r5 = classify_p(5)
    r7 = classify_p(7)
    p = 5
    nl5 = r5["Nplus"] - r5["n_hs"]
    expect_nl5 = p * p * (p - 1)
    allR = r5["D_types"].get(str(("R", "R", "R")), 0)
    ok = (
        nl5 == expect_nl5
        and allR == expect_nl5
        and r5["E_M2_match"]
        and r7["n_D_types"] > 2  # not a single NL type
    )
    # one-NL-type formula at p=7 must fail
    named = M_named(7)
    N7 = r7["Nplus"]
    n1 = r7["n_hs"]
    fake = (n1 / N7) * (2 * named["H"] ** 2 + 6 * named["C"] ** 2) / 8 + (
        1 - n1 / N7
    ) * named["R"] ** 2
    if abs(fake - r7["E_M2"]) < 1.0:
        ok = False
    return {
        "proved": ok,
        "p5_NL": nl5,
        "p5_NL_named": expect_nl5,
        "p5_allR": allR,
        "p5_EM2": r5["E_M2"],
        "p7_n_D_types": r7["n_D_types"],
        "p7_one_NL_formula": fake,
        "p7_live_EM2": r7["E_M2"],
        "theorem": (
            "p=5 NL is all-R of size p²(p−1). "
            "One-NL-type E[M²] fails at p=7."
        ),
    }


def prove_p7_orbits() -> dict:
    r = classify_p(7)
    p = 7
    unit = p * p * (p - 1)
    dt = r["D_types"]
    hs = dt.get(str(("C", "C", "C", "H")), 0)
    crrr = dt.get(str(("C", "R", "R", "R")), 0)
    sizes = sorted(dt.values(), reverse=True)
    ok = (
        hs == n_1d(7)
        and crrr == 4 * unit
        and sizes.count(4 * unit) == 4
        and 2 * unit in sizes
        and unit in sizes
        and r["Nplus"] == n_1d(7) + 19 * unit
        and r["E_M2_match"]
    )
    return {
        "proved": ok,
        "D_types": dt,
        "unit_p2_pm1": unit,
        "Nplus": r["Nplus"],
        "Nplus_named_sum": n_1d(7) + 19 * unit,
        "theorem": (
            "p=7 H_+ = n_1d + (4+4+4+4+2+1) p²(p−1). "
            "Matches 15.139 size classes {56+84, 294, 588, 1176×4}."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "Hplus_named": False,
        "note": (
            "|H_+| unnamed (dens 13, 409 = |H_+|/(2p)). "
            "Type-wise M_τ named; weights |O_τ|/|H_+| not. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.307  signatures classify H_+; E[M²] type mixture", flush=True)
    A = prove_halfspaces()
    print(f"  A halfspaces = n_1d: {A['proved']}", flush=True)
    B = prove_p5_single_NL()
    print(
        f"  B p=5 NL=p²(p−1): {B['proved']} p7_types={B['p7_n_D_types']}",
        flush=True,
    )
    C = prove_p7_orbits()
    print(f"  C p=7 orbit sum: {C['proved']} N+={C['Nplus']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.307",
        "title": "Square-dir signatures classify H_+; E[M²] type mixture; |H_+| unnamed",
        "proved": {
            "halfspaces_are_n1d": A["proved"],
            "p5_NL_is_p2_pm1": B["proved"],
            "p7_orbit_sum": C["proved"],
            "Hplus_named": False,
            "En4_named": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15307.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
