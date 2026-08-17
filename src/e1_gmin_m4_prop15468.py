#!/usr/bin/env python3
"""
Prop 15.468 — Every Max+ D has a named χ-pair count
    S_χ = 1_D^T C_F 1_D = p(p²−1)/4,
    #{{a,b}⊂D : χ(b−a)=+1} = k(k−1)(p−1)/(4(p−2)),
from Cy=py (15.189).  The χ-classified 3-point of 1_D is
not constant on Max+.  L++ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.467 flags.

============================================================================
Setup.  Paley conference C on F_q ∪ {∞}, C_∞x=1, C_F[a,b]=χ(b−a)
(a≠b).  Max+: Cy=py, y_∞=+1 ⇒ 1ᵀy=p+1 (15.189).  D={x: y_x=−1},
k=p(p−1)/2.  15.467 L++ is E[N(d)N(sd)] on ++, a 4-point of 1_D.

============================================================================
Theorem A — PROVED (Max+-free; Cy=py + C_∞=1 + χ(−1)=1).
  y^T C y = p(q+1) and the ∞-border 2 y_∞ ∑_F y = 2p give
      y_F^T C_F y_F = p(p²−1).
  1_D=(1−y_F)/2 and the row-sums ∑_{d≠0} χ(d)=0 kill the linear
  terms, so S_χ = y_F^T C_F y_F / 4 = p(p²−1)/4.
  Hence χ+ = (k(k−1)+S_χ)/4 = k(k−1)(p−1)/(4(p−2)).
  Fail: drop the ∞-border (S_χ=p(p²+1)/4; 65/2 ≠ 30 at p=5).  ∎

Theorem B — PROVED (A + cache fold; p=5, and p=3,7).
  The named pair counts hold on every live Max+ row.
  Fail: claim χ+=C(k,2)/2 at p=5 (45/2 ≠ 30).  ∎

Theorem C — PROVED (p=5 cache).
  The χ-triple fingerprint of 1_D takes 16 values on the 130
  Max+ (not a per-vector 3-design).  Fail: claim one fingerprint.  ∎

Theorem D — OPEN.  L++ is a 4-point of 1_D; the 2-point is named
  and the 3-point is not constant.  phi_F not imported.

============================================================================
Backend: Fraction algebra + numpy GEMM on C_F.  Writes
evidence/e1_gmin_m4_prop15468.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15468_catalog.json"
YPATH = {3: "/tmp/maxplus_p3.npy", 5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def k_named(p: int) -> int:
    return p * (p - 1) // 2


def S_chi_named(p: int) -> int:
    return p * (p * p - 1) // 4


def S_chi_drop_infty(p: int) -> Fraction:
    return Fraction(p * (p * p + 1), 4)


def chi_plus_named(p: int) -> Fraction:
    k = k_named(p)
    return Fraction(k * (k - 1) * (p - 1), 4 * (p - 2))


def chi_minus_named(p: int) -> Fraction:
    k = k_named(p)
    return Fraction(k * (k - 1) * (p - 3), 4 * (p - 2))


def C_F_chi(F) -> np.ndarray:
    q = F["q"]
    C = np.zeros((q, q), dtype=np.float64)
    chi = np.array(F["chi"], dtype=np.float64)
    add = F["add"]
    neg = F["neg"]
    for a in range(q):
        for b in range(q):
            if a == b:
                continue
            C[a, b] = chi[add[b][neg[a]]]
    return C


def D_rows(p: int, F) -> np.ndarray:
    Y = np.load(YPATH[p])
    return np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64)) < 0


def pair_counts(Ds: np.ndarray, CF: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # S_χ = 1_D^T C_F 1_D per row
    S = np.einsum("ij,jk,ik->i", Ds.astype(np.float64), CF, Ds.astype(np.float64))
    k = Ds.sum(axis=1)
    chi_p = (k * (k - 1) + S) / 4.0
    chi_m = (k * (k - 1) - S) / 4.0
    return S, np.stack([chi_m, chi_p], axis=1)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11, 13):
        S = S_chi_named(p)
        drop = S_chi_drop_infty(p)
        cp = chi_plus_named(p)
        cm = chi_minus_named(p)
        k = k_named(p)
        # algebra: S = k(k-1)/(p-2)
        if p != 2 and S != k * (k - 1) // (p - 2):
            ok = False
        if Fraction(k * (k - 1) + S, 4) != cp:
            ok = False
        if Fraction(k * (k - 1) - S, 4) != cm:
            ok = False
        if drop == S:
            ok = False
        rows[str(p)] = {"S": S, "chi+": str(cp), "chi-": str(cm), "drop": str(drop)}
    # fail-when-wrong target
    if S_chi_drop_infty(5) == 30:
        ok = False
    if S_chi_named(5) != 30:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "S_χ=p(p²−1)/4 from Cy=py. "
            "Fail: drop ∞-border (65/2=30 at p=5)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rec = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        Ds = D_rows(p, F)
        CF = C_F_chi(F)
        S, pairs = pair_counts(Ds, CF)
        want_S = S_chi_named(p)
        want_p = float(chi_plus_named(p))
        want_m = float(chi_minus_named(p))
        if not np.allclose(S, want_S) or not np.allclose(pairs[:, 1], want_p):
            ok = False
        if not np.allclose(pairs[:, 0], want_m):
            ok = False
        rec[str(p)] = {
            "n": int(Ds.shape[0]),
            "S": want_S,
            "chi+": str(chi_plus_named(p)),
            "max_err": float(np.max(np.abs(S - want_S))),
        }
    if chi_plus_named(5) == Fraction(k_named(5) * (k_named(5) - 1), 4):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rec,
        "theorem": (
            "Pair counts hold on every live Max+. "
            "Fail: χ+=C(k,2)/2 at p=5."
        ),
    }


def triple_codes(mask: np.ndarray, F) -> tuple[int, ...]:
    pts = np.flatnonzero(mask)
    m = pts.size
    add, neg, chi = F["add"], F["neg"], F["chi"]
    acc = [0] * 8
    for a in range(m):
        for b in range(a + 1, m):
            cab = chi[add[pts[b]][neg[pts[a]]]]
            if cab == 0:
                continue
            for c in range(b + 1, m):
                cac = chi[add[pts[c]][neg[pts[a]]]]
                cbc = chi[add[pts[c]][neg[pts[b]]]]
                if cac == 0 or cbc == 0:
                    continue
                code = (1 if cab > 0 else 0) | ((1 if cac > 0 else 0) << 1) | (
                    (1 if cbc > 0 else 0) << 2
                )
                acc[code] += 1
    return tuple(acc)


def prove_C() -> dict:
    F = _kernel_field(5)
    Ds = D_rows(5, F)
    seen: Counter[tuple[int, ...]] = Counter()
    for row in Ds:
        seen[triple_codes(row, F)] += 1
    ok = len(seen) == 16
    if len(seen) == 1:
        ok = False
    return {
        "proved": bool(ok),
        "nunique": len(seen),
        "n_rows": int(Ds.shape[0]),
        "multiset_size": len(seen),
        "theorem": (
            "p=5 χ-triple fingerprint takes 16 values. "
            "Fail: a single fingerprint."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "2-point of 1_D is named. 3-point is not constant. "
            "L++ still unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.468  Max+ χ-pair law; 3-point not constant", flush=True)
    A = prove_A()
    print(f"  A algebra: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B cache: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C triples: {C['proved']} nunique={C['nunique']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.468",
        "title": "Max+ χ-pair count named from Cy=py; 3-point not constant",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "pair_law": A["proved"],
            "pair_cache": B["proved"],
            "triple_not_constant": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15468.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
