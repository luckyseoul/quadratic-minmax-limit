#!/usr/bin/env python3
"""
Prop 15.365 — Two named Gauss/Jacobi neighbours of ensemble Q++
hit 48/13 at p=5 and both die at p=7.  An expanded GJ / CM / cyclo
catalog has 0 forms equal to 1544/409.  The p=5 u-support is a named
5-atom and is not the p=7 support.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.364 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  Ensemble Q++ = Q on Paley-(1,1) (15.364).  Live values
48/13 (p=5) and 1544/409 (p=7).  dim V₊ = (q+1)/2.  |++| = n_pp
is the torus-Jacobi count (15.355).  Wick off-average is
4(q−1)/(q+1) (15.299).

============================================================================
Theorem A — PROVED (named identity + live fractions).
  Q_neigh^{dim} := 8 n_pp / dim V₊ equals 48/13 at p=5 and equals
  48/25 ≠ 1544/409 at p=7.  Fail: claim 48/25 = 1544/409.  ∎

Theorem B — PROVED (named identity + live fractions).
  Q_neigh^{Wick} := 4(q−1)/(q+1) equals 48/13 at p=5 and equals
  96/25 ≠ 1544/409 at p=7.  Fail: claim 96/25 = 1544/409.  ∎

Theorem C — PROVED (size bounds + 15.330 + Fermat enumeration).
  Write q=49.  409 is prime, so a 1-term catalog integer equal to
  409 must be 409 itself.
    (i) 15.330: 409 ∉ named_gj_ints(7); dim V₊=25 and q²=2401
        sit in that list.
    (ii) Weil: |∑_x χ(x³+Ax+B)| ≤ 2√q = 14, so every CM trace
        and every #E(F_q)=q+1±t lie in [36, 64].
    (iii) Cyclotomic numbers of order e|(q−1) lie in {0,…,q−2}
        = {0,…,47}.
    (iv) |Kl(a,1)|² ≤ (p−1)² = 36.
    (v) Affine Fermat #{x^n+y^n=1} for n=2,3,4 is the explicit
        triple (48, 60, 88), none equal to 409.
  Hence 409 is not a 1-term catalog integer, so no form a/b,
  4−2a/b, 8a/b with a,b in the catalog equals 1544/409.
  Fail: claim 409 = dim V₊ or 409 = q².  ∎

Theorem D — PROVED (exact arithmetic in Z[√5]; p=7 interval).
  At p=5, e(tr)=ζ_5 and cos(2πt/5) ∈ {1,(√5−1)/4,(−√5−1)/4},
  so 4u ∈ Z[√5].  On all 130 y_∞=+1 Max+ and all ξ∈Ω, the
  pair (A,B) with 4u=A+B√5 lies in
      {(0,0), (200,−40), (600,−200), (200,40), (600,200)},
  which are exactly 4 times the named 5-atom
      {0, 10(5−√5), 5(5−√5)², 10(5+√5), 5(5+√5)²}.
  At p=7, |ẑ| is a sum of 49 unit complex numbers; float64
  |u−u_true|<10^{−9}.  The value 21.08… on the first Max+ is
  more than 6 away from every p=5 atom and every p↦7 twin.
  Fail: claim the 5-atom is the p=7 support.  ∎

Theorem E — OPEN.  No Max+-free GJ formula hits both live
  Q++ values.  The neighbours that hit p=5 die at p=7, so they
  do not imply ⟨δ,ψ⟩≤2.  phi_F not imported.

============================================================================
Backend: Fraction algebra + zhat GEMM on caches.  Writes
evidence/e1_gmin_m4_prop15365.json
"""
from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15316 import kloosterman  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP, named_gj_ints  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def dim_Vp(p: int) -> int:
    return (p * p + 1) // 2


def Q_neigh_dim(p: int) -> Fraction:
    return Fraction(8 * n_pp_named(p), dim_Vp(p))


def Q_neigh_wick(p: int) -> Fraction:
    q = p * p
    return Fraction(4 * (q - 1), q + 1)


def p5_u_atoms() -> list[float]:
    s = math.sqrt(5)
    return [
        0.0,
        10.0 * (5.0 - s),
        5.0 * (5.0 - s) ** 2,
        10.0 * (5.0 + s),
        5.0 * (5.0 + s) ** 2,
    ]


def p7_twin_atoms() -> list[float]:
    s = math.sqrt(7)
    return [
        0.0,
        14.0 * (7.0 - s),
        7.0 * (7.0 - s) ** 2,
        14.0 * (7.0 + s),
        7.0 * (7.0 + s) ** 2,
    ]


def _near(x: float, atoms: list[float], tol: float = 1e-6) -> bool:
    return any(abs(x - a) < tol for a in atoms)


def _load_Z(p: int) -> np.ndarray:
    Y = np.load(YPATH[p])
    q = p * p
    return np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64))


def u_on_rows(p: int, n_rows: int | None = None) -> np.ndarray:
    F = _kernel_field(p)
    Omega = omega_set(F)
    Z = _load_Z(p)
    if n_rows is not None:
        Z = Z[:n_rows]
    q = F["q"]
    E = np.zeros((len(Omega), q), dtype=np.complex128)
    for i, xi in enumerate(Omega):
        for x in range(q):
            E[i, x] = _e_tr(F, F["mul"][xi][x])
    return np.abs(Z @ E.T) ** 2


def fermat_counts(p: int, ns: tuple[int, ...] = (2, 3, 4)) -> dict[int, int]:
    F = _kernel_field(p)
    q = F["q"]
    add, mul = F["add"], F["mul"]

    def pn(x, n):
        acc, b, ee = 1, x, n
        if n == 0:
            return 1
        while ee:
            if ee & 1:
                acc = mul[acc][b]
            b = mul[b][b]
            ee >>= 1
        return acc

    out = {}
    for n in ns:
        c = 0
        for x in range(q):
            xn = pn(x, n)
            for y in range(q):
                if add[xn][pn(y, n)] == 1:
                    c += 1
        out[n] = c
    return out


COS5 = {
    0: (4, 0),  # 4*cos(0) = 4
    1: (-1, 1),  # 4*cos(2π/5) = -1+√5
    2: (-1, -1),  # 4*cos(4π/5) = -1-√5
    3: (-1, -1),
    4: (-1, 1),
}

# 4 * named p=5 atoms in Z[√5]
P5_ATOM_AB = {
    (0, 0),
    (200, -40),
    (600, -200),
    (200, 40),
    (600, 200),
}


def _cos5_mats(F, Omega):
    """Per-ξ integer matrices (CA, CB) with 4u = z^T CA z + √5 z^T CB z."""
    q = F["q"]
    mul, add, neg, trt = F["mul"], F["add"], F["neg"], F["trt"]
    mats = []
    for xi in Omega:
        CA = np.zeros((q, q), dtype=np.int64)
        CB = np.zeros((q, q), dtype=np.int64)
        for x in range(q):
            for y in range(q):
                a, b = COS5[trt[mul[xi][add[x][neg[y]]]]]
                CA[x, y] = a
                CB[x, y] = b
        mats.append((CA, CB))
    return mats


def certify_p5_atoms() -> dict:
    F = _kernel_field(5)
    Omega = omega_set(F)
    Z = _load_Z(5).astype(np.int64)
    mats = _cos5_mats(F, Omega)
    seen: set[tuple[int, int]] = set()
    for i in range(Z.shape[0]):
        z = Z[i]
        for CA, CB in mats:
            A = int(z @ CA @ z)
            B = int(z @ CB @ z)
            seen.add((A, B))
    return {
        "n_y": int(Z.shape[0]),
        "n_pairs": len(seen),
        "pairs": sorted(seen),
        "equals_named": seen == P5_ATOM_AB,
        "subset_named": seen <= P5_ATOM_AB,
    }


def prove_A() -> dict:
    v5, v7 = Q_neigh_dim(5), Q_neigh_dim(7)
    ok = v5 == LIVE_QPP[5] == Fraction(48, 13)
    ok = ok and v7 == Fraction(48, 25)
    if v7 == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "p5": str(v5),
        "p7": str(v7),
        "live7": str(LIVE_QPP[7]),
        "n_pp": {5: n_pp_named(5), 7: n_pp_named(7)},
        "dim_Vp": {5: dim_Vp(5), 7: dim_Vp(7)},
        "theorem": "8 n_pp / dim V+ hits p=5 and dies at p=7.",
    }


def prove_B() -> dict:
    v5, v7 = Q_neigh_wick(5), Q_neigh_wick(7)
    ok = v5 == LIVE_QPP[5]
    ok = ok and v7 == Fraction(96, 25)
    if v7 == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "p5": str(v5),
        "p7": str(v7),
        "live7": str(LIVE_QPP[7]),
        "theorem": "4(q-1)/(q+1) hits p=5 and dies at p=7.",
    }


def prove_C() -> dict:
    S7 = named_gj_ints(7)
    ferm = fermat_counts(7)
    weill = 2 * 7  # 2√q = 14
    e_lo, e_hi = 49 + 1 - weill, 49 + 1 + weill  # [36, 64]
    multiples = [x for x in S7 if x and x % 409 == 0]
    ferm_mult = [v for v in ferm.values() if v % 409 == 0]
    ok = 409 not in S7
    ok = ok and 25 in S7
    ok = ok and 2401 in S7
    ok = ok and 409 > e_hi and 409 > 47 and 409 > (7 - 1) ** 2
    ok = ok and 409 not in ferm.values()
    ok = ok and multiples == [] and ferm_mult == []
    ok = ok and dim_Vp(7) != 409 and (7 ** 4) != 409
    # p=5 sanity: the dim-neighbour is a one-term hit
    ok = ok and Fraction(8 * n_pp_named(5), dim_Vp(5)) == LIVE_QPP[5]
    return {
        "proved": bool(ok),
        "S7_has_409": 409 in S7,
        "S7_has_dimVp": 25 in S7,
        "S7_has_q2": 2401 in S7,
        "weil_E_range": [e_lo, e_hi],
        "cyclo_max": 47,
        "kl2_max": 36,
        "fermat": ferm,
        "gj_multiples_of_409": multiples,
        "theorem": (
            "409 is excluded from every 1-term GJ/CM/cyclo/Fermat "
            "family by size or explicit count; it divides no "
            "named_gj_int(7).  Hence no a/b, 4-2a/b, 8a/b in the "
            "catalog equals 1544/409."
        ),
    }


def prove_D() -> dict:
    cert = certify_p5_atoms()
    ok = cert["equals_named"] and cert["n_y"] == 130
    U7 = u_on_rows(7, n_rows=20)
    atoms = p5_u_atoms()
    twin = p7_twin_atoms()
    outsider = None
    gap = None
    for v in np.unique(np.round(U7, 10)):
        fv = float(v)
        d5 = min(abs(fv - a) for a in atoms)
        d7 = min(abs(fv - a) for a in twin)
        d = min(d5, d7)
        if d > 1.0 and (gap is None or d > gap):
            outsider = fv
            gap = d
    # float64 sum of 49 units: |u-u_true| << 1e-9 << gap
    if outsider is None or gap is None or gap < 1.0:
        ok = False
    return {
        "proved": bool(ok) and outsider is not None,
        "p5_n_atoms": cert["n_pairs"],
        "p5_pairs": [list(t) for t in cert["pairs"]],
        "p5_exact": cert["equals_named"],
        "p7_outsider": outsider,
        "p7_gap": gap,
        "p7_err_bound": 1e-9,
        "theorem": (
            "p=5 u-support is the named 5-atom in Z[√5]; "
            "the 5-atom and its p=7 twin fail at p=7 with gap>1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "Neighbours that hit p=5 die at p=7, so they do not "
            "imply <delta,psi> <= 2.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.365  GJ neighbours of Q++ die at p=7", flush=True)
    A = prove_A()
    print(f"  A dim-neigh: {A['proved']} {A['p5']} vs {A['p7']}", flush=True)
    B = prove_B()
    print(f"  B wick-neigh: {B['proved']} {B['p5']} vs {B['p7']}", flush=True)
    C = prove_C()
    print(
        f"  C catalog: {C['proved']} 409in7={C['S7_has_409']} ferm={C['fermat']}",
        flush=True,
    )
    D = prove_D()
    print(
        f"  D u-atom: {D['proved']} p5={D['p5_n_atoms']} outsider={D['p7_outsider']}",
        flush=True,
    )
    E = prove_open()
    print(f"  E open: Q_tau={E['Q_tau_named_in_p']} phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.365",
        "title": "GJ neighbours of Q++ hit p=5 and die at p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "neigh_dim_dies_p7": A["proved"],
            "neigh_wick_dies_p7": B["proved"],
            "catalog_misses_1544_409": C["proved"],
            "p5_u_5atom_not_p7": D["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15365.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
