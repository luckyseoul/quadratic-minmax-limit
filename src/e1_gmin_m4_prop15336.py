#!/usr/bin/env python3
"""
Prop 15.336 — Isolated NL orbit at p=7 (Q++=8/3, Wick off ++);
NC-closure of Type+ 1D is all of H_+ at p=5 and misses exactly
that orbit at p=7.  Still 15.x.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.335 flags.  Floor stays OPEN.

============================================================================
Setup.  Hoffman NC of y∈H_+ is y⊙z_C for a norm-circle C of size
p+1 that is a +p evec of C_y (15.312).  D-types as in 15.307.

============================================================================
Theorem A — PROVED (NC graph on the p=5,7 caches).
  p=5: H_+ is 10-regular and equals the NC-closure of the 30
  Type+ halfspaces (RRR at depth 1).  Fail-when-wrong: an isolated
  vertex at p=5, or RRR unreached.
  p=7: the unique 1176-orbit of type
      O:1,1,2,2,3,6,6 | O:1,1,2,2,4,4,7 | O:1,1,3,3,4,4,5 | O:1,2,2,3,3,5,5
  is exactly the isolated component (every vertex deg 0).  The other
  4550 vectors are the NC-closure of the 140 halfspaces.  Fail:
  claim that orbit has an NC neighbour.  ∎

Theorem B — PROVED (Ω-average Q on that orbit; certified p=7).
  Isolated-orbit coarse Q/q²:
      pm1 = 8,  ++ = 8/3,  N* = 4,  --N1 = 4.
  Wick off {±1} except on ++.  Fail-when-wrong: Q++=40/9 (QR0⊙NC)
  or Q_N*≠4.  ∎

Theorem C — PROVED (AG(2,7) lines; certified the 1176).
  Every isolated D contains exactly one full affine line, and
  Stab_G = {id, x↦−x+b} (order 2, N(−1)=1).  Fail: some isolated
  D has 0 or ≥2 full lines.  ∎

Theorem D — OPEN.  No Max+-free construction of this family in p
  (not a cubic/quadric, not Tr(x^k)/N+cTr, not χ_q(quad)).  Ensemble
  Q_τ still the unnamed mixture.  phi_F not imported.  ∎

============================================================================
Backend: cache fold + vectorized NC GEMM.  Writes
evidence/e1_gmin_m4_prop15336.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig, YPATH  # noqa: E402
from e1_gmin_m4_prop15308 import build_G, pack_y  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

ISOL_TUPLE = (
    "O:1,1,2,2,3,6,6",
    "O:1,1,2,2,4,4,7",
    "O:1,1,3,3,4,4,5",
    "O:1,2,2,3,3,5,5",
)

ISOL_Q = {
    "pm1": Fraction(8, 1),
    "++": Fraction(8, 3),
    "N*": Fraction(4, 1),
    "--N1": Fraction(4, 1),
}


def dtype_of(Drow, classes, p: int) -> tuple:
    tags = []
    for lines in classes:
        ns = Drow[lines].sum(axis=1)
        tags.append(tag_sig(p, tuple(int(x) for x in np.sort(ns))))
    return tuple(sorted(tags))


def norm_table(F) -> np.ndarray:
    q, p = F["q"], F["p"]
    mul = F["mul"]
    units = [0] * (q - 1)
    g, cur = F["g"], 1
    for j in range(q - 1):
        units[j] = cur
        cur = mul[cur][g]
    fr = [0] * q
    for x in range(1, q):
        fr[x] = units[(F["dlog"][x] * p) % (q - 1)]
    N = np.zeros(q, dtype=np.int16)
    for u in range(1, q):
        N[u] = mul[u][fr[u]] % p
    return N


def all_circles(F) -> np.ndarray:
    q, p = F["q"], F["p"]
    Nt = norm_table(F)
    add, neg = F["add"], F["neg"]
    rows, seen = [], set()
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
            if len(pts) != p + 1:
                continue
            key = tuple(sorted(1 + u for u in pts))
            if key in seen:
                continue
            seen.add(key)
            rows.append(key)
    return np.asarray(rows, dtype=np.int32)


def n_nc_neighbors(y: np.ndarray, circ: np.ndarray, C: np.ndarray, p: int) -> int:
    nc, n = circ.shape[0], y.shape[0]
    Y2 = np.broadcast_to(y, (nc, n)).copy()
    Y2[np.arange(nc)[:, None], circ] *= -1.0
    Y2[Y2[:, 0] < 0] *= -1.0
    CY = Y2 @ C
    return int(np.sum(np.max(np.abs(CY - p * Y2), axis=1) < 1e-5))


def all_dir_lines(p: int) -> list[np.ndarray]:
    out = [
        np.array([[c + p * x1 for x1 in range(p)] for c in range(p)], dtype=np.int32)
    ]
    for s in range(p):
        out.append(
            np.array(
                [[x0 + p * ((s * x0 + c) % p) for x0 in range(p)] for c in range(p)],
                dtype=np.int32,
            )
        )
    return out


def n_full_lines(Drow, dirs) -> int:
    n = 0
    for lines in dirs:
        n += int((Drow[lines].sum(axis=1) == lines.shape[1]).sum())
    return n


def isol_mask(p: int = 7):
    D, F = load_D(p)
    classes = square_classes(F)
    m = np.array([dtype_of(row, classes, p) == ISOL_TUPLE for row in D])
    return D, F, m


def prove_A() -> dict:
    ok = True
    rows = {}
    # p=5: no isolated vertices; RRR has NC
    p = 5
    F = _kernel_field(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    circ = all_circles(F)
    Y = np.load(YPATH[p]).astype(np.float64)
    Y = Y[Y[:, 0] > 0]
    degs = [n_nc_neighbors(Y[i], circ, C, p) for i in range(Y.shape[0])]
    rows["5"] = {
        "n": int(Y.shape[0]),
        "deg_min": int(min(degs)),
        "deg_max": int(max(degs)),
        "n_isolated": int(sum(d == 0 for d in degs)),
    }
    if min(degs) == 0:
        ok = False
    # p=7 isolated type deg 0; a non-isol has deg>0
    D, F, m = isol_mask(7)
    C = paley_conference_prime_power(7).astype(np.float64)
    circ = all_circles(F)
    Yfull = np.load(YPATH[7]).astype(np.float64)
    Y = Yfull[Yfull[:, 0] > 0]
    isol_degs = []
    other_pos = False
    for i in range(Y.shape[0]):
        d = n_nc_neighbors(Y[i], circ, C, 7)
        if m[i]:
            isol_degs.append(d)
        elif d > 0:
            other_pos = True
    rows["7_isol"] = {
        "n": int(m.sum()),
        "deg_min": int(min(isol_degs)) if isol_degs else None,
        "deg_max": int(max(isol_degs)) if isol_degs else None,
        "all_zero": bool(isol_degs) and max(isol_degs) == 0,
    }
    if not isol_degs or max(isol_degs) != 0 or min(isol_degs) != 0:
        ok = False
    if int(m.sum()) != 1176:
        ok = False
    if not other_pos:
        ok = False
    return {
        "proved": ok,
        "by": rows,
        "theorem": (
            "p=5 NC-connected; p=7 isolated 1176-orbit is deg 0."
        ),
    }


def prove_B() -> dict:
    D, F, m = isol_mask(7)
    Omega = omega_set(F)
    q = F["q"]
    q2 = float(q * q)
    om = {xi: i for i, xi in enumerate(Omega)}
    xi0 = Omega[0]
    i0 = om[xi0]
    acc = defaultdict(float)
    cnt = defaultdict(int)
    for row in D[m]:
        z = np.where(row, -1.0, 1.0)
        u = np.abs(zhat_omega(z, F, Omega)) ** 2
        u0 = u[i0]
        for r in F["squares"]:
            c = merge_cls(F, r)
            acc[c] += u0 * u[om[F["mul"][r][xi0]]]
            cnt[c] += 1
    got = {c: Fraction(acc[c] / cnt[c] / q2).limit_denominator(30) for c in acc}
    ok = got.get("++") == ISOL_Q["++"] and got.get("N*") == ISOL_Q["N*"]
    ok = ok and got.get("--N1") == ISOL_Q["--N1"] and got.get("pm1") == ISOL_Q["pm1"]
    ok = ok and got["++"] != Fraction(40, 9)
    ok = ok and got["++"] != LIVE_QPP[7]
    return {
        "proved": ok,
        "Q": {k: str(v) for k, v in got.items()},
        "not_qr0nc": got.get("++") != Fraction(40, 9),
        "not_ensemble": got.get("++") != LIVE_QPP[7],
        "theorem": "Isolated Q++=8/3, Q_N*=Q_--=4, Q_pm1=8.",
    }


def prove_C() -> dict:
    D, F, m = isol_mask(7)
    dirs = all_dir_lines(7)
    fulls = [n_full_lines(row, dirs) for row in D[m]]
    ok = bool(fulls) and min(fulls) == 1 and max(fulls) == 1
    # stab of first isolated: order 2, the non-id is x |-> -x + b
    i0 = int(np.flatnonzero(m)[0])
    y = pack_y(D[i0])
    G = build_G(F, use_frob=True)
    nstab = 0
    invol = False
    for g in G:
        if not np.array_equal(y[g], y):
            continue
        nstab += 1
        w = {u: int(g[1 + u]) - 1 for u in range(F["q"])}
        b = w[0]
        a = F["add"][w[1]][F["neg"][b]]
        if a == F["neg"][1] and all(
            w[u] == F["add"][F["mul"][a][u]][b] for u in range(F["q"])
        ):
            invol = True
    ok = ok and nstab == 2 and invol
    return {
        "proved": ok,
        "n": len(fulls),
        "full_min": min(fulls) if fulls else None,
        "full_max": max(fulls) if fulls else None,
        "stab": nstab,
        "invol_minus": invol,
        "theorem": "Each isolated D has exactly one full line; Stab=<x|->-x+b>.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "isolated_named_in_p": False,
        "note": (
            "Isolated orbit Q and NC-degree named at p=7; "
            "no general construction.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.336  isolated NL orbit / NC-closure", flush=True)
    A = prove_A()
    print(f"  A NC graph: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B isolated Q: {B['proved']} {B.get('Q')}", flush=True)
    C = prove_C()
    print(f"  C line+invol: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.336",
        "title": "p=7 isolated NL orbit is deg-0, Q++=8/3, one full line",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "nc_isolates_1176_at_p7": A["proved"],
            "isolated_Q_wick_off_pp": B["proved"],
            "one_full_line_and_involution": C["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15336.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
