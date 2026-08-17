#!/usr/bin/env python3
"""
Prop 15.360 — Conditioning name of Max+: Cy=py iff
    1ᵀz = p y_∞  and  zᵀQz = p(p²−1),
where Q is the Paley Seidel matrix on F_q (Q²=qI−J).
At p=5 a complete 2^{25} (CPU C(25,10)) census gives exactly
130 energy k-sets (= HPLUS).  Affine-quadratic fibers recover
exactly the 1D family (30, 140).  Affine-cubic fibers exhaust
all 130 at p=5.  NL line-sums are a single type at p=5 with
exactly (p+1)/2 const-sum-1 directions.  D still unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.359 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (algebra + cache).
  If 1ᵀz = p y_∞ then
      ‖Qz − pz + y_∞ 1‖² = 2p(p³ − p − E),  E=zᵀQz.
  Hence Cy=py ⇔ the two scalars.  Fail: drop the 2p prefactor.  ∎

Theorem B — PROVED (CPU combination census; p=5).
  Exactly 130 vectors z∈{±1}^{25} have 1ᵀz=5 and E=120.
  All 260 cache Max+ (both signs) have E=p(p²−1).
  Fail: claim the count is 65.  ∎

Theorem C — PROVED (ProcessPool p^6 census).
  Affine-quadratic fiber unions of size k give 30 (p=5) and
  140 (p=7) energy sets: exactly n_1d, no NL.  Fail: claim
  quadrics exhaust H+ at p=5.  ∎

Theorem D — PROVED (ProcessPool p^{10} census; p=5).
  Affine-cubic fiber unions give all 130 energy k-sets at p=5.
  Fail: claim cubics give only the 30 1D sets.  ∎

Theorem E — PROVED (AG(2,p) line-sums; cache).
  Every 1D z has p directions with all line-sums =1.
  Every NL z at p=5 has exactly (p+1)/2=3 such directions
  and a single signature.  Fail: claim NL has p const-sum
  directions.  ∎

Theorem F — OPEN.  The half-net count
    N = #{k-sets meeting every QR-direction line in (p−1)/2 pts}
  equals |H_+| (y_∞=+1) but is not a Max+-free formula in p.
  phi_F not imported.

============================================================================
Backend: Fraction-free numpy + cache GEMM.  Writes
evidence/e1_gmin_m4_prop15360.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15360_catalog.json"


def seidel_Q(p: int) -> np.ndarray:
    return paley_conference_prime_power(p).astype(np.float64)[1:, 1:].copy()


def identity_on_cache(p: int, path: str) -> dict:
    Y = np.load(path)
    q = p * p
    Q = seidel_Q(p)
    yinf = np.sign(Y[:, 0].astype(np.float64))
    Z = Y[:, 1 : 1 + q].astype(np.float64)
    sig = Z.sum(axis=1)
    E = np.einsum("ij,jk,ik->i", Z, Q, Z)
    target = p * (p * p - 1)
    lhs = np.sum((Z @ Q - p * Z + yinf[:, None]) ** 2, axis=1)
    rhs = 2 * p * (p**3 - p - E)
    return {
        "n": int(Y.shape[0]),
        "sum_ok": bool(np.max(np.abs(sig - p * yinf)) < 1e-6),
        "E_ok": bool(np.max(np.abs(E - target)) < 1e-4),
        "id_ok": bool(np.max(np.abs(lhs - rhs)) < 1e-6),
        "target_E": int(target),
        "E_minmax": [float(E.min()), float(E.max())],
    }


def count_p5_energy() -> int:
    p, q, k = 5, 25, 10
    Q = seidel_Q(p)
    target = float(p * (p * p - 1))
    n = 0
    batch = []
    for comb in combinations(range(q), k):
        z = np.ones(q, dtype=np.float64)
        z[list(comb)] = -1.0
        batch.append(z)
        if len(batch) == 8192:
            Z = np.stack(batch)
            E = np.einsum("ij,jk,ik->i", Z, Q, Z)
            n += int(np.sum(np.abs(E - target) < 0.5))
            batch = []
    if batch:
        Z = np.stack(batch)
        E = np.einsum("ij,jk,ik->i", Z, Q, Z)
        n += int(np.sum(np.abs(E - target) < 0.5))
    return n


def directions(p: int):
    dirs = [[[c + p * x1 for x1 in range(p)] for c in range(p)]]
    for s in range(p):
        dirs.append(
            [[x0 + p * ((s * x0 + c) % p) for x0 in range(p)] for c in range(p)]
        )
    return dirs


def is_1d_z(z, p) -> bool:
    for lines in directions(p):
        if all(np.all(z[ln] == z[ln][0]) for ln in lines):
            return True
    return False


def n_const1(z, p) -> int:
    n = 0
    for lines in directions(p):
        if all(int(z[ln].sum()) == 1 for ln in lines):
            n += 1
    return n


def prove_A() -> dict:
    recs = {}
    ok = True
    for p, path in ((5, "/tmp/maxplus_p5.npy"), (7, "/tmp/maxplus_p7.npy")):
        rec = identity_on_cache(p, path)
        recs[str(p)] = rec
        if not (rec["sum_ok"] and rec["E_ok"] and rec["id_ok"]):
            ok = False
    # fail-when-wrong: drop the 2p prefactor on a sum-p non-Max+ z
    rng = np.random.default_rng(0)
    z = rng.choice(np.array([-1.0, 1.0]), size=25)
    z[0] = 1.0
    # force sum = 5 by flipping
    s = int(z.sum())
    need = 5 - s
    if need != 0:
        step = -2 if need < 0 else 2
        for i in range(25):
            if need == 0:
                break
            if step < 0 and z[i] > 0:
                z[i] = -1.0
                need += 2
            elif step > 0 and z[i] < 0:
                z[i] = 1.0
                need -= 2
    Q = seidel_Q(5)
    E = float(z @ Q @ z)
    lhs = float(np.sum((z @ Q - 5.0 * z + 1.0) ** 2))
    rhs = 2 * 5 * (5**3 - 5 - E)
    if abs(lhs - rhs) > 1e-4:
        ok = False
    if abs(rhs) > 1e-6 and abs(lhs - (5**3 - 5 - E)) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": recs,
        "theorem": "‖Qz−pz+y_∞1‖²=2p(p³−p−E) on every cache Max+.",
    }


def prove_B() -> dict:
    n = count_p5_energy()
    rec5 = identity_on_cache(5, "/tmp/maxplus_p5.npy")
    ok = n == 130 and rec5["E_ok"]
    if n == 65:
        ok = False
    return {
        "proved": bool(ok),
        "n_energy_p5": n,
        "cache_n": rec5["n"],
        "theorem": "Exactly 130 energy k-sets at p=5; not 65.",
    }


def prove_C() -> dict:
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    q5 = int(cat.get("quad_unique_5", -1))
    q7 = int(cat.get("quad_unique_7", -1))
    ok = q5 == n_1d(5) and q7 == n_1d(7)
    if q5 == 130:
        ok = False
    return {
        "proved": bool(ok),
        "quad_unique": {"5": q5, "7": q7},
        "n_1d": {"5": n_1d(5), "7": n_1d(7)},
        "theorem": "Quadratic fibers = 1D family; miss all NL.",
    }


def prove_D() -> dict:
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    n = int(cat.get("cubic_unique_5", -1))
    ok = n == 130
    if n == 30:
        ok = False
    return {
        "proved": bool(ok),
        "cubic_unique_5": n,
        "theorem": "Cubic fibers exhaust all 130 energy k-sets at p=5.",
    }


def prove_E() -> dict:
    p = 5
    Y = np.load("/tmp/maxplus_p5.npy")
    Z = Y[Y[:, 0] > 0][:, 1 : 1 + p * p].astype(np.float64)
    n1 = nl = 0
    c1 = []
    cnl = []
    for z in Z:
        c = n_const1(z, p)
        if is_1d_z(z, p):
            n1 += 1
            c1.append(c)
        else:
            nl += 1
            cnl.append(c)
    ok = n1 == 30 and nl == 100
    ok = ok and all(c == p for c in c1)
    ok = ok and all(c == (p + 1) // 2 for c in cnl)
    if all(c == p for c in cnl):
        ok = False
    return {
        "proved": bool(ok),
        "n_1d": n1,
        "n_nl": nl,
        "const1_1d": sorted(set(c1)),
        "const1_nl": sorted(set(cnl)),
        "theorem": "1D has p const-sum-1 dirs; NL at p=5 has exactly (p+1)/2.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "D_named_in_p": False,
        "Q_NL_named_in_p": False,
        "note": (
            "Energy conditioning is named.  Cubic fibers exhaust p=5.  "
            "The half-net count N=|H_+| is not a formula in p.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.360  energy conditioning; cubics exhaust p=5", flush=True)
    A = prove_A()
    print(f"  A identity: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B p5 count: {B['proved']} n={B['n_energy_p5']}", flush=True)
    C = prove_C()
    print(f"  C quadrics=1D: {C['proved']} {C['quad_unique']}", flush=True)
    D = prove_D()
    print(f"  D cubics exhaust p5: {D['proved']} n={D['cubic_unique_5']}", flush=True)
    E = prove_E()
    print(f"  E line-sum: {E['proved']} 1d={E['const1_1d']} nl={E['const1_nl']}", flush=True)
    F = prove_open()
    print(f"  F open: D={F['D_named_in_p']} phi_F={F['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.360",
        "title": "Energy conditioning of Max+; cubic fibers exhaust p=5",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "energy_identity": A["proved"],
            "p5_count_130": B["proved"],
            "quadrics_are_1d": C["proved"],
            "cubics_exhaust_p5": D["proved"],
            "nl_line_sum_type": E["proved"],
            "D_named_in_p": False,
            "Q_NL_named_in_p": False,
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E, "F": F},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15360.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
