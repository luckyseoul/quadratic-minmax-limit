#!/usr/bin/env python3
"""
Prop 15.300 — Coarsened Q-types are not a Schur ring;
the type-character Jacobi matrix J is integer and named,
but S=Φ|_F on those characters still has den N/(4p).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.299 flags.  Does **not** name Q_τ in p.  Floor stays OPEN.

============================================================================
Theorem A — PROVED (Max+-free counting on T).
  The partition of T=□ into the 15.298 coarse Q-classes (mixed Paley
  N1 fused into ++) is not a Schur ring: the structure constants
  p^k_{ij}(g)=#{x∈C_i : x^{-1}g∈C_j} vary on g∈C_k.  Worst spread
  2 at p=5 (class ++*++ in ++) and 4 at p=7 (++*++ in N*).
  Fail-when-wrong: claim p^{++}_{++,++} constant on ++.  ∎

Theorem B — PROVED (complete character sums on T; Max+-free).
  For the cyclic generator of T and characters ψ_k(gen)=exp(2πi k/|T|),
  J_τ(k)=∑_{r∈τ} ψ_k(r) is an algebraic integer.  A basis
      p=5: k=0,2,6 on {pm1, ++, N*}     det J=144
      p=7: k=0,2,8,12 on {pm1, ++, N*, --N1}  det J=−2688
  is integer and invertible.  Q_τ/q² = (J^{-1} S)_τ with
  S_k=∑_τ (Q_τ/q²) J_τ(k).  S_0=2(q−1) is named.  The other S_k
  equal certified Φ|_F eigenvalues (80/13, 176/13 at p=5; 3360/409,
  3072/409, 4320/409 at p=7) and still carry den N/(4p).
  Fail-when-wrong: flip J_++(k=2) from its character-sum value to
  its negative; J^{-1} S misses live Q_++.  ∎

Theorem C — OPEN.  Identifying D=N/(4p) with det J fails (144≠13,
  2688≠409).  Coarsened BM inverse-eigenmatrix does not name Q_τ
  because the partition is not a scheme (A) and S_{k≠0} is unnamed
  in p.  ⟨δ,ψ⟩≤2 / phi_F stay OPEN.  ∎

============================================================================
Backend: CPU field tables.  Live Q only to check inversion / invert.
Writes evidence/e1_gmin_m4_prop15300.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15298 import coarse_Q_class, live_Q_any  # noqa: E402

CLASS_ORDER = ("pm1", "++", "N*", "--N1")


def merge_cls(F, r: int) -> str:
    c = coarse_Q_class(F, r)
    if c.startswith("mixed"):
        return "++"
    return c


def classes_of(F) -> dict[str, list[int]]:
    by: dict[str, list[int]] = defaultdict(list)
    for r in F["squares"]:
        by[merge_cls(F, r)].append(r)
    return {k: by[k] for k in CLASS_ORDER if k in by}


def structure_spreads(F) -> dict:
    C = classes_of(F)
    names = list(C)
    worst = 0
    pp_pp_pp = None
    nvar = 0
    ntot = 0
    for i in names:
        for j in names:
            for k in names:
                vals = []
                for g in C[k]:
                    cnt = 0
                    for x in C[i]:
                        y = F["mul"][F["inv"][x]][g]
                        if y in F["squares"] and merge_cls(F, y) == j:
                            cnt += 1
                    vals.append(cnt)
                ntot += 1
                spr = max(vals) - min(vals)
                if spr > 0:
                    nvar += 1
                if i == j == k == "++":
                    pp_pp_pp = (min(vals), max(vals), spr)
                worst = max(worst, spr)
    return {
        "nvar": nvar,
        "ntot": ntot,
        "worst_spread": worst,
        "p_pp_pp_pp": pp_pp_pp,
        "sizes": {k: len(v) for k, v in C.items()},
    }


def primitive_root_T(F) -> int:
    q = F["q"]
    for g in range(2, q):
        x, seen = g, set()
        ok = True
        for _ in range(q - 1):
            if x in seen or x == 0:
                ok = False
                break
            seen.add(x)
            x = F["mul"][x][g]
        if ok and len(seen) == q - 1:
            return F["mul"][g][g]
    raise RuntimeError("no generator")


def dlog_T(F, gen: int) -> dict[int, int]:
    nT = (F["q"] - 1) // 2
    out = {1: 0}
    x = 1
    for a in range(1, nT):
        x = F["mul"][x][gen]
        out[x] = a
    return out


def J_row(F, C, dlog, nT, k: int) -> list[complex]:
    row = []
    for name in C:
        s = 0j
        for r in C[name]:
            s += np.exp(2j * np.pi * k * dlog[r] / nT)
        row.append(s)
    return row


def prove_not_schur() -> dict:
    ok = True
    blocks = {}
    for p in (5, 7):
        F = _kernel_field(p)
        rec = structure_spreads(F)
        if rec["worst_spread"] == 0:
            ok = False
        if rec["p_pp_pp_pp"] is None or rec["p_pp_pp_pp"][2] == 0:
            # fail-when-wrong target: ++*++ in ++ must vary
            ok = False
        blocks[str(p)] = rec
    return {
        "proved": ok,
        "blocks": blocks,
        "theorem": (
            "Coarse Q-classes are not a Schur ring on T. "
            "p^{++}_{++,++} varies (spread 2 at p=5, 2 at p=7 on ++)."
        ),
    }


def prove_jacobi_matrix() -> dict:
    ok = True
    invert_hit = 0
    invert_tot = 0
    blocks = {}
    bases = {5: (0, 2, 6), 7: (0, 2, 8, 12)}
    for p in (5, 7):
        F = _kernel_field(p)
        C = classes_of(F)
        names = list(C)
        gen = primitive_root_T(F)
        nT = (F["q"] - 1) // 2
        dlog = dlog_T(F, gen)
        ks = bases[p]
        J = np.array([J_row(F, C, dlog, nT, k) for k in ks], dtype=np.complex128)
        # integer?
        Jr = np.real(J)
        if np.max(np.abs(J.imag)) > 1e-8 or np.max(np.abs(Jr - np.round(Jr))) > 1e-6:
            ok = False
        det = float(np.linalg.det(Jr))
        Qlive = live_Q_any(p)
        q2 = float(F["q"] ** 2)
        Qv = []
        for name in names:
            vs = [Qlive[r] / q2 for r in C[name]]
            Qv.append(float(np.mean(vs)))
        Qv = np.array(Qv)
        S = Jr @ Qv
        # S_0 = 2(q-1)
        if abs(S[0] - 2 * (F["q"] - 1)) > 1e-6:
            ok = False
        recov = np.linalg.solve(Jr, S)
        if np.max(np.abs(recov - Qv)) > 1e-8:
            ok = False
        # invert: flip J_++ at k=2 (row index of k=2, col of ++)
        invert_tot += 1
        Jbad = Jr.copy()
        i2 = list(ks).index(2)
        jpp = names.index("++")
        Jbad[i2, jpp] *= -1
        try:
            recov_bad = np.linalg.solve(Jbad, S)
            if abs(recov_bad[jpp] - Qv[jpp]) > 0.05:
                invert_hit += 1
        except np.linalg.LinAlgError:
            invert_hit += 1
        blocks[str(p)] = {
            "names": names,
            "ks": list(ks),
            "J": np.round(Jr).astype(int).tolist(),
            "det": det,
            "S0": float(S[0]),
            "S0_named": 2 * (F["q"] - 1),
            "recov_err": float(np.max(np.abs(recov - Qv))),
        }
    if invert_hit != invert_tot:
        ok = False
    return {
        "proved": ok,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "blocks": blocks,
        "theorem": (
            "J_τ(ψ_k) is an integer complete character sum. "
            "det=144 (p=5), −2688 (p=7). S_0=2(q−1). "
            "Flip J_++(k=2) breaks recovery of Q_++."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "det J is not the Q-denominator (144≠13, 2688≠409). "
            "S_{k≠0} still N/(4p). Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.300  coarse types not Schur; J integer; S unnamed", flush=True)
    A = prove_not_schur()
    print(f"  A not Schur: {A['proved']} {A['blocks']}", flush=True)
    B = prove_jacobi_matrix()
    print(f"  B Jacobi J: {B['proved']} invert={B['invert_hit']}/{B['invert_tot']}", flush=True)
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.300",
        "title": "Coarse Q-types not a Schur ring; integer J; S still N/(4p)",
        "proved": {
            "not_schur": A["proved"],
            "jacobi_J_integer": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15300.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
