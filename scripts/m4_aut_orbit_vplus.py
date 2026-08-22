#!/usr/bin/env python3
"""Boolean 4-point of V_+ on Aut(C)-orbits of 4-sets.

m₄(S)=E[∏_{i∈S} y_i] on Max+ ⊂ V_+ ∩ {±1}^n is Aut-invariant, hence
constant on Aut-orbits of 4-subsets of P¹.  Linear {κ,φ,star} is the
15.597 particular solution, not true m₄ (δ∈E_{4p}^{Aut}).  This script
names the orbit table at p=5 (exact) and p=7 (exact), records
⟨m₄,κ_A⟩ orbit contributions for QVAR and the binding even character,
and tests simple nonlinear formulae in (κ, has_∞, CR).

No flag flip.  Not an identity file.
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
from e1_gmin_cr_classify import canon_cr_fn  # noqa: E402
from e1_gmin_m4_prop15590 import (  # noqa: E402
    MuLab,
    paley_conference,
    signed_generators,
    signed_orbits,
)
from e1_gmin_qvar_box_master import (  # noqa: E402
    A_psi_matrix,
    make_psi,
    permutation_aut_gens,
)


def even_char_psi(p: int, k: int):
    """α_k on F_q^×, even k, 15590 field_ops via make_psi's field."""
    q = p * p
    from e1_gmin_m4_prop15590 import field_ops

    fmul, fadd, fneg, one = field_ops(p)

    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    psi = np.zeros(q, dtype=np.complex128)
    x = one
    for i in range(q - 1):
        psi[x] = np.exp(2j * np.pi * k * i / (q - 1))
        x = fmul(x, gen)
    return psi, fmul, fadd, fneg


def chunk_m4(Y: np.ndarray, S4: np.ndarray, chunk: int = 4096) -> np.ndarray:
    M = len(Y)
    nS = len(S4)
    out = np.empty(nS, dtype=np.float64)
    for a in range(0, nS, chunk):
        b = min(a + chunk, nS)
        s = S4[a:b]
        prod = (
            Y[:, s[:, 0]]
            * Y[:, s[:, 1]]
            * Y[:, s[:, 2]]
            * Y[:, s[:, 3]]
        )
        out[a:b] = prod.mean(axis=0)
    return out


def analyze(p: int) -> dict:
    C = paley_conference(p)
    n = C.shape[0]
    lab = MuLab(p, with_deg6=False)
    Y = lab.Yp.astype(np.float64)
    M = len(Y)
    print(f"\n=== p={p} n={n} |Max+|={M} ===", flush=True)
    S4 = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    nS = len(S4)
    gens = permutation_aut_gens(p, C)
    olab, sg, dead = signed_orbits(S4, gens, n, twist=False)
    uniq = np.unique(olab)
    remap = {int(u): i for i, u in enumerate(uniq)}
    loc = np.array([remap[int(x)] for x in olab], dtype=np.int64)
    nlab = len(uniq)
    sizes = np.bincount(loc, minlength=nlab)
    print(f"  4-sets={nS} Aut-orbits={nlab} dead_labels={len(dead)}", flush=True)

    m4 = chunk_m4(Y, S4)
    i, j, k, l = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    kapC = (
        C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[k, j]
    ).astype(np.int64)
    # φ, star as in 15.597
    Ci = C
    phi = (
        Ci[i, j] * Ci[i, k] * Ci[j, k]
        + Ci[i, j] * Ci[i, l] * Ci[j, l]
        + Ci[i, k] * Ci[i, l] * Ci[k, l]
        + Ci[j, k] * Ci[j, l] * Ci[k, l]
    ).astype(np.int64)
    star = (
        Ci[j, i] * Ci[k, i] * Ci[l, i]
        + Ci[i, j] * Ci[k, j] * Ci[l, j]
        + Ci[i, k] * Ci[j, k] * Ci[l, k]
        + Ci[i, l] * Ci[j, l] * Ci[k, l]
    ).astype(np.int64)
    has_inf = (S4[:, 0] == 0).astype(np.int8)  # combinations are sorted, inf=index 0

    cr_of = canon_cr_fn(p)
    # CR only needed for orbit reps
    # even characters: QVAR k=(q-1)/4 and binding (smallest λ even k)
    q = p * p
    ks = {
        "qvar": (q - 1) // 4,
        "bind": 2 if p == 5 else 8 if p == 7 else 8,
    }
    kapA = {}
    for name, kk in ks.items():
        psi, fmul, fadd, fneg = even_char_psi(p, kk)
        A = A_psi_matrix(p, C, psi, fadd, fneg)
        kapA[name] = (
            A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]
        ).astype(np.float64)

    rows = []
    pair_qvar = 0.0
    pair_bind = 0.0
    # per-orbit aggregates
    print(
        f"  {'oid':>4} {'sz':>7} {'∞':>2} {'κ':>4} {'φ':>5} {'★':>5} {'m4':>18} "
        f"{'n_m4':>4} {'⟨m,κQ⟩':>10} {'⟨m,κb⟩':>10}"
    )
    for oid in range(nlab):
        mask = loc == oid
        sz = int(sizes[oid])
        m4s = m4[mask]
        m4_mean = float(m4s.mean())
        m4_spread = float(np.ptp(m4s))
        kap_u = np.unique(kapC[mask])
        phi_u = np.unique(phi[mask])
        star_u = np.unique(star[mask])
        inf_u = np.unique(has_inf[mask])
        # exact fraction
        fr = Fraction(m4_mean).limit_denominator(M * 20)
        wq = float((m4[mask] * kapA["qvar"][mask]).sum())
        wb = float((m4[mask] * kapA["bind"][mask]).sum())
        pair_qvar += wq
        pair_bind += wb
        kap0 = int(kap_u[0]) if len(kap_u) == 1 else None
        phi0 = int(phi_u[0]) if len(phi_u) == 1 else None
        star0 = int(star_u[0]) if len(star_u) == 1 else None
        inf0 = int(inf_u[0]) if len(inf_u) == 1 else None
        # CR of first set
        S0 = tuple(int(x) for x in S4[mask][0])
        cr = cr_of(S0)
        if oid < 40 or abs(wb) > 1e-6:
            print(
                f"  {oid:4d} {sz:7d} {str(inf0):>2} {str(kap0):>4} {str(phi0):>5} "
                f"{str(star0):>5} {str(fr):>18} {len(set(np.round(m4s,10))):>4} "
                f"{wq:10.4f} {wb:10.4f} cr={cr}"
            )
        rows.append(
            {
                "oid": oid,
                "size": sz,
                "has_inf": inf0,
                "kappa": kap0,
                "phi": phi0,
                "star": star0,
                "m4": str(fr),
                "m4_float": m4_mean,
                "m4_spread": m4_spread,
                "n_m4": int(len(set(np.round(m4s, 10)))),
                "cr": cr,
                "contrib_qvar": wq,
                "contrib_bind": wb,
            }
        )

    # constancy diagnostics
    n_const = sum(1 for r in rows if r["n_m4"] == 1)
    n_kap_split = sum(1 for r in rows if r["kappa"] is None)
    print(f"  orbits with constant m4: {n_const}/{nlab}  κ not Aut-constant: {n_kap_split}")
    print(f"  ⟨m4,κ_QVAR⟩={pair_qvar:.6f}  ⟨m4,κ_bind⟩={pair_bind:.6f}")

    # linear model residual: m4 vs (κ, φ, star) on orbits with all three constant
    good = [r for r in rows if r["kappa"] is not None and r["phi"] is not None and r["star"] is not None]
    if len(good) >= 4:
        X = np.array([[1.0, r["kappa"], r["phi"], r["star"]] for r in good])
        y = np.array([r["m4_float"] for r in good])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coef
        err = float(np.max(np.abs(pred - y)))
        print(f"  lin(1,κ,φ,★) maxerr={err:.4f} coef={np.round(coef, 6)}")
        # 15.597 particular: m4p = ((p²-1)κ - 2φ - 2p star)/D , D=(n-1)(n-6)
        D = (n - 1) * (n - 6)
        part = np.array(
            [((p * p - 1) * r["kappa"] - 2 * r["phi"] - 2 * p * r["star"]) / D for r in good]
        )
        print(f"  15.597 particular vs true m4 maxerr={np.max(np.abs(part-y)):.4f}")

    # (κ, has_inf) as a key: does it determine m4?
    by_ki = defaultdict(set)
    for r in rows:
        by_ki[(r["kappa"], r["has_inf"])].add(round(r["m4_float"], 8))
    split_ki = {k: vs for k, vs in by_ki.items() if len(vs) > 1}
    print(f"  (κ,∞) keys={len(by_ki)} split={len(split_ki)}")
    by_kic = defaultdict(set)
    for r in rows:
        by_kic[(r["kappa"], r["has_inf"], r["cr"])].add(round(r["m4_float"], 8))
    split_kic = {k: vs for k, vs in by_kic.items() if len(vs) > 1}
    print(f"  (κ,∞,CR) keys={len(by_kic)} split={len(split_kic)}")
    by_kips = defaultdict(set)
    for r in rows:
        by_kips[(r["kappa"], r["has_inf"], r["phi"], r["star"])].add(round(r["m4_float"], 8))
    split_kips = {k: vs for k, vs in by_kips.items() if len(vs) > 1}
    print(f"  (κ,∞,φ,★) keys={len(by_kips)} split={len(split_kips)}")

    return {
        "p": p,
        "nS": nS,
        "n_orbits": nlab,
        "n_m4_constant": n_const,
        "pair_qvar": pair_qvar,
        "pair_bind": pair_bind,
        "n_split_kappa_inf": len(split_ki),
        "n_split_kappa_inf_cr": len(split_kic),
        "n_split_kappa_inf_phi_star": len(split_kips),
        "orbits": rows,
    }


def main():
    out = {}
    for p in (5, 7):
        out[str(p)] = analyze(p)
    path = ROOT / "evidence" / "m4_aut_orbit_vplus.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
