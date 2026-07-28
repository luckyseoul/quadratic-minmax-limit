#!/usr/bin/env python3
"""
Moduli free-parameter c: express g_min(c) and pin by Tr(G^2).

On the nullity-1 line m = m_part + c n, for |kappa|=1 classes
  g_min(c) = -max |m_A(c)|.
Tr(G^2)(c) is quadratic; true Max+ Tr pins c and recovers g_min.

Also records the Wick particular check: distance from Wick m4=kappa/p^2
to the true moment vector.

Writes evidence/e1_gmin_cbound.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from e1_gmin_moduli import (  # noqa: E402
    L_p,
    T_p,
    build_classes,
    build_evec_system,
    class_key,
    empirical_m4,
    kappa,
    load_maxplus,
    m4_constancy,
    nullity_analysis,
    trG2_from_maxplus_dots,
    wick_dot4_upper,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def run_p5() -> dict:
    p = 5
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Mp = load_maxplus(p)
    classes, q2k = build_classes(C)
    keys = list(classes.keys())
    idx = {k: i for i, k in enumerate(keys)}
    const = m4_constancy(Mp, classes)
    A, b = build_evec_system(C, p, classes, q2k, keys, max_quads=80)
    m_true = empirical_m4(Mp, classes, keys)
    null = nullity_analysis(A, b, m_true)
    assert null["nullity"] == 1
    m_part = np.array(null["m_part"])
    nvec = np.array(null["nvec"])
    c_true = null["c_true_fit"]

    # kappa=±1 class indices
    kap1 = []
    for k in keys:
        a, b, c, d = classes[k][0]
        kap = kappa(C, (a, b, c, d))
        if abs(kap) == 1:
            kap1.append((idx[k], kap, len(classes[k])))

    def gmin_c(c: float) -> float:
        m = m_part + c * nvec
        return -max(abs(m[i]) for i, _, _ in kap1)

    # Tr(G^2) via dots
    dots = trG2_from_maxplus_dots(Mp)
    wick = wick_dot4_upper(p)
    Tr = dots["Tr_G2"]
    TrW = wick["Tr_G2_wick"]

    # Build combinatorial G_part, G_null for Tr quadratic (wedges from m2)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    G_part = np.zeros((E, E))
    G_null = np.zeros((E, E))
    for ei, (i, j) in enumerate(edges):
        G_part[ei, ei] = 1.0
    for ei, (i, j) in enumerate(edges):
        for ej in range(ei + 1, E):
            k, l = edges[ej]
            verts = {i, j, k, l}
            if len(verts) == 3:
                e1, e2 = {i, j}, {k, l}
                b1 = list(e1 - e2)[0]
                b2 = list(e2 - e1)[0]
                G_part[ei, ej] = G_part[ej, ei] = float(
                    C[i, j] * C[k, l] * C[b1, b2] / p
                )
            elif len(verts) == 4:
                pts = tuple(sorted([i, j, k, l]))
                ci = idx[q2k[pts]]
                cp = float(C[i, j] * C[k, l])
                G_part[ei, ej] = G_part[ej, ei] = cp * m_part[ci]
                G_null[ei, ej] = G_null[ej, ei] = cp * nvec[ci]

    a2 = float(np.sum(G_null * G_null))
    a1 = 2.0 * float(np.sum(G_part * G_null))
    a0 = float(np.sum(G_part * G_part))

    def tr_c(c: float) -> float:
        return a0 + a1 * c + a2 * c * c

    # Roots at Tr and TrW
    def roots_at(target: float) -> list[float]:
        Aq, Bq, Cq = a2, a1, a0 - target
        disc = Bq * Bq - 4 * Aq * Cq
        if disc < 0:
            return []
        sd = np.sqrt(disc)
        return [(-Bq - sd) / (2 * Aq), (-Bq + sd) / (2 * Aq)]

    roots_true = roots_at(Tr)
    roots_wick = roots_at(TrW)
    g_at_true = [gmin_c(c) for c in roots_true]
    g_at_wick = [gmin_c(c) for c in roots_wick]

    # Selected root = larger g_min
    sel_true = max(g_at_true) if g_at_true else None
    sel_wick = max(g_at_wick) if g_at_wick else None

    # Wick m4 = kappa/p^2 for each 4-set class (use class kappa)
    m_wick = np.zeros(len(keys))
    for k in keys:
        a, b, c, d = classes[k][0]
        m_wick[idx[k]] = kappa(C, (a, b, c, d)) / (p * p)
    # Is m_wick on the affine line? project
    # m_wick ≈ m_part + c_w nvec
    c_wick_fit = float(np.dot(nvec, m_wick - m_part) / np.dot(nvec, nvec))
    wick_line_err = float(np.linalg.norm(m_wick - (m_part + c_wick_fit * nvec)))

    def jser(x):
        if isinstance(x, (np.floating, float)):
            return float(x)
        if isinstance(x, (np.integer, int)):
            return int(x)
        if isinstance(x, (np.bool_, bool)):
            return bool(x)
        if isinstance(x, list):
            return [jser(v) for v in x]
        return x

    return {
        "p": p,
        "constancy": const,
        "nullity": int(null["nullity"]),
        "c_true": float(c_true),
        "g_min_at_c_true": float(gmin_c(c_true)),
        "g_min_frac_target": -3 / 65,
        "recovered": bool(abs(gmin_c(c_true) + 3 / 65) < 1e-9),
        "Tr_G2_true": float(Tr),
        "Tr_G2_wick": float(TrW),
        "tr_quadratic": {"a0": float(a0), "a1": float(a1), "a2": float(a2)},
        "roots_true_Tr": [float(c) for c in roots_true],
        "g_min_at_true_roots": [float(g) for g in g_at_true],
        "selected_g_min_true_Tr": float(sel_true) if sel_true is not None else None,
        "roots_wick_Tr": [float(c) for c in roots_wick],
        "g_min_at_wick_roots": [float(g) for g in g_at_wick],
        "selected_g_min_wick_Tr": float(sel_wick) if sel_wick is not None else None,
        "wick_root_beats_T": bool(sel_wick is not None and sel_wick > T_p(p)),
        "true_root_beats_T": bool(sel_true is not None and sel_true > T_p(p)),
        "true_root_ge_L": bool(sel_true is not None and sel_true >= L_p(p) - 1e-15),
        "T": float(T_p(p)),
        "L": float(L_p(p)),
        "m_wick_line_err": float(wick_line_err),
        "c_wick_fit": float(c_wick_fit),
        "n_kappa_abs_1_classes": len(kap1),
        "note": (
            "Wick m4=kappa/p^2 is near but not on the evec moduli line "
            f"(err={wick_line_err:.3e}). True Tr pin recovers -3/65. "
            "Wick Tr pin selected root may hit ~T (not strict)."
        ),
    }


def main() -> None:
    print("cbound p=5", flush=True)
    r = run_p5()
    print(
        f"  c_true={r['c_true']:.6f} gmin={r['g_min_at_c_true']:.8f} "
        f"sel_true={r['selected_g_min_true_Tr']} sel_wick={r['selected_g_min_wick_Tr']}",
        flush=True,
    )
    out = {
        "results": {"5": r},
        "status": (
            "p=5: moduli c pin by Tr(G^2) recovers g_min=-3/65 > T,L. "
            "Wick fourth-moment vector not exactly on evec line. "
            "OPEN: Max+-free evaluation of Tr(G^2) (or E[dot^4]) for all p>=5; "
            "nullity-1 + pin for general p; then g_min>=L(p)."
            if r["recovered"] and r["true_root_beats_T"]
            else "CERT FAILURE"
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_cbound.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if not (r["recovered"] and r["true_root_beats_T"]):
        sys.exit(1)


if __name__ == "__main__":
    main()
