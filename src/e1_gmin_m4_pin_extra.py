#!/usr/bin/env python3
"""
Prop 15.82 follow-up: extra linear pins on type6+CR moduli at p=5,7.

Adds sum_S m4(S)=e4 (Prop 15.73) and denser evec rows; reports nullity.

F17: ProcessPool W≈nproc-2 via refine_moduli_multi.build_classes; atomic JSON.
Writes evidence/e1_gmin_m4_pin_extra.json
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def e4_formula(p: int) -> float:
    """Prop 15.73: ∑_S m4(S) = -p(p-1)(p+1)(p+4)/12."""
    return -p * (p - 1) * (p + 1) * (p + 4) / 12.0


def analyze(p: int, W: int) -> dict:
    from e1_gmin_m4_refine_moduli_multi import build_classes
    from e1_gmin_moduli import (
        build_evec_system,
        nullity_analysis,
        empirical_m4,
        load_maxplus,
        kappa,
    )
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(float)
    Mp = load_maxplus(p)
    print(f"  build type6+cr p={p} W={W}...", flush=True)
    classes, q2k, wall = build_classes(C, p, "type6+cr", W)
    keys = list(classes.keys())
    sizes = np.array([len(classes[k]) for k in keys], dtype=float)
    m_true = empirical_m4(Mp, classes, keys)
    emp_sum = float(np.dot(sizes, m_true))
    e4 = e4_formula(p)
    print(
        f"    ncls={len(keys)} emp_sum_m4={emp_sum:.8f} e4={e4:.8f} "
        f"match={abs(emp_sum - e4) < 1e-6}",
        flush=True,
    )

    dens = 40 if p >= 7 else 60
    dens2 = 200 if p >= 7 else 120
    A, b = build_evec_system(C, p, classes, q2k, keys, max_quads=dens)
    nul0 = nullity_analysis(A, b, m_true)
    print(f"    base nullity={nul0['nullity']} rank={nul0['rank']}", flush=True)

    A1 = np.vstack([A, sizes])
    b1 = np.append(b, e4)
    nul_e4 = nullity_analysis(A1, b1, m_true)
    print(
        f"    +e4 nullity={nul_e4['nullity']} rank={nul_e4['rank']} "
        f"res_part={nul_e4['residual_part']:.2e}",
        flush=True,
    )

    A2, b2 = build_evec_system(C, p, classes, q2k, keys, max_quads=dens2)
    nul_dense = nullity_analysis(A2, b2, m_true)
    print(
        f"    dense nullity={nul_dense['nullity']} rank={nul_dense['rank']}",
        flush=True,
    )

    A3 = np.vstack([A2, sizes])
    b3 = np.append(b2, e4)
    nul_both = nullity_analysis(A3, b3, m_true)
    print(
        f"    dense+e4 nullity={nul_both['nullity']} rank={nul_both['rank']}",
        flush=True,
    )

    # kappa-stratum mean constraints: for each kap value, average m on kap-class
    # is free still — try sum_{κ=±1} size m  and sum_{κ=±3} size m separately
    kappas = np.array([float(kappa(C, classes[k][0])) for k in keys])
    out_rows = []
    out_rhs = []
    for target in (1.0, -1.0, 3.0, -3.0):
        mask = np.abs(kappas - target) < 0.5
        if not np.any(mask):
            continue
        row = sizes * mask.astype(float)
        # empirical sum on that stratum
        rhs = float(np.dot(sizes[mask], m_true[mask]))
        out_rows.append(row)
        out_rhs.append(rhs)
    if out_rows:
        A4 = np.vstack([A2] + out_rows)
        b4 = np.concatenate([b2, np.array(out_rhs)])
        # these are tautological for true m — only constrain the affine space
        # if we use closed forms for stratum sums. Skip as pins of true.
        # Instead: pure C pin sum size*kappa is independent of m
        pass

    # sum kappa is C-only; linear pin on m: none from that.
    # Character of Aut: try mean m on |κ|=1 only equals known?
    # Use sum_{|κ|=1} size * m  from e4 + sum_κ m = known if we know both.
    # ∑ m4 = e4 and ∑ κ known; not enough for each stratum.

    return {
        "p": int(p),
        "mode": "type6+cr",
        "n_classes": len(keys),
        "wall_class_s": float(wall),
        "emp_sum_m4": emp_sum,
        "e4": e4,
        "sum_m4_matches_e4": abs(emp_sum - e4) < 1e-6,
        "base": {
            "nullity": nul0["nullity"],
            "rank": nul0["rank"],
            "residual_part": nul0["residual_part"],
        },
        "plus_e4": {
            "nullity": nul_e4["nullity"],
            "rank": nul_e4["rank"],
            "residual_part": nul_e4["residual_part"],
        },
        "dense": {
            "nullity": nul_dense["nullity"],
            "rank": nul_dense["rank"],
            "residual_part": nul_dense["residual_part"],
        },
        "dense_plus_e4": {
            "nullity": nul_both["nullity"],
            "rank": nul_both["rank"],
            "residual_part": nul_both["residual_part"],
        },
        "wall_s": float(time.perf_counter() - t0),
        "workers": int(W),
    }


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    print(f"pin_extra Prop15.82: W={W}", flush=True)
    out = {
        "prop": "15.82",
        "workers": W,
        "compute": snap,
        "io_policy": "mmap Max+; ProcessPool type6+CR; write_json_atomic",
        "results": {},
        "status": None,
    }
    for p in (5, 7):
        print(f"=== p={p} ===", flush=True)
        out["results"][str(p)] = analyze(p, W)

    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["status"] = (
        f"Extra pins type6+CR: p5 nul base={r5['base']['nullity']} "
        f"+e4={r5['plus_e4']['nullity']} dense+e4={r5['dense_plus_e4']['nullity']}; "
        f"p7 nul base={r7['base']['nullity']} +e4={r7['plus_e4']['nullity']} "
        f"dense+e4={r7['dense_plus_e4']['nullity']}. "
        f"OPEN: p7 multi-param or general pin; lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_pin_extra.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
