#!/usr/bin/env python3
"""
Prop 15.82 (moduli): type6+CR stratification → constant m4 at p=5,7;
evec system + GD/cand on the line (or unique solution).

F17: GPU m4 if needed; ProcessPool for heavy class work; atomic JSON.
Uses refined key (type6, CR) from Prop 15.82 discovery.

Writes evidence/e1_gmin_m4_refine_moduli.json
"""
from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
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


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


# ── refined class key: type6 + CR ───────────────────────────────────────────

_G: dict = {}


def _init_build(C: np.ndarray, p: int) -> None:
    _G["C"] = C
    _G["p"] = p
    _blas1()
    from e1_gmin_moduli import type6
    from e1_gmin_cr_classify import canon_cr_fn

    _G["type6"] = type6
    _G["cr_fn"] = canon_cr_fn(p)


def _shard_build(quads_slice: list) -> dict:
    C = _G["C"]
    t6 = _G["type6"]
    cr = _G["cr_fn"]
    buckets: dict = defaultdict(list)
    for pts in quads_slice:
        pts = tuple(pts)
        key = (t6(C, pts), cr(pts))
        buckets[key].append(pts)
    # return as list of (key, quads) — keys must be picklable
    return [(k, v) for k, v in buckets.items()]


def build_refined_classes(C: np.ndarray, p: int, W: int) -> tuple[dict, dict]:
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nsh = min(W, max(1, len(quads) // 400))
    ss = (len(quads) + nsh - 1) // nsh
    shards = [quads[i : i + ss] for i in range(0, len(quads), ss)]
    classes: dict = defaultdict(list)
    t0 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_build,
        initargs=(C, p),
    ) as ex:
        futs = [ex.submit(_shard_build, sh) for sh in shards]
        for fut in as_completed(futs):
            for k, vs in fut.result():
                classes[k].extend(vs)
    wall = time.perf_counter() - t0
    quad_to_key = {}
    for k, qs in classes.items():
        for q in qs:
            quad_to_key[q] = k
    return dict(classes), quad_to_key, wall


def star_of(C, S, v) -> int:
    pr = 1
    for u in S:
        if u != v:
            pr *= int(C[v, u])
    return pr


def max_star_S1_on_m(C, p, classes, keys, m_vec, quad_to_key, sample=6):
    from e1_gmin_moduli import kappa

    idx = {k: i for i, k in enumerate(keys)}
    inv = 1.0 / (p * p)
    n = C.shape[0]
    max_ss = -1e30
    min_ss = 1e30
    max_abs = 0.0
    n_pos = 0
    n_chk = 0
    for kA, quads in classes.items():
        if abs(kappa(C, quads[0])) != 1:
            continue
        max_abs = max(max_abs, abs(float(m_vec[idx[kA]])))
        for S in quads[:sample]:
            S = tuple(S)
            v = S[0]
            st = star_of(C, S, v)
            others = tuple(x for x in S if x != v)
            S1 = 0.0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                kp = kappa(C, Sp)
                if abs(kp) != 1:
                    continue
                m4p = float(m_vec[idx[quad_to_key[Sp]]])
                S1 += float(C[v, rr]) * (m4p - kp * inv)
            ss = st * S1
            max_ss = max(max_ss, ss)
            min_ss = min(min_ss, ss)
            if ss > 1e-12:
                n_pos += 1
            n_chk += 1
    return {
        "max_star_S1": float(max_ss),
        "min_star_S1": float(min_ss),
        "max_abs_m4_kappa1": float(max_abs),
        "GD": n_pos == 0 and max_ss <= 1e-12,
        "n_positive": n_pos,
        "n_checks": n_chk,
    }


def analyze_p(p: int, W: int) -> dict:
    from e1_gmin_moduli import (
        build_evec_system,
        nullity_analysis,
        empirical_m4,
        m4_constancy,
        load_maxplus,
        kappa,
        trG2_from_maxplus_dots,
    )
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(float)
    Mp = load_maxplus(p)
    print(f"  build type6+CR classes p={p} W={W}...", flush=True)
    classes, quad_to_key, wall_cls = build_refined_classes(C, p, W)
    keys = list(classes.keys())
    const = m4_constancy(Mp, classes, sample=50)
    print(
        f"    n_classes={len(keys)} constant={const['all_constant']} "
        f"({const['n_constant_m4']}/{const['n_classes']}) "
        f"max_std={const['max_std']:.2e} build={wall_cls:.2f}s",
        flush=True,
    )
    if not const["all_constant"]:
        return {
            "p": p,
            "n_classes": len(keys),
            "m4_all_constant": False,
            "constancy": const,
            "wall_s": float(time.perf_counter() - t0),
            "status": "type6+CR failed constancy — abort moduli",
        }

    print(f"  evec system p={p}...", flush=True)
    A, b = build_evec_system(
        C, p, classes, quad_to_key, keys, max_quads=50 if p >= 7 else 60
    )
    m_true = empirical_m4(Mp, classes, keys)
    nul = nullity_analysis(A, b, m_true)
    print(
        f"    rank={nul['rank']} nullity={nul['nullity']} "
        f"fit_err={nul.get('fit_err')}",
        flush=True,
    )

    out: dict = {
        "p": int(p),
        "stratification": "type6+CR",
        "n_classes": len(keys),
        "m4_all_constant": True,
        "constancy": const,
        "wall_class_s": float(wall_cls),
        "nullity": nul["nullity"],
        "rank": nul["rank"],
        "n_var": nul["n_var"],
        "residual_part": nul["residual_part"],
        "M_cand": float(M_cand(p)),
        "wall_s": float(time.perf_counter() - t0),
        "io": "mmap Max+; ProcessPool type6+CR build; write_json_atomic",
        "workers": int(W),
    }

    r_true = max_star_S1_on_m(C, p, classes, keys, m_true, quad_to_key)
    out["at_true_Max+"] = {
        **r_true,
        "m4_le_cand": r_true["max_abs_m4_kappa1"] <= M_cand(p) + 1e-12,
    }
    print(
        f"    true Max+: max|m4|={r_true['max_abs_m4_kappa1']:.10f} "
        f"M_cand={M_cand(p):.10f} GD={r_true['GD']} "
        f"max_star_S1={r_true['max_star_S1']:.8f}",
        flush=True,
    )

    if nul["nullity"] == 1:
        m_part = np.array(nul["m_part"])
        nvec = np.array(nul["nvec"])
        c_true = float(nul["c_true_fit"])
        out["c_true"] = c_true
        out["fit_err"] = float(nul["fit_err"])
        # sample line for affine max star·S1
        cs = np.linspace(c_true - 0.4, c_true + 0.4, 9)
        rows = []
        for c in cs:
            m = m_part + c * nvec
            r = max_star_S1_on_m(C, p, classes, keys, m, quad_to_key, sample=4)
            rows.append(
                {
                    "c": float(c),
                    "max_star_S1": r["max_star_S1"],
                    "max_abs_m4": r["max_abs_m4_kappa1"],
                    "GD": r["GD"],
                }
            )
        X = np.array([[1.0, r["c"]] for r in rows])
        y = np.array([r["max_star_S1"] for r in rows])
        coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        alpha, beta = float(coef[0]), float(coef[1])
        fit_err = float(np.max(np.abs(X @ coef - y)))
        c_GD = -alpha / beta if abs(beta) > 1e-15 else None
        # max star·S1 = α+βc ≤ 0  ⇔  c ≤ c_GD if β>0,  c ≥ c_GD if β<0
        if c_GD is not None and abs(beta) > 1e-15:
            if beta > 0:
                c_on_safe_side = c_true <= c_GD + 1e-12
            else:
                c_on_safe_side = c_true >= c_GD - 1e-12
        else:
            c_on_safe_side = False
        pred_at_true = alpha + beta * c_true
        out["affine"] = {
            "alpha": alpha,
            "beta": beta,
            "fit_maxerr": fit_err,
            "exactly_affine": fit_err < 1e-10,
            "beta_positive": beta > 0,
            "pred_max_star_S1_at_c_true": pred_at_true,
            "pred_matches_true_max": abs(pred_at_true - r_true["max_star_S1"])
            < 1e-9,
        }
        out["c_GD"] = c_GD
        out["c_true_on_safe_side_of_c_GD"] = c_on_safe_side
        out["line_sample"] = rows
        out["moduli_GD_ok"] = bool(
            c_on_safe_side
            and r_true["GD"]
            and out["at_true_Max+"]["m4_le_cand"]
            and out["affine"]["exactly_affine"]
        )
        print(
            f"    nullity1: c*={c_true:.8f} c_GD={c_GD} β={beta:.6f} "
            f"safe_side={c_on_safe_side} pred_S1={pred_at_true:.8f} "
            f"affine_err={fit_err:.2e} moduli_GD_ok={out['moduli_GD_ok']}",
            flush=True,
        )
        tr = trG2_from_maxplus_dots(Mp)
        out["Tr_G2"] = tr["Tr_G2"]
    elif nul["nullity"] == 0:
        m_sol = np.linalg.lstsq(A, b, rcond=1e-10)[0]
        err = float(np.linalg.norm(m_sol - m_true))
        r_sol = max_star_S1_on_m(C, p, classes, keys, m_sol, quad_to_key)
        out["unique_solution"] = {
            "||m_sol-m_true||": err,
            "matches_true": err < 1e-8,
            "max_abs_m4_sol": r_sol["max_abs_m4_kappa1"],
            "GD_sol": r_sol["GD"],
            "m4_le_cand_sol": r_sol["max_abs_m4_kappa1"] <= M_cand(p) + 1e-12,
        }
        out["moduli_GD_ok"] = bool(
            out["unique_solution"]["matches_true"]
            and r_true["GD"]
            and out["at_true_Max+"]["m4_le_cand"]
        )
        print(
            f"    nullity0: ||sol-true||={err:.3e} "
            f"moduli_GD_ok={out['moduli_GD_ok']}",
            flush=True,
        )
    else:
        # nullity ≥ 2: true Max+ GD/cand still load-bearing; multi-param pin OPEN
        out["moduli_GD_ok"] = bool(
            r_true["GD"] and out["at_true_Max+"]["m4_le_cand"]
        )
        out["note"] = (
            f"nullity={nul['nullity']}: multi-parameter moduli; "
            f"true Max+ GD/cand recorded; full c-vector pin OPEN"
        )
        print(
            f"    nullity{nul['nullity']}: true GD/cand; "
            f"moduli_GD_ok(true)={out['moduli_GD_ok']}",
            flush=True,
        )

    out["wall_s"] = float(time.perf_counter() - t0)
    return out


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    print(
        f"refine_moduli Prop15.82: W={W} gpu={snap['gpu'].get('available')}",
        flush=True,
    )
    out: dict = {
        "prop": "15.82",
        "workers": W,
        "compute": snap,
        "stratification": "type6+CR (PGL-complete CR + S4 edge type)",
        "io_policy": (
            "mmap Max+; ProcessPool type6+CR class build; "
            "evec moduli; write_json_atomic"
        ),
        "results": {},
        "status": None,
    }
    for p in (5, 7):
        print(f"=== p={p} ===", flush=True)
        out["results"][str(p)] = analyze_p(p, W)

    r5, r7 = out["results"]["5"], out["results"]["7"]
    both_const = r5.get("m4_all_constant") and r7.get("m4_all_constant")
    both_ok = r5.get("moduli_GD_ok") and r7.get("moduli_GD_ok")
    out["both_m4_constant"] = both_const
    out["both_moduli_GD_ok"] = both_ok
    out["status"] = (
        f"Prop 15.82: type6+CR stratification — m4 constant p=5 "
        f"({r5.get('n_classes')} cls) and p=7 ({r7.get('n_classes')} cls) "
        f"=> {both_const}. "
        f"Moduli GD+cand: p5={r5.get('moduli_GD_ok')} "
        f"(nullity={r5.get('nullity')}), "
        f"p7={r7.get('moduli_GD_ok')} (nullity={r7.get('nullity')}). "
        f"OPEN: prove for all primes p≥5 (or closed c*≤c_GD); lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_refine_moduli.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
