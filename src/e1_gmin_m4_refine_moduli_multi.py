#!/usr/bin/env python3
"""
Prop 15.82 multi-strat moduli: compare constant stratifications for nullity
and GD/cand pin at p=5,7.

F17: ProcessPool W≈nproc-2 for class build; mmap Max+; atomic JSON.
Modes (all m4-constant at p=5,7 per e1_gmin_m4_refine):
  type6+cr, coarse+cr, coarse+cr+kappa

If nullity=1: line pin c* safe-side of c_GD.
If nullity≥2: multi-param — project true Max+ into particular+null; record
  true GD/cand; sample a small grid on null coords around true for GD/cand.

Writes evidence/e1_gmin_m4_refine_moduli_multi.json
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


_G: dict = {}


def _init_build(C: np.ndarray, p: int, mode: str) -> None:
    _G["C"] = C
    _G["p"] = p
    _G["mode"] = mode
    _blas1()
    from e1_gmin_moduli import type6, ext_sum_hist, kappa
    from e1_gmin_cr_classify import canon_cr_fn

    _G["type6"] = type6
    _G["ext_sum_hist"] = ext_sum_hist
    _G["kappa"] = kappa
    _G["cr_fn"] = canon_cr_fn(p)


def _key_of(pts: tuple) -> tuple:
    mode = _G["mode"]
    C = _G["C"]
    t6 = _G["type6"](C, pts)
    ext = _G["ext_sum_hist"](C, pts)
    cr = _G["cr_fn"](pts)
    if mode == "type6+cr":
        return (t6, cr)
    if mode == "coarse+cr":
        return (t6, ext, cr)
    if mode == "coarse+cr+kappa":
        return (t6, ext, cr, _G["kappa"](C, pts))
    raise ValueError(mode)


def _shard_build(quads_slice: list) -> list:
    buckets: dict = defaultdict(list)
    for pts in quads_slice:
        pts = tuple(pts)
        buckets[_key_of(pts)].append(pts)
    return [(k, v) for k, v in buckets.items()]


def build_classes(C: np.ndarray, p: int, mode: str, W: int):
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
        initargs=(C, p, mode),
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


def max_star_S1_on_m(C, p, classes, keys, m_vec, quad_to_key, sample=5):
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
        "n_positive": int(n_pos),
        "n_checks": int(n_chk),
    }


def nullspace_basis(A: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis of ker(A) as columns (n_var, nullity)."""
    u, s, vt = np.linalg.svd(A, full_matrices=True)
    rank = int(np.sum(s > tol * max(1.0, s[0] if len(s) else 1.0)))
    return vt[rank:].T  # (n_var, nullity)


def analyze_mode(p: int, mode: str, W: int) -> dict:
    from e1_gmin_moduli import (
        build_evec_system,
        nullity_analysis,
        empirical_m4,
        m4_constancy,
        load_maxplus,
        kappa,
    )
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(float)
    Mp = load_maxplus(p)
    print(f"  [{mode}] build classes p={p} W={W}...", flush=True)
    classes, quad_to_key, wall_cls = build_classes(C, p, mode, W)
    keys = list(classes.keys())
    const = m4_constancy(Mp, classes, sample=40)
    print(
        f"    n_classes={len(keys)} all_const={const['all_constant']} "
        f"({const['n_constant_m4']}/{const['n_classes']}) "
        f"max_std={const['max_std']:.2e} build={wall_cls:.2f}s",
        flush=True,
    )
    out: dict = {
        "p": int(p),
        "mode": mode,
        "n_classes": len(keys),
        "m4_all_constant": bool(const["all_constant"]),
        "constancy": {
            "n_classes": const["n_classes"],
            "n_constant_m4": const["n_constant_m4"],
            "all_constant": const["all_constant"],
            "max_std": const["max_std"],
        },
        "wall_class_s": float(wall_cls),
        "M_cand": float(M_cand(p)),
        "workers": int(W),
    }
    if not const["all_constant"]:
        out["status"] = "non-constant — skip moduli"
        out["wall_s"] = float(time.perf_counter() - t0)
        return out

    print(f"  [{mode}] evec system p={p} n_cls={len(keys)}...", flush=True)
    A, b = build_evec_system(
        C, p, classes, quad_to_key, keys, max_quads=40 if p >= 7 else 60
    )
    m_true = empirical_m4(Mp, classes, keys)
    nul = nullity_analysis(A, b, m_true)
    print(
        f"    rank={nul['rank']} nullity={nul['nullity']} "
        f"fit_err={nul.get('fit_err')}",
        flush=True,
    )
    out["nullity"] = nul["nullity"]
    out["rank"] = nul["rank"]
    out["n_var"] = nul["n_var"]
    out["residual_part"] = nul["residual_part"]

    r_true = max_star_S1_on_m(C, p, classes, keys, m_true, quad_to_key)
    out["at_true_Max+"] = {
        **r_true,
        "m4_le_cand": r_true["max_abs_m4_kappa1"] <= M_cand(p) + 1e-12,
    }
    print(
        f"    true: max|m4|={r_true['max_abs_m4_kappa1']:.10f} "
        f"GD={r_true['GD']} max_ss={r_true['max_star_S1']:.8f}",
        flush=True,
    )

    if nul["nullity"] == 0:
        out["moduli_GD_ok"] = bool(
            r_true["GD"] and out["at_true_Max+"]["m4_le_cand"]
        )
        out["note"] = "unique solution"
    elif nul["nullity"] == 1:
        m_part = np.array(nul["m_part"])
        nvec = np.array(nul["nvec"])
        c_true = float(nul["c_true_fit"])
        out["c_true"] = c_true
        out["fit_err"] = float(nul["fit_err"])
        cs = np.linspace(c_true - 0.4, c_true + 0.4, 9)
        rows = []
        for c in cs:
            m = m_part + c * nvec
            r = max_star_S1_on_m(
                C, p, classes, keys, m, quad_to_key, sample=4
            )
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
        if c_GD is not None and abs(beta) > 1e-15:
            if beta > 0:
                c_on_safe = c_true <= c_GD + 1e-12
            else:
                c_on_safe = c_true >= c_GD - 1e-12
        else:
            c_on_safe = False
        out["affine"] = {
            "alpha": alpha,
            "beta": beta,
            "fit_maxerr": fit_err,
            "exactly_affine": fit_err < 1e-10,
            "beta_positive": beta > 0,
        }
        out["c_GD"] = c_GD
        out["c_true_on_safe_side_of_c_GD"] = c_on_safe
        out["line_sample"] = rows
        out["moduli_GD_ok"] = bool(
            c_on_safe
            and r_true["GD"]
            and out["at_true_Max+"]["m4_le_cand"]
            and out["affine"]["exactly_affine"]
        )
        print(
            f"    line: c*={c_true:.6f} c_GD={c_GD} "
            f"safe={c_on_safe} moduli_GD_ok={out['moduli_GD_ok']}",
            flush=True,
        )
    else:
        # multi-param: nullspace basis + true coords; sample grid near true
        N = nullspace_basis(A)
        nul_dim = N.shape[1]
        # particular: least-squares
        m_part, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        # coords of true residual in null basis
        res = m_true - m_part
        # N may not be orthonormal after floating; use least-squares coords
        c_true_vec, _, _, _ = np.linalg.lstsq(N, res, rcond=None)
        recon_err = float(np.linalg.norm(m_part + N @ c_true_vec - m_true))
        out["multi"] = {
            "nullity": int(nul_dim),
            "c_true": [float(x) for x in c_true_vec],
            "recon_err": recon_err,
        }
        # grid around true in first 2 null coords (or 1 if only 1 extra)
        kdim = min(2, nul_dim)
        deltas = [-0.15, -0.05, 0.0, 0.05, 0.15]
        samples = []
        worst_abs = 0.0
        any_not_gd = False
        any_over_cand = False
        n_samp = 0
        # product grid on first kdim coords; fix remaining at true
        from itertools import product

        for offs in product(deltas, repeat=kdim):
            c = c_true_vec.copy()
            for i, o in enumerate(offs):
                c[i] = c_true_vec[i] + o
            m = m_part + N @ c
            r = max_star_S1_on_m(
                C, p, classes, keys, m, quad_to_key, sample=3
            )
            n_samp += 1
            worst_abs = max(worst_abs, r["max_abs_m4_kappa1"])
            if not r["GD"]:
                any_not_gd = True
            if r["max_abs_m4_kappa1"] > M_cand(p) + 1e-10:
                any_over_cand = True
            if len(samples) < 12:
                samples.append(
                    {
                        "offs": [float(x) for x in offs],
                        "max_star_S1": r["max_star_S1"],
                        "max_abs_m4": r["max_abs_m4_kappa1"],
                        "GD": r["GD"],
                    }
                )
        out["multi"]["grid"] = {
            "kdim": kdim,
            "deltas": deltas,
            "n_samples": n_samp,
            "worst_max_abs_m4": float(worst_abs),
            "all_GD": not any_not_gd,
            "all_le_cand": not any_over_cand,
            "samples_head": samples,
        }
        # moduli_GD_ok only for true Max+ here; multi-param full pin OPEN
        # unless grid all GD+cand (local evidence, not global proof)
        out["moduli_GD_ok"] = bool(
            r_true["GD"] and out["at_true_Max+"]["m4_le_cand"]
        )
        out["local_multi_pin_ok"] = bool(
            out["multi"]["grid"]["all_GD"]
            and out["multi"]["grid"]["all_le_cand"]
        )
        out["note"] = (
            f"nullity={nul_dim}: multi-param; true GD/cand; "
            f"local grid pin={out['local_multi_pin_ok']}; full pin OPEN"
        )
        print(
            f"    multi nul={nul_dim} recon={recon_err:.2e} "
            f"local_pin={out['local_multi_pin_ok']} "
            f"worst|m4|={worst_abs:.6f}",
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
    modes = ["type6+cr", "coarse+cr", "coarse+cr+kappa"]
    print(
        f"refine_moduli_multi Prop15.82: W={W} modes={modes} "
        f"gpu={snap.get('gpu', {}).get('available')}",
        flush=True,
    )
    out: dict = {
        "prop": "15.82",
        "workers": W,
        "compute": snap,
        "modes": modes,
        "io_policy": (
            "mmap Max+; ProcessPool multi-strat class build; "
            "evec moduli; write_json_atomic"
        ),
        "results": {},
        "status": None,
    }
    for p in (5, 7):
        out["results"][str(p)] = {}
        print(f"=== p={p} ===", flush=True)
        for mode in modes:
            out["results"][str(p)][mode] = analyze_mode(p, mode, W)

    # summary
    summary = {}
    for p in ("5", "7"):
        summary[p] = {}
        for mode, r in out["results"][p].items():
            summary[p][mode] = {
                "n_classes": r.get("n_classes"),
                "constant": r.get("m4_all_constant"),
                "nullity": r.get("nullity"),
                "moduli_GD_ok": r.get("moduli_GD_ok"),
                "local_multi_pin_ok": r.get("local_multi_pin_ok"),
                "c_true_on_safe_side": r.get("c_true_on_safe_side_of_c_GD"),
                "true_le_cand": (r.get("at_true_Max+") or {}).get("m4_le_cand"),
                "true_GD": (r.get("at_true_Max+") or {}).get("GD"),
            }
    out["summary"] = summary
    # best: prefer nullity 0/1 with moduli_GD_ok
    best = {}
    for p in ("5", "7"):
        cands = []
        for mode, r in out["results"][p].items():
            if not r.get("m4_all_constant"):
                continue
            nul = r.get("nullity", 99)
            ok = bool(r.get("moduli_GD_ok"))
            score = (0 if nul <= 1 and ok else 1, nul, -int(ok))
            cands.append((score, mode, r))
        cands.sort(key=lambda x: x[0])
        if cands:
            _, mode, r = cands[0]
            best[p] = {
                "mode": mode,
                "nullity": r.get("nullity"),
                "moduli_GD_ok": r.get("moduli_GD_ok"),
                "n_classes": r.get("n_classes"),
            }
    out["best"] = best
    out["status"] = (
        f"Prop 15.82 multi-strat: best p5={best.get('5')} p7={best.get('7')}. "
        f"OPEN: general p≥5 pin (or closed multi-param); lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_refine_moduli_multi.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
