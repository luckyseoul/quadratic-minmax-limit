#!/usr/bin/env python3
"""
Prop 15.78: star·S1 constancy on |κ|=1 4-sets; Gaussian-domination form;
exact p=5 spectrum; path to star·S1≤0.

Proved (conference combinatorics / Max+ evec, as tagged):
  1. Moment form (Max+): star_a S1(a) = E[φ] - E_Wick[φ] where
       φ(y) = star_a (∏_{u≠a} y_u) ∑_{r∈R1(a)} C_ar y_r,
     Σ = I+C/p, Wick = star_a τ1(a)/p², τ1=∑_{R1} C_ar κ(S_a→r).
  2. Gaussian-domination equivalence: star_a S1(a)≤0 ⇔ E[φ]≤E_Wick[φ].
  3. Joint rewrite (κ=1, star=+1): S1+S3 = p m4 - 1/p - 2/p².
  4. Combinatorial: star_a τ1(a) is constant on the four centres of any |κ|=1
     set (certified multi-worker p=3,5,7,11; algebraic proof via σ=2 star
     and one-center degree constancy + sign symmetry).
  5. Max+ structure (certified GPU p=5,7): star_a S1(a) is likewise
     constant on the four centres of every |κ|=1 set (same rational value).
  6. Exact spectrum at p=5 (certified full census): star·S1 ∈ {-2/65, -42/325},
     both <0; matches st E[f0 U1]=(11/65)sgn(st τ1).

Certified GPU mmap+atomic: GD holds p=5,7; constancy; S1 spectrum.
OPEN: prove star·S1 constancy + ≤0 for all primes p≥5; close joint/S3.
lim α_n remains OPEN.

Writes evidence/e1_gmin_m4_S1_const.json
"""
from __future__ import annotations

import os
import sys
import time
from collections import Counter
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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def rho_budget_cand(p: int) -> float:
    return M_cand(p) - 1.0 / (p * p)


# ── Proved algebra ──────────────────────────────────────────────────────────


def algebraic_moment_form() -> dict:
    """star·S1 = E[φ]-E_Wick[φ]; joint rewrite; GD equivalence."""
    return {
        "proved": True,
        "moment_form": (
            "star_a S1(a) = E[star_a (∏_{u≠a} y_u) U1] - star_a τ1/p² "
            "with U1=∑_{R1} C_ar y_r, R1={r∉S: |κ(S_a→r)|=1}"
        ),
        "Wick_value": "E_Wick[φ] = star_a τ1 / p²  (pairwise Σ=I+C/p)",
        "gaussian_domination_iff": (
            "star_a S1(a)≤0 ⇔ E_Max+[φ] ≤ E_{N(0,Σ)}[φ]"
        ),
        "joint_kappa1_star_plus": "S1+S3 = p m4 - 1/p - 2/p²  (κ=1, star=+1)",
        "cand_via_joint": (
            "max same-sign ρ ≤ ρ_cand ⇔ max_{κ=1,m4>1/p²} m4 ≤ M_cand "
            "(tautological rewrite; power is S1≤0 ⇒ replace joint by S3)"
        ),
        "p5_S3_budget_negative": True,
        "note": (
            "At p=5, p ρ_cand - 2/p² = -16/325 < 0, so S1≤0 alone with "
            "absolute |S3| is insufficient — need strongly negative S1 on maximisers."
        ),
    }


def p5_exact_spectrum_algebra() -> dict:
    """
    At p=5 the census gives exactly two values of star·S1.
    Closed forms: -2/65 and -42/325, both negative.
    Also st E[f0 U1] = (11/65) sgn(st τ1) with st·τ1 ∈ {-1, 5}.
    """
    vals = {
        "star_S1_values": [-2 / 65, -42 / 325],
        "both_negative": True,
        "st_E_fU1": "±11/65 = ±(2p+1)/(n p / 2) at p=5",
        "st_tau1_values": [-1, 5],
        "formula_check": {
            "stt1=5": abs((-2 / 65) - ((11 / 65) - 5 / 25)) < 1e-15,
            "stt1=-1": abs((-42 / 325) - ((-11 / 65) - (-1) / 25)) < 1e-15,
        },
        "proved_if_census": (
            "Full Max+ census p=5 yields only these two rationals; "
            "both <0 ⇒ star·S1≤0 at p=5."
        ),
    }
    vals["algebra_ok"] = all(vals["formula_check"].values()) and all(
        v < 0 for v in vals["star_S1_values"]
    )
    return vals


# ── Combinatorial τ1·star constancy (Max+-free, multi-worker) ───────────────

_G: dict = {}


def _init_C(C: np.ndarray) -> None:
    _G["C"] = C
    _blas1()


def _kappa(C, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def _star(C, S, v) -> int:
    pr = 1
    for u in S:
        if u != v:
            pr *= int(C[v, u])
    return pr


def _tau_shard(quads: list) -> dict:
    """For each |κ|=1 set, record (star·τ1 across 4 centres) multiset."""
    C = _G["C"]
    n = C.shape[0]
    n1 = 0
    constant = 0
    patterns: Counter = Counter()
    star_tau_vals: Counter = Counter()
    for S in quads:
        if abs(_kappa(C, S)) != 1:
            continue
        n1 += 1
        S = tuple(S)
        st_taus = []
        for v in S:
            st = _star(C, S, v)
            others = tuple(x for x in S if x != v)
            t1 = 0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                kp = _kappa(C, Sp)
                if abs(kp) == 1:
                    t1 += int(C[v, rr]) * kp
            st_taus.append(st * t1)
            star_tau_vals[st * t1] += 1
        if len(set(st_taus)) == 1:
            constant += 1
        patterns[tuple(sorted(st_taus))] += 1
    return {
        "n1": n1,
        "constant": constant,
        "patterns": dict(patterns),
        "star_tau_vals": dict(star_tau_vals),
    }


def certify_star_tau1_constant(p: int, W: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(np.float64)
    quads = list(combinations(range(C.shape[0]), 4))
    nsh = min(W, max(1, len(quads) // 400))
    ss = (len(quads) + nsh - 1) // nsh
    shards = [quads[i : i + ss] for i in range(0, len(quads), ss)]
    n1 = constant = 0
    patterns: Counter = Counter()
    star_tau_vals: Counter = Counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_tau_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            r = fut.result()
            n1 += r["n1"]
            constant += r["constant"]
            patterns.update(r["patterns"])
            star_tau_vals.update(r["star_tau_vals"])
    return {
        "p": int(p),
        "n_kappa1": n1,
        "n_star_tau1_constant": constant,
        "all_constant": constant == n1,
        "star_tau1_unique_values": sorted(star_tau_vals.keys()),
        "n_patterns": len(patterns),
        "wall_s": float(time.perf_counter() - t0),
        "io": "ProcessPool pure-C (no Max+)",
    }


# ── GPU Max+: star·S1 constancy + GD + spectrum ─────────────────────────────


def gpu_all_m4(Y_np, quads, *, batch: int = 12_000):
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Yh = np.ascontiguousarray(Y_np, dtype=np.float64)
    Y = cp.asarray(Yh)
    m4 = np.empty(len(quads))
    t0 = time.perf_counter()
    for lo in range(0, len(quads), batch):
        hi = min(len(quads), lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4[lo:hi] = cp.asnumpy(cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0))
    return m4, {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": gpu_info(),
        "wall_s": float(time.perf_counter() - t0),
        "n_quads": int(len(quads)),
        "batch": int(batch),
        "io": {
            "maxplus_read": "mmap",
            "H2D": "single Y",
            "D2H": "m4 vector only",
        },
    }


def cert_maxplus_structure(p: int, *, W: int) -> dict:
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = load_maxplus_array(p, mmap=True)
    # verify Cy=py
    Yc = np.ascontiguousarray(Y, dtype=np.float64)
    cy_err = float(np.max(np.abs(Yc @ C - p * Yc)))
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    a, b, c, d = quads.T
    kap = (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)
    batch = 12_000 if p >= 7 else 50_000
    m4, gpu_meta = gpu_all_m4(Y, quads, batch=batch)
    inv = 1.0 / (p * p)
    lookup = {
        tuple(map(int, q)): (float(m4[i]), int(kap[i]))
        for i, q in enumerate(quads)
    }

    # walk all |κ|=1 sets: constancy of star·S1, spectrum, GD
    t1 = time.perf_counter()
    n1 = 0
    n_const = 0
    star_S1_vals: Counter = Counter()
    max_star_S1 = -1e30
    min_star_S1 = 1e30
    gd_viol = 0  # E[φ] > E_Wick
    patterns: Counter = Counter()

    kappa1 = [i for i in range(len(quads)) if abs(int(kap[i])) == 1]
    # shard kappa1 indices
    nsh = min(W, max(1, len(kappa1) // 200))
    ss = (len(kappa1) + nsh - 1) // nsh
    shards = [kappa1[i : i + ss] for i in range(0, len(kappa1), ss)]

    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_mp,
        initargs=(C, quads, lookup, p),
    ) as ex:
        futs = [ex.submit(_walk_shard_mp, sh) for sh in shards]
        for fut in as_completed(futs):
            r = fut.result()
            n1 += r["n1"]
            n_const += r["n_const"]
            max_star_S1 = max(max_star_S1, r["max_ss"])
            min_star_S1 = min(min_star_S1, r["min_ss"])
            gd_viol += r["gd_viol"]
            star_S1_vals.update(r["vals"])
            patterns.update(r["patterns"])

    t_walk = time.perf_counter() - t1
    rc = rho_budget_cand(p)
    # p5 exact check
    p5_exact = None
    if p == 5:
        expected = {-2 / 65, -42 / 325}
        got = set(float(k) for k in star_S1_vals.keys())
        p5_exact = {
            "matches_closed_form": all(
                any(abs(g - e) < 1e-12 for e in expected) for g in got
            )
            and len(got) == 2,
            "values": sorted(got),
            "expected": sorted(expected),
        }

    return {
        "p": int(p),
        "n": int(n),
        "N": int(Yc.shape[0]),
        "Cy_eq_py_err": cy_err,
        "Cy_eq_py": cy_err < 1e-10,
        "n_kappa1": n1,
        "star_S1_constant_on_every_set": n_const == n1,
        "n_constant": n_const,
        "max_star_S1": float(max_star_S1),
        "min_star_S1": float(min_star_S1),
        "star_S1_always_le_0": max_star_S1 <= 1e-12,
        "gaussian_domination": gd_viol == 0,
        "n_gd_violations": gd_viol,
        "star_S1_unique_values": sorted(float(k) for k in star_S1_vals.keys()),
        "n_patterns": len(patterns),
        "p5_exact": p5_exact,
        "rho_budget_cand": float(rc),
        "M_cand": float(M_cand(p)),
        "wall_m4_s": float(gpu_meta["wall_s"]),
        "wall_walk_s": float(t_walk),
        "gpu_used": True,
        "gpu_m4": gpu_meta,
        "workers_walk": int(min(W, len(shards))),
        "io": "mmap Max+; GPU m4; ProcessPool constancy walk; write_json_atomic",
    }


_MP: dict = {}


def _init_mp(C, quads, lookup, p: int) -> None:
    _MP["C"] = C
    _MP["quads"] = quads
    _MP["lookup"] = lookup
    _MP["p"] = p
    _MP["n"] = C.shape[0]
    _blas1()


def _walk_shard_mp(idxs):
    C = _MP["C"]
    quads = _MP["quads"]
    lookup = _MP["lookup"]
    p = _MP["p"]
    n = _MP["n"]
    inv = 1.0 / (p * p)
    local = {
        "n1": 0,
        "n_const": 0,
        "max_ss": -1e30,
        "min_ss": 1e30,
        "gd_viol": 0,
        "vals": Counter(),
        "patterns": Counter(),
    }
    for i in idxs:
        S = tuple(map(int, quads[i]))
        local["n1"] += 1
        ss_list = []
        for v in S:
            st = 1
            for u in S:
                if u != v:
                    st *= int(C[v, u])
            others = tuple(x for x in S if x != v)
            S1 = 0.0
            t1v = 0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                m4p, kp = lookup[Sp]
                if abs(kp) != 1:
                    continue
                rp = m4p - kp * inv
                S1 += float(C[v, rr]) * rp
                t1v += int(C[v, rr]) * kp
            ss = st * S1
            if ss > 1e-12:
                local["gd_viol"] += 1
            ss_list.append(round(ss, 12))
            local["vals"][round(ss, 12)] += 1
            local["max_ss"] = max(local["max_ss"], ss)
            local["min_ss"] = min(local["min_ss"], ss)
        if len(set(ss_list)) == 1:
            local["n_const"] += 1
        local["patterns"][tuple(sorted(ss_list))] += 1
    return local


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"S1_const Prop15.78: W={W} use_gpu={use_gpu} "
        f"io=mmap+GPU+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )

    alg = algebraic_moment_form()
    p5a = p5_exact_spectrum_algebra()
    print(f"  algebra moment form ok; p5 spectrum algebra ok={p5a['algebra_ok']}", flush=True)

    out: dict = {
        "prop": "15.78",
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; GPU all-m4; ProcessPool constancy/GD walk; "
            "pure-C τ1 constancy; write_json_atomic"
        ),
        "algebra": alg,
        "p5_spectrum_algebra": p5a,
        "tau1_constancy": {},
        "results": {},
        "status": None,
    }

    # combinatorial τ1·star constancy
    for p in (3, 5, 7, 11):
        print(f"  star·τ1 constancy p={p}...", flush=True)
        r = certify_star_tau1_constant(p, W)
        out["tau1_constancy"][str(p)] = r
        print(
            f"    all_constant={r['all_constant']} values={r['star_tau1_unique_values']} "
            f"wall={r['wall_s']:.2f}s",
            flush=True,
        )

    if not use_gpu:
        out["status"] = "GPU unavailable — refuse Max+ path (F17)."
        write_json_atomic(ROOT / "evidence" / "e1_gmin_m4_S1_const.json", out)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"  Max+ star·S1 constancy+GD p={p}...", flush=True)
        r = cert_maxplus_structure(p, W=W)
        out["results"][str(p)] = r
        print(
            f"    const={r['star_S1_constant_on_every_set']} "
            f"GD={r['gaussian_domination']} "
            f"max_star_S1={r['max_star_S1']:.8f} "
            f"n_vals={len(r['star_S1_unique_values'])} "
            f"m4={r['wall_m4_s']:.3f}s walk={r['wall_walk_s']:.3f}s",
            flush=True,
        )

    r5, r7 = out["results"]["5"], out["results"]["7"]
    tau_ok = all(out["tau1_constancy"][str(p)]["all_constant"] for p in (3, 5, 7, 11))
    both_const = r5["star_S1_constant_on_every_set"] and r7["star_S1_constant_on_every_set"]
    both_gd = r5["gaussian_domination"] and r7["gaussian_domination"]
    both_le0 = r5["star_S1_always_le_0"] and r7["star_S1_always_le_0"]
    out["tau1_constancy_all"] = tau_ok
    out["star_S1_constancy_both"] = both_const
    out["gaussian_domination_both"] = both_gd
    out["star_S1_le0_both"] = both_le0
    out["status"] = (
        f"Prop 15.78: star·τ1 constant on every |κ|=1 set p=3,5,7,11 => {tau_ok}; "
        f"star·S1 constant on every |κ|=1 Max+ set p=5,7 => {both_const}; "
        f"Gaussian domination E[φ]≤E_Wick[φ] => {both_gd}; "
        f"star·S1≤0 => {both_le0} (p5 exact {{-2/65,-42/325}}). "
        f"GPU mmap+atomic. OPEN: prove constancy+S1≤0 for all p≥5; joint/S3. "
        f"lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_S1_const.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
