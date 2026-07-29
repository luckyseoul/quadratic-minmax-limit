#!/usr/bin/env python3
"""
Prop 15.72: resolvent-gain structure from n1,n3,d1,d3 (toward |m4|≤L).

Proved (Max+-free conference / Paley combinatorics):
  1. On |κ|=3, Tκ/κ ∈ {±8} (from Prop 15.68: Tκ∈{±24}, |κ|=3).
  2. Given n1,n3 (Prop 15.71) and constant d3=p²-5 on |κ|=1 (Prop 15.68),
     handshaking forces reverse degrees on |κ|=3:
       d1^(3) = n1·d3 / n3 = 3(p²-1),   d3^(3) = p²-9.
  3. Target algebra: same-sign |ρ|≤(p-4)/(2p²) ⇔ |m4|≤L_abs on |κ|=1;
     gain ≤(p-4)/48 with source amp 24/p² is sufficient.
  4. Separate vanishing (certified p=3,5,7): on every |κ|=1 4-set,
       ∑_{ext→|κ|=1} C_vr κ(S') = 0 and ∑_{ext→|κ|=3} C_vr κ(S') = 0
     (strictly stronger than Tκ=0, which is only the sum of the two).

Certified multi-worker:
  - reverse degrees constant on |κ|=3 at p=3,5,7
  - separate vanishing full census p=3,5,7
  - type6 Max+-free (4pI-T)ρ=Tκ/p² solve: max|m| and same-sign |ρ| vs budgets
  - empirical Max+ gain p=5,7 (mmap) when caches exist

OPEN: prove gain ≤(p-4)/48 (or |m4|≤M_mid/L) for all primes p≥5.
Writes evidence/e1_gmin_m4_resolvent_gain.json (atomic).
"""
from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations, permutations
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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def rho_budget(p: int) -> float:
    return (p - 4) / (2.0 * p * p)


def gain_budget(p: int) -> float:
    return (p - 4) / 48.0


def source_amp(p: int) -> float:
    return 24.0 / (p * p)


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    )


def star_sum(C, S) -> int:
    s = 0
    for v in S:
        pr = 1
        for u in S:
            if u != v:
                pr *= int(C[v, u])
        s += pr
    return s


def algebra_reverse_degrees() -> dict:
    """Proved: d1^(3)=3(p²-1), d3^(3)=p²-9 from n1,n3,d3 via handshaking."""
    from e1_gmin_m4_stratum import (
        d3_formula,
        n1_formula,
        n3_formula,
        n_of_p,
    )

    rows = {}
    ok_all = True
    for p in range(3, 60, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        n = n_of_p(p)
        n1, n3 = n1_formula(p), n3_formula(p)
        d3 = d3_formula(p)  # p²-5 = n-6
        # handshaking
        d1_3 = n1 * d3 // n3
        d3_3 = 4 * (n - 4) - d1_3
        pred_d1_3 = 3 * (p * p - 1)
        pred_d3_3 = p * p - 9
        match = d1_3 == pred_d1_3 and d3_3 == pred_d3_3
        # also check exact rational identity without integer div
        exact = abs(n1 * d3 / n3 - pred_d1_3) < 1e-12
        ok = match and exact and d1_3 + d3_3 == 4 * (n - 4)
        ok_all = ok_all and ok
        rows[str(p)] = {
            "n1": n1,
            "n3": n3,
            "d3_from_kappa1": d3,
            "d1_from_kappa3": d1_3,
            "d3_from_kappa3": d3_3,
            "pred_d1_3": pred_d1_3,
            "pred_d3_3": pred_d3_3,
            "match": match,
            "ok": ok,
        }
    return {
        "identity_d1_from_kappa3": "n1*d3/n3 = 3(p^2-1)",
        "identity_d3_from_kappa3": "4(n-4)-3(p^2-1) = p^2-9",
        "hypothesis": "d3=p^2-5 constant on every |κ|=1 4-set (Prop 15.68 / 15.71)",
        "by_p": rows,
        "proved_under_hypothesis": ok_all,
        "note": (
            "Handshaking e13=n1*d3=n3*d1^(3) is unconditional once d3 is constant "
            "on |κ|=1. Constancy certified for Paley p=3,5,7,11."
        ),
    }


def algebra_gain_targets() -> dict:
    rows = {}
    for p in range(5, 60, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        rb, gb, sa = rho_budget(p), gain_budget(p), source_amp(p)
        # sufficient: same-sign |ρ| ≤ gb * sa = (p-4)/48 * 24/p² = (p-4)/(2p²) = rb
        rows[str(p)] = {
            "rho_budget": rb,
            "gain_budget": gb,
            "source_amp_24_over_p2": sa,
            "gain_times_source_eq_rho_budget": abs(gb * sa - rb) < 1e-15,
            "L_abs": L_abs(p),
            "M_mid": M_mid(p),
            "M_cand": M_cand(p),
            "wick": 1.0 / (p * p),
            "wick_plus_rho_budget": 1.0 / (p * p) + rb,
            "equals_L": abs(1.0 / (p * p) + rb - L_abs(p)) < 1e-15,
        }
    return {
        "sufficient_condition": (
            "same-sign |ρ|≤(p-4)/(2p²) on |κ|=1 ⇒ |m4|≤L_abs; "
            "equivalently resolvent gain from |κ|=3 source amp 24/p² into same-sign "
            "|κ|=1 is ≤(p-4)/48"
        ),
        "by_p": rows,
        "proved_algebra": all(r["equals_L"] and r["gain_times_source_eq_rho_budget"] for r in rows.values()),
    }


def tkappa_over_kappa_identity() -> dict:
    """K4: on |κ|=3, Tκ/κ ∈ {±8}."""
    from itertools import product

    edges = list(combinations(range(4), 2))
    ratios = set()
    for bits in product([-1, 1], repeat=6):
        e = dict(zip(edges, bits))

        def E(i, j):
            return e[tuple(sorted((i, j)))]

        kap = E(0, 1) * E(2, 3) + E(0, 2) * E(1, 3) + E(0, 3) * E(1, 2)
        star = 0
        for v in range(4):
            pr = 1
            for u in range(4):
                if u != v:
                    pr *= E(v, u)
            star += pr
        T = -6 * star
        if abs(kap) == 3:
            ratios.add(T // kap)
    return {
        "proved": ratios <= {8, -8} and ratios == {8, -8},
        "Tkappa_over_kappa_on_abs3": sorted(ratios),
        "method": "64 K4 edge labelings; Tκ=-6*star (Prop 15.68)",
    }


# ---------------------------------------------------------------------------
# ProcessPool workers
# ---------------------------------------------------------------------------
_G: dict = {}


def _init_C(C: np.ndarray) -> None:
    _G["C"] = C
    _blas1()


def _vanish_shard(quads: list) -> dict:
    """Per |κ|=1 set: sk1, sk3, d1, d3."""
    C = _G["C"]
    n = C.shape[0]
    n1 = 0
    sk1_nonzero = sk3_nonzero = 0
    d1s, d3s = set(), set()
    for S in quads:
        k = kappa(C, S)
        if abs(k) != 1:
            continue
        n1 += 1
        Sset = set(S)
        sk1 = sk3 = 0.0
        d1 = d3 = 0
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                kp = kappa(C, Sp)
                w = float(C[v, r])
                if abs(kp) == 1:
                    d1 += 1
                    sk1 += w * kp
                elif abs(kp) == 3:
                    d3 += 1
                    sk3 += w * kp
        d1s.add(d1)
        d3s.add(d3)
        if abs(sk1) > 1e-9:
            sk1_nonzero += 1
        if abs(sk3) > 1e-9:
            sk3_nonzero += 1
    return {
        "n1": n1,
        "sk1_nonzero": sk1_nonzero,
        "sk3_nonzero": sk3_nonzero,
        "d1s": list(d1s),
        "d3s": list(d3s),
    }


def _k3_deg_shard(quads: list) -> dict:
    C = _G["C"]
    n = C.shape[0]
    n3 = 0
    d1s, d3s = set(), set()
    tkr = set()
    for S in quads:
        k = kappa(C, S)
        if abs(k) != 3:
            continue
        n3 += 1
        Sset = set(S)
        d1 = d3 = 0
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                ak = abs(kappa(C, Sp))
                if ak == 1:
                    d1 += 1
                elif ak == 3:
                    d3 += 1
        d1s.add(d1)
        d3s.add(d3)
        T = -6 * star_sum(C, S)
        tkr.add(T // k)
    return {
        "n3": n3,
        "d1s": list(d1s),
        "d3s": list(d3s),
        "Tkappa_over_kappa": list(tkr),
    }


def certify_vanishing_and_degrees(p: int, W: int) -> dict:
    _blas1()
    from e1_gmin_m4_stratum import d1_formula, d3_formula, n1_formula, n3_formula
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    n_shards = min(W, max(1, len(quads) // 400))
    ss = (len(quads) + n_shards - 1) // n_shards
    shards = [quads[i : i + ss] for i in range(0, len(quads), ss)]

    n1 = sk1_nz = sk3_nz = 0
    d1s: set[int] = set()
    d3s: set[int] = set()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_vanish_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            r = fut.result()
            n1 += r["n1"]
            sk1_nz += r["sk1_nonzero"]
            sk3_nz += r["sk3_nonzero"]
            d1s.update(r["d1s"])
            d3s.update(r["d3s"])

    n3 = 0
    d1_3: set[int] = set()
    d3_3: set[int] = set()
    tkr: set[int] = set()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_k3_deg_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            r = fut.result()
            n3 += r["n3"]
            d1_3.update(r["d1s"])
            d3_3.update(r["d3s"])
            tkr.update(r["Tkappa_over_kappa"])

    pred_d1_3 = 3 * (p * p - 1)
    pred_d3_3 = p * p - 9
    return {
        "p": int(p),
        "n": int(n),
        "wall_s": float(time.perf_counter() - t0),
        "n1_census": n1,
        "n1_formula": n1_formula(p),
        "n3_census": n3,
        "n3_formula": n3_formula(p),
        "separate_vanishing_sk1": sk1_nz == 0,
        "separate_vanishing_sk3": sk3_nz == 0,
        "sk1_nonzero_count": sk1_nz,
        "sk3_nonzero_count": sk3_nz,
        "d1_kappa1_unique": sorted(d1s),
        "d3_kappa1_unique": sorted(d3s),
        "d1_kappa1_match": d1s == {d1_formula(p)},
        "d3_kappa1_match": d3s == {d3_formula(p)},
        "d1_kappa3_unique": sorted(d1_3),
        "d3_kappa3_unique": sorted(d3_3),
        "d1_kappa3_match": d1_3 == {pred_d1_3},
        "d3_kappa3_match": d3_3 == {pred_d3_3},
        "Tkappa_over_kappa_unique": sorted(tkr),
        "Tkappa_over_kappa_match_pm8": tkr <= {8, -8} and tkr == {8, -8},
    }


def type6_canon(C, S):
    a, b, c, d = S
    raw = {
        (0, 1): int(C[a, b]),
        (0, 2): int(C[a, c]),
        (0, 3): int(C[a, d]),
        (1, 2): int(C[b, c]),
        (1, 3): int(C[b, d]),
        (2, 3): int(C[c, d]),
    }
    best = None
    for perm in permutations(range(4)):
        tup = []
        for i, j in combinations(range(4), 2):
            oi, oj = perm[i], perm[j]
            tup.append(raw[tuple(sorted((oi, oj)))])
        t = tuple(tup)
        if best is None or t < best:
            best = t
    return best


def _cls_shard(quads: list) -> dict:
    """Build type6 class map for a shard of 4-sets."""
    C = _G["C"]
    _blas1()
    out: dict = defaultdict(list)
    for S in quads:
        out[type6_canon(C, S)].append(S)
    return dict(out)


def _type6_row(pack: tuple) -> tuple[int, np.ndarray]:
    """Average T row for one type6 class (sample of quads)."""
    C = _G["C"]
    n = C.shape[0]
    i, take, kidx_items = pack
    kidx = dict(kidx_items)
    nk = len(kidx)
    acc = np.zeros(nk)
    for S in take:
        Sset = set(S)
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                acc[kidx[type6_canon(C, Sp)]] += float(C[v, r])
    return i, acc / max(len(take), 1)


def type6_resolvent_solve(p: int, W: int) -> dict:
    """Max+-free type6-averaged (4pI-T)ρ = Tκ/p²."""
    _blas1()
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    classes: dict = defaultdict(list)
    # class build: multi-worker shards
    quads = list(combinations(range(n), 4))
    n_shards = min(W, max(1, len(quads) // 300))
    ss = (len(quads) + n_shards - 1) // n_shards
    shards = [quads[i : i + ss] for i in range(0, len(quads), ss)]

    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_cls_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            for k, lst in fut.result().items():
                classes[k].extend(lst)

    keys = sorted(classes.keys())
    nk = len(keys)
    kidx = {k: i for i, k in enumerate(keys)}
    kap = np.array([float(kappa(C, classes[k][0])) for k in keys])
    Tkap = np.array([-6.0 * star_sum(C, classes[k][0]) for k in keys])

    # T matrix rows in parallel
    packs = []
    rng = np.random.default_rng(p)
    for i, kA in enumerate(keys):
        quads_c = classes[kA]
        if len(quads_c) <= 150 or p <= 5:
            take = quads_c
        else:
            take = [
                quads_c[j]
                for j in rng.choice(len(quads_c), 150, replace=False)
            ]
        packs.append((i, take, list(kidx.items())))

    Tmat = np.zeros((nk, nk))
    with ProcessPoolExecutor(
        max_workers=min(W, max(2, len(packs))),
        initializer=_init_C,
        initargs=(C,),
    ) as ex:
        futs = [ex.submit(_type6_row, pack) for pack in packs]
        for fut in as_completed(futs):
            i, row = fut.result()
            Tmat[i] = row

    A = 4 * p * np.eye(nk) - Tmat
    b = Tkap / (p * p)
    rho, *_ = np.linalg.lstsq(A, b, rcond=None)
    m4 = kap / (p * p) + rho
    resid = float(np.max(np.abs(A @ rho - b)))
    evals = np.linalg.eigvals(Tmat).real
    idx1 = np.where(np.abs(kap) == 1)[0]
    max_m1 = float(np.max(np.abs(m4[idx1]))) if len(idx1) else None
    max_rho = float(np.max(np.abs(rho[idx1]))) if len(idx1) else None
    same = [
        abs(float(rho[i]))
        for i in idx1
        if float(rho[i]) * float(kap[i]) > 0
    ]
    max_same = float(max(same)) if same else 0.0
    gb = gain_budget(p) if p > 4 else None
    rb = rho_budget(p) if p > 4 else None
    sa = source_amp(p)
    eff_gain = max_same / sa if sa > 0 else None

    return {
        "p": int(p),
        "n_type6": nk,
        "wall_s": float(time.perf_counter() - t0),
        "lam_max_T_type6": float(np.max(evals)),
        "lam_min_T_type6": float(np.min(evals)),
        "gap_4p": float(4 * p - np.max(evals)),
        "lstsq_resid": resid,
        "max_abs_m4_type6_kappa1": max_m1,
        "max_abs_rho_kappa1": max_rho,
        "max_abs_rho_same_sign": max_same,
        "effective_gain": eff_gain,
        "gain_budget": gb,
        "rho_budget": rb,
        "source_amp": sa,
        "L_abs": L_abs(p),
        "M_mid": M_mid(p),
        "M_cand": M_cand(p),
        "type6_m4_le_L": max_m1 is not None and max_m1 <= L_abs(p) + 1e-9,
        "type6_m4_le_mid": max_m1 is not None and max_m1 <= M_mid(p) + 1e-9,
        "type6_m4_le_cand": max_m1 is not None and max_m1 <= M_cand(p) + 1e-9,
        "same_rho_le_budget": rb is None or max_same <= rb + 1e-9,
        "gain_le_budget": gb is None or (eff_gain is not None and eff_gain <= gb + 1e-9),
        "io": "Max+-free type6; ProcessPool T rows",
    }


def empirical_maxplus_gain(p: int) -> dict | None:
    """True Max+ same-sign |ρ| / (24/p²) via mmap."""
    try:
        from io_atomic import load_maxplus_array
        from minmax_quadratic import paley_conference_prime_power
    except Exception:
        return None
    try:
        Y = np.ascontiguousarray(load_maxplus_array(p, mmap=True), dtype=float)
    except FileNotFoundError:
        return None
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    max_m = max_rho = max_same = 0.0
    n1 = 0
    # full for p=5; stride for p=7
    step = 1 if p <= 5 else 15
    for S in list(combinations(range(n), 4))[::step]:
        k = kappa(C, S)
        if abs(k) != 1:
            continue
        n1 += 1
        a, b, c, d = S
        m4 = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        rho = m4 - k / (p * p)
        max_m = max(max_m, abs(m4))
        max_rho = max(max_rho, abs(rho))
        if rho * k > 0:
            max_same = max(max_same, abs(rho))
    sa = source_amp(p)
    return {
        "p": int(p),
        "n_kappa1_sampled": n1,
        "stride": step,
        "max_abs_m4": max_m,
        "max_abs_rho": max_rho,
        "max_abs_rho_same_sign": max_same,
        "effective_gain": max_same / sa,
        "gain_budget": gain_budget(p),
        "rho_budget": rho_budget(p),
        "L_abs": L_abs(p),
        "M_mid": M_mid(p),
        "m4_le_L": max_m <= L_abs(p) + 1e-9,
        "m4_le_mid": max_m <= M_mid(p) + 1e-9,
        "same_le_budget": max_same <= rho_budget(p) + 1e-9,
        "gain_le_budget": max_same / sa <= gain_budget(p) + 1e-9,
        "io": "mmap Max+",
    }


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    print(f"resolvent_gain Prop15.72: W={W}", flush=True)

    tk = tkappa_over_kappa_identity()
    print(f"  Tκ/κ on |κ|=3: {tk}", flush=True)
    rev = algebra_reverse_degrees()
    print(f"  reverse degrees algebra: {rev['proved_under_hypothesis']}", flush=True)
    tgt = algebra_gain_targets()
    print(f"  gain target algebra: {tgt['proved_algebra']}", flush=True)

    certs = {}
    for p in (3, 5, 7):
        print(f"  vanish+degrees p={p}...", flush=True)
        certs[str(p)] = certify_vanishing_and_degrees(p, W)
        r = certs[str(p)]
        print(
            f"    vanish sk1={r['separate_vanishing_sk1']} sk3={r['separate_vanishing_sk3']} "
            f"d3_k3={r['d3_kappa3_match']} Tκ/κ={r['Tkappa_over_kappa_match_pm8']} "
            f"wall={r['wall_s']:.2f}s",
            flush=True,
        )

    type6 = {}
    for p in (3, 5, 7):
        print(f"  type6 resolvent p={p}...", flush=True)
        type6[str(p)] = type6_resolvent_solve(p, W)
        r = type6[str(p)]
        print(
            f"    max|m4|={r['max_abs_m4_type6_kappa1']} le_L={r['type6_m4_le_L']} "
            f"same_ρ={r['max_abs_rho_same_sign']} gain={r['effective_gain']} "
            f"le_bud={r['gain_le_budget']} wall={r['wall_s']:.2f}s",
            flush=True,
        )

    emp = {}
    for p in (5, 7):
        print(f"  empirical Max+ gain p={p}...", flush=True)
        e = empirical_maxplus_gain(p)
        if e:
            emp[str(p)] = e
            print(
                f"    max|m4|={e['max_abs_m4']:.8f} same_ρ={e['max_abs_rho_same_sign']:.8f} "
                f"gain={e['effective_gain']:.8f} le_bud={e['gain_le_budget']}",
                flush=True,
            )

    certs_ok = all(
        certs[k]["separate_vanishing_sk1"]
        and certs[k]["separate_vanishing_sk3"]
        and certs[k]["d1_kappa3_match"]
        and certs[k]["d3_kappa3_match"]
        for k in certs
    )
    type6_L_ok = all(type6[k]["type6_m4_le_L"] for k in type6 if int(k) >= 5)

    out = {
        "workers": W,
        "prop": "15.72",
        "tkappa_over_kappa": tk,
        "reverse_degrees_algebra": rev,
        "gain_targets": tgt,
        "certs_vanish_degrees": certs,
        "type6_resolvent": type6,
        "empirical_maxplus_gain": emp,
        "certs_ok": certs_ok,
        "type6_le_L_for_p_ge_5": type6_L_ok,
        "io": "write_json_atomic; mmap Max+ for empirical; ProcessPool type6/vanish",
        "status": (
            f"Prop 15.72: reverse degrees d1^(3)=3(p^2-1), d3^(3)=p^2-9 proved from "
            f"n1,n3,d3 handshaking; Tκ/κ∈{{±8}} on |κ|=3 proved; gain⇔L algebra proved; "
            f"separate κ-weighted vanishing + reverse deg certs p=3,5,7 ({certs_ok}); "
            f"type6 Max+-free resolvent le_L p≥5 ({type6_L_ok}). "
            f"OPEN: prove resolvent gain ≤(p-4)/48 (or |m4|≤M_mid/L) for all primes p≥5 "
            f"without per-p Max+/type6 census. lim α_n remains OPEN."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_resolvent_gain.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
