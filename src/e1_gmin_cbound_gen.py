#!/usr/bin/env python3
"""
Generalize moduli c-pin (Prop 15.53–15.54 / 15.70) to p=5,7.

For each p:
  - parallel build of refined classes (ProcessPool over quad shards)
  - evec system (parallel over classes) → nullity-1 line m = m_part + c n
  - g_min(c) = -max |m| on |κ|=1 classes
  - Tr(G²) from Max+ dots via CORRECT all-pairs identity
      Tr = (1/4)(E[dot^4] - 2n E[dot²] + n²)
    with E over ALL ordered pairs incl. diagonal (matches trG2_from_maxplus_dots;
    off-diagonal-only was a known bug → empty roots at p=5)
  - pin c by Tr; check g_min ≥ -M_mid and ≥ L

F17: multi-worker class/evec; GPU for dense Gram Tr; mmap Max+; atomic evidence.
Writes evidence/e1_gmin_cbound_gen.json (+ optional /tmp class caches).
"""
from __future__ import annotations

import os
import pickle
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


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


# ---------------------------------------------------------------------------
# ProcessPool workers (module-level; file script only — stdin ProcessPool fails)
# ---------------------------------------------------------------------------

_G_C: np.ndarray | None = None
_G_KEYS: list | None = None
_G_IDX: dict | None = None
_G_Q2K: dict | None = None
_G_P: int | None = None
_G_N: int | None = None


def _init_C(C: np.ndarray) -> None:
    global _G_C
    _G_C = C
    _blas1()


def _init_evec(C: np.ndarray, p: int, keys: list, q2k: dict) -> None:
    global _G_C, _G_KEYS, _G_IDX, _G_Q2K, _G_P, _G_N
    _G_C = C
    _G_KEYS = keys
    _G_IDX = {k: i for i, k in enumerate(keys)}
    _G_Q2K = q2k
    _G_P = p
    _G_N = C.shape[0]
    _blas1()


def _shard_class_keys(quads: list[tuple]) -> list[tuple]:
    """Return list of (key, quad) for a shard of 4-sets."""
    from e1_gmin_moduli import class_key

    assert _G_C is not None
    C = _G_C
    out = []
    for q in quads:
        out.append((class_key(C, q), q))
    return out


def _evec_row_pack(pack: tuple) -> tuple[np.ndarray, float]:
    """One class → one averaged evec equation row. State from _init_evec."""
    from e1_gmin_moduli import m2_val, reduce_moment

    C = _G_C
    assert C is not None and _G_IDX is not None and _G_Q2K is not None
    p, n = int(_G_P), int(_G_N)
    kA, sample = pack
    idx = _G_IDX
    nv = len(_G_KEYS)  # type: ignore[arg-type]
    coef = np.zeros(nv)
    rhs = 0.0
    count = 0
    for quad in sample:
        pts = list(quad)
        for pos in range(4):
            i = pts[pos]
            others = [pts[t] for t in range(4) if t != pos]
            coef[idx[kA]] += p
            for r in range(n):
                if r == i:
                    continue
                red = reduce_moment([r] + others)
                Cr = float(C[i, r])
                if red[0] == "1":
                    rhs += Cr
                elif red[0] == "m2":
                    rhs += Cr * m2_val(C, p, red[1], red[2])
                elif red[0] == "m4":
                    coef[idx[_G_Q2K[red[1]]]] -= Cr
            count += 1
    return coef / count, rhs / count


def build_classes_parallel(C: np.ndarray, W: int) -> tuple[dict, dict]:
    """Shard C(n,4) over ProcessPool workers; merge class maps."""
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nq = len(quads)
    # ~W shards; avoid empty shards
    n_shards = min(W, max(1, nq // 200))
    shard_size = (nq + n_shards - 1) // n_shards
    shards = [quads[i : i + shard_size] for i in range(0, nq, shard_size)]

    classes: dict = defaultdict(list)
    quad_to_key: dict = {}
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_shard_class_keys, sh) for sh in shards]
        for fut in as_completed(futs):
            for k, q in fut.result():
                classes[k].append(q)
                quad_to_key[q] = k
    return dict(classes), quad_to_key


def build_evec_parallel(
    C: np.ndarray,
    p: int,
    classes: dict,
    q2k: dict,
    keys: list,
    *,
    max_quads: int,
    W: int,
) -> tuple[np.ndarray, np.ndarray]:
    n = C.shape[0]
    packs = []
    for kA, quads in classes.items():
        if len(quads) <= max_quads:
            sample = quads
        else:
            step = max(1, len(quads) // max_quads)
            sample = quads[::step][:max_quads]
        packs.append((kA, sample))

    rows = [None] * len(packs)
    bvec = [0.0] * len(packs)
    with ProcessPoolExecutor(
        max_workers=min(W, max(2, len(packs))),
        initializer=_init_evec,
        initargs=(C, p, keys, q2k),
    ) as ex:
        futs = {ex.submit(_evec_row_pack, pack): i for i, pack in enumerate(packs)}
        for fut in as_completed(futs):
            i = futs[fut]
            coef, rhs = fut.result()
            rows[i] = coef
            bvec[i] = rhs
    return np.array(rows), np.array(bvec)


def trG2_gpu_allpairs(Mp: np.ndarray) -> dict:
    """
    Tr(G²) matching e1_gmin_moduli.trG2_from_maxplus_dots exactly.

    Uses ALL ordered pairs (incl. diagonal). Off-diagonal-only is WRONG
    (p=5: ~1407 vs true 1808 → empty quadratic roots).
    """
    N, n = Mp.shape
    try:
        import cupy as cp

        Yg = cp.asarray(np.ascontiguousarray(Mp, dtype=np.float64))
        G = Yg @ Yg.T
        E_dot2 = float(cp.mean(G * G).get())
        E_dot4 = float(cp.mean(G ** 4).get())
        # K_ab = ((dot)^2 - n)/2; Tr = mean(K^2)
        K = (G * G - n) / 2.0
        Tr = float(cp.mean(K * K).get())
        via = "gpu_allpairs"
    except Exception as e:
        # CPU chunked (same identity as trG2_from_maxplus_dots)
        from e1_gmin_moduli import trG2_from_maxplus_dots

        d = trG2_from_maxplus_dots(Mp)
        return {
            "Tr_G2": d["Tr_G2"],
            "E_dot2": d["E_dot2"],
            "E_dot4": d["E_dot4"],
            "TrG2_from_moments": d["TrG2_from_moments"],
            "tr_via": f"cpu_chunked({e.__class__.__name__})",
        }
    Tr_mom = 0.25 * (E_dot4 - 2 * n * E_dot2 + n * n)
    return {
        "Tr_G2": Tr,
        "E_dot2": E_dot2,
        "E_dot4": E_dot4,
        "TrG2_from_moments": Tr_mom,
        "identity_match": abs(Tr - Tr_mom) < 1e-6,
        "tr_via": via,
    }


def cache_path(p: int, name: str) -> Path:
    return Path(f"/tmp/e1_cbound_p{p}_{name}")


def save_pickle_atomic(path: Path, obj) -> None:
    from io_atomic import write_bytes_atomic

    write_bytes_atomic(path, pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))


def load_pickle(path: Path):
    with open(path, "rb") as f:
        return pickle.load(f)


def run_p(p: int, W: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_moduli import (
        L_p,
        T_p,
        empirical_m4,
        kappa,
        m4_constancy,
        nullity_analysis,
        wick_dot4_upper,
    )
    from io_atomic import load_maxplus_array, write_npy_atomic
    from minmax_quadratic import paley_conference_prime_power

    t_all = time.perf_counter()
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    # mmap Max+ (Wieferich-style shared page cache)
    Mp = np.ascontiguousarray(load_maxplus_array(p, mmap=True), dtype=float)

    # --- classes (parallel) or cache ---
    cls_cache = cache_path(p, "classes.pkl")
    t0 = time.perf_counter()
    if cls_cache.is_file():
        classes, q2k = load_pickle(cls_cache)
        class_via = f"cache:{cls_cache}"
    else:
        classes, q2k = build_classes_parallel(C, W)
        save_pickle_atomic(cls_cache, (classes, q2k))
        class_via = f"parallel_W={W}"
    wall_class = time.perf_counter() - t0
    keys = list(classes.keys())
    idx = {k: i for i, k in enumerate(keys)}

    const = m4_constancy(Mp, classes, sample=30 if p >= 7 else 40)

    # --- evec system (parallel over classes) ---
    max_q = 40 if p >= 7 else 80
    t0 = time.perf_counter()
    evec_cache = cache_path(p, "evec.npz")
    if evec_cache.is_file():
        z = np.load(evec_cache)
        A, b = z["A"], z["b"]
        # keys order must match — recompute if n_var mismatch
        if A.shape[1] != len(keys):
            A, b = build_evec_parallel(C, p, classes, q2k, keys, max_quads=max_q, W=W)
            write_npy_atomic(str(evec_cache).replace(".npz", "_A.npy"), A)  # keep both
            np.savez(str(evec_cache) + ".tmp", A=A, b=b)
            os.replace(str(evec_cache) + ".tmp", evec_cache)
            evec_via = f"rebuild_parallel_W={W}"
        else:
            evec_via = f"cache:{evec_cache}"
    else:
        A, b = build_evec_parallel(C, p, classes, q2k, keys, max_quads=max_q, W=W)
        # atomic-ish: write tmp then replace
        tmp = str(evec_cache) + ".tmp.npz"
        np.savez(tmp, A=A, b=b)
        os.replace(tmp, evec_cache)
        evec_via = f"parallel_W={W}"
    wall_evec = time.perf_counter() - t0

    m_true = empirical_m4(Mp, classes, keys)
    null = nullity_analysis(A, b, m_true)

    # --- Tr(G²) correct all-pairs (GPU preferred) ---
    t0 = time.perf_counter()
    tr_info = trG2_gpu_allpairs(Mp)
    wall_tr = time.perf_counter() - t0
    Tr = float(tr_info["Tr_G2"])

    wick = wick_dot4_upper(p)
    TrW = float(wick["Tr_G2_wick"])
    mid = M_mid(p)

    base = {
        "p": int(p),
        "n": int(n),
        "N": int(Mp.shape[0]),
        "n_classes": len(keys),
        "io": {
            "maxplus": "mmap load_maxplus_array",
            "class_cache": str(cls_cache),
            "class_via": class_via,
            "evec_via": evec_via,
            "evidence": "write_json_atomic",
        },
        "wall_class_s": float(wall_class),
        "wall_evec_s": float(wall_evec),
        "wall_tr_s": float(wall_tr),
        "wall_total_s": float(time.perf_counter() - t_all),
        "constancy_all": bool(const.get("all_constant")),
        "constancy_max_std": float(const.get("max_std", float("nan"))),
        "nullity": int(null["nullity"]),
        "Tr_G2": float(Tr),
        "Tr_G2_from_moments": float(tr_info.get("TrG2_from_moments", Tr)),
        "Tr_identity_match": bool(tr_info.get("identity_match", True)),
        "E_dot2": float(tr_info["E_dot2"]),
        "E_dot4": float(tr_info["E_dot4"]),
        "Tr_G2_wick": float(TrW),
        "tr_via": tr_info["tr_via"],
        "M_mid": mid,
        "L": float(L_p(p)),
        "T": float(T_p(p)),
        # empirical max |m4| on class reps with |kappa|=1 (even if non-constant)
        "emp_max_abs_m4_kappa1": None,
        "selected_g_min_true_Tr": None,
        "selected_ge_neg_mid": False,
        "selected_ge_L": False,
        "selected_gt_T": False,
        "recovered_p5_target": False,
        "pin_status": "ok",
    }

    # Empirical bound from class representatives (always available)
    kap1_abs = []
    for k in keys:
        a, b_, c, d = classes[k][0]
        kap = kappa(C, (a, b_, c, d))
        if abs(kap) == 1:
            kap1_abs.append(abs(float(m_true[idx[k]])))
    if kap1_abs:
        base["emp_max_abs_m4_kappa1"] = float(max(kap1_abs))
        base["emp_g_min"] = -float(max(kap1_abs))
        base["emp_ge_neg_mid"] = base["emp_g_min"] >= -mid - 1e-12
        base["emp_ge_L"] = base["emp_g_min"] >= float(L_p(p)) - 1e-12

    if "m_part" not in null or "nvec" not in null:
        base["pin_status"] = (
            f"skip_pin: nullity={null['nullity']} constancy={const.get('all_constant')} "
            f"(need nullity≥1 and m_part/nvec; p=7 coarse classes not m4-constant)"
        )
        base["wall_total_s"] = float(time.perf_counter() - t_all)
        return base

    m_part = np.array(null["m_part"])
    nvec = np.array(null["nvec"])
    c_true = null.get("c_true_fit")

    kap1 = []
    for k in keys:
        a, b_, c, d = classes[k][0]
        kap = kappa(C, (a, b_, c, d))
        if abs(kap) == 1:
            kap1.append(idx[k])

    def gmin_c(c: float) -> float:
        m = m_part + c * nvec
        return -max(abs(m[i]) for i in kap1) if kap1 else float("nan")

    # Tr quadratic coefficients via combinatorial G_part, G_null
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
                cpv = float(C[i, j] * C[k, l])
                G_part[ei, ej] = G_part[ej, ei] = cpv * m_part[ci]
                G_null[ei, ej] = G_null[ej, ei] = cpv * nvec[ci]

    a2 = float(np.sum(G_null * G_null))
    a1 = 2.0 * float(np.sum(G_part * G_null))
    a0 = float(np.sum(G_part * G_part))

    def roots_at(target: float) -> list[float]:
        Aq, Bq, Cq = a2, a1, a0 - target
        disc = Bq * Bq - 4 * Aq * Cq
        if disc < 0 or abs(Aq) < 1e-30:
            return []
        sd = float(np.sqrt(disc))
        return [(-Bq - sd) / (2 * Aq), (-Bq + sd) / (2 * Aq)]

    roots_true = roots_at(Tr)
    roots_wick = roots_at(TrW)
    g_true = [gmin_c(c) for c in roots_true]
    g_wick = [gmin_c(c) for c in roots_wick]
    sel_true = max(g_true) if g_true else None
    sel_wick = max(g_wick) if g_wick else None

    cs = np.linspace(-2, 2, 801)
    gs = np.array([gmin_c(float(c)) for c in cs])
    hits = []
    for i in range(len(cs) - 1):
        if (gs[i] + mid) * (gs[i + 1] + mid) <= 0:
            hits.append(float(0.5 * (cs[i] + cs[i + 1])))

    base.update(
        {
            "c_true_fit": float(c_true) if c_true is not None else None,
            "g_min_at_c_true": float(gmin_c(c_true)) if c_true is not None else None,
            "tr_quadratic": {"a0": a0, "a1": a1, "a2": a2},
            "roots_true_Tr": roots_true,
            "g_min_at_true_roots": g_true,
            "selected_g_min_true_Tr": sel_true,
            "selected_g_min_wick_Tr": sel_wick,
            "selected_ge_neg_mid": sel_true is not None
            and sel_true >= -mid - 1e-12,
            "selected_ge_L": sel_true is not None
            and sel_true >= float(L_p(p)) - 1e-12,
            "selected_gt_T": sel_true is not None and sel_true > float(T_p(p)),
            "c_where_gmin_crosses_neg_mid": hits,
            "max_gmin_on_scan": float(np.max(gs)),
            "min_gmin_on_scan": float(np.min(gs)),
            "recovered_p5_target": (
                p == 5
                and sel_true is not None
                and abs(sel_true + 3.0 / 65.0) < 1e-9
            ),
            "wall_total_s": float(time.perf_counter() - t_all),
        }
    )
    return base


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot
    from io_atomic import write_json_atomic
    from workers import require_workers, assert_not_single_core_thrash

    W = require_workers()
    snap = compute_snapshot()
    print(
        f"cbound_gen: W={W} gpu={snap['gpu'].get('available')} "
        f"io=mmap+atomic (Wieferich-style)",
        flush=True,
    )

    out: dict = {
        "workers": W,
        "compute": snap,
        "io_policy": (
            "mmap Max+ shared across workers; pickle class caches via write_bytes_atomic; "
            "evidence via write_json_atomic (tmp+fsync+replace)"
        ),
        "results": {},
    }
    # Sequential over p (each uses full W for classes + one GPU context for Tr)
    for p in (5, 7):
        print(f"  run_p({p}) W={W}...", flush=True)
        out["results"][str(p)] = run_p(p, W)
        r = out["results"][str(p)]
        print(
            f"    nullity={r['nullity']} sel_g={r['selected_g_min_true_Tr']} "
            f"ge_mid={r['selected_ge_neg_mid']} ge_L={r['selected_ge_L']} "
            f"Tr={r['Tr_G2']:.4f} via={r['tr_via']} "
            f"class={r['wall_class_s']:.1f}s evec={r['wall_evec_s']:.1f}s "
            f"tr={r['wall_tr_s']:.2f}s total={r['wall_total_s']:.1f}s",
            flush=True,
        )
        # F17 live check between primes
        try:
            assert_not_single_core_thrash()
        except Exception as e:
            print(f"  F17 warn: {e}", flush=True)

    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["status"] = (
        f"Moduli pin p=5: g_min={r5['selected_g_min_true_Tr']} "
        f"ge_mid={r5['selected_ge_neg_mid']} recovered_3/65={r5.get('recovered_p5_target')}; "
        f"p=7: pin_status={r7.get('pin_status')} emp_g={r7.get('emp_g_min')} "
        f"emp_ge_mid={r7.get('emp_ge_neg_mid')} nullity={r7.get('nullity')} "
        f"constancy={r7.get('constancy_all')}. "
        f"OPEN: Max+-free Tr(G²)/E[dot^4] + constant-m4 classes for general p≥5."
    )
    path = ROOT / "evidence" / "e1_gmin_cbound_gen.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
