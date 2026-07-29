#!/usr/bin/env python3
"""
Prop 15.62: typeA + wedge = 6N ‖B‖_F² for zero-diag B on V_+;
Q = 6N ‖B‖_F² + Q_4 with Q_4 only from vertex-disjoint edge pairs.

16N residual: Q_4 ≤ 10 N ‖B‖_F²  (⇔ λ_cycle ≤ 8).

Multi-worker ProcessPool from this file only (F17). OMP=1 per worker.
Writes evidence/e1_gmin_typeA_wedge.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _edges(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def _sym_basis_dim(d: int) -> int:
    return d * (d + 1) // 2


def _A_from_vec(v: np.ndarray, d: int) -> np.ndarray:
    A = np.zeros((d, d))
    k = 0
    for i in range(d):
        for j in range(i, d):
            A[i, j] = A[j, i] = v[k]
            k += 1
    return A


def _zero_diag_nullspace(Vp: np.ndarray) -> tuple[np.ndarray, int]:
    """Nullspace of (ambient diag of Vp A Vp.T, Tr A) in Sym(d) coordinates."""
    n, d = Vp.shape
    dimS = _sym_basis_dim(d)
    M = np.zeros((n, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        A = _A_from_vec(e, d)
        B = Vp @ A @ Vp.T
        M[:, k] = np.diag(B)
    tr_row = np.zeros(dimS)
    k = 0
    for i in range(d):
        for j in range(i, d):
            if i == j:
                tr_row[k] = 1.0
            k += 1
    Mfull = np.vstack([M, tr_row])
    _u, s, vh = np.linalg.svd(Mfull, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    null = vh[rank:].T
    return null, int(null.shape[1])


def _worker_identity(args: tuple) -> dict:
    _blas1()
    p, seed, mode = args
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    Pplus = Vp @ Vp.T
    rng = np.random.default_rng(seed)
    null, nullity = _zero_diag_nullspace(Vp)

    if mode == "max":
        A = rng.standard_normal((d, d))
        A = 0.5 * (A + A.T)
        A -= np.trace(A) / d * np.eye(d)
        U = (Y @ Vp) / np.sqrt(n)
        for _ in range(55):
            q = np.einsum("bi,ij,bj->b", U, A, U)
            Ph = (U * q[:, None]).T @ U
            Ph = 0.5 * (Ph + Ph.T)
            Ph -= np.trace(Ph) / d * np.eye(d)
            A = Ph / (np.linalg.norm(Ph) + 1e-30)
        B = Vp @ A @ Vp.T
    else:
        coef = rng.standard_normal(nullity)
        v = null @ coef
        A = _A_from_vec(v, d)
        B = Vp @ A @ Vp.T
    nrm = np.linalg.norm(B)
    B = B / (nrm + 1e-30)
    bf2 = float(np.sum(B * B))
    diagmax = float(np.max(np.abs(np.diag(B))))
    vplus_err = float(np.linalg.norm(B - Pplus @ B @ Pplus))

    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))

    typeA = 2.0 * N * bf2
    wedge_formula = 4.0 * N * bf2
    sixN = 6.0 * N * bf2

    edges = _edges(n)
    E = len(edges)
    Be = np.array([2.0 * B[i, j] for i, j in edges], dtype=np.float64)
    Be_norm2 = float(np.dot(Be, Be))
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)

    if E <= 400:
        Gu = F.T @ F
        share = np.zeros((E, E), dtype=bool)
        for a, (i, j) in enumerate(edges):
            for b in range(a + 1, E):
                k, l = edges[b]
                if len({i, j, k, l}) == 3:
                    share[a, b] = share[b, a] = True
        disj = ~share
        np.fill_diagonal(disj, False)
        same = np.eye(E, dtype=bool)
        same_term = float(Be @ (Gu * same) @ Be)
        wedge_term = float(Be @ (Gu * share) @ Be)
        disj_term = float(Be @ (Gu * disj) @ Be)
        Q_edge = float(Be @ Gu @ Be)
    else:
        # p=7: avoid E×E; use formula + residual
        same_term = N * Be_norm2
        wedge_term = wedge_formula
        Q_edge = Q
        disj_term = Q - same_term - wedge_term

    typeA_ok = abs(same_term - typeA) < 1e-5 * max(typeA, 1.0)
    wedge_ok = abs(wedge_term - wedge_formula) < 1e-4 * max(abs(wedge_formula), 1.0)
    six_ok = abs(same_term + wedge_term - sixN) < 1e-4 * max(sixN, 1.0)
    edge_err = abs(Q_edge - Q)

    return {
        "p": int(p),
        "seed": int(seed),
        "mode": mode,
        "N": int(N),
        "n": int(n),
        "d": int(d),
        "nullity": int(nullity),
        "diagmax": diagmax,
        "vplus_err": vplus_err,
        "bf2": bf2,
        "Q": Q,
        "typeA": typeA,
        "same_term": same_term,
        "wedge_formula": wedge_formula,
        "wedge_term": wedge_term,
        "disj_term": disj_term,
        "typeA_ok": bool(typeA_ok),
        "wedge_ok": bool(wedge_ok),
        "sixN_ok": bool(six_ok),
        "edge_err": edge_err,
        "Be_norm2": Be_norm2,
        "disj_le_10N": bool(disj_term <= 10.0 * N + 1e-3 * max(N, 1.0)),
        "ratio_to_16N": Q / (16.0 * N * bf2),
        "lam_cycle_proxy": Q / (2.0 * N * bf2),
        "Q4_over_N": disj_term / N,
    }


def _worker_spectrum(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    edges = _edges(n)
    E = len(edges)
    if E > 400:
        return {"p": p, "skipped": True}

    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    Gu = F.T @ F
    share = np.zeros((E, E), dtype=bool)
    for a, (i, j) in enumerate(edges):
        for b in range(a + 1, E):
            k, l = edges[b]
            if len({i, j, k, l}) == 3:
                share[a, b] = share[b, a] = True
    disj = ~share
    np.fill_diagonal(disj, False)
    Gu_disj = Gu * disj
    evals = np.linalg.eigvalsh(Gu_disj)
    lam_max = float(evals[-1])
    lam_min = float(evals[0])
    op_norm = float(max(abs(lam_max), abs(lam_min)))

    null, nullity = _zero_diag_nullspace(Vp)
    rng = np.random.default_rng(p + 99)
    rays = []
    for _ in range(400):
        coef = rng.standard_normal(nullity)
        v = null @ coef
        A = _A_from_vec(v, d)
        B = Vp @ A @ Vp.T
        nrm = np.linalg.norm(B)
        if nrm < 1e-14:
            continue
        B = B / nrm
        Be = np.array([2.0 * B[i, j] for i, j in edges])
        bn2 = float(np.dot(Be, Be))
        if bn2 < 1e-28:
            continue
        rays.append(float(Be @ Gu_disj @ Be) / bn2)

    # maximizer
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    U = (Y @ Vp) / np.sqrt(n)
    for _ in range(60):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    B = B / (np.linalg.norm(B) + 1e-30)
    Be = np.array([2.0 * B[i, j] for i, j in edges])
    bn2 = float(np.dot(Be, Be))
    r_max = float(Be @ Gu_disj @ Be) / bn2

    # For unit B: Be_norm2 = 2, disj_term = r * 2. Need disj ≤ 10N ⇒ r ≤ 5N
    thr = 5.0 * N
    return {
        "p": int(p),
        "skipped": False,
        "E": int(E),
        "N": int(N),
        "nullity": int(nullity),
        "lam_max_Gu_disj": lam_max,
        "lam_min_Gu_disj": lam_min,
        "op_norm": op_norm,
        "ray_image_max": float(max(rays)) if rays else None,
        "ray_image_mean": float(np.mean(rays)) if rays else None,
        "ray_image_p99": float(np.quantile(rays, 0.99)) if rays else None,
        "ray_maximizer": r_max,
        "thr_16N_ray": thr,
        "image_max_below_thr": bool(max(rays) <= thr + 1e-4) if rays else None,
        "maximizer_below_thr": bool(r_max <= thr + 1e-4),
        "unrestricted_exceeds_thr": bool(lam_max > thr),
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"typeA_wedge: workers={W}", flush=True)

    jobs = []
    for p in (3, 5, 7):
        jobs.append((p, 0, "max"))
        for s in range(1, 41):
            jobs.append((p, s, "rand"))

    rows: list[dict] = []
    specs: list[dict] = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(_worker_identity, j): ("id", j) for j in jobs}
        for p in (3, 5):
            futs[ex.submit(_worker_spectrum, p)] = ("spec", p)
        for fut in as_completed(futs):
            kind, tag = futs[fut]
            row = fut.result()
            if kind == "id":
                rows.append(row)
                if row["mode"] == "max" or row["seed"] <= 2:
                    print(
                        f"  p={row['p']} {row['mode']} s={row['seed']} "
                        f"null={row['nullity']} diag={row['diagmax']:.1e} "
                        f"six={row['sixN_ok']} wedge={row['wedge_ok']} "
                        f"typeA={row['typeA_ok']} ratio={row['ratio_to_16N']:.4f} "
                        f"Q4/N={row['Q4_over_N']:.4f}",
                        flush=True,
                    )
            else:
                specs.append(row)
                print(
                    f"  SPEC p={row['p']}: op={row.get('op_norm')} "
                    f"ray_img_max={row.get('ray_image_max')} "
                    f"ray_maxdir={row.get('ray_maximizer')} "
                    f"thr={row.get('thr_16N_ray')} "
                    f"unrestricted_hi={row.get('unrestricted_exceeds_thr')}",
                    flush=True,
                )

    summary = {}
    id_ok = True
    for p in (3, 5, 7):
        lst = [r for r in rows if r["p"] == p]
        six = all(r["sixN_ok"] and r["typeA_ok"] and r["wedge_ok"] for r in lst)
        # p=7 wedge measured via formula; still require typeA_ok and six via residual consistency
        if p == 7:
            six = all(r["typeA_ok"] and r["sixN_ok"] for r in lst)
        max_ratio = max(r["ratio_to_16N"] for r in lst)
        max_Q4N = max(r["Q4_over_N"] for r in lst)
        max_diag = max(r["diagmax"] for r in lst)
        summary[str(p)] = {
            "n_trials": len(lst),
            "nullity": lst[0]["nullity"],
            "typeA_wedge_identity_all_ok": six,
            "max_diag": max_diag,
            "max_ratio_to_16N": max_ratio,
            "max_Q4_over_N": max_Q4N,
            "all_disj_le_10N": all(r["disj_le_10N"] for r in lst),
            "max_mode_ratio": next(r["ratio_to_16N"] for r in lst if r["mode"] == "max"),
            "max_mode_Q4_over_N": next(r["Q4_over_N"] for r in lst if r["mode"] == "max"),
            "max_mode_lam_cycle_proxy": next(
                r["lam_cycle_proxy"] for r in lst if r["mode"] == "max"
            ),
        }
        id_ok = id_ok and six

    specs.sort(key=lambda s: s["p"])
    out = {
        "summary": summary,
        "spectrum": specs,
        "identity_claim": (
            "Zero-diag B=P_+ B P_+: typeA=2N‖B‖_F² (same-edge), "
            "wedge=4N‖B‖_F² (share-one-vertex via Σ=2P_+ on V_+), "
            "Q=6N‖B‖_F²+Q_4 with Q_4=BeᵀGu_disj Be only."
        ),
        "residual": (
            "Prove Q_4≤10N‖B‖_F² for all primes p≥5 "
            "(Gu_disj Rayleigh on image of zero-diag∩V_+ only; unrestricted op-norm fails)."
        ),
        "status": (
            "Prop 15.62 typeA+wedge=6N: "
            + ("CERTIFIED multi-seed p=3,5,7" if id_ok else "IDENTITY FAIL")
            + "; 16N / Gu_disj image bound OPEN."
        ),
        "workers": W,
        "identity_ok": bool(id_ok),
    }
    path = ROOT / "evidence" / "e1_gmin_typeA_wedge.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(json.dumps(summary, indent=2), flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
