#!/usr/bin/env python3
"""Compute Phi spectra for the dominant p=11 profile strata.

The script rebuilds exact pair Grams for k<=5, obtains k=6 by subtraction
from the verified total pair Gram, constructs an orthonormal edge-coordinate
basis of Z, and evaluates ``Phi = 4 A.T (G/N) A``.  It tests whether the
dominant strata individually cross the floor after the p=7 top-stratum value
5.6888... fell short.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import cupy as cp
import numpy as np

sys.path.insert(0, "/home/nick/quadratic-minmax-limit/src")
from minmax_quadratic import paley_conference_prime_power  # type: ignore


P = 11
N = P * P + 1
COUNTS = {"k1": 2772, "k3": 24200, "k4": 58080, "k5": 1306800, "k6": 36065260}


def accumulate_pair_gram(
    rows: np.ndarray, iu0: cp.ndarray, iu1: cp.ndarray, label: str
) -> np.ndarray:
    npair = int(iu0.size)
    acc = np.zeros((npair, npair), dtype=np.int64)
    chunk = int(os.environ.get("PAIR_GRAM_CHUNK", "50000"))
    nch = (len(rows) + chunk - 1) // chunk
    started = time.time()
    for ci, lo in enumerate(range(0, len(rows), chunk), start=1):
        yg = cp.asarray(rows[lo : lo + chunk], dtype=cp.int8)
        qg = (yg[:, iu0] * yg[:, iu1]).astype(cp.float32)
        gram = qg.T @ qg
        acc += cp.asnumpy(gram).astype(np.int64)
        del yg, qg, gram
        if ci % 10 == 0 or ci == nch:
            cp.get_default_memory_pool().free_all_blocks()
            print(f"{label}: {ci}/{nch}, elapsed={time.time()-started:.1f}s", flush=True)
    return acc


def z_edge_basis() -> np.ndarray:
    """Return A[pair,t]=B_t[i,j] for an orthonormal Frobenius basis of Z."""
    c = paley_conference_prime_power(P).astype(np.float64)
    eig, vec = np.linalg.eigh(c)
    vp = vec[:, eig > 1e-8]
    d = vp.shape[1]
    su = np.triu_indices(d)
    off = su[0] != su[1]

    diag_constraints = np.empty((N, len(su[0])), dtype=np.float64)
    for x in range(N):
        vx = vp[x]
        row = (vx[:, None] * vx[None, :])[su]
        row[off] *= np.sqrt(2.0)
        diag_constraints[x] = row
    _u, singular, vt = np.linalg.svd(diag_constraints, full_matrices=True)
    rank = int(np.sum(singular > 1e-9))
    kernel = vt[rank:].T

    iu = np.triu_indices(N, 1)
    vi = vp[iu[0]]
    vj = vp[iu[1]]
    edge_functionals = vi[:, su[0]] * vj[:, su[1]]
    edge_functionals[:, off] += vi[:, su[1][off]] * vj[:, su[0][off]]
    edge_functionals[:, off] /= np.sqrt(2.0)
    edge_basis = edge_functionals @ kernel

    gram = 2.0 * edge_basis.T @ edge_basis
    orth_error = float(np.max(np.abs(gram - np.eye(gram.shape[0]))))
    if orth_error > 2e-8:
        raise RuntimeError(f"Z edge basis orthogonality error {orth_error}")
    print(f"Z basis: shape={edge_basis.shape}, orth_error={orth_error:.3e}", flush=True)
    return edge_basis


def phi_spectrum(gram: np.ndarray, count: int, edge_basis: np.ndarray) -> dict[str, object]:
    ag = cp.asarray(edge_basis, dtype=cp.float32)
    mg = cp.asarray(gram, dtype=cp.float32) / np.float32(count)
    phi_gpu = np.float32(4.0) * (ag.T @ (mg @ ag))
    phi = cp.asnumpy(phi_gpu).astype(np.float64)
    del ag, mg, phi_gpu
    cp.get_default_memory_pool().free_all_blocks()
    phi = (phi + phi.T) / 2
    eig = np.linalg.eigvalsh(phi)
    return {
        "count_eps_plus": count,
        "lambda_min": float(eig[0]),
        "lambda_max": float(eig[-1]),
        "lambda_min_ge_6": bool(eig[0] >= 6 - 2e-5),
        "trace": float(np.trace(phi)),
        "trace2": float(np.sum(phi * phi)),
    }


def main() -> None:
    root = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    yall = np.load(root / "maxplus_p11_eps1.npy", mmap_mode="r")
    total = np.asarray(np.load(root / "G_pairmoment_p11.npy", mmap_mode="r"), dtype=np.int64)
    if yall.shape != (sum(COUNTS.values()), N):
        raise RuntimeError(f"unexpected Max+ shape {yall.shape}")

    iu = np.triu_indices(N, 1)
    iu0 = cp.asarray(iu[0].astype(np.int32))
    iu1 = cp.asarray(iu[1].astype(np.int32))
    grams: dict[str, np.ndarray] = {}
    lo = 0
    for label in ("k1", "k3", "k4", "k5"):
        hi = lo + COUNTS[label]
        grams[label] = accumulate_pair_gram(yall[lo:hi], iu0, iu1, label)
        lo = hi
    grams["k6"] = total - sum(grams.values())
    if not np.array_equal(sum(grams.values()), total):
        raise RuntimeError("stratum Gram reconstruction failed")

    edge_basis = z_edge_basis()
    report = {
        "p": P,
        "dimZ": edge_basis.shape[1],
        "spectra": {
            "k4": phi_spectrum(grams["k4"], COUNTS["k4"], edge_basis),
            "k5": phi_spectrum(grams["k5"], COUNTS["k5"], edge_basis),
            "k6": phi_spectrum(grams["k6"], COUNTS["k6"], edge_basis),
            "all": phi_spectrum(total, len(yall), edge_basis),
        },
    }
    output = root / "stratum_floor_p11.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
