#!/usr/bin/env python3
"""Best-response magnitude identity on Paley vs plus-I beaters.

For a maximizer with sign(Ax)=x, w_i=|(Ax)_i|, Phi=(1/2) sum w_i,
sum w_i^2 = x^T A^2 x. Conference uniquely has A^2=(n-1)I so
||Ax||_2 is constant on the cube. The Phi defect vs n*sqrt(n-1)/2
is exactly the L1/L2 flatness defect of w.

This script records that identity on the known beaters (serial: n<=26
MITM/local search, tiny matrices). Not a census.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (  # noqa: E402
    form_Q,
    paley_conference_matrix,
    paley_conference_prime_power,
    halfspace_boolean_vector,
    phi,
    phi_mitm,
    phi_local,
)


def lift(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1
    return b


def maximizer(A: np.ndarray, seed: int = 0) -> np.ndarray:
    n = A.shape[0]
    if n <= 20:
        npat = 1 << (n - 1)
        best_q, best_x = -1.0, None
        for xb in range(npat):
            x = np.ones(n, dtype=np.float64)
            for i in range(n - 1):
                x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
            q = abs(form_Q(A, x))
            if q > best_q:
                best_q, best_x = q, x.copy()
        return best_x
    rng = np.random.default_rng(seed)
    best_q, best_x = -1.0, None
    for s in range(200):
        x = rng.choice([-1.0, 1.0], size=n)
        Ax = A @ x
        qfull = float(x @ Ax)
        cur = abs(qfull) / 2.0
        improved = True
        while improved:
            improved = False
            for i in range(n):
                delta = -4.0 * x[i] * Ax[i]
                new_abs = abs(qfull + delta) / 2.0
                if new_abs > cur + 1e-12:
                    xi = x[i]
                    x[i] = -xi
                    Ax -= 2.0 * xi * A[:, i]
                    qfull += delta
                    cur = new_abs
                    improved = True
        if cur > best_q:
            best_q, best_x = cur, x.copy()
    return best_x


def stats(name: str, A: np.ndarray, x: np.ndarray) -> dict:
    n = A.shape[0]
    Ax = A @ x
    q = float(x @ Ax)
    if q < 0:
        x = -x
        Ax = -Ax
        q = -q
    w = Ax * x  # |(Ax)_i| iff x = sign(Ax)
    sign_ok = bool(np.all(w >= -1e-9))
    Phi = 0.5 * q
    if sign_ok:
        Phi = 0.5 * float(np.sum(w))
    s2 = float(x @ (A @ A) @ x)
    mu = np.sqrt(n - 1.0)
    mean_w = float(np.mean(w))
    rms_w = float(np.sqrt(s2 / n))
    var_w = float(np.mean((w - mean_w) ** 2))
    flatness = mean_w / rms_w if rms_w else 0.0
    Delta = A @ A - (n - 1) * np.eye(n)
    Dx = Delta @ x
    conference_phi = 0.5 * n * mu
    rec = {
        "name": name,
        "n": n,
        "Phi": Phi,
        "conference_phi": conference_phi,
        "gap": conference_phi - Phi,
        "op": float(np.max(np.abs(np.linalg.eigvalsh(A)))),
        "sign_aligned": sign_ok,
        "mean_w": mean_w,
        "rms_w": rms_w,
        "sqrt_n_1": mu,
        "var_w": var_w,
        "flatness_rho_proxy": flatness,
        "sum_w": float(np.sum(w)),
        "sum_w2": s2,
        "xT_Delta_x": float(x @ Dx),
        "||Delta||_F": float(np.linalg.norm(Delta, "fro")),
        "||Delta x||": float(np.linalg.norm(Dx)),
        "max_|w_i-mu|": float(np.max(np.abs(w - mu))),
        "w_min": float(np.min(w)),
        "w_max": float(np.max(w)),
    }
    return rec


def main() -> None:
    rows = []

    C10 = paley_conference_prime_power(3)
    x10 = halfspace_boolean_vector(3)
    rows.append(stats("Paley10", C10, x10))
    rows.append(stats("plusI_C5", lift(cycle5()), maximizer(lift(cycle5()))))

    C14 = paley_conference_matrix(13)
    rows.append(stats("Paley14", C14, maximizer(C14)))

    C18 = paley_conference_matrix(17)
    rows.append(stats("Paley18", C18, maximizer(C18)))

    C26 = paley_conference_prime_power(5)
    x26 = halfspace_boolean_vector(5)
    rows.append(stats("Paley26", C26, x26))

    n26 = ROOT / "evidence/ns_port_n26_undercut_20260912/stage01.npz"
    A26 = np.load(n26)["A"].astype(np.float64)
    # exact maximizer via mitm-scale local; confirm Phi=61
    xB = maximizer(A26, seed=1)
    rec = stats("plusI_n26", A26, xB)
    rec["phi_mitm"] = float(phi_mitm(A26))
    rec["Q_found"] = float(abs(form_Q(A26, xB)))
    rows.append(rec)

    out = ROOT / "evidence/beater_magnitude_20260919"
    out.mkdir(parents=True, exist_ok=True)
    (out / "magnitude.json").write_text(json.dumps(rows, indent=2))
    for r in rows:
        print(
            f"{r['name']:12s} n={r['n']:3d} Phi={r['Phi']:7.2f} "
            f"gap={r['gap']:7.2f} mean_w={r['mean_w']:.3f} "
            f"rms={r['rms_w']:.3f} mu={r['sqrt_n_1']:.3f} "
            f"var_w={r['var_w']:.3f} flat={r['flatness_rho_proxy']:.4f} "
            f"||D||F={r['||Delta||_F']:.2f} ||Dx||={r['||Delta x||']:.2f}"
        )


if __name__ == "__main__":
    main()
