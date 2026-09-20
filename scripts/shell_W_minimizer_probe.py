#!/usr/bin/env python3
"""Probe W(A) vs 2√2 Φ(A) on recorded exact small-n minimizers.

Not a census and not a proof. Architecture diagnostic for the astra
exchange lemma: among Φ-minimizers, does the energy-shell max in
NOTE_2026-09-05_GAUSSIAN_ENERGY_SHELL_UPPER (3) already sit under
2√2 Φ, and which (h,k) pair attains it?

Global Φ-minimizers cannot drop Φ by any Seidel flip (Φ(A')≥m_n).
If W>2√2 Φ on every recorded minimizer, the sufficient bound is the
wrong target and the split (8) on central shells is the live gap.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from minmax_quadratic import form_Q, paley_conference_matrix  # noqa: E402

KAPPA = 2.0 / math.pi


def cycle5() -> np.ndarray:
    b = np.ones((5, 5)) - np.eye(5)
    for i in range(5):
        b[i, (i + 1) % 5] = b[(i + 1) % 5, i] = -1.0
    return b


def plus_i(b: np.ndarray) -> np.ndarray:
    d = b + np.eye(b.shape[0])
    return np.block([[b, d], [d, -b]])


def all_states(n: int) -> np.ndarray:
    m = 1 << (n - 1)
    X = np.ones((m, n), dtype=np.float64)
    for i in range(n - 1):
        bit = (np.arange(m) >> i) & 1
        X[:, i + 1] = np.where(bit, 1.0, -1.0)
    return X


def shells(A: np.ndarray) -> dict:
    n = A.shape[0]
    X = all_states(n)
    Q = 0.5 * np.einsum("ia,ab,ib->i", X, A, X)
    # Q(-x)=Q(x); x_0=+1 already
    vals = {}
    for q, x in zip(Q, X):
        h = int(round(q))
        vals.setdefault(h, []).append(x)
    out = {h: np.stack(xs, axis=0) for h, xs in vals.items()}
    return out, float(np.max(np.abs(Q)))


def midpoint_params(A: np.ndarray):
    eig = np.linalg.eigvalsh(A)
    a = float(eig[-1])
    b = float(-eig[0])
    if a <= 0 or b <= 0:
        a = max(a, 1e-9)
        b = max(b, 1e-9)
    alpha = 0.5 * (a - b)
    mu = 0.5 * (a * a + b * b)
    return a, b, alpha, mu


def T_of(n, h, a, b, alpha, mu):
    v = n - 2.0 * KAPPA * alpha * h / mu
    t = KAPPA * (2.0 * h - alpha * n) / mu
    return v, t


def width_mc(S: np.ndarray, v: float, t: float, A: np.ndarray, rng, n_samp: int) -> float:
    """E max_{x in S} g^T x, g ~ N(0, vI + t A)."""
    n = A.shape[0]
    T = v * np.eye(n) + t * A
    # PSD clip for numerical midpoint error
    evals, evecs = np.linalg.eigh(T)
    evals = np.clip(evals, 0.0, None)
    cholish = evecs * np.sqrt(evals)
    G = rng.standard_normal((n_samp, n)) @ cholish.T
    return float(np.mean(np.max(G @ S.T, axis=1)))


def sample_Z(A: np.ndarray, n_samp: int, rng) -> np.ndarray:
    """Z with vec(Z) ~ N(0, C), C = I + κ (A⊗A)/μ at spectral midpoint α=0."""
    n = A.shape[0]
    a, b, alpha, mu = midpoint_params(A)
    # C = I + κ H/μ, H = A⊗A - α(A⊗I+I⊗A)
    I = np.eye(n)
    H = np.kron(A, A) - alpha * (np.kron(A, I) + np.kron(I, A))
    C = np.eye(n * n) + KAPPA * H / mu
    evals, evecs = np.linalg.eigh(C)
    evals = np.clip(evals, 0.0, None)
    fac = evecs * np.sqrt(evals)
    G = rng.standard_normal((n_samp, n * n)) @ fac.T
    return G.reshape(n_samp, n, n)


def EM_of(A: np.ndarray, n_samp: int = 4000, seed: int = 1) -> dict:
    """Monte Carlo E max_{x,y} |Q(x)-Q(y)+x^T Z y| with x_0=y_0=+1 doubled by signs."""
    n = A.shape[0]
    X = all_states(n)
    Q = 0.5 * np.einsum("ia,ab,ib->i", X, A, X)
    rng = np.random.default_rng(seed)
    Zs = sample_Z(A, n_samp, rng)
    # For each sample: max_{i,j} |Q_i - Q_j + X_i Z X_j|
    # XZ : (n_samp, m, n) = Z @ X.T -> (n_samp, n, m); then X @ that
    Xt = X.T
    maxes = np.empty(n_samp)
    for s in range(n_samp):
        cross = X @ (Zs[s] @ Xt)
        M = np.abs(Q[:, None] - Q[None, :] + cross)
        maxes[s] = float(np.max(M))
    return {
        "E_M": round(float(np.mean(maxes)), 4),
        "E_M_std": round(float(np.std(maxes) / math.sqrt(n_samp)), 4),
        "E_M_minus_target": None,  # filled by caller
        "n_samp_M": n_samp,
    }


def W_of(A: np.ndarray, n_samp: int = 8000, seed: int = 0) -> dict:
    sh, phi = shells(A)
    n = A.shape[0]
    a, b, alpha, mu = midpoint_params(A)
    rng = np.random.default_rng(seed)
    hs = sorted(sh)
    widths = {}
    for h in hs:
        for k in hs:
            v, t = T_of(n, k, a, b, alpha, mu)
            widths[(h, k)] = width_mc(sh[h], v, t, A, rng, n_samp)
    best = None
    best_val = -1e18
    rows = []
    for h in hs:
        for k in hs:
            val = abs(h - k) + widths[(h, k)] + widths[(k, h)]
            rows.append(
                {
                    "h": h,
                    "k": k,
                    "abs_hk": abs(h - k),
                    "w_h_Tk": round(widths[(h, k)], 4),
                    "w_k_Th": round(widths[(k, h)], 4),
                    "shell": round(val, 4),
                    "n_h": int(sh[h].shape[0]),
                    "n_k": int(sh[k].shape[0]),
                }
            )
            if val > best_val:
                best_val = val
                best = rows[-1]
    target = 2.0 * math.sqrt(2.0) * phi
    return {
        "n": n,
        "Phi": phi,
        "target_2sqrt2_Phi": round(target, 4),
        "W": round(best_val, 4),
        "W_minus_target": round(best_val - target, 4),
        "W_over_n32": round(best_val / (n ** 1.5), 4),
        "target_over_n32": round(target / (n ** 1.5), 4),
        "attained": best,
        "shell_sizes": {str(h): int(sh[h].shape[0]) for h in hs},
        "spectrum_ab": [round(a, 4), round(b, 4)],
        "midpoint_alpha_mu": [round(alpha, 4), round(mu, 4)],
        "top_pairs": sorted(rows, key=lambda r: -r["shell"])[:8],
    }


def main() -> None:
    mats = []
    c5 = cycle5()
    mats.append(("C5_m5", c5))
    mats.append(("Paley6_m6", paley_conference_matrix(5)))
    mats.append(("plusI_C5_m10", plus_i(c5)))
    b13 = np.loadtxt(ROOT / "evidence" / "coherent_BI_lift_20260919" / "B13.txt")
    mats.append(("B13_m13", b13))
    recs = []
    for name, A in mats:
        print("BEGIN", name, "n=", A.shape[0], flush=True)
        rec = W_of(A)
        ns_m = 300 if A.shape[0] >= 13 else 4000
        em = EM_of(A, n_samp=ns_m)
        em["E_M_minus_target"] = round(em["E_M"] - rec["target_2sqrt2_Phi"], 4)
        rec.update(em)
        rec["name"] = name
        recs.append(rec)
        print(
            f"  Phi={rec['Phi']} W={rec['W']} target={rec['target_2sqrt2_Phi']} "
            f"W-tgt={rec['W_minus_target']} E_M={rec['E_M']} E_M-tgt={rec['E_M_minus_target']} "
            f"attained={rec['attained']}",
            flush=True,
        )
    out = {
        "note": "Energy-shell W vs 2√2 Φ on recorded exact minimizers. Not a theorem.",
        "records": recs,
    }
    path = ROOT / "evidence" / "case_b_gamma_p5" / "shell_W_minimizer_probe.json"
    # keep with campaign receipts but this is the Gaussian-shell probe
    path = ROOT / "evidence" / "shell_W_minimizer_probe_20260920.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps({r["name"]: {k: r[k] for k in r if k != "top_pairs"} for r in recs}, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
