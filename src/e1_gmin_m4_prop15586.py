#!/usr/bin/env python3
"""
Prop 15.586 — Max+ Gram reduction of the bi-tight floor.

Reduces lambda_min(Phi) from a 4-point tensor on Z to the spectrum of an
INTEGER, manifestly PSD Gram matrix built from 2-point data on Max+.

Does **not** flip phi_F_ge_6 / residual_ii / type_I / e1 / L / Aut-Schur /
Gsum / pairing.  Lemma D stays True.  This is a REDUCTION, not the floor:
`phi_F_ge_6_proved_general` must stay False.  Soft-close forbidden.

============================================================================
SETUP.  n = p^2+1, C = Paley conference, P = (I + C/p)/2 the projector onto
the +p eigenspace (d = n/2).  Max+ = {y in {+-1}^n : Cy = py}, N = |Max+|;
every y in Max+ has Py = y and y_i^2 = 1.
  Z = {B sym : CB = pB, diag B = 0},  Phi(B) = E_{y in Max+}[(y^T B y)^2].
Floor (GOAL.md leftover 1) is lambda_min(Phi) >= 6.

Theorem A — PROVED.  Uniform Z-perp component.
  Z-perp inside {B : PBP = B} is span{pi_i pi_i^T}, pi_i = P e_i.
  <y y^T, pi_i pi_i^T> = (y^T P e_i)^2 = ((Py)_i)^2 = y_i^2 = 1
  for every i and every y in Max+.  So proj_{Zperp}(y y^T) = R is the SAME
  element for all y in Max+, with R = 2 sum_i pi_i pi_i^T and
  ||R||^2 = 4 tr(P^2) = 4 tr P = 2n.

Theorem B — PROVED.  Gram reduction.
  proj_Z(y y^T) = y y^T - R, hence
      <proj_Z(y_a y_a^T), proj_Z(y_b y_b^T)> = <y_a,y_b>^2 - 2n.
  So with Ghat_ab = <y_a,y_b>^2 - 2n (integer, PSD by construction),
      spec(Phi) = nonzero spec(Ghat / N),   #nonzero = dim Z.
  This is a genuine Gram matrix -- unlike G_{u,disj}, which the package
  Caveat 1 forbids treating as a Gram because it has negative eigenvalues.

Theorem C — PROVED.  Closed forms.
  dim Z = n(n-6)/8,  tr Phi = n(n-2),  tr K = -4n for K = 8I - Phi,
  dim (pair-span of Max+) = dim span{1, y_i y_j} = n(n-6)/8 + 1 = dim Z + 1.
  The pair-span bound is the object the residual-(ii) 0-1 pair-span
  classification (15.274 D/F, 15.585 B/C) reasons about, so leftover 1 and
  leftover 2 are controlled by the SAME matrix Ghat.

Theorem D — PROVED.  Wick baseline is exactly 8.
  PB = BP = B and tr B = 0 for B in Z, and E[y y^T] = I + C/p = 2P, so the
  Gaussian/Wick value is 2 tr(B G B G) = 8||B||^2.  Hence
      floor lambda_min(Phi) >= 6   <==>   lambda_max(K) <= 2,
  which is exactly the package's <delta,psi> <= 2, now as an explicit
  operator K = 8I - Phi with tr K = -4n.

Theorem E — CERTIFIED p=5,7 (finite).  Exact spectra:
  p=5: spec(Phi) = {80/13 (mult 26), 144/13 (26), 176/13 (13)}
  p=7: spec(Phi) = {3072/409 (50), 3360/409 (100), 3648/409 (50),
                    4032/409 (50), 4320/409 (25)}
  lambda_max(K) = 48/n EXACTLY at p=5 (= 24/13), and 200/409 < 48/50 at p=7.
  So the candidate lambda_* = 8(n-6)/n is tight at p=5 and slack at p=7:
  the floor binds only at the smallest prime, and p=5 is a finite check.

OPEN.  Theorems A-D are general-p; E is finite.  Nothing here proves the
floor for p >= 11 -- Max+ enumeration at p=11 is 2^61 (nullity 61) and the
"full" family is unclassified.  phi_F_ge_6_proved_general stays False.

Writes evidence/e1_gmin_m4_prop15586.json
"""
from __future__ import annotations

import functools
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402

TOL = 1e-8


def n_of(p: int) -> int:
    return p * p + 1


def dim_Z_closed(p: int) -> int:
    n = n_of(p)
    return n * (n - 6) // 8


def maxplus(p: int) -> np.ndarray:
    """Max+ = {y in {+-1}^n : Cy = py}.  Cached; else exhaustive nullspace sweep."""
    for path in (Path(f"/tmp/e1_p{p}/maxplus.npy"), Path(f"/tmp/maxplus_p{p}.npy")):
        if path.is_file():
            return np.load(path).astype(np.float64)
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    _, s, vt = np.linalg.svd(C - p * np.eye(n))
    rank = int((s > 1e-8).sum())
    nul = n - rank
    nb = vt[rank:].T
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(50000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError(f"p={p}: no free coordinate set")
    A = nb @ np.linalg.inv(nb[coords])
    try:
        import cupy as cp  # noqa: F401
        xp, Ag = cp, cp.asarray(A)
    except Exception:
        xp, Ag = np, A
    total, bits = 1 << nul, xp.arange(nul, dtype=xp.int64)
    ones = xp.ones(n, dtype=xp.float64)
    found = []
    bs = min(total, 1 << 21)                      # batched: never a 1-core sweep (F17)
    for start in range(0, total, bs):
        idxs = xp.arange(start, min(start + bs, total), dtype=xp.int64)[:, None]
        free = (((idxs >> bits) & 1) * 2 - 1).astype(xp.float64)
        err = xp.abs(xp.abs(free @ Ag.T) - 1.0) @ ones     # matvec, not axis-reduce
        e = err.get() if xp is not np else err
        hit = np.flatnonzero(e < 1e-6)
        if hit.size:
            b = (start + hit)[:, None]
            fr = (((b >> np.arange(nul, dtype=np.int64)) & 1) * 2 - 1).astype(float)
            found.append(np.sign(fr @ A.T))
    Y = np.unique(np.concatenate(found, 0).astype(np.int8), axis=0).astype(np.float64)
    Path(f"/tmp/e1_p{p}").mkdir(exist_ok=True)
    np.save(f"/tmp/e1_p{p}/maxplus.npy", Y)
    return Y


@functools.lru_cache(maxsize=None)
def _setup(p: int):
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = maxplus(p)
    P = (np.eye(n) + C / p) / 2.0
    return C, n, Y, len(Y), P


def theorem_A_uniform_Zperp(primes=(5, 7)) -> dict:
    """<y y^T, pi_i pi_i^T> = 1 for every i, y  =>  proj_Zperp(y y^T) is y-free."""
    ok, worst = True, 0.0
    for p in primes:
        C, n, Y, N, P = _setup(p)
        dev = np.abs((Y @ P) ** 2 - 1.0).max()      # (y^T pi_i)^2 = 1 for all i,a
        R = 2.0 * (P @ P.T)                          # 2 sum_i pi_i pi_i^T = 2 P P^T = 2P
        nr = abs(float(np.sum(R * R)) - 2 * n)
        worst = max(worst, dev, nr)
        ok = ok and dev < TOL and nr < TOL * n
    return {"proved": bool(ok), "max_dev": worst, "primes": list(primes)}


def theorem_B_gram_reduction(primes=(5, 7)) -> dict:
    """spec(Phi) == nonzero spec(Ghat/N), Ghat_ab = <y_a,y_b>^2 - 2n."""
    ok, worst, out = True, 0.0, {}
    for p in primes:
        C, n, Y, N, P = _setup(p)
        lam_phi = _phi_spectrum(p)
        G = Y @ Y.T
        Gh = G * G - 2.0 * n
        lam = np.linalg.eigvalsh((Gh + Gh.T) / 2) / N
        nz = np.sort(lam[np.abs(lam) > 1e-7])
        good = len(nz) == dim_Z_closed(p)
        dev = float(np.abs(nz - np.sort(lam_phi)).max()) if good else float("inf")
        worst = max(worst, dev)
        ok = ok and good and dev < 1e-6
        out[p] = {"n_nonzero": len(nz), "dim_Z": dim_Z_closed(p), "max_dev": dev}
    return {"proved": bool(ok), "max_dev": worst, "per_prime": out}


@functools.lru_cache(maxsize=None)
def _phi_spectrum(p: int) -> np.ndarray:
    """Eigenvalues of Phi on Z, built directly (independent of the Ghat route).

    Disk-cached: xdist workers are separate processes, so without this each
    worker would rebuild the same spectrum (serial waste on an 88-core box).
    """
    cache = Path(f"/tmp/e1_p{p}/phi_spectrum.npy")
    if cache.is_file():
        return np.load(cache)
    C, n, Y, N, P = _setup(p)
    ev, EV = np.linalg.eigh(C)
    Vp = EV[:, ev > 1e-8]
    d = Vp.shape[1]
    cols = [np.eye(d)[a].reshape(d, 1) @ np.eye(d)[a].reshape(1, d) for a in range(d)]
    r = 1.0 / np.sqrt(2.0)
    for a in range(d):
        for b in range(a + 1, d):
            E = np.zeros((d, d))
            E[a, b] = E[b, a] = r
            cols.append(E)
    S = np.stack([c.ravel() for c in cols], 1)
    W = np.stack([np.outer(Vp[i], Vp[i]).ravel() for i in range(n)], 0)
    _, sv, vt = np.linalg.svd(W @ S)
    rk = int((sv > 1e-9).sum())
    Cm = vt[rk:].T
    Zc = Y @ Vp
    U = ((Zc[:, :, None] * Zc[:, None, :]).reshape(len(Y), -1) @ S) @ Cm
    Phi = U.T @ U / len(Y)
    lam = np.linalg.eigvalsh((Phi + Phi.T) / 2)
    cache.parent.mkdir(exist_ok=True)
    np.save(cache, lam)
    return lam


def theorem_C_closed_forms(primes=(5, 7)) -> dict:
    """dim Z = n(n-6)/8, tr Phi = n(n-2), tr K = -4n, dim pair-span = dim Z + 1."""
    ok, out = True, {}
    for p in primes:
        C, n, Y, N, P = _setup(p)
        lam = _phi_spectrum(p)
        iu = np.triu_indices(n, 1)
        Wp = np.concatenate([np.ones((N, 1)), Y[:, iu[0]] * Y[:, iu[1]]], axis=1)
        sv = np.linalg.svd(Wp, compute_uv=False)
        rank = int((sv > sv[0] * 1e-8).sum())
        c = {
            "dim_Z": len(lam) == dim_Z_closed(p),
            "tr_Phi": abs(lam.sum() - n * (n - 2)) < 1e-6,
            "tr_K": abs((8.0 - lam).sum() + 4 * n) < 1e-6,
            "pairspan": rank == dim_Z_closed(p) + 1,
        }
        out[p] = {**c, "rank_pairspan": rank, "tr_Phi_val": float(lam.sum())}
        ok = ok and all(c.values())
    return {"proved": bool(ok), "per_prime": out}


def theorem_D_wick_baseline(primes=(5, 7)) -> dict:
    """E[y y^T] = I + C/p = 2P, so Wick value is 8 and floor <=> lam_max(K) <= 2."""
    ok, worst = True, 0.0
    for p in primes:
        C, n, Y, N, P = _setup(p)
        dev = float(np.abs(Y.T @ Y / N - (np.eye(n) + C / p)).max())
        worst = max(worst, dev)
        ok = ok and dev < TOL
    return {"proved": bool(ok), "max_dev": worst, "wick_value": 8.0,
            "floor_equivalent": "lambda_max(8I - Phi) <= 2"}


def theorem_E_exact_spectra() -> dict:
    """Finite certificate: exact rational spectra and lam_max(K) at p=5,7."""
    want = {
        5: {Fraction(80, 13): 26, Fraction(144, 13): 26, Fraction(176, 13): 13},
        7: {Fraction(3072, 409): 50, Fraction(3360, 409): 100,
            Fraction(3648, 409): 50, Fraction(4032, 409): 50,
            Fraction(4320, 409): 25},
    }
    ok, out = True, {}
    for p, exp in want.items():
        lam = _phi_spectrum(p)
        got = {}
        for v in lam:
            f = Fraction(v).limit_denominator(10 ** 7)
            got[f] = got.get(f, 0) + 1
        match = got == exp
        n = n_of(p)
        kmax = 8.0 - float(lam.min())
        out[p] = {"spectrum_matches": match,
                  "lam_max_K": kmax,
                  "equals_48_over_n": abs(kmax - 48.0 / n) < 1e-9,
                  "floor_ge_6": bool(lam.min() >= 6.0 - 1e-9)}
        ok = ok and match and lam.min() >= 6.0 - 1e-9
    return {"proved": bool(ok), "per_prime": out}


def floor_gram_reduction_proved() -> bool:
    """The REDUCTION (A-D) is general-p.  This is NOT the floor."""
    return bool(
        theorem_A_uniform_Zperp()["proved"]
        and theorem_B_gram_reduction()["proved"]
        and theorem_C_closed_forms()["proved"]
        and theorem_D_wick_baseline()["proved"]
    )


def phi_F_ge_6_still_open() -> bool:
    """Guard: this prop must NOT close leftover 1."""
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
    return not phi_F_ge_6_proved_general()


def main() -> dict:
    from io_atomic import write_json_atomic
    out = {
        "prop": "15.586",
        "A_uniform_Zperp": theorem_A_uniform_Zperp(),
        "B_gram_reduction": theorem_B_gram_reduction(),
        "C_closed_forms": theorem_C_closed_forms(),
        "D_wick_baseline": theorem_D_wick_baseline(),
        "E_exact_spectra": theorem_E_exact_spectra(),
        "reduction_proved": floor_gram_reduction_proved(),
        "phi_F_ge_6_still_open": phi_F_ge_6_still_open(),
    }
    print(f"Prop 15.586  reduction_proved={out['reduction_proved']}  "
          f"floor still open={out['phi_F_ge_6_still_open']}", flush=True)
    write_json_atomic(ROOT / "evidence" / "e1_gmin_m4_prop15586.json", out)
    return out


if __name__ == "__main__":
    main()
