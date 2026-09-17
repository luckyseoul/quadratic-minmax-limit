#!/usr/bin/env python3
"""Independent numerical verification of the 2026-09-15 explicit lower-bound estimates.

Target: ``evidence/NOTE_2026-09-15_EXPLICIT_ORIGINAL_LOWER_BOUND.md``.

Scope. Corroboration only; no all-orders statement is proved or refuted here.
The all-orders argument is analytic; this script re-derives, for an explicit
complete symmetric zero-diagonal signing ``A``, every object in the note's
chain and checks its displayed estimates on explicit families:

  * ``M=A/sqrt(n)``, ``S=sign(M)`` (kernel `+1`), ``|M|``, ``ell=tr|M|/n``,
    ``q=(n-1)/n``, ``d=q+1-2*ell``;
  * the rounding frame ``R`` of Section 2, the exact sign-Gaussian
    covariance ``C_X=(2/pi)*arcsin(R)`` (entrywise), ``v_i=(S C_X S)_ii``
    and ``E F_i^2=(M C_X M)_ii``;
  * estimates (4), (5), (6), (7), the diagonal-variance estimate, and the
    field-variance estimate -- these are deterministic in the signings and
    are checked exactly up to float64 rounding;
  * the Gaussianization estimate (8) and the update bound (9) -- checked by
    Monte Carlo, since ``E|F_i|`` has no elementary exact form (for (9) the
    exact ``Phi(A)`` is additionally used on exhaustive small orders).

Families: all signings of order ``n <= 7`` (exhaustive), Paley conference
signings of order ``q+1`` for prime ``q = 1 mod 4``, conference signings
with forced local flips, uniform-random signings, and the all-ones
off-diagonal signing. Self-contained: numpy + stdlib only, so it can be
replayed unchanged on other mesh nodes.

Usage examples:
  python3 verify_estimates.py --exhaustive 6 --out result_estimates_ex6.json
  python3 verify_estimates.py --families --mc-samples 200000 --out result_families.json
  python3 verify_estimates.py --exhaustive 7 --workers 80 --out result_estimates_ex7.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import socket
import sys
import time
from pathlib import Path

import numpy as np

KAPPA = 2.0 / math.pi
TOL = 1e-9


# --------------------------------------------------------------------------
# signings
# --------------------------------------------------------------------------
def num_entries(n: int) -> int:
    return n * (n - 1) // 2


def decode_signing(n: int, idx: int) -> np.ndarray:
    """Signing from a bit pattern over the strict upper triangle (row-major)."""
    i, j = np.triu_indices(n, 1)
    bits = (np.asarray(idx, dtype=np.int64) >> np.arange(num_entries(n))) & 1
    vals = np.where(bits == 1, 1.0, -1.0)
    A = np.zeros((n, n))
    A[i, j] = vals
    A[j, i] = vals
    return A


def paley_conference(q: int) -> np.ndarray:
    """Symmetric conference signing C of order n=q+1 for prime q = 1 mod 4.

    Index 0 is the point at infinity. Self-contained reimplementation
    (validated against C^2 = q I and against the repo builder where present).
    """
    if q % 4 != 1:
        raise ValueError("need q = 1 mod 4")
    for d in range(2, int(math.isqrt(q)) + 1):
        if q % d == 0:
            raise ValueError(f"{q} is not prime (use the repo GF builder for prime powers)")
    leg = np.zeros(q, dtype=np.float64)
    inner = np.zeros(q, dtype=np.float64)
    for u in range(1, q):
        chi = pow(u, (q - 1) // 2, q)
        leg[u] = 1.0 if chi == 1 else -1.0
    n = q + 1
    C = np.zeros((n, n))
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    inner = np.zeros((q, q))
    for u in range(q):
        for v in range(q):
            inner[u, v] = leg[(u - v) % q]
    C[1:, 1:] = inner
    return C


def all_ones_offdiag(n: int) -> np.ndarray:
    A = np.ones((n, n))
    np.fill_diagonal(A, 0.0)
    return A


# --------------------------------------------------------------------------
# per-signing analysis
# --------------------------------------------------------------------------
def build_frame(S: np.ndarray) -> tuple[np.ndarray, list, list]:
    """The Section 2 rounding frame R (correlation matrix, 0 <= R <= 4I)."""
    n = S.shape[0]
    diagS = np.diag(S)
    I_idx = [i for i in range(n) if abs(diagS[i]) <= 0.5]
    B_idx = [i for i in range(n) if abs(diagS[i]) > 0.5]
    R = np.eye(n)
    P = np.eye(n) + S
    if I_idx:
        u = 1.0 / np.sqrt(1.0 + diagS[I_idx])
        R[np.ix_(I_idx, I_idx)] = P[np.ix_(I_idx, I_idx)] * np.outer(u, u)
    return R, I_idx, B_idx


def monte_carlo_f(M: np.ndarray, R: np.ndarray, samples: int, seed: int,
                  chunk: int = 50000) -> tuple[float, float]:
    """Estimate f=(1/n) sum_i E|F_i|, F=MX, X=sign(G), G~N(0,R)."""
    n = M.shape[0]
    lam, Q = np.linalg.eigh(R)
    lam = np.clip(lam, 0.0, None)
    L = Q * np.sqrt(lam)
    rng = np.random.default_rng(seed)
    total = 0.0
    total_sq = 0.0
    count = 0
    remaining = samples
    while remaining > 0:
        m = min(chunk, remaining)
        remaining -= m
        Z = rng.standard_normal((m, n))
        G = Z @ L.T
        X = np.sign(G)
        X[X == 0] = 1.0
        F = X @ M.T
        row = np.abs(F).mean(axis=1)
        total += float(row.sum())
        total_sq += float((row * row).sum())
        count += m
    f = total / count
    var = max(total_sq / count - f * f, 0.0)
    return f, math.sqrt(var / count)


def phi_exact(A: np.ndarray) -> float:
    """Exact max|Q_A(x)| over the full cube (small n only)."""
    n = A.shape[0]
    m = 1 << n
    idx = np.arange(m, dtype=np.uint32)
    X = np.empty((m, n))
    for j in range(n):
        X[:, j] = 1.0 - 2.0 * ((idx >> j) & 1)
    Q = 0.5 * np.einsum("ij,jk,ik->i", X, A, X)
    return float(np.abs(Q).max())


def analyze(A: np.ndarray, mc_samples: int = 0, seed: int = 12345,
            update_check: bool = False) -> dict:
    """All deterministic estimate checks for one signing, plus optional MC."""
    n = A.shape[0]
    M = A / math.sqrt(n)
    lam, Q = np.linalg.eigh(M)
    sgn = np.where(lam >= 0.0, 1.0, -1.0)
    S = (Q * sgn) @ Q.T
    ell = float(np.abs(lam).sum() / n)
    q = (n - 1.0) / n
    d = q + 1.0 - 2.0 * ell
    out: dict = {"n": n, "d": d, "ell": ell}

    # identity: ||M-S||_F^2 = n d  (Section 2)
    out["frob_identity_err"] = abs(float(np.sum((M - S) ** 2)) - n * d)
    # sum_i S_ii^2 <= n d  (used in the proof of (5))
    diag_sq = float(np.sum(np.diag(S) ** 2))
    out["sdiag_over_nd"] = (diag_sq / (n * d)) if d > 0 else float("nan")
    out["sdiag_le_nd"] = bool(diag_sq <= n * d + 1e-9)

    R, I_idx, B_idx = build_frame(S)
    eigR = np.linalg.eigvalsh(R)
    # Tolerance note: entries of R can equal +/-1 exactly, so the entrywise
    # (2/pi) arcsin evaluation loses ~1e-8 (arcsin derivative blows up at 1);
    # the PSD/diagonal checks below therefore use a 1e-6 float64 tolerance.
    TOL_PSD = 1e-6
    checks = {}
    # (4) R is a correlation matrix with R <= 4I
    checks["c4"] = {
        "diag_err": float(np.abs(np.diag(R) - 1.0).max()),
        "psd_slack": float(-eigR[0]),
        "upper_slack": float(4.0 - eigR[-1]),
        "ok": bool(np.abs(np.diag(R) - 1.0).max() <= TOL_PSD
                   and -eigR[0] <= TOL_PSD and 4.0 - eigR[-1] >= -TOL_PSD),
    }
    # (5) ||R-(I+S)||_F <= 4 sqrt(n d)
    lhs5 = float(np.linalg.norm(R - (np.eye(n) + S)))
    rhs5 = 4.0 * math.sqrt(n * d)
    checks["c5"] = {"lhs": lhs5, "rhs": rhs5,
                    "ratio": (lhs5 / rhs5) if rhs5 > 0 else None,
                    "ok": bool(lhs5 <= rhs5 + 1e-9 * max(1.0, rhs5))}

    # exact sign-Gaussian covariance
    C = (2.0 / math.pi) * np.arcsin(np.clip(R, -1.0, 1.0))
    checks["cX"] = {
        "diag_err": float(np.abs(np.diag(C) - 1.0).max()),
        "psd_slack": float(-np.linalg.eigvalsh(C)[0]),
        "ok": bool(np.abs(np.diag(C) - 1.0).max() <= TOL_PSD
                   and -np.linalg.eigvalsh(C)[0] <= TOL_PSD),
    }
    # (6) e := E Q_M(X)/n >= kappa ell/2 - 2 kappa sqrt(d) - 2(1-kappa)/sqrt(n)
    e_val = float(np.trace(M @ C) / (2.0 * n))
    rhs6 = KAPPA * ell / 2.0 - 2.0 * KAPPA * math.sqrt(d) - 2.0 * (1.0 - KAPPA) / math.sqrt(n)
    checks["c6"] = {"e": e_val, "rhs": rhs6,
                    "ratio": (rhs6 / e_val) if (rhs6 > 0 and e_val > 0) else None,
                    "active": bool(rhs6 > 0.0),
                    "ok": bool(e_val >= rhs6 - 1e-9)}
    # (7) ||C_X-(I+kappa S)||_F
    lhs7 = float(np.linalg.norm(C - (np.eye(n) + KAPPA * S)))
    a_n = 1.0 + 1.0 / math.sqrt(n) + 1.0 / n
    rhs7 = (4.0 * KAPPA + 5.0 * (1.0 - KAPPA) * a_n) * math.sqrt(n * d) \
        + (1.0 - KAPPA) / math.sqrt(n)
    checks["c7"] = {"lhs": lhs7, "rhs": rhs7, "ratio": lhs7 / rhs7,
                    "ok": bool(lhs7 <= rhs7 + 1e-9 * max(1.0, rhs7))}
    # diagonal variance estimate
    v = np.diag(S @ C @ S)
    lhs_v = float(np.abs(v - 1.0).mean())
    rhs_v = (5.0 + 5.0 * (1.0 - KAPPA) * (1.0 / math.sqrt(n) + 1.0 / n)) * math.sqrt(d) \
        + (1.0 - KAPPA) / n
    checks["cv"] = {"lhs": lhs_v, "rhs": rhs_v,
                    "ratio": (lhs_v / rhs_v) if rhs_v > 0 else None,
                    "ok": bool(lhs_v <= rhs_v + 1e-9 * max(1.0, rhs_v))}
    # field-variance estimate
    F2 = np.diag(M @ C @ M)
    lhs_fv = float(np.abs(np.sqrt(F2) - np.sqrt(np.maximum(v, 0.0))).mean())
    rhs_fv = 2.0 * math.sqrt(d)
    checks["cfv"] = {"lhs": lhs_fv, "rhs": rhs_fv,
                     "ratio": (lhs_fv / rhs_fv) if rhs_fv > 0 else None,
                     "ok": bool(lhs_fv <= rhs_fv + 1e-9 * max(1.0, rhs_fv))}
    # b_n of (8)
    b_n = (7.0 + 5.0 * (1.0 - KAPPA) * (1.0 / math.sqrt(n) + 1.0 / n)) * math.sqrt(d) \
        + (1.0 - KAPPA) / n
    extras = {"b_n": b_n, "sqrt_kappa": math.sqrt(KAPPA),
              "rhs8": math.sqrt(KAPPA) * (1.0 - b_n),
              "rhs8_7sqrtd": math.sqrt(KAPPA) * (1.0 - 7.0 * math.sqrt(d))}
    out["checks"] = checks
    out["extras"] = extras

    if mc_samples > 0:
        f_est, f_se = monte_carlo_f(M, R, mc_samples, seed)
        extras["f_est"] = f_est
        extras["f_se"] = f_se
        extras["e4_empirical"] = max(0.0, extras["rhs8"] - f_est)
        extras["f_over_rhs8"] = f_est / extras["rhs8"]
        extras["mc_samples"] = mc_samples
    if update_check:
        alpha = phi_exact(A) / n ** 1.5
        rng = np.random.default_rng(seed + 1)
        lam2, Q2 = np.linalg.eigh(R)
        lam2 = np.clip(lam2, 0.0, None)
        L2 = Q2 * np.sqrt(lam2)
        G = rng.standard_normal((256, n)) @ L2.T
        X = np.sign(G)
        X[X == 0] = 1.0
        Fm = X @ M.T
        Y = np.sign(Fm)
        Y[Y == 0] = 1.0
        p = 0.1
        Z = (1.0 - p) * X + p * Y
        Qz = 0.5 * np.einsum("ij,jk,ik->i", Z, M, Z)
        extras["update_alpha"] = alpha
        # the Boolean bound used in the note is |Q_M(z)| <= n alpha, i.e.
        # <= Phi(A)/sqrt(n) for every z in [-1,1]^n
        extras["update_max_abs_ratio"] = float(np.abs(Qz).max() / (n * alpha))
        # exact expansion check of (9) on a random subset of states
        e_xi = (1.0 - p) ** 2 * 0.5 * np.einsum("ij,jk,ik->i", X, M, X) \
            + p * (1.0 - p) * np.abs(Fm).sum(axis=1) \
            + p ** 2 * 0.5 * np.einsum("ij,jk,ik->i", Y, M, Y)
        extras["update_expansion_err"] = float(np.abs(e_xi - Qz).max())
    return out


# --------------------------------------------------------------------------
# aggregators
# --------------------------------------------------------------------------
def merge_check_agg(agg: dict, res: dict, where: dict) -> None:
    """Track max ratio / violations per check for one analyzed signing."""
    for key, val in res["checks"].items():
        slot = agg.setdefault(key, {"max_ratio": None, "where": None,
                                    "violations": 0, "inactive": 0})
        if not val.get("ok", False):
            slot["violations"] += 1
        ratio = val.get("ratio")
        if ratio is None:
            slot["inactive"] += 1
            continue
        if slot["max_ratio"] is None or ratio > slot["max_ratio"]:
            slot["max_ratio"] = ratio
            slot["where"] = dict(where)
            slot["where"].update({"d": res["d"]})
    if not res["sdiag_le_nd"]:
        agg.setdefault("sdiag", {"violations": 0, "max_ratio": None,
                                 "where": None, "inactive": 0})["violations"] += 1
    if res["frob_identity_err"] > 1e-8:
        agg.setdefault("frob_identity", {"violations": 0, "max_ratio": None,
                                         "where": None, "inactive": 0})["violations"] += 1
    ex = res.get("extras", {})
    if "update_max_abs_ratio" in ex:
        slot = agg.setdefault("update9", {"max_ratio": None, "where": None,
                                          "violations": 0, "inactive": 0})
        r = ex["update_max_abs_ratio"]
        if slot["max_ratio"] is None or r > slot["max_ratio"]:
            slot["max_ratio"] = r
            slot["where"] = dict(where)
        if r > 1.0 + 1e-9:
            slot["violations"] += 1


def _exhaustive_worker(payload):
    n, lo, hi, mc_samples, update_check, seed = payload
    agg: dict = {}
    count = 0
    for idx in range(lo, hi):
        A = decode_signing(n, idx)
        res = analyze(A, mc_samples=0, seed=seed,
                      update_check=update_check)
        merge_check_agg(agg, res, {"idx": idx})
        count += 1
    return {"count": count, "agg": agg}


# --------------------------------------------------------------------------
# modes
# --------------------------------------------------------------------------
def run_exhaustive(n: int, workers: int, mc_samples: int,
                   update_check: bool = False) -> dict:
    from concurrent.futures import ProcessPoolExecutor
    total = 1 << num_entries(n)
    if update_check and n > 5:
        raise ValueError("exhaustive update checks are only enabled for n <= 5")
    chunk = max(1, total // (workers * 8))
    ranges = [(lo, min(lo + chunk, total)) for lo in range(0, total, chunk)]
    t0 = time.time()
    agg: dict = {}
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for part in pool.map(_exhaustive_worker,
                             [(n, lo, hi, 0, update_check, 0)
                              for lo, hi in ranges]):
            done += part["count"]
            for key, slot in part["agg"].items():
                cur = agg.setdefault(key, {"max_ratio": None, "where": None,
                                           "violations": 0, "inactive": 0})
                cur["violations"] += slot.get("violations", 0)
                cur["inactive"] += slot.get("inactive", 0)
                if slot.get("max_ratio") is not None and (
                        cur["max_ratio"] is None or slot["max_ratio"] > cur["max_ratio"]):
                    cur["max_ratio"] = slot["max_ratio"]
                    cur["where"] = slot.get("where")
    return {"mode": "exhaustive", "n": n, "signings": done,
            "checks": agg, "seconds": round(time.time() - t0, 2)}


def run_families(mc_samples: int, seed: int, update_check: bool = False) -> dict:
    results = {}
    agg: dict = {}
    # conferences (prime q = 1 mod 4)
    primes = [q for q in range(5, 160) if q % 4 == 1 and
              all(q % d for d in range(2, int(math.isqrt(q)) + 1))]
    for q in primes:
        C = paley_conference(q)
        n = q + 1
        gram = C @ C
        gram_err = float(np.abs(gram - q * np.eye(n)).max())
        # MC is cheap for n <= 62; deterministic for all
        samples = mc_samples if n <= 62 else 0
        res = analyze(C, mc_samples=samples, seed=seed + q,
                      update_check=(update_check and n <= 18))
        key = f"paley_q{q}"
        results[key] = res
        results[key]["gram_err"] = gram_err
        merge_check_agg(agg, res, {"family": key})
    # conference + forced local flips
    for q in (13, 29, 61):
        C = paley_conference(q)
        n = q + 1
        rng = np.random.default_rng(seed + q)
        for k in (1, 3, 10):
            for rep in range(3):
                A = C.copy()
                iu = np.triu_indices(n, 1)
                pick = rng.choice(len(iu[0]), size=k, replace=False)
                for p in pick:
                    i, j = int(iu[0][p]), int(iu[1][p])
                    A[i, j] = -A[i, j]
                    A[j, i] = A[i, j]
                res = analyze(A, mc_samples=mc_samples // 5, seed=seed + 7 * q + rep)
                key = f"paley_q{q}_flip{k}_{rep}"
                results[key] = res
                merge_check_agg(agg, res, {"family": key})
    # random signings and the all-ones signing
    for n in (8, 16, 32, 64, 128, 256):
        for rep in range(3):
            rng = np.random.default_rng(seed + 1000 * n + rep)
            iu = np.triu_indices(n, 1)
            vals = rng.choice([-1.0, 1.0], size=len(iu[0]))
            A = np.zeros((n, n))
            A[iu] = vals
            A = A + A.T
            res = analyze(A, mc_samples=mc_samples if n <= 64 else mc_samples // 5,
                          seed=seed + n + rep,
                          update_check=(update_check and n <= 16))
            key = f"random_n{n}_{rep}"
            results[key] = res
            merge_check_agg(agg, res, {"family": key})
    for n in (8, 32, 128):
        A = all_ones_offdiag(n)
        res = analyze(A, mc_samples=mc_samples // 5, seed=seed + n)
        key = f"ones_n{n}"
        results[key] = res
        merge_check_agg(agg, res, {"family": key})
    return {"mode": "families", "checks": agg, "results": results}


def run_hillclimb(n: int, restarts: int, steps: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(n, 1)
    best = {"ratio": -1.0, "A": None}
    for rep in range(restarts):
        vals = rng.choice([-1.0, 1.0], size=len(iu[0]))
        A = np.zeros((n, n))
        A[iu] = vals
        A = A + A.T

        def ratio_of(A):
            M = A / math.sqrt(n)
            lam, Q = np.linalg.eigh(M)
            S = (Q * np.where(lam >= 0, 1.0, -1.0)) @ Q.T
            ell = float(np.abs(lam).sum() / n)
            d = (n - 1.0) / n + 1.0 - 2.0 * ell
            if d <= 1e-12:
                return -1.0, 0.0
            R, _, _ = build_frame(S)
            lhs = float(np.linalg.norm(R - (np.eye(n) + S)))
            return lhs / (4.0 * math.sqrt(n * d)), d

        cur, _ = ratio_of(A)
        for _ in range(steps):
            improved = False
            order = rng.permutation(len(iu[0]))
            for p in order:
                i, j = iu[0][p], iu[1][p]
                A[i, j] = -A[i, j]
                A[j, i] = A[i, j]
                r, _ = ratio_of(A)
                if r > cur + 1e-12:
                    cur = r
                    improved = True
                    break
                A[i, j] = -A[i, j]
                A[j, i] = A[i, j]
            if not improved:
                break
        if cur > best["ratio"]:
            best = {"ratio": cur, "A": A.copy()}
    return {"mode": "hillclimb", "n": n, "restarts": restarts, "steps": steps,
            "best_ratio_c5": best["ratio"]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exhaustive", type=int, default=0)
    ap.add_argument("--families", action="store_true")
    ap.add_argument("--update-check", action="store_true",
                    help="also check the update bound (9) where Phi(A) is enumerable")
    ap.add_argument("--hillclimb", type=int, default=0,
                    help="order n for the exploratory ratio search")
    ap.add_argument("--restarts", type=int, default=40)
    ap.add_argument("--steps", type=int, default=40)
    ap.add_argument("--mc-samples", type=int, default=200000)
    ap.add_argument("--workers", type=int, default=80)
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    out: dict = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "kappa": KAPPA,
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runs": [],
    }
    if args.exhaustive:
        out["runs"].append(run_exhaustive(args.exhaustive, args.workers,
                                          args.mc_samples,
                                          update_check=args.update_check))
    if args.families:
        out["runs"].append(run_families(args.mc_samples, args.seed,
                                        update_check=args.update_check))
    if args.hillclimb:
        out["runs"].append(run_hillclimb(args.hillclimb, args.restarts,
                                         args.steps, args.seed))
    out["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    text = json.dumps(out, indent=1, sort_keys=True, default=str)
    if args.out:
        args.out.write_text(text + "\n")
        print(f"wrote {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
