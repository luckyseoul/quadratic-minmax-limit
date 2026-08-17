#!/usr/bin/env python3
"""
Prop 15.298 — Full-ensemble Q is the half-Gauss image of L_χ.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.297 flags.  Does **not** name the three off values in p.
Floor stays OPEN.

============================================================================
Setup.  q=p², k=|D|=p(p−1)/2, N(δ)=|D∩(D−δ)|, u=|ẑ|² on Ω,
Q(r)=E[u(ξ)u(rξ)] for r∈T=□.  L_χ(s)=E[N(d)N(s d)] for any d
with χ(d)=χ.  Half-Gauss (β≠0):
    S_χ(β)=(-1 + χ(β) G)/2,   G=G(χ,ξ)=±p on Ω,
    S_χ(0)=(q−1)/2.
μ_+=p(p−1)/4, μ_-=p(p−3)/4.  Ensemble is only ×□-invariant.

============================================================================
Theorem A — PROVED (Fourier inversion of N; Max+-free).
  |ẑ(ξ)|² = k + ∑_{δ≠0} (p(2−p)+4 N(δ)) e_ξ(δ)  (15.279 T),
  so Q(r)/16 equals
      k² + k(G_+ + G_r)
      + ∑_{s≠0} ∑_{χ=±1} L_χ(s) S_χ(1+r s)
  with G_+ = ∑_{δ≠0} E[N(δ)] e_ξ(δ) the named Gauss image of μ_±.
  L_χ(s) is constant on: s∈□, the 15.290 Paley×norm type of s;
  s∉□, the pair (Paley(s), N_{q/p}(s)∈F_p).
  Fail-when-wrong: replace L_χ(s) by μ_+ μ_- for every nonsquare s
  (the naive 15.279 W extension).  That residual is ~430 q⁰ at p=5.  ∎

Theorem B — CERTIFIED p=3,5,7 (not a general orbit proof).
  Q takes at most three values on T\\{±1}:
      Q_++  on Paley ++ (sub and ++ N1 merge),
      Q_N*  on every N* class (Paley types merge),
      Q_{--1} on Paley -- ∩ N1.
  At p=5 the mixed Paley N1 class equals Q_++ and --N1 is empty.
  At p=3 the only off type is --N1 and Q=0.
  Fail-when-wrong: claim Q(++ sub)≠Q(++ N1) at p=7, or that the
  three N* Paley classes split.  ∎

Theorem C — OPEN.  The three off values are not yet a Gauss/Jacobi
  formula in p (denominators 13, 409 = N_+/(2p), N_+ unnamed).
  Conference+row-sum 4-linear system on Aut_∞ 4-set types has
  nullity 1 (p=5) and 5 (p=7); it does not name the 4-point.
  ⟨δ,ψ⟩≤2 / phi_F stay OPEN.  ∎

============================================================================
Backend: CPU, vectorized N on the p=5,7 caches (one pass).  p=3
optional.  Serial outer: two primes, inherently a single cache fold.
Writes evidence/e1_gmin_m4_prop15298.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import (  # noqa: E402
    field_norm,
    live_Q_on_T,
    norm_type,
    omega_set,
    paley_pair,
    type_key,
    zhat_omega,
)

YPATH = {3: "/tmp/maxplus_p3.npy", 5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def mu_plus(p: int) -> Fraction:
    return Fraction(p * (p - 1), 4)


def mu_minus(p: int) -> Fraction:
    return Fraction(p * (p - 3), 4)


def ns_key(F, s: int) -> tuple:
    if F["chi"][s] == 1:
        return type_key(F, s)
    return ("ns", paley_pair(F, s), int(field_norm(F, s)))


def coarse_Q_class(F, r: int) -> str:
    if r in (1, F["neg1"]):
        return "pm1"
    nt = norm_type(F, r)
    if nt == "N*":
        return "N*"
    a, b = paley_pair(F, r)
    if a == 1 and b == 1:
        return "++"
    if a == -1 and b == -1 and nt == "N1":
        return "--N1"
    return f"mixed:{a}:{b}:{nt}"


def load_N(p: int, F) -> np.ndarray:
    Y = np.load(YPATH[p])
    Ds = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64)) < 0
    add = np.array(F["add"], dtype=np.int32)
    return np.stack([(Ds & Ds[:, add[:, d]]).sum(axis=1) for d in range(F["q"])], axis=1)


def live_Q_any(p: int) -> dict[int, float]:
    if p in (5, 7):
        return live_Q_on_T(p)
    F = _kernel_field(p)
    Y = np.load(YPATH[p])
    Z = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + F["q"]].astype(np.float64))
    Omega = omega_set(F)
    U = np.stack([np.abs(zhat_omega(Z[i], F, Omega)) ** 2 for i in range(Z.shape[0])])
    xi0 = Omega[0]
    om_idx = {eta: i for i, eta in enumerate(Omega)}
    u0 = U[:, om_idx[xi0]]
    out = {}
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * U[:, om_idx[eta]]))
    return out


def half_gauss(F, beta: int, chi_d: int, Gchi, q: int):
    if beta == 0:
        return (q - 1) / 2
    return (-1 + chi_d * F["chi"][beta] * Gchi) / 2


def build_L(N: np.ndarray, F) -> dict:
    q = F["q"]
    ENN = (N.astype(np.float64).T @ N) / N.shape[0]
    L = {}
    for d in range(1, q):
        for s in range(1, q):
            key = (int(F["chi"][d]), ns_key(F, s))
            if key not in L:
                L[key] = float(ENN[d, F["mul"][s][d]])
    return L


def predict_Q(p: int, F, L: dict, EN: np.ndarray, naive_ns: bool = False) -> dict[int, float]:
    q = F["q"]
    k = p * (p - 1) / 2.0
    Omega = omega_set(F)
    xi = Omega[0]
    Gchi = sum(F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q))
    mu_p = float(mu_plus(p))
    mu_m = float(mu_minus(p))
    mu_cross = mu_p * mu_m

    def G_of(scale):
        acc = 0j
        for d in range(1, q):
            acc += EN[d] * _e_tr(F, F["mul"][F["mul"][xi][scale]][d] if scale != 1 else F["mul"][xi][d])
        return acc

    Gp = G_of(1)
    out = {}
    for r in F["squares"]:
        Gr = G_of(r)
        acc = k * k + k * (Gp + Gr)
        for s in range(1, q):
            beta = F["add"][1][F["mul"][r][s]]
            for cd in (1, -1):
                if naive_ns and F["chi"][s] == -1:
                    Lval = mu_cross
                else:
                    Lval = L[(cd, ns_key(F, s))]
                acc += Lval * half_gauss(F, beta, cd, Gchi, q)
        out[int(r)] = float((16 * acc).real)
    return out


def prove_gauss_identity(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    ok = True
    naive_must_fail = False
    blocks = []
    for p in primes:
        F = _kernel_field(p)
        q2 = float(F["q"] ** 2)
        N = load_N(p, F)
        EN = N.mean(axis=0)
        L = build_L(N, F)
        Q = live_Q_any(p)
        pred = predict_Q(p, F, L, EN, naive_ns=False)
        naive = predict_Q(p, F, L, EN, naive_ns=True)
        max_err = max(abs(pred[r] - Q[r]) for r in Q)
        naive_err = max(abs(naive[r] - Q[r]) for r in Q)
        if max_err / q2 > 1e-9:
            ok = False
        if naive_err / q2 > 1e-3:
            naive_must_fail = True
        blocks.append(
            {
                "p": p,
                "Nplus": int(N.shape[0]),
                "n_L_keys": len(L),
                "max_err_over_q2": max_err / q2,
                "naive_ns_err_over_q2": naive_err / q2,
            }
        )
    if not naive_must_fail:
        ok = False
    return {
        "proved": ok and naive_must_fail,
        "blocks": blocks,
        "theorem": (
            "Q(r)/16 = k² + k(G_++G_r) + ∑_s ∑_χ L_χ(s) S_χ(1+rs). "
            "μ_+μ_- on all nonsquare s fails."
        ),
    }


def prove_Q_coarsening(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [3, 5, 7]
    ok = True
    merge_seen = {"pp_merge": False, "Nstar_merge": False}
    blocks = []
    for p in primes:
        if p not in YPATH or not Path(YPATH[p]).exists():
            continue
        F = _kernel_field(p)
        Q = live_Q_any(p)
        q2 = float(F["q"] ** 2)
        by_c: dict[str, list[float]] = defaultdict(list)
        by_t: dict[tuple, list[float]] = defaultdict(list)
        for r, Qr in Q.items():
            by_c[coarse_Q_class(F, r)].append(Qr / q2)
            by_t[type_key(F, r)].append(Qr / q2)
        spreads = {c: float(max(vs) - min(vs)) for c, vs in by_c.items()}
        if max(spreads.values()) > 1e-6:
            ok = False
        means = {c: str(Fraction(float(np.mean(vs))).limit_denominator(2000)) for c, vs in by_c.items()}
        if p == 7:
            tpp_sub = by_t.get((1, 1, "sub"), [])
            tpp_n1 = by_t.get((1, 1, "N1"), [])
            if tpp_sub and tpp_n1 and abs(tpp_sub[0] - tpp_n1[0]) / q2 < 1e-6:
                merge_seen["pp_merge"] = True
            nstar = [by_t[t][0] for t in by_t if len(t) == 3 and t[2] == "N*"]
            if len(nstar) >= 2 and max(nstar) - min(nstar) < 1e-6 * q2:
                merge_seen["Nstar_merge"] = True
        blocks.append({"p": p, "spreads": spreads, "means": means, "n_coarse": len(by_c)})
    if not (merge_seen["pp_merge"] and merge_seen["Nstar_merge"]):
        ok = False
    return {
        "proved": ok,
        "merges": merge_seen,
        "blocks": blocks,
        "theorem": (
            "At p=3,5,7, Q is constant on {pm1, ++, N*, --N1, mixed}. "
            "++ sub merges with ++ N1; all N* Paley classes merge."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "note": (
            "Gauss image names Q in terms of L_χ, not in p. "
            "Three off values still unnamed. Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.298  Q = half-Gauss image of L_χ; 3 off-classes certified", flush=True)
    A = prove_gauss_identity()
    print(f"  A Gauss identity: {A['proved']} {A['blocks']}", flush=True)
    B = prove_Q_coarsening()
    print(f"  B Q coarsening: {B['proved']} {B['merges']}", flush=True)
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.298",
        "title": "Full-ensemble Q is the half-Gauss image of L_χ",
        "proved": {
            "gauss_identity": A["proved"],
            "Q_coarsening_certified": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15298.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
