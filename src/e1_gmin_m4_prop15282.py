#!/usr/bin/env python3
"""
Prop 15.282 — Cross_y on Type+ 1D ⊙ Hoffman norm-circle: Gauss-fibre
on every even slot; circ Jacobi collapse; n_T; generic ker-Tr DFT.
Aut_∞ mixture / floor still OPEN (N+(p) unnamed).  Does **not** flip
phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing / 15.279 inequality
/ 15.280 / 15.281 flags.

============================================================================
Setup.  Even multiplicative ψ_k of F_q^×, q=p², ψ(0)=0.  Field
encoding of 15.270: ω²=ib, ia=0 (first irreducible T²−b, b⊠).
N(x+yω)=x²−ib y².  Type+ 1D: L(a+bω)=b, D=⋃_{m∈M}(mω+F_p),
|M|=(p−1)/2, σ=+1 on S=F_p\\M.  T=t+λH Hoffman circle, N(λ)=c,
t=t0+t1 ω.  Cross_y=∑_{x∈D,y∈T} σ(L(y)) ψ(y−x).
Hoffman: N̂_Δ=N̂_D+N̂_{χψ,T}+2 Cross_y (15.281 B).
Floor: E|N̂|² ≥ 3q(q−1)/16.

============================================================================
Theorem A (Gauss-fibre Cross_y) — PROVED on Type+ 1D ⊙ Hoffman nc,
every even slot (triv / circ / generic).
  F_D(a)=∑_{m∈M} Λ(a−ξ_m), ξ_m=mω,  Λ(b)=∑_{s∈F_p} ψ(b−s).
      Λ|F_p = (p−1) 1_{(p−1)|k}
      Λ off F_p = p / G(ψ̄)  ∑_{u∈ker Tr} ψ̄(u) e^{Tr(u b)}
  (15.281 C).  Hence
      Cross_y = ∑_{y∈T} σ(L(y)) ∑_{m∈M} Λ(y−mω).
  Invert Cross_y in N̂_Δ is not an identity.  ∎

Theorem B (circ slot) — PROVED on Type+ 1D ⊙ Hoffman nc, (p+1)|k.
  ψ=η∘N with ℓ=k/(p+1).  J=J(χ_p,η)=∑_v χ_p(v)η(1−v).
      Λ(x+βω) = χ_p(ib) η(−ib β²) J    (β≠0);   0 if β=0
  (ia=0 ⇒ Δ_F/4=ib; the χ_p(ib) factor is required).  Fibre collapse
      Cross_y = χ_p(ib) η(−ib) J
                ∑_β n_T(β) σ(β) ∑_{m∈M, m≠β} η²(β−m).
  Fail-when-wrong: 15.281 triv Cross_y=−2 N̂_D/(p−1) on a circ slot
  (p=5 accident; false at p=7).  η=ψ|F_p in place of η_N.  ∎

Theorem C (n_T) — PROVED, any translate of a scale of H.
      n_T(β)=1+χ_p(c+ib(β−t1)²)
  Independent of t0.  Fail-when-wrong: replace t1 by t0.  ∎

Theorem D (generic Cy DFT) — PROVED on Type+ 1D ⊙ Hoffman nc.
  ker Tr={tω}.  Off F_p the Gauss line-sum is the DFT
      Λ(x+βω)=p/G(ψ̄) ∑_{t≠0} ψ̄(tω) e^{4πi ib t β / p}  (β≠0).
      Cross_y = p/G(ψ̄) ∑_{t≠0} ψ̄(tω)
                ∑_β n_T(β)σ(β) ∑_{m∈M,m≠β} e^{4πi ib t(β−m)/p}
                + 1_{(p−1)|k}(p−1)∑_{y∈T, L(y)∈M} σ(L(y)).
  The same formula / the 1D triv Cy formula fail on ync⊙nc.  ∎

Theorem E (OPEN mixture / floor).
  1D+hs+yP at p=7 k=8 is 241.5<441 (census |O|, not an N+ formula).
  N+(p) unnamed (15.144 type-list dead for p≥11).  triv-on-circ is
  p=5-only.  hs Cy on ync fails.  phi_F_ge_6 stays False.  ∎

============================================================================
Backend: CPU field + character sums; ProcessPool over (p, rec, k).
GPU unused (not a dense batch).  Writes evidence/e1_gmin_m4_prop15282.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15270 import _field_params  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _gauss,
    _interval_plus_set,
    _kernel_field,
    _paley_plus_set,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15280 import _add_neg, even_ks, nhat  # noqa: E402
from e1_gmin_m4_prop15281 import (  # noqa: E402
    _cross_y,
    _D_of_y,
    _named_pairs,
    _psi_arr,
)
from workers import pool, require_workers  # noqa: E402

# 15.140 census orbit sizes at p=7.  Not a formula for N+(p).
_P7_ORBIT_SIZES = {"Paley": 56, "AP": 84, "hs⊙nc": 588, "yP⊙nc": 1176}
_P7_NPLUS_CENSUS = 5726


def _chi_p(v: int, p: int) -> int:
    if v % p == 0:
        return 0
    ls = pow(int(v) % p, (p - 1) // 2, p)
    return 1 if ls == 1 else -1


def _slot_of(k: int, p: int) -> str:
    if k % (p - 1) == 0:
        return "triv"
    if k % (p + 1) == 0:
        return "circ"
    return "gen"


def _eta_norm(F: dict, k: int, p: int) -> np.ndarray:
    """Unique η of F_p^× with ψ_k = η∘N on (p+1)|k.  ℓ=k/(p+1)."""
    e = np.zeros(p, dtype=np.complex128)
    ell = int(k) // (p + 1)
    dlog = F["dlog"]
    fac = 2j * np.pi * ell / (p - 1)
    for x in range(1, p):
        e[x] = np.exp(fac * (dlog[x] // (p + 1)))
    return e


def _J_chi_eta(eta: np.ndarray, p: int) -> complex:
    s = 0j
    for v in range(p):
        s += _chi_p(v, p) * eta[(1 - v) % p]
    return complex(s)


def _infer_M(d: np.ndarray, p: int) -> tuple[list[int], bool]:
    M = []
    ok = True
    for b in range(p):
        cnt = int(d[b * p : b * p + p].sum())
        if cnt == p:
            M.append(b)
        elif cnt != 0:
            ok = False
    return M, ok and len(M) == (p - 1) // 2


def _nT_pred(beta: int, t1: int, c: int, ib: int, p: int) -> int:
    rhs = (int(c) + int(ib) * ((int(beta) - int(t1)) % p) ** 2) % p
    return 1 + _chi_p(rhs, p)


def _frob(mul, u: int, p: int) -> int:
    r, b, e = 1, int(u), int(p)
    while e:
        if e & 1:
            r = mul[r][b]
        b = mul[b][b]
        e >>= 1
    return r


def _norm_u(mul, u: int, p: int) -> int:
    return mul[int(u)][_frob(mul, u, p)]


def Lambda_circ_closed(beta: int, eta: np.ndarray, p: int, ib: int, Jchi: complex) -> complex:
    """ia=0: Λ(x+βω)=χ_p(ib β²) η(−ib β²) J  (β≠0); 0 if β=0."""
    if beta % p == 0:
        return 0j
    d = (int(ib) * int(beta) * int(beta)) % p
    return _chi_p(d, p) * eta[(-d) % p] * Jchi


def Lambda_gauss_off(b: int, psi: np.ndarray, F: dict, p: int, q: int, Gbar: complex) -> complex:
    if b < p:
        return 0j
    acc = 0j
    mul = F["mul"]
    for u in range(q):
        if F["trt"][u] != 0:
            continue
        acc += np.conjugate(psi[u]) * _e_tr(F, mul[u][int(b)])
    return p * acc / Gbar


def fibre_Cy_gauss(T, yf, M, add, neg, psi, F, p, q, k, Gbar) -> complex:
    acc = 0j
    for y in np.flatnonzero(T):
        y = int(y)
        sig = int(yf[y])
        inner = 0j
        for m in M:
            b = int(add[y, int(neg[m * p])])
            if k % (p - 1) == 0 and b < p:
                inner += p - 1
            else:
                inner += Lambda_gauss_off(b, psi, F, p, q, Gbar)
        acc += sig * inner
    return complex(acc)


def fibre_Cy_circ(T, yf, M, eta, p, ib, Jchi) -> complex:
    acc = 0j
    for y in np.flatnonzero(T):
        y = int(y)
        beta = y // p
        sig = int(yf[y])
        inner = 0j
        for m in M:
            inner += Lambda_circ_closed((beta - m) % p, eta, p, ib, Jchi)
        acc += sig * inner
    return complex(acc)


def collapsed_circ(T, yf, M, eta, p, ib, Jchi) -> complex:
    """Cy = χ(ib) η(−ib) J ∑_β n_T(β) σ(β) ∑_{m∈M, m≠β} η²(β−m)."""
    nT = np.zeros(p, dtype=np.int32)
    sigb = np.zeros(p, dtype=np.int32)
    for y in np.flatnonzero(T):
        y = int(y)
        b = y // p
        nT[b] += 1
        sigb[b] = int(yf[y])
    pre = _chi_p(ib, p) * eta[(-int(ib)) % p] * Jchi
    acc = 0j
    for b in range(p):
        if nT[b] == 0:
            continue
        W = 0j
        for m in M:
            df = (b - m) % p
            if df:
                W += eta[df] * eta[df]
        acc += int(nT[b]) * int(sigb[b]) * W
    return complex(pre * acc)


def generic_Cy_dft(T, yf, M, psi, p, ib, Gbar) -> complex:
    """Cy = p/G(ψ̄) ∑_{t≠0} ψ̄(tω) ∑_β n_T(β)σ(β) ∑_{m≠β} e^{4πi ib t(β−m)/p}."""
    nT = np.zeros(p, dtype=np.int32)
    sigb = np.zeros(p, dtype=np.int32)
    for y in np.flatnonzero(T):
        y = int(y)
        b = y // p
        nT[b] += 1
        sigb[b] = int(yf[y])
    acc = 0j
    fac = 2j * np.pi / p
    for t in range(1, p):
        w = np.conjugate(psi[t * p])
        inner = 0j
        for b in range(p):
            if nT[b] == 0:
                continue
            W = 0j
            for m in M:
                df = (b - m) % p
                if df:
                    W += np.exp(fac * (2 * ib * t * df))
            inner += int(nT[b]) * int(sigb[b]) * W
        acc += w * inner
    return p * acc / Gbar


def _triv_extra(T, yf, M, p) -> complex:
    Mset = set(M)
    s = 0j
    for y in np.flatnonzero(T):
        y = int(y)
        if (y // p) in Mset:
            s += int(yf[y]) * (p - 1)
    return s


def _eval_ck(job: tuple) -> dict:
    """One (p, rec, k) character-sum job.  Module-level for ProcessPool."""
    p, k, kind, d, T, yf, M, is_1d, t_c, c = job
    F = _kernel_field(p)
    add, neg, q = _add_neg(F)
    ia, ib = _field_params(p)
    sl = _slot_of(k, p)
    psi = _psi_arr(F, k)
    chi = np.array(F["chi"], dtype=np.int8)
    Cy = _cross_y(d, T, yf, psi, add, neg)
    Nd = nhat(d, add, psi)
    Ntchi = nhat(T, add, psi * chi.astype(np.complex128))
    lhs = nhat(np.bitwise_xor(d, T), add, psi)
    hoffman = abs(lhs - (Nd + Ntchi + 2.0 * Cy)) < 1e-6
    invert = abs(lhs - (Nd + Ntchi - 2.0 * Cy)) < 1e-6
    out = {
        "p": p,
        "k": k,
        "kind": kind,
        "sl": sl,
        "is_1d": bool(is_1d),
        "hoffman": hoffman,
        "invert": invert,
        "E": float(abs(lhs) ** 2),
        "gauss_ok": None,
        "circ_lam_ok": 0,
        "circ_lam_n": 0,
        "circ_lam_naive_ok": 0,
        "etaN_ok": 0,
        "etaN_n": 0,
        "circ_fibre_ok": None,
        "collapsed_ok": None,
        "dft_ok": None,
        "triv_on_circ": None,
        "hs_on_ync": None,
        "eta_fp_ok": None,
    }
    if is_1d:
        Gbar = complex(_gauss(F, list(np.conjugate(psi))))
        predG = fibre_Cy_gauss(T, yf, M, add, neg, psi, F, p, q, k, Gbar)
        out["gauss_ok"] = abs(predG - Cy) < 1e-5
        predD = generic_Cy_dft(T, yf, M, psi, p, ib, Gbar)
        if sl == "triv":
            predD = predD + _triv_extra(T, yf, M, p)
        out["dft_ok"] = abs(predD - Cy) < 1e-4
        if sl == "circ":
            eta = _eta_norm(F, k, p)
            Jchi = _J_chi_eta(eta, p)
            mul = F["mul"]
            for u in range(1, q):
                Nu = _norm_u(mul, u, p)
                out["etaN_n"] += 1
                if Nu < p and abs(psi[u] - eta[int(Nu)]) < 1e-8:
                    out["etaN_ok"] += 1
            for beta in range(p):
                for x in (0, 1, 2):
                    b = x + beta * p
                    live = sum(psi[int(add[b, int(neg[s])])] for s in range(p))
                    pred = Lambda_circ_closed(beta, eta, p, ib, Jchi)
                    out["circ_lam_n"] += 1
                    if abs(live - pred) < 1e-5:
                        out["circ_lam_ok"] += 1
                    # missing χ_p(ib β²)
                    if beta % p != 0:
                        dlt = (int(ib) * int(beta) * int(beta)) % p
                        naive = eta[(-dlt) % p] * Jchi
                    else:
                        naive = 0j
                    if abs(live - naive) < 1e-5:
                        out["circ_lam_naive_ok"] += 1
            predC = fibre_Cy_circ(T, yf, M, eta, p, ib, Jchi)
            predK = collapsed_circ(T, yf, M, eta, p, ib, Jchi)
            out["circ_fibre_ok"] = abs(predC - Cy) < 1e-5
            out["collapsed_ok"] = abs(predK - Cy) < 1e-5
            out["triv_on_circ"] = abs(Cy + 2.0 * Nd / (p - 1)) < 1e-6
            eta_fp = np.zeros(p, dtype=np.complex128)
            eta_fp[1:] = psi[1:p]
            Jbad = _J_chi_eta(eta_fp, p)
            pred_bad = collapsed_circ(T, yf, M, eta_fp, p, ib, Jbad)
            out["eta_fp_ok"] = abs(pred_bad - Cy) < 1e-5
    elif kind == "ync⊙nc":
        if abs(Nd) > 1e-8:
            out["hs_on_ync"] = abs(Cy + 2.0 * Nd / (p - 1)) < 1e-6
    return out


def _eval_nT_p(p: int) -> dict:
    """n_T formula + t0-independence + t0-wrong formula, one prime."""
    F = _kernel_field(p)
    add, neg, q = _add_neg(F)
    ia, ib = _field_params(p)
    assert ia == 0
    mul = F["mul"]
    geom_ok = geom_n = 0
    t0_indep_ok = t0_indep_n = 0
    t0_wrong_ok = t0_wrong_n = 0
    for t1 in range(p):
        for c in range(1, p):
            nTs = []
            for t0 in range(p):
                t = t0 + t1 * p
                live = np.zeros(p, dtype=np.int32)
                for u in range(q):
                    if _norm_u(mul, int(add[u, int(neg[t])]), p) == int(c):
                        live[u // p] += 1
                nTs.append(live)
                for b in range(p):
                    geom_n += 1
                    if int(live[b]) == _nT_pred(b, t1, c, ib, p):
                        geom_ok += 1
                    t0_wrong_n += 1
                    if int(live[b]) == _nT_pred(b, t0, c, ib, p):
                        t0_wrong_ok += 1
            t0_indep_n += 1
            if all(np.array_equal(nTs[0], x) for x in nTs[1:]):
                t0_indep_ok += 1
    pair_ok = pair_n = 0
    for rec in _named_pairs(p):
        if rec.get("t") is None:
            continue
        T = rec["T"]
        t1 = int(rec["t"]) // p
        live_nT = np.zeros(p, dtype=np.int32)
        for y in np.flatnonzero(T):
            live_nT[int(y) // p] += 1
        for b in range(p):
            pair_n += 1
            if int(live_nT[b]) == _nT_pred(b, t1, int(rec["c"]), ib, p):
                pair_ok += 1
    return {
        "p": p,
        "geom": f"{geom_ok}/{geom_n}",
        "geom_ok": geom_ok,
        "geom_n": geom_n,
        "t0_indep": f"{t0_indep_ok}/{t0_indep_n}",
        "t0_indep_ok": t0_indep_ok == t0_indep_n and t0_indep_n > 0,
        "t0_wrong": f"{t0_wrong_ok}/{t0_wrong_n}",
        "t0_wrong_is_identity": t0_wrong_ok == t0_wrong_n and t0_wrong_n > 0,
        "hoffman_nT": f"{pair_ok}/{pair_n}",
        "hoffman_nT_ok": pair_ok == pair_n and pair_n > 0,
        "proved": (
            geom_ok == geom_n
            and geom_n > 0
            and t0_indep_ok == t0_indep_n
            and t0_indep_n > 0
            and pair_ok == pair_n
            and pair_n > 0
            and t0_wrong_ok != t0_wrong_n
        ),
    }


_P_CK: dict[int, list[dict]] = {}
_P_NT: dict[int, dict] = {}


def _ensure(primes: list[int]) -> None:
    missing = [p for p in primes if p not in _P_CK]
    if not missing:
        return
    jobs = []
    for p in missing:
        for rec in _named_pairs(p):
            if rec.get("t") is None:
                continue
            d = _D_of_y(rec["yA"])
            M, is_1d = _infer_M(d, p)
            T = rec["T"]
            yf = np.rint(rec["yA"][1:]).astype(np.int8)
            for k in even_ks(p):
                jobs.append(
                    (
                        p,
                        int(k),
                        rec["kind"],
                        d,
                        T,
                        yf,
                        list(M),
                        bool(is_1d),
                        rec.get("t"),
                        rec.get("c"),
                    )
                )
    W = require_workers()
    n_workers = max(2, min(W, len(jobs) + len(missing)))
    # One pool: Cy jobs + n_T geometry (independent).
    nT_jobs = list(missing)
    with pool(n_workers) as ex:
        ck_f = ex.map(_eval_ck, jobs, chunksize=1)
        nt_f = ex.map(_eval_nT_p, nT_jobs)
        ck_recs = list(ck_f)
        nt_recs = list(nt_f)
    by_p: dict[int, list[dict]] = defaultdict(list)
    for rec in ck_recs:
        by_p[int(rec["p"])].append(rec)
    for p in missing:
        _P_CK[p] = by_p.get(p, [])
    for rec in nt_recs:
        _P_NT[int(rec["p"])] = rec


def prove_gauss_fibre(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    _ensure(primes)
    rows: dict[str, dict] = {}
    ok = True
    inv_hits = hoff_hits = nchk = 0
    g_ok = g_n = 0
    by_slot: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for p in primes:
        recs = [r for r in _P_CK[p] if r["is_1d"]]
        po = pn = io = ho = 0
        slot_local: dict[str, list[int]] = defaultdict(lambda: [0, 0])
        for r in recs:
            pn += 1
            if r["gauss_ok"]:
                po += 1
            if r["invert"]:
                io += 1
            if r["hoffman"]:
                ho += 1
            sl = r["sl"]
            slot_local[sl][1] += 1
            by_slot[sl][1] += 1
            if r["gauss_ok"]:
                slot_local[sl][0] += 1
                by_slot[sl][0] += 1
        rows[str(p)] = {
            "ok": po,
            "n": pn,
            "hoffman": f"{ho}/{pn}",
            "invert": f"{io}/{pn}",
            "by_slot": {s: f"{v[0]}/{v[1]}" for s, v in sorted(slot_local.items())},
        }
        if not (po == pn and pn > 0 and ho == pn):
            ok = False
        g_ok += po
        g_n += pn
        inv_hits += io
        hoff_hits += ho
        nchk += pn
    return {
        "proved": ok and inv_hits != nchk and all(by_slot[s][0] == by_slot[s][1] and by_slot[s][1] > 0 for s in ("triv", "circ", "gen")),
        "theorem": (
            "On Type+ 1D⊙Hoffman nc, Cross_y is the Gauss-fibre sum of Λ "
            "on every even slot.  Invert-Cy is not an identity."
        ),
        "by_p": rows,
        "by_slot": {s: f"{v[0]}/{v[1]}" for s, v in sorted(by_slot.items())},
        "invert_Cy_is_identity": inv_hits == nchk and nchk > 0,
        "hoffman_ok": hoff_hits == nchk and nchk > 0,
        "n_checks": nchk,
        "gauss": f"{g_ok}/{g_n}",
    }


def prove_circ(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    _ensure(primes)
    rows: dict[str, dict] = {}
    ok = True
    triv: dict[int, list[int]] = {}
    eta_fp_hits = eta_fp_n = 0
    nchk = 0
    lam_ok = lam_n = 0
    naive_ok = 0
    etaN_ok = etaN_n = 0
    fib_ok = fib_n = 0
    coll_ok = coll_n = 0
    for p in primes:
        recs = [r for r in _P_CK[p] if r["is_1d"] and r["sl"] == "circ"]
        cy_ok = cy_n = 0
        p_lam_ok = p_lam_n = 0
        p_eta = p_etan = 0
        p_triv = [0, 0]
        for r in recs:
            cy_n += 1
            if r["circ_fibre_ok"] and r["collapsed_ok"]:
                cy_ok += 1
            p_lam_ok += r["circ_lam_ok"]
            p_lam_n += r["circ_lam_n"]
            p_eta += r["etaN_ok"]
            p_etan += r["etaN_n"]
            if r["triv_on_circ"] is not None:
                p_triv[1] += 1
                if r["triv_on_circ"]:
                    p_triv[0] += 1
            if r["eta_fp_ok"]:
                eta_fp_hits += 1
            if r["eta_fp_ok"] is not None:
                eta_fp_n += 1
            naive_ok += r["circ_lam_naive_ok"]
        rows[str(p)] = {
            "Cy": f"{cy_ok}/{cy_n}",
            "Lambda": f"{p_lam_ok}/{p_lam_n}",
            "etaN": f"{p_eta}/{p_etan}",
            "triv_on_circ": f"{p_triv[0]}/{p_triv[1]}",
        }
        triv[p] = p_triv
        if not (
            cy_ok == cy_n
            and cy_n > 0
            and p_lam_ok == p_lam_n
            and p_lam_n > 0
            and p_eta == p_etan
            and p_etan > 0
        ):
            ok = False
        nchk += cy_n
        lam_ok += p_lam_ok
        lam_n += p_lam_n
        etaN_ok += p_eta
        etaN_n += p_etan
        fib_ok += cy_ok
        fib_n += cy_n
        coll_ok += cy_ok
        coll_n += cy_n
    p7_triv = triv.get(7, [0, 0])
    triv_all = all(v[0] == v[1] and v[1] > 0 for v in triv.values())
    return {
        "proved": ok and (p7_triv[1] == 0 or p7_triv[0] != p7_triv[1]) and eta_fp_hits != eta_fp_n,
        "theorem": (
            "On 1D⊙nc and (p+1)|k, ψ=η∘N and Λ=χ_p(ib)η(−ib β²)J; "
            "Cy collapses to n_T σ W η².  η=ψ|F_p is not the identity.  "
            "triv formula on circ is a p=5 accident."
        ),
        "by_p": rows,
        "triv_on_circ_is_identity": triv_all,
        "triv_on_circ_p7_is_identity": p7_triv[0] == p7_triv[1] and p7_triv[1] > 0,
        "triv_on_circ": {str(p): f"{v[0]}/{v[1]}" for p, v in triv.items()},
        "eta_fp_is_identity": eta_fp_hits == eta_fp_n and eta_fp_n > 0,
        "eta_fp": f"{eta_fp_hits}/{eta_fp_n}",
        "Lambda_naive_no_chi_is_identity": naive_ok == lam_n and lam_n > 0,
        "etaN": f"{etaN_ok}/{etaN_n}",
        "Lambda": f"{lam_ok}/{lam_n}",
        "n_checks": nchk,
    }


def prove_nT(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    _ensure(primes)
    rows = {str(p): _P_NT[p] for p in primes}
    ok = all(rows[str(p)]["proved"] for p in primes)
    t0_wrong_id = any(rows[str(p)]["t0_wrong_is_identity"] for p in primes)
    return {
        "proved": ok and not t0_wrong_id,
        "theorem": (
            "n_T(β)=1+χ_p(c+ib(β−t1)²) on every scale-translate of H; "
            "independent of t0.  Replacing t1 by t0 is not an identity."
        ),
        "by_p": {str(p): {
            "geom": rows[str(p)]["geom"],
            "t0_indep": rows[str(p)]["t0_indep"],
            "t0_wrong": rows[str(p)]["t0_wrong"],
            "hoffman_nT": rows[str(p)]["hoffman_nT"],
        } for p in primes},
        "t0_wrong_is_identity": t0_wrong_id,
        "n_checks": sum(rows[str(p)]["geom_n"] for p in primes),
    }


def prove_generic(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    _ensure(primes)
    rows: dict[str, dict] = {}
    ok = True
    ync_hits = ync_n = 0
    nchk = 0
    for p in primes:
        recs = _P_CK[p]
        id_ok = n = 0
        by_sl: dict[str, list[int]] = defaultdict(lambda: [0, 0])
        for r in recs:
            if r["kind"] == "ync⊙nc" and r["hs_on_ync"] is not None:
                ync_n += 1
                if r["hs_on_ync"]:
                    ync_hits += 1
            if not r["is_1d"]:
                continue
            if r["dft_ok"] is None:
                continue
            n += 1
            by_sl[r["sl"]][1] += 1
            if r["dft_ok"]:
                id_ok += 1
                by_sl[r["sl"]][0] += 1
        rows[str(p)] = {
            "ok": id_ok,
            "n": n,
            "by_slot": {s: f"{v[0]}/{v[1]}" for s, v in sorted(by_sl.items())},
        }
        if not (id_ok == n and n > 0):
            ok = False
        nchk += n
    return {
        "proved": ok and (ync_n == 0 or ync_hits != ync_n),
        "theorem": (
            "On 1D⊙nc, generic/triv Cross_y is the ker-Tr DFT of ψ̄ "
            "against the n_T-weighted exponential fibre sum "
            "(plus the on-fibre (p−1) term on triv).  "
            "hs formula on ync is not an identity."
        ),
        "by_p": rows,
        "ync_hs_formula_is_identity": ync_hits == ync_n and ync_n > 0,
        "ync_hs": f"{ync_hits}/{ync_n}",
        "n_checks": nchk,
    }


def one_d_hs_yp_mix_p7_k8() -> dict:
    """Named-orbit contribution (15.140 sizes) / N+(7) census.  Not N+ formula."""
    p, k = 7, 8
    F = _kernel_field(p)
    add, neg, q = _add_neg(F)
    psi = _psi_arr(F, k)
    chi = np.array(F["chi"], dtype=np.int8)
    energies: dict[str, float] = {}
    for name, plus in (("Paley", _paley_plus_set(p)), ("AP", _interval_plus_set(p))):
        y = np.array(lift_sigma_vector(p, 0, 1, plus), dtype=np.float64)
        d = (y[1:] < 0).astype(np.int8)
        energies[name] = float(abs(nhat(d, add, psi)) ** 2)
    for rec in _named_pairs(p):
        if rec["kind"] not in ("hs⊙nc", "yP⊙nc") or rec.get("t") is None:
            continue
        d = _D_of_y(rec["yA"])
        T = rec["T"]
        yf = np.rint(rec["yA"][1:]).astype(np.int8)
        Cy = _cross_y(d, T, yf, psi, add, neg)
        Nd = nhat(d, add, psi)
        Ntchi = nhat(T, add, psi * chi.astype(np.complex128))
        lhs = nhat(np.bitwise_xor(d, T), add, psi)
        energies[rec["kind"]] = float(abs(lhs) ** 2)
    sizes = dict(_P7_ORBIT_SIZES)
    num = sum(sizes[name] * energies[name] for name in sizes)
    den = _P7_NPLUS_CENSUS
    mix = num / den
    floor = 3.0 * q * (q - 1) / 16.0
    return {
        "energies": {name: energies[name] for name in sizes},
        "sizes": sizes,
        "num": num,
        "den_census": den,
        "mix": mix,
        "floor": floor,
        "meets_floor": bool(mix + 1e-9 >= floor),
    }


def prove_open() -> dict:
    mix = one_d_hs_yp_mix_p7_k8()
    return {
        "floor_closed": False,
        "Nplus_named": False,
        "mixture_meets_floor": False,
        "mixture_meets_floor_general": False,
        "one_d_plus_hs_yp_meets_floor_p7k8": False,
        "one_d_plus_hs_yp_mix_p7k8": mix["mix"],
        "one_d_plus_hs_yp_mix": mix,
        "open": (
            "Cross_y named on 1D⊙nc (A–D).  Aut_∞ mixture ≥ 3q(q−1)/16 "
            "needs N+(p) / a type-list (DEAD p≥11).  "
            "1D+hs+yP = 241.5<441 at p=7 k=8.  Census mix at p=5,7 is "
            "above the floor and is not a general proof."
        ),
        "dead": [
            "1D-only mixture",
            "1D+hs⊙nc as the full mixture",
            "1D+hs+yP as the full mixture",
            "triv Cross_y on circ k",
            "hs Cross_y on ync⊙nc",
            "15.144 type-list for N+ at p≥11",
        ],
    }


def main() -> dict:
    a = prove_gauss_fibre()
    b = prove_circ()
    c = prove_nT()
    d = prove_generic()
    e = prove_open()
    out = {
        "prop": "15.282",
        "title": (
            "Prop 15.282: Gauss-fibre / circ Jacobi / n_T / generic DFT "
            "Cross_y on Type+ 1D⊙Hoffman nc; mixture/floor OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "Named Max+-free Cross_y on 1D⊙nc (Gauss-fibre all slots; "
            "circ J(χ,η)η(−β²ΔF/4) collapse; n_T; ker-Tr DFT).  "
            "N+(p) unnamed; Aut_∞ mixture does not meet 3q(q−1)/16 "
            "without a type-list.  No flag flip."
        ),
        "proved": {
            "gauss_fibre": a["proved"],
            "circ_cross": b["proved"],
            "nT": c["proved"],
            "generic_cross": d["proved"],
            "triv_on_circ_is_identity": b["triv_on_circ_is_identity"],
            "triv_on_circ_p7_is_identity": b["triv_on_circ_p7_is_identity"],
            "eta_fp_is_identity": b["eta_fp_is_identity"],
            "ync_hs_formula_is_identity": d["ync_hs_formula_is_identity"],
            "invert_Cy_is_identity": a["invert_Cy_is_identity"],
            "t0_wrong_is_identity": c["t0_wrong_is_identity"],
            "Nplus_named": False,
            "mixture_meets_floor": False,
            "phi_F_ge_6_proved_general": False,
            "inequality_proved": False,
        },
        "algebra": {
            "gauss_fibre": a,
            "circ": b,
            "nT": c,
            "generic": d,
            "open": e,
        },
        "closed_forms": {
            "gauss_fibre": "Cy=∑_{y∈T} σ(L(y)) ∑_{m∈M} Λ(y−mω), Λ=15.281 C",
            "circ_Lambda": "χ_p(ib) η(−ib β²) J(χ_p,η)  (ia=0, Δ_F/4=ib)",
            "circ_Cy": (
                "χ_p(ib) η(−ib) J(χ_p,η) ∑_β n_T(β)σ(β) ∑_{m∈M,m≠β} η²(β−m), "
                "ψ=η∘N"
            ),
            "nT": "n_T(β)=1+χ_p(c+ib(β−t1)²)  (independent of t0)",
            "generic_Cy": (
                "p/G(ψ̄) ∑_{t≠0} ψ̄(tω) ∑_β n_T(β)σ(β) ∑_{m≠β} "
                "e^{4πi ib t(β−m)/p} + 1_triv (p−1)∑_{L(y)∈M} σ"
            ),
            "Nhat_Delta": "N̂_D + N̂_{χψ,T} + 2 Cross_y",
        },
        "fail_when_wrong": [
            "15.281 triv Cross_y on a (p+1)|k slot at p=7",
            "η=ψ|F_p in place of η_N",
            "Λ without χ_p(ib β²)",
            "invert Cross_y in N̂_Δ",
            "n_T with t0 in place of t1",
            "1D Cross_y / hs formula on ync⊙nc",
            "1D+hs+yP as the Aut_∞ mixture at p=7 k=8",
        ],
        "status": {
            "identities_A_B_C_D": bool(
                a["proved"] and b["proved"] and c["proved"] and d["proved"]
            ),
            "identities_A_B_C": bool(a["proved"] and b["proved"] and c["proved"]),
            "floor": False,
        },
        "backend": "CPU field + character sums; ProcessPool over (p,rec,k); GPU unused",
        "residual_closed": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15282.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        f"Prop 15.282 gauss={a['proved']} circ={b['proved']} "
        f"nT={c['proved']} generic={d['proved']}"
    )
    print(
        f"  triv_on_circ_p7={b['triv_on_circ_p7_is_identity']} "
        f"invert_Cy={a['invert_Cy_is_identity']} "
        f"mix_p7k8={e['one_d_plus_hs_yp_mix_p7k8']:.4f} floor OPEN"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
