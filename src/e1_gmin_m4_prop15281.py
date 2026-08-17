#!/usr/bin/env python3
"""
Prop 15.281 — Cross_ψ on Hoffman pairs as a Jacobi / Gauss / Hoffman
sum; |N̂_Δ|² on 1D⊙nc at triv slots.  Aut_∞ mixture / floor still OPEN
(N+(p) unnamed).  Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur /
Gsum / pairing / 15.279 inequality / 15.280 flags.

============================================================================
Setup.  Even multiplicative ψ_k of F_q^×, q=p², ψ(0)=0.
N̂_S(ψ)=∑_{x,y∈S} ψ(y−x).  Cross_ψ(A,B)=∑_{x∈A,y∈B} ψ(y−x).
Cross_y(D,T)=∑_{x∈D,y∈T} y_A(y) ψ(y−x).
T=t+λH a norm circle, H=ker N, |H|=p+1, N(λ)=c.
Hoffman: χ(b−a)=y_A(a)y_A(b) on T; I=D∩T, |I|=|T|/2 (15.280 C).
1_H(x)=(1/(p−1))∑_{ℓ=0}^{p−2} (η_ℓ∘N)(x) on F_q^×.
J(φ,ψ)=∑_v φ(v)ψ(1−v).  F_D^χ(a)=∑_{x∈D} χ(a−x).
η_ℓ∘N = ψ_{ℓ(p+1)}.  Floor: E|N̂|² ≥ 3q(q−1)/16.

============================================================================
Theorem A (Jacobi Cross on a circle) — PROVED, any D, circle T.
  S(α)=∑_{z∈H} ψ(z−α).  For α≠0,
      S(α)=ψ(α)/(p−1) ∑_ℓ (η_ℓ∘N)(α) J(η_ℓ∘N, ψ).
  Cross(D,T)=ψ(λ) ∑_{x∈D} S((x−t)/λ), field division.
  Hence
      Cross_ψ(D, t+λH)
        = 1_{t∈D}(p+1) 1_{(p+1)|k} ψ(λ)
          + (1/(p−1)) ∑_ℓ J(η_ℓ∘N, ψ) η_ℓ(c^{-1}) F_D^{ψ·(η_ℓ∘N)}(t).
  η_ℓ(c^{-1})=ψ_ℓ(c^{-1})  (dlog_q/(q−1), not /(p−1)).
  Fail-when-wrong: invert the Cross sign in 15.280 A; use dlog/(p−1)
  for η_ℓ(c^{-1}).  ∎

Theorem B (Hoffman N̂_Δ / N̂_I / Cross(I,J)) — PROVED on Hoffman pairs.
  χψ is even and shares the (p+1)|k slot with ψ, so |N̂_{χψ,T}|²
  is 15.280 B.  S_1=∑_{z∈H} ψ(z−1),  J_H=∑_{w∈H\\{1}} χ(w−1)ψ(w),
  μ=∑_{x∈T} y_A(x) ψ(x−t).
      Cross(I,J)=(N̂_T − N̂_{χψ,T})/4
      μ = y_A(a0) ψ(a0−t) [1 + χ_p(c) J_H]   (any a0∈T)
      N̂_I=(N̂_T + N̂_{χψ,T} − 2 S_1 μ)/4
      N̂_{DΔT} = N̂_D + N̂_{χψ,T} + 2 Cross_y(D,T).
  The last is the E_KJ expansion with Cross evaluated: N̂_K+N̂_J+2 Cross(K,J)
  collapses by Hoffman.  Invert Cross_y is not an identity.
  Half-formula N̂_I=(N̂_T+N̂_{χψ,T})/4 holds iff S_1 μ=0 (not automatic).
  Fail-when-wrong: invert Cross_y; apply B on a non-Hoffman T.  ∎

Theorem C (1D Gauss F_D and triv |N̂_Δ|²) — PROVED on Type+ 1D ⊙ nc.
  D=⋃_{m∈M} (ξ_m+F_p) with |M|=(p−1)/2, 0∉M (L=b fibres).
  F_D(a)=∑_{m∈M} Λ(a−ξ_m),  Λ(b)=∑_{s∈F_p} ψ(b−s).
      Λ|F_p = (p−1) 1_{(p−1)|k}
      Λ off F_p = p / G(ψ̄)  ∑_{u∈ker Tr} ψ̄(u) e^{Tr(u b)}.
  On (p−1)|k (and 1D ⊙ Hoffman circle):
      Cross_y = −2 N̂_D /(p−1)
      |N̂_Δ|² = p²(p+1)²(p−5)²/16
  (N̂_{χψ,T}=0 on these k).  The same Cross_y formula fails on ync⊙nc
  and on size-12.  p=5 accident Cross=N̂_D/|D| fails on yP⊙nc at p=7,11
  and on p=7 hs off triv.  Fail-when-wrong: apply the 1D triv formula
  to ync; apply Cross=N̂_D/|D| to yP⊙nc; invert (p+1)|k on N̂_T.  ∎

Theorem D (OPEN mixture / floor).
  E|N̂|²=∑_O (|O|/N+(p)) E_O.  1D E_O and 1D⊙nc triv E_O are named.
  General-k Cross_y on 1D is the Jacobi+Gauss of C; on ync it is the
  same Jacobi of F_{D_A}+F_{T_0}−2 F_{I_0} (algebra).  Weights |O|/N+
  stay unnamed (N+=6,130,5726 at p=3,5,7; 409 prime; 15.144 type-list
  dead for p≥11).  1D+hs⊙nc still 224<441 at p=7 k=8.  phi_F_ge_6
  stays False.  ∎

============================================================================
Backend: CPU field + character sums.  GPU unused.
Writes evidence/e1_gmin_m4_prop15281.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15139 import (  # noqa: E402
    _field_ops,
    affine_halfspace,
    double_switch,
    norm_circle_z_on,
    seidel_switch,
)
from e1_gmin_m4_prop15279 import _e_tr, _gauss, _kernel_field, _psi_vec  # noqa: E402
from e1_gmin_m4_prop15280 import (  # noqa: E402
    _add_neg,
    cross,
    even_ks,
    nhat,
    unit_circle,
)
from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power  # noqa: E402


def _psi_arr(F: dict, k: int) -> np.ndarray:
    q = F["q"]
    k = int(k) % (q - 1)
    return np.array(_psi_vec(F, k), dtype=np.complex128)


def _D_of_y(y: np.ndarray) -> np.ndarray:
    return (np.rint(y[1:]) < 0).astype(np.int8)


def _cross_y(d, t, yf, psi, add, neg) -> complex:
    ia = np.flatnonzero(d)
    ib = np.flatnonzero(t)
    if ia.size == 0 or ib.size == 0:
        return 0j
    diffs = add[np.ix_(ib, neg[ia])]
    w = yf[ib].astype(np.complex128)
    return complex((w[:, None] * psi[diffs]).sum())


def _F_at(d, a, psi, add, neg) -> complex:
    idx = np.flatnonzero(d)
    if idx.size == 0:
        return 0j
    return complex(psi[add[int(a), neg[idx]]].sum())


def _jacobi(phi, psi, add, neg, q) -> complex:
    one = 1
    s = 0j
    for v in range(q):
        s += phi[v] * psi[int(add[one, int(neg[v])])]
    return complex(s)


def _pred_cross_jacobi(d, t_c, c, lam, psi, F, add, neg, p, q, k) -> complex:
    slot_term = 0j
    if int(d[int(t_c)]) == 1 and k % (p + 1) == 0:
        slot_term = (p + 1) * psi[int(lam)]
    inv_c = F["inv"][int(c)]
    acc = 0j
    for ell in range(p - 1):
        phi = _psi_arr(F, ell * (p + 1))
        J = _jacobi(phi, psi, add, neg, q)
        eta = _psi_arr(F, ell)[int(inv_c)] if ell else (1.0 + 0j)
        acc += J * eta * _F_at(d, int(t_c), psi * phi, add, neg)
    return slot_term + acc / (p - 1)


def _pred_cross_jacobi_bad_eta(d, t_c, c, lam, psi, F, add, neg, p, q, k) -> complex:
    """Fail-when-wrong: η_ℓ(c^{-1}) via dlog_q/(p−1)."""
    slot_term = 0j
    if int(d[int(t_c)]) == 1 and k % (p + 1) == 0:
        slot_term = (p + 1) * psi[int(lam)]
    inv_c = F["inv"][int(c)]
    dlog = F["dlog"]
    acc = 0j
    for ell in range(p - 1):
        phi = _psi_arr(F, ell * (p + 1))
        J = _jacobi(phi, psi, add, neg, q)
        if ell == 0:
            eta = 1.0 + 0j
        else:
            eta = np.exp(2j * np.pi * ell * dlog[int(inv_c)] / (p - 1))
        acc += J * eta * _F_at(d, int(t_c), psi * phi, add, neg)
    return slot_term + acc / (p - 1)


def _pick_lambda(F, p, c: int) -> int:
    q = p * p
    mul = F["mul"]
    # N(u)=u * u^p ; u^p via repeated square
    def frob(u: int) -> int:
        r, b, e = 1, u, p
        while e:
            if e & 1:
                r = mul[r][b]
            b = mul[b][b]
            e >>= 1
        return r

    for u in range(1, q):
        if mul[u][frob(u)] == int(c):
            return u
    raise RuntimeError(c)


def _S1(H, psi, add, neg) -> complex:
    s = 0j
    for z in H:
        s += psi[int(add[int(z), int(neg[1])])]
    return complex(s)


def _JH(H, psi, chi, add, neg) -> complex:
    s = 0j
    for w in H:
        if int(w) == 1:
            continue
        wm1 = int(add[int(w), int(neg[1])])
        s += int(chi[wm1]) * psi[int(w)]
    return complex(s)


def _named_pairs(p: int) -> list[dict]:
    """ystar (hs⊙nc) and, at p≥7, one yP⊙nc + one ync⊙nc if they exist."""
    q = p * p
    hs = halfspace_boolean_vector(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, addf, negf, frob, qq = _field_ops(p)
    out: list[dict] = []
    rec = construct_ystar(p)
    if rec.get("found"):
        T = np.zeros(q, dtype=np.int8)
        for s in rec["S"]:
            T[s - 1] = 1
        out.append(
            {
                "kind": "hs⊙nc",
                "yA": hs,
                "y": rec["y"],
                "T": T,
                "t": rec["t"],
                "c": rec["c"],
                "hoffman": True,
            }
        )
        if p >= 7:
            Cnc = seidel_switch(C, rec["y"])
            ncs = norm_circle_z_on(Cnc, p, mul, addf, negf, frob, q)
            for nc in ncs[:1]:
                y = double_switch(rec["y"], nc["z"], C, p)
                if y is None:
                    continue
                T2 = np.zeros(q, dtype=np.int8)
                for s in nc["S"]:
                    T2[s - 1] = 1
                out.append(
                    {
                        "kind": "ync⊙nc",
                        "yA": rec["y"],
                        "y": y,
                        "T": T2,
                        "t": nc["t"],
                        "c": nc["c"],
                        "hoffman": True,
                    }
                )
    if p >= 7:
        qr = {(x * x) % p for x in range(1, p)}
        yP = affine_halfspace(p, (0, 1), qr | {0})
        Cp = seidel_switch(C, yP)
        for nc in norm_circle_z_on(Cp, p, mul, addf, negf, frob, q)[:1]:
            y = double_switch(yP, nc["z"], C, p)
            if y is None:
                continue
            T = np.zeros(q, dtype=np.int8)
            for s in nc["S"]:
                T[s - 1] = 1
            out.append(
                {
                    "kind": "yP⊙nc",
                    "yA": yP,
                    "y": y,
                    "T": T,
                    "t": nc["t"],
                    "c": nc["c"],
                    "hoffman": True,
                }
            )
    return out


def prove_jacobi_cross(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = {}
    ok = True
    bad_eta_hits = 0
    nchk = 0
    for p in primes:
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        pairs = _named_pairs(p)
        ks = even_ks(p)
        id_ok = bad_ok = n = 0
        for rec in pairs:
            if rec.get("t") is None:
                continue
            d = _D_of_y(rec["yA"])
            T = rec["T"]
            lam = _pick_lambda(F, p, int(rec["c"]))
            for k in ks:
                psi = _psi_arr(F, k)
                live = cross(d, T, psi, add, neg)
                pred = _pred_cross_jacobi(
                    d, rec["t"], rec["c"], lam, psi, F, add, neg, p, q, k
                )
                pred_bad = _pred_cross_jacobi_bad_eta(
                    d, rec["t"], rec["c"], lam, psi, F, add, neg, p, q, k
                )
                n += 1
                if abs(live - pred) < 1e-5:
                    id_ok += 1
                if abs(live - pred_bad) < 1e-5:
                    bad_ok += 1
        rows[str(p)] = {"ok": id_ok, "n": n, "bad_eta_ok": bad_ok}
        if id_ok != n or n == 0:
            ok = False
        bad_eta_hits += bad_ok
        nchk += n
    return {
        "proved": ok,
        "theorem": (
            "Cross(D, t+λH) is the Jacobi combination of F_D^{ψ·(η_ℓ∘N)}(t). "
            "η via dlog/(p−1) is not the identity."
        ),
        "by_p": rows,
        "bad_eta_is_identity": bad_eta_hits == nchk and nchk > 0,
        "n_checks": nchk,
    }


def prove_hoffman_delta(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = {}
    ok = True
    inv_hits = 0
    nchk = 0
    half_ok = half_n = 0
    for p in primes:
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        chi = np.array(F["chi"], dtype=np.int8)
        H = unit_circle(F)
        ks = even_ks(p)
        d_ok = ij_ok = mu_ok = ni_ok = n = 0
        inv_ok = 0
        for rec in _named_pairs(p):
            d = _D_of_y(rec["yA"])
            T = rec["T"]
            yf = np.rint(rec["yA"][1:]).astype(np.int8)
            I = d & T
            J = T & (1 - d)
            t_c = rec.get("t")
            Tidx = np.flatnonzero(T)
            a0 = int(Tidx[0])
            c = rec.get("c")
            for k in ks:
                psi = _psi_arr(F, k)
                psi_chi = psi * chi.astype(np.complex128)
                Nd = nhat(d, add, psi)
                Nt = nhat(T, add, psi)
                Ntchi = nhat(T, add, psi_chi)
                Ni = nhat(I, add, psi)
                Xij = cross(I, J, psi, add, neg)
                Cy = _cross_y(d, T, yf, psi, add, neg)
                lhs = nhat(np.bitwise_xor(d, T), add, psi)
                n += 1
                if abs(Xij - 0.25 * (Nt - Ntchi)) < 1e-6:
                    ij_ok += 1
                if abs(lhs - (Nd + Ntchi + 2.0 * Cy)) < 1e-6:
                    d_ok += 1
                if abs(lhs - (Nd + Ntchi - 2.0 * Cy)) < 1e-6:
                    inv_ok += 1
                if t_c is not None and c is not None:
                    s1 = _S1(H, psi, add, neg)
                    mu = 0j
                    ib = np.flatnonzero(T)
                    mu = complex(
                        (yf[ib].astype(np.complex128) * psi[add[ib, int(neg[int(t_c)])]]).sum()
                    )
                    cmod = int(c) % p
                    ls = pow(cmod, (p - 1) // 2, p) if cmod else 0
                    chi_p_c = 0 if ls == 0 else (1 if ls == 1 else -1)
                    pred_mu = complex(
                        int(yf[a0])
                        * psi[int(add[a0, int(neg[int(t_c)])])]
                        * (1.0 + chi_p_c * _JH(H, psi, chi, add, neg))
                    )
                    if abs(pred_mu - mu) < 1e-6:
                        mu_ok += 1
                    if abs(Ni - 0.25 * (Nt + Ntchi - 2.0 * s1 * mu)) < 1e-6:
                        ni_ok += 1
                    half_n += 1
                    if abs(Ni - 0.25 * (Nt + Ntchi)) < 1e-6:
                        half_ok += 1
        rows[str(p)] = {
            "delta": f"{d_ok}/{n}",
            "Cross_IJ": f"{ij_ok}/{n}",
            "mu": f"{mu_ok}/{n}",
            "Nhat_I": f"{ni_ok}/{n}",
            "invert_Cy": f"{inv_ok}/{n}",
        }
        if not (d_ok == n and ij_ok == n and mu_ok == n and ni_ok == n and n > 0):
            ok = False
        inv_hits += inv_ok
        nchk += n
    # non-Hoffman leak of B at p=5
    p = 5
    F = _kernel_field(p)
    add, neg, q = _add_neg(F)
    chi = np.array(F["chi"], dtype=np.int8)
    mul, addf, negf, frob, qq = _field_ops(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    hs = halfspace_boolean_vector(p)
    Cp = seidel_switch(C, hs)
    d = _D_of_y(hs)
    yf = np.rint(hs[1:]).astype(np.int8)
    leaked = False
    n_non = 0
    psi = _psi_arr(F, 2)
    psi_chi = psi * chi.astype(np.complex128)
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if mul(addf(u, negf(t)), frob(addf(u, negf(t)))) == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if all(Cp[S[i], S[j]] > 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
                continue
            n_non += 1
            T = np.zeros(q, dtype=np.int8)
            T[pts] = 1
            I = d & T
            J = T & (1 - d)
            Xij = cross(I, J, psi, add, neg)
            Nt = nhat(T, add, psi)
            Ntchi = nhat(T, add, psi_chi)
            Cy = _cross_y(d, T, yf, psi, add, neg)
            lhs = nhat(np.bitwise_xor(d, T), add, psi)
            Nd = nhat(d, add, psi)
            if abs(Xij - 0.25 * (Nt - Ntchi)) > 1e-6 and abs(
                lhs - (Nd + Ntchi + 2.0 * Cy)
            ) > 1e-6:
                leaked = True
                break
        if leaked:
            break
    return {
        "proved": ok and leaked,
        "theorem": (
            "Hoffman: N̂_Δ=N̂_D+N̂_{χψ,T}+2 Cross_y, Cross(I,J)=(N̂T−N̂χT)/4, "
            "N̂_I via S_1 μ, μ via J_H.  Invert-Cy and non-Hoffman fail."
        ),
        "by_p": rows,
        "invert_Cy_is_identity": inv_hits == nchk and nchk > 0,
        "half_formula_always": half_ok == half_n and half_n > 0,
        "half_ok": half_ok,
        "half_n": half_n,
        "nonhoffman_breaks_B_p5": leaked,
        "n_nonhoffman_scanned": n_non,
        "n_checks": nchk,
    }


def prove_triv_1d(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = {}
    ok = True
    yp_ratio_hits = 0
    yp_n = 0
    ync_triv_hits = 0
    ync_triv_n = 0
    for p in primes:
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        pred_E = p * p * (p + 1) ** 2 * (p - 5) ** 2 / 16.0
        cy_ok = e_ok = n_1d = 0
        ratio_ok = ratio_n = 0
        g_ok = g_n = 0
        # Gauss Λ off F_p on a few k
        for k in even_ks(p)[:4]:
            psi = _psi_arr(F, k)
            Gbar = complex(_gauss(F, list(np.conjugate(psi))))
            for b in (p, p + 1, p + 2):
                live = sum(psi[int(add[b, int(neg[s])])] for s in range(p))
                acc = 0j
                for u in range(q):
                    if F["trt"][u] != 0:
                        continue
                    acc += np.conjugate(psi[u]) * _e_tr(F, F["mul"][u][b])
                pred = p * acc / Gbar
                g_n += 1
                if abs(pred - live) < 1e-5:
                    g_ok += 1
        for rec in _named_pairs(p):
            d = _D_of_y(rec["yA"])
            T = rec["T"]
            yf = np.rint(rec["yA"][1:]).astype(np.int8)
            den = float(d.sum()) or 1.0
            oned = rec["kind"] in ("hs⊙nc", "yP⊙nc")
            for k in even_ks(p):
                psi = _psi_arr(F, k)
                Nd = nhat(d, add, psi)
                Cy = _cross_y(d, T, yf, psi, add, neg)
                Xdt = cross(d, T, psi, add, neg)
                lhs = nhat(np.bitwise_xor(d, T), add, psi)
                if rec["kind"] == "yP⊙nc":
                    yp_n += 1
                    if abs(Xdt - Nd / den) < 1e-6:
                        yp_ratio_hits += 1
                if rec["kind"] == "hs⊙nc":
                    ratio_n += 1
                    if abs(Xdt - Nd / den) < 1e-6:
                        ratio_ok += 1
                if k % (p - 1) != 0 or abs(Nd) < 1e-8:
                    continue
                want_cy = -2.0 / (p - 1) * Nd
                if oned:
                    n_1d += 1
                    if abs(Cy - want_cy) < 1e-6:
                        cy_ok += 1
                    if abs(abs(lhs) ** 2 - pred_E) < 1e-3:
                        e_ok += 1
                else:
                    ync_triv_n += 1
                    if abs(Cy - want_cy) < 1e-6:
                        ync_triv_hits += 1
        rows[str(p)] = {
            "Cy_triv": f"{cy_ok}/{n_1d}",
            "energy_triv": f"{e_ok}/{n_1d}",
            "pred_E": pred_E,
            "hs_ratio_Nd_over_D": f"{ratio_ok}/{ratio_n}",
            "Lambda_gauss": f"{g_ok}/{g_n}",
        }
        if not (cy_ok == n_1d and e_ok == n_1d and n_1d > 0 and g_ok == g_n):
            ok = False
    return {
        "proved": ok
        and (yp_n == 0 or yp_ratio_hits != yp_n)
        and (ync_triv_n == 0 or ync_triv_hits != ync_triv_n),
        "theorem": (
            "1D F_D is a Gauss line-sum; triv Cross_y=−2 N̂_D/(p−1); "
            "|N̂_Δ|²=p²(p+1)²(p−5)²/16.  Formula fails on ync; "
            "Cross=N̂_D/|D| fails on yP."
        ),
        "by_p": rows,
        "yp_ratio_is_identity": yp_ratio_hits == yp_n and yp_n > 0,
        "yp_ratio": f"{yp_ratio_hits}/{yp_n}",
        "ync_triv_formula_is_identity": ync_triv_hits == ync_triv_n and ync_triv_n > 0,
        "ync_triv": f"{ync_triv_hits}/{ync_triv_n}",
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "Nplus_named": False,
        "mixture_meets_floor": False,
        "one_d_plus_hs_meets_floor_p7k8": False,
        "open": (
            "Cross is named (A–C).  Aut_∞ mixture still needs N+(p) and "
            "general-p orbit masses.  1D+hs⊙nc ≈224<441 at p=7 k=8."
        ),
        "dead": [
            "linear Cross in (N̂_D,N̂_T,N̂_I) for all types",
            "1D-only mixture",
            "1D+hs⊙nc as the full mixture",
            "Cross=N̂_D/|D| for general 1D types",
            "N̂_I=(N̂_T+N̂_{χψ})/4 on every Hoffman type",
            "15.144 type-list for N+ at p≥11",
        ],
    }


def main() -> dict:
    a = prove_jacobi_cross()
    b = prove_hoffman_delta()
    c = prove_triv_1d()
    d = prove_open()
    out = {
        "prop": "15.281",
        "title": (
            "Prop 15.281: Cross_ψ Jacobi/Hoffman/Gauss on Hoffman pairs; "
            "1D triv |N̂_Δ|²; mixture/floor OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "Named Max+-free Cross (A) and Hoffman N̂_Δ (B) and 1D triv energy (C). "
            "N+(p) unnamed; Aut_∞ mixture does not meet 3q(q−1)/16 without it. "
            "No flag flip."
        ),
        "proved": {
            "jacobi_cross": a["proved"],
            "hoffman_delta": b["proved"],
            "triv_1d_energy": c["proved"],
            "bad_eta_is_identity": a["bad_eta_is_identity"],
            "invert_Cy_is_identity": b["invert_Cy_is_identity"],
            "yp_ratio_is_identity": c["yp_ratio_is_identity"],
            "ync_triv_formula_is_identity": c["ync_triv_formula_is_identity"],
            "Nplus_named": False,
            "mixture_meets_floor": False,
            "phi_F_ge_6_proved_general": False,
            "inequality_proved": False,
        },
        "algebra": {"jacobi": a, "hoffman": b, "triv_1d": c, "open": d},
        "closed_forms": {
            "Cross": (
                "1_{t∈D}(p+1)1_{(p+1)|k} ψ(λ) + "
                "(1/(p−1))∑_ℓ J(η_ℓ∘N,ψ) η_ℓ(c^{-1}) F_D^{ψ(η_ℓ∘N)}(t)"
            ),
            "Nhat_Delta": "N̂_D + N̂_{χψ,T} + 2 Cross_y",
            "Cross_IJ": "(N̂_T − N̂_{χψ,T})/4",
            "Nhat_I": "(N̂_T + N̂_{χψ,T} − 2 S_1 μ)/4",
            "triv_1d_energy": "p²(p+1)²(p−5)²/16",
        },
        "fail_when_wrong": [
            "η_ℓ(c^{-1}) via dlog_q/(p−1)",
            "invert Cross_y in N̂_Δ",
            "Hoffman B on non-Hoffman T",
            "Cross=N̂_D/|D| on yP⊙nc",
            "1D triv Cross_y formula on ync⊙nc",
            "invert (p+1)|k vanishing of N̂_T (15.280 B)",
        ],
        "status": {
            "identities_A_B_C": bool(a["proved"] and b["proved"] and c["proved"]),
            "floor": False,
        },
        "backend": "CPU field + character sums; GPU unused",
        "residual_closed": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15281.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        f"Prop 15.281 jacobi={a['proved']} hoffman={b['proved']} triv1d={c['proved']}"
    )
    print(
        f"  invert_Cy={b['invert_Cy_is_identity']} "
        f"yp_ratio={c['yp_ratio_is_identity']} floor OPEN"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
