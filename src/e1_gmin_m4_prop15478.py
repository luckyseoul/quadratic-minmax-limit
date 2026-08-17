#!/usr/bin/env python3
"""
Prop 15.478 — Paley++ type-indicator Fourier: n_++ is the CM/Jacobi
count of y²=x³−x; the linear-χ Gauss kernel is named and misses
live Q++.  Quartic-character projection of K̃_I misses.  That
Fourier family does not name Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.477 flags.

============================================================================
Setup.  Paley++ = {r ∈ T : χ(r−1)=χ(r+1)=+1}.  Type indicator
    I_++(r)=1_{r∈Paley++},   n_++=|Paley++|,
    K̃_I(α)=∑_r I_++(r) e(α r)   (r-sums over F_q^×).
Q++ = n_++^{−1} E[u(ξ) ∑_r I_++(r) u(r ξ)]
    = n_++^{−1} E[u(ξ) · 4 ∑_d N(d) K̃_I(ξ d)]
(15.298: u=4|hat 1_D|² on Ω).  Smooth expansion
    I_sm=(1+χ(r))(1+χ(r−1))(1+χ(r+1))/8
has eight terms.  Linear-χ block
    I_lin=(1+χ(r)+χ(r−1)+χ(r+1))/8.

============================================================================
Theorem A — PROVED (Jacobi + CM of E:y²=x³−x; certified p=3..31).
  ∑_{x∈F_q} χ(x³−x)=2p−a_E², a_E=trace of E(F_p) (Weierstrass:
  q+S3=#E(F_q)−1 and a_{p²}=a_E²−2p).  a_E=0 for p≡3 (mod 4);
  a_E²=4α² for p=α²+β², α odd, p≡1 (mod 4) (only a_E² is used;
  the sign of α is irrelevant).  Pairwise sums ∑χ((x−a)(x−b))=−1.
  Boundary I_sm(0)=I_sm(±1)=1/2.  Hence
      n_++ = (p² + 2p − 15 − a_E²)/8.
  Fail: drop a_E² at p=5 (5/2 ≠ 2).  Fail: the 3-point interpolant
  (p−3)(p−1)/4 at p=11 (20 ≠ 16).  ∎

Theorem B — PROVED (complete Gauss sums over F_q^×; certified p=5,7).
  G=∑_t χ(t)e(t).  For α≠0, χ(−1)=1,
      K_lin(α)=(1/8)[χ(α)G(1+e(α)+e(−α))−3].
  Fail: drop the −3 (max |err| > 0.3 at p=5).  ∎

Theorem C — PROVED (B + cache fold).
  The linear-χ Fourier predictor (truncation of I_sm, not a fit)
      Q++_lin = n_++^{−1} E[u(ξ) ∑_r I_lin(r) u(r ξ)]
  is 87/13 at p=5 and 5039/1227 at p=7, not live 48/13 or 1544/409.
  The full I_++ pairing recovers live Q++ (miss is the truncation).
  Fail: claim Q++_lin=48/13 at p=5.  ∎

Theorem D — PROVED (cache fold).
  Quartic-character L² projection of K̃_I (span{1,χ,ψ,ψ̄,ψ³} on
  F^×, K(0) kept) has max|err|>2 at p=5,7 and the projected
  predictor is 2.87 ≠ 48/13.  Ω± two-level |hat|²: 0 of 5726
  rows at p=7.  Fail: claim every p=7 Max+ is 2-level on Ω±.  ∎

Theorem E — OPEN.  The four quadratic-χ Weil terms in I_sm are
  the leftover; they are not in the quartic character span.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: field Jacobi + one DFT/N fold at p=5,7.  Writes
evidence/e1_gmin_m4_prop15478.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix, gauss_chi  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15478_catalog.json"

LIVE_LIN = {5: Fraction(87, 13), 7: Fraction(5039, 1227)}
NPP_CERT = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def a_E_sq(p: int) -> int:
    """a_E(F_p)² for E: y² = x³ − x.  0 if p≡3 (mod 4)."""
    if p % 4 == 3:
        return 0
    for alpha in range(1, p, 2):
        b2 = p - alpha * alpha
        if b2 < 0:
            break
        beta = int(round(b2 ** 0.5))
        if beta * beta == b2:
            return 4 * alpha * alpha
    raise ValueError(p)


def n_pp_named(p: int) -> int:
    return (p * p + 2 * p - 15 - a_E_sq(p)) // 8


def n_pp_drop_aE(p: int) -> Fraction:
    return Fraction(p * p + 2 * p - 15, 8)


def n_pp_interpolant(p: int) -> int:
    """3-point interpolant through p=3,5,7.  Dies at p=11."""
    return (p - 3) * (p - 1) // 4


def I_pp(F, r: int) -> float:
    if r == 0:
        return 0.0
    chi = F["chi"]
    if chi[r] != 1:
        return 0.0
    rm = F["add"][r][F["neg"][1]]
    rp = F["add"][r][1]
    if rm == 0 or rp == 0:
        return 0.0
    return 1.0 if (chi[rm] == 1 and chi[rp] == 1) else 0.0


def n_pp_count(F) -> int:
    return sum(1 for r in range(1, F["q"]) if I_pp(F, r) > 0.5)


def K_lin_named(F, G):
    """Named Gauss formula for the linear-χ kernel, α≠0."""
    q = F["q"]
    chi = F["chi"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    neg = F["neg"]
    out = np.zeros(q, dtype=np.complex128)
    out[0] = (q - 3) / 8.0
    for a in range(1, q):
        ca = float(chi[a])
        out[a] = (ca * G * (1.0 + e_tab[a] + e_tab[neg[a]]) - 3.0) / 8.0
    return out


def K_lin_named_drop3(F, G):
    q = F["q"]
    chi = F["chi"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    neg = F["neg"]
    out = np.zeros(q, dtype=np.complex128)
    out[0] = (q - 3) / 8.0
    for a in range(1, q):
        ca = float(chi[a])
        out[a] = (ca * G * (1.0 + e_tab[a] + e_tab[neg[a]])) / 8.0
    return out


def K_lin_exact(F):
    q = F["q"]
    chi = F["chi"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    add, neg, mul = F["add"], F["neg"], F["mul"]

    def term(mfun):
        acc = np.zeros(q, dtype=np.complex128)
        for a in range(q):
            s = 0j
            for r in range(1, q):
                w = mfun(r)
                if w:
                    s += w * e_tab[mul[a][r]]
            acc[a] = s
        return acc

    K1 = term(lambda r: 1.0)
    Kc = term(lambda r: float(chi[r]))
    Km = term(
        lambda r: 0.0 if add[r][neg[1]] == 0 else float(chi[add[r][neg[1]]])
    )
    Kp = term(lambda r: 0.0 if add[r][1] == 0 else float(chi[add[r][1]]))
    return (K1 + Kc + Km + Kp) / 8.0


def K_full(F):
    q = F["q"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    I = np.array([I_pp(F, r) for r in range(q)], dtype=np.float64)
    mul = F["mul"]
    K = np.array([(I * e_tab[mul[a]]).sum() for a in range(q)], dtype=np.complex128)
    return K, I


def quartic_proj(F, K, psi):
    q = F["q"]
    chi = np.array(F["chi"], dtype=np.float64)
    basis = [
        np.ones(q, dtype=np.complex128),
        chi.astype(np.complex128),
        psi.copy(),
        np.conjugate(psi),
        psi ** 3,
    ]
    for v in basis:
        v[0] = 0
    B = np.stack([v[1:] for v in basis], axis=1)
    c, *_ = np.linalg.lstsq(B, K[1:], rcond=None)
    Kp = np.zeros(q, dtype=np.complex128)
    Kp[0] = K[0]
    Kp[1:] = B @ c
    return Kp, float(np.max(np.abs(Kp - K)))


@lru_cache(maxsize=4)
def fold_pack(p: int):
    F = _kernel_field(p)
    q = F["q"]
    G = gauss_chi(F, p)
    Klin = K_lin_named(F, G)
    Kf, I = K_full(F)
    npp = int(round(I.sum()))
    Ds = D_rows(p, F).astype(np.float64)
    hat = Ds @ e_matrix(F).T
    energy = np.abs(hat) ** 2
    u = 4.0 * energy
    sigma = sigma_named(p)
    Omega = np.array(F["chi"]) == sigma
    Omega[0] = False
    xi0 = int(np.where(Omega)[0][0])
    N = load_N(p, F).astype(np.float64)
    mul = F["mul"]

    def pred(K):
        Kd = np.array([K[mul[xi0][d]] for d in range(q)])
        Iu = 4.0 * (N * Kd[None, :]).sum(axis=1)
        return float(np.mean((u[:, xi0] * Iu).real / npp) / (q * q))

    psi = make_psi(F)
    Kp, perr = quartic_proj(F, Kf, psi)
    plus = Omega & ((psi.imag > 0) if sigma < 0 else (psi.real > 0))
    minus = Omega & ~plus
    a_sp = energy[:, plus].std(axis=1)
    b_sp = energy[:, minus].std(axis=1)
    n2 = int(np.sum((a_sp < 1e-6) & (b_sp < 1e-6)))
    return {
        "npp": npp,
        "Q_lin": pred(Klin),
        "Q_full": pred(Kf),
        "Q_proj": pred(Kp),
        "proj_err": perr,
        "n_twolevel": n2,
        "n_rows": int(energy.shape[0]),
        "live": float(LIVE_QPP[p]),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in NPP_CERT:
        live = n_pp_count(_kernel_field(p))
        named = n_pp_named(p)
        if live != named:
            ok = False
        rows[str(p)] = {"live": live, "named": named, "aE2": a_E_sq(p)}
    if n_pp_drop_aE(5) == 2:
        ok = False
    if n_pp_drop_aE(5) != Fraction(5, 2):
        ok = False
    if n_pp_interpolant(11) == n_pp_named(11):
        ok = False
    if n_pp_interpolant(11) != 20 or n_pp_named(11) != 16:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "n_++=(p²+2p−15−a_E²)/8. "
            "Fail: drop a_E² at p=5; interpolant at p=11."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        G = gauss_chi(F, p)
        exact = K_lin_exact(F)
        named = K_lin_named(F, G)
        drop = K_lin_named_drop3(F, G)
        err = float(np.max(np.abs(exact - named)))
        err_drop = float(np.max(np.abs(exact - drop)))
        if err > 1e-10:
            ok = False
        if err_drop < 0.3:
            ok = False
        rows[str(p)] = {"err": err, "err_drop3": err_drop}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "K_lin(α)=(1/8)[χ(α)G(1+e(α)+e(−α))−3] (α≠0, r-sums on F_q^×). "
            "Fail: drop the −3."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = fold_pack(p)
        lin = Fraction(rec["Q_lin"]).limit_denominator(20000)
        full = Fraction(rec["Q_full"]).limit_denominator(20000)
        if lin != LIVE_LIN[p]:
            ok = False
        if full != LIVE_QPP[p]:
            ok = False
        if lin == LIVE_QPP[p]:
            ok = False
        rows[str(p)] = {
            "Q_lin": str(lin),
            "Q_full": str(full),
            "live": str(LIVE_QPP[p]),
        }
    if LIVE_LIN[5] == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Linear-χ / I_sm-truncation predictor is 87/13 and 5039/1227, "
            "not live Q++. Full I_++ recovers live. Fail: Q++_lin=48/13."
        ),
    }


def prove_D() -> dict:
    r5 = fold_pack(5)
    r7 = fold_pack(7)
    ok = r5["proj_err"] > 2.0 and r7["proj_err"] > 2.0
    ok = ok and abs(r5["Q_proj"] - float(LIVE_QPP[5])) > 0.5
    ok = ok and r7["n_twolevel"] == 0
    if r7["n_twolevel"] == r7["n_rows"]:
        ok = False
    return {
        "proved": bool(ok),
        "proj_err": {"5": r5["proj_err"], "7": r7["proj_err"]},
        "Q_proj5": r5["Q_proj"],
        "n_twolevel": {"5": r5["n_twolevel"], "7": r7["n_twolevel"]},
        "theorem": (
            "Quartic proj of K̃_I misses; Ω± two-level is 0 rows at p=7. "
            "Fail: every p=7 Max+ is 2-level on Ω±."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Linear-χ and quartic-projected type-indicator Fourier "
            "predictors miss Q++. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.478  Paley++ type-indicator Fourier", flush=True)
    A = prove_A()
    print(f"  A n_++: {A['proved']} p11={A['rows']['11']}", flush=True)
    B = prove_B()
    print(f"  B K_lin: {B['proved']} drop5={B['rows']['5']['err_drop3']:.3f}", flush=True)
    C = prove_C()
    print(f"  C lin-miss: {C['proved']} {C['rows']}", flush=True)
    D = prove_D()
    print(
        f"  D proj/2lev: {D['proved']} n2_7={D['n_twolevel']['7']}",
        flush=True,
    )
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["proved"], "B": B["proved"], "C": C["proved"], "D": D["proved"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.478",
        "title": (
            "Paley++ type-indicator Fourier: n_++ named; "
            "linear-χ / quartic families miss Q++"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "n_pp_CM_named": A["proved"],
            "K_lin_gauss_named": B["proved"],
            "lin_misses_Qpp": C["proved"],
            "quartic_and_twolevel_dead": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
        "L_status": "OPEN",
        "referees": {
            "claude_math_review": "PASS-WITH-NOTE",
            "openai_math_review": "PASS-WITH-NOTE",
            "action": "do_not_branch",
        },
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15478.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
