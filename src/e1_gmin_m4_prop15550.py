#!/usr/bin/env python3
"""
Prop 15.550 — Max+-free Kloosterman name of
    S(λ)=Σ_v χ(v²−1) e_λ(v):
    S(λ)= χ(−1) G₂(1)/G₁ · Kl(1, λ²/4) = Kl(1, λ²/4),
Kl(a,b)=Σ_{x≠0} e_1(a x + b x^{-1}).  Fail ±G(χ), χ(λ)G,
Kl(1,λ), Kl(1,λ²), Kl(1,−λ²/4).  K_λ^{all} is the W-rewrite
W=q^{-1}[p J(α)+Σ_Ω ŷ(η) J(α−η)] with J=e S, not a p-integer.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / Lemma D / 15.279–15.549 flags.
Does **not** import phi_F.  Does **not** overwrite 15.522–15.549.

============================================================================
Setup.  F_q, q=p², χ quadratic, e_1(x)=exp(2πi Tr(x)/p).
G₁=Σ_{x≠0} χ(x) e_1(x), G₂(α)=Σ_v e_1(α v²).  χ(z)=G₁^{-1}
Σ_t χ(t) e_1(−t z).  15.549: K_λ=K_λ^{all}−N q(3q−5);
K_λ^{all}=Σ_{a,b,c,d} χ((d−a)(d−b)(c−a)(c−b)) T_{abcd}
e_ξ(a−b) e_{rξ}(c−d).  J(a,b;γ)=Σ_t χ((t−a)(t−b)) e_γ(t).
W(a,b;α)=Σ_t y_t χ((t−a)(t−b)) e_α(t).  Cy=py ⇒ ŷ(0)=p and
ŷ vanishes off {0}∪Ω, Ω={ξ: G(χ,ξ)=p}.

============================================================================
Theorem A — PROVED (Gauss inversion + complete square; p=5..19).
  Inner Gauss Σ_v e_1(−t v²+λ v)=e_1(λ²/(4t)) G₂(−t) and
  G₂(−t)=χ(−t) G₂(1) give
      S(λ)= χ(−1) G₂(1)/G₁ · Kl(1, λ²/4).
  G₂(1)=G₁ by Σ_v e(v²)=Σ_x(1+χ(x))e(x)=G₁ (no ±p needed)
  and χ(−1)=1 on F_{p²}, so the prefactor is +1 and S=Kl.
  The same cancellation actually gives S=Kl on every odd q.
  S(0)=Kl(1,0)=Σ_{x≠0}e_1(x)=−1, not q−1 (that is Kl(0,0)).
  J(a,a;γ)=−e_γ(a) only for γ≠0; J(a,a;0)=q−1.  Fail:
  Kl(1, −λ²/4).  ∎

Theorem B — PROVED (fail-when-wrong dies; p=5,7,11,13,17,19).
  S≢G(χ), S≢−G(χ), S≢χ(λ)G(χ), S≢Kl(1,λ), S≢Kl(1,λ²).
  Each max |err| ≥ 3.  ∎

Theorem C — OPEN.  For a≠b,
    J(a,b;γ)=e_γ((a+b)/2) S(γ(a−b)/2)=e_γ((a+b)/2) Kl(1, γ²(a−b)²/16).
  Fourier inversion + Cy=py:
    W(a,b;α)=q^{-1}[ p J(α) + Σ_{η∈Ω} ŷ(η) J(α−η) ].
  Then K_λ^{all}=Σ_y Σ_{a,b} y_a y_b e_ξ(a−b) W(rξ) W(−rξ)
  reconstructs 249050 at p=5.  The ŷ=0 term is
  (p/q)² Σ G(a,b) e_ξ(a−b) J(rξ)J(−rξ)=650, not the bulk.
  Fail the frequency α+η (p=5 max |W err|≥30).  K_all/N equals
  24905/13 (p=5) and 3003413/409 (p=7): Max+ denominators,
  not a Gauss/Jacobi integer in p.  16pA / Q_τ unnamed.

============================================================================
Backend: field tables O(q²); p=5 W-rewrite on the Max+ cache.
Writes evidence/e1_gmin_m4_prop15550.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field

EV = ROOT / "evidence" / "e1_gmin_m4_prop15550.json"
YP5 = Path("/tmp/maxplus_p5.npy")

PRIMES = (5, 7, 11, 13, 17, 19)


def _half(F) -> int:
    add = F["add"]
    for t in range(F["q"]):
        if add[t][t] == 1:
            return t
    raise RuntimeError("no 1/2")


def _e_tab(F) -> np.ndarray:
    return np.array([_e_tr(F, t) for t in range(F["q"])], dtype=np.complex128)


def kloosterman(F, a: int, b: int) -> complex:
    q = F["q"]
    add, mul, inv = F["add"], F["mul"], F["inv"]
    e = _e_tab(F)
    acc = 0j
    for x in range(1, q):
        acc += e[add[mul[a][x]][mul[b][inv[x]]]]
    return acc


def S_of(F, lam: int) -> complex:
    q = F["q"]
    chi, add, mul, neg = F["chi"], F["add"], F["mul"], F["neg"]
    e = _e_tab(F)
    acc = 0j
    for v in range(q):
        acc += chi[add[mul[v][v]][neg[1]]] * e[mul[lam][v]]
    return acc


def S_named(F, lam: int) -> complex:
    mul = F["mul"]
    inv4 = mul[_half(F)][_half(F)]
    return kloosterman(F, 1, mul[mul[lam][lam]][inv4])


def S_named_pref(F, lam: int) -> complex:
    return pref_S(F) * S_named(F, lam)


def S_wrong_sign(F, lam: int) -> complex:
    """Fail-when-wrong: Kl(1, −λ²/4)."""
    mul, neg = F["mul"], F["neg"]
    inv4 = mul[_half(F)][_half(F)]
    return kloosterman(F, 1, neg[mul[mul[lam][lam]][inv4]])


def S_wrong_drop4(F, lam: int) -> complex:
    """Fail-when-wrong: drop the 1/4, Kl(1, λ²)."""
    return kloosterman(F, 1, F["mul"][lam][lam])


def S_wrong_lam(F, lam: int) -> complex:
    """Fail-when-wrong: Kl(1, λ)."""
    return kloosterman(F, 1, lam)


def Gchi1(F) -> complex:
    e = _e_tab(F)
    return sum(F["chi"][x] * e[x] for x in range(1, F["q"]))


def G2(F, alpha: int) -> complex:
    e = _e_tab(F)
    mul = F["mul"]
    return sum(e[mul[alpha][mul[v][v]]] for v in range(F["q"]))


def pref_S(F) -> complex:
    return F["chi"][F["neg"][1]] * G2(F, 1) / Gchi1(F)


def omega_set_local(F) -> list[int]:
    q, p = F["q"], F["p"]
    om = []
    for xi in range(1, q):
        chat = sum(
            F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)
        )
        if abs(chat.real - p) < 1e-6 and abs(chat.imag) < 1e-6:
            om.append(xi)
    return om


def pp_r(F) -> int:
    add, neg, chi = F["add"], F["neg"], F["chi"]
    for r in F["squares"]:
        rm, rp = add[r][neg[1]], add[r][1]
        if rm and rp and chi[rm] == 1 and chi[rp] == 1:
            return int(r)
    raise RuntimeError("no ++ r")


def J_raw(F, a: int, b: int, gam: int) -> complex:
    q = F["q"]
    chi, add, mul, neg = F["chi"], F["add"], F["mul"], F["neg"]
    e = _e_tab(F)
    acc = 0j
    for t in range(q):
        acc += chi[mul[add[t][neg[a]]][add[t][neg[b]]]] * e[mul[gam][t]]
    return acc


def J_named(F, a: int, b: int, gam: int) -> complex:
    """e_γ((a+b)/2) S(γ(a−b)/2) off-diagonal; diagonal Gauss."""
    if a == b:
        phase = _e_tr(F, F["mul"][gam][a])
        return phase * ((F["q"] - 1) if gam == 0 else -1)
    h = _half(F)
    mid = F["mul"][F["add"][a][b]][h]
    arg = F["mul"][gam][F["mul"][F["add"][a][F["neg"][b]]][h]]
    return _e_tr(F, F["mul"][gam][mid]) * S_of(F, arg)


def load_Yp5(q: int) -> np.ndarray:
    raw = np.load(YP5)
    block = raw[raw[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)
    return np.where(np.sign(block) < 0, -1.0, 1.0)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in PRIMES:
        F = _kernel_field(p)
        q = F["q"]
        g1, g2 = Gchi1(F), G2(F, 1)
        pref = pref_S(F)
        chi_m1 = int(F["chi"][F["neg"][1]])
        s0, k0 = S_of(F, 0), kloosterman(F, 1, 0)
        err_K = err_pref = 0.0
        for lam in range(q):
            got = S_of(F, lam)
            named = S_named(F, lam)
            err_K = max(err_K, abs(got - named))
            err_pref = max(err_pref, abs(got - pref * named))
        row_ok = (
            err_K < 1e-8
            and err_pref < 1e-8
            and abs(g1 - g2) < 1e-8
            and abs(pref - 1) < 1e-8
            and chi_m1 == 1
            and abs(s0 + 1) < 1e-8
            and abs(k0 + 1) < 1e-8
        )
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "q": q,
            "G1_re": float(np.real(g1)),
            "G2_re": float(np.real(g2)),
            "chi_m1": chi_m1,
            "pref_re": float(np.real(pref)),
            "S0_re": float(np.real(s0)),
            "maxerr_S_eq_K": float(err_K),
            "maxerr_S_eq_pref_K": float(err_pref),
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "S(λ)=χ(−1)G₂(1)/G₁·Kl(1,λ²/4)=Kl(1,λ²/4). "
            "Fail Kl(1,−λ²/4)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in PRIMES:
        F = _kernel_field(p)
        q = F["q"]
        g1 = Gchi1(F)
        err = {
            "G": 0.0,
            "mG": 0.0,
            "chiG": 0.0,
            "drop4": 0.0,
            "lam": 0.0,
            "sign": 0.0,
        }
        for lam in range(q):
            got = S_of(F, lam)
            err["G"] = max(err["G"], abs(got - g1))
            err["mG"] = max(err["mG"], abs(got + g1))
            if lam != 0:
                err["chiG"] = max(err["chiG"], abs(got - F["chi"][lam] * g1))
            err["drop4"] = max(err["drop4"], abs(got - S_wrong_drop4(F, lam)))
            err["lam"] = max(err["lam"], abs(got - S_wrong_lam(F, lam)))
            err["sign"] = max(err["sign"], abs(got - S_wrong_sign(F, lam)))
        row_ok = all(v > 1.0 for v in err.values())
        if not row_ok:
            ok = False
        rows[str(p)] = {k: float(v) for k, v in err.items()}
        rows[str(p)]["ok"] = bool(row_ok)
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "S≢±G(χ), S≢χ(λ)G, S≢Kl(1,λ), S≢Kl(1,λ²), "
            "S≢Kl(1,−λ²/4)."
        ),
    }


def check_J_sample(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    Om = omega_set_local(F)
    ns = next(x for x in range(1, q) if F["chi"][x] == -1)
    maxerr = 0.0
    n = 0
    for a in (0, 1, 2):
        for b in range(min(q, 8)):
            for gam in (0, 1, Om[0], ns):
                maxerr = max(maxerr, abs(J_raw(F, a, b, gam) - J_named(F, a, b, gam)))
                n += 1
    return {"p": p, "n": n, "maxerr": float(maxerr), "ok": bool(maxerr < 1e-8)}


def check_rewrite_p5() -> dict:
    p = 5
    F = _kernel_field(p)
    q = F["q"]
    Y = load_Yp5(q)
    N = Y.shape[0]
    Omega = omega_set_local(F)
    xi = int(Omega[0])
    eta = int(F["mul"][pp_r(F)][xi])
    add = np.array(F["add"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    mul = np.array(F["mul"], dtype=np.int32)
    chi = np.array(F["chi"], dtype=np.float64)
    e_tab = np.array(
        [_e_tab(F)[mul[al]] for al in range(q)], dtype=np.complex128
    )
    ta = add[:, neg]
    prod = mul[ta[:, :, None], ta[:, None, :]]
    Q = chi[prod].transpose(1, 2, 0)
    yhat = {int(om): Y @ e_tab[om] for om in Omega}
    off = [
        float(np.max(np.abs(Y @ e_tab[al])))
        for al in range(1, q)
        if al not in set(Omega)
    ]

    def J_at(gam: int) -> np.ndarray:
        return (Q * e_tab[gam][None, None, :]).sum(axis=2)

    def W_direct(alpha: int) -> np.ndarray:
        kern = Q * e_tab[alpha][None, None, :]
        return (Y @ kern.reshape(q * q, q).T).reshape(N, q, q)

    def W_pred(alpha: int, plus: bool) -> np.ndarray:
        acc = float(p) * J_at(alpha)
        for om in Omega:
            gam = int(F["add"][alpha][om if plus else F["neg"][om]])
            acc = acc + yhat[om][:, None, None] * J_at(gam)[None, :, :]
        return acc / q

    Wr = W_direct(eta)
    Wm = W_direct(int(F["neg"][eta]))
    err_minus = float(np.max(np.abs(Wr - W_pred(eta, False))))
    err_plus = float(np.max(np.abs(Wr - W_pred(eta, True))))
    e_xi = e_tab[xi]
    Kall = 0j
    for a in range(q):
        for b in range(q):
            phase = e_xi[a] * np.conj(e_xi[b])
            Kall += phase * np.sum(Y[:, a] * Y[:, b] * Wr[:, a, b] * Wm[:, a, b])
    Jr, Jm = J_at(eta), J_at(int(F["neg"][eta]))
    G = Y.T @ Y
    term0 = 0j
    for a in range(q):
        for b in range(q):
            term0 += G[a, b] * e_xi[a] * np.conj(e_xi[b]) * Jr[a, b] * Jm[a, b]
    term0 *= (p * p) / (q * q)
    return {
        "N": int(N),
        "yhat0_err": float(np.max(np.abs(Y.sum(axis=1) - p))),
        "max_off_Omega": float(max(off) if off else 0.0),
        "err_W_aminus": err_minus,
        "err_W_aplus": err_plus,
        "K_all": float(Kall.real),
        "term0": float(term0.real),
        "match_249050": bool(abs(Kall.real - 249050) < 1e-4),
        "term0_not_bulk": bool(abs(term0.real - 249050) > 1.0),
        "plus_fails": bool(err_plus > 1.0),
        "minus_ok": bool(err_minus < 1e-8),
        "p5_Kall_over_N": "24905/13",
        "p7_Kall_over_N": "3003413/409",
    }


def prove_open() -> dict:
    Jrows = {str(p): check_J_sample(p) for p in (5, 7, 11, 13)}
    rew = check_rewrite_p5()
    j_ok = all(r["ok"] for r in Jrows.values())
    rewrite_ok = (
        rew["match_249050"]
        and rew["term0_not_bulk"]
        and rew["plus_fails"]
        and rew["minus_ok"]
        and rew["yhat0_err"] < 1e-8
        and rew["max_off_Omega"] < 1e-8
    )
    return {
        "proved": False,
        "K_all_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "J_named_ok": bool(j_ok),
        "rewrite_ok": bool(rewrite_ok),
        "J": Jrows,
        "rewrite": rew,
        "note": (
            "W=q^{-1}[pJ+Σ_Ω ŷ J] and J=e S reconstruct K_λ^{all} "
            "(p=5: 249050) but K_all/N has Max+ denominators 13,409, "
            "not a Gauss/Jacobi integer in p. 16pA / Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.550  S(λ)=Kl(1, λ²/4); K_λ^{all} rewrite open", flush=True)
    A = prove_A()
    print(f"  A S=Kl: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B fail-when-wrong: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C K_all named: {C['K_all_named_in_p']}", flush=True)
    out = {
        "prop": "15.550",
        "title": "S(λ)=Kl(1, λ²/4); K_λ^{all} is a W-rewrite not a p-law",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "S_is_kloosterman": A["proved"],
            "S_fail_when_wrong": B["proved"],
            "K_all_named_in_p": False,
            "NUM_SUM_16pA_general": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": "S(λ)=χ(−1)G2(1)/G1·Kl(1,λ²/4)=Kl(1,λ²/4)",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
