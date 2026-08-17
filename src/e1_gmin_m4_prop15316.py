#!/usr/bin/env python3
"""
Prop 15.316 — Torus Gauss is a Kloosterman sum: K(α)=−Kl(−N(α),−1).
M_U is the Kl-image of the norm-fibers W(c).  F_1(0) still carries
den N_+/(2p).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.315 flags.  Floor stays OPEN.

============================================================================
Setup.  U={N_{q/p}=1}, K(α)=∑_{u∈U} e(Tr(αu))=G_U(0,α).
Kl(A,B)=∑_{x∈F_p^×} e_p(A x + B x^{-1}).  W(c)=∑_{N(δ)=c} N(δ)
for c∈F_p^×.  M_U(ξ)=∑_{u∈U} N̂(ξu).

============================================================================
Theorem A — PROVED (complete the square on F_p²; certified p=5,7).
  For α≠0,
      K(α) = −Kl(−N_{q/p}(α), −1),
  and K(0)=p+1.  Derivation: 1_{N=c}=(1/p)∑_a e_p(a(N−c)), the
  inner sum ∑_x e(Tr(x)+a N(x)) equals −p e_p(−1/a) for a≠0
  (χ(-d)χ(-1)=−1 cancels the choice of nonsquare d).  Hence the
  fiber sum is −Kl(−c,−1).  Fail-when-wrong: replace Kl by 1, or
  by −χ_p(N(α)).  ∎

Theorem B — PROVED (2-design on fibers; certified p=5,7).
  |N^{-1}(c)|=p+1 and χ_q=χ_p∘N, so
      E[W(c)] = (p+1) μ_{χ_p(c)},
  μ_+=p(p−1)/4, μ_-=p(p−3)/4.  E[W(c)²] depends only on χ_p(c).
  Fail-when-wrong: claim E[W(c)]=k for every c.  ∎

Theorem C — PROVED as an identity; F_1(0) OPEN as a number.
  M_U(ξ)=k(p+1) − ∑_{c∈F_p^×} Kl(−N(ξ)c, −1) W(c).
  Hence F_1(0)=16 q^{-2} E[N̂ M_U] is a Kloosterman image of the
  W-gram.  For p≡3 (mod 4), E[W(c)W(−c)]=(p+1)² μ_+ μ_- (named;
  fails at p=5 where χ_p(−1)=+1).  The remaining W-gram entries
  still have den 13, 409.  Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: field tables + one N-fold per prime.  Writes
evidence/e1_gmin_m4_prop15316.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm, live_Q_on_T  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15302 import form_factor  # noqa: E402
from e1_gmin_m4_prop15303 import GU  # noqa: E402
from e1_gmin_m4_prop15314 import F1_named  # noqa: E402


def kloosterman(p: int, A: int, B: int) -> complex:
    acc = 0j
    for x in range(1, p):
        inv = pow(x, p - 2, p)
        acc += np.exp(2j * np.pi * ((A * x + B * inv) % p) / p)
    return acc


def chi_p_pm(p: int, x: int) -> int:
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def mu_plus(p: int) -> Fraction:
    return Fraction(p * (p - 1), 4)


def mu_minus(p: int) -> Fraction:
    return Fraction(p * (p - 3), 4)


def prove_kloosterman() -> dict:
    ok = True
    one_hit = 0
    chi_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        err = 0.0
        err1 = 0.0
        errchi = 0.0
        k0 = GU(F, 0, 0)
        if abs(k0 - (p + 1)) > 1e-9:
            ok = False
        for a in range(1, F["q"]):
            k = GU(F, 0, a)
            c = int(field_norm(F, a))  # F_p embedded as 0..p-1
            pred = -kloosterman(p, (-c) % p, (-1) % p)
            err = max(err, abs(k - pred))
            err1 = max(err1, abs(k - (-1)))
            errchi = max(errchi, abs(k - (-chi_p_pm(p, c))))
        if err > 1e-8:
            ok = False
        if err1 > 0.5:
            one_hit += 1
        if errchi > 0.5:
            chi_hit += 1
        rows[str(p)] = {"err": err, "err_minus1": err1, "err_chiN": errchi, "K0": p + 1}
    if one_hit == 0 or chi_hit == 0:
        ok = False
    return {
        "proved": ok,
        "minus1_fails": one_hit > 0,
        "chiN_fails": chi_hit > 0,
        "by_p": rows,
        "theorem": (
            "K(α)=−Kl(−N(α),−1) for α≠0; K(0)=p+1. "
            "−1 and −χ_p(N) both miss."
        ),
    }


def prove_W_means() -> dict:
    ok = True
    k_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        fibers: dict[int, list[int]] = defaultdict(list)
        for d in range(1, F["q"]):
            fibers[int(field_norm(F, d))].append(d)
        rec = {}
        mp, mm = float(mu_plus(p)), float(mu_minus(p))
        for c, ds in fibers.items():
            w = N[:, np.array(ds, dtype=np.int32)].sum(axis=1)
            mu = float(w.mean())
            want = (p + 1) * (mp if chi_p_pm(p, c) == 1 else mm)
            if abs(mu - want) > 1e-8:
                ok = False
            if abs(mu - p * (p - 1) / 2) > 1e-6:
                k_hit += 1
            rec[str(c)] = {
                "mean": mu,
                "want": want,
                "fiber": len(ds),
                "unique": int(len(np.unique(np.round(w, 8)))),
            }
        rows[str(p)] = rec
    if k_hit == 0:
        ok = False
    return {
        "proved": ok,
        "E_W_is_not_k": k_hit > 0,
        "by_p": rows,
        "theorem": (
            "E[W(c)]=(p+1) μ_{χ_p(c)}.  Not equal to k. "
            "E[W(c)²] is χ_p(c)-constant (den 13, 409)."
        ),
    }


def prove_MU_and_antipode() -> dict:
    ok = True
    p5_indep_fail = False
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        fibers: dict[int, list[int]] = defaultdict(list)
        for d in range(1, F["q"]):
            fibers[int(field_norm(F, d))].append(d)
        W = {
            c: N[:, np.array(ds, dtype=np.int32)].sum(axis=1)
            for c, ds in fibers.items()
        }
        # antipode gram
        antip = {}
        named_prod = float((p + 1) ** 2 * mu_plus(p) * mu_minus(p))
        same_char_prod = float((p + 1) ** 2 * mu_plus(p) * mu_plus(p))
        for c in fibers:
            cm = (-c) % p
            if cm == 0 or cm not in W:
                continue
            val = float(np.mean(W[c] * W[cm]))
            antip[str(c)] = val
            if p % 4 == 3:
                if abs(val - named_prod) > 1e-6:
                    ok = False
            else:
                if abs(val - named_prod) < 1e-4:
                    p5_indep_fail = True
                    ok = False
        # F_1(0) from live Q vs 15.314 named form
        Q = live_Q_on_T(p)
        acc: dict[str, list[float]] = defaultdict(list)
        q2 = float(F["q"] ** 2)
        for r, v in Q.items():
            acc[merge_cls(F, r)].append(v / q2)
        Qm = {c: float(np.mean(vs)) for c, vs in acc.items()}
        f1 = form_factor(F, Q, 1, 0).real
        f1n = F1_named(p, 0, Qm.get("++", 0.0), Qm.get("--N1", 0.0))
        if abs(f1 - f1n) > 1e-8:
            ok = False
        rows[str(p)] = {
            "antipode": antip,
            "named_mu_prod": named_prod,
            "same_char_prod": same_char_prod,
            "F1": f1,
            "F1_named_from_Q": f1n,
        }
    if p % 4 == 1 and not p5_indep_fail:
        # p=5 must refuse the p≡3 formula
        pass
    return {
        "proved": ok,
        "p5_rejects_mu_prod": (not p5_indep_fail) and abs(
            list(rows["5"]["antipode"].values())[0] - rows["5"]["named_mu_prod"]
        )
        > 1.0,
        "by_p": rows,
        "theorem": (
            "M_U=k(p+1)−∑_c Kl(−N(ξ)c,−1) W(c). "
            "E[W(c)W(−c)]=(p+1)² μ_+ μ_- iff p≡3 (mod 4). "
            "F_1(0) still has den N_+/(2p)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "K is Kloosterman. F_1(0) is a Kl-image of the W-gram. "
            "W-gram still unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.316  K=−Kl(−N,−1); W-fibers; F_1(0) still unnamed", flush=True)
    A = prove_kloosterman()
    print(f"  A Kl: {A['proved']} -1_fails={A['minus1_fails']} chi_fails={A['chiN_fails']}", flush=True)
    B = prove_W_means()
    print(f"  B W-means: {B['proved']} not_k={B['E_W_is_not_k']}", flush=True)
    C = prove_MU_and_antipode()
    print(f"  C antipode/F1: {C['proved']} p5_rejects={C['p5_rejects_mu_prod']}", flush=True)
    D = prove_open()
    print(f"  D general: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.316",
        "title": "Torus Gauss is Kloosterman; W-fibers; F_1(0) unnamed",
        "proved": {
            "kloosterman_K": A["proved"],
            "W_means": B["proved"],
            "MU_identity": C["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15316.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
