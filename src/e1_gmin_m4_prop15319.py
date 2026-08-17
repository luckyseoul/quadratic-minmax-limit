#!/usr/bin/env python3
"""
Prop 15.319 — Ω is the Hasse–Davenport χ-coset; C(1) is named;
∑_a Kl(a)Kl(sa) is elementary; A_□ is the named-Kl image of the
U-form factors.  A_□ / Q_τ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.318 flags.  Floor stays OPEN.

============================================================================
Setup.  q=p².  Ω={ξ: χ̂(ξ)=p}.  K(α)=∑_{u∈U} e(Tr(αu))=−Kl(−N(α),−1)
(15.316).  Kl(n):=Kl(1,n)=∑_{x≠0} e(x+n/x) (real).  C(r)=∑_{ξ∈Ω} K(ξ)K(rξ)
depends only on s=N(r)∈F_p^{×□}.  W(1)=∑_{N(δ)=1} N(δ), A_□=E[W(1)²].
ε=−χ_p(−1).  N̂=u/4 on Ω.

============================================================================
Theorem A — PROVED (Hasse–Davenport + quadratic Gauss; certified p=5,7,11).
  χ̂(ξ)=G_q χ_q(ξ) for ξ≠0, and G_q=−G_p(χ_p)²=−χ_p(−1)p.  Hence
      Ω={ξ: χ_q(ξ)=−χ_p(−1)}
  (nonsquares if p≡1 (mod 4), squares if p≡3).  Fail-when-wrong:
  claim Ω is always the squares (false at p=5).  ∎

Theorem B — PROVED (complete exponential sums on F_p; certified p=5,7,11).
  ∑_{a≠0} Kl(a)Kl(sa) = p²−p−1 if s=1, and =−(p+1) if s≠1.
  Off-diagonal contribution is p−2 (s=1) and −2 (s≠1); plus the
  y=−sx diagonal this is the displayed total.
  ∑_{a≠0} χ_p(a) Kl(a)² = −p  (J(χ,χ)=−χ(−1), not −χ(−1)p, because
  χ² is trivial so the usual Jacobi formula does not apply).
  Hence C(1)=(p+1)(p²−p−1−ε p)/2 is named in p.  Fail-when-wrong:
  C(1)=|Ω|; or C(s) constant for all square s≠1 (false at p=11:
  values 192 and −336).  ∎

Theorem C — PROVED as a Fourier identity; A_□ OPEN as a number
  (certified p=5,7, err=0).
  W(1)=q^{-1}[(p+1)k² + ∑_{ξ∈Ω} K(ξ) N̂(ξ)], and ∑_Ω K=|Ω|
  (average K on Ω is 1).  Expanding E[W²] and using
  E[N̂(ξ)N̂(rξ)]=Q(r)/16 on Ω gives
      A_□ = (p+1)² k⁴/q² + (p+1) k² |Ω|/q
            + (1/16) ∑_{c∈F_p^{×□}} C(c) F_{√c}(0).
  At p=5,7 this is α+β F_1(0)+γ F_*(0) with named Kl coefficients
  (C takes two values).  At p≥11, C takes ≥3 values, so the third
  term is not a 2-class form.  No lin/lin or named-feature formula
  hits 12360/13 and 2907408/409.  Q_τ unnamed.  phi_F not imported.
  Fail-when-wrong: drop the third term (err = third ≠ 0).  ∎

============================================================================
Backend: F_p Kloosterman (Max+-free) + one N/Q fold per prime for C.
Writes evidence/e1_gmin_m4_prop15319.json
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
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm, live_Q_on_T, omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15302 import form_factor  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm, kloosterman  # noqa: E402

_ID_CACHE: dict[int, dict] = {}


def omega_chi_named(p: int) -> int:
    return -chi_p_pm(p, p - 1)


def kl1(p: int, n: int) -> complex:
    return kloosterman(p, 1, n % p)


def kl_conv_full_named(p: int, s: int) -> int:
    if s % p == 1:
        return p * p - p - 1
    return -(p + 1)


def chi_kl2_named(p: int) -> int:
    return -p


def C1_named(p: int) -> int:
    eps = omega_chi_named(p)
    return (p + 1) * (p * p - p - 1 + eps * chi_kl2_named(p)) // 2


def C_from_kl(p: int, s: int) -> float:
    """C(s)=(p+1) ∑_{χ(a)=ε} Kl(a)Kl(sa)."""
    eps = omega_chi_named(p)
    acc = 0j
    for a in range(1, p):
        if chi_p_pm(p, a) != eps:
            continue
        acc += kl1(p, a) * kl1(p, (s * a) % p)
    return float(((p + 1) * acc).real)


def sqrt_fp(p: int, c: int) -> int:
    for x in range(1, p):
        if (x * x) % p == c % p:
            return x
    raise ValueError(f"no sqrt of {c} mod {p}")


def live_kl_moments(p: int) -> dict:
    s2 = 0j
    sx = 0j
    conv = {}
    for a in range(1, p):
        ka = kl1(p, a)
        s2 += ka * ka
        sx += chi_p_pm(p, a) * ka * ka
    for s in range(1, p):
        acc = 0j
        for a in range(1, p):
            acc += kl1(p, a) * kl1(p, (s * a) % p)
        conv[s] = float(acc.real)
    return {
        "sum_kl2": float(s2.real),
        "sum_chi_kl2": float(sx.real),
        "conv": conv,
    }


def prove_omega() -> dict:
    ok = True
    always_square_hit = 0
    rows = {}
    for p in (5, 7, 11):
        F = _kernel_field(p)
        Omega = omega_set(F)
        want = omega_chi_named(p)
        chis = [int(F["chi"][xi]) for xi in Omega]
        match = all(c == want for c in chis) and len(Omega) == (F["q"] - 1) // 2
        if not match:
            ok = False
        if all(c == 1 for c in chis):
            always_square_hit += 1
        rows[str(p)] = {
            "want": want,
            "Om": len(Omega),
            "all_want": match,
            "all_squares": all(c == 1 for c in chis),
        }
    if always_square_hit == 3:
        ok = False
    if rows["5"]["all_squares"]:
        ok = False
    return {
        "proved": ok,
        "always_squares_false": not rows["5"]["all_squares"],
        "by_p": rows,
        "theorem": "Ω={χ_q=−χ_p(−1)} via G_q=−χ_p(−1)p.",
    }


def prove_kl_conv() -> dict:
    ok = True
    cstar_const_hit = 0
    rows = {}
    for p in (5, 7, 11):
        mom = live_kl_moments(p)
        if abs(mom["sum_kl2"] - (p * p - p - 1)) > 1e-6:
            ok = False
        if abs(mom["sum_chi_kl2"] - chi_kl2_named(p)) > 1e-6:
            ok = False
        conv_err = 0.0
        for s, live in mom["conv"].items():
            conv_err = max(conv_err, abs(live - kl_conv_full_named(p, s)))
        if conv_err > 1e-6:
            ok = False
        cvals = {c: C_from_kl(p, c) for c in range(1, p) if chi_p_pm(p, c) == 1}
        others = [v for c, v in cvals.items() if c != 1]
        spread = (max(others) - min(others)) if others else 0.0
        if spread < 1.0:
            cstar_const_hit += 1
        c1_live = cvals[1]
        c1_om = (p * p - 1) // 2  # |Ω|
        rows[str(p)] = {
            "sum_kl2": mom["sum_kl2"],
            "sum_chi_kl2": mom["sum_chi_kl2"],
            "conv_err": conv_err,
            "C1_live": c1_live,
            "C1_named": C1_named(p),
            "C1_vs_Om": c1_om,
            "C_off1": others,
            "Cstar_spread": spread,
        }
        if abs(c1_live - C1_named(p)) > 1e-6:
            ok = False
        if abs(c1_live - c1_om) < 1e-6:
            ok = False
    if cstar_const_hit == 3:
        ok = False
    if rows["11"]["Cstar_spread"] < 1.0:
        ok = False
    return {
        "proved": ok,
        "C1_eq_Omega_false": abs(rows["5"]["C1_live"] - 12) > 1.0,
        "Cstar_constant_false": rows["11"]["Cstar_spread"] > 1.0,
        "by_p": rows,
        "theorem": (
            "∑ Kl(a)Kl(sa) is p²−p−1 or −(p+1); ∑ χ Kl²=−p; "
            "C(1) named; C off 1 not constant."
        ),
    }


def _identity_one(p: int) -> dict:
    if p in _ID_CACHE:
        return _ID_CACHE[p]
    F = _kernel_field(p)
    q = F["q"]
    k = p * (p - 1) / 2.0
    Q = live_Q_on_T(p)
    N = load_N(p, F).astype(np.float64)
    Uidx = [d for d in range(1, q) if int(field_norm(F, d)) == 1]
    W = N[:, np.array(Uidx, dtype=np.int32)].sum(axis=1)
    live = float(np.mean(W**2))
    term0 = (p + 1) ** 2 * k**4 / float(q * q)
    Om = (q - 1) // 2
    term1 = (p + 1) * k**2 * Om / float(q)
    third = 0.0
    Cmap = {}
    for c in range(1, p):
        if chi_p_pm(p, c) != 1:
            continue
        Cc = C_from_kl(p, c)
        Cmap[c] = Cc
        scale = sqrt_fp(p, c)
        Fc = form_factor(F, Q, scale, 0).real
        third += Cc * Fc / 16.0
    pred = term0 + term1 + third
    hist = defaultdict(int)
    for x in W:
        hist[int(round(float(x)))] += 1
    out = {
        "p": p,
        "live": live,
        "pred": pred,
        "err": pred - live,
        "term0": term0,
        "term1": term1,
        "third": third,
        "C": Cmap,
        "W_hist": dict(sorted(hist.items())),
        "n_rows": int(N.shape[0]),
        "W_mean": float(np.mean(W)),
        "W_var": float(np.var(W)),
    }
    _ID_CACHE[p] = out
    return out


def prove_identity() -> dict:
    ok = True
    drop_hit = 0
    rows = {}
    for p in (5, 7):
        rec = _identity_one(p)
        if abs(rec["err"]) > 1e-7:
            ok = False
        dropped = rec["term0"] + rec["term1"]
        if abs(dropped - rec["live"]) > 1.0:
            drop_hit += 1
        rows[str(p)] = {
            "live": rec["live"],
            "pred": rec["pred"],
            "err": rec["err"],
            "third": rec["third"],
            "W_hist": rec["W_hist"],
            "nuniq_W": len(rec["W_hist"]),
        }
    if drop_hit == 0:
        ok = False
    if rows["5"]["nuniq_W"] <= 2:
        ok = False
    return {
        "proved": ok,
        "drop_third_fails": drop_hit > 0,
        "W_not_two_valued": rows["5"]["nuniq_W"] > 2,
        "by_p": rows,
        "theorem": (
            "A_□ = named(p) + (1/16)∑_{c□} C(c) F_{√c}(0). "
            "Drop third fails. W(1) is not 2-valued."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "A_□ is a named-Kl linear image of the U-form factors. "
            "C(s) for s≠1 still carries the unnamed χ-Kloosterman "
            "convolution.  Q_τ / phi_F stay OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.319  Ω / C(1) / Kl-conv named; A_□ = Kl-image of F_c", flush=True)
    A = prove_omega()
    print(f"  A Omega: {A['proved']} always_sq={not A['always_squares_false']}", flush=True)
    B = prove_kl_conv()
    print(
        f"  B Kl-conv/C1: {B['proved']} C1=Om={not B['C1_eq_Omega_false']} "
        f"Cstar_const={not B['Cstar_constant_false']}",
        flush=True,
    )
    C = prove_identity()
    print(
        f"  C identity: {C['proved']} drop3={C['drop_third_fails']} "
        f"W_spread={C['W_not_two_valued']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D A_□/Q: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.319",
        "title": "Ω and C(1) named; A_□ is Kl-linear in F_c; Q_τ open",
        "proved": {
            "omega_HD": A["proved"],
            "kl_conv_and_C1": B["proved"],
            "Asq_identity": C["proved"],
            "A_square_named": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15319.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
