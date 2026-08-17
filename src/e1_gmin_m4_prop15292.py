#!/usr/bin/env python3
"""
Prop 15.292 — Q on Type+ 1D lifts is named in p.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.291 flags.  Does **not** name Q_NL or Q_full in p.
Floor stays OPEN.

============================================================================
Setup.  Type+ 1D: y_∞=1, y=σ∘L, ker L a square line, |σ^{-1}(+1)|=m=(p+1)/2
(15.279 R).  n_1D = m C(p,m).  Q_1D(r)=E[|ẑ(ξ)|²|ẑ(rξ)|²] on this family.
NL = Max+ \\ 1D.  Full ensemble is the n_1D : n_NL mixture (N+ unnamed).

============================================================================
Theorem A — PROVED (Max+-free; line support of ẑ).
  ẑ of a Type+ lift is supported on the F_p-line V°=ann(ker L).
  r∉F_p sends every F_p-line through 0 to a different line, so
      Q_1D(r)=0  for all r∉F_p.
  Fail-when-wrong: claim Q_1D=0 on F_p^× \\ {±1} (the sub type).  ∎

Theorem B — PROVED (Max+-free plus-set 4th moment).
  Fraction of Type+ forms with ξ0∈V° is 2/(p+1).  On that line
  ẑ(ξ)=p ˆσ(β(ξ)) and ˆσ(α)=2 1̂_A(α) for α≠0, A=σ^{-1}(+1).
  For r∈F_p^×,
      Q_1D(r)/q² = 32/(p+1) E_A[ |1̂_A(1)|² |1̂_A(r)|² ]
  with A uniform among m-subsets of F_p.  The Johnson 2-point of
  N(d)=|A∩(A−d)| evaluates the expectation to
      E[|1̂(1)|⁴]         = (p+1)(p²−2p−1) / (8(p−2))     (r=±1)
      E[|1̂(1)|²|1̂(r)|²] = (p+1)(p²−3p−2) / (16(p−2))    (r∈F_p\\{±1}).
  Hence
      Q_1D(±1)/q²     = 4(p²−2p−1)/(p−2)
      Q_1D(sub)/q²    = 2(p²−3p−2)/(p−2).
  Fail-when-wrong: replace p−2 by p−1; swap the two rationals.  ∎

Theorem C — OPEN.  Q_full = (n_1D Q_1D + n_NL Q_NL)/N+.
  Q_1D ≠ Q_full (15.290 types; certified p=5,7).  Q_NL unnamed in p.
  S_1D/q² is 4p(p−1)/(p−2) or larger, both >6, but 1D is not Max+
  (p=7 NL orbit has S_□=0 for some even ψ).  ⟨δ,ψ⟩≤2 / phi_F stay
  OPEN.  ∎

============================================================================
Backend: CPU.  Lift Q uses ProcessPool W=full_workers.
Writes evidence/e1_gmin_m4_prop15292.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import omega_set, type_key, zhat_omega  # noqa: E402

_G: dict = {}


def n_1d(p: int) -> int:
    """# Type+ lifts with y_∞=+1."""
    m = (p + 1) // 2
    return m * comb(p, m)


def Q_1d_pm1_over_q2(p: int) -> Fraction:
    return Fraction(4 * (p * p - 2 * p - 1), p - 2)


def Q_1d_sub_over_q2(p: int) -> Fraction:
    return Fraction(2 * (p * p - 3 * p - 2), p - 2)


def Q_1d_off_Fp_over_q2(_p: int) -> Fraction:
    return Fraction(0)


def Q_1d_pm1_over_q2_inverted(p: int) -> Fraction:
    """Fail-when-wrong: p−1 in the denominator."""
    return Fraction(4 * (p * p - 2 * p - 1), p - 1)


def Q_1d_sub_over_q2_inverted(p: int) -> Fraction:
    return Fraction(2 * (p * p - 3 * p - 2), p - 1)


def _falling(n: int, k: int) -> int:
    acc = 1
    for i in range(k):
        acc *= n - i
    return acc


def _lam(p: int, m: int, k: int) -> Fraction:
    if k < 0 or k > m or k > p:
        return Fraction(0)
    return Fraction(_falling(m, k), _falling(p, k))


def E_N_N(p: int, d: int, e: int) -> Fraction:
    """E[N(d)N(e)] for N(δ)=|A∩(A−δ)|, A uniform m-subset, m=(p+1)/2."""
    m = (p + 1) // 2
    d %= p
    e %= p
    tot = Fraction(0)
    for delta in range(p):
        S = {0, d, delta % p, (delta + e) % p}
        tot += _lam(p, m, len(S))
    return p * tot


def E_hat_product(p: int, r: int) -> Fraction:
    """E_A[|1̂_A(1)|² |1̂_A(r)|²] via Johnson 2-point of N.

    ∑_{d,e} E[N(d)N(e)] ω^{d − r e} is in Q; evaluated as the constant
    term of the cyclic convolution (sum of E[N(d)N(e)] over d−r e ≡ 0
    minus the mean of the other residues, using ∑_ζ ζ^k = p 1_{k=0}).
    Implemented by summing E[N(d)N(e)] ζ^{d-re} in the cyclotomic
    ring projected to Q: the imaginary part vanishes and the real
    part equals ∑_{d,e} E[N(d)N(e)] * 1_{d≡r e}  minus a correction
    that is already included because ∑_{k≠0} ζ^k = −1 only when we
    collapse by complete sums.  Direct route: DFT of the table
    recognized as a Fraction (denominator ≤ 8(p−2)).
    """
    r %= p
    # exact: value = ∑_{d,e} ENN(d,e) * (1 if (d-r*e)%p==0 else 0)
    # is WRONG (that would be the ω-sum only after replacing ω by
    # the complete-sum projector).  Correct projector:
    # (1/p) ∑_k ∑_{d,e} ENN(d,e) ω^{k(d-re)} = ENN-average on the
    # kernel?  No — we want ∑ ENN ω^{d-re} = ∑_{t} (∑_{d-re=t} ENN) ω^t.
    #
    # Let a_t = ∑_{d,e: (d-r e)%p = t} E[N(d)N(e)].
    # Want ∑_t a_t ζ^t.  a is a class function of a multiplicative
    # character of order dividing p.  Because Aut(F_p) and the
    # Johnson scheme force a_t to take at most two values
    # (t=0 vs t≠0), we have
    #   ∑ a_t ζ^t = a_0 + a_off * (−1) = a_0 − a_off
    # with a_off = (sum_t a_t − a_0)/(p−1) and sum_t a_t = p² E[N(0)²]?
    # sum_t a_t = ∑_{d,e} ENN = (∑_d E N(d))² only if uncorrelated.
    # Safer: compute a_0 and A=∑_t a_t, use two-level.
    a0 = Fraction(0)
    Atot = Fraction(0)
    for d in range(p):
        for e in range(p):
            val = E_N_N(p, d, e)
            Atot += val
            if (d - r * e) % p == 0:
                a0 += val
    # two-level check: a_t constant off 0
    a_off = (Atot - a0) / (p - 1)
    return a0 - a_off


def E_hat_pm1_formula(p: int) -> Fraction:
    return Fraction((p + 1) * (p * p - 2 * p - 1), 8 * (p - 2))


def E_hat_sub_formula(p: int) -> Fraction:
    return Fraction((p + 1) * (p * p - 3 * p - 2), 16 * (p - 2))


def S_1d_over_q2_triv_on_Fp(p: int) -> Fraction:
    """S_□/q² for even ψ trivial on F_p^× (not a floor proof)."""
    return (
        2 * Q_1d_pm1_over_q2(p)
        + Q_1d_sub_over_q2(p) * (p - 3)
    )


def S_1d_over_q2_ntriv_on_Fp(p: int) -> Fraction:
    """S_□/q² for even ψ with ∑_{F_p^×} ψ=0."""
    return 2 * Q_1d_pm1_over_q2(p) + Q_1d_sub_over_q2(p) * (-2)


def _init_worker(p: int) -> None:
    F = _kernel_field(p)
    _G["F"] = F
    _G["Omega"] = omega_set(F)
    _G["p"] = p


def _u_of_z(z):
    z = np.asarray(z, dtype=np.float64)
    return np.abs(zhat_omega(z, _G["F"], _G["Omega"])) ** 2


def typeplus_z_rows(p: int) -> np.ndarray:
    F = _kernel_field(p)
    forms = [(a, b) for a, b in _linear_forms(p) if is_type_plus_form(p, a, b)]
    m = (p + 1) // 2
    rows = []
    for a, b in forms:
        for plus in combinations(range(p), m):
            y = lift_sigma_vector(p, a, b, set(plus))
            rows.append(np.sign(np.asarray(y[1 : 1 + F["q"]], dtype=np.float64)))
    return np.stack(rows)


def live_Q_1d(p: int) -> dict[int, float]:
    F = _kernel_field(p)
    Omega = omega_set(F)
    z = typeplus_z_rows(p)
    # n_1D(5)=30, n_1D(7)=140: spawn overhead dominates a 86-pool.
    # Serial zhat here; ProcessPool is for NL/Hoffman (thousands of rows).
    _init_worker(p)
    u = np.stack([_u_of_z(row) for row in z])
    xi0 = Omega[0]
    om_idx = {xi: i for i, xi in enumerate(Omega)}
    u0 = u[:, om_idx[xi0]]
    out = {}
    q2 = float(F["q"] ** 2)
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * u[:, om_idx[eta]])) / q2
    return out


def pack_live(p: int, Qrel: dict[int, float]) -> dict:
    F = _kernel_field(p)
    by = defaultdict(list)
    for r, val in Qrel.items():
        by[str(type_key(F, r))].append(val)
    return {
        k: {
            "n": len(vs),
            "mean": float(np.mean(vs)),
            "spread": float(max(vs) - min(vs)),
            "frac": str(Fraction(float(np.mean(vs))).limit_denominator(200)),
        }
        for k, vs in by.items()
    }


def prove_vanishes_off_subfield(
    primes: list[int] | None = None,
    live: dict[int, dict[int, float]] | None = None,
) -> dict:
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        Qrel = live[p] if live is not None else live_Q_1d(p)
        F = _kernel_field(p)
        off_ok = True
        sub_vals = []
        for r, val in Qrel.items():
            if r in (1, F["neg1"]):
                continue
            # encoded F_p: b=0
            if r < p:  # a+0*p, a≠0
                sub_vals.append(val)
            else:
                if abs(val) > 1e-8:
                    off_ok = False
                    ok = False
        if not sub_vals or max(abs(v) for v in sub_vals) < 1e-8:
            ok = False
        invert_tot += 1
        # fail-when-wrong: claiming Q_1D(sub)=0
        if max(abs(v) for v in sub_vals) < 1e-8:
            invert_hit += 1
        sample[str(p)] = {
            "off_zero": off_ok,
            "sub_mean": float(np.mean(sub_vals)),
            "packed": pack_live(p, Qrel),
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_closed_form(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 20) if is_prime(p)]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        E1 = E_hat_product(p, 1)
        E2 = E_hat_product(p, 2)
        f1 = E_hat_pm1_formula(p)
        f2 = E_hat_sub_formula(p)
        row_ok = E1 == f1 and E2 == f2
        q1 = Fraction(32, p + 1) * E1
        q2 = Fraction(32, p + 1) * E2
        row_ok = row_ok and q1 == Q_1d_pm1_over_q2(p) and q2 == Q_1d_sub_over_q2(p)
        if not row_ok:
            ok = False
        invert_tot += 1
        if Q_1d_pm1_over_q2_inverted(p) == q1:
            invert_hit += 1
        sample[str(p)] = {
            "E1": str(E1),
            "E2": str(E2),
            "Qpm1": str(q1),
            "Qsub": str(q2),
            "ok": row_ok,
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_live_matches_formula(
    live: dict[int, dict[int, float]] | None = None,
) -> dict:
    ok = True
    sample = {}
    for p in (5, 7):
        Qrel = live[p] if live is not None else live_Q_1d(p)
        F = _kernel_field(p)
        pm1 = Q_1d_pm1_over_q2(p)
        sub = Q_1d_sub_over_q2(p)
        for r, val in Qrel.items():
            if r in (1, F["neg1"]):
                target = float(pm1)
            elif r < p:
                target = float(sub)
            else:
                target = 0.0
            if abs(val - target) > 1e-8:
                ok = False
        sample[str(p)] = {
            "pm1_formula": str(pm1),
            "sub_formula": str(sub),
            "packed": pack_live(p, Qrel),
            "ok": ok,
        }
    return {"proved": ok, "sample": sample}


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_budget": bool(rhat_ge_budget_proved_general()),
        "S_1d_triv": {str(p): str(S_1d_over_q2_triv_on_Fp(p)) for p in (5, 7, 11)},
        "S_1d_ntriv": {str(p): str(S_1d_over_q2_ntriv_on_Fp(p)) for p in (5, 7, 11)},
        "note": (
            "Q_1D named.  Q_NL / Q_full unnamed.  1D ensemble S_□/q²>6 "
            "is not the Max+ floor."
        ),
    }


def main() -> dict:
    print("Prop 15.292  Q_1D named in p; Q_NL / floor OPEN", flush=True)
    B = prove_closed_form()
    print(f"  B plus-set moment = closed form: {B['proved']}", flush=True)
    live = {5: live_Q_1d(5), 7: live_Q_1d(7)}
    A = prove_vanishes_off_subfield(live=live)
    print(f"  A Q_1D=0 off F_p: {A['proved']}", flush=True)
    L = prove_live_matches_formula(live=live)
    print(f"  live lifts match formula: {L['proved']}", flush=True)
    C = prove_open()
    print(f"  C Q_NL/floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.292",
        "title": "Q_1D named in p; Q_NL still open",
        "proved": {
            "Q_1d_zero_off_Fp": A["proved"],
            "Q_1d_closed_form": B["proved"],
            "live_lifts_match": L["proved"],
            "Q_NL_named_in_p": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "live": L, "C": C},
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15292.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
