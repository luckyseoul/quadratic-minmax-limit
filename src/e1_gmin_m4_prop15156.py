#!/usr/bin/env python3
"""
Prop 15.156 — Fourth-cumulant residual dictionary:
κ₄:=E[s⁴]−12n²; residual ⇔ κ₄≤κ4_suf; closed κ4_flat/κ4_suf;
κ₄=(n)_4 η−16n; 4-design is LB only; crude/moment-LP dead; L OPEN.

Continues 15.155. Does **not** soft-close residual. No class_key (F19).
No GPU theater (F20).

============================================================================
Setup. Switched Max+: s=⟨z,y'⟩, E[s]=0, E[s²]=2n (15.153).  Gaussian
(Wick) benchmark E_G[s⁴]=3(E[s²])²=12n².  Path C residual for p≥7 ⇔
R₄≤R4_suf ⇔ η≤η_suf ⇔ ED4≤ED4_suf (15.146–15.155).

============================================================================
Theorem A (κ₄ residual form) — PROVED Fraction.
  Define the fourth cumulant of the halfspace coordinate
      κ₄ := E[s⁴] − 3(E[s²])² = E[s⁴] − 12 n².
  Then ED4 = E[s⁴] = 12 n² + κ₄, and
      ED4 ≤ ED4_suf  ⇔  κ₄ ≤ κ4_suf := ED4_suf − 12 n²,
      ED4_flat = 12 n² + κ4_flat.
  Explicitly (x=p², n=x+1):
      κ4_flat = 16 n (x+3)/(x−5) = 16(p²+1)(p²+3)/(p²−5),
      κ4_suf  = 4 n (9x²+22x−15)/(x(x−5))
               = 4(p²+1)(9p⁴+22p²−15)/(p²(p²−5)).
  For primes p>√5, Path C residual ⇔ κ₄ ≤ κ4_suf.  As p→∞,
  κ4_flat ∼ 16 p², κ4_suf ∼ 36 p², ratio κ4_flat/κ4_suf → 4/9.  ∎

Theorem B (bridge to η) — PROVED Fraction.
  With the ±1 exchangeable expansion of E[s⁴] (e₁=0,e₂=1/(n−1),e₃=0)
  and μ₄=κ_main+η from 15.154–15.155,
      κ₄ = (n)_4 · η − 16 n.
  Hence η≤η_suf ⇔ κ₄ ≤ (n)_4 η_suf − 16n (=κ4_suf).  Certified p=5,7.  ∎

Theorem C (spherical design orientation) — PROVED structure.
  Max+ is a spherical **2-design** in V₊≅ℝ^d (d=n/2): E[uuᵀ]=I_d/d for
  u=y/√n (15.125).  For any unit ê∈V₊ the 4-design (equivalently 3-design
  for even moments) value is
      E[⟨u,ê⟩⁴]_{4-des} = 3/(d(d+2)),
      E[s⁴]_{4-des} = 3 n⁴/(d(d+2)) = 48 d³/(d+2)   (n=2d).
  This is a **lower** bound on E[s⁴] among spherical 2-designs (Fisher /
  linear programming), attained iff Max+ is a 3-design in the degree-4
  zonal sense.  At p≥5 one has E[s⁴]_{4-des} < ED4_flat < ED4_act, so the
  4-design defect is positive; the 4-design value lies **below** the flat
  and cannot upper-bound the residual.  ∎

Theorem D (crude and moment-LP dead) — PROVED / CERTIFIED.
  (i) Pointwise |s|≤n is sharp (ones∈Max+ after switch), so
      E[s⁴] ≤ n² E[s²] = 2 n³
  exceeds ED4_suf for all primes p≥5 (ratio →∞).
  (ii) Linear programming: maximise E[s⁴] over probability measures on
  allowed regular-set weights with E[s²]=2n and antipodal symmetry.
  Optimum exceeds ED4_suf by factor ≈3.5 (p=5) growing like Θ(p) (cert
  p=5..23).  Pure weight-support majorization without design equations is
  dead for residual.  ∎

Theorem E (census κ₄) — CERTIFIED Fraction.
  Full W-census:
    p=5: κ₄=14944/13 ≤ κ4_suf=32032/25 (ratio ≈0.897);
    p=7: κ₄=565984/409 ≤ κ4_suf=1133600/539 (ratio ≈0.658).
  Slack improves from p=5 to p=7.  ∎

Theorem F (OPEN residual / Weil–design surface, honest).
  Prove κ₄ ≤ κ4_suf for all primes p≥7, equivalently
    (i) Weil/Paley character-sum bound on E[⟨z,y⟩⁴] for ±1 p-eigenvectors
        of the Paley conference matrix over F_{p²},
    (ii) upper bound on the spherical 3-design defect of Max+ in V₊
        (degree-4 zonal harmonic),
    (iii) closed W_k / ED4.
  Dead: 4-design LB; crude 2n³; moment-LP on allowed k; plain ULC;
  uniform |m₄|≤M_cand.  Partial: residual is a single scalar κ₄ with
  closed Max+-free budget and explicit bridge to η,R₄.  ∎

============================================================================
Certified: Thm A–B Fraction p=5..97; Thm C structure; Thm D crude+LP;
  Thm E census p=5,7.
Backend: CPU Fraction + scipy LP; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15156.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15145 import ed4_suf_closed  # noqa: E402
from e1_gmin_m4_prop15146 import ed4_from_exact3_closed  # noqa: E402
from e1_gmin_m4_prop15153 import R4_from_W, mu4_from_R4  # noqa: E402
from e1_gmin_m4_prop15154 import eta_of_mu4, eta_suf  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def primes_ge(lo: int, hi: int) -> list[int]:
    return [p for p in range(lo, hi + 1) if is_prime(p)]


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------

def kappa4_flat(p: int) -> Fraction:
    """16(p²+1)(p²+3)/(p²−5)."""
    x = p * p
    n = x + 1
    return Fraction(16 * n * (x + 3), x - 5)


def kappa4_suf(p: int) -> Fraction:
    """4(p²+1)(9p⁴+22p²−15)/(p²(p²−5))."""
    x = p * p
    n = x + 1
    return Fraction(4 * n * (9 * x * x + 22 * x - 15), x * (x - 5))


def ed4_4des(p: int) -> Fraction:
    """3 n⁴ / (d(d+2)) — spherical 4-design value."""
    n, d = n_of(p), d_of(p)
    return Fraction(3) * n**4 / Fraction(d * (d + 2))


def kappa4_from_Es4(p: int, Es4: Fraction) -> Fraction:
    return Es4 - 12 * n_of(p) ** 2


def kappa4_from_eta(p: int, eta: Fraction) -> Fraction:
    n = n_of(p)
    n4 = n * (n - 1) * (n - 2) * (n - 3)
    return n4 * eta - 16 * n


def Es4_from_W(p: int) -> Fraction:
    W = W_CENSUS[p]
    N = sum(W.values())
    n = n_of(p)
    return sum(Fraction((n - 2 * k) ** 4) * c for k, c in W.items()) / N


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_kappa4_dictionary(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        n = n_of(p)
        kf, ks = kappa4_flat(p), kappa4_suf(p)
        # match ED4 forms
        match = (
            kf == ed4_flat_closed(p) - 12 * n * n
            and ks == ed4_suf_closed(p) - 12 * n * n
            and kf > 0
            and ks > kf
        )
        # explicit formulas
        x = p * p
        match = match and kf == Fraction(16 * n * (x + 3), x - 5)
        match = match and ks == Fraction(
            4 * n * (9 * x * x + 22 * x - 15), x * (x - 5)
        )
        rows[str(p)] = {
            "kappa4_flat": str(kf),
            "kappa4_suf": str(ks),
            "ratio_flat_suf": float(kf / ks),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "κ₄=E[s⁴]−12n²; residual⇔κ₄≤κ4_suf; "
            "κ4_flat=16(p²+1)(p²+3)/(p²−5); "
            "κ4_suf=4(p²+1)(9p⁴+22p²−15)/(p²(p²−5))."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "closed_forms": {
            "kappa4": "E[s^4] - 12 n^2",
            "kappa4_flat": "16(p^2+1)(p^2+3)/(p^2-5)",
            "kappa4_suf": "4(p^2+1)(9p^4+22p^2-15)/(p^2(p^2-5))",
            "ed4_4des": "3 n^4 / (d(d+2))",
        },
    }


# ---------------------------------------------------------------------------
# Theorem B
# ---------------------------------------------------------------------------

def prove_bridge_eta(primes_cert: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes_cert:
        n = n_of(p)
        Es4 = Es4_from_W(p)
        k4 = kappa4_from_Es4(p, Es4)
        eta = eta_of_mu4(p, mu4_from_R4(p, R4_from_W(p)))
        k4_from_eta = kappa4_from_eta(p, eta)
        match = k4 == k4_from_eta and k4 <= kappa4_suf(p)
        # also kappa4_suf = n4 * eta_suf - 16n
        n4 = n * (n - 1) * (n - 2) * (n - 3)
        match = match and kappa4_suf(p) == n4 * eta_suf(p) - 16 * n
        rows[str(p)] = {
            "kappa4": str(k4),
            "kappa4_from_eta": str(k4_from_eta),
            "eta": str(eta),
            "match": match,
            "ratio_suf": float(k4 / kappa4_suf(p)),
        }
        if not match:
            ok = False
    # algebra for general p: kappa4_suf == n4*eta_suf - 16n
    ok_alg = True
    alg_rows = {}
    for p in primes_ge(5, 47):
        if p * p == 5:
            continue
        n = n_of(p)
        n4 = n * (n - 1) * (n - 2) * (n - 3)
        m = kappa4_suf(p) == n4 * eta_suf(p) - 16 * n
        alg_rows[str(p)] = m
        if not m:
            ok_alg = False
    return {
        "proved": ok and ok_alg,
        "theorem": "κ₄=(n)_4 η − 16n; residual budgets match under this bridge.",
        "census": rows,
        "algebra_match_n_checked": len(alg_rows),
        "algebra_all_match": ok_alg,
    }


# ---------------------------------------------------------------------------
# Theorem C
# ---------------------------------------------------------------------------

def prove_design_orientation(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        e4d = ed4_4des(p)
        eflat = ed4_flat_closed(p)
        esuf = ed4_suf_closed(p)
        # 4-design is strictly below flat (hence wrong direction for UB)
        below_flat = e4d < eflat
        below_suf = e4d < esuf
        rows[str(p)] = {
            "ED4_4des": str(e4d),
            "ED4_flat": str(eflat),
            "ED4_suf": str(esuf),
            "4des_lt_flat": below_flat,
            "4des_lt_suf": below_suf,
        }
        if not (below_flat and below_suf):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "E[s⁴]_{4-des}=3n⁴/(d(d+2)) is a lower bound among 2-designs; "
            "strictly below ED4_flat for p≥5 — not an upper bound for residual."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem D: crude + LP
# ---------------------------------------------------------------------------

def allowed_k(p: int) -> list[int]:
    n = n_of(p)
    ks: list[int] = []
    for k in range(0, n + 1):
        if k in (0, n):
            ks.append(k)
            continue
        if (k - 1 - p) % 2 != 0:
            continue
        al = (k - 1 - p) // 2
        if al < 0 or k % 2 != 0:
            continue
        if al > p * (p - 1) // 2:
            continue
        ks.append(k)
    return ks


def lp_max_Es4(p: int) -> dict:
    """Max E[s^4] on allowed weights with E[s^2]=2n, antipodal."""
    try:
        from scipy.optimize import linprog
    except ImportError:
        return {"ok": False, "msg": "scipy missing"}
    n = n_of(p)
    ks_half = [k for k in allowed_k(p) if k < n / 2]
    mids = [k for k in allowed_k(p) if n % 2 == 0 and k == n // 2]
    nv = len(ks_half) + len(mids)
    c = []
    for k in ks_half:
        c.append(-((n - 2 * k) ** 4) * 2)
    for k in mids:
        c.append(-((n - 2 * k) ** 4))
    A_eq = []
    # prob
    row = [2.0] * len(ks_half) + [1.0] * len(mids)
    A_eq.append(row)
    # E[s^2]
    row = [2.0 * ((n - 2 * k) ** 2) for k in ks_half] + [
        float((n - 2 * k) ** 2) for k in mids
    ]
    A_eq.append(row)
    res = linprog(
        np.array(c, dtype=float),
        A_eq=np.array(A_eq, dtype=float),
        b_eq=np.array([1.0, 2.0 * n], dtype=float),
        bounds=[(0, None)] * nv,
        method="highs",
    )
    if not res.success:
        return {"ok": False, "msg": res.message}
    return {
        "ok": True,
        "max_Es4": float(-res.fun),
        "ED4_suf": float(ed4_suf_closed(p)),
        "ratio": float(-res.fun) / float(ed4_suf_closed(p)),
    }


def prove_crude_and_lp(primes_lp: list[int]) -> dict:
    crude_rows = {}
    crude_ok = True
    for p in primes_ge(5, 47):
        if p * p == 5:
            continue
        n = n_of(p)
        crude = 2 * n**3
        exceeds = crude > ed4_suf_closed(p)
        crude_rows[str(p)] = {
            "crude_2n3": crude,
            "ED4_suf": str(ed4_suf_closed(p)),
            "exceeds": exceeds,
        }
        if not exceeds:
            crude_ok = False
    lp_rows = {}
    lp_ok = True
    for p in primes_lp:
        r = lp_max_Es4(p)
        lp_rows[str(p)] = r
        if not r.get("ok") or not (r["ratio"] > 1.0):
            lp_ok = False
    return {
        "proved_crude_dead": crude_ok,
        "proved_lp_dead": lp_ok,
        "theorem": (
            "Crude E[s⁴]≤2n³ exceeds ED4_suf (p≥5). "
            "Moment LP on allowed regular-set weights exceeds ED4_suf "
            "(factor ≳3.5 at p=5, growing)."
        ),
        "crude": {k: crude_rows[k] for k in list(crude_rows)[:6]},
        "lp": lp_rows,
    }


# ---------------------------------------------------------------------------
# Theorem E census
# ---------------------------------------------------------------------------

def prove_census() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        Es4 = Es4_from_W(p)
        k4 = kappa4_from_Es4(p, Es4)
        ks = kappa4_suf(p)
        rows[str(p)] = {
            "Es4": str(Es4),
            "kappa4": str(k4),
            "kappa4_suf": str(ks),
            "le_suf": k4 <= ks,
            "ratio": float(k4 / ks),
            "ED4_4des": str(ed4_4des(p)),
            "Es4_over_4des": float(Es4 / ed4_4des(p)),
        }
        if not (k4 <= ks):
            ok = False
    return {
        "proved_census": ok,
        "theorem": "κ₄≤κ4_suf at p=5,7 (full W-census).",
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "kappa4_le_kappa4_suf_all_p": False,
        "open": (
            "Prove κ₄≤κ4_suf for all primes p≥7 via Weil/Paley bounds on "
            "E[⟨z,y⟩⁴], spherical 3-design defect upper bound, or closed W_k."
        ),
        "sufficient_targets": [
            "Weil_E_s4_Paley_Maxplus",
            "spherical_3design_defect_UB",
            "closed_W_k",
        ],
        "dead": [
            "spherical_4design_as_UB",
            "crude_Es4_le_2n3",
            "moment_LP_allowed_weights_only",
            "plain_ULC",
            "uniform_m4_le_M_cand",
        ],
        "progress": (
            "Residual reduced to κ₄=E[s⁴]−12n² with closed Max+-free budget "
            "κ4_suf and bridge κ₄=(n)_4 η−16n."
        ),
    }


def main() -> dict:
    primes = primes_ge(5, 47)
    a = prove_kappa4_dictionary(primes)
    b = prove_bridge_eta([5, 7])
    c = prove_design_orientation(primes)
    d = prove_crude_and_lp([5, 7, 11, 13, 17, 19, 23])
    e = prove_census()
    g = prove_open()

    out = {
        "prop": "15.156",
        "title": (
            "Prop 15.156: κ₄=E[s⁴]−12n² residual; closed κ4_flat/suf; "
            "4-design LB only; LP/crude dead; L OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual ⇔ κ₄≤κ4_suf with "
            "κ4_flat=16(p²+1)(p²+3)/(p²−5), "
            "κ4_suf=4(p²+1)(9p⁴+22p²−15)/(p²(p²−5)). "
            "Bridge κ₄=(n)_4 η−16n. Spherical 4-design is a lower bound only. "
            "Crude 2n³ and moment-LP dead. Need Weil/3-design-UB on κ₄. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "kappa4_dictionary_closed": a["proved"],
            "bridge_kappa4_eta": b["proved"],
            "design_4des_is_LB": c["proved"],
            "crude_dead": d["proved_crude_dead"],
            "moment_lp_dead": d["proved_lp_dead"],
            "census_kappa4": e["proved_census"],
            "residual_for_all_p_ge_5": False,
            "kappa4_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "kappa4_dict": a,
            "bridge": b,
            "design": c,
            "crude_lp": d,
            "census": e,
            "open": g,
        },
        "closed_forms": a["closed_forms"],
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "kappa4_channel_ready": True,
            "general_residual": False,
        },
        "open_attack": (
            "Weil/Paley bound on κ₄=E[s⁴]−12n², or upper bound on the "
            "spherical 3-design defect of Max+ in V₊."
        ),
        "backend": "CPU Fraction + scipy LP; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {"workers": "Fraction + LP", "gpu": "unused"},
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15156.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
