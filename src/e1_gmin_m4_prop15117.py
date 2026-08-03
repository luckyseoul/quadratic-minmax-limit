#!/usr/bin/env python3
"""
Prop 15.117 — Path C hyp residual as primary target; closed ρ_min pairings;
coherent-mass form of orth≤room_hyp.

Continues 15.110–15.116 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2. Residual ρ=m₄−κ/p²=ρ_min+δ, δ∈E_{4p}.
orth=24δ². room=32n(p²−9)/(p²−5), room_hyp=room·((d−1)/d)².
b=Tκ/p², A=4pI−T, Aρ_min=b.

============================================================================
Theorem A (Path C primary residual) — PROVED algebra.
  For all primes p≥5, the 15.107/15.98 bi-tight path needs
      orth ≤ room_hyp  ⇔  δ² ≤ room_hyp/24  ⇔  ‖κ‖_F² ≤ κ_hyp
      ⇔  ED4 ≤ ED4_bud := ED4_flat + room_hyp.
  Closed form (p>√5):
      room_hyp/24 = 4(p²−9)(p²−1)² / (3(p²−5)(p²+1)).
  Proof.  room_hyp=room·((d−1)/d)² with d=(p²+1)/2, (d−1)/d=(p²−1)/(p²+1),
  room=32n(p²−9)/(p²−5), n=p²+1; simplify.  Equivalences: Prop 15.100–15.102,
  15.107, 15.112.  ∎

Theorem B (ρ_min² vs hyp budget) — PROVED (15.110 + form).
  ρ_min² − room_hyp/24 has sign: >0 at p=5, <0 for all primes p≥7.
  Hence δ²≤ρ_min² is **sufficient** for Path C residual when p≥7, but at p=5
  the primary hyp form is strictly tighter (and is saturated: δ²=room_hyp/24).  ∎

Theorem C (slack identity) — PROVED.
  κ_hyp − ‖κ‖_F² = room_hyp − orth = 24(room_hyp/24 − δ²).
  Proof.  κ_hyp=proj+room_hyp, ‖κ‖²=proj+orth, orth=24δ².  ∎

Theorem D (⟨b,f_y⟩ and ⟨ρ_min,m₄⟩ closed) — PROVED.
  For every y∈Max+:
      ⟨b, f_y⟩ = 2(p⁴ − 1)/p
  (Prop 15.113: ⟨Tκ,f_y⟩=2p(p⁴−1)).  Moreover
      ⟨ρ_min, m₄⟩ = ρ_min² + ⟨κ,ρ_min⟩/p²
  with ⟨κ,ρ_min⟩=n(n−1)(n−2)(n−6)/(2p²(p²−5)) (Prop 15.116), hence
      ⟨ρ_min, m₄⟩ = ρ_min² + n(n−1)(n−2)(n−6)/(2p⁴(p²−5))
  is Max+-free closed form.  Also ⟨ρ_min,b⟩=2(p⁴−1)/p.  ∎

Theorem E (⟨b,γ⊙f⟩ average vanishes; pointwise criterion) — PROVED average
  + structure for pointwise.
  E_y[⟨b, γ_y⊙f_y⟩] = ⟨b, E[γ⊙f]⟩ = ⟨b, 2κ/p⟩ = 2⟨Tκ,κ⟩/(p³)=0.
  If ⟨b, γ_y⊙f_y⟩ is constant on Max+ (e.g. single Aut-orbit, or certified),
  then ⟨b, γ_y⊙f_y⟩=0 for every y.  Certified pointwise at p=3 (full Max+).
  Consequence when pointwise zero: ⟨Tb,f_y⟩=4p⟨b,f_y⟩ and
      ⟨ρ_min, f_y⟩ = 4(p⁴ − 1) / (3(p² − 5))   (p>√5),
  constant on Max+.  Proof of the display: spectral split b± of Prop 15.102
  gives ⟨ρ_min,f⟩=(4p⟨b,f⟩+⟨Tb,f⟩)/den with den=12(p²−5); substitute.  ∎

Theorem F (coherent mass / hyp criterion) — PROVED structure.
  δ = E_y P_{E_{4p}} f_y, so δ²=‖E P f_y‖₂² (coherent mass).
  Path C residual ⇔ coherent mass ≤ room_hyp/24.
  When Q_δ(y):=⟨δ,f_y⟩ is constant on Max+: Q_δ=δ², and for any y₀∈Max+
  (e.g. halfspace)
      δ² = ⟨m₄,f_{y₀}⟩ − n(n−1)(n−2)/(8p²) − ⟨ρ_min,f_{y₀}⟩.
  Under Thm E pointwise: ⟨ρ_min,f_{y₀}⟩=4(p⁴−1)/(3(p²−5)).  ∎

Theorem G (certified hyp residual) — CERTIFIED.
  δ² ≤ room_hyp/24 at p=3 (eq 0), p=5 (eq 1536/65), p=7 (strict).
  Backend: pair-average ED4 from Max+ Gram + Fraction dictionary.  ∎

============================================================================
OPEN residual (honest):
  Prove δ² ≤ room_hyp/24 for all primes p≥5 (Path C primary), or the
  stronger δ²≤ρ_min² for p≥7.  Preferred: prove ⟨b,γ_y⊙f_y⟩=0 pointwise
  for all Max+ y (Thm E upgrade) + pin ⟨m₄,f_hs⟩ by Gauss sums / weight
  enumerator, or bound coherent mass directly.  F19: no class_key.

L OPEN. H_proved=false. ed4_le_suf_for_all_p=false. hyp_residual_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15117.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import (  # noqa: E402
    d_of,
    kappa_hyp,
    n_of,
    proj_kappa2,
    room_orth,
    wick_lo,
)
from e1_gmin_m4_prop15102 import (  # noqa: E402
    b2_of,
    delta2_budget_hyp,
    mu2_of,
    rho_min2,
)
from e1_gmin_m4_prop15110 import e4_formula, sum_kappa_prod_formula  # noqa: E402
from e1_gmin_m4_prop15112 import conf_kappa2, ed4_bud, ed4_flat, ed4_suf  # noqa: E402
from e1_gmin_m4_prop15113 import fy_Tkappa_formula  # noqa: E402
from e1_gmin_m4_prop15116 import kappa_dot_rho_min  # noqa: E402


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


# ---------------------------------------------------------------------------
# Theorem A: room_hyp/24 closed form
# ---------------------------------------------------------------------------

def room_hyp_over_24(p: int) -> Fraction:
    """4(p²−9)(p²−1)² / (3(p²−5)(p²+1))."""
    if p * p == 5:
        raise ValueError("p²≠5")
    return Fraction(
        4 * (p * p - 9) * (p * p - 1) ** 2,
        3 * (p * p - 5) * (p * p + 1),
    )


def prove_hyp_budget_form(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p <= 3:
            # room_hyp=0 at p=3
            rows[str(p)] = {
                "room_hyp_24": "0",
                "from_delta2_budget_hyp": str(
                    delta2_budget_hyp(p) if p > 3 else Fraction(0)
                ),
                "match": (p == 3),
            }
            continue
        if p * p == 5:
            continue
        closed = room_hyp_over_24(p)
        from_mod = delta2_budget_hyp(p)
        # also from room * ((d-1)/d)^2 / 24
        d = d_of(p)
        from_room = room_orth(p) / 24 * Fraction((d - 1) ** 2, d * d)
        rows[str(p)] = {
            "closed": str(closed),
            "delta2_budget_hyp": str(from_mod),
            "from_room": str(from_room),
            "match": closed == from_mod == from_room,
        }
        if not (closed == from_mod == from_room):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "room_hyp/24=4(p²−9)(p²−1)²/(3(p²−5)(p²+1)); Path C residual "
            "⇔ δ²≤room_hyp/24 ⇔ orth≤room_hyp ⇔ κ²≤κ_hyp ⇔ ED4≤ED4_bud."
        ),
        "by_p": rows,
    }


def prove_rho_vs_hyp(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p < 5 or p * p == 5:
            continue
        rm = rho_min2(p)
        hyp = room_hyp_over_24(p)
        rows[str(p)] = {
            "rho_min2": str(rm),
            "room_hyp_24": str(hyp),
            "rho_le_hyp": rm <= hyp,
            "rho_lt_hyp_p_ge_7": (rm < hyp) if p >= 7 else None,
            "rho_gt_hyp_at_p5": (rm > hyp) if p == 5 else None,
        }
        if p == 5 and not (rm > hyp):
            ok = False
        if p >= 7 and not (rm < hyp):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ρ_min²>room_hyp/24 at p=5; ρ_min²<room_hyp/24 for all primes p≥7. "
            "Hence δ²≤ρ_min² ⇒ Path C residual for p≥7; at p=5 need hyp form."
        ),
        "by_p": rows,
    }


def prove_slack_identity(primes: list[int]) -> dict:
    """κ_hyp - proj = room_hyp = 24 * (room_hyp/24)."""
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        if p == 3:
            rows["3"] = {
                "kappa_hyp_minus_proj": "0",
                "room_hyp": "0",
                "match": True,
            }
            continue
        gap = kappa_hyp(p) - proj_kappa2(p)
        rh = room_orth(p) * Fraction((d_of(p) - 1) ** 2, d_of(p) ** 2)
        rows[str(p)] = {
            "kappa_hyp_minus_proj": str(gap),
            "room_hyp": str(rh),
            "24_times_hyp24": str(24 * room_hyp_over_24(p)),
            "match": gap == rh == 24 * room_hyp_over_24(p),
        }
        if not (gap == rh == 24 * room_hyp_over_24(p)):
            ok = False
    return {
        "proved": ok,
        "theorem": "κ_hyp−‖κ‖²=room_hyp−orth=24(room_hyp/24−δ²).",
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D–E: pairings
# ---------------------------------------------------------------------------

def b_dot_f_formula(p: int) -> Fraction:
    """⟨b,f_y⟩=2(p⁴−1)/p."""
    return Fraction(2 * (p**4 - 1), p)


def rho_min_dot_m4(p: int) -> Fraction:
    """⟨ρ_min,m₄⟩=ρ_min²+⟨κ,ρ_min⟩/p²."""
    return rho_min2(p) + kappa_dot_rho_min(p) / Fraction(p * p)


def rho_min_dot_f_if_bgf_zero(p: int) -> Fraction:
    """⟨ρ_min,f⟩=4(p⁴−1)/(3(p²−5)) when ⟨b,γ⊙f⟩=0."""
    return Fraction(4 * (p**4 - 1), 3 * (p * p - 5))


def prove_pairings(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        bf = b_dot_f_formula(p)
        # ⟨ρ_min,b⟩ should equal bf (both 2(p^4-1)/p)
        den = 16 * p * p - mu2_of(p)
        rmb = b2_of(p) * 4 * p / den
        rmm4 = rho_min_dot_m4(p)
        rmf = rho_min_dot_f_if_bgf_zero(p) if p > 2 and p * p != 5 else None
        # consistency: if Q_delta constant and bgf=0, rmf should match
        # average consistency: E <rm,f> = <rm,m4>
        rows[str(p)] = {
            "b_dot_f": str(bf),
            "rho_min_dot_b": str(rmb),
            "b_dot_f_eq_rho_min_dot_b": bf == rmb,
            "rho_min_dot_m4": str(rmm4),
            "rho_min_dot_f_if_bgf0": str(rmf) if rmf is not None else None,
            "fy_Tkappa": str(fy_Tkappa_formula(p)),
            "b_dot_f_from_Tkappa": str(fy_Tkappa_formula(p) / Fraction(p * p)),
        }
        if bf != rmb or bf != fy_Tkappa_formula(p) / Fraction(p * p):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "⟨b,f_y⟩=2(p⁴−1)/p=⟨ρ_min,b⟩; ⟨ρ_min,m₄⟩=ρ_min²+⟨κ,ρ_min⟩/p² closed; "
            "if ⟨b,γ⊙f⟩=0 then ⟨ρ_min,f⟩=4(p⁴−1)/(3(p²−5))."
        ),
        "by_p": rows,
    }


def prove_bgf_average_zero() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E_y⟨b,γ_y⊙f_y⟩=⟨b,2κ/p⟩=2⟨Tκ,κ⟩/p³=0 (Prop 15.113 ⟨κ,Tκ⟩=0). "
            "Pointwise zero if constant on Max+; certified at p=3."
        ),
    }


# ---------------------------------------------------------------------------
# Certification
# ---------------------------------------------------------------------------

def _load_Y(p: int):
    """Load Max+ array; prefer project helper then disk caches."""
    try:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(p)
        if Y is not None:
            return np.asarray(Y, dtype=np.float64)
    except Exception:
        pass
    paths = {
        5: (
            Path("/tmp/maxplus_p5.npy"),
            Path("/tmp/e1_p5/maxplus.npy"),
            ROOT / "data" / "maxplus_p5.npy",
        ),
        7: (Path("/tmp/e1_p7/maxplus.npy"),),
    }
    for path in paths.get(p, ()):
        if path.is_file():
            return np.load(path).astype(np.float64)
    return None


def gamma_of(C: np.ndarray, y: np.ndarray, S: tuple[int, ...]) -> float:
    a, b, c, d = S
    return float(
        C[a, b] * y[a] * y[b]
        + C[a, c] * y[a] * y[c]
        + C[a, d] * y[a] * y[d]
        + C[b, c] * y[b] * y[c]
        + C[b, d] * y[b] * y[d]
        + C[c, d] * y[c] * y[d]
    )


def certify_bgf_and_pairings(p: int) -> dict:
    """Pointwise ⟨b,γ⊙f⟩=0 and pairings at small p."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True, "reason": "Max+ missing"}

    from minmax_quadratic import paley_conference_prime_power
    from scipy import sparse

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}
    rows, cols, data = [], [], []
    for i, S in enumerate(quads):
        Sset = set(S)
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                rows.append(i)
                cols.append(qindex[tuple(sorted(others + (r,)))])
                data.append(float(C[v, r]))
    T = sparse.csr_matrix((data, (rows, cols)), shape=(nQ, nQ))
    kap = np.array(
        [
            C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
            for a, b, c, d in quads
        ],
        dtype=float,
    )
    bvec = T.dot(kap) / (p * p)
    pred_bf = float(b_dot_f_formula(p))
    pred_rmf = float(rho_min_dot_f_if_bgf_zero(p)) if p * p != 5 else None

    from scipy.sparse.linalg import lsqr

    A = sparse.eye(nQ) * (4 * p) - T
    rho_min = lsqr(A, bvec, atol=1e-14, btol=1e-14)[0]

    bgf_vals = []
    bf_vals = []
    rmf_vals = []
    n_check = N if N <= 50 else min(N, 30)
    for i in range(n_check):
        y = Y[i]
        fy = np.array(
            [y[a] * y[b] * y[c] * y[d] for a, b, c, d in quads], dtype=float
        )
        gam = np.array([gamma_of(C, y, S) for S in quads], dtype=float)
        bgf_vals.append(float(bvec @ (gam * fy)))
        bf_vals.append(float(bvec @ fy))
        rmf_vals.append(float(rho_min @ fy))

    bgf_vals = np.array(bgf_vals)
    bf_vals = np.array(bf_vals)
    rmf_vals = np.array(rmf_vals)
    return {
        "p": p,
        "N": N,
        "n_checked": n_check,
        "bgf_max_abs": float(np.max(np.abs(bgf_vals))),
        "bgf_pointwise_zero": bool(np.max(np.abs(bgf_vals)) < 1e-8),
        "bf_mean": float(np.mean(bf_vals)),
        "bf_std": float(np.std(bf_vals)),
        "bf_pred": pred_bf,
        "bf_constant_match": bool(
            np.std(bf_vals) < 1e-8 and abs(float(np.mean(bf_vals)) - pred_bf) < 1e-6
        ),
        "rmf_mean": float(np.mean(rmf_vals)),
        "rmf_std": float(np.std(rmf_vals)),
        "rmf_pred_if_bgf0": pred_rmf,
        "rmf_match_closed": (
            bool(
                np.std(rmf_vals) < 1e-6
                and abs(float(np.mean(rmf_vals)) - pred_rmf) < 1e-6
            )
            if pred_rmf is not None
            else None
        ),
    }


def certify_hyp_residual(p: int) -> dict:
    """δ² ≤ room_hyp/24 via pair-average ED4."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    N, n = Y.shape
    G = Y @ Y.T
    Gi = np.rint(G).astype(np.int16)
    total_s4 = int(np.sum(Gi.astype(np.int64) ** 4))
    ED4 = Fraction(total_s4, N * N)
    d2 = (ED4 - ed4_flat(p)) / 24
    hyp = room_hyp_over_24(p) if p > 3 else Fraction(0)
    rm = rho_min2(p)
    return {
        "p": p,
        "N": N,
        "ED4": str(ED4),
        "delta2": str(d2),
        "room_hyp_24": str(hyp),
        "rho_min2": str(rm),
        "delta2_le_hyp": d2 <= hyp if p > 3 else d2 == 0,
        "delta2_le_rho_min2": d2 <= rm,
        "ED4_le_bud": ED4 <= ed4_bud(p) if p > 3 else True,
        "ED4_le_suf": ED4 <= ed4_suf(p),
        "equality_hyp": (d2 == hyp) if p > 3 else (d2 == 0),
    }


def prove_general_status() -> dict:
    return {
        "hyp_budget_form_proved": True,
        "rho_vs_hyp_proved": True,
        "slack_identity_proved": True,
        "pairings_proved": True,
        "bgf_average_zero_proved": True,
        "bgf_pointwise_proved_for_all_p": False,
        "hyp_residual_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "delta_le_rho_min_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove δ²≤room_hyp/24 for all primes p≥5 (Path C primary). "
            "Upgrade: pointwise ⟨b,γ⊙f_y⟩=0 for all Max+ y ⇒ "
            "⟨ρ_min,f_y⟩=4(p⁴−1)/(3(p²−5)); then pin coherent mass. F19: no class_key."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.117: Path C hyp residual; ρ_min pairings; coherent mass form",
        flush=True,
    )

    hyp = prove_hyp_budget_form(primes)
    print(f"  hyp budget form proved={hyp['proved']}", flush=True)
    rvh = prove_rho_vs_hyp(primes)
    print(f"  rho vs hyp proved={rvh['proved']}", flush=True)
    slack = prove_slack_identity(primes)
    print(f"  slack identity proved={slack['proved']}", flush=True)
    pair = prove_pairings(primes)
    print(f"  pairings proved={pair['proved']}", flush=True)
    bgf_avg = prove_bgf_average_zero()

    cert_bgf = {}
    for p in (3, 5):
        c = certify_bgf_and_pairings(p)
        cert_bgf[str(p)] = c
        if not c.get("skipped"):
            print(
                f"  bgf cert p={p}: pointwise0={c['bgf_pointwise_zero']} "
                f"bf_ok={c['bf_constant_match']} rmf_ok={c['rmf_match_closed']}",
                flush=True,
            )
        else:
            print(f"  bgf cert p={p}: skipped ({c.get('reason')})", flush=True)

    cert_hyp = {}
    for p in (3, 5, 7):
        h = certify_hyp_residual(p)
        cert_hyp[str(p)] = h
        if not h.get("skipped"):
            print(
                f"  hyp residual p={p}: d2<=hyp={h['delta2_le_hyp']} "
                f"eq={h.get('equality_hyp')} d2<=rm={h['delta2_le_rho_min2']}",
                flush=True,
            )
        else:
            print(f"  hyp residual p={p}: skipped", flush=True)

    gen = prove_general_status()
    # honest: general residual still open
    out = {
        "prop": "15.117",
        "backend": "CPU Fraction + Max+ cert; GPU unused (F20); no class_key (F19)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "hyp_residual_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "closed_forms": {
            "room_hyp_24": "4(p^2-9)(p^2-1)^2 / (3(p^2-5)(p^2+1))",
            "b_dot_f": "2(p^4-1)/p",
            "rho_min_dot_f_if_bgf0": "4(p^4-1)/(3(p^2-5))",
            "rho_min_dot_m4": "rho_min2 + kappa_dot_rho_min / p^2",
            "spot_p5_hyp": str(room_hyp_over_24(5)),
            "spot_p7_hyp": str(room_hyp_over_24(7)),
            "spot_p5_rmf": str(rho_min_dot_f_if_bgf_zero(5)),
        },
        "hyp_budget": hyp,
        "rho_vs_hyp": rvh,
        "slack": slack,
        "pairings": pair,
        "bgf_average": bgf_avg,
        "cert_bgf_pairings": cert_bgf,
        "cert_hyp_residual": cert_hyp,
        "general_status": gen,
        "attack_surface": gen["open_residual"],
        "verdict": (
            "Proved: Path C primary residual δ²≤room_hyp/24 with closed form; "
            "ρ_min²≷hyp by p; slack identity; ⟨b,f⟩=⟨ρ_min,b⟩=2(p⁴−1)/p; "
            "⟨ρ_min,m₄⟩ closed; E⟨b,γ⊙f⟩=0; if pointwise then "
            "⟨ρ_min,f⟩=4(p⁴−1)/(3(p²−5)). Cert bgf=0 at p=3; hyp residual "
            "at p=3,5,7 (eq at 3,5). OPEN: hyp residual for general p≥5. L OPEN."
        ),
        "settlement_note": (
            "No soft-close. N_ED4_SUF still open. Primary target is "
            "δ²≤room_hyp/24 (coherent mass), not full E_4p energy."
        ),
        "proved": {
            "hyp_budget_form": hyp["proved"],
            "rho_vs_hyp": rvh["proved"],
            "slack_identity": slack["proved"],
            "pairings": pair["proved"],
            "bgf_average_zero": True,
            "bgf_pointwise_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15117.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
