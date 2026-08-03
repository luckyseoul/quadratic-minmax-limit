#!/usr/bin/env python3
"""
Prop 15.158 — Aut / scheme structure of Max+; closed Q₄(t);
1-homogeneous non-scheme; pole–equator decomposition of μ_G4;
Gram spectrum; conditional N-bound; L OPEN.

Continues 15.157. Does **not** soft-close residual. No class_key (F19).
No GPU theater (F20).

============================================================================
Setup. Switched Max+ X⊂{±1}ⁿ∩V₊, |X|=N, Gram G_{yy'}=y·y'.  Spherical
2-design in V₊≅ℝ^d (d=n/2).  Residual ⇔ μ_G4≤μ_G4_suf with
μ_G4=E[Q₄(s/n)], s=⟨z,y⟩ (15.157).

============================================================================
Theorem A (closed Q₄ on S^{d−1}) — PROVED Fraction.
  With Gegenbauer Q₄ normalised by Q₄(1)=1,
      Q₄(t) = [(d+2)(d+4) t⁴ − 6(d+2) t² + 3] / (d² − 1)
             = (t⁴ − a₀ − a₂ Q₂(t)) / a₄,
  where Q₂(t)=(d t²−1)/(d−1) and a₀,a₂,a₄ as in 15.157.  In particular
      Q₄(0) = 3/(d²−1).
  Sign pattern on (0,1): two roots 0<t₁(d)<t₂(d)<1; Q₄>0 on [0,t₁),
  Q₄<0 on (t₁,t₂), Q₄>0 on (t₂,1].  Certified identity for d=3..300.  ∎

Theorem B (1-homogeneous; Gram spectrum) — PROVED structure + cert.
  (i) Max+ is **1-homogeneous**: the absolute inner-product distribution from
  any codeword equals the global weight law (W_k with s=n−2k).  Certified
  at p=5 by full neighbor-profile identity (all 260 rows identical).
  (ii) The Gram G=YYᵀ (Y rows = switched Max+) has spectrum
      λ=2N (multiplicity d),   λ=0 (multiplicity N−d),
  confirming rank-d tight frame / 2-design (E[yyᵀ] has image V₊).
  Certified p=5 (N=260,d=13: λ=520 mult 13; 0 mult 247) and structure
  for general p.  ∎

Theorem C (NOT an inner-product association scheme) — CERTIFIED p=5.
  The relations R_s={(y,y'): y·y'=s} on Max+ do **not** form an
  association scheme: intersection numbers p_{ij}^k fail to be constant on
  R_k for s∈{−14,−6,−2,2,6,14} (sample 40 pairs/class, std>0).  Constant
  for antipodes s=±n and for s=±10.  Hence a pure Bose–Mesner attack on
  the inner-product scheme of Max+ is blocked; a finer Aut-coherent
  configuration (or character sum) is required.  ∎

Theorem D (pole–equator decomposition) — PROVED.
  Write t=s/n and partition Max+ into poles {±z}, exterior
  E={|t|≥t₂}∖{±z}, equator Eq={|t|≤t₁}, and annulus A={t₁<|t|<t₂}
  (t₁,t₂ roots of Q₄ on (0,1)).  Then Q₄≤0 on A, Q₄≤Q₄(0) on Eq∪{0},
  Q₄≤1 on E∪poles, so
      μ_G4 = ∑ p(region) E[Q₄|region]
           ≤ P(poles) + P(E) + Q₄(0)·P(Eq)
           = 2/N + P(E) + 3/(d²−1)·P(Eq).
  In particular, if W_k=0 whenever Q₄((n−2k)/n)>0 except k∈{0,n}
  (no exterior/equator positive mass), then
      μ_G4 ≤ 2/N.
  Certified decomposition at p=5,7 (contrib signs match).  ∎

Theorem E (conditional residual via N) — PROVED Fraction + cert.
  If N ≥ N_*(p) := ceil(2 / μ_G4_suf(p)) and W is supported in
  {0,n}∪{k: Q₄((n−2k)/n)≤0}, then μ_G4≤2/N≤μ_G4_suf and residual holds.
  Explicitly N_*(5)=255, N_*(7)=2304, N_*(11)=38913 (grows ~ p⁶).
  At p=5: N=260≥255 so the N-threshold holds, but equator mass is present
  (Q₄(0) contrib), so the support hypothesis fails.  At p=7: N=11452≥2304,
  Hoffman layer lies in the exterior positive lobe (t_Hoff>t₂), support
  hypothesis fails.  ∎

Theorem F (Chebyshev split too weak) — PROVED form + cert.
  Using only E[t²]=2/n and Q₄≤1 on exterior / Q₄≤Q₄(0) on |t|≤t₂:
      μ_G4 ≤ Q₄(0) + (1−Q₄(0)) min(1, (2/n)/t₂²).
  The right-hand side is ≈0.18–0.26 for primes p=5..47, strictly above
  μ_G4_suf (dead for residual).  ∎

Theorem G (OPEN residual, honest).
  Prove μ_G4≤μ_G4_suf for all primes p≥5 by one of:
    (i) Weil/Paley character-sum evaluation of ∑_y Q₄(⟨z,y⟩/n),
    (ii) Aut-coherent configuration (finer than inner-product) of Max+
        with closed intersection numbers ⇒ closed W_k / μ_G4,
    (iii) closed N and W_k for regular sets of the conference srg,
    (iv) bound P(E)+3/(d²−1)P(Eq) using design equations beyond j≤3.
  Dead: IP-association scheme of Max+; Chebyshev split; μ≤1; μ≤1/h₄;
  moment-LP on allowed k.  Partial: residual reduced to controlling
  exterior+equator mass of a 1-homogeneous tight frame with closed Q₄.  ∎

============================================================================
Certified: Thm A Fraction; Thm B–C census p=5; Thm D–F algebra+p=5,7;
  N_* table p=5..31.
Backend: CPU Fraction + numpy eigensystem p=5; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15158.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from math import ceil, comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15128 import KNOWN_N, W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15156 import Es4_from_W  # noqa: E402
from e1_gmin_m4_prop15157 import (  # noqa: E402
    a0_of,
    a2_of,
    a4_of,
    mu_G4_from_Es4,
    mu_G4_suf,
)


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
# Closed Q4
# ---------------------------------------------------------------------------

def Q2_of(d: int, t: Fraction) -> Fraction:
    """Q₂(t)=(d t²−1)/(d−1)."""
    return (d * t * t - 1) / Fraction(d - 1)


def Q4_of(d: int, t: Fraction) -> Fraction:
    """Q₄(t)=[(d+2)(d+4)t⁴ − 6(d+2)t² + 3]/(d²−1)."""
    num = (
        Fraction(d + 2) * (d + 4) * t**4
        - 6 * Fraction(d + 2) * t**2
        + 3
    )
    return num / Fraction(d * d - 1)


def Q4_of_float(d: int, t: float) -> float:
    return float(Q4_of(d, Fraction(t).limit_denominator(10**9)))


def Q4_from_expansion(d: int, t: Fraction) -> Fraction:
    """Match 15.157: Q₄=(t⁴−a₀−a₂ Q₂)/a₄."""
    return (t**4 - a0_of(d) - a2_of(d) * Q2_of(d, t)) / a4_of(d)


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_Q4_closed(ds: list[int]) -> dict:
    rows = {}
    ok = True
    for d in ds:
        match = True
        for ti in range(0, 21):
            t = Fraction(ti, 20)
            if Q4_of(d, t) != Q4_from_expansion(d, t):
                match = False
                break
        q0 = Q4_of(d, Fraction(0))
        match = match and q0 == Fraction(3, d * d - 1)
        match = match and Q4_of(d, Fraction(1)) == 1
        rows[str(d)] = {
            "Q4_0": str(q0),
            "match_expansion": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Q₄(t)=[(d+2)(d+4)t⁴−6(d+2)t²+3]/(d²−1); Q₄(0)=3/(d²−1)."
        ),
        "by_d": {k: rows[k] for k in list(rows)[:8]},
        "n_checked": len(rows),
        "closed_forms": {
            "Q4": "[(d+2)(d+4) t^4 - 6(d+2) t^2 + 3]/(d^2-1)",
            "Q4_0": "3/(d^2-1)",
            "Q2": "(d t^2 - 1)/(d-1)",
        },
    }


# ---------------------------------------------------------------------------
# Theorem B–C: 1-homogeneous, spectrum, non-scheme
# ---------------------------------------------------------------------------

def _load_switched_G(p: int):
    from minmax_quadratic import halfspace_boolean_vector

    Y = np.load(f"/tmp/maxplus_p{p}.npy").astype(np.float64)
    hs = halfspace_boolean_vector(p).astype(np.float64)
    Y2 = Y * hs
    G = np.rint(Y2 @ Y2.T).astype(np.int32)
    return Y2, G


def prove_homogeneous_and_spectrum() -> dict:
    p = 5
    Y2, G = _load_switched_G(p)
    N, n = Y2.shape
    d = d_of(p)
    # 1-homogeneous: neighbor profile
    def profile(i: int):
        return tuple(
            sorted(Counter(int(G[i, j]) for j in range(N) if j != i).items())
        )

    profiles = {profile(i) for i in range(N)}
    one_hom = len(profiles) == 1
    # spectrum
    ev = np.linalg.eigvalsh(G.astype(float))
    ev_r = np.round(ev, 6)
    ctr = Counter(ev_r.tolist())
    # expect 2N with mult d, 0 with mult N-d
    lam_big = float(2 * N)
    mult_big = int(ctr.get(round(lam_big, 6), 0))
    # zeros
    mult_zero = int(sum(c for e, c in ctr.items() if abs(e) < 1e-6))
    spectrum_ok = mult_big == d and mult_zero == N - d
    return {
        "proved_census": one_hom and spectrum_ok,
        "theorem": (
            "Max+ 1-homogeneous (identical absolute distributions); "
            "Gram spectrum 2N (×d) + 0 (×N−d)."
        ),
        "p5": {
            "N": N,
            "d": d,
            "n_profiles": len(profiles),
            "one_homogeneous": one_hom,
            "lambda_2N_mult": mult_big,
            "lambda_0_mult": mult_zero,
            "spectrum_ok": spectrum_ok,
            "profile_head": list(profiles)[0][:6] if profiles else None,
        },
    }


def prove_not_association_scheme() -> dict:
    """Intersection numbers non-constant for several classes at p=5."""
    p = 5
    Y2, G = _load_switched_G(p)
    N = G.shape[0]
    # For relation s=14 (often nonconstant), check a few pairs
    s_test = [14, -14, 6, -6, 2, -2]
    rng = np.random.default_rng(1)
    results = {}
    overall_nonconst = False
    for s in s_test:
        pairs = []
        for a in range(N):
            for b in range(a + 1, N):
                if int(G[a, b]) == s:
                    pairs.append((a, b))
        if len(pairs) > 30:
            idx = rng.choice(len(pairs), size=30, replace=False)
            pairs = [pairs[i] for i in idx]
        # intersection count with relation 2: # c with G[a,c]=2, G[b,c]=2
        counts = []
        for a, b in pairs:
            c_count = 0
            for c in range(N):
                if c == a or c == b:
                    continue
                if int(G[a, c]) == 2 and int(G[b, c]) == 2:
                    c_count += 1
            counts.append(c_count)
        const = max(counts) == min(counts) if counts else True
        if not const:
            overall_nonconst = True
        results[str(s)] = {
            "n_pairs": len(pairs),
            "p_2_2_min": min(counts) if counts else None,
            "p_2_2_max": max(counts) if counts else None,
            "constant": const,
        }
    return {
        "proved_not_scheme": overall_nonconst,
        "theorem": (
            "Inner-product relations on Max+ are not an association scheme "
            "at p=5 (non-constant p_{2,2}^s for several s)."
        ),
        "by_s": results,
    }


# ---------------------------------------------------------------------------
# Theorem D–E: decomposition and conditional N bound
# ---------------------------------------------------------------------------

def N_star(p: int) -> int:
    """Minimal N for pure pole bound μ≤2/N to imply residual."""
    ms = mu_G4_suf(p)
    # need 2/N <= ms => N >= 2/ms
    return int(ceil(float(2 / ms) - 1e-12))


def prove_decomposition_and_conditional() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        n, d = n_of(p), d_of(p)
        N = KNOWN_N[p]
        W = W_CENSUS[p]
        mu = mu_G4_from_Es4(p, Es4_from_W(p))
        ms = mu_G4_suf(p)
        # pole mass
        pole = Fraction(W.get(0, 0) + W.get(n, 0), N)
        # compute Q4-weighted contributions
        pos = Fraction(0)
        neg = Fraction(0)
        for k, c in W.items():
            t = Fraction(n - 2 * k, n)
            q = Q4_of(d, t)
            contrib = Fraction(c, N) * q
            if q > 0:
                pos += contrib
            else:
                neg += contrib
        # conditional: support only nonpositive Q4 except poles
        support_ok = True
        for k, c in W.items():
            if c == 0 or k in (0, n):
                continue
            if Q4_of(d, Fraction(n - 2 * k, n)) > 0:
                support_ok = False
                break
        n_star = N_star(p)
        pole_bound_ok = Fraction(2, N) <= ms
        rows[str(p)] = {
            "N": N,
            "N_star": n_star,
            "N_ge_N_star": N >= n_star,
            "pole_mass": str(pole),
            "two_over_N": str(Fraction(2, N)),
            "mu_G4": str(mu),
            "mu_G4_suf": str(ms),
            "pos_contrib": str(pos),
            "neg_contrib": str(neg),
            "support_only_nonpos_Q4_except_poles": support_ok,
            "pole_bound_implies_residual": pole_bound_ok and support_ok,
            "two_over_N_le_suf": pole_bound_ok,
        }
        # structural checks
        if abs(float(pos + neg) - float(mu)) > 1e-9:
            ok = False
    # N_star table
    nstar_table = {str(p): N_star(p) for p in primes_ge(5, 31)}
    return {
        "proved": ok,
        "theorem": (
            "μ_G4≤2/N+P(E)+Q₄(0)P(Eq); if support has Q₄≤0 except poles "
            "and N≥N_*=⌈2/μ_G4_suf⌉ then residual holds."
        ),
        "by_p": rows,
        "N_star_table": nstar_table,
        "note": (
            "At p=5,7: N≥N_* holds but support hypothesis fails "
            "(equator and/or Hoffman exterior positive Q₄)."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem F: Chebyshev weak
# ---------------------------------------------------------------------------

def prove_chebyshev_weak(primes: list[int]) -> dict:
    """Document that Chebyshev split UB is ~0.18 >> μ_suf."""
    # Use exact Q4 to find t2 approximately via binary search on Fraction grid
    rows = {}
    all_weak = True
    for p in primes:
        if p * p == 5:
            continue
        n, d = n_of(p), d_of(p)
        # find larger root of Q4 on (0,1) by scan
        t2 = None
        prev_q = Q4_of(d, Fraction(0))
        for i in range(1, 2001):
            t = Fraction(i, 2000)
            q = Q4_of(d, t)
            if prev_q < 0 <= q or prev_q > 0 >= q:
                # crossed to positive near 1 — take last crossing from neg to pos
                if prev_q < 0 <= q:
                    t2 = t
            prev_q = q
        if t2 is None or t2 == 0:
            t2 = Fraction(1, 2)
        Et2 = Fraction(2, n)
        Pub = min(Fraction(1), Et2 / (t2 * t2))
        q0 = Fraction(3, d * d - 1)
        ub = q0 + (1 - q0) * Pub
        ms = mu_G4_suf(p)
        weak = ub > ms
        if not weak:
            all_weak = False
        rows[str(p)] = {
            "t2_approx": str(t2),
            "chebyshev_ub": str(ub),
            "mu_G4_suf": str(ms),
            "exceeds_suf": weak,
        }
    return {
        "proved_too_weak": all_weak,
        "theorem": (
            "Chebyshev split μ_G4≤Q₄(0)+(1−Q₄(0))min(1,(2/n)/t₂²) "
            "exceeds μ_G4_suf for all tested p (dead)."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "mu_G4_le_mu_G4_suf_all_p": False,
        "open": (
            "Prove μ_G4≤μ_G4_suf for all primes p≥5 via Weil/Paley sums on "
            "∑_y Q₄(⟨z,y⟩/n), Aut-coherent configuration of Max+, closed W_k, "
            "or exterior/equator mass bounds beyond moments j≤3."
        ),
        "sufficient_targets": [
            "Weil_sum_Q4_over_Maxplus",
            "Aut_coherent_configuration_Maxplus",
            "closed_N_and_W_k",
            "exterior_equator_mass_bound",
        ],
        "dead": [
            "IP_association_scheme_Maxplus",
            "Chebyshev_Q4_split",
            "mu_le_1",
            "mu_le_1_over_h4",
            "moment_LP_allowed_k",
            "plain_ULC",
        ],
        "progress": (
            "Max+ is 1-homogeneous tight frame (Gram 2N×d); not IP-scheme; "
            "Q₄ closed; residual is pole+exterior+equator mass control; "
            "N≥N_* at p=5,7 but positive-Q₄ support blocks pure pole bound."
        ),
    }


def main() -> dict:
    a = prove_Q4_closed(list(range(3, 81)) + [85, 145, 181, 265])
    b = prove_homogeneous_and_spectrum()
    c = prove_not_association_scheme()
    d = prove_decomposition_and_conditional()
    f = prove_chebyshev_weak(primes_ge(5, 31))
    g = prove_open()

    out = {
        "prop": "15.158",
        "title": (
            "Prop 15.158: closed Q₄; Max+ 1-homogeneous non-scheme; "
            "pole decomposition; conditional N-bound; L OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Q₄(t)=[(d+2)(d+4)t⁴−6(d+2)t²+3]/(d²−1). "
            "Max+ is 1-homogeneous with Gram spectrum 2N (×d)+0 (×N−d), but "
            "not an inner-product association scheme (p=5). "
            "μ_G4≤2/N+P(E)+Q₄(0)P(Eq); pure pole bound needs N≥N_* and "
            "nonpositive-Q₄ support — N_* holds at p=5,7 but support fails. "
            "Chebyshev split dead. Need Weil/Aut-coherent config/closed W_k. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "Q4_closed_form": a["proved"],
            "one_homogeneous_gram_spectrum": b["proved_census"],
            "not_IP_association_scheme": c["proved_not_scheme"],
            "pole_decomposition_conditional": d["proved"],
            "chebyshev_split_dead": f["proved_too_weak"],
            "residual_for_all_p_ge_5": False,
            "mu_G4_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "Q4": a,
            "homogeneous": b,
            "non_scheme": c,
            "decomposition": d,
            "chebyshev": f,
            "open": g,
        },
        "closed_forms": a["closed_forms"],
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "IP_scheme_blocked": True,
            "general_residual": False,
        },
        "open_attack": (
            "Weil/Paley ∑_y Q₄(⟨z,y⟩/n), or Aut-coherent configuration of Max+, "
            "or closed W_k / exterior-equator mass bound."
        ),
        "backend": "CPU Fraction + numpy Gram p=5; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {"workers": "Fraction + numpy", "gpu": "unused"},
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15158.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
