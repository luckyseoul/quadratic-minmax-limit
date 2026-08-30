#!/usr/bin/env python3
"""
Prop 15.595 — The δ-hierarchy: which GOAL units one bound ‖δ‖² ≤ c·n closes,
and which it provably cannot.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SCOPE NOTE (read first).  The corrected ``e1_closed_general()`` is False.
The public claim is gated on four current GOAL units: required bi-tight
levels 2/3 (True by 15.720), residual_ii_k_ge_4p (False),
type_I_multilevel (False), and lemma_D (True).  The delta hierarchy below
belongs to an older spectral route; 15.720 bypassed its ``leftover 1`` as an
acceptance gate.  This proposition remains useful dependency information but
does not identify the current theorem blockers.

SETUP (15.594):  δ = m₄⁺ − m₄_part ∈ ker(4pI − T), V = 24‖δ‖².

PROVED — the hierarchy of ‖δ‖² requirements (all exact, data-free)
  A. Three distinct GOAL/chain units are bounds on the SAME ‖δ‖²:

     | unit                                   | requires ‖δ‖² ≤        | ~limit |
     |----------------------------------------|------------------------|--------|
     | leftover 1  (λ_min(Φ) ≥ 6)             | n(λ̄−6)²/48             | n/12   |
     | leftover 3  (3A+B > 0)                 | c₃(p)·n/24             | ~2.9n  |
     | residual-(i) R≤2p (15.217 delta_room)  | (p²−1)n(3p²−47)/(24(p²−5)) | ~n²/8 |

     so  ‖δ‖² ≤ n/12  ⟹  all three, the binding one being leftover 1.
     Verified numerically: the ordering
        n(λ̄−6)²/48  <  c₃n/24  <  delta_room_for_R(p)
     holds at every tested prime (5..47).
  B. Measured ‖δ‖²/n = 0.9089 (p=5), 0.2085 (p=7), 0.01941 (p=11).
     Against the binding threshold ≈ 1/12: FAILS at p=5,7 (both handed to
     census — 15.275 L and the 109/2863 census) and HOLDS at p=11 with
     4.3× margin.  Three points, recorded, NOT extrapolated.

PROVED (negative) — leftover 2 does NOT reduce to ‖δ‖²
  C. Leftover 2 lives on size-k subsets, k=4p.  Any four-point statement
     there is Σ_{S⊆G} m₄(S); the δ-part is bounded only by Cauchy–Schwarz
     over the C(k,4) subsets:
        |Σ_{S⊆G} δ(S)| ≤ √C(k,4)·‖δ‖ ,
     while the signal scale is ≈ C(k,4)/n.  The ratio is
        error/signal ≈ n^{3/2}/(√12·√C(k,4)) ≈ p/11.4 → ∞,
     crossing 1 at p ≈ 11.  So even the optimal ‖δ‖² ≤ n/12 gives an error
     bar that EXCEEDS the entire four-point signal for every p ≥ 11 and
     diverges.  **Do not attempt leftover 2 through an L² bound on δ.**
     Leftover 2 remains genuinely separate (Max−, Walsh/minus-slice route).

HISTORICAL CONSEQUENCE
  D. In the retired spectral formulation, leftovers 1 and 3 shared the
     delta inequality while leftover 2 was independent.  That historical
     route had two roots.  It is not the current acceptance ledger: after
     15.720, the two live roots are the multi-level Type-I gate and the
     non-Walsh residual-(ii) gate at even k>=4p.

Writes evidence/e1_gmin_m4_prop15595.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from math import comb, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    return p * p + 1


def lambda_bar(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def req_leftover1(p: int) -> Fraction:
    """||delta||^2 <= n(lbar-6)^2/48."""
    return Fraction(n_of(p), 48) * (lambda_bar(p) - 6) ** 2


def req_leftover3(p: int, c3: Fraction) -> Fraction:
    """||delta||^2 <= c3 n / 24."""
    return c3 * n_of(p) / 24


def req_residual_i(p: int) -> Fraction:
    """15.217 delta_room_for_R: (p^2-1)(p^2+1)(3p^2-47)/(24(p^2-5))."""
    return Fraction((p - 1) * (p + 1) * (p * p + 1) * (3 * p * p - 47),
                    24 * (p * p - 5))


def hierarchy(p: int, c3: Fraction = Fraction(1550, 100)) -> dict:
    r1, r3, ri = req_leftover1(p), req_leftover3(p, c3), req_residual_i(p)
    n = n_of(p)
    return {"p": p, "n": n,
            "req_leftover1": r1, "req_leftover3": r3, "req_residual_i": ri,
            "req_leftover1_per_n": float(r1 / n),
            "req_leftover3_per_n": float(r3 / n),
            "req_residual_i_per_n": float(ri / n),
            "ordered": r1 < r3 < ri,
            "binding": "leftover1"}


def leftover2_error_vs_signal(p: int, c: Fraction = Fraction(1, 12)) -> dict:
    """Cauchy-Schwarz error bar on sum_{S in G} delta(S) vs the signal scale."""
    n, k = n_of(p), 4 * p
    Ck4 = comb(k, 4)
    err = sqrt(Ck4) * sqrt(float(c) * n)
    signal = Ck4 / n
    return {"p": p, "k": k, "C_k_4": Ck4, "error_bar": err,
            "signal_scale": signal, "ratio": err / signal,
            "useless": err >= signal}


def closes_with(c: Fraction, primes=(5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)) -> dict:
    """Given a hypothetical ||delta||^2 <= c*n, which units close at which primes."""
    out = {}
    for p in primes:
        n = n_of(p)
        out[p] = {"leftover1": c * n <= req_leftover1(p),
                  "leftover3": c * n <= req_leftover3(p, Fraction(1550, 100)),
                  "residual_i": c * n <= req_residual_i(p)}
    return out


def main():
    t0 = time.time()
    out = {"prop": "15.595",
           "title": "delta-hierarchy: what one bound closes, and what it provably cannot"}
    primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    out["hierarchy"] = {}
    ordered_all = True
    for p in primes:
        h = hierarchy(p)
        ordered_all &= h["ordered"]
        out["hierarchy"][p] = {k: (str(v) if isinstance(v, Fraction) else v)
                               for k, v in h.items()}
    out["hierarchy_ordered_all_primes"] = ordered_all
    out["leftover2_no_reduction"] = {
        p: leftover2_error_vs_signal(p) for p in (5, 7, 11, 13, 17, 23)}
    out["leftover2_useless_from_p11"] = all(
        leftover2_error_vs_signal(p)["useless"] for p in (11, 13, 17, 23))
    out["measured_delta_sq_per_n"] = {5: 0.9089, 7: 0.2085, 11: 0.01941}
    out["closes_with_n_over_12"] = {
        str(k): v for k, v in closes_with(Fraction(1, 12)).items()}
    out["historical_delta_route_roots"] = [
        "R1 (retired): ||P_{E_4p} m4+||^2 <= n/12",
        "R2 (still related to live residual ii): Max- minus-slice/Walsh",
    ]
    out["current_open_gates"] = [
        "non-Walsh residual (ii), even k>=4p",
        "multi-level Type I",
    ]
    out["scope_note"] = (
        "e1_closed_general()=False and agrees with the current four-unit gate; "
        "the delta inequality is no longer an acceptance unit after 15.720."
    )
    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    (ROOT / "evidence" / "e1_gmin_m4_prop15595.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.595  historical delta hierarchy; current gates unchanged")
    print(f"  hierarchy r1 < r3 < residual_i at all primes 5..47: {ordered_all}")
    for p in (5, 11, 23, 47):
        h = out["hierarchy"][p]
        print(f"    p={p:3d}: r1={h['req_leftover1_per_n']:.4f} n  "
              f"r3={h['req_leftover3_per_n']:.3f} n  ri={h['req_residual_i_per_n']:.2f} n")
    print(f"  leftover 2 unreachable by an L2 delta bound from p=11 on: "
          f"{out['leftover2_useless_from_p11']}")
    for p in (5, 11, 23):
        e = out["leftover2_no_reduction"][p]
        print(f"    p={p:3d}: error/signal = {e['ratio']:.2f}")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()
