#!/usr/bin/env python3
"""
Prop 15.595 — The δ-hierarchy: which GOAL units one bound ‖δ‖² ≤ c·n closes,
and which it provably cannot.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SCOPE NOTE (read first).  `e1_closed_general()` returns True, but the repo's
own rule (e1_main_chain_status) records this as "a separate wiring fact":
the public claim is gated on FOUR GOAL units —
    phi_F_ge_6 (leftover 1) · residual_ii_k_ge_4p (leftover 2)
    · type_I_multilevel (leftover 3) · lemma_D
of which lemma_D is already True and the other three are False.  Nothing
below changes any of that; this proposition only maps the dependencies.

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

CONSEQUENCE
  D. Of the four GOAL units: lemma_D is True; leftovers 1 and 3 share one
     open inequality; leftover 2 is independent of it.  The E(1) leftover
     program therefore has exactly TWO open roots, not three:
        (R1)  ‖P_{E_{4p}} m₄⁺‖² ≤ n/12      [closes leftovers 1 and 3]
        (R2)  the Max− minus-slice/Walsh statement  [leftover 2 alone]

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
    out["open_roots"] = [
        "R1: ||P_{E_4p} m4+||^2 <= n/12  -> closes leftovers 1 and 3",
        "R2: Max- minus-slice/Walsh       -> leftover 2 alone",
    ]
    out["scope_note"] = ("e1_closed_general()=True is a separate wiring fact "
                         "(e1_main_chain_status); the gate is the four GOAL "
                         "units, three of which are False.")
    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    (ROOT / "evidence" / "e1_gmin_m4_prop15595.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.595  delta-hierarchy; two open roots, not three")
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
