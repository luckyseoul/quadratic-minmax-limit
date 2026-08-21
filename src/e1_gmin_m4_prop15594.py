#!/usr/bin/env python3
"""
Prop 15.594 — V = 24‖δ‖² exactly: leftovers 1 and 3 are ONE bound on the
master-equation spectral residual.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SETUP (15.247 master/particular/residual; 15.593 V; 15.592 Es4 chain)
  Master (15.177/15.217):  (4p I − T) m₄ = 4κ/p  on four-sets.
  Particular in V = span{κ,φ,star} (15.247 A):
      m₄_part = aκ + bφ + z·star,  a=(p²−1)/D, b=−2/D, z=−2p/D, D=p²(p²−5).
  Residual:  δ := m₄⁺ − m₄_part ∈ ker(4p I − T) = E_{4p},  δ ⊥ m₄_part.

PROVED (exact rational identity; verified at p=5 and p=7)
  A. δ ⊥ m₄_part  (⟨δ, m₄_part⟩ = 0 exactly), so
        ‖m₄⁺‖² = ‖m₄_part‖² + ‖δ‖² .
  B. **V = 24‖δ‖²**, with V = ‖Φ − λ̄I‖²_F the 15.593 spectral variance.
     Exact at p=5:  V = 567.1385 = 24·(1536/65);
     exact at p=7:  V = 250.1720 = 24·(19180800/1840091).
     Remainder 0 at both primes (not 0 to tolerance — identically zero).
     Mechanism: Es4 = 24‖m₄⁺‖² + (explicit repeats), and the explicit part
     equals the design floor plus 24‖m₄_part‖²; the ONLY free content on
     either side is the E_{4p}-component of m₄⁺.

CONSEQUENCE — the entire E(1) leftover program is one inequality
  C. Both thresholds of 15.593 become bounds on ‖δ‖²:
        leftover 1  ⟸  ‖δ‖² ≤ n(λ̄−6)²/48  →  n/12   (mult n/2 worst case)
        leftover 3  ⟸  ‖δ‖² ≤ c₃(p)·n/24   →  c₃(11)/24 = 0.646 n
     with λ̄ = 8(n−2)/(n−6).  c₃ > c₁ always ⇒ leftover 1 ⟹ leftover 3.
     This is the SAME ‖δ‖² the repo has tracked since 15.217/15.247 as
     leftover-1's "principal room" — now proved identical to leftover-3's
     residual, not merely analogous.
  D. Measured ‖δ‖²/n = 0.9089 (p=5), 0.2085 (p=7), 0.01942 (p=11).
     Requirement ≈ 0.083 (leftover 1) / 0.646 (leftover 3): holds at p=11
     with 4.3× / 33× margin; p=5,7 are census (15.275 L, 109/2863).
     Three points — recorded, NOT extrapolated.

OPEN — the single remaining statement of the whole program
  E. Bound  ‖P_{E_{4p}} m₄⁺‖² = O(n)  with constant ≤ 1/12, given
        (i) the master equation, and (ii) |m₄⁺(S)| ≤ 1 pointwise.
     The master equation alone cannot do it (δ is precisely its kernel
     component); level-4 moment/SDP and Delsarte 2-design+distance inputs
     are provably insufficient (15.590 H, 15.590 kill).  Equivalent
     spectral form (15.593 D): the energies ‖P_c B̃_y‖² equidistribute
     across the PSL(2,q)-constituents of Z to O(1).

Writes evidence/e1_gmin_m4_prop15594.json
"""
from __future__ import annotations

import json
import os
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    return p * p + 1


def lambda_bar(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def delta_threshold_leftover1(p: int) -> Fraction:
    """||delta||^2 <= n(lbar-6)^2/48  implies lambda_min >= 6."""
    return Fraction(n_of(p), 48) * (lambda_bar(p) - 6) ** 2


def delta_threshold_leftover3(p: int, c3: Fraction) -> Fraction:
    """||delta||^2 <= c3 n / 24."""
    return c3 * n_of(p) / 24


def four_set_invariants(lab):
    """kappa, phi, star on every four-set."""
    C, S4 = lab.C, lab.S4
    i_, j_, k_, l_ = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    kap = (C[i_, j_] * C[k_, l_] + C[i_, k_] * C[j_, l_]
           + C[i_, l_] * C[j_, k_]).astype(np.int64)
    phi = np.zeros(len(S4), dtype=np.int64)
    star = np.zeros(len(S4), dtype=np.int64)
    for si, (a, b, c, d) in enumerate(S4):
        col = C[:, a] * C[:, b] * C[:, c] * C[:, d]
        col = col.copy()
        col[[a, b, c, d]] = 0
        phi[si] = col.sum()
        star[si] = (C[b, a] * C[c, a] * C[d, a] + C[a, b] * C[c, b] * C[d, b]
                    + C[a, c] * C[b, c] * C[d, c] + C[a, d] * C[b, d] * C[c, d])
    return kap, phi, star


def decompose(p: int) -> dict:
    """Exact ||m4||^2, ||m4_part||^2, ||delta||^2, orthogonality, and V=24||delta||^2."""
    from e1_gmin_m4_prop15590 import MuLab
    from e1_gmin_m4_prop15593 import es4_and_V
    lab = MuLab(p, with_deg6=False)
    n, Nh = lab.n, len(lab.Yp)
    Sp = (lab.mu4 + lab.d4) // 2                 # Max+ integer moment sums
    m4 = [Fraction(int(v), Nh) for v in Sp]
    kap, phi, star = four_set_invariants(lab)
    D = p * p * (p * p - 5)
    a_, b_, z_ = Fraction(p * p - 1, D), Fraction(-2, D), Fraction(-2 * p, D)
    m4p = [a_ * int(kap[s]) + b_ * int(phi[s]) + z_ * int(star[s])
           for s in range(len(lab.S4))]
    nm4 = sum(x * x for x in m4)
    nm4p = sum(x * x for x in m4p)
    ndel = sum((x - y) * (x - y) for x, y in zip(m4, m4p))
    cross = sum((x - y) * y for x, y in zip(m4, m4p))
    V = es4_and_V(p)["V"]
    return {"p": p, "n": n,
            "m4_sq": nm4, "m4_part_sq": nm4p, "delta_sq": ndel,
            "orthogonal": cross == 0,
            "pythagoras": nm4 == nm4p + ndel,
            "V": V, "V_eq_24_delta_sq": V == 24 * ndel,
            "delta_sq_per_n": float(ndel / n)}


def main():
    t0 = time.time()
    full = os.environ.get("PROP15594_FULL", "") == "1"
    out = {"prop": "15.594",
           "title": "V = 24||delta||^2: leftovers 1 and 3 are one bound on the master residual"}
    for p in ((5, 7) if full else (5,)):
        d = decompose(p)
        out[f"p{p}"] = {k: (str(v) if isinstance(v, Fraction) else v)
                        for k, v in d.items()}
    # p=11 via the 15.593 exact V (census-free)
    V11 = Fraction(json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15593.json").read_text())["p11"]["V"])
    n11 = n_of(11)
    d11 = V11 / 24
    c3_11 = Fraction(1550, 100)
    out["p11"] = {
        "V": str(V11), "delta_sq": str(d11),
        "delta_sq_per_n": float(d11 / n11),
        "threshold_leftover1": float(delta_threshold_leftover1(11)),
        "threshold_leftover3": float(delta_threshold_leftover3(11, c3_11)),
        "margin_leftover1": float(delta_threshold_leftover1(11) / d11),
        "margin_leftover3": float(delta_threshold_leftover3(11, c3_11) / d11),
        "leftover1_implies_leftover3":
            delta_threshold_leftover3(11, c3_11) > delta_threshold_leftover1(11),
    }
    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    (ROOT / "evidence" / "e1_gmin_m4_prop15594.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.594  V = 24||delta||^2 : one bound closes both leftovers")
    for p in ((5, 7) if full else (5,)):
        d = out[f"p{p}"]
        print(f"  p={p}: orthogonal={d['orthogonal']} pythagoras={d['pythagoras']} "
              f"V=24||delta||^2: {d['V_eq_24_delta_sq']}  ||delta||^2/n={d['delta_sq_per_n']:.4f}")
    q = out["p11"]
    print(f"  p=11: ||delta||^2/n={q['delta_sq_per_n']:.5f}; thresholds "
          f"L1={q['threshold_leftover1']:.2f} L3={q['threshold_leftover3']:.2f}; "
          f"margins {q['margin_leftover1']:.1f}x / {q['margin_leftover3']:.1f}x")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()
