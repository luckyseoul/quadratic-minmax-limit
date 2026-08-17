#!/usr/bin/env python3
"""
Prop 15.338 — Named k-vector for p≡3 that specialises to 15.337;
the interpolating pair-rule α is p=7-only Max+.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.337 flags.  Floor stays OPEN.  Still 15.x.

============================================================================
Setup.  Centrally symmetric |D|=p(p−1)/2 through 0 exists iff
p≡3 (mod 4).  Lines through 0 carry (p−1)/2 antipodal pairs.
k(ℓ)=#pairs of D on ℓ.  t² labels pairs on a finite-slope line
by squares □⊂F_p.  The predicate “t²−α is a square (incl. 0)”
has three levels (certified p=7,11,19):
    α=0 ⇒ k=(p−1)/2,   α∈□ ⇒ k=(p+1)/4,   α∉□∪{0} ⇒ k=(p−3)/4.

============================================================================
Theorem A — PROVED (Max+-free arithmetic; p≡3 (mod 4), p≥7).
  One full line, one empty line,
      n_high=(p−3)/4  lines with k=(p+1)/4,
      n_low=(3p−1)/4  lines with k=(p−3)/4
  give
      ∑k = (p(p−1)−2)/4,   |D|=p(p−1)/2.
  At p=7 this is exactly 15.337’s (3,2,1⁵,0).  Fail-when-wrong:
  claim n_high=1 at p=11 (that is the p=7 value; live formula is 2).  ∎

Theorem B — PROVED (character-sum census; p=7,11,19).
  On □, #{u∈□ : u−α ∈ □∪{0}} equals (p−1)/2, (p+1)/4, (p−3)/4
  according as α=0, α∈□, α nonsquare.  Fail: claim α∈□ gives k=1
  at p=7 (it gives 2).  ∎

Theorem C — PROVED as a negative (certified p=7,11,19,23).
  The unique deg≤5 interpolant of the p=7 isolated α-table,
      α(s)=−2−s+s²(1+s)^{p−4}   (F_p),
  with full x-axis and empty y-axis, is Max+ and isolated at p=7,
  and is not Max+ at p=11,19,23 (sizes 59,185,265 ≠ 55,171,253).
  Fail-when-wrong: claim the interpolant is Max+ at p=11.  ∎

Theorem D — OPEN.  FLT / deg≤2 / consecutive-square constant-α
  rules give 0 Max+ at p=11.  Pair-selection on the named k-vector
  is not yet a Gauss formula that stays Max+ for p≡3.  Q_τ / phi_F
  not imported.  ∎

============================================================================
Backend: F_p arithmetic + one Paley C per prime.  Writes
evidence/e1_gmin_m4_prop15338.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402

# pair_rule helpers live in scratch for hunts; reimplement the tiny bits here
# so the shipped unit does not import scratch.


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def n_high(p: int) -> int:
    return (p - 3) // 4


def k_high(p: int) -> int:
    return (p + 1) // 4


def k_low(p: int) -> int:
    return (p - 3) // 4


def k_full(p: int) -> int:
    return (p - 1) // 2


def n_low(p: int) -> int:
    return (3 * p - 1) // 4


def sum_k_named(p: int) -> int:
    return (
        k_full(p)
        + 0
        + n_high(p) * k_high(p)
        + n_low(p) * k_low(p)
    )


def size_from_k(p: int) -> int:
    return 1 + 2 * sum_k_named(p)


def k_hist_p7() -> dict[int, int]:
    return {3: 1, 2: 1, 1: 5, 0: 1}


def squares(p: int) -> set[int]:
    return {pow(i, 2, p) for i in range(1, p)}


def k_of_alpha(p: int, a: int) -> int:
    sq = squares(p)
    return sum(1 for u in sq if chi_p(p, (u - a) % p) >= 0)


def alpha_interpolant(p: int, s: int) -> int:
    """deg-5 interpolant of the p=7 isolated table, written in p."""
    return ((-2 - s) + (s * s) * pow((1 + s) % p, p - 4, p)) % p


def pair_ts(p: int) -> list[int]:
    seen, out = set(), []
    for t in range(1, p):
        if t in seen:
            continue
        seen.add(t)
        seen.add((-t) % p)
        out.append(min(t, (-t) % p))
    return out


def build_D_interpolant(p: int) -> np.ndarray:
    ts = pair_ts(p)
    mask = np.zeros(p * p, dtype=bool)
    mask[0] = True
    for t in range(1, p):
        mask[t] = True
    for s in range(1, p):
        a = alpha_interpolant(p, s)
        for t in ts:
            u = (t * t) % p
            if chi_p(p, (u - a) % p) >= 0:
                x, y = t % p, (s * t) % p
                mask[x + p * y] = True
                mask[(-x) % p + p * ((-y) % p)] = True
    return mask


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in range(7, 80):
        if not is_prime(p) or p % 4 != 3:
            continue
        sk = sum_k_named(p)
        want = (p * (p - 1) - 2) // 4
        sz = size_from_k(p)
        if sk != want or sz != p * (p - 1) // 2:
            ok = False
        rows[str(p)] = {
            "n_high": n_high(p),
            "n_low": n_low(p),
            "k_full": k_full(p),
            "k_high": k_high(p),
            "k_low": k_low(p),
            "sum_k": sk,
            "size": sz,
        }
    # specialise p=7
    ok = ok and n_high(7) == 1 and k_high(7) == 2 and k_low(7) == 1
    ok = ok and k_full(7) == 3 and n_low(7) == 5
    # fail-when-wrong
    if n_high(11) == 1:
        ok = False
    ok = ok and n_high(11) == 2
    return {
        "proved": bool(ok),
        "p7_hist": {
            "full": k_full(7),
            "high": [n_high(7), k_high(7)],
            "low": [n_low(7), k_low(7)],
            "empty": 0,
        },
        "matches_15337": n_high(7) == 1 and k_high(7) == 2 and n_low(7) == 5,
        "n_high_p11": n_high(11),
        "by_p": {k: rows[k] for k in list(rows)[:8]},
        "theorem": "Named k-vector specialises to (3,2,1^5,0) at p=7.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (7, 11, 19):
        sq = squares(p)
        ns = [i for i in range(1, p) if i not in sq]
        k0 = k_of_alpha(p, 0)
        kh = {k_of_alpha(p, a) for a in sq}
        kl = {k_of_alpha(p, a) for a in ns}
        if k0 != k_full(p) or kh != {k_high(p)} or kl != {k_low(p)}:
            ok = False
        rows[str(p)] = {"k0": k0, "k_square": sorted(kh), "k_ns": sorted(kl)}
    # fail: α square gives k=1 at p=7
    if k_of_alpha(7, 1) == 1:
        ok = False
    ok = ok and k_of_alpha(7, 1) == 2
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "t²−α square realises the three named k-levels.",
    }


def prove_C() -> dict:
    from minmax_quadratic import paley_conference_prime_power

    # p=7 table recovered
    tab7 = {s: alpha_interpolant(7, s) for s in range(1, 7)}
    expect = {1: 5, 2: 6, 3: 4, 4: 6, 5: 3, 6: 6}
    ok = tab7 == expect
    rows = {}
    for p in (7, 11, 19, 23):
        mask = build_D_interpolant(p)
        C = paley_conference_prime_power(p).astype(np.float64)
        y = np.ones(p * p + 1, dtype=np.float64)
        y[1:] = np.where(mask, -1.0, 1.0)
        ev = float(np.max(np.abs(C @ y - p * y)))
        n = int(mask.sum())
        want = p * (p - 1) // 2
        rows[str(p)] = {
            "n": n,
            "want": want,
            "ev": ev,
            "maxplus": ev < 1e-4,
            "size_ok": n == want,
        }
        if p == 7:
            ok = ok and ev < 1e-4 and n == 21
        else:
            ok = ok and ev >= 1e-4
            if ev < 1e-4:
                ok = False
    # fail-when-wrong: Max+ at p=11
    if rows["11"]["maxplus"]:
        ok = False
    return {
        "proved": bool(ok),
        "alpha7": tab7,
        "by_p": rows,
        "theorem": "Interpolant α is isolated-Max+ at p=7 only.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "pair_rule_Maxplus_general": False,
        "note": (
            "k-vector named for p≡3.  Pair-selection that stays Max+ "
            "for p=11,19 is still OPEN.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.338  named k-vector; interpolant p=7-only", flush=True)
    A = prove_A()
    print(f"  A k-vector: {A['proved']} n_high(11)={A['n_high_p11']}", flush=True)
    B = prove_B()
    print(f"  B t²−α levels: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C interpolant p=7-only: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.338",
        "title": "Named k-vector for p≡3; interpolant α is p=7-only Max+",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "kvec_specialises_15337": A["proved"],
            "t2_alpha_three_levels": B["proved"],
            "interpolant_p7_only": C["proved"],
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15338.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
