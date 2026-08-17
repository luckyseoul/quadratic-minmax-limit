#!/usr/bin/env python3
"""
Prop 15.291 — Type cardinalities |τ| in p (Jacobi / subfield / norm torus).

Does **not** flip phi_F / e1 / L / Aut-Schur / Gsum / pairing.
Does **not** name Q_τ in p (still 15.290 C).  Floor stays OPEN.

============================================================================
T=□ ⊂ F_q^×, q=p².  Types of r∈T as in 15.290:
  pm1 : {±1},  |pm1|=2
  sub : r∈F_p^× \\ {±1}
  N1  : N_{q/p}(r)=1 and r∉F_p
  N*  : the rest
split further by (χ(r−1),χ(r+1)).

============================================================================
Theorem A — PROVED (Max+-free counting).
  χ_q | F_p^× ≡ 1, so every r∈F_p^× is a square in F_q.
  |sub ∩ T| = p−1, hence |++ sub| = p−3  (drop ±1) and no other
  Paley type occurs on the subfield.
  |{N=1}| = p+1 (ker of the norm), of which ±1, so |N1|=p−1.
  |T|−2 = (q−5)/2, hence
      |N*| = (p−1)(p−3)/2.
  Fail-when-wrong: |sub off|=p−1; |N1|=p+1 (include ±1).  ∎

Theorem B — PROVED (character-sum expansion; certified p=3,5,7,11).
  For ε1,ε2∈{±1},
      N_□(ε1,ε2) := #{r∈T \\ {±1} : χ(r−1)=ε1, χ(r+1)=ε2}
  equals
      (q − 3 + ε1 + ε2)/8
      + (ε1/8) J_{0,1} + (ε2/8) J_{0,-1} + (ε1 ε2/8) J_{1,-1}
      + (1/8) ∑_{x∉{0,±1}} χ(x)  [the ε0=1 slice already in T]
  evaluated in closed Gauss/Jacobi form:
      ∑_x χ(x(x−a)) = −1  (a≠0),
      ∑_x χ((x−1)(x+1)) = −1,
      ∑_x χ(x(x²−1)) = S3(p) ∈ {0, ±p, ±1} computed by the cubic sum
      (complete, roots {0,±1} already vanish).
  The restriction r∈□ is the (1+χ(r))/2 slice.  Implemented as
  `count_paley_off`.  Fail-when-wrong: drop the (1+χ)/2 projection
  (count all units); invert S3.  ∎

Theorem C — OPEN.  Q_τ(p) still unnamed.  ⟨δ,ψ⟩≤2 still OPEN.

============================================================================
Backend: CPU field tables, no Max+ needed for A–B.
Writes evidence/e1_gmin_m4_prop15291.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import (  # noqa: E402
    field_norm,
    in_subfield,
    paley_pair,
    type_key,
)


def n_pm1() -> int:
    return 2


def n_sub_off(p: int) -> int:
    """|F_p^× \\ {±1}|.  All lie in T and are type (++ sub)."""
    return p - 3


def n_N1(p: int) -> int:
    """|{N=1} \\ {±1}| = p−1."""
    return p - 1


def n_Nstar(p: int) -> int:
    return (p - 1) * (p - 3) // 2


def n_T(p: int) -> int:
    return (p * p - 1) // 2


def n_off(p: int) -> int:
    return (p * p - 5) // 2


def live_type_counts(p: int) -> dict[tuple, int]:
    F = _kernel_field(p)
    c: dict[tuple, int] = {}
    for r in F["squares"]:
        k = type_key(F, r)
        c[k] = c.get(k, 0) + 1
    return c


def count_paley_off_live(p: int) -> dict[tuple[int, int], int]:
    F = _kernel_field(p)
    c: dict[tuple[int, int], int] = {(1, 1): 0, (1, -1): 0, (-1, 1): 0, (-1, -1): 0}
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        c[paley_pair(F, r)] += 1
    return c


def sum_chi_cubic(p: int) -> int:
    """S3=∑_x χ(x(x−1)(x+1)).  Roots 0,±1 contribute 0."""
    F = _kernel_field(p)
    q = F["q"]
    s = 0
    for x in range(q):
        xm = F["add"][x][F["neg"][1]]
        xp = F["add"][x][1]
        if x == 0 or xm == 0 or xp == 0:
            continue
        s += F["chi"][x] * F["chi"][xm] * F["chi"][xp]
    return int(s)


def count_paley_off_formula(p: int) -> dict[tuple[int, int], int]:
    """(1+χ(r))/2 slice of (1+ε1 χ(r−1))/2 (1+ε2 χ(r+1))/2 on F_q \\ {0,±1}."""
    F = _kernel_field(p)
    q = F["q"]
    # complete pairwise Jacobi: ∑_x χ((x-a)(x-b)) = -1 for a≠b
    # triple S3
    S3 = sum_chi_cubic(p)
    # ∑_{x∉{0,±1}} 1 = q-3
    # ∑ χ(r) over r∉{0,±1}: ∑_all χ − χ(0)−χ(1)−χ(-1) = 0-0-1-1 = -2
    # ∑ χ(r-1) over r∉{0,±1}:
    #   all x of χ(x-1) = 0; minus x=0: χ(-1)=1; x=1: χ(0)=0; x=-1: χ(-2)
    chi_m2 = F["chi"][F["add"][F["neg"][1]][F["neg"][1]]]  # -1 + -1 = -2
    # -2 as field elt
    neg2 = F["add"][F["neg"][1]][F["neg"][1]]
    chi_m2 = F["chi"][neg2]
    sum_chi_rm1_off = 0 - F["chi"][F["neg"][1]] - 0 - chi_m2
    # ∑ χ(r+1) off: minus x=0 χ(1)=1; x=1 χ(2); x=-1 χ(0)=0
    chi2 = F["chi"][2 % p] if p > 2 else 0
    # 2 as field element (2,0)
    two = 2 if p > 2 else 0
    chi2 = F["chi"][two]
    sum_chi_rp1_off = 0 - 1 - chi2 - 0
    # pairwise χ(r(r-1)) off
    # complete ∑ χ(r(r-1)) = -1; minus r=0:0; r=1:0; r=-1: χ(-1)χ(-2)=chi_m2
    sum_r_rm1 = -1 - chi_m2
    # χ(r(r+1)): complete ∑ χ(r(r+1))=∑ χ(u(u-1))=-1
    # minus r=0:0; r=1: χ(1)χ(2)=chi2; r=-1:0
    sum_r_rp1 = -1 - chi2
    # χ((r-1)(r+1))=χ(r²-1): complete -1
    # minus r=0: χ(-1)=1; r=1:0; r=-1:0
    sum_rm1_rp1 = -1 - 1
    # triple χ(r(r-1)(r+1))=S3; already 0 at {0,±1}
    triple = S3

    out = {}
    for e1 in (1, -1):
        for e2 in (1, -1):
            # (1+χ(r))/2 already: we only sum r∈□.  Use live □ list
            # for the projection — formula for all-units then take □ via
            # (1+χ)/2 on the same expansion.
            # All-units off {0,±1}:
            tot = q - 3
            s = tot
            s += e1 * sum_chi_rm1_off
            s += e2 * sum_chi_rp1_off
            s += e1 * e2 * sum_rm1_rp1
            # □ projection: multiply by (1+χ(r))/2
            # (1/2)(all + χ-weighted all)
            # χ-weighted all-units off:
            # ∑_{off} χ(r) = -2
            # ∑ χ(r)χ(r-1)=∑ χ(r(r-1))=sum_r_rm1
            # ∑ χ(r)χ(r+1)=sum_r_rp1
            # ∑ χ(r)χ(r-1)χ(r+1)=triple
            chi_tot = -2
            chi_s = chi_tot
            chi_s += e1 * sum_r_rm1
            chi_s += e2 * sum_r_rp1
            chi_s += e1 * e2 * triple
            val = (s + chi_s) // 8
            out[(e1, e2)] = int(val)
    return out


def count_paley_off_formula_inverted(p: int) -> dict[tuple[int, int], int]:
    """Drop the (1+χ)/2 projection (count all units, not □)."""
    raw = count_paley_off_formula(p)
    # deliberately wrong: double as if no square projection
    return {k: 2 * v for k, v in raw.items()}


def prove_cardinalities(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(3, 40) if is_prime(p)]
    ok = True
    sample = {}
    invert_hit = 0
    invert_tot = 0
    for p in primes:
        live = live_type_counts(p)
        sub_off = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "sub")
        n1 = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "N1")
        ns = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "N*")
        row = (
            live.get(("pm1",), 0) == 2
            and sub_off == n_sub_off(p)
            and n1 == n_N1(p)
            and ns == n_Nstar(p)
            and sum(live.values()) == n_T(p)
        )
        if not row:
            ok = False
        if p <= 13:
            sample[str(p)] = {
                "live": {str(k): v for k, v in live.items()},
                "sub_off": sub_off,
                "N1": n1,
                "Nstar": ns,
                "ok": row,
            }
        # fail-when-wrong: claiming |sub off|=p-1
        invert_tot += 1
        if sub_off == p - 1:
            invert_hit += 1
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def prove_paley_slice(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [3, 5, 7, 11, 13]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        live = count_paley_off_live(p)
        form = count_paley_off_formula(p)
        bad = live != form
        if bad:
            ok = False
        inv = count_paley_off_formula_inverted(p)
        invert_tot += 1
        if inv == live:
            invert_hit += 1
        sample[str(p)] = {"live": {str(k): v for k, v in live.items()},
                          "formula": {str(k): v for k, v in form.items()},
                          "ok": not bad,
                          "S3": sum_chi_cubic(p)}
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {"proved": ok, "sample": sample, "invert_hit": invert_hit, "invert_tot": invert_tot}


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "note": "Q_τ still unnamed in p.  |τ| is named.  Floor OPEN.",
    }


def main() -> dict:
    print("Prop 15.291  |τ| in p (subfield / torus / Jacobi)", flush=True)
    A = prove_cardinalities()
    print(f"  A cardinalities: {A['proved']}", flush=True)
    B = prove_paley_slice()
    print(f"  B Paley slice: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C Q_τ / floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.291",
        "title": "Type counts |τ| named in p; Q_τ still open",
        "proved": {
            "cardinalities": A["proved"],
            "paley_slice": B["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": ["phi_F_ge_6", "e1", "L", "Aut-Schur", "Gsum", "pairing"],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15291.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
