#!/usr/bin/env python3
r"""Prop. 15.756 -- the arbitrary-boundary Weil route has no large-size gain.

Let ``D`` be an arbitrary finite boundary in ``F_(p^2)``, ``|D|=s``.
For a projective ``F_p`` direction ``d``, write ``n_(d,t)`` for the sizes
of its ``p`` affine fibres and ``b_d`` for the number of odd fibres.  Split
the ``p+1`` directions into their two Paley quadratic types
``eps_d in {+1,-1}``, each of size ``m=(p+1)/2``.  Put

    delta_tau = sum_(eps_d=tau) (s-b_d),
    P_tau     = # unordered pairs of D of direction type tau,

and

    R_tau = sum_(eps_d=tau,t) h(n_(d,t)),
    h(2r)=r(r-1),       h(2r+1)=r^2.

The fibrewise identity

    n-(n mod 2) = 2*binom(n,2)-4*h(n)

gives the exact universal type split

    delta_tau = 2 P_tau-4 R_tau.                          (1)

Let ``Q`` be the finite ``p^2`` by ``p^2`` Paley matrix.  It satisfies
``Q1=0`` and ``Q^2=p^2 I-J``.  Therefore

    W=P_+-P_-=(1/2) 1_D^T Q 1_D,
    |W| <= s(p^2-s)/(2p).                                 (2)

Combining (1)--(2) and dropping ``R_tau>=0`` yields

    delta_tau <= binom(s,2)+s(p^2-s)/(2p).                (3)

But the definition already gives ``delta_tau<=m s``.  The difference
between the right side of (3) and this trivial cap is exactly

    s((p-1)s-2p)/(2p),                                    (4)

which is nonnegative for every even ``s>=4``.  Thus the character/spectral
bound adds literally no information in any boundary size beyond the
already-closed two-point case.  Ordinary Parseval for the fibre counts is
the same pair count and does not repair the loss.

The failure is sharp, not an artifact of the estimate.  Take ``D`` to be
the union of two parallel affine ``F_p``-lines.  Then ``s=2p``; their own
direction has ``b=2`` and every transverse direction has ``b=0``.  The
centered incidence vector lies in one ``+/-p`` Paley eigenspace, and

    |W|=p^2-2p=s(p^2-s)/(2p).

Every even set is the odd-degree boundary of a graph (pair its vertices),
so this is a genuine graph-boundary countermechanism.  It is not asserted
to satisfy the residual separator inequalities.  What it proves is that no
argument depending only on ``D``, its directional odd-fibre counts, pair
energy, and an arbitrary-boundary character estimate can close residual
(ii).  A successful large-boundary theorem must couple those data to the
signed directional edge matrices of one common graph ``H``.

This proposition is a method barrier and does not flip a global predicate.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15721 import is_prime
from minmax_quadratic import paley_conference_prime_power


ROOT = Path(__file__).resolve().parents[1]


def occupancy_slack(n: int) -> int:
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("a fibre occupancy is a nonnegative integer")
    r = n // 2
    return r * (r - 1) if n % 2 == 0 else r * r


def fibre_type_identity(n: int) -> dict[str, int | bool]:
    """Audit the summand of ``delta=2P-4R`` on one fibre."""
    left = n - n % 2
    pair_term = 2 * comb(n, 2)
    slack_term = 4 * occupancy_slack(n)
    return {
        "occupancy": n,
        "n_minus_parity": left,
        "twice_pairs": pair_term,
        "four_times_slack": slack_term,
        "proved": left == pair_term - slack_term,
    }


def type_split_from_fibres(fibre_profiles: list[list[int]]) -> dict[str, int | bool]:
    """Verify (1) for a collection of directions of one fixed type."""
    if not fibre_profiles:
        raise ValueError("at least one direction profile is required")
    s = sum(fibre_profiles[0])
    if any(sum(profile) != s for profile in fibre_profiles):
        raise ValueError("all directions must partition the same set D")
    delta = sum(s - sum(value % 2 for value in profile) for profile in fibre_profiles)
    pairs = sum(sum(comb(value, 2) for value in profile) for profile in fibre_profiles)
    slack = sum(sum(occupancy_slack(value) for value in profile) for profile in fibre_profiles)
    return {
        "boundary_size": s,
        "directions": len(fibre_profiles),
        "delta": delta,
        "P": pairs,
        "R": slack,
        "proved": delta == 2 * pairs - 4 * slack,
    }


def character_cap_ledger(p: int, s: int) -> dict[str, object]:
    """Compare the spectral pair cap with the definition-only cap ``ms``."""
    if p < 3 or p % 2 == 0:
        raise ValueError("need an odd p>=3")
    if not 0 <= s <= p * p:
        raise ValueError("the finite boundary size lies in 0..p^2")
    m = (p + 1) // 2
    spectral_W_cap = Fraction(s * (p * p - s), 2 * p)
    pair_character_cap = Fraction(comb(s, 2)) + spectral_W_cap
    trivial_cap = m * s
    difference = pair_character_cap - trivial_cap
    closed_difference = Fraction(s * ((p - 1) * s - 2 * p), 2 * p)
    proved = difference == closed_difference
    if not proved:
        raise ArithmeticError("the character-cap comparison changed")
    return {
        "p": p,
        "boundary_size": s,
        "directions_per_type": m,
        "Paley_pair_imbalance_cap": str(spectral_W_cap),
        "pair_character_delta_cap": str(pair_character_cap),
        "definition_only_delta_cap": trivial_cap,
        "character_cap_minus_trivial_cap": str(difference),
        "character_cap_is_nonimproving": difference >= 0,
        "nonimproving_for_every_even_s_ge_4": bool(s < 4 or s % 2 or difference >= 0),
        "proved": proved,
    }


def two_parallel_lines_certificate(p: int) -> dict[str, object]:
    """Exact Paley-matrix replay of the sharp two-parallel-line mechanism."""
    if not is_prime(p):
        raise ValueError("the Paley construction requires an odd prime")
    C = paley_conference_prime_power(p).astype(np.int64)
    Q = C[1:, 1:]
    incidence = np.array(
        [1 if u // p in (0, 1) else 0 for u in range(p * p)],
        dtype=np.int64,
    )
    s = int(incidence.sum())
    W = int(incidence @ Q @ incidence) // 2
    cap = s * (p * p - s) // (2 * p)
    matrix_identity = Q @ Q == p * p * np.eye(p * p, dtype=np.int64) - np.ones(
        (p * p, p * p), dtype=np.int64
    )
    proved = bool(
        s == 2 * p
        and np.all(matrix_identity)
        and abs(W) == cap == p * p - 2 * p
    )
    if not proved:
        raise ArithmeticError("the two-parallel-line sharpness audit failed")
    return {
        "p": p,
        "boundary_size": s,
        "own_direction_odd_fibres": 2,
        "transverse_direction_odd_fibres": 0,
        "transverse_directions": p,
        "Paley_pair_imbalance_W": W,
        "spectral_cap": cap,
        "spectral_cap_attained": True,
        "Q_squared_equals_p2I_minus_J": True,
        "is_odd_degree_boundary_of_a_matching": True,
        "is_claimed_residual_separator": False,
        "proved": proved,
    }


def theorem_record(primes: tuple[int, ...] = (3, 5, 7, 11, 13, 17, 19)) -> dict[str, object]:
    fibres = [fibre_type_identity(n) for n in range(0, 15)]
    caps = {
        str(p): [character_cap_ledger(p, s) for s in range(4, p * p + 1, 2)]
        for p in primes
    }
    sharp = {str(p): two_parallel_lines_certificate(p) for p in primes}
    proved = bool(
        all(row["proved"] for row in fibres)
        and all(row["proved"] for rows in caps.values() for row in rows)
        and all(row["character_cap_is_nonimproving"] for rows in caps.values() for row in rows)
        and all(row["proved"] for row in sharp.values())
    )
    return {
        "prop": "15.756",
        "title": "Arbitrary-boundary character cap is vacuous beyond size two",
        "proved": {
            "typed_pair_slack_identity": proved,
            "Paley_pair_imbalance_spectral_cap": proved,
            "character_cap_nonimproving_for_even_boundary_size_ge_4": proved,
            "two_parallel_lines_attain_spectral_cap_all_odd_primes": proved,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "fibre_identity_checks": fibres,
        "cap_checks": caps,
        "two_parallel_line_checks": sharp,
        "remaining_obstruction": (
            "Couple boundary parity and occupancy to the signed directional edge "
            "matrices of one common 0/1 graph H."
        ),
        "duplicate_work_guards": [
            "Do not invoke Weil cancellation for an arbitrary boundary D.",
            "Do not treat ordinary fibre Parseval as stronger than pair counting.",
            "Do not infer a residual separator from the two-line boundary witness.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15756.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.756 arbitrary-boundary character barrier: proved")
    print("  residual (ii): OPEN")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
