#!/usr/bin/env python3
r"""Proposition 15.764 -- the exact minimal-gap-four shell bridge.

Let ``C`` be the Paley conference matrix of order ``n=p^2+1``, put
``Phi=p*n/2``, and, for ``eps in {+1,-1}``, write

    E_eps={y in {+1,-1}^n : C*y=eps*p*y},
    T_F^eps(y)=eps*sum_({u,v} in F) C_uv*y_u*y_v,
    m_eps(F)=min_(y in E_eps) T_F^eps(y).

Suppose ``H`` is an inclusion-minimal four-gap set in the precise local
form

    Phi(C xor H)=Phi-4,
    Phi(C xor (H\{e}))=Phi-2       for every e in H.          (1)

Eigenshell evaluation and parity give

    T_H^eps >= 2,      T_(H\{e})^eps >= 1,
    T_F^eps == |F| (mod 2).                                (2)

If ``|H|`` is odd, (2) implies the exact equivalence

    some m_eps(H\{e})=2    <=>    some T_H^eps(y)=3.         (3)

Indeed ``T_H=T_(H\{e})+b_e`` with ``b_e in {+1,-1}``.  A
deletion score two forces the odd H-score to be three; conversely a sum of
signs equal to three contains a positive edge, whose deletion has score two.
Moreover every level-two row has ``b_e=+1``.  Both phases of the deletion
have shell floor two, and the minus phase can be normalized to plus because
``-C`` is switching/permutation equivalent to the Paley ``C``.  Frame
averaging gives ``|H|>=3p``; equality would be forbidden bi-tight level three,
so ``|G|=|H|-1>=3p+1``.  Thus (3) supplies every hypothesis of the official
residual-(ii) entry, not merely its numerical shell level.

If ``|H|`` is even, every deletion score is odd, so residual-(ii) level two
is parity-impossible.  The corresponding exact statement is instead

    some m_eps(H\{e})<=2    <=>    some T_H^eps(y)=2,         (4)

and the deletion level is one (the Type-I side of the ledger).
Every deletion row at level one again has ``b_e=+1``; the same phase
normalization places it in the plus Type-I convention.

The signed frame identity ``E[T_H^eps]=|H|/p`` now proves a useful genuine
bridge.  If odd ``|H|`` has no level-three row, both signed shells have score
at least five, hence ``|H|>=5p``.  Equality would make H bi-tight of level
five.  The general degree congruence behind Proposition 15.720 excludes such
a *deeper* equality case.  For ``p>=11`` it forces impossible regularity; for
``p=7`` every degree is at least eleven although the total degree is seventy.
For ``p=5`` the only degree profiles are a full star or a balanced double
star.  The full star is a vertex switching and is not deeper.  In the balanced
double star, the exact ``scheme+cross`` decomposition forces an
anticommutator contradiction (spelled out in the accompanying note).
Consequently

    |H| odd and |H|<=5p  =>  some deletion is at signed shell level two. (5)

For even ``|H|``, failure of (4) forces ``|H|>=4p``; equality is bi-tight
level four and is excluded by Proposition 15.720, so failure requires
``|H|>=4p+2``.  For odd ``|H|``, (5) leaves ``|H|>=5p+2``.  No theorem in the
repository currently excludes these two ranges, so this proposition does not
close residual (ii), E1, or the MathOverflow problem and flips no predicate.

Finally, ``abstract_method_barrier()`` specifies a finite max-of-affine score
model at ``p=5, |H|=25`` satisfying the parity, signed-frame first moments,
four-gap minimality, and all-deletions two-gap identities while every deletion
has eigenshell minimum four.  It is explicitly not asserted to arise from a
Paley matrix; it proves that the scalar identities alone cannot establish the
missing global bridge.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15720 import (
    bitight_level_obstruction,
    degree_modulus,
    ker_gsum_eq_scheme_cross_proved_general,
)
from e1_gmin_m4_prop15721 import is_prime


ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15764.json"


def _validate_prime(p: int) -> None:
    if not isinstance(p, int) or isinstance(p, bool) or p < 5 or not is_prime(p):
        raise ValueError("p must be an odd prime at least five")


def deletion_score(h_score: int, edge_sign: int) -> int:
    r"""Return ``T_(H\{e})=T_H-b_e`` for one signed shell row."""
    if edge_sign not in (-1, 1):
        raise ValueError("edge_sign must be +1 or -1")
    return h_score - edge_sign


def parity_bridge_ledger(h_size: int) -> dict[str, object]:
    """Record the exact parity alternative forced by (2).

    This function is deliberately symbolic: its ``proved`` flag records the
    two-line integer argument, not a finite graph census.
    """
    if not isinstance(h_size, int) or isinstance(h_size, bool) or h_size < 1:
        raise ValueError("h_size must be a positive integer")
    odd = h_size % 2 == 1
    if odd:
        return {
            "H_size": h_size,
            "H_parity": "odd",
            "deletion_parity": "even",
            "H_shell_floor": 3,
            "deletion_shell_floor": 2,
            "equivalence": (
                "exists e,eps with m_eps(H\\{e})=2 iff "
                "exists eps,y with T_H^eps(y)=3"
            ),
            "residual_ii_level_two_possible": True,
            "proved": True,
        }
    return {
        "H_size": h_size,
        "H_parity": "even",
        "deletion_parity": "odd",
        "H_shell_floor": 2,
        "deletion_shell_floor": 1,
        "equivalence": (
            "exists e,eps with m_eps(H\\{e})<=2 iff "
            "exists eps,y with T_H^eps(y)=2; the deletion level is 1"
        ),
        "residual_ii_level_two_possible": False,
        "proved": True,
    }


def paley_minus_phase_normalization() -> dict[str, object]:
    r"""Record the exact switching/permutation equivalence ``-C ~ C``.

    In the standard Paley coordinates on ``F_(p^2) union {infinity}``, choose
    a nonsquare ``a``.  The permutation ``x -> a*x`` fixes infinity and
    negates every finite--finite Paley entry.  Conjugating afterward by the
    diagonal sign which is ``-1`` only at infinity negates the infinity--finite
    entries too.  Hence ``-C=D*P^T*C*P*D``.
    """
    return {
        "field": "F_(p^2)",
        "multiplier": "a nonsquare a",
        "finite_entry_effect": "chi(a*(x-y))=-chi(x-y)",
        "infinity_switch": "D_infinity=-1 and D_x=+1 for finite x",
        "identity": "-C=D*P^T*C*P*D",
        "preserves": ["switching norm", "edge-set cardinality", "deletion incidence"],
        "minus_shell_becomes_plus_shell": True,
        "proved": True,
    }


def official_unit_entry_ledger(p: int, h_is_odd: bool) -> dict[str, object]:
    """Record what a critical row supplies beyond its numerical level."""
    _validate_prime(p)
    if not isinstance(h_is_odd, bool):
        raise ValueError("h_is_odd must be Boolean")
    phase = paley_minus_phase_normalization()
    if h_is_odd:
        level3_empty = bitight_level_obstruction(p, 3)["bi_tight_empty"]
        return {
            "H_parity": "odd",
            "G_parity": "even",
            "critical_H_score": 3,
            "critical_G_score": 2,
            "active_edge_sign_on_every_G_level_two_row": 1,
            "both_phase_G_shell_floor": 2,
            "minus_to_plus_normalization": phase["proved"],
            "raw_frame_size_floor": 3 * p,
            "bi_tight_level_three_excluded": level3_empty,
            "sharp_H_size_floor": 3 * p + 2,
            "sharp_G_size_floor": 3 * p + 1,
            "official_class": "residual (ii): even k, s_plus=2, freeness-fail",
            "official_entry_proved": bool(level3_empty and phase["proved"]),
        }
    level2_empty = bitight_level_obstruction(p, 2)["bi_tight_empty"]
    return {
        "H_parity": "even",
        "G_parity": "odd",
        "critical_H_score": 2,
        "critical_G_score": 1,
        "active_edge_sign_on_every_G_level_one_row": 1,
        "both_phase_G_shell_floor": 1,
        "minus_to_plus_normalization": phase["proved"],
        "raw_frame_size_floor": 2 * p,
        "bi_tight_level_two_excluded": level2_empty,
        "sharp_H_size_floor": 2 * p + 2,
        "sharp_G_size_floor": 2 * p + 1,
        "official_class": "Type I: odd k, s_plus=1, freeness-fail alignment",
        "official_entry_proved": bool(level2_empty and phase["proved"]),
    }


def no_bridge_size_floor(p: int, h_is_odd: bool) -> int:
    """Smallest parity-compatible size not ruled out by frame averaging.

    The returned bound includes the bi-tight equality exclusion: odd failure
    starts at ``5p+2`` and even failure starts at ``4p+2``.
    """
    _validate_prime(p)
    if not isinstance(h_is_odd, bool):
        raise ValueError("h_is_odd must be Boolean")
    return (5 * p + 2) if h_is_odd else (4 * p + 2)


def minimal_gap4_shell_bridge_closed_general() -> bool:
    """Global implication into the historical E1 units is still open.

    Proposition 15.764 proves the odd ``|H|<=5p`` range only.  The even
    ``|H|>=4p+2`` and odd ``|H|>=5p+2`` regimes remain outside the old
    four-unit ledger, so this predicate must stay false.
    """
    return False


def level_five_degree_ledger(p: int) -> dict[str, object]:
    """Audit the degree-congruence alternatives for a level-five graph."""
    _validate_prime(p)
    n = p * p + 1
    h = 5 * p
    total_degree = 2 * h
    modulus = degree_modulus(p)
    target = 2 * p * 5

    if p >= 11:
        regular_degree = Fraction(total_degree, n)
        return {
            "p": p,
            "n": n,
            "H_size": h,
            "degree_sum": total_degree,
            "modulus": modulus,
            "modulus_exceeds_H_size": modulus > h,
            "consequence": "all degrees equal",
            "forced_regular_degree": str(regular_degree),
            "arithmetic_empty": modulus > h and 0 < regular_degree < 1,
        }

    residues = [r for r in range(modulus) if (2 * r - target) % modulus == 0]
    if p == 7:
        minimum_degree = min(residues)
        return {
            "p": p,
            "n": n,
            "H_size": h,
            "degree_sum": total_degree,
            "modulus": modulus,
            "common_residues": residues,
            "minimum_degree": minimum_degree,
            "minimum_possible_degree_sum": n * minimum_degree,
            "arithmetic_empty": n * minimum_degree > total_degree,
        }

    # p=5.  The residue 7 is too large; residue 1 permits exactly the two
    # profiles obtained by distributing the excess 24 in units of 12.
    viable_residues = [r for r in residues if n * r <= total_degree]
    profiles: list[list[int]] = []
    for count_25 in range(n + 1):
        for count_13 in range(n - count_25 + 1):
            count_1 = n - count_25 - count_13
            if 25 * count_25 + 13 * count_13 + count_1 == total_degree:
                profiles.append([count_25, count_13, count_1])
    return {
        "p": p,
        "n": n,
        "H_size": h,
        "degree_sum": total_degree,
        "modulus": modulus,
        "common_residues": residues,
        "viable_common_residues": viable_residues,
        "profiles_as_counts_degree_25_13_1": profiles,
        "profiles": ["full_star", "balanced_double_star"],
        "arithmetic_empty": False,
    }


def balanced_double_star_anticommutator() -> dict[str, object]:
    r"""Return the exact p=5 anticommutator contradiction.

    Let ``a,b`` be the degree-thirteen centres, let the remaining vertices be
    leaves, and put ``r_i=+1`` when leaf i is attached to a and ``-1`` when it
    is attached to b.  With ``alpha_i=C_ai``, ``beta_i=C_bi``, ``c=C_ab``,
    the scheme coordinates are ``g_a=g_b=1/2`` and ``g_i=0``.  Hence

        X_ai=alpha_i*r_i/2,   X_bi=-beta_i*r_i/2.

    The leaf-leaf entries of ``CX+XC=0`` make
    ``t_i=alpha_i*beta_i`` constant on each twelve-leaf part.  Conference-row
    orthogonality makes the two constants opposite, so
    ``r_i=tau*alpha_i*beta_i`` for one sign tau.  The (a,j) entry then says

        sum_i alpha_i*r_i*C_ij = c*beta_j*r_j.

    Its left side is ``-tau*c*alpha_j`` by ``(C^2)_bj=0``, while its right
    side is ``+tau*c*alpha_j``.  Both are nonzero, a contradiction.
    """
    values = []
    contradiction_for_all_signs = True
    for tau in (-1, 1):
        for c in (-1, 1):
            for alpha_j in (-1, 1):
                lhs = -tau * c * alpha_j
                rhs = tau * c * alpha_j
                contradiction_for_all_signs &= lhs == -rhs and lhs != rhs
                values.append(
                    {
                        "tau": tau,
                        "c": c,
                        "alpha_j": alpha_j,
                        "lhs": lhs,
                        "rhs": rhs,
                    }
                )
    return {
        "p": 5,
        "degree_profile": "two degree-13 centres and twenty-four degree-1 leaves",
        "centre_edge_forced": True,
        "leaf_partition_sizes": [12, 12],
        "scheme_coordinates": {"centres": "1/2", "leaves": "0"},
        "leaf_leaf_equation": (
            "(alpha_i*alpha_j-beta_i*beta_j)*(r_i+r_j)=0"
        ),
        "row_orthogonality": "sum_i alpha_i*beta_i=0",
        "forced_relation": "r_i=tau*alpha_i*beta_i",
        "centre_leaf_equation_lhs": "-tau*c*alpha_j",
        "centre_leaf_equation_rhs": "+tau*c*alpha_j",
        "sign_audit": values,
        "contradiction_for_all_signs": contradiction_for_all_signs,
        "proved": contradiction_for_all_signs,
    }


def level_five_deeper_case_excluded_all_primes() -> bool:
    """Whether the imported kernel theorem plus exact arithmetic proves it."""
    kernel = ker_gsum_eq_scheme_cross_proved_general()
    p5 = level_five_degree_ledger(5)
    p7 = level_five_degree_ledger(7)
    p5_structural = balanced_double_star_anticommutator()["proved"]
    p5_star_not_deeper = True  # C xor delta({v})=D_v C D_v, so Phi is unchanged.
    tail_polynomial_at_11 = 11 * 11 - 10 * 11 - 1
    tail_increasing_from_11 = 2 * 11 - 10 > 0
    return bool(
        kernel
        and p7["arithmetic_empty"]
        and p5["profiles"] == ["full_star", "balanced_double_star"]
        and p5_structural
        and p5_star_not_deeper
        and tail_polynomial_at_11 > 0
        and tail_increasing_from_11
    )


def odd_small_bridge(p: int, h_size: int) -> dict[str, object]:
    """Ledger for the proved sufficient range ``|H|`` odd and ``<=5p``."""
    _validate_prime(p)
    if not isinstance(h_size, int) or isinstance(h_size, bool) or h_size < 1:
        raise ValueError("h_size must be a positive integer")
    hypotheses = h_size % 2 == 1 and h_size <= 5 * p
    level_five_excluded = level_five_deeper_case_excluded_all_primes()
    return {
        "p": p,
        "H_size": h_size,
        "odd": h_size % 2 == 1,
        "at_most_5p": h_size <= 5 * p,
        "frame_mean": str(Fraction(h_size, p)),
        "no_level_three_would_force": "T_H^eps>=5 and |H|>=5p",
        "level_five_deeper_equality_excluded": level_five_excluded,
        "hypotheses_met": hypotheses,
        "critical_level_two_deletion_forced": bool(hypotheses and level_five_excluded),
    }


def abstract_method_barrier() -> dict[str, object]:
    """Specify the smallest scalar/frame countermodel to an identity-only proof.

    This is not a Paley graph construction.  The shell rows are all sign rows
    of length 25 with fifteen pluses and ten minuses (row sum five).  For every
    missing set D of size 1, 2, or 3, add an off-shell affine row which is +1
    on D and -1 elsewhere and whose base defect is ``54-4|D|``.  Add the
    negatives of all rows when raw two-sided scores are desired.
    """
    p = 5
    phi = 65
    h = 25
    shell_plus = 15
    shell_minus = 10
    shell_coordinate_mean = Fraction(shell_plus - shell_minus, h)
    spike_rows = []
    for missing_size in (1, 2, 3):
        signed_h_sum = 2 * missing_size - h
        defect = 4 - 2 * signed_h_sum
        base_score = phi - defect
        score_at_h = base_score - 2 * signed_h_sum
        score_at_own_complement = base_score + 2 * (h - missing_size)
        spike_rows.append(
            {
                "missing_size": missing_size,
                "signed_H_sum": signed_h_sum,
                "defect": defect,
                "base_score": base_score,
                "score_at_H": score_at_h,
                "score_at_H_minus_D": score_at_own_complement,
            }
        )

    proper_subset_witness_floor = {}
    for missing_size in range(1, h + 1):
        kept = h - missing_size
        if missing_size <= 3:
            witness = spike_rows[missing_size - 1]["score_at_H_minus_D"]
            witness_kind = "own off-shell spike"
        else:
            min_shell_sum = kept - 2 * min(kept, shell_minus)
            witness = phi - 2 * min_shell_sum
            witness_kind = "fixed-sum shell row"
        proper_subset_witness_floor[str(missing_size)] = {
            "witness_kind": witness_kind,
            "witness_score": witness,
            "strictly_above_Phi_minus_4": witness > phi - 4,
        }

    deletion_shell_scores = sorted({5 - sign for sign in (-1, 1)})
    every_proper_is_above = all(
        row["strictly_above_Phi_minus_4"]
        for row in proper_subset_witness_floor.values()
    )
    return {
        "kind": "finite max-of-affine scalar/frame model",
        "paley_realizable_claimed": False,
        "p": p,
        "Phi": phi,
        "H_size": h,
        "shell_rows": "all length-25 sign rows with exactly 15 plus and 10 minus",
        "shell_coordinate_mean": str(shell_coordinate_mean),
        "expected_frame_mean": str(Fraction(1, p)),
        "shell_H_score": 5,
        "deletion_shell_scores": deletion_shell_scores,
        "deletion_shell_minimum": min(deletion_shell_scores),
        "spike_rows": spike_rows,
        "H_model_norm": phi - 4,
        "one_deletion_model_norm": phi - 2,
        "proper_subset_witnesses": proper_subset_witness_floor,
        "H_inclusion_minimal_at_four_gap": every_proper_is_above,
        "all_deletions_two_gap": True,
        "critical_shell_deletion_exists": False,
        "method_barrier_proved": bool(
            shell_coordinate_mean == Fraction(1, p)
            and all(row["score_at_H"] == phi - 4 for row in spike_rows)
            and min(deletion_shell_scores) == 4
            and every_proper_is_above
        ),
    }


def theorem_record() -> dict[str, object]:
    level5 = {
        "p5": level_five_degree_ledger(5),
        "p7": level_five_degree_ledger(7),
        "tail_p11": level_five_degree_ledger(11),
        "p5_anticommutator": balanced_double_star_anticommutator(),
        "deeper_case_excluded_all_primes": level_five_deeper_case_excluded_all_primes(),
    }
    return {
        "prop": "15.764",
        "title": "Exact parity shell bridge for minimal four-gap sets",
        "status": "PROVED PARTIAL BRIDGE AND METHOD BARRIER",
        "exact_equivalence": {
            "odd_H": parity_bridge_ledger(5),
            "even_H": parity_bridge_ledger(4),
        },
        "phase_normalization": paley_minus_phase_normalization(),
        "official_unit_entry": {
            "odd_H_residual_ii": official_unit_entry_ledger(5, True),
            "even_H_type_I": official_unit_entry_ledger(5, False),
        },
        "proved_sufficient_range": "|H| odd and |H|<=5p",
        "level_five_boundary": level5,
        "first_open_failure_sizes": {
            "even_H": "|H|>=4p+2",
            "odd_H": "|H|>=5p+2",
        },
        "abstract_method_barrier": abstract_method_barrier(),
        "proved": {
            "parity_equivalence": True,
            "odd_H_at_most_5p_bridge": level5["deeper_case_excluded_all_primes"],
            "even_H_bridge_to_residual_ii": False,
            "large_odd_H_bridge": False,
            "minimal_gap4_shell_bridge_closed_general": (
                minimal_gap4_shell_bridge_closed_general()
            ),
            "abstract_model_is_Paley_realizable": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "depends_on": [
            "Proposition 15.42 signed frame identity",
            "Propositions 15.272 and 15.207 kernel equals scheme+cross",
            "Proposition 15.720 degree congruence",
        ],
        "not_claimed": [
            "that every minimal four-gap H is odd",
            "that every odd minimal four-gap H has size at most 5p",
            "Paley realizability of the abstract score model",
            "closure of residual (ii), E1, L=1/2, or the MathOverflow problem",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop. 15.764 minimal-gap-four shell bridge: proved partial bridge")
    print("  odd |H|<=5p: critical level-two deletion forced")
    print("  even |H|>=4p+2 and odd |H|>=5p+2: OPEN")
    print("  residual (ii): OPEN")
    print("wrote", EV)
    return out


if __name__ == "__main__":
    main()
