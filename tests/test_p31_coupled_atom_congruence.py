from __future__ import annotations

from collections import Counter

import pytest

from e1_gmin_m4_p31_coupled_atom_congruence import (
    FIXED_DIRECTION,
    TARGET_PROFILE,
    compact_atom_vector,
    construct_coupled_ledger_graph,
    coupled_atom_congruence_certificate,
    signed_compact_lattice_decomposition,
)


def test_compact_lattice_constructor_replays_even_degree_vectors() -> None:
    target = {(0, 1): 3, (0, 2): -1, (1, 2): -1}
    atoms = signed_compact_lattice_decomposition(target, label_count=4)
    replay: Counter[tuple[int, int]] = Counter()
    for atom, multiplicity in atoms.items():
        for edge, coefficient in compact_atom_vector(atom).items():
            replay[edge] += multiplicity * coefficient
            if not replay[edge]:
                del replay[edge]
    assert dict(replay) == target
    assert sum(atoms.values()) == -sum(target.values())


def test_compact_lattice_constructor_rejects_odd_degree() -> None:
    with pytest.raises(ValueError, match="even signed degrees"):
        signed_compact_lattice_decomposition({(0, 1): 1}, label_count=3)


def test_constructed_graph_has_exact_hard_fixed_ledger() -> None:
    edges, centers, counts_before_fill = construct_coupled_ledger_graph()
    assert len(edges) == 479
    assert len(centers) == 16
    assert len(counts_before_fill) == 32
    assert max(counts_before_fill) <= min(TARGET_PROFILE)
    assert FIXED_DIRECTION == 1


def test_coupled_congruence_certificate_is_pinned() -> None:
    record = coupled_atom_congruence_certificate()
    assert record["proved"] is True
    assert record["parallel_profile"] == TARGET_PROFILE
    assert record["graph_sha256"] == (
        "36aea8d59a4131042de02a999a1f36070cc9d69150c19f17771543d45e46d116"
    )
    assert record["signed_atom_decompositions_sha256"] == (
        "42745425e873264598bbec521a83e744db1db81ef0d80a39a32bfcecb144bee7"
    )
    assert record["row_count"] == 32
    assert record["all_rows_in_required_signed_integer_atom_lattice"] is True
    assert record["pure_lattice_congruence_cut_found"] is False
    assert record["nonnegative_atom_decomposition_constructed"] is False
    assert record["minimum_compact_coefficient"] == -17
    assert record["maximum_compact_coefficient"] == 5
    assert record["total_signed_compact_support"] == 16_048
    assert record["total_signed_compact_l1"] == 25_451
