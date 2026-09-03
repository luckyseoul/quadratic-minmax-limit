from fractions import Fraction

from e1_gmin_m4_prop15756 import (
    character_cap_ledger,
    fibre_type_identity,
    theorem_record,
    two_parallel_lines_certificate,
    type_split_from_fibres,
)


def test_fibre_identity_and_type_sum():
    assert all(fibre_type_identity(n)["proved"] for n in range(30))
    row = type_split_from_fibres([[0, 1, 2, 3, 4], [2, 2, 2, 2, 2]])
    assert row["proved"]


def test_character_cap_is_vacuous_for_even_size_at_least_four():
    for p in (3, 5, 7, 11, 13):
        for s in range(4, p * p + 1, 2):
            row = character_cap_ledger(p, s)
            assert row["proved"]
            assert row["character_cap_is_nonimproving"]
            assert Fraction(row["character_cap_minus_trivial_cap"]) >= 0


def test_two_parallel_lines_are_sharp():
    for p in (3, 5, 7, 11, 13):
        row = two_parallel_lines_certificate(p)
        assert row["proved"]
        assert row["spectral_cap_attained"]
        assert not row["is_claimed_residual_separator"]


def test_theorem_record_keeps_residual_open():
    out = theorem_record((3, 5, 7, 11))
    assert out["proved"]["character_cap_nonimproving_for_even_boundary_size_ge_4"]
    assert not out["proved"]["residual_ii_closed"]
    assert out["L_status"] == "OPEN"
