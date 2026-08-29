from fractions import Fraction

from e1_gmin_m4_prop15690 import (
    dilation_energy_identity,
    p7_psd_autocorrelation_obstruction,
    theorem_record,
)


def test_dilation_energy_is_exactly_strong_r1_at_threshold():
    row = dilation_energy_identity(49, Fraction(1, 12))
    assert row["V_over_n"] == 2
    assert row["S_K"] == 48
    assert row["strong_R1_equivalence"] is True


def test_psd_autocorrelation_relaxation_has_explicit_violation():
    row = p7_psd_autocorrelation_obstruction()
    assert row["violation_factor"] == Fraction(23_040, 11)
    assert row["is_actual_uniform_full_Max_plus_ensemble"] is False
    assert theorem_record()["closes_R1"] is False
