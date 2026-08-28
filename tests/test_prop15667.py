import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15667 import (
    ARTIFACT_SHA256,
    CRT_PRODUCT,
    ENDPOINTS_E120,
    ENDPOINTS_E800,
    MAX_COUNT_BOUND,
    MAX_FOURTH_MOMENT_BOUND,
    MAX_SECOND_MOMENT_BOUND,
    P11_EARLY_PROFILE,
    TRACE_PIVOT_EXPONENTS,
    endpoint_plateau_audit,
    raw_trace_from_profile_moment,
    theorem_record,
    tight_frame_second_moment_holds,
)


ROOT = Path(__file__).resolve().parents[1]


def test_prop15667_profile_crt_and_early_moment_identities():
    assert CRT_PRODUCT > max(
        MAX_COUNT_BOUND, MAX_SECOND_MOMENT_BOUND, MAX_FOURTH_MOMENT_BOUND
    )
    expected_raw = {
        11: Fraction(0),
        20: Fraction(89792, 11),
        24: Fraction(7076),
        27: Fraction(538752),
    }
    for exponent, count, moment2, moment4, _harmonic, raw in P11_EARLY_PROFILE:
        assert tight_frame_second_moment_holds(11, exponent, count, moment2)
        assert raw_trace_from_profile_moment(11, exponent, count, moment4) == Fraction(raw)
        assert Fraction(raw) == expected_raw[exponent]


def test_prop15667_modular_ranks_and_honest_scope():
    row = theorem_record()
    proved = row["proved"]
    assert row["finite_profile"]["glue_dual_codewords"] == 11**10
    assert row["finite_profile"]["translation_scalar_representatives"] == 21_437_340
    assert row["finite_profile"]["weighted_histogram_tuples"] == 2_558_543
    assert row["finite_profile"]["nonempty_tight_frame_checks"] == 51
    modular = row["modular_reconstruction"]
    assert modular["scalar_space_dimension"] == modular["scalar_prefix_rank"] == 41
    assert modular["trace_affine_dimension"] == modular["trace_prefix_rank"] == 32
    assert tuple(modular["trace_pivot_exponents"]) == TRACE_PIVOT_EXPONENTS
    assert modular["scalar_held_out_coefficients_matched"] == 32
    assert modular["trace_held_out_coefficients_matched"] == 28
    assert proved["aggregate_trace_conservation_closes_R1"] is False
    assert proved["R1"] is False
    assert proved["global_QVAR"] is False
    assert proved["type_I_multilevel"] is False
    assert proved["limit_exists"] is False
    assert row["L_status"] == "OPEN"


def test_prop15667_endpoint_plateau_is_exactly_seven_of_eight():
    audit = endpoint_plateau_audit()
    assert audit["equal_endpoint_count"] == 7
    assert audit["total_endpoint_count"] == 8
    assert audit["sole_late_contraction"] == "circle-low-Weil maximum"
    assert ENDPOINTS_E120["circle-low-Weil"]["minimum"] == ENDPOINTS_E800[
        "circle-low-Weil"
    ]["minimum"]
    assert Fraction(ENDPOINTS_E800["circle-low-Weil"]["maximum"]) < Fraction(
        ENDPOINTS_E120["circle-low-Weil"]["maximum"]
    )


def test_prop15667_artifact_hashes_and_committed_record_are_pinned():
    assert ARTIFACT_SHA256
    assert all(len(value) == 64 for value in ARTIFACT_SHA256.values())
    committed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15667.json").read_text()
    )
    assert committed == theorem_record()
