import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15668 import (
    AGGREGATE_TRACE_TARGET,
    ARTIFACT_SHA256,
    BROAD_DIMENSIONS,
    BROAD_TARGETS,
    CRT_PRODUCT,
    ENDPOINTS_E800,
    MAX_COUNT_BOUND,
    MAX_PROFILE_LEGENDRE_FOURTH_BOUND,
    MAX_SQUARED_EXCESS_BOUND,
    P11_DELTA_SQ,
    PIVOT_EXPONENTS,
    finite_p11_r1_audit,
    phi_interval,
    phi_scale,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_prop15668_marked_crt_and_broad_target_conservation():
    assert CRT_PRODUCT > max(
        MAX_COUNT_BOUND,
        MAX_SQUARED_EXCESS_BOUND,
        MAX_PROFILE_LEGENDRE_FOURTH_BOUND,
    )
    assert sum(BROAD_DIMENSIONS.values()) == 1769
    assert sum(
        (
            BROAD_DIMENSIONS[channel] * BROAD_TARGETS[channel]
            for channel in BROAD_DIMENSIONS
        ),
        Fraction(),
    ) == AGGREGATE_TRACE_TARGET


def test_prop15668_rank_holdouts_and_exact_endpoint_route_limit():
    row = theorem_record()
    modular = row["modular_reconstruction"]
    assert modular["affine_dimension_per_channel"] == 32
    assert modular["first_full_rank_exponent"] == 92
    assert tuple(modular["pivot_exponents"]) == PIVOT_EXPONENTS
    assert modular["held_out_coefficients_matched_per_channel"] == 28
    assert modular["reconstructed_through"] == 800
    assert phi_scale() > 0
    assert set(ENDPOINTS_E800) == set(
        row["exact_broad_conservation_lp"]["intervals"]
    )
    for case in ENDPOINTS_E800:
        lower, upper = phi_interval(case)
        assert lower < upper
        assert lower < 6
    assert row["proved"]["broad_channel_cone_closes_R1"] is False


def test_prop15668_finite_p11_strong_r1_is_exact_and_scoped():
    audit = finite_p11_r1_audit()
    assert P11_DELTA_SQ == Fraction(1_382_747_375_360, 583_792_784_981)
    assert Fraction(audit["strong_margin"]) > 0
    assert Fraction(audit["exact_R1_margin"]) > 0
    assert audit["strong_R1_holds"] is True
    assert audit["exact_R1_holds"] is True
    assert audit["source"] == "evidence/e1_gmin_r1_principal_pge11.json"
    assert audit["source_sha256"] == ARTIFACT_SHA256["finite_p11_r1_record"]
    row = theorem_record()
    assert row["proved"]["finite_p11_strong_R1_from_independent_full_census"]
    assert row["proved"]["general_R1"] is False
    assert row["proved"]["global_QVAR"] is False
    assert row["proved"]["residual_ii"] is False
    assert row["proved"]["limit_exists"] is False
    assert row["L_status"] == "OPEN"


def test_prop15668_artifact_hashes_and_committed_record_are_pinned():
    assert ARTIFACT_SHA256
    assert all(len(value) == 64 for value in ARTIFACT_SHA256.values())
    committed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15668.json").read_text()
    )
    assert committed == theorem_record()
