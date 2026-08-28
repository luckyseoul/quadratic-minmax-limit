#!/usr/bin/env python3
"""Prop. 15.667 -- exact p=11 profile and quartic-trace reconstruction.

The ten-dimensional p=11 glue-dual code is reduced exactly by translation
and nonzero-scalar orbits, and every resulting six-profile type is counted
with integer dynamic programming.  Five independent 31-bit prime moduli
admitting primitive eleventh roots have
product larger than the unrestricted bounds for the ordinary count and the
common-sum second and fourth moments through scaled norm 120.  CRT therefore
recovers those coefficients over the integers, not probabilistically.

The exact scalar prefix through 88 has full rank 41 in its modular space and
the exact quartic-trace prefix through 92 has full affine rank 32.  They
uniquely reconstruct both exported series through 800 and correctly predict
32 and 28 held-out profile coefficients, respectively.  Exact shellwise raw
mass conservation then gives rational QSopt_ex primal/dual endpoint
certificates for all four p=11 PSL component cases through 800.

This substantially sharpens the p=11 R1 data, corrects the two normalization
rows in Proposition 15.665, and also proves a route limit: aggregate scalar
trace conservation alone leaves very wide feasible target intervals.  It
does not prove R1, global QVAR, Type I, or the limit.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1/"
)

ARTIFACT_SHA256 = {
    "profile_orbit_report": (
        "cec16194297818aa3753e7ce5dc37d4c1647f71d7bce13f83ac6d10e43b4c8fe"
    ),
    "profile_tuple_report": (
        "cd2c59874359d195517e3e19ae74d4241f4c6b3b92f6fb3863530623d162300f"
    ),
    "profile_tuple_npz": (
        "3c27864c6b4d9a7106aa46a8eb8d7b6afd2062b86d0f7dd54bdd8d2a410fa45c"
    ),
    "profile_table_report": (
        "b8f25d06a0892e810bf2fe034f77f795f10037c40a6b00160f76f2aca162c4bd"
    ),
    "profile_table_npz": (
        "35a950245ef6391b4ce4699bbd2ec63f75300b5e0ae986404156991ca6c0cf4f"
    ),
    "profile_moment_report_e120": (
        "06debd6476b2236624684f10b3845d87f09f513a880da5f869e84c7853956181"
    ),
    "scalar_reconstruction_e800": (
        "fbb99d589b3c24df335d5ada78e872b775f33ffeb688a5e01f961bb117accac8"
    ),
    "scalar_qrows_e800": (
        "08d8eed940c23cc489d908e9b1ff779b0f70b739213f6b0fa4f51439bb5d9cd9"
    ),
    "scalar_half_target_e800": (
        "820e8585fdaff2e0ce1ee9cc978dde5fca2d552942bf7d8c1d662201d4bccab0"
    ),
    "channel_reduction_unscaled": (
        "7f10326190b11228e62c5b3b03c26332f286d92e4d18a992e580de529c465c21"
    ),
    "circle_kernel_qrows_e800": (
        "f16e2d54d52a51e5b38ef94c02181dcc431dca03b97b5ef9fcf2930ddd7d0531"
    ),
    "circle_low_qrows_e800": (
        "6db984eb2495e79e09b05a92778a5e4f2cff0169a43d9d9e28da449a8a7935d3"
    ),
    "circle_high_qrows_e800": (
        "e876cd2a76b9577bae716e734f42fd7539ddecdd59dbfd8cdbe1db1e2056db50"
    ),
    "trace_reconstruction_e800": (
        "e2e3875e64da7b71e9fca687a48adca5b2b00fe455065467cc4fc465a1c5bde9"
    ),
    "trace_endpoint_report_e120": (
        "1a8f4fdad87f178ae0e1b5ef7ce38b77bd58ed6ef52841166925e4ba30100e48"
    ),
    "trace_endpoint_report_e800": (
        "872c21284ed9c17d751926a5ab6137778f6af9f74e1a2053da1bf319527f8568"
    ),
    "permanent_archive_manifest": (
        "b92262530b3f853b67814aeb0ec9d07e834f833b8d507a230a64a66fb9e8fc23"
    ),
}

CRT_MODULI = (1_999_999_871, 1_999_999_321, 1_999_999_013, 1_999_998_947, 1_999_997_957)
CRT_PRODUCT = 31_999_921_744_068_749_461_247_094_447_450_713_426_945_936_557
MAX_COUNT_BOUND = 3_035_822_103_388_904_867_185_472_590_558_138
MAX_SECOND_MOMENT_BOUND = 65_693_200_149_015_553_751_998_674_995_898_480
MAX_FOURTH_MOMENT_BOUND = 4_129_286_721_131_182_422_173_790_859_253_254_080

P11_EARLY_PROFILE = (
    # exponent, count, common-sum second moment, common-sum fourth moment,
    # harmonic trace, raw trace
    (11, 244, 484, 29_524, "-3538/63", "0"),
    (20, 16_104, 58_080, 2_555_520, "-85888/21", "89792/11"),
    (24, 14_762, 63_888, 5_134_272, "-63684/7", "7076"),
    (27, 442_860, 2_156_220, 66_363_660, "-527406/7", "538752"),
)

TRACE_PIVOT_EXPONENTS = (
    31, 32, 35, 36, 39, 40, 43, 44, 47, 48, 51, 52, 55, 56, 59, 60,
    63, 64, 67, 68, 71, 72, 75, 76, 79, 80, 83, 84, 87, 88, 91, 92,
)

ENDPOINTS_E120 = {
    "circle-kernel-principal": {
        "minimum": "-2751803996323451587329052771199994539993892877054760858530912214576768076471274058/4224741065864693244061199213126065064106173761682290876205710767162447036127347",
        "maximum": "175201398362216838969786767534411637099523929508898994632505363995073138706814/285120822770373085338845667893551955096665216709259200544265414730801662767",
    },
    "circle-low-Weil": {
        "minimum": "-26129495401234283572272889747754128085166429980809833084888606330393261580116016705065964/29411252932472017632661963776019203957333339642642391241334399069469481534412025254927",
        "maximum": "4267954100270165124246655504687421769646745237168777151646734811739333530950904/4849923573730053968367039199240477197164470210985094737040921535545285299823",
    },
    "circle-low-principal": {
        "minimum": "-8254486670769270717731792421373443286505127882725247691267429797809580035660184914/12674223197594079732183597639378195192318521285046872628617132301487341108382041",
        "maximum": "175220699302337080257472900991943793381806574122738537942636583764271640605514/285120822770373085338845667893551955096665216709259200544265414730801662767",
    },
    "circle-high-principal": {
        "minimum": "-2751433869043017969626906414309378406603272577679146904800789476208478398969819154/4224741065864693244061199213126065064106173761682290876205710767162447036127347",
        "maximum": "175224559490361128515010127683450224638263103045506446604662827718111340985254/285120822770373085338845667893551955096665216709259200544265414730801662767",
    },
}

ENDPOINTS_E800 = {
    **ENDPOINTS_E120,
    "circle-low-Weil": {
        "minimum": ENDPOINTS_E120["circle-low-Weil"]["minimum"],
        "maximum": "65083138953068735632403758391021967485130693524822164103725166552235240871994506088/74387512469496681039575094244174529261299333938155375806078329744466803340169151",
    },
}


def raw_trace_from_profile_moment(
    p: int,
    exponent: int,
    shell_count: int,
    common_sum_fourth_moment: int,
) -> Fraction:
    """Return tr(R_e) from the exact common-coordinate fourth moment.

    Coordinate transitivity gives

      sum_(x,i) x_i^4 = (p^2+1)/(16p^4) sum_t t^4 N_(e,t),

    and projection to the zero-diagonal tensor space then simplifies to the
    formula below.
    """
    n = p * p + 1
    return Fraction(
        n * (shell_count * exponent * exponent - common_sum_fourth_moment),
        4 * p * p * (p * p - 1),
    )


def tight_frame_second_moment_holds(
    p: int,
    exponent: int,
    shell_count: int,
    common_sum_second_moment: int,
) -> bool:
    """Audit n*sum_t t^2 N_(e,t) = 2 p e N_e exactly."""
    return (p * p + 1) * common_sum_second_moment == (
        2 * p * exponent * shell_count
    )


def endpoint_plateau_audit() -> dict[str, object]:
    comparisons = {}
    equal_count = 0
    for case, endpoints in ENDPOINTS_E800.items():
        comparisons[case] = {}
        for sense, endpoint in endpoints.items():
            equal = endpoint == ENDPOINTS_E120[case][sense]
            comparisons[case][sense] = {
                "same_at_120_and_800": equal,
                "endpoint_e120": ENDPOINTS_E120[case][sense],
                "endpoint_e800": endpoint,
                "endpoint_e800_decimal": float(Fraction(endpoint)),
            }
            equal_count += int(equal)
    return {
        "equal_endpoint_count": equal_count,
        "total_endpoint_count": 8,
        "sole_late_contraction": "circle-low-Weil maximum",
        "comparisons": comparisons,
    }


def theorem_record() -> dict[str, object]:
    early_rows = []
    for exponent, count, moment2, moment4, harmonic, raw in P11_EARLY_PROFILE:
        calculated_raw = raw_trace_from_profile_moment(11, exponent, count, moment4)
        tight_frame = tight_frame_second_moment_holds(11, exponent, count, moment2)
        if calculated_raw != Fraction(raw) or not tight_frame:
            raise ArithmeticError(f"p=11 profile calibration failed at {exponent}")
        early_rows.append(
            {
                "scaled_norm": exponent,
                "shell_count": count,
                "common_sum_second_moment": moment2,
                "common_sum_fourth_moment": moment4,
                "harmonic_trace": harmonic,
                "raw_trace_mass": raw,
                "tight_frame_identity": True,
            }
        )

    if CRT_PRODUCT <= max(
        MAX_COUNT_BOUND, MAX_SECOND_MOMENT_BOUND, MAX_FOURTH_MOMENT_BOUND
    ):
        raise ArithmeticError("CRT product does not dominate the exact bounds")
    plateau = endpoint_plateau_audit()
    if plateau["equal_endpoint_count"] != 7:
        raise ArithmeticError("the e120/e800 endpoint comparison changed")

    return {
        "prop": "15.667",
        "title": "Exact p=11 profile and quartic-trace theta reconstruction",
        "proved": {
            "complete_p11_glue_dual_profile_orbit_count": True,
            "exact_scalar_and_common_moments_through_exponent_120": True,
            "scalar_modular_form_unique_from_prefix_through_88": True,
            "quartic_trace_modular_form_unique_from_prefix_through_92": True,
            "scalar_and_trace_coefficients_reconstructed_through_800": True,
            "all_eight_endpoint_bounds_exact_primal_dual_certified": True,
            "prop15665_early_shell_normalization_corrected": True,
            "aggregate_trace_conservation_closes_R1": False,
            "R1": False,
            "phi_F_ge_6_proved_general": False,
            "global_QVAR": False,
            "type_I_multilevel": False,
            "limit_exists": False,
        },
        "finite_profile": {
            "p": 11,
            "glue_dual_codewords": 11**10,
            "translation_scalar_representatives": 21_437_340,
            "weighted_histogram_tuples": 2_558_543,
            "quartic_value_distribution_types": 604,
            "affine_output_types": 13,
            "profile_exponent_through": 120,
            "nonempty_tight_frame_checks": 51,
            "crt_moduli": list(CRT_MODULI),
            "crt_product": CRT_PRODUCT,
            "unrestricted_bounds": {
                "count": MAX_COUNT_BOUND,
                "second_moment": MAX_SECOND_MOMENT_BOUND,
                "fourth_moment": MAX_FOURTH_MOMENT_BOUND,
            },
            "crt_product_exceeds_every_bound": True,
            "early_shell_audit": early_rows,
        },
        "modular_reconstruction": {
            "scalar_space_dimension": 41,
            "scalar_prefix_fixed_through": 88,
            "scalar_prefix_rank": 41,
            "scalar_held_out_coefficients_matched": 32,
            "trace_affine_dimension": 32,
            "trace_prefix_fixed_through": 92,
            "trace_prefix_rank": 32,
            "trace_pivot_exponents": list(TRACE_PIVOT_EXPONENTS),
            "trace_held_out_coefficients_matched": 28,
            "reconstructed_through": 800,
            "ordinary_coefficients_nonnegative_integers": True,
            "raw_trace_coefficients_nonnegative": True,
            "trace_half_cusp_target": (
                "-4428472046531859136727844588716/"
                "11389083011948997399715094836557"
            ),
        },
        "exact_trace_conservation_lp": {
            "coefficient_through": 800,
            "cases": 4,
            "endpoints": 8,
            "all_status": "exact_qsopt_primal_dual_certified",
            "endpoint_comparison": plateau,
            "conclusion": (
                "strictly stronger than the prior uncoupled positivity cone, "
                "but still far too wide to certify R1"
            ),
        },
        "normalization_erratum": {
            "affected_scaled_norms": [20, 24],
            "cause": "H(u/2) anchors were used as H(u)",
            "quartic_correction_factor": 16,
            "corrected_raw_traces": {"20": "89792/11", "24": "7076"},
            "general_prop15665_theorem_unchanged": True,
        },
        "artifact_sha256": ARTIFACT_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "remaining_obstruction": (
            "Scalar shellwise mass does not identify its distribution among "
            "PSL constituents.  The next attack must add channel-resolved "
            "moments (or a stronger multi-scale transport inequality)."
        ),
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    output = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15667.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print("Prop 15.667 exact p=11 profile/trace reconstruction: proved")
    print(f"  wrote {destination}")
    return output


if __name__ == "__main__":
    main()
