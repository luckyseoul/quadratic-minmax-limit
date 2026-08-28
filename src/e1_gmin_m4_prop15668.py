#!/usr/bin/env python3
"""Prop. 15.668 -- exact p=11 broad-channel theta reconstruction.

The marked p=11 profile census refines Proposition 15.667 by retaining the
fourth Legendre-convolution statistic needed to split every raw shell
operator among the kernel, low, and high eigenspaces of the square-circle
operator.  Five-modulus CRT gives those three masses exactly through
exponent 120.  Their 32-dimensional affine modular spaces all reach full
rank at exponent 92, so the three series and their half-cusp targets are
uniquely reconstructed through exponent 800; 28 coefficients per channel
are held out and reproduced exactly.

Separate shell-mass and transformed-target conservation in those three
channels gives eight rational QSopt_ex endpoint certificates.  This is a
strict refinement of aggregate trace conservation, but the certified cone
still contains principal target values whose corresponding Phi eigenvalue
is below six.  Thus the broad-channel route does not prove R1.

Independently, the previously computed full p=11 census gives
``||delta||^2 < n/12`` exactly.  That is a genuine finite theorem and is
recorded here to keep it distinct from the failed cone relaxation; it is
not an all-prime R1 proof.  No general gate is flipped.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-28-r1-broad-channel/"
)

ARTIFACT_SHA256 = {
    "finite_p11_r1_record": (
        "bc78840eadf843db041cc601f949b0f305a1a8027b3920bc95d37029e616a6f4"
    ),
    "profile_types_report": (
        "2896ea32b9c93e0cfeb623745b66fdfdf73ce0675569c8234da818e8967fe918"
    ),
    "profile_types_npz": (
        "3e2375e8f5b938b3c2c1f5054966c46e187cf3f1d948c58b946d37b4e9380b7f"
    ),
    "channel_tuples_report": (
        "abeb850686e8136c87748f7ab618d2aa8888532c97fb503b6fff010245d226a0"
    ),
    "channel_tuples_npz": (
        "e873061507543195c266ef27a5dbd0f853d38a4c9cbcf0773878054be6a845d6"
    ),
    "profile_tables_report": (
        "e324cfce6269bd9882140bed46ebd482859e0c97a0b4713df4d7239ae6521390"
    ),
    "profile_tables_npz": (
        "5fbbf88cf811d32e1261c62006e63f39bdf399b79f418b5bda4f2146162b1165"
    ),
    "channel_moments_e120": (
        "c8a43925bb3b5800fa7f875a454b477bcd2499354e2fefe925b4bfe1595404fd"
    ),
    "broad_reconstruction_e800": (
        "17b732ec3e23161aec5b47b8cbec516157bc5285968df9d38df70507be1d98aa"
    ),
    "broad_endpoint_report_e800": (
        "981ec5d76c98ed26491391422b1c9da919f276d0c4b7f9befa182136bbabd743"
    ),
    "permanent_archive_manifest": (
        "d1ef69b9af7007c0d2f09a3a5ea8a014cde62d9ed6109175cf4a6496d06b3f07"
    ),
}

CRT_PRODUCT = 31_999_921_744_068_749_461_247_094_447_450_713_426_945_936_557
MAX_COUNT_BOUND = 3_035_822_103_388_904_867_185_472_590_558_138
MAX_SQUARED_EXCESS_BOUND = 16_393_439_358_300_086_282_801_551_989_013_945_200
MAX_PROFILE_LEGENDRE_FOURTH_BOUND = (
    106_481_782_073_508_797_440_446_372_773_921_289_512_628
)

BROAD_DIMENSIONS = {
    "circle-kernel": 1220,
    "circle-low": 305,
    "circle-high": 244,
}
BROAD_TARGETS = {
    "circle-kernel": Fraction(
        -10_463_154_194_187_058_501_821_423_212,
        56_945_415_059_744_986_998_575_474_182_785,
    ),
    "circle-low": Fraction(
        -19_210_249_628_300_203_741_452_825_212,
        56_945_415_059_744_986_998_575_474_182_785,
    ),
    "circle-high": Fraction(
        -2_883_758_999_278_296_860_880_307_324,
        11_389_083_011_948_997_399_715_094_836_557,
    ),
}
AGGREGATE_TRACE_TARGET = Fraction(
    -4_428_472_046_531_859_136_727_844_588_716,
    11_389_083_011_948_997_399_715_094_836_557,
)

PIVOT_EXPONENTS = (
    31, 32, 35, 36, 39, 40, 43, 44, 47, 48, 51, 52, 55, 56, 59, 60,
    63, 64, 67, 68, 71, 72, 75, 76, 79, 80, 83, 84, 87, 88, 91, 92,
)

ENDPOINTS_E800 = {
    "circle-kernel-principal": {
        "minimum": Fraction(
            -20393214511994032102477715910967308502211567924873803328450969019796376791929093752,
            38997734455663404666676299368603589516164284351106413958321849410422191158270821,
        ),
        "maximum": Fraction(
            626109633353046861967860955031899963605092893773873611780782307993937793339,
            1231302858699550783534704421148631966771030827344055951972256735692927850,
        ),
    },
    "circle-low-Weil": {
        "minimum": Fraction(
            -614637209198103599449909593052575420613248482308606575340734290927431054795553554,
            1607293310935398445413194957995739623334893743024074035602997588152140810613225,
        ),
        "maximum": Fraction(
            3497856335795515178722001254323006595838887219230673902358130345020609405768,
            8901422367137672027534852695997733404019350295080661376130321391096381825,
        ),
    },
    "circle-low-principal": {
        "minimum": Fraction(
            -38171263947463143453483304508320244077437222634283532132953798378885337472912,
            148454955807949385415452383886835932549038394718977206437499988565535609887,
        ),
        "maximum": Fraction(
            48879250999250215194622957387645997857946746508476527911086247554517678763574,
            181956604590918344685436362505964389642665247765210976024426908142370657345,
        ),
    },
    "circle-high-principal": {
        "minimum": Fraction(
            -149998907810299899467846909035777716383932183836983531367583527535586143249,
            684211297483694608038023106534426480852365916696807470560172413472759193,
        ),
        "maximum": Fraction(
            149998561320510606268132354344421762525301832767317650482392871095164379897,
            684211297483694608038023106534426480852365916696807470560172413472759193,
        ),
    },
}

CERTIFICATE_SIZES = {
    "circle-kernel-principal": (160, 3091),
    "circle-low-Weil": (128, 2705),
    "circle-low-principal": (160, 3091),
    "circle-high-principal": (160, 3091),
}

P11_DELTA_SQ = Fraction(1_382_747_375_360, 583_792_784_981)


def phi_scale() -> Fraction:
    """Poisson conversion from a harmonic target to its Phi scalar."""
    p = 11
    n = p * p + 1
    zdim = n * (n - 6) // 8
    lbar = Fraction(8 * (n - 2), n - 6)
    spherical = Fraction(8 * n, n + 4)
    return -Fraction(zdim) * (lbar - spherical) / AGGREGATE_TRACE_TARGET


def phi_interval(case: str) -> tuple[Fraction, Fraction]:
    """Map one certified harmonic-target interval to a Phi interval."""
    endpoints = ENDPOINTS_E800[case]
    spherical = Fraction(8 * 122, 126)
    scale = phi_scale()
    return (
        spherical - scale * endpoints["maximum"],
        spherical - scale * endpoints["minimum"],
    )


def finite_p11_r1_audit() -> dict[str, object]:
    """Exact finite-p=11 R1 margins from the independent full census."""
    n = 122
    lbar = Fraction(240, 29)
    strong_threshold = Fraction(n, 12)
    exact_threshold = Fraction(n, 48) * (lbar - 6) ** 2
    return {
        "delta_squared": str(P11_DELTA_SQ),
        "strong_n_over_12_threshold": str(strong_threshold),
        "strong_margin": str(strong_threshold - P11_DELTA_SQ),
        "exact_R1_threshold": str(exact_threshold),
        "exact_R1_margin": str(exact_threshold - P11_DELTA_SQ),
        "strong_R1_holds": P11_DELTA_SQ <= strong_threshold,
        "exact_R1_holds": P11_DELTA_SQ <= exact_threshold,
        "source": "evidence/e1_gmin_r1_principal_pge11.json",
        "source_sha256": ARTIFACT_SHA256["finite_p11_r1_record"],
        "scope": "independent exact full Max+ census; not the broad-channel cone",
    }


def theorem_record() -> dict[str, object]:
    if CRT_PRODUCT <= max(
        MAX_COUNT_BOUND,
        MAX_SQUARED_EXCESS_BOUND,
        MAX_PROFILE_LEGENDRE_FOURTH_BOUND,
    ):
        raise ArithmeticError("CRT product does not dominate all marked bounds")
    if sum(BROAD_DIMENSIONS.values()) != 1769:
        raise ArithmeticError("broad dimensions no longer sum to dim Z")
    weighted_target = sum(
        (
            BROAD_DIMENSIONS[channel] * BROAD_TARGETS[channel]
            for channel in BROAD_DIMENSIONS
        ),
        Fraction(),
    )
    if weighted_target != AGGREGATE_TRACE_TARGET:
        raise ArithmeticError("broad transformed targets do not conserve trace")
    if phi_scale() <= 0:
        raise ArithmeticError("unexpected Poisson target orientation")

    intervals = {
        case: {
            "harmonic_minimum": str(endpoints["minimum"]),
            "harmonic_maximum": str(endpoints["maximum"]),
            "harmonic_minimum_decimal": float(endpoints["minimum"]),
            "harmonic_maximum_decimal": float(endpoints["maximum"]),
            "Phi_minimum": str(phi_interval(case)[0]),
            "Phi_maximum": str(phi_interval(case)[1]),
            "Phi_minimum_decimal": float(phi_interval(case)[0]),
            "Phi_maximum_decimal": float(phi_interval(case)[1]),
            "contains_sub_six_Phi_value": phi_interval(case)[0] < 6,
            "variables": CERTIFICATE_SIZES[case][0],
            "constraints_verified_per_endpoint": CERTIFICATE_SIZES[case][1],
            "dual_stationarity_equations_verified_per_endpoint": (
                CERTIFICATE_SIZES[case][0]
            ),
        }
        for case, endpoints in ENDPOINTS_E800.items()
    }
    finite = finite_p11_r1_audit()
    if finite["strong_R1_holds"] is not True:
        raise ArithmeticError("independent finite p=11 R1 audit changed")

    return {
        "prop": "15.668",
        "title": "Exact p=11 broad-channel theta reconstruction and route limit",
        "proved": {
            "exact_marked_profile_moments_through_exponent_120": True,
            "three_broad_theta_series_unique_from_prefix_through_92": True,
            "three_broad_theta_series_reconstructed_through_800": True,
            "twenty_eight_held_out_coefficients_match_per_channel": True,
            "all_broad_raw_masses_nonnegative_through_800": True,
            "broad_masses_and_targets_conserve_aggregate_trace": True,
            "all_eight_broad_endpoints_exact_primal_dual_certified": True,
            "broad_channel_cone_closes_R1": False,
            "finite_p11_strong_R1_from_independent_full_census": True,
            "general_R1": False,
            "phi_F_ge_6_proved_general": False,
            "global_QVAR": False,
            "type_I_multilevel": False,
            "residual_ii": False,
            "limit_exists": False,
        },
        "marked_profile_census": {
            "quartic_coefficient_tuples": 11**4,
            "domain_affine_profile_types": 1007,
            "domain_output_affine_dynamic_programs": 20,
            "translation_scalar_representatives": 21_437_340,
            "unique_sorted_channel_profile_tuples": 2_584_901,
            "dual_codewords_reconstructed": 11**10,
            "cpu_gpu_entries_audited": 268,
            "crt_moduli": 5,
            "crt_product": CRT_PRODUCT,
            "crt_product_exceeds_every_bound": True,
            "ordinary_profile_coefficients_matched": True,
            "classified_shell_channel_masses_matched": True,
        },
        "modular_reconstruction": {
            "broad_dimensions": BROAD_DIMENSIONS,
            "broad_half_cusp_targets": {
                key: str(value) for key, value in BROAD_TARGETS.items()
            },
            "aggregate_trace_target": str(AGGREGATE_TRACE_TARGET),
            "affine_dimension_per_channel": 32,
            "prefix_fixed_through": 92,
            "first_full_rank_exponent": 92,
            "pivot_exponents": list(PIVOT_EXPONENTS),
            "profile_available_through": 120,
            "held_out_coefficients_matched_per_channel": 28,
            "reconstructed_through": 800,
        },
        "exact_broad_conservation_lp": {
            "coefficient_through": 800,
            "cases": 4,
            "endpoints": 8,
            "all_status": "exact_qsopt_primal_dual_certified",
            "poisson_scale": str(phi_scale()),
            "intervals": intervals,
            "route_limit": (
                "The exact broad-channel cone still contains certified "
                "principal target points mapping below Phi=6; finer "
                "character-resolved information or a uniform inequality is required."
            ),
        },
        "finite_p11_R1": finite,
        "artifact_sha256": ARTIFACT_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "literature_and_oeis_scope": (
            "Weighted-theta literature confirms the modular framework; targeted "
            "exact-number OEIS searches found no matching sequence and supply no theorem."
        ),
        "remaining_obstruction": (
            "Prove R1 and mixed-character QVAR uniformly for all primes p>=11; "
            "the finite p=11 census cannot be promoted to a p-law."
        ),
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    output = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15668.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print("Prop 15.668 exact p=11 broad-channel reconstruction: proved")
    print(f"  wrote {destination}")
    return output


if __name__ == "__main__":
    main()
