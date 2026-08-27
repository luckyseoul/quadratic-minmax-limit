#!/usr/bin/env python3
"""Prop. 15.665 -- conserved positive quartic mass on every R1 dual shell.

For a complete equal-norm shell X in the Paley dual lattice and

    Z = {W = W^T : PWP = W, diag(W) = 0},

put b_x = projection_Z(xx^T) and R_X = sum_x b_x tensor b_x.  Then R_X is
positive semidefinite.  The degree-four harmonic shell operator is

    A_X = R_X - rho_X I,   rho_X = 2 |X| r^2/[d(d+2)],

where r=||x||^2.  Every PSL constituent of Z is multiplicity-free, so if
q_{X,c} is the scalar of R_X on a constituent of dimension m_c, then

    q_{X,c} >= 0,   sum_c m_c q_{X,c} = tau_X,
    0 <= q_{X,c} <= tau_X/m_c.

The trace tau_X is the coefficient of one explicit scalar weighted theta
series.  This is new nonlinear channel coupling, but it does not by itself
prove R1; the transformed target still needs a certified theta inequality.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15631 import harmonic_min_shell_sum
from e1_gmin_m4_prop15634 import second_shadow_spectrum
from e1_gmin_m4_prop15635 import third_pair_harmonic_coefficient
from e1_gmin_m4_prop15640 import harmonic_spectrum

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ("circle-kernel", "circle-low", "circle-high")


def dimensions(p: int) -> tuple[int, int, int]:
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    return n, d, zdim


def constituent_partition(p: int) -> list[dict[str, int | str]]:
    """Multiplicity-free real PSL constituents grouped by circle channel."""
    n, d, _zdim = dimensions(p)
    rows: list[dict[str, int | str]] = [
        {
            "channel": "circle-kernel",
            "kind": "principal",
            "dimension": n,
            "count": (p - 1) * (p - 3) // 8,
        }
    ]
    if p % 4 == 1:
        rows.extend(
            [
                {
                    "channel": "circle-low",
                    "kind": "principal",
                    "dimension": n,
                    "count": (p - 1) // 4,
                },
                {
                    "channel": "circle-high",
                    "kind": "Weil",
                    "dimension": d,
                    "count": 1,
                },
                {
                    "channel": "circle-high",
                    "kind": "principal",
                    "dimension": n,
                    "count": (p - 5) // 4,
                },
            ]
        )
    else:
        rows.extend(
            [
                {
                    "channel": "circle-low",
                    "kind": "Weil",
                    "dimension": d,
                    "count": 1,
                },
                {
                    "channel": "circle-low",
                    "kind": "principal",
                    "dimension": n,
                    "count": (p - 3) // 4,
                },
                {
                    "channel": "circle-high",
                    "kind": "principal",
                    "dimension": n,
                    "count": (p - 3) // 4,
                },
            ]
        )
    rows = [row for row in rows if int(row["count"]) > 0]
    if sum(int(row["dimension"]) * int(row["count"]) for row in rows) != dimensions(p)[2]:
        raise ArithmeticError("constituent dimensions do not sum to dim Z")
    return rows


def broad_channel_dimensions(p: int) -> dict[str, int]:
    output = {channel: 0 for channel in CHANNELS}
    for row in constituent_partition(p):
        output[str(row["channel"])] += int(row["dimension"]) * int(row["count"])
    return output


def diagonal_gram_inverse(p: int) -> dict[str, Fraction]:
    """K=P hadamard P and K^-1=a I+b J on the ambient diagonal space."""
    return {
        "I_coefficient": Fraction(4 * p * p, p * p - 1),
        "J_coefficient": -Fraction(2, p * p - 1),
    }


def projected_rank_one_norm_coefficients(p: int) -> dict[str, Fraction]:
    """||projection_Z(xx^T)||^2 = a r^2 + b sum_i x_i^4."""
    return {
        "radius_fourth_coefficient": Fraction(p * p + 1, p * p - 1),
        "coordinate_fourth_coefficient": -Fraction(4 * p * p, p * p - 1),
    }


def radial_shift(p: int, shell_count: int, radius_squared: Fraction) -> Fraction:
    _n, d, _zdim = dimensions(p)
    return Fraction(2 * shell_count, d * (d + 2)) * radius_squared**2


def trace_harmonic_radial_correction(p: int) -> Fraction:
    """Coefficient of r^2 added to ||projection_Z(xx^T)||^2 in H_trace."""
    _n, d, zdim = dimensions(p)
    return -Fraction(4 * zdim, d * (d + 4)) + Fraction(
        2 * zdim, (d + 2) * (d + 4)
    )


def trace_harmonic_zonal_identity(p: int) -> dict[str, Fraction]:
    """Express H_trace as the orbit sum of one coordinate-zonal quartic.

    For ell_i(x)=x_i on range(P), ||ell_i||^2=P_ii=1/2.  Thus

      Z_i = ell_i^4 - 3*r*ell_i^2/(d+4)
            + 3*r^2/[4(d+2)(d+4)]

    is harmonic.  Coordinate transitivity and sum_i ell_i^2=r give
    H_trace = orbit_factor * Z_i at the level of theta series.
    """
    n, d, _zdim = dimensions(p)
    coordinate_coefficient = projected_rank_one_norm_coefficients(p)[
        "coordinate_fourth_coefficient"
    ]
    radial_coefficient = (
        projected_rank_one_norm_coefficients(p)["radius_fourth_coefficient"]
        + trace_harmonic_radial_correction(p)
    )
    expected_radial = -coordinate_coefficient * Fraction(3, 2 * (d + 2))
    if radial_coefficient != expected_radial:
        raise ArithmeticError("trace polynomial is not the claimed harmonic orbit sum")
    return {
        "coordinate_fourth_coefficient": coordinate_coefficient,
        "radius_fourth_coefficient": radial_coefficient,
        "one_coordinate_theta_factor": coordinate_coefficient * n,
    }


def p11_early_shell_audit() -> list[dict[str, object]]:
    p = 11
    n, d, zdim = dimensions(p)
    channel_dims = broad_channel_dimensions(p)
    counts = {
        p: 2 * n,
        2 * (p - 1): p * (p + 1) * n,
        2 * (p + 1): p * p * n,
        3 * p - 6: p * p * (p - 1) * (p + 7) * n // 6,
    }
    spectra = {
        p: {channel: harmonic_min_shell_sum(p) for channel in CHANNELS},
        2 * (p - 1): {
            str(row["channel"]): Fraction(row["eigenvalue"])
            for row in second_shadow_spectrum(p)
        },
        2 * (p + 1): {
            channel: third_pair_harmonic_coefficient(p) for channel in CHANNELS
        },
        3 * p - 6: {
            str(row["channel"]): Fraction(row["eigenvalue"])
            for row in harmonic_spectrum(p)
        },
    }
    expected_tau = {
        11: Fraction(0),
        20: Fraction(923784, 77),
        24: Fraction(436943, 28),
        27: Fraction(538752),
    }
    rows = []
    for exponent in sorted(counts):
        trace_harmonic = sum(
            (
                Fraction(channel_dims[channel]) * spectra[exponent][channel]
                for channel in CHANNELS
            ),
            Fraction(),
        )
        rho = radial_shift(p, counts[exponent], Fraction(exponent, 2 * p))
        tau = trace_harmonic + zdim * rho
        if tau != expected_tau[exponent] or tau < 0:
            raise ArithmeticError(f"p=11 shell trace mismatch at {exponent}")
        raw_by_channel = {
            channel: spectra[exponent][channel] + rho for channel in CHANNELS
        }
        if any(value < 0 for value in raw_by_channel.values()):
            raise ArithmeticError(f"negative raw channel at exponent {exponent}")
        if sum(
            (Fraction(channel_dims[channel]) * raw_by_channel[channel] for channel in CHANNELS),
            Fraction(),
        ) != tau:
            raise ArithmeticError(f"channel conservation mismatch at {exponent}")
        rows.append(
            {
                "scaled_norm": exponent,
                "shell_count": counts[exponent],
                "harmonic_trace": str(trace_harmonic),
                "radial_shift": str(rho),
                "raw_trace_mass": str(tau),
                "raw_channel_eigenvalues": {
                    channel: str(value) for channel, value in raw_by_channel.items()
                },
            }
        )
    return rows


def theorem_record() -> dict[str, object]:
    partition_checks = {
        str(p): {
            "Z_dimension": dimensions(p)[2],
            "constituents": constituent_partition(p),
            "broad_channel_dimensions": broad_channel_dimensions(p),
        }
        for p in (5, 7, 11, 13, 17, 19)
    }
    return {
        "prop": "15.665",
        "title": "Conserved positive quartic mass on every R1 dual shell",
        "proved": {
            "raw_shell_operator_positive_semidefinite": True,
            "harmonic_shell_is_raw_minus_radial_scalar": True,
            "multiplicity_weighted_raw_mass_is_conserved": True,
            "channel_upper_bound_tau_over_dimension": True,
            "closed_trace_harmonic_polynomial": True,
            "R1": False,
            "phi_F_ge_6_proved_general": False,
            "global_QVAR": False,
        },
        "formulas": {
            "Z_dimension": "n(n-6)/8",
            "raw_radial_shift": "2*N*r^2/[d(d+2)] where r=||x||^2",
            "projected_rank_one_norm": (
                "r^2 - 4p^2/(p^2-1)*sum_i x_i^4 + 2r^2/(p^2-1)"
            ),
            "channel_bound": "0 <= q_(s,c) <= tau_s/dim(c)",
            "trace_zonal_theta_factor": "-4p^2(p^2+1)/(p^2-1)",
        },
        "partition_checks": partition_checks,
        "p11_early_shell_audit": p11_early_shell_audit(),
        "remaining_obstruction": (
            "A certified modular or multi-scale theta inequality must still "
            "transport these shellwise positive masses to the odd-coset target."
        ),
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    output = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15665.json"
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print("Prop 15.665 conserved positive shell mass: proved")
    print(f"  wrote {destination}")
    return output


if __name__ == "__main__":
    main()
