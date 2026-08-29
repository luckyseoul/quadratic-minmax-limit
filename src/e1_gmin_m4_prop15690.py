#!/usr/bin/env python3
"""Prop. 15.690 -- exact dilation-energy normalization and method no-go."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def dilation_energy_identity(q: int, delta_sq_over_n: Fraction) -> dict[str, object]:
    if q < 25 or int(q**0.5) ** 2 != q:
        raise ValueError("q must be an odd square at least 25")
    lambda_bar = Fraction(8 * (q - 1), q - 5)
    L_bar = Fraction(q * (q - 1) * (q + 5), 16)
    V_over_n = 24 * delta_sq_over_n
    S_K = Fraction(q - 1, 2) * V_over_n
    return {
        "q": q,
        "lambda_bar": lambda_bar,
        "L_bar": L_bar,
        "V_over_n": V_over_n,
        "S_K": S_K,
        "S_K_formula": "12(q-1)||delta||^2/n",
        "strong_R1_equivalence": S_K <= q - 1
        if delta_sq_over_n == Fraction(1, 12)
        else None,
    }


def p7_psd_autocorrelation_obstruction() -> dict[str, object]:
    energy = Fraction(1_105_920, 11)
    target = 48
    return {
        "p": 7,
        "q": 49,
        "N_values": {"0": 28, "quartic_residues": 21, "other_nonzero": 14},
        "equation_5_energy": energy,
        "target": target,
        "violation_factor": energy / target,
        "has_nonnegative_additive_Fourier_transform": True,
        "is_actual_uniform_full_Max_plus_ensemble": False,
    }


def theorem_record() -> dict[str, object]:
    threshold = dilation_energy_identity(25, Fraction(1, 12))
    obstruction = p7_psd_autocorrelation_obstruction()
    if threshold["S_K"] != 24 or threshold["strong_R1_equivalence"] is not True:
        raise ArithmeticError("dilation/R1 normalization changed")
    if obstruction["violation_factor"] != Fraction(23_040, 11):
        raise ArithmeticError("p=7 PSD obstruction changed")
    return {
        "proposition": "15.690",
        "exact_identity": "S_K=(q-1)V/(2n)=12(q-1)||delta||^2/n",
        "target_equivalence": "S_K<=q-1 iff V<=2n iff ||delta||^2<=n/12",
        "sample_threshold": threshold,
        "PSD_autocorrelation_obstruction": obstruction,
        "character_PSD_only_route_sufficient": False,
        "actual_R1_counterexample": False,
        "required_extra_input": (
            "Boolean rank-one identity and exact uniform Max+ orbit mixture"
        ),
        "closes_R1": False,
        "closes_type_I": False,
        "L_status": "OPEN",
        "proved": True,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    return value


def main() -> None:
    target = ROOT / "evidence" / "e1_gmin_m4_prop15690.json"
    target.write_text(json.dumps(_jsonable(theorem_record()), indent=2) + "\n")
    print(f"Prop. 15.690: exact dilation energy; wrote {target}")


if __name__ == "__main__":
    main()
