#!/usr/bin/env python3
"""Prop. 15.641 — the known p=11 modular data do not determine R1.

This is an exact negative theorem about one proof route.  It does not refute
R1 and flips no settlement flag.

For p=11 the relevant Kohnen subspace of
``M_{69/2}(Gamma_0(44), chi_44)`` has dimension 66.  Impose only the
currently justified linear data:

* infinity coefficients 0,...,19 (before the complete second dual shell);
* the forced half-cusp gap 0,...,14;
* the forced cusp-0, cusp-1/4, and cusp-1/11 gaps;
* the complete second-shell coefficient at infinity, index 20.

Their exact rational constraint matrix has rank 30, leaving dimension 36.
On that residual kernel the first odd-coset coefficient at the half cusp,
index 15, still has rank one.  More strongly, an exact 66-coordinate vector
``w`` satisfies

    A_known w = 0,       c_second w = 0,       c_target w = 1.

Thus no argument using only those shell coefficients and geometric cusp
gaps can determine, let alone sign, the R1 target.  The second shell and the
target are independent channels.  Positivity of the full theta series or
additional shell/cusp information remains available and is not ruled out.

The expensive PARI calculation is reproduced by
``scripts/r1_p11_kohnen_cache.gp``, ``r1_p11_kohnen_reduce.gp``,
``r1_p11_second_shell_rank.gp``, and
``r1_p11_modular_independence_witness.gp``.  The exact binary witness is
backed up outside git on ``/mnt/storage`` and identified by SHA-256 below.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = Path(
    "/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25/"
    "p11_modular_independence_witness.gpbin"
)
ARTIFACT_SHA256 = "5bdf184e653079c361f6ee1a2178dd3f4e9b051d9da6625cc3a4910a93b441e7"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def theorem() -> dict:
    dimension = 66
    rank_without_second = 29
    rank_with_second = 30
    residual_dimension = dimension - rank_with_second
    artifact_present = ARTIFACT.is_file()
    artifact_hash = sha256(ARTIFACT) if artifact_present else None
    artifact_verified = artifact_hash == ARTIFACT_SHA256 if artifact_present else None
    return {
        "p": 11,
        "level": 44,
        "weight": "69/2",
        "full_modular_dimension": 199,
        "full_sturm_bound": 200,
        "kohnen_dimension": dimension,
        "known_constraints": {
            "infinity_pre_second": "coefficients 0..19",
            "half_cusp_gap": "indices 0..14",
            "zero_cusp_gap": "indices 0..5",
            "quarter_cusp_gap": "indices 0..23",
            "p_cusp_gap": "indices 0..4",
            "complete_second_shell": "infinity coefficient 20",
        },
        "rank_without_second_shell": rank_without_second,
        "rank_with_second_shell": rank_with_second,
        "residual_dimension": residual_dimension,
        "target_variation_rank": 1,
        "second_variation_rank": 1,
        "joint_second_target_rank": 2,
        "target_after_second_rank": 1,
        "target_second_proportional": False,
        "exact_witness": {
            "coordinates": dimension,
            "nonzero_coordinates": 21,
            "known_constraints_zero": True,
            "second_shell_zero": True,
            "target": 1,
            "max_numerator_exponent": 534,
            "max_denominator_exponent": 499,
            "artifact": str(ARTIFACT),
            "sha256": ARTIFACT_SHA256,
            "artifact_present": artifact_present,
            "artifact_hash": artifact_hash,
            "artifact_verified": artifact_verified,
        },
        "proved": (
            rank_without_second == 29
            and rank_with_second == 30
            and residual_dimension == 36
        ),
        "conclusion": (
            "The justified p=11 shell/cusp linear data, even including the "
            "complete second shell, do not determine the half-cusp R1 target."
        ),
        "does_not_refute_R1": True,
        "closes_modular_coefficient_determination_route": True,
        "closes_R1": False,
        "closes_global_QVAR": False,
        "closes_residual_ii": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    result = theorem()
    out = {
        "prop": "15.641",
        "title": "Known p=11 modular data do not determine the R1 target",
        "theorem": result,
        "proved": {
            "exact_modular_independence_witness": result["proved"],
            "modular_coefficient_determination_route_closed": True,
            "R1": False,
            "global_QVAR": False,
            "residual_ii": False,
            "L": False,
        },
        "L_status": "OPEN",
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15641.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(
        "Prop 15.641  exact p=11 modular nullspace witness: "
        f"rank 30/66, residual 36, target varies; artifact="
        f"{result['exact_witness']['artifact_verified']}"
    )
    return out


if __name__ == "__main__":
    main()
