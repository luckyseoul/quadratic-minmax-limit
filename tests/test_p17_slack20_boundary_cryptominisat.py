import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p17_slack20_boundary_cryptominisat.py"
SPEC = importlib.util.spec_from_file_location("p17_slack20_boundary_cms", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_exact_census_and_signature_deduplication():
    profiles = MODULE.slack20_profiles()
    manifest = MODULE.signature_manifest()

    assert len(profiles) == 193
    assert [row["census_index"] for row in profiles] == list(range(1364, 1557))
    assert all(row["pair_slack"] == 20 for row in profiles)
    assert all(
        2 <= row["phase_profiles_b"]["0"].get(0, 0) <= 6
        for row in profiles
    )
    assert manifest["profile_count"] == 193
    assert manifest["signature_count"] == 184
    assert manifest["multiplicity_histogram"] == {"1": 175, "2": 9}
    assert sorted(
        index
        for signature in manifest["signatures"]
        for index in signature["profile_indices"]
    ) == list(range(193))


def test_normalization_and_residual_reflection_accounting():
    manifest = MODULE.signature_manifest()
    normalization = manifest["normalization"]
    reflection = manifest["residual_direction_reflection"]

    assert normalization == {
        "mode": "phase-zero-b0-pair",
        "selected_points": [0, 1],
        "canonical_direction": [0, 1],
        "phase": 0,
        "b": 0,
        "c_H": -1,
        "lossless_for_every_profile": True,
    }
    assert reflection["phase_zero_fixed"] == [17]
    assert reflection["phase_one_fixed"] == [0]
    assert reflection["raw_assignment_count"] == 1_971_382
    assert reflection["orbit_count"] == 985_730
    assert reflection["fixed_assignment_count"] == 78
    assert min(
        row["direction_assignment_count_after_normalization"]
        for row in manifest["signatures"]
    ) == 280
    assert max(
        row["direction_assignment_count_after_normalization"]
        for row in manifest["signatures"]
    ) == 60_480


def test_radon_geometry_has_expected_native_xor_dimensions():
    geometry = MODULE.radon_geometry()
    assert geometry["normalized_index"] == 17
    assert geometry["directions"][17] == (0, 1)
    assert geometry["point_variables"] == 289
    assert geometry["line_parity_variables"] == 306
    assert geometry["native_xor_constraints"] == 595
    assert sum(row[2] == 0 for row in geometry["records"]) == 9
    assert sum(row[2] == 1 for row in geometry["records"]) == 9


def test_profile_and_signature_cli_resolve_to_same_case():
    manifest = MODULE.signature_manifest()
    duplicate = next(
        row for row in manifest["signatures"] if row["multiplicity"] == 2
    )
    by_signature = MODULE.resolve_signature(
        signature_index=duplicate["signature_index"]
    )
    first = MODULE.resolve_signature(profile_index=duplicate["profile_indices"][0])
    second = MODULE.resolve_signature(profile_index=duplicate["profile_indices"][1])
    assert by_signature is first or by_signature == first
    assert first == second

    with pytest.raises(ValueError):
        MODULE.resolve_signature()
    with pytest.raises(ValueError):
        MODULE.resolve_signature(signature_index=0, profile_index=0)
    with pytest.raises(ValueError):
        MODULE.resolve_signature(signature_index=184)
    with pytest.raises(ValueError):
        MODULE.resolve_signature(profile_index=193)


def test_atomic_json_output_and_manifest_cli(tmp_path):
    output = tmp_path / "nested" / "manifest.json"
    MODULE.main(["--list-signatures", "--output", str(output)])
    payload = json.loads(output.read_text())
    assert payload["profile_count"] == 193
    assert payload["signature_count"] == 184
    assert not list(output.parent.glob(".manifest.json.*.tmp"))


def test_auditor_rejects_a_normalized_set_with_wrong_profile():
    signature = MODULE.resolve_signature(signature_index=0)
    audit = MODULE.audit_boundary(range(16), signature["phase_profiles_b"])
    assert audit["point_set_valid"] is True
    assert audit["normalization_valid"] is True
    assert audit["inverse_radon_valid"] is True
    assert audit["phase_histograms_valid"] is False
    assert audit["valid"] is False
