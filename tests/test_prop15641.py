from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import e1_gmin_m4_prop15641 as prop  # noqa: E402


def test_exact_rank_and_scope() -> None:
    result = prop.theorem()
    assert result["kohnen_dimension"] == 66
    assert result["rank_without_second_shell"] == 29
    assert result["rank_with_second_shell"] == 30
    assert result["residual_dimension"] == 36
    assert result["joint_second_target_rank"] == 2
    assert result["target_after_second_rank"] == 1
    assert result["target_second_proportional"] is False
    assert result["exact_witness"]["known_constraints_zero"] is True
    assert result["exact_witness"]["second_shell_zero"] is True
    assert result["exact_witness"]["target"] == 1
    assert result["does_not_refute_R1"] is True
    assert result["closes_R1"] is False


def test_reproduction_scripts_present() -> None:
    for name in (
        "r1_p11_kohnen_cache.gp",
        "r1_p11_kohnen_reduce.gp",
        "r1_p11_second_shell_rank.gp",
        "r1_p11_modular_independence_witness.gp",
        "r1_p11_verify_modular_witness.gp",
    ):
        assert (ROOT / "scripts" / name).is_file()


def test_backed_up_artifact_when_mounted() -> None:
    if prop.ARTIFACT.is_file():
        assert prop.sha256(prop.ARTIFACT) == prop.ARTIFACT_SHA256
