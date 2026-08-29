"""Regression tests for finite-field subtraction in edge-lift models."""
from __future__ import annotations

import importlib.util
from pathlib import Path

from minmax_quadratic import paley_conference_prime_power


ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _assert_conference_data_matches_canonical(script_name: str, p: int) -> None:
    module = _load_script(script_name)
    edges, signs = module.conference_data()
    canonical = paley_conference_prime_power(p)
    assert len(edges) == len(signs) == (p * p + 1) * (p * p) // 2
    for (left, right), sign in zip(edges, signs):
        assert sign == int(canonical[left, right]), (left, right)


def test_p17_edge_lift_uses_componentwise_finite_field_subtraction():
    _assert_conference_data_matches_canonical(
        "p17_slack0_conic_edge_means_cpsat", 17
    )


def test_p19_b16_edge_lift_uses_componentwise_finite_field_subtraction():
    _assert_conference_data_matches_canonical("p19_slack20_b16_lift_cpsat", 19)


def test_p19_all_b2_edge_lift_uses_componentwise_finite_field_subtraction():
    _assert_conference_data_matches_canonical("p19_slack20_allb2_lift_cpsat", 19)
