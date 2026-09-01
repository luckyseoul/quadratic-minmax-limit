"""Exact, cache-free base certificates for Proposition 15.750."""
from __future__ import annotations

import builtins
import json
from pathlib import Path

import numpy as np
import pytest

import e1_type_i_small_prime_exact as exact
from e1_type_i_small_prime_exact import (
    EXPECTED_SHA256,
    ExactCertificateError,
    type_I_badcase_small_primes_exact,
    verify_type_i_badcase_farkas,
)


ROOT = Path(__file__).resolve().parents[1]


def test_p5_exact_farkas_identity_with_cardinality_side_condition() -> None:
    row = verify_type_i_badcase_farkas(5)
    assert row["variables"] == 325
    assert row["rows"] == 231
    assert row["distinct_eigenvectors"] == 97
    assert row["row_kinds"] == {
        "lower": 133,
        "minus_bad": 36,
        "plus_ge": 22,
        "plus_le": 39,
        "sum_le": 1,
    }
    assert row["exact_lhs_nonzeros"] == 0
    assert row["exact_rhs"] == -144
    assert row["sha256"] == EXPECTED_SHA256[5]
    assert row["uses_upper_bounds"] is False
    assert row["proved"] is True


def test_p7_strengthened_full_cone_farkas_identity() -> None:
    row = verify_type_i_badcase_farkas(7)
    assert row["variables"] == 1225
    assert row["rows"] == 1226
    assert row["distinct_eigenvectors"] == 534
    assert row["row_kinds"] == {
        "lower": 692,
        "minus_bad": 258,
        "plus_ge": 121,
        "plus_le": 155,
    }
    assert row["exact_lhs_nonzeros"] == 0
    assert row["exact_rhs"] < 0
    assert row["sha256"] == EXPECTED_SHA256[7]
    assert row["uses_upper_bounds"] is False
    assert row["proved"] is True


def test_base_verification_is_independent_of_scipy_and_npy_caches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    exact._verify_default.cache_clear()

    def no_npy_cache(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("an eigenshell cache was read")

    monkeypatch.setattr(np, "load", no_npy_cache)
    real_import = builtins.__import__

    def no_scipy(
        name: str,
        globals: dict[str, object] | None = None,
        locals: dict[str, object] | None = None,
        fromlist: tuple[str, ...] = (),
        level: int = 0,
    ) -> object:
        if name == "scipy" or name.startswith("scipy."):
            raise AssertionError("SciPy was imported")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", no_scipy)
    assert type_I_badcase_small_primes_exact() is True


def test_corrupted_multiplier_is_rejected(tmp_path: Path) -> None:
    source = ROOT / "evidence" / "e1_type_i_badcase_farkas_p5.json"
    payload = json.loads(source.read_text())
    payload["rows"][0]["weight"] += 1
    corrupted = tmp_path / "bad-weight.json"
    corrupted.write_text(json.dumps(payload))
    with pytest.raises(ExactCertificateError, match=r"A\^T lambda"):
        verify_type_i_badcase_farkas(5, corrupted)


def test_corrupted_boolean_eigenvector_is_rejected(tmp_path: Path) -> None:
    source = ROOT / "evidence" / "e1_type_i_badcase_farkas_p7.json"
    payload = json.loads(source.read_text())
    eigen_row = next(row for row in payload["rows"] if "y" in row)
    eigen_row["y"][0] *= -1
    corrupted = tmp_path / "bad-eigenvector.json"
    corrupted.write_text(json.dumps(payload))
    with pytest.raises(ExactCertificateError, match="Boolean eigenvector"):
        verify_type_i_badcase_farkas(7, corrupted)


def test_certificate_domain_is_exactly_the_two_base_primes() -> None:
    with pytest.raises(ValueError, match="only for p=5,7"):
        verify_type_i_badcase_farkas(11)
