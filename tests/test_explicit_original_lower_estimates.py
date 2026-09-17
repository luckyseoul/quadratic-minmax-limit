"""Checks for the 2026-09-15 explicit lower-bound verification artifacts.

These are cheap deterministic predicates: the Paley builder, the exhaustive
order-4 estimate checks (including the update bound), the exact conference
frame identity, the exact-rational certificates of the sharpened constant,
and the consistency of ``result.json`` with the files it pins.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EVID = ROOT / "evidence" / "explicit_original_lower_20260915"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


verify = _load("explicit_lower_verify_estimates", EVID / "verify_estimates.py")
sharpen = _load("explicit_lower_sharpen_constant", EVID / "sharpen_constant.py")


def test_paley_builder_is_a_conference_signing():
    for q in (5, 13, 29, 101):
        C = verify.paley_conference(q)
        assert np.abs(C @ C - q * np.eye(q + 1)).max() == 0.0
        assert np.all(np.diag(C) == 0)
        iu = np.triu_indices(q + 1, 1)
        assert np.all(np.abs(C[iu]) == 1)


def test_exhaustive_order_four_estimates_hold():
    agg: dict = {}
    for idx in range(1 << verify.num_entries(4)):
        res = verify.analyze(verify.decode_signing(4, idx), update_check=True)
        verify.merge_check_agg(agg, res, {"idx": idx})
    assert agg, "no checks ran"
    for key, slot in agg.items():
        assert slot["violations"] == 0, key


def test_conference_frame_is_exact_and_estimates_active():
    C = verify.paley_conference(13)
    res = verify.analyze(C)
    # For exact conference signings R = I + S identically (float64: ~2e-15).
    assert res["checks"]["c5"]["lhs"] < 1e-12
    assert res["checks"]["c5"]["ok"]
    assert res["checks"]["c4"]["ok"]
    assert res["checks"]["cX"]["ok"]
    # The (6) estimate is active near saturation and holds with slack.
    assert res["checks"]["c6"]["rhs"] > 0.0
    assert res["checks"]["c6"]["ok"]


def test_sharpen_certificates_and_ceiling():
    for eps in (Fraction(1, 10 ** 6), Fraction(4, 10 ** 6)):
        cert = sharpen.certify(eps)
        assert cert["contradiction"] and cert["d_bound_ok"]
    assert not sharpen.certify(Fraction(9, 2 * 10 ** 6))["contradiction"]
    assert sharpen.note_exact_reproduction()["matches_published"]
    ceiling = Fraction(sharpen.ceiling()["ceiling_eps"])
    assert 4 * 10 ** -6 < float(ceiling) < 4.6e-6
    assert sharpen.certify(ceiling)["contradiction"]


def test_result_receipt_is_consistent_with_pinned_files():
    receipt = json.loads((EVID / "result.json").read_text())
    for key in ("note", "review", "readme"):
        pinned = receipt[key]
        digest = hashlib.sha256((ROOT / pinned["path"]).read_bytes()).hexdigest()
        assert digest == pinned["sha256"], key
    for name, digest in receipt["scripts"].items():
        assert hashlib.sha256((EVID / name).read_bytes()).hexdigest() == digest
    assert receipt["violations_total"] == 0
    assert receipt["convergence_proved"] is False
