"""Route-neutral registry/wiring tests; synthetic entries are NOT MO proofs."""
from dataclasses import replace
from hashlib import sha256
import importlib

import pytest

import original_mo_status as registry
from original_mo_status import Conclusion, ReviewedCompletionProof


def _synthetic_entry(tmp_path, conclusion, suffix="one", limit_value=None):
    """Create pinned plumbing fixtures, explicitly never production evidence."""
    evidence = tmp_path / "evidence"
    evidence.mkdir(exist_ok=True)
    theorem = f"SYNTHETIC TEST ONLY: no mathematical theorem. Fixture {suffix}.\n"
    review = f"SYNTHETIC TEST ONLY: no mathematical review. Fixture {suffix}.\n"
    theorem_path = f"evidence/synthetic-{suffix}-theorem.md"
    review_path = f"evidence/synthetic-{suffix}-review.md"
    (tmp_path / theorem_path).write_text(theorem, encoding="utf-8")
    (tmp_path / review_path).write_text(review, encoding="utf-8")
    return ReviewedCompletionProof(
        proof_id=f"synthetic-test-only-{suffix}",
        conclusion=conclusion,
        theorem="Synthetic fixture; NOT a proof of MathOverflow 413935",
        theorem_path=theorem_path,
        theorem_sha256=sha256(theorem.encode()).hexdigest(),
        review_path=review_path,
        review_sha256=sha256(review.encode()).hexdigest(),
        limit_value=limit_value,
    )


def _install_entries(monkeypatch, tmp_path, *entries):
    monkeypatch.setattr(registry, "ROOT", tmp_path)
    monkeypatch.setattr(registry, "_reviewed_completion_entries", lambda: entries)


def _synthetic_docs(monkeypatch, tmp_path, declaration):
    import e1_main_chain_status as chain

    monkeypatch.setattr(chain, "ROOT", tmp_path)
    for name in ("HANDOFF.md", "STATUS.md", "solution.md"):
        (tmp_path / name).write_text(declaration + "\n", encoding="utf-8")
    package = tmp_path / "evidence" / "share"
    package.mkdir(parents=True, exist_ok=True)
    (package / "denseness_path_package.md").write_text(
        "Historical optional Paley route: open.\n", encoding="utf-8"
    )


def _stub_paley(monkeypatch, closed):
    """Avoid executing historical theorem checkers or artifact-writing mains."""
    import e1_main_chain_status as chain

    monkeypatch.setattr(chain, "four_e1_units_closed", lambda: {"closed": closed})
    monkeypatch.setattr(chain, "e1_closed_general", lambda: closed)
    monkeypatch.setattr(chain, "gsum_disj_lb_proved_general", lambda: False)
    monkeypatch.setattr(
        chain, "run_bitight_chain",
        lambda: {"bi_tight_required_levels_empty_for_all_p_ge_5": closed},
    )
    monkeypatch.setattr(chain, "prop15168_main", lambda: {"proved": {}})
    monkeypatch.setattr(chain, "e1_residual_open", lambda: {"open": [] if closed else ["Paley"]})


def test_production_registry_is_empty_and_route_neutral():
    assert registry._reviewed_completion_entries() == ()
    status = registry.original_mo_status()
    assert status["registry_valid"] is True
    assert status["status"] == status["limit_status"] == "OPEN"
    assert status["problem_settled"] is False
    assert status["existence_proved"] is False
    assert status["nonexistence_proved"] is False
    assert status["value_proved"] is False
    assert status["limit_value"] is None
    assert status["reviewed_completion_proofs"] == []
    assert status["required_optional_routes"] == []


def test_status_has_no_caller_supplied_proof_boolean():
    with pytest.raises(TypeError):
        registry.original_mo_status(proved=True)


@pytest.mark.parametrize(
    "conclusion,value,result,existence,nonexistence",
    [
        (Conclusion.EXISTENCE, None, "EXISTENCE_PROVED", True, False),
        (Conclusion.VALUE, "1/2", "VALUE_PROVED", True, False),
        (Conclusion.NONEXISTENCE, None, "NONEXISTENCE_PROVED", False, True),
    ],
)
def test_distinct_reviewed_conclusions(monkeypatch, tmp_path, conclusion, value,
                                      result, existence, nonexistence):
    entry = _synthetic_entry(tmp_path, conclusion, limit_value=value)
    _install_entries(monkeypatch, tmp_path, entry)
    status = registry.original_mo_status()
    assert status["registry_valid"] is True
    assert status["status"] == result
    assert status["problem_settled"] is True
    assert status["existence_proved"] is existence
    assert status["nonexistence_proved"] is nonexistence
    assert status["value_proved"] is (value is not None)
    assert status["limit_value"] == value
    assert status["reviewed_completion_proofs"][0]["proof_id"] == entry.proof_id


@pytest.mark.parametrize("bad", [True, {"proved": True}, "PROVED"])
def test_unreviewed_flags_cannot_enter_registry(monkeypatch, tmp_path, bad):
    _install_entries(monkeypatch, tmp_path, bad)
    status = registry.original_mo_status()
    assert status["registry_valid"] is False
    assert status["problem_settled"] is False
    assert status["registry_errors"]


@pytest.mark.parametrize(
    "changes",
    [
        {"proof_id": ""},
        {"problem": "a finite signing, not the original all-orders sequence"},
        {"conclusion": "existence"},
        {"theorem": ""},
        {"limit_value": "1/2"},
        {"review_path": "evidence/missing-review.md"},
        {"theorem_path": "../outside.md"},
        {"review_sha256": "not-a-hash"},
        {"theorem_sha256": "0" * 64},
    ],
)
def test_incomplete_or_stale_entries_fail_closed(monkeypatch, tmp_path, changes):
    entry = _synthetic_entry(tmp_path, Conclusion.EXISTENCE)
    _install_entries(monkeypatch, tmp_path, replace(entry, **changes))
    status = registry.original_mo_status()
    assert status["registry_valid"] is False
    assert status["problem_settled"] is False
    assert status["limit_status"] == "OPEN"


def test_value_entry_requires_value(monkeypatch, tmp_path):
    entry = _synthetic_entry(tmp_path, Conclusion.VALUE)
    _install_entries(monkeypatch, tmp_path, entry)
    assert registry.original_mo_status()["registry_valid"] is False


def test_review_must_be_separate_pinned_artifact(monkeypatch, tmp_path):
    entry = _synthetic_entry(tmp_path, Conclusion.EXISTENCE)
    entry = replace(entry, review_path=entry.theorem_path, review_sha256=entry.theorem_sha256)
    _install_entries(monkeypatch, tmp_path, entry)
    assert registry.original_mo_status()["registry_valid"] is False


@pytest.mark.parametrize("case", ["duplicate_id", "existence_conflict", "value_conflict"])
def test_conflicting_registry_fails_closed(monkeypatch, tmp_path, case):
    if case == "value_conflict":
        first = _synthetic_entry(tmp_path, Conclusion.VALUE, "one", "1/2")
        second = _synthetic_entry(tmp_path, Conclusion.VALUE, "two", "1/3")
    else:
        first = _synthetic_entry(tmp_path, Conclusion.EXISTENCE, "one")
        second = (
            first if case == "duplicate_id"
            else _synthetic_entry(tmp_path, Conclusion.NONEXISTENCE, "two")
        )
    _install_entries(monkeypatch, tmp_path, first, second)
    status = registry.original_mo_status()
    assert status["registry_valid"] is False
    assert status["problem_settled"] is False


@pytest.mark.parametrize("declaration", ["L OPEN", "**L CLOSED.**"])
def test_paley_flags_alone_cannot_close_original_problem(monkeypatch, tmp_path, declaration):
    import e1_main_chain_status as chain

    _install_entries(monkeypatch, tmp_path)
    _synthetic_docs(monkeypatch, tmp_path, declaration)
    _stub_paley(monkeypatch, True)
    out = chain.run_main_chain()
    assert out["optional_route_diagnostics"]["paley_gap_two"]["closed"] is True
    assert out["E1_closed"] is True  # Preserved route-local diagnostic.
    assert out["Main_closed"] is False
    assert out["L_closed"] is False
    assert out["L_status"] == "OPEN"
    assert out["writeup_L_closed"] is False
    assert out["docs"]["docs_ok"] is (declaration == "L OPEN")


@pytest.mark.parametrize(
    "conclusion,value,declaration",
    [
        (Conclusion.EXISTENCE, None, "Existence PROVED; value unidentified."),
        (Conclusion.VALUE, "1/2", "L=1/2 CLOSED"),
        (Conclusion.NONEXISTENCE, None, "Nonexistence PROVED"),
        (Conclusion.NONEXISTENCE, None, "Non-existence CLOSED"),
        (Conclusion.NONEXISTENCE, None, "Non existence PROVED"),
    ],
)
def test_open_paley_route_cannot_veto_reviewed_global_proof(
    monkeypatch, tmp_path, conclusion, value, declaration,
):
    import e1_main_chain_status as chain

    entry = _synthetic_entry(tmp_path, conclusion, limit_value=value)
    _install_entries(monkeypatch, tmp_path, entry)
    _synthetic_docs(monkeypatch, tmp_path, declaration)
    _stub_paley(monkeypatch, False)
    out = chain.run_main_chain()
    assert out["optional_route_diagnostics"]["paley_gap_two"]["closed"] is False
    assert out["E1_closed"] is False
    assert out["Main_closed"] is True
    assert out["L_closed"] is (conclusion is not Conclusion.NONEXISTENCE)
    assert out["L_status"] == ("NONEXISTENT" if conclusion is Conclusion.NONEXISTENCE else "CLOSED")
    assert out["L_value"] == value
    assert out["value_proved"] is (value is not None)
    assert out["docs"]["docs_ok"] is True
    assert out["docs"]["global_claim_hits"] == []


def test_existence_proof_does_not_license_value_claim(monkeypatch, tmp_path):
    import e1_main_chain_status as chain

    _install_entries(monkeypatch, tmp_path, _synthetic_entry(tmp_path, Conclusion.EXISTENCE))
    _synthetic_docs(monkeypatch, tmp_path, "L=1/2 CLOSED")
    _stub_paley(monkeypatch, False)
    docs = chain.check_docs_L_status()
    assert docs["original_mo"]["existence_proved"] is True
    assert docs["original_mo"]["value_proved"] is False
    assert docs["docs_ok"] is False
    assert docs["global_claim_hits"]


@pytest.mark.parametrize("conclusion", [None, Conclusion.EXISTENCE, Conclusion.NONEXISTENCE])
def test_three_legacy_live_status_wrappers_follow_registry(monkeypatch, tmp_path, conclusion):
    import e1_main_chain_status as chain

    entries = () if conclusion is None else (_synthetic_entry(tmp_path, conclusion),)
    _install_entries(monkeypatch, tmp_path, *entries)
    monkeypatch.setattr(chain, "four_e1_units_closed", lambda: {"closed": conclusion is None})
    expected = registry.original_mo_status()["limit_status"]
    for name in (
        "e1_gmin_leftover1_qvar_principal", "e1_gmin_global_qvar", "e1_gmin_qvar_k_ge_7",
    ):
        assert importlib.import_module(name).live_L_status() == expected


@pytest.mark.parametrize("number", [15168, 15169, 15170, 15171])
def test_legacy_implication_alias_is_explicitly_route_local(number):
    module = importlib.import_module(f"e1_gmin_m4_prop{number}")
    for e1, bitight in ((False, False), (False, True), (True, False), (True, True)):
        out = module.main_L_from_e1(e1, bitight)
        assert out["scope"] == "optional_paley_route"
        assert out["sufficient_not_necessary"] is True
        assert out["paley_route_closed"] is (e1 and bitight)
        assert out["E1"] is e1 and out["bi_tight"] is bitight
        assert out["L_closed"] is (e1 and bitight)  # Existing route-local API.
        assert out["L_status"] == ("CLOSED" if e1 and bitight else "OPEN")
        assert "not a necessity claim" in out["rule"]
    assert registry.original_mo_status()["problem_settled"] is False
