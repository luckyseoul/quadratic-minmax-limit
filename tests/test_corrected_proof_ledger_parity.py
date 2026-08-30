"""Cross-artifact guards for the corrected all-finite proof ledger."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"
CORRECTED_CHAIN = (
    "2503 -> 2219 -> 1744 -> 1481 -> 1368 -> 1228 -> 1215 -> "
    "1213 -> 1020 -> 869 -> 321 -> 19 -> 14 -> 0"
)
CHAIN_STEPS = (
    ("15700", 2503, 284, 2219),
    ("15701", 2219, 475, 1744),
    ("15702", 1744, 263, 1481),
    ("15703", 1481, 113, 1368),
    ("15704", 1368, 140, 1228),
    ("15705", 1228, 13, 1215),
    ("15706", 1215, 2, 1213),
    ("15707", 1213, 193, 1020),
    ("15708", 1020, 151, 869),
    ("15709", 869, 548, 321),
    ("15710", 321, 302, 19),
    ("15711", 19, 5, 14),
    ("15712", 14, 14, 0),
)


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def test_corrected_json_chain_is_contiguous_and_historical_payloads_are_quarantined():
    previous_after = None
    for suffix, before, excluded, after in CHAIN_STEPS:
        canonical = _load(EVIDENCE / f"e1_gmin_m4_prop{suffix}.json")
        assert canonical["record_status"] == "CORRECTED_REPLAY_SUMMARY"
        assert canonical["profile_count_before"] == before
        assert canonical["profiles_excluded_here"] == excluded
        assert canonical["profile_count_after"] == after
        assert canonical["corrected_chain"] == CORRECTED_CHAIN
        assert canonical["boundary_gate_superseded_by"] == "15.721"
        if previous_after is not None:
            assert before == previous_after
        previous_after = after

        historical = ROOT / canonical["historical_payload"]
        assert historical.exists()
        assert _load(historical)["record_status"] == (
            "HISTORICAL_PRE_15.723_PAYLOAD"
        )


def test_partial_15705_block_is_closed_only_at_15709():
    step_15705 = _load(EVIDENCE / "e1_gmin_m4_prop15705.json")
    step_15709 = _load(EVIDENCE / "e1_gmin_m4_prop15709.json")
    assert step_15705["proof_status"] == "PARTIAL"
    assert step_15705["stage_status"] == "OPEN_AT_THIS_STAGE"
    assert step_15705["remaining_slack_sixteen_profiles"] == 74
    assert step_15705["remaining_slack_sixteen_status_here"] == "OPEN_UNTIL_15.709"
    assert step_15709["historical_slack_sixteen_rows_received_from_15_705"] == 74
    assert step_15709["historical_slack_sixteen_rows_excluded_here"] == 74


def test_canonical_documents_do_not_soft_close_retracted_reductions():
    documents = (
        ROOT / "STATUS.md",
        ROOT / "HANDOFF.md",
        ROOT / "README.md",
        ROOT / "solution.md",
        EVIDENCE / "PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    for document in documents:
        text = document.read_text()
        assert "OPEN_RETRACTED_REDUCTION" in text, document
        assert "15.678" in text, document
        assert "15.684" in text, document
        assert "15.721" in text, document

    status = (ROOT / "STATUS.md").read_text()
    solution = (ROOT / "solution.md").read_text()
    assert "**15.678 PROVED" not in status
    assert "**15.684 PROVED" not in status
    assert (
        "## Proposition 15.678 — exceptional p=17 first survivor is impossible"
        not in solution
    )
    assert (
        "## Proposition 15.684 — the p=23 next endpoint reduces to 203 exact profiles"
        not in solution
    )


def test_corrected_chain_is_identical_in_canonical_narratives():
    for relative in (
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    ):
        normalized = " ".join((ROOT / relative).read_text().split())
        assert CORRECTED_CHAIN in normalized, relative

    solution = (ROOT / "solution.md").read_text()
    assert "2503\\to2219\\to1744\\to1481\\to1368" in solution
    assert "15.705 is partial" in solution
    assert "leaves 74 slack-sixteen rows" in solution


def test_global_api_prose_separates_bounded_true_from_global_false():
    status = (ROOT / "STATUS.md").read_text()
    handoff = (ROOT / "HANDOFF.md").read_text()
    solution = (ROOT / "solution.md").read_text()
    audit = (EVIDENCE / "PROPOSITION_DEDUP_AUDIT_2026-08-30.md").read_text()

    combined = "\n".join((status, handoff, solution, audit))
    assert "e1_bounded_residual_split_closed()" in combined
    assert "residual_ii_bounded_even_k_le_4p_minus_2_closed=True" in status
    assert "`residual_ii_full_closed=False`" in solution
    assert "`e1_closed_general()` is now the corrected global predicate" in audit

    assert "`e1_closed_general()` can be `True`" not in audit
    assert "`residual_ii_full_closed` = affine" not in status
    assert "Residual (ii) later closed by ND" not in solution
    assert "`e1_closed_general` True is old incomplete wiring" not in handoff
    assert "Live `e1` is still the old AND" not in handoff


def test_active_evidence_never_serializes_the_obsolete_global_e1_true():
    """Current certificates must agree with the authoritative open gate."""
    stale_fields = (
        '"e1": true',
        '"e1_closed_general": true',
        '"E1": true',
        '"E1_closed": true',
        '"E1_closed_general": true',
        '"L_closed": true',
        '"L_status": "CLOSED"',
    )
    offenders: list[str] = []
    for path in sorted(EVIDENCE.glob("e1*.json")):
        if "historical" in path.name or "retracted" in path.name:
            continue
        text = path.read_text()
        if any(field in text for field in stale_fields):
            offenders.append(path.name)
    assert offenders == []
