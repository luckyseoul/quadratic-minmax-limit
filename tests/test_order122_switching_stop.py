import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEARCH = ROOT / "evidence" / "order122_switching"


def test_order122_switching_sources_match_the_pinned_hashes():
    expected = {
        "gm121_scan.py": "a9001052233f6b9d0aef67de1e58263e660dda234608bb126b342c2861872df8",
        "gm_oa121.py": "4a2d33a5cdff484a9a54d5faa25b54961f67b57200854fda06494e791ba4fa24",
        "peisert121_exact.py": "2ae21c03bbe1421fe019ae714fc3e895c2b8e4544a50358044927a1f5d2e639c",
        "pn_wqh3_search.cpp": "999fb0b45d18487adfc5712134d4a45bf541c5e25b770beada399a41524d4d64",
        "pn_wqh4_linesubsets.cpp": "747fcc05e6e92499d585eef33dd433ab94cdf654d1db1f9b22f986ec54af2b51",
        "pn_wqh_lines.cpp": "6585819a7c9642b4da0989784cb7e39bcfc7d8dafa4c310c64091d0bb5af3e5c",
        "wqh121_general.py": "1df403a9be9d088439cf632da03f03cccc4f66cdd71e1bf8ff3169b30b105a19",
    }
    observed = {
        name: hashlib.sha256((SEARCH / name).read_bytes()).hexdigest()
        for name in expected
    }
    assert observed == expected


def test_stopping_note_records_exact_scope_and_nonclosure():
    text = (ROOT / "evidence" / "NOTE_2026-09-03_ORDER122_SWITCHING_STOP.md").read_text()
    flat = " ".join(text.split())
    for row in (
        "NONE tested=7988904 wqh_nontrivial=0",
        "The total is 6,534,660 exactly checked candidate pairs.",
        "mask=63  exact_candidates=8  NONE",
        "mask=95  exact_candidates=0  NONE",
        "mask=111 exact_candidates=0  NONE",
        "mask=119 exact_candidates=1  NONE",
    ):
        assert row in text
    assert "entire line-supported WQH `(4,4)` ansatz is empty" in flat
    assert "All 40 cells returned `INFEASIBLE`" in text
    assert "do not classify" in flat
    assert "remain OPEN" in flat
