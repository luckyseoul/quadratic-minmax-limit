"""Drive e1_main_chain_status: L OPEN, docs_ok, no soft-close."""
from __future__ import annotations

from e1_main_chain_status import check_docs_L_status, run_main_chain


def test_main_chain_L_open_and_docs_ok():
    out = run_main_chain()
    assert out["L_status"] == "OPEN"
    assert out["L_closed"] is False
    assert out["E1_closed"] is False
    assert out["gsum_disj_lb_proved_general"] is False
    docs = out["docs"]
    assert docs["soft_close_detected"] is False
    assert docs["HANDOFF_shows_L_OPEN"] is True
    assert docs["docs_ok"] is True


def test_soft_close_detector_flags_bare_L_CLOSED():
    """Unit-level: detector pattern must catch bare '**L CLOSED.**' style."""
    import re

    pat = re.compile(r"\*\*L CLOSED\.\*\*|\*\*L CLOSED\*\*|E\(1\)\s*\(CLOSED", re.I)
    assert pat.search("**E(1) (CLOSED, 15.168–171):** ... **L CLOSED.**")
    assert not pat.search("**L OPEN.** Denseness path blocked")
