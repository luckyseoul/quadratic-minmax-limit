"""Drive e1_main_chain_status: L OPEN, docs_ok, no soft-close."""
from __future__ import annotations

from e1_main_chain_status import check_docs_L_status, run_main_chain


def test_main_chain_L_open_and_docs_ok():
    out = run_main_chain()
    assert out["L_status"] == "OPEN"
    assert out["L_closed"] is False
    assert out["gsum_disj_lb_proved_general"] is False
    assert out["writeup_L_closed"] is False
    docs = out["docs"]
    assert docs["soft_close_detected"] is False
    assert docs["overclaim_detected"] is False
    assert docs["HANDOFF_shows_L_OPEN"] is True
    assert docs["docs_ok"] is True
    units = docs["four_e1_units"]
    assert units["phi_F_ge_6"] is False
    assert units["residual_ii_k_ge_4p"] is False
    assert units["type_I_multilevel"] is False
    assert units["closed"] is False


def test_solution_does_not_assert_limit_theorem():
    from pathlib import Path

    text = (Path(__file__).resolve().parents[1] / "solution.md").read_text(
        encoding="utf-8", errors="replace"
    )
    head = text[:4000]
    assert "**Main Theorem (limit).**" not in head
    assert "E(1) and" in head and "are not complete" in head
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text(
        encoding="utf-8", errors="replace"
    )
    assert "blocked**\nby residual **(i)** only" not in readme
    assert "L=\\lim_n\\alpha_n$ is **OPEN**" in readme or "is **OPEN**" in readme


def test_soft_close_detector_flags_bare_L_CLOSED():
    """Unit-level: detector pattern must catch bare '**L CLOSED.**' style."""
    import re

    pat = re.compile(r"\*\*L CLOSED\.\*\*|\*\*L CLOSED\*\*|E\(1\)\s*\(CLOSED|E1_closed_general=true\. L CLOSED", re.I)
    assert pat.search("**E(1) (CLOSED, 15.168–171):** ... **L CLOSED.**")
    assert pat.search("E1_closed_general=true. L CLOSED (via E1∧bi-tight).")
    assert not pat.search("**L OPEN.** Denseness path blocked")


def test_props_15170_171_body_not_soft_closed():
    """solution.md Props 15.170–171 must not claim residual/E1/L CLOSED while hinge open."""
    from pathlib import Path
    import re

    sol = Path(__file__).resolve().parents[1] / "solution.md"
    text = sol.read_text(encoding="utf-8", errors="replace")
    # extract from Prop 15.170 onward
    i = text.find("## Prop 15.170")
    assert i >= 0
    tail = text[i:]
    # forbidden soft-close phrases (post-retraction)
    bad = [
        r"E1_closed_general=true\. L CLOSED",
        r"\*\*E\(1\) closed\*\*",
        r"Closes residual \(i\) of E\(1\) for all primes",
        r"Closes residual \(ii\) of E\(1\) for all primes",
        r"association-scheme min \$-12",
        r"disj: association-scheme min",
    ]
    for pat in bad:
        assert not re.search(pat, tail), f"soft-close residue matched: {pat}"
    assert "L OPEN" in tail or "OPEN for general" in tail or "not complete" in text[:4000]
    assert "gsum_disj_lb_proved_general" in tail or "NOT proved for general" in tail
