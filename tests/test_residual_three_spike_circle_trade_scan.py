import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from residual_three_spike_circle_trade_scan import scan


def test_circle_trade_scan_detects_known_positive_mismatch():
    row = scan(11, 0, 11 * 11)
    assert row["circle_count_checked"] == 605
    assert row["switchable_circle_count"] == 4
    assert row["mismatch_histogram"] == {"0": 3, "2": 1}
    assert row["positive_mismatch_exists"]
    assert row["all_hits_pair_rechecked"]


def test_circle_trade_scan_keeps_zero_mismatch_distinct():
    row = scan(13, 0, 13 * 13)
    assert row["circle_count_checked"] == 1014
    assert row["switchable_circle_count"] == 8
    assert row["mismatch_histogram"] == {"0": 8}
    assert not row["positive_mismatch_exists"]
    assert row["all_hits_pair_rechecked"]
