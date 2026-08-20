import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_small_exceptional_shell_probe_certificate():
    data = json.loads(
        (ROOT / "evidence" / "exceptional_projection_shell_probe.json").read_text()
    )
    assert data["pointwise_exceptional_floor_false"] is True
    assert data["nonzero_shell_ge_3n_false"] is True
    assert data["small_prime_full_shells"]["7"]["shells"] == [
        {"norm_sq": "0", "count": 2352},
        {"norm_sq": "928/3", "count": 8400},
        {"norm_sq": "608", "count": 700},
    ]
    witness = data["p11_nonzero_below_3n_witness"]
    assert Fraction(witness["norm_sq"]) == Fraction(4304, 15)
    assert Fraction(witness["norm_sq"]) < witness["target_3n"]


def test_full_p11_exceptional_shell_certificate():
    data = json.loads(
        (
            ROOT
            / "evidence"
            / "maxplus_p11"
            / "exceptional_projection_shells_p11_xpu.json"
        ).read_text()
    )
    histogram = {Fraction(value): count for value, count in data["histogram"].items()}
    assert data["row_count_eps_plus"] == 37_457_112
    assert len(histogram) == 37
    assert data["zero_count"] == 0
    assert min(histogram) == Fraction(4304, 15)
    assert histogram[Fraction(4304, 15)] == 442_860
    assert sum(histogram.values()) == data["row_count_eps_plus"]
    mean = sum(value * count for value, count in histogram.items()) / sum(
        histogram.values()
    )
    assert mean == Fraction(74_989_008, 141_883)
    assert mean > data["target_3n"]

    assert set(data["strata"]) == {"k1", "k3", "k4", "k5", "k6"}
    assert sum(part["count"] for part in data["strata"].values()) == data[
        "row_count_eps_plus"
    ]
    for part in data["strata"].values():
        assert sum(part["histogram"].values()) == part["count"]
