import hashlib
import json
from pathlib import Path

from e1_gmin_m4_p31_equi_zero68_mitm import (
    p31_equianharmonic_alignment_census,
    p31_equianharmonic_zero68_mitm_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "evidence" / "p31_equianharmonic_zero68_mitm.cpp"
MANIFEST = (
    ROOT / "evidence" / "p31_equianharmonic_zero68_mitm_manifest.json"
)
SOURCE_SHA256 = "14e23138797b8bde9edbbd447c69ee735ec0a85600145b2ecf0988e490c59520"


def test_complete_atom_alignment_census_forces_total_deficit_three():
    row = p31_equianharmonic_alignment_census()
    assert row["ae_atom_count"] == 4_495
    assert row["compact_atom_count"] == 13_485
    assert row["ae_score_counts"] == {
        -3: 9,
        -2: 1,
        -1: 702,
        0: 3_071,
        1: 702,
        2: 1,
        3: 9,
    }
    assert row["compact_score_counts"] == {
        -2: 111,
        -1: 2_133,
        0: 8_997,
        1: 2_133,
        2: 111,
    }
    assert row["maximum_compact_score"] == 2
    assert row["maximal_ae_cycle_count"] == 9
    assert row["maximal_ae_cycles_are_edge_disjoint"]
    assert row["maximal_ae_target_edges_covered"] == 27
    assert row["broken_target_edges"] == 2
    assert row["maximal_compact_antipodal_off_edge_count"] == 5
    assert row["maximal_compact_supported_off_orbit_count"] == 51


def test_exact_partition_totals_exclude_simultaneous_degree_six_and_eight():
    row = p31_equianharmonic_zero68_mitm_certificate()
    assert row["maximum_total_score"] == 32
    assert row["required_target_score"] == 29
    assert row["total_deficit"] == 3
    assert row["deficit_partitions"] == ((3,), (2, 1), (1, 1, 1))
    assert row["partitions"]["deficit_3"]["completions"] == 13_528_344
    assert (
        row["partitions"]["deficit_2_plus_1"]["exceptional_pairs"]
        == 20_697_666
    )
    assert (
        row["partitions"]["three_deficit_1_with_ae"][
            "exceptional_multisets"
        ]
        == 2_278_045
    )
    all_compact = row["partitions"]["three_compact_deficit_1"]
    assert all_compact["all_unordered_multisets"] == 1_619_689_995
    assert all_compact["unsupported_projection_3sum_hits"] == 2_027_542
    assert all_compact["completions"] == 108_480_057
    assert row["maximal_completions_tested"] == 230_314_710
    assert row["edge_hits"] == 17_076
    assert row["zero_degree_six_and_eight_hits"] == 0
    assert row["p31_b7_k11_constant_conic_zero68_fiber"] == "UNSAT"
    assert row["proved"]
    assert not row["residual_ii_closed"]
    assert row["L_status"] == "OPEN"


def test_machine_manifest_and_exhaustive_source_are_pinned():
    manifest = json.loads(MANIFEST.read_text())
    assert manifest["status"] == "UNSAT"
    assert manifest["scope"].startswith("p=31, b=7, k=11")
    assert manifest["totals"] == {
        "maximal_completions_tested": 230_314_710,
        "edge_hits": 17_076,
        "zero_6_8_hits": 0,
    }
    assert not manifest["global_claims"]["residual_ii_closed"]
    assert not manifest["global_claims"]["L_closed"]
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert manifest["reproduction"]["source_sha256"] == SOURCE_SHA256
