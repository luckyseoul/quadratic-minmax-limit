"""Regression checks for the G(23^2,2) deduplication artifact."""

import json
from pathlib import Path

from sporadic_peisert529_exact import (
    CANONICAL_CERTIFICATE_SHA256,
    CONSTRUCTED_CSV_SHA256,
    PUBLISHED_CSV_SHA256,
    exact_report,
)


def test_sporadic_peisert529_exact_constructor_and_witnesses() -> None:
    report = exact_report()
    assert report["status"] == "regularizable_linear_OA_PN_not_residual_ii_closure"
    assert report["common_neighbor_histogram"] == {"131": 264, "132": 264}
    assert report["witnesses"]["+23"]["int16_sha256"] == (
        "19706232e3f181513356b8515cb3897f79b3827d36975be0549ea704a9c260b7"
    )
    assert report["witnesses"]["-23"]["int16_sha256"] == (
        "1d9873cc4139bb507a16837e4bfeebac36906bd19c2270f90ca01f2076efc62e"
    )
    assert report["published_csv"]["sha256"] == PUBLISHED_CSV_SHA256
    assert report["constructed_csv_sha256"] == CONSTRUCTED_CSV_SHA256
    assert report["canonical_certificate_sha256"] == CANONICAL_CERTIFICATE_SHA256
    recorded = json.loads(
        (Path(__file__).parents[1] / "evidence" / "sporadic_peisert529_exact.json")
        .read_text()
    )
    assert recorded == report
