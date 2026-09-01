"""Exact leaf certificates for the ``p=5,7`` Type-I bad-case boxes.

The floating HiGHS runs in Propositions 15.408 and 15.410 were useful for
discovering these supports, but they are not theorem dependencies.  This
module regenerates the Paley conference matrix, checks every stored Boolean
eigenvector, and verifies an integer Farkas contradiction.  It deliberately
imports no proposition module, optimizer, or eigenshell cache.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

from minmax_quadratic import paley_conference_prime_power


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATES = {
    5: ROOT / "evidence" / "e1_type_i_badcase_farkas_p5.json",
    7: ROOT / "evidence" / "e1_type_i_badcase_farkas_p7.json",
}
EXPECTED_SHA256 = {
    5: "f3e8c8a0f85fcaf95bc5b0556eced1a3699735e3bc11b78285bb4d25abd8008a",
    7: "40a6e5156817a421dcc7debe75a55af5749ffd6cdcb12cda78ead75a3d0cc8db",
}
EXPECTED_ROW_KINDS = {
    5: {
        "lower": 133,
        "minus_bad": 36,
        "plus_ge": 22,
        "plus_le": 39,
        "sum_le": 1,
    },
    7: {
        "lower": 692,
        "minus_bad": 258,
        "plus_ge": 121,
        "plus_le": 155,
    },
}


class ExactCertificateError(ArithmeticError):
    """A checked premise or integer certificate identity failed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ExactCertificateError(message)


def _integer_conference(p: int) -> np.ndarray:
    raw = np.asarray(paley_conference_prime_power(p))
    C = raw.astype(np.int64)
    n = p * p + 1
    _require(np.array_equal(raw, C), "conference matrix is not integral")
    _require(C.shape == (n, n), "conference matrix shape changed")
    _require(np.array_equal(C, C.T), "conference matrix is not symmetric")
    _require(np.count_nonzero(np.diag(C)) == 0, "conference diagonal changed")
    _require(
        set(np.unique(C)).issubset({-1, 0, 1}),
        "conference alphabet changed",
    )
    _require(
        np.array_equal(C @ C, p * p * np.eye(n, dtype=np.int64)),
        "conference square identity changed",
    )
    return C


def _verify(
    p: int,
    path: Path,
    expected_digest: str | None,
) -> dict[str, object]:
    try:
        payload = path.read_bytes()
        cert = json.loads(payload)
    except (OSError, json.JSONDecodeError) as exc:
        raise ExactCertificateError(
            f"cannot load exact certificate {path}: {exc}"
        ) from exc

    digest = hashlib.sha256(payload).hexdigest()
    if expected_digest is not None:
        _require(digest == expected_digest, "certificate digest changed")
    _require(type(p) is int and p in CERTIFICATES, "only p=5,7 are certified")
    _require(cert.get("p") == p, "certificate prime mismatch")
    _require(
        cert.get("eliminated_edge_index") is None,
        "tracked certificate must check the full edge cone",
    )
    rows = cert.get("rows")
    _require(isinstance(rows, list), "certificate rows are missing")
    _require(cert.get("row_count") == len(rows), "certificate row count changed")

    C = _integer_conference(p)
    n = len(C)
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    _require(edges[0] == (0, 1), "distinguished-edge ordering changed")
    edge_a = np.asarray([a for a, _b in edges], dtype=np.int64)
    edge_b = np.asarray([b for _a, b in edges], dtype=np.int64)
    edge_c = np.asarray([C[a, b] for a, b in edges], dtype=np.int64)

    lhs = np.zeros(len(edges), dtype=object)
    rhs = 0
    kinds: Counter[str] = Counter()
    eigenvectors: set[tuple[int, ...]] = set()

    for row_number, record in enumerate(rows):
        _require(isinstance(record, dict), f"row {row_number} is not an object")
        kind = record.get("kind")
        weight = record.get("weight")
        _require(
            type(weight) is int and weight > 0,
            f"row {row_number} has invalid weight",
        )
        kinds[kind] += 1

        feature: np.ndarray | None = None
        if kind in {"plus_le", "plus_ge", "minus_bad"}:
            raw_y = record.get("y")
            _require(
                isinstance(raw_y, list) and len(raw_y) == n,
                f"row {row_number} has invalid eigenvector length",
            )
            _require(
                all(type(value) is int and value in (-1, 1) for value in raw_y),
                f"row {row_number} eigenvector is not Boolean",
            )
            y_tuple = tuple(raw_y)
            y = np.asarray(y_tuple, dtype=np.int64)
            eigenvectors.add(y_tuple)
            eigenvalue = p if kind.startswith("plus_") else -p
            _require(
                np.array_equal(C @ y, eigenvalue * y),
                f"row {row_number} is not a {eigenvalue:+d} Boolean eigenvector",
            )
            feature = edge_c * y[edge_a] * y[edge_b]

        if kind == "plus_le":
            _require(feature is not None, "plus row has no feature vector")
            A = feature.astype(object)
            b = int(3 - 2 * feature[0])
        elif kind == "plus_ge":
            _require(feature is not None, "plus row has no feature vector")
            A = (-feature).astype(object)
            b = int(-(3 - 2 * feature[0]))
        elif kind == "sum_le":
            A = np.ones(len(edges), dtype=object)
            b = 3 * p - 2
        elif kind == "minus_bad":
            _require(feature is not None, "minus row has no feature vector")
            A = feature.astype(object)
            b = min(-1, int(-3 * feature[0]))
        elif kind == "lower":
            index = record.get("index")
            _require(
                type(index) is int and 0 <= index < len(edges),
                f"row {row_number} has invalid lower-bound index",
            )
            A = np.zeros(len(edges), dtype=object)
            A[index] = -1
            b = 0
        else:
            raise ExactCertificateError(
                f"row {row_number} has unknown kind {kind!r}"
            )

        lhs += weight * A
        rhs += weight * b

    bad_coordinates = [i for i, value in enumerate(lhs) if value != 0]
    _require(
        not bad_coordinates,
        f"A^T lambda is nonzero at coordinates {bad_coordinates[:10]}",
    )
    _require(
        type(cert.get("rhs")) is int and rhs == cert["rhs"],
        "certificate RHS changed",
    )
    _require(rhs < 0, "Farkas RHS is not strictly negative")
    row_kinds = dict(sorted(kinds.items()))
    _require(row_kinds == EXPECTED_ROW_KINDS[p], "certificate support changed")
    return {
        "proved": True,
        "p": p,
        "variables": len(edges),
        "rows": len(rows),
        "distinct_eigenvectors": len(eigenvectors),
        "row_kinds": row_kinds,
        "exact_lhs_nonzeros": len(bad_coordinates),
        "exact_rhs": rhs,
        "sha256": digest,
        "uses_scipy": False,
        "uses_eigenshell_cache": False,
        "uses_upper_bounds": False,
    }


@lru_cache(maxsize=2)
def _verify_default(p: int) -> dict[str, object]:
    return _verify(p, CERTIFICATES[p], EXPECTED_SHA256[p])


def verify_type_i_badcase_farkas(
    p: int,
    certificate_path: str | Path | None = None,
) -> dict[str, object]:
    """Verify and expose one exact small-prime bad-case exclusion."""
    if type(p) is not int or p not in CERTIFICATES:
        raise ValueError("exact Type-I base certificates exist only for p=5,7")
    if certificate_path is None:
        return dict(_verify_default(p))
    return _verify(p, Path(certificate_path), None)


def type_I_badcase_small_primes_exact() -> bool:
    """Return true exactly when both tracked integer identities verify."""
    return all(verify_type_i_badcase_farkas(p)["proved"] for p in (5, 7))


def type_i_badcase_small_primes_exact() -> bool:
    """PEP-8 alias for :func:`type_I_badcase_small_primes_exact`."""
    return type_I_badcase_small_primes_exact()
