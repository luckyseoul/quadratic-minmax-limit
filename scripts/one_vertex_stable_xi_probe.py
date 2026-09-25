#!/usr/bin/env python3
"""Reproduce the stable-spectrum/Xi diagnostic for a finite signing.

Default input is the exact order-15 minimizer used in
NOTE_2026-09-25_COORDINATE_VARIATION_ONE_VERTEX.md.

No third-party packages are required.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


def q_and_fields(A: list[list[int]], x: list[int]) -> tuple[int, list[int]]:
    n = len(A)
    fields = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    q = sum(x[i] * fields[i] for i in range(n)) // 2
    return q, fields


def projective_states(n: int):
    # x and -x define the same antipodal class; fix x_0=+1.
    for mask in range(1 << (n - 1)):
        x = [1] * n
        for j in range(1, n):
            if (mask >> (j - 1)) & 1:
                x[j] = -1
        yield x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "matrix",
        nargs="?",
        default="evidence/K15_exact_minimizer_20260911.json",
    )
    args = ap.parse_args()

    data = json.loads(Path(args.matrix).read_text())
    A = data["A"]
    F = int(data["Phi"])
    n = len(A)

    hist: Counter[int] = Counter()
    phi_check = 0
    stable_count = 0

    for x in projective_states(n):
        q, fields = q_and_fields(A, x)
        phi_check = max(phi_check, abs(q))
        w = [x[i] * fields[i] for i in range(n)]

        energies: list[int] = []
        if all(v >= 0 for v in w):
            energies.append(q)
        if all(v <= 0 for v in w):
            energies.append(-q)

        if energies:
            E = max(energies)
            D = F - E
            hist[D] += 1
            stable_count += 1

    if phi_check != F:
        raise SystemExit(f"stored Phi={F}, enumerated Phi={phi_check}")

    r = F * ((1.0 + 1.0 / n) ** 1.5 - 1.0)
    xi = sum(
        count / (D + r) ** 2
        for D, count in hist.items()
        if D + r <= n
    )
    threshold = 1.0 / (18.0 * math.pi)

    out = {
        "matrix": args.matrix,
        "n": n,
        "Phi": F,
        "stable_antipodal_classes": stable_count,
        "deficit_histogram": {str(k): hist[k] for k in sorted(hist)},
        "neutral_increment": r,
        "Xi_neutral": xi,
        "komlos_threshold": threshold,
        "Xi_over_threshold": xi / threshold,
    }

    # Regression for the pinned default witness.
    if n == 15 and F == 27 and Path(args.matrix).name == "K15_exact_minimizer_20260911.json":
        expected = {0: 66, 2: 66, 4: 39, 6: 117, 12: 13, 14: 39}
        assert dict(sorted(hist.items())) == expected
        assert abs(r - 2.7445120988729617) < 1e-12
        assert abs(xi - 14.141428337053116) < 1e-12

    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
