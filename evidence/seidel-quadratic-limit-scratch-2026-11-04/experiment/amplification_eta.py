#!/usr/bin/env python3
"""CORE §7 amplification eta census on exact F(n) + F(16) branches."""
from __future__ import annotations

import json
from pathlib import Path

F = [0, 1, 3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21, 27]  # n=1..15


def H(f: float) -> float:
    return f ** (2.0 / 3.0)


def eta_double(n: int, f_n: int, f_2n: int) -> float:
    # H(2n) <= 2 H(n) + 2n eta  => eta >= max(0, (H(2n)-2H(n))/(2n))
    return max(0.0, (H(f_2n) - 2 * H(f_n)) / (2 * n))


def eta_triple(n: int, f_n: int, f_3n: int) -> float:
    return max(0.0, (H(f_3n) - 3 * H(f_n)) / (3 * n))


def main() -> None:
    rows = []
    for n in range(2, 8):
        if 2 * n <= 15:
            e2 = eta_double(n, F[n - 1], F[2 * n - 1])
            rows.append({"n": n, "kind": "double", "m": 2 * n, "F_n": F[n - 1], "F_m": F[2 * n - 1], "eta": e2, "sub": e2 == 0.0})
        if 3 * n <= 15:
            e3 = eta_triple(n, F[n - 1], F[3 * n - 1])
            rows.append({"n": n, "kind": "triple", "m": 3 * n, "F_n": F[n - 1], "F_m": F[3 * n - 1], "eta": e3, "sub": e3 == 0.0})

    branches = {}
    for f16 in (28, 30):
        e = eta_double(8, 10, f16)
        branches[str(f16)] = {
            "F16": f16,
            "eta8_double": e,
            "sub_neutral_H": e == 0.0,
            "H8": H(10),
            "H16": H(f16),
            "two_H8": 2 * H(10),
        }

    # Sum of eta* over exact double/triple rows (crude finite E)
    etas = [r["eta"] for r in rows]
    payload = {
        "hypothesis": "CORE §7 eta vanishes on several exact doublings; F(16)=28 makes eta_8=0 while F(16)=30 forces eta_8>0",
        "exact_multipliers": rows,
        "f16_branches": branches,
        "num_exact_zero_eta": sum(1 for r in rows if r["sub"]),
        "max_exact_eta": max(etas) if etas else None,
        "prediction": {
            "some_zero_eta": any(r["sub"] for r in rows),
            "f28_eta8_zero": branches["28"]["sub_neutral_H"],
            "f30_eta8_positive": branches["30"]["eta8_double"] > 0,
        },
    }
    payload["holds"] = (
        payload["prediction"]["some_zero_eta"]
        and payload["prediction"]["f28_eta8_zero"]
        and payload["prediction"]["f30_eta8_positive"]
    )
    out = Path(__file__).resolve().parents[1] / "results" / "experiment6_amplification_eta.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload["prediction"], indent=2))
    print("holds=", payload["holds"])
    print("f16 branches", json.dumps(branches, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
