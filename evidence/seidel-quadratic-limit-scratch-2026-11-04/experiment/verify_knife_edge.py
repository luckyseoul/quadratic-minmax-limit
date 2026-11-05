#!/usr/bin/env python3
"""Machine check of the order-16 amplification knife-edge theorem."""
from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    # Cubed form of H(28) < 2 H(8) < H(30) with H(n)=F(n)^{2/3}, F(8)=10.
    assert 28**2 < 8 * 10**2 < 30**2
    assert 784 < 800 < 900

    def H(f: int) -> float:
        return f ** (2.0 / 3.0)

    two_H8 = 2 * H(10)
    assert H(28) < two_H8 < H(30)

    # Dichotomy consequences
    assert H(28) <= two_H8  # F(16)=28 => doubling OK
    assert not (H(30) <= two_H8)  # F(16)=30 => doubling fails

    eta30 = (H(30) - two_H8) / 16.0
    assert eta30 > 0.023 and eta30 < 0.024

    payload = {
        "theorem": "H(16)<=2*H(8) iff F(16)=28, given F(16) in {28,30} and F(8)=10",
        "cubed_comparison": [28**2, 8 * 10**2, 30**2],
        "two_H8": two_H8,
        "H28": H(28),
        "H30": H(30),
        "eta8_if_F16_30": eta30,
        "verified": True,
    }
    out = Path(__file__).resolve().parents[1] / "results" / "experiment7_knife_edge.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("KNIFE_EDGE_VERIFIED", json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
