#!/usr/bin/env python3
"""Search nsq-circle flips of named z that stay Max-."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15613 import named_z, _finv  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def search(p):
    z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
    Cmat = paley_conference_prime_power(p)
    nsq_b = [e for e in range(1, q) if chi(e) == -1]
    # unique directions
    used = set()
    dirs = []
    for b in nsq_b:
        if b in used:
            continue
        dirs.append(b)
        for t in range(1, p):
            used.add(mul(t, b))
    n_ok = 0
    examples = []
    zf = z.astype(np.float64)
    for b in dirs:
        for s in range(p):  # parallels: s + F_p b, skip those through 0
            L = [add(s, mul(t, b)) for t in range(p)]
            if 0 in L:
                continue
            IL = [_finv(mul, q, x) for x in L]
            negIL = {((p - x % p) % p) + ((p - x // p) % p) * p for x in IL}
            for a in range(1, min(q, p * 3)):  # cap a for p=11
                if a in negIL:
                    continue
                circle = [a] + [add(a, x) for x in IL]
                if 0 in circle:
                    continue
                z2 = z.copy()
                for x in circle:
                    z2[1 + x] *= -1
                if np.allclose(Cmat @ z2.astype(np.float64), -p * z2.astype(np.float64), atol=1e-5):
                    n_ok += 1
                    if len(examples) < 3:
                        examples.append({"b": int(b), "s": int(s), "a": int(a), "chi_b": int(chi(b))})
                    break
            if n_ok >= 3:
                break
        if n_ok >= 3:
            break
    print(f"p={p} nsq_dirs={len(dirs)} n_ok_found={n_ok} ex={examples}", flush=True)
    return n_ok, examples


def main():
    for p in (5, 7):
        search(p)


if __name__ == "__main__":
    main()
