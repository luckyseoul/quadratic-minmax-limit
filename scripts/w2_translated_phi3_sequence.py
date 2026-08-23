#!/usr/bin/env python3
"""Full Phi_3 translation sequence from four affine-line pullbacks.

For the translated antipodal boundary pair w_a, consecutive halfspaces
differ only on the two affine lines

    ell(x) = (p-1)/2-a,  -1-a.

The switched-pole signs are independent of a, so w_{a+1}+w_a is exactly
the xor of the pullbacks of those lines under the two antipodal poles.
At Phi_3 one square/nonsquare Fourier component of the cyclic generator is
zero and the other is nonzero.  One direct w_0 evaluation plus line-residue
tables therefore recovers every content c_a mod Phi_3 in O(p^2) time.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_gamma  # noqa: E402
from w2_pole_fourier_fast import named_z_without_conference  # noqa: E402
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    apply_pole,
    is_prime,
    pole_action_data,
)
from w2_translated_phi3_boundary_scan import (  # noqa: E402
    gf4_mul,
    gf4_star,
    phi3_residue,
)


GF4_POWERS = np.asarray([1, 2, 3], dtype=np.uint8)


def record(p: int) -> dict:
    t0 = time.time()
    _z, _bits, q, mul, _add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, ib = field_ctx(p)
    if q2 != q:
        raise AssertionError("field contexts disagree")

    sig = next(value for value in range(1, q) if chi(value) == -1)
    sinv = _finv(mul, q, sig)
    field_points = np.arange(q, dtype=np.int64)
    point_a = field_points % p
    point_b = field_points // p
    sinv_a, sinv_b = sinv % p, sinv // p
    ell = (
        sinv_a * point_b
        + sinv_b * point_a
        + sinv_b * point_b * ia
    ) % p

    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = (
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
        )

    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    ncoord = (q - 1) // 2
    square = np.empty(ncoord, dtype=np.int64)
    point = 1
    for j in range(ncoord):
        square[j] = point
        point = mul(gen, point)
    nonsquare = np.fromiter(
        (mul(omega, point) for point in square),
        dtype=np.int64,
        count=ncoord,
    )

    half = (p - 1) // 2
    gamma, qg, _mulg, _b = named_gamma(p)
    if qg != q:
        raise AssertionError("gamma field disagrees")
    gamma_residues = [
        phi3_residue(gamma, square),
        phi3_residue(gamma, nonsquare),
    ]
    active = [i for i, value in enumerate(gamma_residues) if value]
    if len(active) != 1:
        raise AssertionError(
            f"Phi_3 generator should have one live component: "
            f"p={p}, residues={gamma_residues}"
        )
    active_component = active[0]
    gamma_inverse = gf4_star(gamma_residues[active_component])

    # Residue of each pullback 1_{ell(pi_t(x))=r}, split by the two
    # multiplicative point-orbits.  The possible infinity correction is an
    # all-ones vector, whose Phi_3 residue is zero because 3|ncoord.
    line_residues = np.zeros((2, 2, p), dtype=np.uint8)
    shifted = ell
    z0 = np.empty(q + 1, dtype=np.int8)
    z0[0] = -1
    z0[1:] = np.where(shifted <= half, 1, -1).astype(np.int8)
    bits0 = ((1 - z0) // 2).astype(np.uint8)
    pair0 = None
    for action_index, u in enumerate((half, p - half)):
        pole_t = lam * pow(u, p - 2, p) % p
        perm, switch = pole_action_data(
            p, pole_t, ia, ib, inv_fp, legendre
        )
        _in_u, wfn = apply_pole(z0, bits0, perm, switch)
        pair0 = wfn.copy() if pair0 is None else pair0 ^ wfn
        for component, orbit in enumerate((square, nonsquare)):
            sources = perm[1 + orbit]
            finite = sources != 0
            levels = ell[sources[finite] - 1]
            exponent_classes = np.arange(ncoord, dtype=np.int64)[finite] % 3
            for exponent_class, value in enumerate(GF4_POWERS):
                parity = np.bincount(
                    levels[exponent_classes == exponent_class], minlength=p
                ) & 1
                line_residues[action_index, component] ^= (
                    parity.astype(np.uint8) * value
                )
        del perm, switch, wfn

    initial_residues = [
        phi3_residue(pair0, square),
        phi3_residue(pair0, nonsquare),
    ]
    residues = initial_residues.copy()
    content_values = []
    inactive_component_nonzero = 0
    for a in range(p):
        if residues[1 - active_component] != 0:
            inactive_component_nonzero += 1
        content_values.append(
            gf4_mul(residues[active_component], gamma_inverse)
        )
        r_enter = (half - a) % p
        r_leave = (p - 1 - a) % p
        for component in range(2):
            delta = 0
            for action_index in range(2):
                delta ^= int(
                    line_residues[action_index, component, r_enter]
                )
                delta ^= int(
                    line_residues[action_index, component, r_leave]
                )
            residues[component] ^= delta

    if residues != initial_residues:
        raise AssertionError(f"translation recurrence did not close at p={p}")
    if inactive_component_nonzero:
        raise AssertionError(
            f"inactive Phi_3 component nonzero {inactive_component_nonzero} "
            f"times at p={p}"
        )
    reflection_mismatches = sum(
        content_values[(half - a) % p] != content_values[a]
        for a in range(p)
    )
    if reflection_mismatches:
        raise AssertionError(
            f"c_(m-a)=c_a failed {reflection_mismatches} times at p={p}"
        )

    valid = content_values[1:half]
    representatives = content_values[1 : half // 2 + 1]
    zero_positions = [a for a, value in enumerate(content_values) if value == 0]
    valid_zero_positions = [a for a in range(1, half) if content_values[a] == 0]
    return {
        "p": p,
        "boundary_u": half,
        "gamma_residues": gamma_residues,
        "active_component": "square" if active_component == 0 else "nonsquare",
        "content_value_counts": {
            str(value): content_values.count(value) for value in range(4)
        },
        "n_zeros_all_translations": len(zero_positions),
        "n_zeros_valid_translations": len(valid_zero_positions),
        "n_valid_translations": len(valid),
        "n_good_valid_translations": sum(bool(value) for value in valid),
        "all_valid_bad": bool(valid) and not any(valid),
        "reflection_c_m_minus_a_equals_c_a": True,
        "mod12_5_total_zero_formula": (
            len(zero_positions) == (2 * p - 7) // 3
            if p % 12 == 5
            else None
        ),
        "first_good_representative_a": next(
            (a for a, value in enumerate(representatives, 1) if value),
            None,
        ),
        "zero_positions": zero_positions,
        "content_values": content_values,
        "elapsed_seconds": time.time() - t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument(
        "--omit-sequences",
        action="store_true",
        help="drop zero_positions/content_values from the JSON output",
    )
    args = ap.parse_args()
    primes = [p for p in range(max(5, args.start), args.stop + 1) if is_prime(p)]
    t0 = time.time()
    if args.workers == 1:
        rows = list(map(record, primes))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(record, primes, chunksize=1))
    if args.omit_sequences:
        for row in rows:
            row.pop("zero_positions")
            row.pop("content_values")
    out = {
        "range": [args.start, args.stop],
        "n_primes": len(rows),
        "all_valid_bad": [row["p"] for row in rows if row["all_valid_bad"]],
        "mod12_5_total_zero_formula_failures": [
            row["p"]
            for row in rows
            if row["mod12_5_total_zero_formula"] is False
        ],
        "max_zero_fraction_all_translations": max(
            (row["n_zeros_all_translations"] / row["p"] for row in rows),
            default=None,
        ),
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"p={args.start}..{args.stop} primes={len(rows)} "
        f"all-valid-bad={out['all_valid_bad']} "
        f"max-zero-fraction={out['max_zero_fraction_all_translations']:.6f} "
        f"seconds={out['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
