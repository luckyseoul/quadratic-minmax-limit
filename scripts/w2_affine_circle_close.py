#!/usr/bin/env python3
"""Exact certificates for the affine-halfspace close of the W2 circle route.

For a subset T of F_p with |T|=(p+1)/2, let h_T(infinity)=1 and

    h_T(a+b*omega) = +1 iff b is in T.

When ker(b)=F_p is a square direction, C h_T = p h_T for every such T.
After a nonsquare dilation sigma and the infinity sign switch,

    z_T(infinity)=-1,  z_T(sigma*u)=h_T(u)

is a Max-minus vector.  If 0 is in T, all z_T have the same restriction
to the standard nonsquare circle {infinity} union sigma*F_p.

For any two points sigma*u, sigma*v outside that circle, choose T so that
the membership signs of their transverse coordinates b,d have product
chi(u-v).  Then the pair is in the U slice.  Such a T always exists for
p>=5.  Flipping z_T on the circle supplies the circle incidence word as a
difference of two points in U.

This script checks the construction with exact integer arithmetic, exhausts
all affine subsets through p=11, and exhausts the elementary subset choice
for all outside pairs at the requested primes.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def sub(p: int, u: int, v: int) -> int:
    return (u % p - v % p) % p + ((u // p - v // p) % p) * p


def affine_halfspace(p: int, chosen: set[int]) -> np.ndarray:
    if len(chosen) != (p + 1) // 2:
        raise ValueError("chosen subset has the wrong size")
    h = np.ones(p * p + 1, dtype=np.int8)
    for u in range(p * p):
        h[1 + u] = 1 if u // p in chosen else -1
    return h


def choose_transverse_subset(
    p: int, b: int, d: int, target_product: int
) -> set[int]:
    """Choose T containing zero with s_T(b)s_T(d)=target_product."""
    if p < 5 or b == 0 or d == 0 or target_product not in (-1, 1):
        raise ValueError("outside-pair construction requires p>=5 and b,d!=0")
    k = (p + 1) // 2
    include = {0}
    exclude: set[int] = set()
    if b == d:
        if target_product != 1:
            raise ValueError("equal transverse coordinates force product +1")
        include.add(b)
    elif target_product == 1:
        include.update((b, d))
    else:
        include.add(b)
        exclude.add(d)
    for x in range(p):
        if len(include) == k:
            break
        if x not in include and x not in exclude:
            include.add(x)
    if len(include) != k:
        raise RuntimeError("failed to extend the prescribed memberships")

    def sign(x: int) -> int:
        return 1 if x in include else -1

    if 0 not in include or sign(b) * sign(d) != target_product:
        raise RuntimeError("constructed subset failed its sign prescription")
    return include


def standard_completion(p: int, chosen: set[int]):
    q, mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    sigma = next(x for x in range(1, q) if chi(x) == -1)
    sinv = pow_field(mul, q, sigma, q - 2)
    h = affine_halfspace(p, chosen)
    z = np.empty(q + 1, dtype=np.int8)
    z[0] = -1
    for x in range(q):
        z[1 + x] = h[1 + mul(sinv, x)]
    circle = np.zeros(q + 1, dtype=np.uint8)
    circle[0] = 1
    for t in range(p):
        circle[1 + mul(sigma, t)] = 1
    return h, z, circle, sigma


def pow_field(mul, q: int, x: int, exponent: int) -> int:
    out, base = 1, x
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        base = mul(base, base)
        exponent >>= 1
    return out


def exhaustive_affine_subsets(p: int, C: np.ndarray) -> dict:
    k = (p + 1) // 2
    tested = failures = 0
    for values in combinations(range(p), k):
        h = affine_halfspace(p, set(values))
        tested += 1
        if not np.array_equal(C @ h.astype(np.int64), p * h.astype(np.int64)):
            failures += 1
    return {"tested": tested, "failures": failures}


def exhaustive_outside_pairs(p: int) -> dict:
    q, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    outside = [u for u in range(q) if u // p != 0]
    tested = equal_transverse = failures = 0
    for ai, u in enumerate(outside):
        b = u // p
        for v in outside[ai + 1 :]:
            d = v // p
            target = chi(sub(p, u, v))
            tested += 1
            equal_transverse += b == d
            try:
                chosen = choose_transverse_subset(p, b, d, target)
            except (ValueError, RuntimeError):
                failures += 1
                continue
            sb = 1 if b in chosen else -1
            sd = 1 if d in chosen else -1
            if sb * sd != target:
                failures += 1
    return {
        "tested": tested,
        "equal_transverse": equal_transverse,
        "selection_failures": failures,
    }


def exact_frobenius_witness(p: int, C: np.ndarray) -> dict:
    q, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    sigma = next(x for x in range(1, q) if chi(x) == -1)
    sinv = pow_field(mul, q, sigma, q - 2)
    circle_vertices = {0} | {1 + mul(sigma, t) for t in range(p)}
    pair = None
    for i in range(1, q + 1):
        if i in circle_vertices:
            continue
        candidate_u = mul(sinv, i - 1)
        candidate_v = frob(candidate_u)
        j = 1 + mul(sigma, candidate_v)
        if j not in circle_vertices and i < j:
            pair = (i, j, candidate_u, candidate_v)
            break
    if pair is None:
        raise RuntimeError("failed to find a Frobenius-conjugate outside pair")
    i, j, u, v = pair
    b, d = u // p, v // p
    target = chi(sub(p, u, v))
    chosen = choose_transverse_subset(p, b, d, target)
    h, z, circle, sigma2 = standard_completion(p, chosen)
    if sigma2 != sigma:
        raise RuntimeError("standard nonsquare direction mismatch")
    zflip = z.copy()
    zflip[circle.astype(bool)] *= -1
    h_exact = np.array_equal(C @ h.astype(np.int64), p * h.astype(np.int64))
    z_exact = np.array_equal(C @ z.astype(np.int64), -p * z.astype(np.int64))
    flip_exact = np.array_equal(
        C @ zflip.astype(np.int64), -p * zflip.astype(np.int64)
    )
    pair_u = int(C[i, j]) * int(z[i]) * int(z[j]) == -1
    pair_u_flip = int(C[i, j]) * int(zflip[i]) * int(zflip[j]) == -1
    diff = ((1 - z) // 2) ^ ((1 - zflip) // 2)
    support_exact = np.array_equal(diff.astype(np.uint8), circle)
    return {
        "frobenius_pair": [int(i), int(j)],
        "transverse_coordinates": [int(b), int(d)],
        "target_product": int(target),
        "chosen_subset": sorted(chosen),
        "affine_plus_exact": bool(h_exact),
        "completion_minus_exact": bool(z_exact),
        "circle_flip_minus_exact": bool(flip_exact),
        "pair_in_U_before_and_after": bool(pair_u and pair_u_flip),
        "difference_support_is_circle": bool(support_exact),
        "witness_bits_hex": np.packbits(
            ((1 - z) // 2).astype(np.uint8), bitorder="little"
        ).tobytes().hex(),
    }


def direct_p3_slice_rank() -> dict:
    """The p=3 exception is finite; W2 is vacuous and the slice has rank 4."""
    from e1_gmin_m4_prop15406 import gf2_rref, load_minus

    Y, C = load_minus(3)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    in_u = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1] == -1
    bits = ((1 - Y[in_u]) // 2).astype(np.uint8)
    directions = bits ^ bits[0]
    rank = int(gf2_rref(directions.copy())[2])
    return {
        "U_points": int(len(bits)),
        "direction_rank": rank,
        "target_rank": 4,
        "full_slice": rank == 4,
        "W2_vacuous": True,
    }


def run(primes: list[int]) -> dict:
    t0 = time.time()
    rows = {}
    all_ok = True
    for p in primes:
        C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
        row = {
            "all_affine_subsets": (
                exhaustive_affine_subsets(p, C) if p <= 11 else None
            ),
            "all_outside_pairs": exhaustive_outside_pairs(p),
            "frobenius_witness": exact_frobenius_witness(p, C),
        }
        subset_ok = row["all_affine_subsets"] is None or (
            row["all_affine_subsets"]["failures"] == 0
        )
        pair_ok = row["all_outside_pairs"]["selection_failures"] == 0
        witness_ok = all(
            row["frobenius_witness"][key]
            for key in (
                "affine_plus_exact",
                "completion_minus_exact",
                "circle_flip_minus_exact",
                "pair_in_U_before_and_after",
                "difference_support_is_circle",
            )
        )
        row["all_checks"] = bool(subset_ok and pair_ok and witness_ok)
        all_ok = all_ok and row["all_checks"]
        rows[str(p)] = row
        print(json.dumps({"p": p, **row}, indent=2), flush=True)
    return {
        "title": "W2 affine-halfspace circle completion close",
        "proved_for_p_ge_5": True,
        "p3_direct": direct_p3_slice_rank(),
        "all_checks": bool(all_ok),
        "rows": rows,
        "seconds": round(time.time() - t0, 3),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", type=int, nargs="+", default=[5, 7, 11, 19])
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = run(args.primes)
    print(json.dumps(out, indent=2), flush=True)
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
