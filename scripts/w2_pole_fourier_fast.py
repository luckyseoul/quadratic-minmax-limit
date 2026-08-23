#!/usr/bin/env python3
"""Fast factor-orbit test for the normalized switched pole family.

For w=c(D)gamma and a root zeta of an irreducible factor f of g,

    sum_j zeta^j w(gen^j) = c(zeta) sum_j zeta^j gamma(gen^j).

The second factor is nonzero because gamma is cyclic.  Thus f divides the
content c exactly when the multiplicative Fourier sum on the left vanishes.
This avoids constructing and inverting the full gamma-orbit coordinate map.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_gamma  # noqa: E402
from minmax_quadratic import halfspace_boolean_vector  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


def named_z_without_conference(p: int, tau: int = 0):
    """Construct 15.613's named z without allocating the conference matrix."""
    q, mul, add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    h = np.sign(halfspace_boolean_vector(p)).astype(np.int8)
    sig = next(e for e in range(1, q) if chi(e) == -1)
    sinv = _finv(mul, q, sig)
    neg_tau = ((p - tau % p) % p) + ((p - tau // p) % p) * p
    z = np.zeros(q + 1, dtype=np.int8)
    z[0] = np.int8(-h[0])
    for x in range(q):
        z[1 + x] = h[1 + add(mul(sinv, x), neg_tau)]
    return z, ((1 - z) // 2).astype(np.uint8), q, mul, add, chi, sinv // p


def mul_x_mod(value: int, modulus: int, degree: int) -> int:
    """Multiply a polynomial residue by X in F2[X]/modulus."""
    value <<= 1
    if value & (1 << degree):
        value ^= modulus
    return value & ((1 << degree) - 1)


def fourier_factor_residues(
    wfn: np.ndarray,
    orbit_points: list[int],
    factors: list[list[int]],
) -> list[int]:
    """Multiplicative Fourier residues, one in each factor field."""
    residues = []
    for i, fac in enumerate(factors):
        degree = len(fac) - 1
        modulus = sum(int(v) << j for j, v in enumerate(fac))
        power = 1
        acc = 0
        for point in orbit_points:
            if wfn[point]:
                acc ^= power
            power = mul_x_mod(power, modulus, degree)
        residues.append(acc)
    return residues


def fourier_factor_mask(
    wfn: np.ndarray,
    square_points: list[int],
    nonsquare_points: list[int],
    factors: list[list[int]],
    use_square: list[bool],
) -> int:
    """Return factors whose selected nonzero-generator projection vanishes."""
    square = fourier_factor_residues(wfn, square_points, factors)
    nonsquare = fourier_factor_residues(wfn, nonsquare_points, factors)
    mask = 0
    for i, choose_square in enumerate(use_square):
        if (square[i] if choose_square else nonsquare[i]) == 0:
            mask |= 1 << i
    return mask


def switched_wfn(p: int, t: int, z, bits, q, mul, add, chi):
    mat = (1, 0, t, p - 1)
    perm = _mobius_perm(p, *mat)
    if not np.array_equal(perm[perm], np.arange(q + 1)):
        raise AssertionError(f"not an involution at p={p}, t={t}")
    y = np.empty_like(z)
    for j in range(q + 1):
        src = int(perm[j])
        if j == 0:
            sw = chi(t)
        else:
            sw = chi(add(mul(t, j - 1), p - 1))
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    ybits = ((1 - y) // 2).astype(np.uint8)
    in_u = int(ybits[0]) == 1 and int(ybits[1]) == 0
    diff = (bits ^ ybits) & 1
    wfn = diff[1 : 1 + q].copy()
    if diff[0]:
        wfn ^= 1
    return in_u, wfn


def run(p: int, offsets: list[int], antipodal_pairs: bool = False) -> dict:
    # Keep factorization optional for the factorization-free scanners that
    # import the field/pole helpers from this module (notably on NUKA).
    from w2_ramanujan_mask_spectrum import aut_factor_orbits, factor_g_irreducible

    t0 = time.time()
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    g, factors = factor_g_irreducible(p)
    factor_orbits = aut_factor_orbits(p, g, factors)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    ncoord = (q - 1) // 2
    square_points = []
    point = 1
    for _ in range(ncoord):
        square_points.append(point)
        point = mul(gen, point)
    if point != 1 or len(set(square_points)) != ncoord:
        raise AssertionError(f"bad square cycle at p={p}")
    nonsquare_points = [mul(omega, point) for point in square_points]

    gamma, _q2, _mul2, _b = named_gamma(p)
    gamma_square = fourier_factor_residues(gamma, square_points, factors)
    gamma_nonsquare = fourier_factor_residues(gamma, nonsquare_points, factors)
    if any(a == 0 and b == 0 for a, b in zip(gamma_square, gamma_nonsquare)):
        raise AssertionError("gamma Fourier projection vanished on both point-orbits")
    use_square = [value != 0 for value in gamma_square]

    lo = (p + 1) // 2
    us = [lo + offset for offset in offsets]
    if len(set(us)) != len(us) or not all(lo <= u < p for u in us):
        raise ValueError("offsets must select distinct points in the upper interval")
    rows = []
    for u in us:
        t = lam * pow(u, p - 2, p) % p
        in_u, wfn = switched_wfn(p, t, z, bits, q, mul, add, chi)
        if not in_u:
            raise AssertionError(f"normalized upper point missed U at p={p}, u={u}")
        mask = fourier_factor_mask(
            wfn, square_points, nonsquare_points, factors, use_square
        )
        row = {"u": u, "t": t, "factor_mask": mask}
        if antipodal_pairs:
            opposite_u = p - u
            opposite_t = lam * pow(opposite_u, p - 2, p) % p
            opposite_in_u, opposite_wfn = switched_wfn(
                p, opposite_t, z, bits, q, mul, add, chi
            )
            if opposite_in_u:
                raise AssertionError(
                    f"both antipodal pole parameters entered U at p={p}, u={u}"
                )
            pair_mask = fourier_factor_mask(
                wfn ^ opposite_wfn,
                square_points,
                nonsquare_points,
                factors,
                use_square,
            )
            row.update(
                {
                    "opposite_u": opposite_u,
                    "opposite_t": opposite_t,
                    "antipodal_factor_mask": pair_mask,
                    "antipodal_unit_content": pair_mask == 0,
                }
            )
        rows.append(row)

    orbit_masks = [sum(1 << i for i in orbit) for orbit in factor_orbits]
    orbit_bad_counts = [
        sum((row["factor_mask"] & omask) == omask for row in rows)
        for omask in orbit_masks
    ]
    for row in rows:
        row["unit_content"] = row["factor_mask"] == 0
        row["valid_for_all_factor_orbits"] = all(
            (row["factor_mask"] & omask) != omask for omask in orbit_masks
        )
    result = {
        "p": p,
        "pole_lambda": lam,
        "u_offsets": offsets,
        "factor_degrees": [len(fac) - 1 for fac in factors],
        "factor_orbits": factor_orbits,
        "factor_orbit_degrees": [
            [len(factors[i]) - 1 for i in orbit] for orbit in factor_orbits
        ],
        "rows": rows,
        "orbit_product_bad_counts": orbit_bad_counts,
        "collective_W2_on_rows": all(count < len(rows) for count in orbit_bad_counts),
        "n_single_aut_orbit_witnesses": sum(
            row["valid_for_all_factor_orbits"] for row in rows
        ),
        "elapsed_seconds": time.time() - t0,
    }
    if antipodal_pairs:
        result["antipodal_pairs"] = {
            "n_pairs": len(rows),
            "n_unit_content": sum(row["antipodal_unit_content"] for row in rows),
            "all_unit_content": all(
                row["antipodal_unit_content"] for row in rows
            ),
        }
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--offsets", default="0,1")
    ap.add_argument("--antipodal-pairs", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.offsets == "all":
        offsets = list(range((args.p - 1) // 2))
    else:
        offsets = [int(value) for value in args.offsets.split(",")]
    result = run(args.p, offsets, antipodal_pairs=args.antipodal_pairs)
    print(
        f"p={args.p} rows={len(result['rows'])} factors={len(result['factor_degrees'])} "
        f"orbits={len(result['factor_orbits'])} "
        f"collective={result['collective_W2_on_rows']} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )
    output = args.output or Path(f"/tmp/w2_pole_fourier_fast_p{args.p}.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
