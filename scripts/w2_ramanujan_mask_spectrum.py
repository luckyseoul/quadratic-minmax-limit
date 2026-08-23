#!/usr/bin/env python3
"""Joint coprimality spectrum for the switched split-involution W2 class.

This is a lab script, not a numbered proposition.  It keeps the complete
gcd mask against the irreducible factors of

    g = (X^m + 1)/(X + 1),  m = oddpart((p^2 - 1)/2),

instead of only the per-factor marginals.  The mask histogram gives the
exact polynomial-Moebius/Ramanujan expansion of the unit-content count and
therefore retains all correlations between bad factors.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15612 import _f2_divmod  # noqa: E402
from e1_gmin_m4_prop15613 import _finv, named_gamma, named_z  # noqa: E402
from e1_gmin_m4_prop15616 import _g_factors  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from walsh_linecode_rank import _mobius_perm  # noqa: E402


CTX: dict[str, object] = {}


def split_involutions(p: int) -> list[tuple[int, int, int, int]]:
    """One representative per projective pair +/-M."""
    out: list[tuple[int, int, int, int]] = []
    seen: set[tuple[int, int, int, int]] = set()
    for a in range(p):
        for b in range(p):
            need = (1 - a * a) % p
            if b == 0:
                if need:
                    continue
                cs = range(p)
            else:
                cs = [(need * pow(b, p - 2, p)) % p]
            for c in cs:
                raw = (a, b, c, (-a) % p)
                neg = tuple((-x) % p for x in raw)
                key = min(raw, neg)
                if key not in seen:
                    seen.add(key)
                    out.append(key)
    return out


def factor_g_irreducible(p: int) -> tuple[list[int], list[list[int]]]:
    """Factor g completely over F2; coefficients are low degree first."""
    g, _blocks = _g_factors(p)
    x = sp.symbols("x")
    poly = sp.Poly(sum(int(c) * x**i for i, c in enumerate(g)), x, modulus=2)
    _unit, raw = sp.factor_list(poly, modulus=2)
    factors: list[list[int]] = []
    for fac, exponent in raw:
        if exponent != 1:
            raise AssertionError(f"g is not square-free at p={p}: {fac}^{exponent}")
        coeff = [int(v) & 1 for v in reversed(fac.all_coeffs())]
        factors.append(coeff)
    factors.sort(key=lambda f: (len(f), f))
    prod = [1]
    for fac in factors:
        nxt = [0] * (len(prod) + len(fac) - 1)
        for i, a in enumerate(prod):
            if a:
                for j, b in enumerate(fac):
                    if b:
                        nxt[i + j] ^= 1
        prod = nxt
    if prod != list(map(int, g)):
        raise AssertionError(f"factorization did not reconstruct g at p={p}")
    return list(map(int, g)), factors


def substitute_power_mod(poly: list[int], power: int, modulus: list[int]) -> list[int]:
    out = [0] * ((len(poly) - 1) * power + 1)
    for i, coeff in enumerate(poly):
        if coeff:
            out[i * power] = 1
    return _f2_divmod(out, modulus)[1]


def aut_factor_orbits(p: int, g: list[int], factors: list[list[int]]) -> list[list[int]]:
    """Orbits under X->X^p and X->X^{-1} in F2[X]/g."""
    m = len(g)
    neighbours: dict[int, set[int]] = {i: set() for i in range(len(factors))}
    for i, fac in enumerate(factors):
        for power in (p, m - 1):
            transformed = substitute_power_mod(fac, power, g)
            hits = [j for j, other in enumerate(factors) if _poly_gcd(transformed, other) != [1]]
            if len(hits) != 1:
                raise AssertionError(
                    f"factor transform at p={p}, i={i}, power={power} hits {hits}"
                )
            j = hits[0]
            neighbours[i].add(j)
            neighbours[j].add(i)
    orbits: list[list[int]] = []
    unseen = set(range(len(factors)))
    while unseen:
        root = min(unseen)
        todo = [root]
        orbit = set()
        while todo:
            i = todo.pop()
            if i in orbit:
                continue
            orbit.add(i)
            todo.extend(neighbours[i] - orbit)
        unseen -= orbit
        orbits.append(sorted(orbit))
    return orbits


def coordinate_solver(p: int):
    """Precompute c from w=sum c_k D^k gamma using independent rows."""
    gamma, q, mul, _b = named_gamma(p)
    ncoord = (q - 1) // 2
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    cols = [gamma.copy()]
    cur = gamma.copy()
    for _ in range(ncoord - 1):
        cur = _dil_fn(cur, mul, gen, q)
        cols.append(cur.copy())
    matrix = np.stack(cols, axis=1).astype(np.uint8)
    _r, pivot_rows, rank = gf2_rref(matrix.T)
    if rank != ncoord:
        raise AssertionError(f"gamma orbit rank {rank} != {ncoord} at p={p}")
    pivot_rows = pivot_rows[:ncoord]
    square = matrix[pivot_rows, :]
    aug = np.concatenate([square, np.eye(ncoord, dtype=np.uint8)], axis=1)
    reduced, pivots, inv_rank = gf2_rref(aug)
    if inv_rank != ncoord or pivots[:ncoord] != list(range(ncoord)):
        raise AssertionError(f"coordinate minor singular at p={p}")
    inverse = reduced[:, ncoord:].copy()
    if not np.array_equal((inverse.astype(np.int32) @ square.astype(np.int32)) & 1,
                          np.eye(ncoord, dtype=np.int32)):
        raise AssertionError(f"bad coordinate inverse at p={p}")
    return np.asarray(pivot_rows, dtype=np.int32), inverse, matrix, mul, gen


def setup(p: int) -> None:
    z, bits, eigen, in_u, q, mul, add, chi, sig = named_z(p)
    if not eigen or not in_u:
        raise AssertionError(f"named z failed at p={p}")
    g, factors = factor_g_irreducible(p)
    factor_orbits = aut_factor_orbits(p, g, factors)
    pivot_rows, inverse, matrix, mul2, gen = coordinate_solver(p)
    if mul2 is not mul:
        # field_ctx returns fresh closures; equality is not semantically useful.
        pass
    CTX.clear()
    CTX.update(
        p=p,
        q=q,
        z=z,
        bits=bits,
        mul=mul,
        add=add,
        chi=chi,
        g=g,
        factors=factors,
        factor_orbits=factor_orbits,
        pole_lambda=_finv(mul, q, sig) // p,
        pivot_rows=pivot_rows,
        inverse=inverse,
        basis_matrix=matrix,
        gen=gen,
    )


def content_row(mat: tuple[int, int, int, int], require_in_u: bool = True):
    """Return content data, optionally retaining rows outside the U slice."""
    p = int(CTX["p"])
    q = int(CTX["q"])
    z = CTX["z"]
    bits = CTX["bits"]
    mul = CTX["mul"]
    add = CTX["add"]
    chi = CTX["chi"]
    pivot_rows = CTX["pivot_rows"]
    inverse = CTX["inverse"]
    basis_matrix = CTX["basis_matrix"]
    factors = CTX["factors"]
    factor_orbits = CTX["factor_orbits"]
    a, b, c, d0 = mat
    perm = _mobius_perm(p, a, b, c, d0)
    if not np.array_equal(perm[perm], np.arange(q + 1)):
        raise AssertionError(f"split matrix is not an involution: {mat}")
    y = np.empty_like(z)
    for j in range(q + 1):
        src = int(perm[j])
        if j == 0:
            sw = chi(c) if c else 1
        else:
            lin = add(mul(c, j - 1), d0)
            sw = chi(lin)
            if sw == 0:
                sw = 1
        y[j] = np.int8(int(sw) * int(z[src]))
    ybits = ((1 - y) // 2).astype(np.uint8)
    if require_in_u and not (int(ybits[0]) == 1 and int(ybits[1]) == 0):
        return None
    diff = (bits ^ ybits) & 1
    wfn = diff[1 : 1 + q].copy()
    if diff[0]:
        wfn ^= 1
    cvec = (inverse @ wfn[pivot_rows]) & 1
    recon = (basis_matrix @ cvec) & 1
    if not np.array_equal(recon, wfn):
        raise AssertionError(f"coordinate reconstruction failed for {mat}")
    coeff = list(map(int, cvec))
    content_int = sum(int(v) << i for i, v in enumerate(coeff))
    mask = 0
    for i, fac in enumerate(factors):
        if _poly_gcd(coeff, fac) != [1]:
            mask |= 1 << i
    gcd_poly = _poly_gcd(coeff, CTX["g"])
    full_data = None
    if mask == (1 << len(factors)) - 1:
        quotient, remainder = _f2_divmod(coeff, CTX["g"])
        if remainder != [0]:
            raise AssertionError(f"full mask without g divisibility for {mat}")
        quotient_int = sum(int(v) << i for i, v in enumerate(quotient))
        full_data = {
            "quotient_degree": len(quotient) - 1,
            "quotient_weight": sum(map(int, quotient)),
            "quotient_at_1": sum(map(int, quotient)) & 1,
            "quotient_hex": hex(quotient_int),
        }
    return mat, mask, len(gcd_poly) - 1, int(cvec.sum()), full_data, content_int


def class_row(mat: tuple[int, int, int, int]):
    """ProcessPool entry point: retain only the fixed U pair-slice."""
    return content_row(mat, require_in_u=True)


def summarize(p: int, rows, elapsed: float, n_candidates: int, family: str) -> dict:
    factors = CTX["factors"]
    factor_orbits = CTX["factor_orbits"]
    hist = Counter(mask for _mat, mask, _gd, _wt, _full, _ci in rows)
    gcd_degree_hist = Counter(gd for _mat, _mask, gd, _wt, _full, _ci in rows)
    bad_count_hist = Counter(mask.bit_count() for mask in hist.elements())
    n = len(rows)
    r = len(factors)
    ie_by_order = []
    running = 0
    for k in range(r + 1):
        unsigned = sum(count * math.comb(mask.bit_count(), k) for mask, count in hist.items())
        signed = unsigned if k % 2 == 0 else -unsigned
        running += signed
        ie_by_order.append(
            {
                "order": k,
                "unsigned_intersection_sum": unsigned,
                "signed_contribution": signed,
                "partial_sum": running,
            }
        )
    n_unit = hist.get(0, 0)
    if running != n_unit:
        raise AssertionError(f"Moebius total {running} != unit count {n_unit}")
    model = math.prod(1.0 - 2.0 ** (-(len(f) - 1)) for f in factors)
    orbit_masks = [sum(1 << i for i in orbit) for orbit in factor_orbits]
    orbit_bad_counts = [
        sum(count for mask, count in hist.items() if (mask & omask) == omask)
        for omask in orbit_masks
    ]
    orbit_bad_matrices = [
        [list(mat) for mat, mask, _gd, _wt, _full, _ci in rows
         if (mask & omask) == omask]
        for omask in orbit_masks
    ]
    aut_unit_rows = [
        mat
        for mat, mask, _gd, _wt, _full, _ci in rows
        if all((mask & omask) != omask for omask in orbit_masks)
    ]
    xor_content_int = 0
    for _mat, _mask, _gd, _wt, _full, content_int in rows:
        xor_content_int ^= content_int
    ncoord = int(CTX["inverse"].shape[0])
    xor_coeff = [(xor_content_int >> i) & 1 for i in range(ncoord)]
    xor_gcd = _poly_gcd(xor_coeff, CTX["g"])
    xor_factor_mask = 0
    for i, fac in enumerate(factors):
        if _poly_gcd(xor_coeff, fac) != [1]:
            xor_factor_mask |= 1 << i
    xor_bad_orbits = [
        i for i, omask in enumerate(orbit_masks)
        if (xor_factor_mask & omask) == omask
    ]
    mask_examples: dict[str, list[list[int]]] = {}
    full_g_examples = []
    stratum_rows: dict[str, list[tuple]] = {
        "affine_c0": [],
        "pole_b0_c_nonzero": [],
        "generic_bc_nonzero": [],
    }
    for mat, mask, _gd, content_weight, full_data, _content_int in rows:
        key = str(mask)
        if len(mask_examples.setdefault(key, [])) < 8:
            mask_examples[key].append(list(mat))
        if full_data is not None:
            full_g_examples.append(
                {"matrix": list(mat), "content_weight": content_weight, **full_data}
            )
        if mat[2] == 0:
            stratum_rows["affine_c0"].append((mat, mask))
        elif mat[1] == 0:
            stratum_rows["pole_b0_c_nonzero"].append((mat, mask))
        else:
            stratum_rows["generic_bc_nonzero"].append((mat, mask))
    strata = {}
    for name, subset in stratum_rows.items():
        sh = Counter(mask.bit_count() for _mat, mask in subset)
        su = sum(mask == 0 for _mat, mask in subset)
        strata[name] = {
            "n_inU": len(subset),
            "n_unit": su,
            "unit_rate": su / len(subset) if subset else None,
            "bad_factor_count_histogram": {str(k): v for k, v in sorted(sh.items())},
        }
    return {
        "p": p,
        "family": family,
        "n_candidates": n_candidates,
        "n_inU": n,
        "n_unit": n_unit,
        "unit_rate": n_unit / n,
        "random_unit_model": model,
        "model_ratio": (n_unit / n) / model,
        "factor_degrees": [len(f) - 1 for f in factors],
        "factor_orbits": factor_orbits,
        "factor_orbit_degrees": [
            [len(factors[i]) - 1 for i in orbit] for orbit in factor_orbits
        ],
        "orbit_product_bad_counts": orbit_bad_counts,
        "orbit_product_bad_matrices": orbit_bad_matrices if family != "class" else None,
        "collective_W2_on_family": all(v < n for v in orbit_bad_counts),
        "n_single_aut_orbit_witnesses": len(aut_unit_rows),
        "single_aut_orbit_witness_matrices": [list(mat) for mat in aut_unit_rows[:100]],
        "xor_all_family_content": {
            "weight": xor_content_int.bit_count(),
            "gcd_degree_with_g": len(xor_gcd) - 1,
            "coprime_to_g": xor_gcd == [1],
            "factor_mask": xor_factor_mask,
            "bad_factor_orbits": xor_bad_orbits,
            "valid_for_all_factor_orbits": not xor_bad_orbits,
            "content_hex": hex(xor_content_int),
        },
        "pole_lambda": int(CTX["pole_lambda"]) if family == "pole" else None,
        "n_irreducible_factors": r,
        "mask_histogram": {str(mask): count for mask, count in sorted(hist.items())},
        "mask_examples": mask_examples,
        "unit_matrices": [
            list(mat) for mat, mask, _gd, _wt, _full, _ci in rows if mask == 0
        ],
        "compact_rows": [
            {
                "matrix": list(mat),
                "u_parameter": (
                    int(CTX["pole_lambda"])
                    * pow(mat[2], p - 2, p)
                    % p
                    if family == "pole" else None
                ),
                "mask": mask,
                "gcd_degree": gd,
                "content_weight": wt,
                "content_hex": hex(content_int),
            }
            for mat, mask, gd, wt, _full, content_int in rows
        ] if family != "class" else None,
        "full_g_examples": full_g_examples,
        "matrix_strata": strata,
        "bad_factor_count_histogram": {str(k): v for k, v in sorted(bad_count_hist.items())},
        "gcd_degree_histogram": {str(k): v for k, v in sorted(gcd_degree_hist.items())},
        "inclusion_exclusion_by_order": ie_by_order,
        "elapsed_seconds": elapsed,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--family", choices=("class", "pole", "affine"), default="class")
    ap.add_argument(
        "--pole-u-offsets",
        type=str,
        help="comma-separated offsets from u=(p+1)/2; pole family only",
    )
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    t0 = time.time()
    print(f"setup p={args.p}", flush=True)
    setup(args.p)
    if args.family == "class":
        mats = split_involutions(args.p)
    elif args.family == "pole":
        if args.pole_u_offsets:
            offsets = [int(v) for v in args.pole_u_offsets.split(",")]
            lo = (args.p + 1) // 2
            us = [lo + offset for offset in offsets]
            if len(set(us)) != len(us) or not all(lo <= u < args.p for u in us):
                ap.error("pole u offsets must select distinct points in the upper interval")
            lam = int(CTX["pole_lambda"])
            mats = [
                (1, 0, lam * pow(u, args.p - 2, args.p) % args.p, args.p - 1)
                for u in us
            ]
        else:
            mats = [(1, 0, c, args.p - 1) for c in range(1, args.p)]
    else:
        mats = [(1, b, 0, args.p - 1) for b in range(args.p)]
    print(
        f"p={args.p} class={len(mats)} factors="
        f"{[len(f)-1 for f in CTX['factors']]} workers={args.workers}",
        flush=True,
    )
    rows = []
    if args.workers == 1:
        iterator = map(class_row, mats)
        for row in iterator:
            if row is not None:
                rows.append(row)
    else:
        # Linux fork workers inherit the read-only precomputed coordinate map.
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            for i, row in enumerate(ex.map(class_row, mats, chunksize=4), 1):
                if row is not None:
                    rows.append(row)
                if i % 200 == 0:
                    print(f"  {i}/{len(mats)} inU={len(rows)}", flush=True)
    result = summarize(args.p, rows, time.time() - t0, len(mats), args.family)
    print(
        f"p={args.p} inU={result['n_inU']} unit={result['n_unit']} "
        f"rate={result['unit_rate']:.4f} bad-K={result['bad_factor_count_histogram']}",
        flush=True,
    )
    output = args.output or Path(f"/tmp/w2_ramanujan_mask_p{args.p}.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
