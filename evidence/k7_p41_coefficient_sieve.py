#!/usr/bin/env python3
"""p=41 k=7 coupled sieve, Boolean reconstruction, and |Z_ψ|² census.

The unique energy partition is 7×b_min=T=210.  After translation kills
the degree-4 kernel, types are depressed quintics.  Singer orbits of
every 7-subset are sieved by the cubic/quadratic/linear kernels; survivors
are Boolean-lifted (k=6 endpoint test) and scored by:

  Z_ψ = ∑_{d≠0} ψ(d) N(d),   N(d)=|D ∩ (D−d)|,

which is the same bilinear form as the k=6 quartic kernel (ψ = i^{ord}).
a_L = ∑_s h_L(s)² is checked in 2pℤ and against 2p·30.

Backend: ProcessPool over leadings and over Singer orbits.  GPU unused
(F_p 4×4 + lookups, not a dense batch).  Serial only inside one orbit.

This module is the proving/refuting unit imported by
`qvar_k_ge_7_proved_general`: the (41,7) stratum clears QVAR iff the
computed E|Z_ψ|² meets 3q(q−1)/16.
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations, product
from math import prod
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k5_p29_coefficient_sieve import (  # noqa: E402
    cyclic_direction_permutation,
    field_ctx,
    homogeneous_matrix,
    kernel_modp,
    square_directions,
    subset_orbits,
)

P = 41
K = 7
GMIN = 30
REPO = HERE.parent


def inverse_modp(matrix: np.ndarray, p: int) -> np.ndarray:
    n = matrix.shape[0]
    work = np.concatenate(
        [np.asarray(matrix, dtype=np.int64) % p, np.eye(n, dtype=np.int64)],
        axis=1,
    )
    for column in range(n):
        pivot = next(row for row in range(column, n) if work[row, column] % p)
        work[[column, pivot]] = work[[pivot, column]]
        work[column] *= pow(int(work[column, column]), p - 2, p)
        work[column] %= p
        for row in range(n):
            if row != column and work[row, column]:
                work[row] -= work[row, column] * work[column]
                work[row] %= p
    return work[:, n:] % p


def _n_workers(default: int) -> int:
    env = os.environ.get("K7_P41_WORKERS")
    if env:
        return max(1, int(env))
    n = os.cpu_count() or 4
    cap = min(default, max(1, n - 2))
    if os.environ.get("PYTEST_CURRENT_TEST"):
        return min(cap, 40)
    return cap


def types_for_leading(args: tuple[int, int, int]) -> tuple[int, list]:
    p, leading, target = args
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    centered = np.where(s <= midpoint, s, s - p)
    shifted = np.stack([centered[(s + constant) % p] for constant in range(p)])
    shifted_sq = shifted * shifted
    cubic = np.repeat(np.arange(p, dtype=np.int64), p * p)
    quadratic = np.tile(np.repeat(np.arange(p, dtype=np.int64), p), p)
    linear = np.tile(np.arange(p, dtype=np.int64), p * p)
    row_count = p**3
    rows = np.repeat(np.arange(row_count), p)
    s2, s3, s5 = s * s, s**3, s**5
    values = (
        leading * s5[None, :]
        + cubic[:, None] * s3[None, :]
        + quadratic[:, None] * s2[None, :]
        + linear[:, None] * s[None, :]
    ) % p
    counts = np.zeros((row_count, p), dtype=np.int16)
    np.add.at(counts, (rows, values.ravel()), 1)
    standard_sum = counts @ shifted.T
    standard_energy = counts @ shifted_sq.T
    endpoint_counts = np.stack(
        [counts[:, (midpoint - constant) % p] for constant in range(p)],
        axis=1,
    )
    replacements = standard_sum // p
    valid = (
        (standard_sum % p == 0)
        & (replacements >= 0)
        & (replacements <= endpoint_counts)
    )
    energy = (standard_energy + p * replacements) // (2 * p)
    hit = valid & (energy == target)
    recs = []
    if np.any(hit):
        ri, constants = np.nonzero(hit)
        recs = [
            (int(cubic[row]), int(quadratic[row]), int(linear[row]), int(constant))
            for row, constant in zip(ri, constants)
        ]
    return leading, recs


def enumerate_types(p: int, target: int) -> dict[int, list]:
    workers = _n_workers(p - 1)
    jobs = [(p, leading, target) for leading in range(1, p)]
    out: dict[int, list] = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for leading, recs in pool.map(types_for_leading, jobs, chunksize=1):
            if recs:
                out[int(leading)] = recs
    return out


def index_types(types_by_leading: dict[int, list]) -> dict[int, dict]:
    out = {}
    for leading, recs in types_by_leading.items():
        by_c = defaultdict(list)
        for rec in recs:
            by_c[rec[0]].append(rec)
        out[leading] = dict(by_c)
    return out


def sieve_one_orbit(args):
    subset, orbit_size, forms, types_by_leading, types_by_cubic, p = args
    selected = [forms[i] for i in subset]
    top = kernel_modp(homogeneous_matrix(selected, 5, p), p)
    if len(top) != 1 or np.any(top[0] == 0):
        return (tuple(subset), int(orbit_size), 0, [])
    kappa = top[0]
    cubic_m = homogeneous_matrix(selected, 3, p)
    quad_m = homogeneous_matrix(selected, 2, p)
    lin_m = homogeneous_matrix(selected, 1, p)
    hits = []
    n_coeff = 0
    for scalar in range(1, p):
        leading = [(scalar * int(kappa[j])) % p for j in range(7)]
        groups = [types_by_leading.get(leading[j], []) for j in range(7)]
        if any(len(g) == 0 for g in groups):
            continue
        order = sorted(range(7), key=lambda j: len(groups[j]))
        a0, a1, a2 = order[0], order[1], order[2]
        b = order[3:]
        try:
            inv = inverse_modp(cubic_m[:, b], p)
        except StopIteration:
            return (tuple(subset), int(orbit_size), 0, [])
        g0, g1, g2 = groups[a0], groups[a1], groups[a2]
        col_a0, col_a1, col_a2 = cubic_m[:, a0], cubic_m[:, a1], cubic_m[:, a2]
        for t0 in g0:
            for t1 in g1:
                for t2 in g2:
                    rhs = -(col_a0 * t0[0] + col_a1 * t1[0] + col_a2 * t2[0]) % p
                    c_b = (inv @ rhs) % p
                    pools = []
                    ok = True
                    for idx, bj in enumerate(b):
                        pool = types_by_cubic.get(leading[bj], {}).get(
                            int(c_b[idx]), []
                        )
                        if not pool:
                            ok = False
                            break
                        pools.append(pool)
                    if not ok:
                        continue
                    for u0, u1, u2, u3 in product(*pools):
                        rec = [None] * 7
                        rec[a0], rec[a1], rec[a2] = t0, t1, t2
                        rec[b[0]], rec[b[1]], rec[b[2]], rec[b[3]] = u0, u1, u2, u3
                        cvec = np.array([r[0] for r in rec], dtype=np.int64)
                        dvec = np.array([r[1] for r in rec], dtype=np.int64)
                        evec = np.array([r[2] for r in rec], dtype=np.int64)
                        if np.any((cvec @ cubic_m.T) % p):
                            continue
                        if np.any((dvec @ quad_m.T) % p):
                            continue
                        if np.any((evec @ lin_m.T) % p):
                            continue
                        n_coeff += 1
                        hits.append((leading, rec))
    return (tuple(subset), int(orbit_size), n_coeff, hits)


def official_zpsi_from_D(negative: np.ndarray, p: int) -> complex:
    """Z_ψ = ∑_{d≠0} ψ(d) N(d), N(d)=|D ∩ (D−d)|, F_q in the (1,p) basis."""
    q, multiply, _character, _trace = field_ctx(p)
    D = np.flatnonzero(negative)
    n = np.zeros(q, dtype=np.int64)
    for x in D:
        xi, xj = int(x) % p, int(x) // p
        for y in D:
            yi, yj = int(y) % p, int(y) // p
            d = (xi - yi) % p + ((xj - yj) % p) * p
            if d:
                n[d] += 1
    real, imag = np.zeros(q, dtype=np.int8), np.zeros(q, dtype=np.int8)
    units = ((1, 0), (0, 1), (-1, 0), (0, -1))
    value = 1
    from k5_p29_coefficient_sieve import primitive_element

    generator = primitive_element(p)
    for exponent in range(q - 1):
        real[value], imag[value] = units[exponent % 4]
        value = multiply(value, generator)
    z = 0j
    for d in range(1, q):
        z += int(n[d]) * (int(real[d]) + 1j * int(imag[d]))
    return z


def boolean_and_moment(
    subset,
    orbit_size,
    forms,
    coordinates,
    coeff_hits,
    kernel_real,
    kernel_imag,
    p,
    gmin,
):
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    selected_coordinates = [coordinates[i] for i in subset]
    two_p = 2 * p
    expected_a = two_p * gmin
    hist = Counter()
    n_bool = 0
    n_a_ok = 0
    n_bool_a_eq = 0
    n_official_match = 0
    n_checked_official = 0
    naive_signed = []
    first_Dmask = None
    q, multiply, character, _trace = field_ctx(p)
    # g_L: annihilator direction as F_q element from the linear form
    g_vals = []
    for form in [forms[i] for i in subset]:
        cx, cy = int(form[0]), int(form[1])
        g_vals.append(cy + ((-cx) % p) * p)
    for leading, rec in coeff_hits:
        polynomial = np.zeros((7, p), dtype=np.int64)
        for j in range(7):
            a, (c, d, e, f) = leading[j], rec[j]
            polynomial[j] = (a * s**5 + c * s**3 + d * s**2 + e * s + f) % p
        centered = np.where(polynomial <= midpoint, polynomial, polynomial - p)
        replacements = np.sum(centered, axis=1) // p
        endpoint_choices = []
        skip = False
        for j in range(7):
            locs = np.where(polynomial[j] == midpoint)[0]
            r = int(replacements[j])
            if r < 0 or r > len(locs):
                skip = True
                break
            endpoint_choices.append(list(combinations(locs, r)))
        if skip:
            continue
        for endpoint_sets in product(*endpoint_choices):
            profiles = centered.copy()
            for j, endpoint_set in enumerate(endpoint_sets):
                if endpoint_set:
                    profiles[j, list(endpoint_set)] -= p
            a_vec = np.sum(profiles * profiles, axis=1)
            if np.any(a_vec % two_p):
                continue
            n_a_ok += 1
            point_sum = sum(
                profiles[j][selected_coordinates[j]] for j in range(7)
            )
            if not np.all(
                (point_sum == midpoint) | (point_sum == -midpoint - 1)
            ):
                continue
            negative = (point_sum == -midpoint - 1).astype(np.int64)
            real = int(negative @ kernel_real @ negative)
            imag = int(negative @ kernel_imag @ negative)
            val = real * real + imag * imag
            hist[val] += orbit_size
            n_bool += orbit_size
            if first_Dmask is None:
                first_Dmask = negative.copy()
            if np.all(a_vec == expected_a) and not np.any(a_vec % two_p):
                n_bool_a_eq += orbit_size
            if n_checked_official < 3:
                z_off = official_zpsi_from_D(negative, p)
                n_checked_official += 1
                if abs(z_off.real - real) < 1e-6 and abs(z_off.imag - imag) < 1e-6:
                    n_official_match += 1
            if np.all(a_vec == expected_a):
                signed = 0j
                for j, g in enumerate(g_vals):
                    # ψ(g) via quartic table on g; p≡1 so this is NOT Z_ψ
                    signed += a_vec[j] * (
                        kernel_real[g, 0] + 1j * kernel_imag[g, 0]
                    )
                naive_signed.append(signed)
    return {
        "n_boolean": n_bool,
        "hist": hist,
        "n_a_divisible": n_a_ok,
        "n_boolean_a_L_eq_2pgmin": n_bool_a_eq,
        "n_official_checked": n_checked_official,
        "n_official_match": n_official_match,
        "expected_a": expected_a,
        "naive_signed_sample": naive_signed[:3],
        "first_Dmask": first_Dmask,
    }


_CACHE: dict | None = None


def scan_p41_k7(*, write_json: bool = False, recompute: bool | None = None) -> dict:
    """Full coupled sieve + Boolean + |Z_ψ|².

    Cached per process.  If a JSON dump exists and recompute is false,
    load it (leftover dumps / cheap qvar imports).  The AF test sets
    K7_P41_RECOMPUTE=1 so pytest actually runs the sieve.
    """
    global _CACHE
    if recompute is None:
        recompute = os.environ.get("K7_P41_RECOMPUTE") == "1"
    if _CACHE is not None and not recompute:
        return _CACHE
    path = HERE / "k7_p41_coefficient_sieve.json"
    if not recompute and path.is_file():
        _CACHE = json.loads(path.read_text())
        return _CACHE
    p, target = P, GMIN
    T = (p * p - 1) // 8
    q = p * p
    threshold = 3 * q * (q - 1) // 16
    types_by_leading = enumerate_types(p, target)
    types_by_cubic = index_types(types_by_leading)
    n_types = sum(len(v) for v in types_by_leading.values())
    square = square_directions(p)
    coordinates = [c for c, _ in square]
    forms = [f for _, f in square]
    perm = cyclic_direction_permutation(forms, p)
    orbits = subset_orbits(len(forms), K, perm)
    jobs = [
        (subset, orbit_size, forms, types_by_leading, types_by_cubic, p)
        for subset, orbit_size in orbits
    ]
    workers = _n_workers(80)
    hit_orbits = []
    n_coeff = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for subset, orbit_size, n_c, hits in pool.map(
            sieve_one_orbit, jobs, chunksize=1
        ):
            n_coeff += n_c * orbit_size
            if n_c:
                hit_orbits.append(
                    {
                        "subset": list(subset),
                        "orbit_size": orbit_size,
                        "n_coeff_rep": n_c,
                        "hits": hits,
                    }
                )
    kernel_real, kernel_imag = quartic_kernel(p)
    hist = Counter()
    n_bool = 0
    n_official_checked = n_official_match = 0
    n_a_divisible = 0
    n_bool_a_eq = 0
    first_Dmask = None
    expected_a = 2 * p * target
    for orb in hit_orbits:
        part = boolean_and_moment(
            tuple(orb["subset"]),
            orb["orbit_size"],
            forms,
            coordinates,
            orb["hits"],
            kernel_real,
            kernel_imag,
            p,
            target,
        )
        hist.update(part["hist"])
        n_bool += part["n_boolean"]
        n_official_checked += part["n_official_checked"]
        n_official_match += part["n_official_match"]
        n_a_divisible += part["n_a_divisible"]
        n_bool_a_eq += part["n_boolean_a_L_eq_2pgmin"]
        if first_Dmask is None:
            first_Dmask = part.get("first_Dmask")
        orb.pop("hits", None)
        orb["n_boolean_weighted"] = part["n_boolean"]
    moment = (
        sum(value * count for value, count in hist.items()) / n_bool
        if n_bool
        else None
    )
    clears = bool(n_bool == 0 or (moment is not None and moment >= threshold))
    report = {
        "p": p,
        "k": K,
        "normalized_total_T": T,
        "minimum_quintic_b": target,
        "unique_energy_partition": [target] * K,
        "n_square_directions": len(forms),
        "n_direction_subsets": sum(s for _, s in orbits),
        "n_cyclic_subset_orbits": len(orbits),
        "n_min_energy_leadings": len(types_by_leading),
        "min_energy_leadings": sorted(types_by_leading),
        "n_types_at_min": n_types,
        "n_types_by_leading": {
            str(L): len(v) for L, v in sorted(types_by_leading.items())
        },
        "total_coefficient_candidates": n_coeff,
        "n_hit_orbits": len(hit_orbits),
        "hit_orbits": [
            {
                "subset": o["subset"],
                "orbit_size": o["orbit_size"],
                "n_coeff_rep": o["n_coeff_rep"],
                "n_boolean_weighted": o["n_boolean_weighted"],
            }
            for o in hit_orbits
        ],
        "boolean_representatives_mod_translation": n_bool,
        "abs_Zpsi_sq_histogram": {str(v): c for v, c in sorted(hist.items())},
        "E_abs_Zpsi_sq": 0 if moment == 0 else moment,
        "QVAR_threshold": threshold,
        "clears_QVAR": clears,
        "k7_empty": n_bool == 0,
        "stratum_qvar": bool(clears and n_bool > 0) or (n_bool == 0),
        "a_L_expected": expected_a,
        "a_L_in_2pZ_on_endpoint_tuples": n_a_divisible > 0,
        "boolean_a_L_eq_2p_gmin": n_bool_a_eq,
        "boolean_a_L_all_match": n_bool > 0 and n_bool_a_eq == n_bool,
        "official_Z_checked": n_official_checked,
        "official_Z_matches_kernel": (
            n_official_checked > 0 and n_official_match == n_official_checked
        ),
        "p_mod_4": p % 4,
        "integer_signed_sum_is_Z": p % 4 == 3,
        "maxplus_Cy_eq_py": False,
        "maxplus_residual": None,
    }
    if first_Dmask is not None:
        sys.path.insert(0, str(REPO / "src"))
        from minmax_quadratic import paley_conference_prime_power

        C = paley_conference_prime_power(p)
        y = np.ones(q + 1, dtype=np.float64)
        y[1:] = np.where(first_Dmask.astype(bool), -1.0, 1.0)
        resid = float(np.max(np.abs(C @ y - p * y)))
        report["maxplus_residual"] = resid
        report["maxplus_Cy_eq_py"] = resid < 1e-8
    # Vacuous empty strata clear; nonempty with E below floor do not.
    if n_bool > 0:
        report["stratum_qvar"] = bool(moment >= threshold)
    if write_json or recompute:
        slim = dict(report)
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(slim, indent=2) + "\n")
        tmp.replace(path)
    _CACHE = report
    return report


def stratum_clears_qvar() -> bool:
    return bool(scan_p41_k7()["stratum_qvar"])


if __name__ == "__main__":
    out = scan_p41_k7(write_json=True)
    print(json.dumps({k: v for k, v in out.items() if k != "hit_orbits"}, indent=2))
    print("hit_orbits", out["hit_orbits"], flush=True)
    print(
        f"E={out['E_abs_Zpsi_sq']} thr={out['QVAR_threshold']} "
        f"clears={out['stratum_qvar']} empty={out['k7_empty']} "
        f"official_match={out['official_Z_matches_kernel']}",
        flush=True,
    )
