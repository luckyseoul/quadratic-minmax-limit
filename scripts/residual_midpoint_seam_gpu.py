#!/usr/bin/env python3
"""GPU probes for the uncollapsed midpoint seam in residual (ii).

This is discovery code, not theorem evidence.  It uses the *anti-diagonal*
coefficient profile

    A_L(a) = sum_{s<t, (s+t)/2=a} W_L(s,t),

never the endpoint row sum.  If D_L(a) records parallel edges by midpoint
line, then G_L(a)=eta_L(D_L(a)+A_L(a)) must be an integral affine point-
Radon transform.  Equivalently

    sum_L G_L(Lx) == T (mod p)                 for every x in F_p^2.

The ``anneal`` command searches this new necessary condition for the compact
all-prime atom templates of Proposition 15.758.  The ``identity`` command
checks the coefficientwise seam on an actual random simple graph at a large
held-out prime.  Neither command is a finite-prime census or a graph lift.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import time
from pathlib import Path

import numpy as np


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def legendre_table(p: int) -> np.ndarray:
    out = np.zeros(p, dtype=np.int32)
    for a in range(1, p):
        out[a] = 1 if pow(a, (p - 1) // 2, p) == 1 else -1
    return out


def smallest_nonsquare(p: int) -> int:
    legendre = legendre_table(p)
    return next(a for a in range(2, p) if legendre[a] == -1)


def directions_and_signs(p: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Return functionals (1,lambda), (0,1) and Paley line signs."""
    nonsquare = smallest_nonsquare(p)
    legendre = legendre_table(p)
    directions = np.array([(1, lam) for lam in range(p)] + [(0, 1)], dtype=np.int32)
    signs = []
    for left, right in directions:
        # Kernel vector (right,-left); norm x^2-d*y^2.
        norm = (int(right) ** 2 - nonsquare * int(left) ** 2) % p
        signs.append(int(legendre[norm]))
    signs_array = np.array(signs, dtype=np.int32)
    if np.count_nonzero(signs_array == 1) != (p + 1) // 2:
        raise ArithmeticError("projective Paley signs did not split evenly")
    return directions, signs_array, nonsquare


def balanced(total: int, count: int, baseline: int = 0) -> list[int]:
    q, rem = divmod(total - baseline * count, count)
    if q < 0:
        raise ValueError("balanced total lies below baseline")
    return [baseline + q + int(i < rem) for i in range(count)]


def add_edge_midpoint(row: np.ndarray, p: int, u: int, v: int, coefficient: int) -> None:
    row[((u + v) * pow(2, -1, p)) % p] += coefficient


def add_star(row: np.ndarray, p: int, centre: int, coefficient: int) -> None:
    for other in range(p):
        if other != centre:
            add_edge_midpoint(row, p, centre, other, coefficient)


def random_triple(rng: random.Random, p: int) -> tuple[int, int, int]:
    a, b, c = rng.sample(range(p), 3)
    return a, b, c


def add_triangle(row: np.ndarray, p: int, triple: tuple[int, int, int], kind: str) -> None:
    a, b, c = triple
    if kind == "plus":
        for u, v in ((a, b), (a, c), (b, c)):
            add_edge_midpoint(row, p, u, v, 1)
    elif kind == "minus":
        add_edge_midpoint(row, p, a, b, 1)
        add_edge_midpoint(row, p, a, c, -1)
        add_edge_midpoint(row, p, b, c, -1)
    else:
        raise ValueError(kind)


def add_omitted_pair(row: np.ndarray, p: int, a: int, b: int) -> None:
    add_star(row, p, a, -1)
    add_star(row, p, b, -1)
    add_edge_midpoint(row, p, a, b, 1)


def compact_template(p: int, branch: str, endpoint: str, seed: int) -> dict[str, object]:
    if not is_prime(p) or p < 29:
        raise ValueError("probe requires a held-out prime p>=29")
    directions, eta, nonsquare = directions_and_signs(p)
    hard = np.flatnonzero(eta == 1).tolist()
    opposite = np.flatnonzero(eta == -1).tolist()
    m = (p + 1) // 2
    rng = random.Random(seed)
    anti = np.zeros((p + 1, p), dtype=np.int32)
    parallel = np.zeros(p + 1, dtype=np.int32)

    if branch == "B":
        if p % 4 != 1:
            raise ValueError("branch B requires p=1 mod 4")
        r = (p - 1) // 4
        lower, upper = 2 * r * r - 5 * r, 4 * r * r - 6 * r - 3
        t = lower if endpoint == "lower" else upper
        excesses = balanced(t, m)
        q_values = balanced(6 * r + t, m, r)
        signed_total = p + 4
        for direction, excess in zip(hard, excesses):
            add_star(anti[direction], p, rng.randrange(p), 1)
            for _ in range(excess):
                add_triangle(anti[direction], p, random_triple(rng, p), "minus")
            parallel[direction] = 5 + excess
        for direction, q_value in zip(opposite, q_values):
            a, b = rng.sample(range(p), 2)
            add_omitted_pair(anti[direction], p, a, b)
            for _ in range(r - 2):
                add_triangle(anti[direction], p, random_triple(rng, p), "plus")
            for _ in range(q_value - r):
                add_triangle(anti[direction], p, random_triple(rng, p), "minus")
            parallel[direction] = q_value
    elif branch == "C":
        if p % 4 != 3:
            raise ValueError("branch C requires p=3 mod 4")
        r = (p - 3) // 4
        lower, upper = 2 * r * r - 4 * r - 2, 4 * r * r - 2 * r - 5
        t = lower if endpoint == "lower" else upper
        excesses = balanced(t + 1, m)
        q_values = balanced(10 * r + 6 + t, m, r + 2)
        signed_total = 4 - p
        for direction, excess in zip(hard, excesses):
            add_star(anti[direction], p, rng.randrange(p), -1)
            for _ in range(excess):
                add_triangle(anti[direction], p, random_triple(rng, p), "minus")
            parallel[direction] = 3 + excess
        for direction, q_value in zip(opposite, q_values):
            for _ in range(r - 1):
                add_triangle(anti[direction], p, random_triple(rng, p), "plus")
            for _ in range(q_value - r - 2):
                add_triangle(anti[direction], p, random_triple(rng, p), "minus")
            parallel[direction] = q_value
    else:
        raise ValueError("branch must be B or C")

    expected_offdiag = eta * signed_total - parallel
    actual_offdiag = anti.sum(axis=1)
    if not np.array_equal(actual_offdiag, expected_offdiag):
        raise ArithmeticError("compact atom anti-diagonal totals changed")
    edge_count = 4 * p + 2 * t + 1
    if int(parallel.sum()) != edge_count:
        raise ArithmeticError("parallel counts do not sum to |H|")

    particle_directions = np.repeat(np.arange(p + 1, dtype=np.int32), parallel)
    points = np.array([(x, y) for x in range(p) for y in range(p)], dtype=np.int32)
    labels = (directions[:, 0, None] * points[:, 0] + directions[:, 1, None] * points[:, 1]) % p
    base = (np.sum(eta[:, None] * anti[np.arange(p + 1)[:, None], labels], axis=0) - signed_total) % p
    payload = {
        "p": p,
        "branch": branch,
        "endpoint": endpoint,
        "seed": seed,
        "nonsquare": nonsquare,
        "eta": eta.tolist(),
        "parallel": parallel.tolist(),
        "anti": anti.tolist(),
        "signed_total": signed_total,
        "t": t,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return {
        **payload,
        "directions": directions,
        "particle_directions": particle_directions,
        "base": base.astype(np.int32),
        "edge_count": edge_count,
        "template_sha256": hashlib.sha256(canonical).hexdigest(),
    }


ANNEAL_SOURCE = r'''
extern "C" __global__
void init_residual(const int *base, const int *particle_dir, const int *eta,
                   const int *state, int *residual, int chains, int points,
                   int total, int p) {
    int index = blockDim.x * blockIdx.x + threadIdx.x;
    int limit = chains * points;
    if (index >= limit) return;
    int chain = index / points;
    int point = index - chain * points;
    int x = point / p;
    int y = point - x * p;
    int value = base[point];
    for (int k = 0; k < total; ++k) {
        int d = particle_dir[k];
        int label = d == p ? y : (x + d * y) % p;
        value += eta[d] * (label == state[chain * total + k]);
    }
    value %= p;
    if (value < 0) value += p;
    residual[index] = value;
}

__device__ __forceinline__ unsigned int xs32(unsigned int &x) {
    x ^= x << 13; x ^= x >> 17; x ^= x << 5; return x;
}
__device__ __forceinline__ int centre(int x, int p) {
    return x <= p/2 ? x : x-p;
}
__device__ __forceinline__ int line_point(int d, int a, int q, int p) {
    if (d == p) return q*p + a;
    int x = a - (d*q)%p; if (x < 0) x += p;
    return x*p + q;
}

extern "C" __global__
void anneal(const int *particle_dir, const int *eta, int *state, int *residual,
            long long *objective, long long *best_seen, int chains, int points,
            int total, int p, int steps, unsigned int seed, float temp0) {
    int chain = blockDim.x * blockIdx.x + threadIdx.x;
    if (chain >= chains) return;
    unsigned int rng = seed ^ (0x9e3779b9U * (chain + 1));
    long long obj = objective[chain], best = obj;
    int roff = chain * points, soff = chain * total;
    for (int step = 0; step < steps && best != 0; ++step) {
        int k = (int)(xs32(rng) % (unsigned int)total);
        int olda = state[soff+k];
        int newa = (int)(xs32(rng) % (unsigned int)(p-1));
        if (newa >= olda) ++newa;
        int d = particle_dir[k], sgn = eta[d];
        long long change = 0;
        for (int q=0; q<p; ++q) {
            int idx = roff + line_point(d, olda, q, p);
            int before = residual[idx], after = before - sgn;
            after %= p; if (after < 0) after += p;
            int cb=centre(before,p), ca=centre(after,p);
            change += (long long)ca*ca-(long long)cb*cb;
        }
        for (int q=0; q<p; ++q) {
            int idx = roff + line_point(d, newa, q, p);
            int before = residual[idx], after = before + sgn;
            after %= p; if (after < 0) after += p;
            int cb=centre(before,p), ca=centre(after,p);
            change += (long long)ca*ca-(long long)cb*cb;
        }
        float temp = temp0 * (1.0f - (float)step/(float)steps) + 0.05f;
        float u = (float)(xs32(rng) & 0x00ffffffU) / 16777216.0f;
        bool accept = change <= 0 || u < expf(-(float)change/temp);
        if (!accept) continue;
        for (int q=0; q<p; ++q) {
            int idx = roff + line_point(d, olda, q, p);
            int after = (residual[idx] - sgn) % p;
            if (after < 0) after += p; residual[idx] = after;
        }
        for (int q=0; q<p; ++q) {
            int idx = roff + line_point(d, newa, q, p);
            int after = (residual[idx] + sgn) % p;
            if (after < 0) after += p; residual[idx] = after;
        }
        state[soff+k] = newa;
        obj += change; if (obj < best) best = obj;
    }
    objective[chain] = obj; best_seen[chain] = best;
}
'''


def centred_square_sum(residue: np.ndarray, p: int) -> int:
    centred = np.where(residue <= p // 2, residue, residue - p)
    return int(np.dot(centred.astype(np.int64), centred.astype(np.int64)))


def anneal(args: argparse.Namespace) -> dict[str, object]:
    import cupy as cp

    template = compact_template(args.p, args.branch, args.endpoint, args.template_seed)
    p = args.p
    chains = args.chains
    total = len(template["particle_directions"])
    cp.random.seed(args.search_seed)
    state = cp.random.randint(0, p, size=(chains, total), dtype=cp.int32)
    base = cp.asarray(template["base"], dtype=cp.int32)
    particle_dir = cp.asarray(template["particle_directions"], dtype=cp.int32)
    eta = cp.asarray(template["eta"], dtype=cp.int32)
    residual = cp.empty((chains, p * p), dtype=cp.int32)
    init = cp.RawKernel(ANNEAL_SOURCE, "init_residual")
    kernel = cp.RawKernel(ANNEAL_SOURCE, "anneal")
    block = 256
    init(((chains * p * p + block - 1) // block,), (block,),
         (base, particle_dir, eta, state, residual, chains, p*p, total, p))
    centred = cp.where(residual <= p // 2, residual, residual - p).astype(cp.int64)
    objective = cp.sum(centred * centred, axis=1, dtype=cp.int64)
    initial_min = int(cp.min(objective).get())
    best_seen = cp.empty_like(objective)
    start = time.time()
    kernel(((chains + block - 1)//block,), (block,),
           (particle_dir, eta, state, residual, objective, best_seen,
            chains, p*p, total, p, args.steps, args.search_seed, args.temperature))
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.time() - start
    final_index = int(cp.argmin(objective).get())
    final_objective = int(objective[final_index].get())
    best_value = int(cp.min(best_seen).get())
    final_state = cp.asnumpy(state[final_index])

    exact = np.array(template["base"], copy=True)
    directions = np.asarray(template["directions"])
    points = np.array([(x, y) for x in range(p) for y in range(p)], dtype=np.int32)
    labels = (directions[:, 0, None] * points[:, 0] + directions[:, 1, None] * points[:, 1]) % p
    for atom, direction in enumerate(template["particle_directions"]):
        exact[labels[int(direction)] == int(final_state[atom])] += int(template["eta"][int(direction)])
    exact %= p
    exact_objective = centred_square_sum(exact, p)
    if exact_objective != final_objective:
        raise ArithmeticError("CPU exact replay disagrees with GPU objective")
    device = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())
    device_name = device["name"].decode() if isinstance(device["name"], bytes) else str(device["name"])
    return {
        "command": "anneal",
        "classification": "randomized necessary-condition probe; not theorem evidence",
        "host": platform.node(),
        "architecture": platform.machine(),
        "gpu_backend": "cupy",
        "gpu_device": device_name,
        "p": p,
        "branch": args.branch,
        "endpoint": args.endpoint,
        "t": template["t"],
        "edge_count": template["edge_count"],
        "signed_total": template["signed_total"],
        "template_seed": args.template_seed,
        "search_seed": args.search_seed,
        "template_sha256": template["template_sha256"],
        "chains": chains,
        "steps_per_chain": args.steps,
        "initial_minimum_centered_residue_square": initial_min,
        "minimum_seen_centered_residue_square": best_value,
        "final_replay_centered_residue_square": exact_objective,
        "final_nonzero_point_residues": int(np.count_nonzero(exact)),
        "integral_point_radon_completion_found": exact_objective == 0,
        "elapsed_seconds": elapsed,
        "changed_premise": "uncollapsed anti-diagonal midpoint profiles and hidden diagonal counts",
    }


def identity(args: argparse.Namespace) -> dict[str, object]:
    import cupy as cp

    p = args.p
    if not is_prime(p) or p < 29:
        raise ValueError("identity probe requires a prime p>=29")
    directions, eta, nonsquare = directions_and_signs(p)
    legendre = legendre_table(p)
    rng = random.Random(args.template_seed)
    edge_set: set[tuple[int, int]] = set()
    point_count = p * p
    while len(edge_set) < args.edges:
        u, v = rng.sample(range(point_count), 2)
        edge_set.add((min(u, v), max(u, v)))
    edges = np.asarray(sorted(edge_set), dtype=np.int32)
    ux, uy = edges[:, 0] // p, edges[:, 0] % p
    vx, vy = edges[:, 1] // p, edges[:, 1] % p
    dx, dy = (vx - ux) % p, (vy - uy) % p
    tau = legendre[(dx * dx - nonsquare * dy * dy) % p]
    inv2 = pow(2, -1, p)
    mx, my = ((ux + vx) * inv2) % p, ((uy + vy) * inv2) % p

    # GPU scatter builds the two halves independently.
    c_tau = cp.asarray(tau)
    c_mindex = cp.asarray(mx * p + my)
    g = cp.zeros(point_count, dtype=cp.int64)
    cp.add.at(g, c_mindex, c_tau)
    diagonal = cp.zeros((p + 1, p), dtype=cp.int64)
    anti = cp.zeros((p + 1, p), dtype=cp.int64)
    endpoint = cp.zeros((p + 1, p), dtype=cp.int64)
    c_ux, c_uy, c_vx, c_vy = map(cp.asarray, (ux, uy, vx, vy))
    for index, (left, right) in enumerate(directions):
        su = (int(left) * c_ux + int(right) * c_uy) % p
        sv = (int(left) * c_vx + int(right) * c_vy) % p
        midpoint_label = ((su + sv) * inv2) % p
        mask = su == sv
        cp.add.at(diagonal[index], su[mask], 1)
        cp.add.at(anti[index], midpoint_label[~mask], int(eta[index]) * c_tau[~mask])
        coefficient = int(eta[index]) * c_tau[~mask]
        cp.add.at(endpoint[index], su[~mask], coefficient)
        cp.add.at(endpoint[index], sv[~mask], coefficient)
    G = cp.asarray(eta[:, None]) * (diagonal + anti)
    direct = cp.zeros_like(G)
    points_x = cp.repeat(cp.arange(p, dtype=cp.int32), p)
    points_y = cp.tile(cp.arange(p, dtype=cp.int32), p)
    for index, (left, right) in enumerate(directions):
        label = (int(left) * points_x + int(right) * points_y) % p
        cp.add.at(direct[index], label, g)
    discrepancy = cp.asnumpy(G - direct)
    totals = cp.asnumpy(cp.sum(G, axis=1))
    reconstructed = cp.zeros(point_count, dtype=cp.int64)
    for index, (left, right) in enumerate(directions):
        label = (int(left) * points_x + int(right) * points_y) % p
        reconstructed += G[index, label]
    reconstructed = (reconstructed - int(cp.sum(g).get())) // p
    reconstruction_error = cp.asnumpy(reconstructed - g)

    # The endpoint and midpoint halves cancel the unknown parallel-edge
    # diagonal exactly.  If d_z is the signed endpoint degree function and
    # kappa=d_z-2g, then H=eta*(E-2A)=Radon(kappa).
    signed_degree = cp.zeros(point_count, dtype=cp.int64)
    cp.add.at(signed_degree, cp.asarray(edges[:, 0]), c_tau)
    cp.add.at(signed_degree, cp.asarray(edges[:, 1]), c_tau)
    kappa = signed_degree - 2 * g
    H = cp.asarray(eta[:, None]) * (endpoint - 2 * anti)
    direct_kappa = cp.zeros_like(H)
    for index, (left, right) in enumerate(directions):
        label = (int(left) * points_x + int(right) * points_y) % p
        cp.add.at(direct_kappa[index], label, kappa)
    curvature_discrepancy = cp.asnumpy(H - direct_kappa)
    reconstructed_kappa_numerator = cp.zeros(point_count, dtype=cp.int64)
    for index, (left, right) in enumerate(directions):
        label = (int(left) * points_x + int(right) * points_y) % p
        reconstructed_kappa_numerator += H[index, label]
    curvature_inversion_remainder = cp.asnumpy(reconstructed_kappa_numerator % p)
    reconstructed_kappa = reconstructed_kappa_numerator // p
    curvature_inversion_error = cp.asnumpy(reconstructed_kappa - kappa)
    kappa_host = cp.asnumpy(kappa)
    signed_degree_host = cp.asnumpy(signed_degree)
    curvature_energy = int(cp.sum(H * H).get())
    kappa_energy = int(cp.sum(kappa * kappa).get())
    boundary_size = int(np.count_nonzero(signed_degree_host % 2))
    cp.cuda.runtime.deviceSynchronize()
    device = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())
    device_name = device["name"].decode() if isinstance(device["name"], bytes) else str(device["name"])
    edge_hash = hashlib.sha256(edges.tobytes()).hexdigest()
    return {
        "command": "identity",
        "classification": "exact coefficientwise identity validation; theorem still requires symbolic proof",
        "host": platform.node(),
        "architecture": platform.machine(),
        "gpu_backend": "cupy",
        "gpu_device": device_name,
        "p": p,
        "random_simple_edges": len(edges),
        "edge_list_sha256": edge_hash,
        "template_seed": args.template_seed,
        "maximum_anti_diagonal_identity_error": int(np.max(np.abs(discrepancy))),
        "directional_totals_all_equal": bool(np.all(totals == totals[0])),
        "common_total": int(totals[0]),
        "maximum_exact_inversion_error": int(np.max(np.abs(reconstruction_error))),
        "anti_diagonal_seam_validated": bool(not np.any(discrepancy) and not np.any(reconstruction_error)),
        "formula": "G_L(a)=eta_L[D_L(a)+sum_(s+t=2a)W_L(s,t)]",
        "inverse": "p*g(x)=sum_L G_L(Lx)-T",
        "maximum_endpoint_midpoint_curvature_error": int(np.max(np.abs(curvature_discrepancy))),
        "maximum_curvature_inversion_error": int(np.max(np.abs(curvature_inversion_error))),
        "nonzero_curvature_inversion_remainders_mod_p": int(np.count_nonzero(curvature_inversion_remainder)),
        "curvature_total_zero": bool(int(kappa_host.sum()) == 0),
        "curvature_parity_equals_graph_boundary": bool(np.array_equal(kappa_host % 2, signed_degree_host % 2)),
        "graph_boundary_size": boundary_size,
        "curvature_energy": curvature_energy,
        "p_times_point_curvature_energy": p * kappa_energy,
        "parseval_exact": bool(curvature_energy == p * kappa_energy),
        "energy_dominates_p_times_boundary": bool(curvature_energy >= p * boundary_size),
        "endpoint_midpoint_curvature_seam_validated": bool(
            not np.any(curvature_discrepancy)
            and not np.any(curvature_inversion_error)
            and not np.any(curvature_inversion_remainder)
            and int(kappa_host.sum()) == 0
            and np.array_equal(kappa_host % 2, signed_degree_host % 2)
            and curvature_energy == p * kappa_energy
            and curvature_energy >= p * boundary_size
        ),
        "curvature_formula": "H_L=eta_L*(E_L-2*A_L)=Radon(kappa)_L, kappa=d_z-2*g",
        "curvature_inverse": "p*kappa(x)=sum_L H_L(Lx)",
        "curvature_parseval": "sum_(L,a) H_L(a)^2=p*sum_x kappa(x)^2>=p*|boundary|",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    anneal_parser = sub.add_parser("anneal")
    anneal_parser.add_argument("--p", type=int, required=True)
    anneal_parser.add_argument("--branch", choices=("B", "C"), required=True)
    anneal_parser.add_argument("--endpoint", choices=("lower", "upper"), default="lower")
    anneal_parser.add_argument("--template-seed", type=int, default=15762)
    anneal_parser.add_argument("--search-seed", type=int, default=1)
    anneal_parser.add_argument("--chains", type=int, default=1024)
    anneal_parser.add_argument("--steps", type=int, default=20000)
    anneal_parser.add_argument("--temperature", type=float, default=8.0)
    anneal_parser.add_argument("--output", type=Path, required=True)
    identity_parser = sub.add_parser("identity")
    identity_parser.add_argument("--p", type=int, default=127)
    identity_parser.add_argument("--edges", type=int, default=8192)
    identity_parser.add_argument("--template-seed", type=int, default=15762)
    identity_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = anneal(args) if args.command == "anneal" else identity(args)
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
