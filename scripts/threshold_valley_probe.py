#!/usr/bin/env python3
"""Bounded n=6 generic threshold-valley countermodel probe.

The GPU computes every complete-signing score against all 32 antipodal
Boolean representatives; BOTH energy signs remain in every norm/active set.
At most 128 deterministic reference signings are searched. This is not a
Paley counterexample, a conference construction, or a limit assertion.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import random
import sys
import time


N = 6
EDGES = tuple((i, j) for i in range(N) for j in range(i + 1, N))
EDGE_COUNT = len(EDGES)
SIGNING_COUNT = 1 << EDGE_COUNT
STATE_COUNT = 1 << (N - 1)
MAX_REFERENCES = 128
DEFAULT_SEED = 20260904


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def reference_masks(seed: int = DEFAULT_SEED, count: int = MAX_REFERENCES) -> list[int]:
    require(type(seed) is int, "seed must be an integer")
    require(type(count) is int and 1 <= count <= MAX_REFERENCES,
            "reference count must be between 1 and 128; expansion is not authorized")
    return random.Random(seed).sample(range(SIGNING_COUNT), count)


def check_mask(mask: int) -> None:
    require(type(mask) is int and 0 <= mask < SIGNING_COUNT, "invalid n=6 signing mask")


def state_vector(index: int) -> tuple[int, ...]:
    require(type(index) is int and 0 <= index < STATE_COUNT, "invalid Boolean state index")
    return (1,) + tuple(1 - 2 * ((index >> (i - 1)) & 1) for i in range(1, N))


def signing_matrix(mask: int) -> list[list[int]]:
    check_mask(mask)
    matrix = [[0] * N for _ in range(N)]
    for bit, (i, j) in enumerate(EDGES):
        matrix[i][j] = matrix[j][i] = 1 - 2 * ((mask >> bit) & 1)
    return matrix


def integer_scores(mask: int) -> tuple[int, ...]:
    """Independent plain-integer score implementation for small witnesses."""
    matrix = signing_matrix(mask)
    return tuple(sum(matrix[i][j] * x[i] * x[j] for i, j in EDGES)
                 for x in (state_vector(index) for index in range(STATE_COUNT)))


def active_rows(scores: tuple[int, ...]) -> list[tuple[int, int]]:
    norm = max(abs(value) for value in scores)
    return [(index, 1 if value > 0 else -1)
            for index, value in enumerate(scores) if abs(value) == norm]


def rational_valley_lambda(a_scores: tuple[int, ...], c_scores: tuple[int, ...]) -> Fraction:
    """Strict descent on every signed A-active row gives an explicit interval."""
    require(len(a_scores) == len(c_scores) == STATE_COUNT, "need all 32 state scores")
    norm = max(abs(value) for value in a_scores)
    require(all(sign * c_scores[index] < norm for index, sign in active_rows(a_scores)),
            "not every signed A-active row descends")
    endpoint = Fraction(1)
    for a_value, c_value in zip(a_scores, c_scores):
        for sign in (-1, 1):
            intercept = sign * a_value
            slope = sign * (c_value - a_value)
            if slope > 0:
                require(intercept < norm, "a positive-slope active row prevents a valley")
                endpoint = min(endpoint, Fraction(norm - intercept, slope))
    result = endpoint / 2
    require(0 < result < 1, "valley parameter is not interior")
    require(max(abs((1 - result) * a + result * c)
                for a, c in zip(a_scores, c_scores)) < norm,
            "rational parameter failed the full signed norm check")
    return result


def restoration_masks(a_mask: int, c_mask: int):
    """All nonempty restorations from A toward C, including C exactly once."""
    check_mask(a_mask)
    check_mask(c_mask)
    disagreement = a_mask ^ c_mask
    subset = disagreement
    while subset:
        yield subset, a_mask ^ subset
        subset = (subset - 1) & disagreement


def verify_witness(witness: dict[str, object]) -> dict[str, object]:
    """Exact CPU verification; no GPU table or search verdict is trusted."""
    a_mask, c_mask = witness["A_mask"], witness["C_mask"]
    check_mask(a_mask)
    check_mask(c_mask)
    a_scores, c_scores = integer_scores(a_mask), integer_scores(c_mask)
    a_norm = max(map(abs, a_scores))
    c_norm = max(map(abs, c_scores))
    disagreement = a_mask ^ c_mask
    parameter = Fraction(witness["lambda"])
    interior_scores = tuple((1 - parameter) * a + parameter * c
                            for a, c in zip(a_scores, c_scores))
    interior_norm = max(map(abs, interior_scores))
    a_active, c_active = active_rows(a_scores), active_rows(c_scores)
    restored_norms = Counter()
    restored_sizes = Counter()
    failing_restorations = []
    for subset, mask in restoration_masks(a_mask, c_mask):
        norm = max(map(abs, integer_scores(mask)))
        restored_norms[norm] += 1
        restored_sizes[subset.bit_count()] += 1
        if norm < a_norm + 2 and len(failing_restorations) < 8:
            failing_restorations.append({"subset_mask": subset, "restored_mask": mask,
                                         "norm": norm})
    checks = {
        "n_is_six": witness.get("n") == N,
        "nonempty_disagreement": disagreement != 0,
        "disagreement_mask_matches": witness.get("H_mask") == disagreement,
        "complete_A_matrix_matches": witness.get("A_matrix") == signing_matrix(a_mask),
        "complete_C_matrix_matches": witness.get("C_matrix") == signing_matrix(c_mask),
        "Phi_A_equals_Phi_C_minus_four": a_norm == c_norm - 4,
        "reported_M_matches": witness.get("M") == a_norm,
        "all_signed_A_active_rows_descend": all(
            sign * c_scores[index] < a_norm for index, sign in a_active),
        "all_signed_C_maximizers_drop_by_at_least_four": all(
            sign * a_scores[index] <= a_norm - 4 for index, sign in c_active),
        "all_nonempty_restorations_checked": sum(restored_norms.values()) == (1 << disagreement.bit_count()) - 1,
        "all_nonempty_restorations_have_norm_at_least_M_plus_two": not failing_restorations,
        "lambda_is_strictly_interior": 0 < parameter < 1,
        "interior_full_norm_strictly_below_M": interior_norm < a_norm,
    }
    return {"verifier": "plain Python integer and Fraction arithmetic; no GPU table",
            "host": platform.node(), "architecture": platform.machine(),
            "python": platform.python_version(), "A_mask": a_mask, "C_mask": c_mask,
            "H_mask": disagreement, "H_size": disagreement.bit_count(),
            "Phi_A": a_norm, "Phi_C": c_norm, "lambda": str(parameter),
            "interior_norm": str(interior_norm),
            "A_scores_all_32_states": list(a_scores),
            "C_scores_all_32_states": list(c_scores),
            "interior_scores_all_32_states": [str(value) for value in interior_scores],
            "A_active_signed_rows": [
                {"state_index": index, "state": list(state_vector(index)), "sign": sign,
                 "signed_A_score": sign * a_scores[index], "signed_C_score": sign * c_scores[index]}
                for index, sign in a_active],
            "C_active_signed_rows": [
                {"state_index": index, "state": list(state_vector(index)), "sign": sign,
                 "signed_A_score": sign * a_scores[index], "signed_C_score": sign * c_scores[index]}
                for index, sign in c_active],
            "restoration_count": sum(restored_norms.values()),
            "restoration_norm_histogram": {str(key): value for key, value in sorted(restored_norms.items())},
            "restoration_subset_size_histogram": {str(key): value for key, value in sorted(restored_sizes.items())},
            "minimum_restoration_norm": min(restored_norms, default=None),
            "failing_restorations": failing_restorations,
            "coverage": {"antipodal_representatives": 32, "signed_rows": 64,
                         "global_state_negation_only_quotiented": True,
                         "both_energy_signs_retained": True},
            "checks": checks, "verified": all(checks.values())}


def labelled_switch_classification(delta_mask: int) -> dict[str, object]:
    """Vertex switching and optional global sign, with vertex labels fixed."""
    check_mask(delta_mask)
    cut_states = []
    complemented_cut_states = []
    for state_index in range(STATE_COUNT):
        state = state_vector(state_index)
        cut = sum(1 << bit for bit, (i, j) in enumerate(EDGES) if state[i] != state[j])
        if delta_mask == cut:
            cut_states.append(state_index)
        if delta_mask == (SIGNING_COUNT - 1) ^ cut:
            complemented_cut_states.append(state_index)
    return {"vertex_labels_fixed": True, "vertex_permutations_included": False,
            "is_cut": bool(cut_states), "is_complemented_cut": bool(complemented_cut_states),
            "cut_state_indices": cut_states,
            "complemented_cut_state_indices": complemented_cut_states,
            "same_labelled_signed_switch_orbit": bool(cut_states or complemented_cut_states)}


def conference_structure(mask: int) -> dict[str, object]:
    """Exact one-matrix check, not a conference signing search."""
    matrix = signing_matrix(mask)
    square = [[sum(matrix[i][k] * matrix[k][j] for k in range(N))
               for j in range(N)] for i in range(N)]
    return {"mask": mask, "matrix_square": square,
            "square_equals_five_identity": all(
                square[i][j] == (5 if i == j else 0) for i in range(N) for j in range(N))}


def verify_repair(repair: dict[str, object]) -> dict[str, object]:
    """Replay ONE near miss and its minimum restoration using integer arithmetic."""
    exact = verify_witness(repair)
    a_mask, c_mask, subset = repair["A_mask"], repair["C_mask"], repair["D_mask"]
    check_mask(subset)
    target = exact["Phi_A"]
    low = []
    for other_subset, mask in restoration_masks(a_mask, c_mask):
        norm = max(map(abs, integer_scores(mask)))
        if norm <= target:
            low.append((other_subset.bit_count(), other_subset, mask, norm))
    low.sort()
    require(bool(low), "near miss contains no integral repair")
    minimum_size, minimum_subset, restored_mask, restored_norm = low[0]
    orbit = labelled_switch_classification(subset)
    checks = {key: value for key, value in exact["checks"].items()
              if key != "all_nonempty_restorations_have_norm_at_least_M_plus_two"}
    checks.update({
        "original_all_subset_gate_fails": not exact["checks"][
            "all_nonempty_restorations_have_norm_at_least_M_plus_two"],
        "D_is_nonempty_subset_of_H": subset != 0 and subset & ~(a_mask ^ c_mask) == 0,
        "D_is_first_minimum_size_low_restoration": subset == minimum_subset,
        "D_size_matches": repair.get("D_size") == minimum_size,
        "D_edges_match": repair.get("D_edges") == [
            list(pair) for bit, pair in enumerate(EDGES) if subset & (1 << bit)],
        "restored_mask_matches": repair.get("restored_A_mask") == restored_mask,
        "restored_matrix_matches": repair.get("restored_A_matrix") == signing_matrix(restored_mask),
        "restored_norm_matches": repair.get("restored_Phi_A") == restored_norm,
        "every_single_restoration_fails_to_repair": all(
            max(map(abs, integer_scores(a_mask ^ (1 << bit)))) >= target + 2
            for bit in range(EDGE_COUNT) if (a_mask ^ c_mask) & (1 << bit)),
        "reported_labelled_switch_classification_matches": repair.get("labelled_switch") == orbit,
    })
    low_sizes = Counter(row[0] for row in low)
    return {"verifier": "plain Python integer and Fraction arithmetic; no GPU table",
            "host": platform.node(), "architecture": platform.machine(),
            "python": platform.python_version(), "A_mask": a_mask, "C_mask": c_mask,
            "D_mask": subset, "minimum_repair_size": minimum_size,
            "minimum_size_repair_count": low_sizes[minimum_size],
            "minimum_repair_tie_break": "ascending D mask after minimum cardinality",
            "low_restoration_count": len(low),
            "low_restoration_size_histogram": {str(key): value for key, value in sorted(low_sizes.items())},
            "restored_A_mask": restored_mask, "restored_Phi_A": restored_norm,
            "restored_A_scores_all_32_states": list(integer_scores(restored_mask)),
            "labelled_switch": orbit,
            "conference_checks": {"A": conference_structure(a_mask),
                                  "restored_A": conference_structure(restored_mask)},
            "original_near_miss": exact,
            "checks": checks, "verified": all(checks.values())}


GPU_SOURCE = r'''
extern "C" __global__ void complete_scores(short* scores) {
    const unsigned index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= 32768u * 32u) return;
    const unsigned signing = index >> 5;
    const unsigned state = index & 31u;
    int total = 0;
    unsigned bit = 0;
    for (unsigned i = 0; i < 6; ++i) {
        const int xi = i == 0 ? 1 : 1 - 2 * (int)((state >> (i - 1)) & 1u);
        for (unsigned j = i + 1; j < 6; ++j, ++bit) {
            const int xj = 1 - 2 * (int)((state >> (j - 1)) & 1u);
            const int edge = 1 - 2 * (int)((signing >> bit) & 1u);
            total += edge * xi * xj;
        }
    }
    scores[index] = (short)total;
}
'''


def gpu_near_miss_repair(record: dict[str, object]) -> dict[str, object]:
    """Authorized new invariant on the FIRST recorded singles-passing reference.

    Reconstruct the fixed n=6 score table once; do not add references or weaken
    the original stronger C-active M-4 filter. This is not a repeated search.
    """
    import cupy as cp

    require(record.get("classification") == "BOUNDED INCONCLUSIVE PROBE"
            and record.get("found") is False and record.get("witness") is None,
            "repair replay requires the original bounded inconclusive record")
    trace = next(row for row in record["reference_trace"]
                 if row["all_single_restoration_candidates"] > 0)
    c_mask = trace["C_mask"]
    started = time.monotonic()
    kernel = cp.RawKernel(GPU_SOURCE, "complete_scores")
    scores = cp.empty((SIGNING_COUNT, STATE_COUNT), dtype=cp.int16)
    kernel(((scores.size + 255) // 256,), (256,), (scores,))
    cp.cuda.runtime.deviceSynchronize()
    table_hash = hashlib.sha256(cp.asnumpy(scores).tobytes()).hexdigest()
    require(table_hash == record["gpu"]["score_table_sha256"],
            "reconstructed score table differs from original search")
    norms = cp.max(cp.abs(scores), axis=1)
    masks = cp.arange(SIGNING_COUNT, dtype=cp.int32)
    c_scores = scores[c_mask]
    c_norm = int(norms[c_mask].item())
    target = c_norm - 4
    norm_ok = norms == target
    a_active = cp.abs(scores) == target
    signed_c_at_a = cp.where(scores > 0, c_scores, -c_scores)
    valley_ok = cp.all(~a_active | (signed_c_at_a < target), axis=1)
    c_active = cp.abs(c_scores) == c_norm
    c_drop_ok = cp.all(~c_active | (scores * cp.sign(c_scores) <= target - 4), axis=1)
    eligible = norm_ok & valley_ok & c_drop_ok
    disagreements = masks ^ c_mask
    single_ok = eligible.copy()
    for bit in range(EDGE_COUNT):
        edge = 1 << bit
        single_ok &= ((disagreements & edge) == 0) | (norms[masks ^ edge] >= target + 2)
    candidate_masks = cp.asnumpy(cp.flatnonzero(single_ok)).tolist()
    observed = {"norm_candidates": int(cp.count_nonzero(norm_ok).item()),
                "active_row_candidates": int(cp.count_nonzero(eligible).item()),
                "all_single_restoration_candidates": len(candidate_masks)}
    require(all(observed[key] == trace[key] for key in observed),
            "first-reference filters do not reproduce original counts")
    a_mask = int(candidate_masks[0])
    h_mask = a_mask ^ c_mask
    restorations = list(restoration_masks(a_mask, c_mask))
    restoration_norms = cp.asnumpy(norms[cp.asarray([mask for _, mask in restorations])]).tolist()
    low = sorted((subset.bit_count(), subset, mask, int(norm))
                 for (subset, mask), norm in zip(restorations, restoration_norms) if norm <= target)
    require(bool(low), "selected near miss unexpectedly passes all-subset gate")
    d_size, d_mask, restored_mask, restored_norm = low[0]
    parameter = rational_valley_lambda(integer_scores(a_mask), integer_scores(c_mask))
    repair = {"classification": "ONE NEAR-MISS MINIMUM INTEGRAL REPAIR; NOT A COUNTERMODEL",
              "n": N, "A_mask": a_mask, "C_mask": c_mask, "H_mask": h_mask,
              "H_size": h_mask.bit_count(), "M": target,
              "A_matrix": signing_matrix(a_mask), "C_matrix": signing_matrix(c_mask),
              "lambda": str(parameter), "D_mask": d_mask, "D_size": d_size,
              "D_edges": [list(pair) for bit, pair in enumerate(EDGES) if d_mask & (1 << bit)],
              "restored_A_mask": restored_mask, "restored_Phi_A": restored_norm,
              "restored_A_matrix": signing_matrix(restored_mask),
              "labelled_switch": labelled_switch_classification(d_mask),
              "selection": {"reference_index": trace["reference_index"],
                            "rule": "first recorded C with singles; first A by ascending mask; minimum (D size,D mask)",
                            "C_active_filter": "sign(C_score)*A_score <= M-4 (stronger than universal M-2 floor)",
                            "original_counts_reproduced": observed,
                            "new_references_examined": 0,
                            "selected_recorded_references_reconstructed": 1},
              "gpu": {"host": platform.node(), "score_table_sha256": table_hash,
                      "score_table_reconstructed_once": True,
                      "restorations_checked": len(restorations),
                      "low_restoration_count": len(low),
                      "elapsed_seconds": round(time.monotonic() - started, 6)},
              "repair_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    exact = verify_repair(repair)
    require(exact["verified"], "GPU repair failed plain-integer and Fraction verification")
    repair["soulkiller_exact_verification"] = exact
    return repair


def gpu_search(seed: int = DEFAULT_SEED, count: int = MAX_REFERENCES) -> dict[str, object]:
    """Run the authorized finite probe; the complete score table is GPU-only."""
    import cupy as cp

    refs = reference_masks(seed, count)
    started = time.monotonic()
    kernel = cp.RawKernel(GPU_SOURCE, "complete_scores")
    scores = cp.empty((SIGNING_COUNT, STATE_COUNT), dtype=cp.int16)
    kernel(((scores.size + 255) // 256,), (256,), (scores,))
    cp.cuda.runtime.deviceSynchronize()
    norms = cp.max(cp.abs(scores), axis=1)
    masks = cp.arange(SIGNING_COUNT, dtype=cp.int32)
    require(bool(cp.all(norms % 2 == 1).item()), "integer norm parity check failed")
    require(int(norms[0].item()) == EDGE_COUNT, "all-positive score normalization failed")
    # This is a spot check, not an independent host recomputation of the table.
    require(cp.asnumpy(scores[0]).tolist() == list(integer_scores(0)),
            "GPU/CPU all-positive score normalization disagrees")
    table_bytes = cp.asnumpy(scores).tobytes()
    traces = []
    witness = None
    for position, c_mask in enumerate(refs):
        c_scores = scores[c_mask]
        c_norm = int(norms[c_mask].item())
        target = c_norm - 4
        norm_ok = norms == target
        a_active = cp.abs(scores) == target
        signed_c_at_a = cp.where(scores > 0, c_scores, -c_scores)
        valley_ok = cp.all(~a_active | (signed_c_at_a < target), axis=1)
        c_active = cp.abs(c_scores) == c_norm
        c_drop_ok = cp.all(~c_active | (scores * cp.sign(c_scores) <= target - 4), axis=1)
        eligible = norm_ok & valley_ok & c_drop_ok
        eligible_count = int(cp.count_nonzero(eligible).item())
        disagreement = masks ^ c_mask
        single_ok = eligible.copy()
        for bit in range(EDGE_COUNT):
            edge = 1 << bit
            single_ok &= ((disagreement & edge) == 0) | (norms[masks ^ edge] >= target + 2)
        singles_count = int(cp.count_nonzero(single_ok).item())
        trace = {"reference_index": position, "C_mask": c_mask, "Phi_C": c_norm,
                 "target_M": target, "norm_candidates": int(cp.count_nonzero(norm_ok).item()),
                 "active_row_candidates": eligible_count,
                 "all_single_restoration_candidates": singles_count,
                 "all_subset_candidates": 0, "all_subset_transform_run": False}
        if singles_count:
            # Integer subset-zeta transform: count ALL X subset H with
            # Phi(C xor X)<=M. Exactly one means A is the unique low vertex
            # in that entire cube; every nonempty restoration then has >=M+2.
            low_subset_count = (norms[masks ^ c_mask] <= target).astype(cp.int32)
            for bit in range(EDGE_COUNT):
                width = 1 << bit
                blocks = low_subset_count.reshape(-1, 2 * width)
                blocks[:, width:] += blocks[:, :width]
            exact_ok = single_ok & (low_subset_count[disagreement] == 1)
            exact_masks = cp.asnumpy(cp.flatnonzero(exact_ok))
            trace["all_subset_transform_run"] = True
            trace["all_subset_candidates"] = int(len(exact_masks))
            if len(exact_masks):
                a_mask = int(exact_masks[0])
                a_cpu = tuple(int(value) for value in cp.asnumpy(scores[a_mask]))
                c_cpu = tuple(int(value) for value in cp.asnumpy(c_scores))
                parameter = rational_valley_lambda(a_cpu, c_cpu)
                h_mask = a_mask ^ c_mask
                witness = {"n": N, "A_mask": a_mask, "C_mask": c_mask,
                           "A_mask_hex": hex(a_mask), "C_mask_hex": hex(c_mask),
                           "H_mask": h_mask, "H_mask_hex": hex(h_mask),
                           "H_edges": [list(pair) for bit, pair in enumerate(EDGES) if h_mask & (1 << bit)],
                           "H_size": h_mask.bit_count(), "M": target,
                           "A_matrix": signing_matrix(a_mask), "C_matrix": signing_matrix(c_mask),
                           "lambda": str(parameter),
                           "gpu_all_subset_low_vertex_count": int(low_subset_count[h_mask].item())}
                exact = verify_witness(witness)
                require(exact["verified"], "GPU candidate failed independent exact verification")
                witness["soulkiller_exact_verification"] = exact
        traces.append(trace)
        print(json.dumps({"progress": trace, "found": witness is not None}), flush=True)
        if witness is not None:
            break
    device = cp.cuda.runtime.getDeviceProperties(0)
    name = device["name"]
    if isinstance(name, bytes):
        name = name.decode()
    return {"probe": "threshold valley versus all-subset integral restoration",
            "classification": "GENERIC SIGNING COUNTERMODEL" if witness else "BOUNDED INCONCLUSIVE PROBE",
            "n": N, "edge_count": EDGE_COUNT, "seed": seed,
            "reference_limit": count, "reference_masks_planned": refs,
            "references_examined": len(traces), "reference_trace": traces,
            "found": witness is not None, "witness": witness,
            "gpu": {"host": platform.node(), "device": name,
                    "compute_capability": [device["major"], device["minor"]],
                    "score_dtype": "int16, integer accumulation",
                    "complete_signing_count": SIGNING_COUNT,
                    "antipodal_states_per_signing": STATE_COUNT,
                    "signed_rows_per_signing": 2 * STATE_COUNT,
                    "score_table_sha256": hashlib.sha256(table_bytes).hexdigest(),
                    "all_subset_method": "exact integer subset-zeta low-vertex count",
                    "elapsed_seconds": round(time.monotonic() - started, 6)},
            "normalization": {"Q": "sum_(i<j) A_ij*x_i*x_j", "Phi": "max_x abs(Q)",
                              "mask_bit_one_means_edge": -1, "edge_order": [list(pair) for pair in EDGES],
                              "state_representative": "x_0=+1; state bits set x_1,...,x_5 to -1",
                              "both_energy_signs_retained": True},
            "independent_node_verifications": [],
            "paley_counterexample_claimed": False, "conference_matrix_claimed": False,
            "residual_ii_closed": False, "limit_closed": False,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    search = sub.add_parser("search")
    search.add_argument("--seed", type=int, default=DEFAULT_SEED)
    search.add_argument("--references", type=int, default=MAX_REFERENCES)
    verify = sub.add_parser("verify")
    verify.add_argument("record", type=Path)
    repair = sub.add_parser("repair")
    repair.add_argument("record", type=Path)
    verify_repair_parser = sub.add_parser("verify-repair")
    verify_repair_parser.add_argument("record", type=Path)
    args = parser.parse_args()
    if args.mode == "search":
        result = gpu_search(args.seed, args.references)
        print("RESULT_JSON=" + json.dumps(result, sort_keys=True), flush=True)
        return 0
    record = json.loads(args.record.read_text())
    if args.mode == "repair":
        result = gpu_near_miss_repair(record)
        print("REPAIR_JSON=" + json.dumps(result, sort_keys=True), flush=True)
        return 0
    if args.mode == "verify-repair":
        result = verify_repair(record.get("near_miss_repair", record))
        print(json.dumps(result, sort_keys=True), flush=True)
        return 0 if result["verified"] else 1
    witness = record.get("witness", record)
    require(isinstance(witness, dict), "record contains no witness to verify")
    result = verify_witness(witness)
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    sys.exit(main())
