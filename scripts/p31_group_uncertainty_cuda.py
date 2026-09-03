#!/usr/bin/env python3
"""V100 heuristic for the open even-support grouped uncertainty case.

For p=31, search even nonzero words f on the 480 antipodal point classes
for ``silent_direction_count > wt(f)``.  Odd weight has a separate symbolic
proof and is intentionally excluded.  Every returned state is reconstructed
on the CPU; this is discovery evidence, never an infeasibility certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_symmetric_halved_mod2 import (  # noqa: E402
    _antipodal_classes,
    _directions,
    _evaluate,
)


P = 31
H = 15
D = 32
N = 480
WORDS = 8


CUDA_SOURCE = r'''
extern "C" {

__device__ __forceinline__ unsigned int xs32(unsigned int &x) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x ? x : 0x9e3779b9U;
}

__global__ void hunt_even(
    const signed char *point_map,
    const unsigned int seed,
    const unsigned int iterations,
    int *margin_out,
    unsigned short *weight_out,
    unsigned short *silent_out,
    unsigned long long *state_out,
    const int chains) {
  const int gid = blockDim.x * blockIdx.x + threadIdx.x;
  if (gid >= chains) return;
  unsigned int rng = seed ^ (0x9e3779b9U * (unsigned int)(gid + 1));
  unsigned short masks[32];
  unsigned long long state[8];
  unsigned long long best_state[8];
  #pragma unroll
  for (int a = 0; a < 32; ++a) masks[a] = 0;
  #pragma unroll
  for (int w = 0; w < 8; ++w) state[w] = 0ULL;

  unsigned short weight = 0;
  const unsigned int mode = ((unsigned int)gid) & 3U;
  if (mode < 2U) {
    int first = (int)(xs32(rng) % 480U);
    int second = (int)(xs32(rng) % 479U);
    if (second >= first) ++second;
    state[first >> 6] ^= 1ULL << (first & 63);
    state[second >> 6] ^= 1ULL << (second & 63);
    weight = 2;
    if (mode == 1U) {
      #pragma unroll
      for (int j = 0; j < 2; ++j) {
        const int point = (int)(xs32(rng) % 480U);
        const unsigned long long bit = 1ULL << (point & 63);
        const bool was_set = (state[point >> 6] & bit) != 0;
        state[point >> 6] ^= bit;
        weight = (unsigned short)(weight + (was_set ? -1 : 1));
      }
    }
  } else {
    const int count = mode == 2U ? 16 : 30;
    for (int j = 0; j < count; ++j) {
      const int point = (int)(xs32(rng) % 480U);
      const unsigned long long bit = 1ULL << (point & 63);
      const bool was_set = (state[point >> 6] & bit) != 0;
      state[point >> 6] ^= bit;
      weight = (unsigned short)(weight + (was_set ? -1 : 1));
    }
  }
  if (weight == 0) {
    state[0] = 3ULL;
    weight = 2;
  }

  for (int point = 0; point < 480; ++point) {
    if ((state[point >> 6] & (1ULL << (point & 63))) == 0) continue;
    #pragma unroll
    for (int a = 0; a < 32; ++a) {
      const int block = (int)point_map[point * 32 + a];
      if (block >= 0) masks[a] ^= (unsigned short)(1U << block);
    }
  }
  unsigned short silent = 0;
  #pragma unroll
  for (int a = 0; a < 32; ++a) silent += masks[a] == 0;
  int best_margin = (int)silent - (int)weight;
  unsigned short best_weight = weight;
  unsigned short best_silent = silent;
  #pragma unroll
  for (int w = 0; w < 8; ++w) best_state[w] = state[w];

  for (unsigned int iteration = 0; iteration < iterations; ++iteration) {
    const int first = (int)(xs32(rng) % 480U);
    int second = (int)(xs32(rng) % 479U);
    if (second >= first) ++second;
    const unsigned long long bit1 = 1ULL << (first & 63);
    const unsigned long long bit2 = 1ULL << (second & 63);
    const bool remove1 = (state[first >> 6] & bit1) != 0;
    const bool remove2 = (state[second >> 6] & bit2) != 0;
    const unsigned short new_weight = (unsigned short)(
        weight + (remove1 ? -1 : 1) + (remove2 ? -1 : 1));
    if (new_weight == 0) continue;

    unsigned short new_silent = silent;
    #pragma unroll
    for (int a = 0; a < 32; ++a) {
      const unsigned short old_mask = masks[a];
      unsigned short new_mask = old_mask;
      const int block1 = (int)point_map[first * 32 + a];
      const int block2 = (int)point_map[second * 32 + a];
      if (block1 >= 0) new_mask ^= (unsigned short)(1U << block1);
      if (block2 >= 0) new_mask ^= (unsigned short)(1U << block2);
      if (old_mask == 0 && new_mask != 0) --new_silent;
      if (old_mask != 0 && new_mask == 0) ++new_silent;
    }
    const int old_score = (int)weight - (int)silent;
    const int new_score = (int)new_weight - (int)new_silent;
    const int uphill = new_score - old_score;
    bool accept = uphill <= 0;
    if (!accept) {
      const unsigned int phase = iteration & 2047U;
      const unsigned int shift = min(15, 3 + uphill + (phase > 1535U ? 0 : 2));
      accept = (xs32(rng) & ((1U << shift) - 1U)) == 0U;
    }
    if (!accept) continue;

    state[first >> 6] ^= bit1;
    state[second >> 6] ^= bit2;
    #pragma unroll
    for (int a = 0; a < 32; ++a) {
      const int block1 = (int)point_map[first * 32 + a];
      const int block2 = (int)point_map[second * 32 + a];
      if (block1 >= 0) masks[a] ^= (unsigned short)(1U << block1);
      if (block2 >= 0) masks[a] ^= (unsigned short)(1U << block2);
    }
    weight = new_weight;
    silent = new_silent;
    const int margin = (int)silent - (int)weight;
    if (margin > best_margin ||
        (margin == best_margin && weight < best_weight)) {
      best_margin = margin;
      best_weight = weight;
      best_silent = silent;
      #pragma unroll
      for (int w = 0; w < 8; ++w) best_state[w] = state[w];
    }
  }
  margin_out[gid] = best_margin;
  weight_out[gid] = best_weight;
  silent_out[gid] = best_silent;
  #pragma unroll
  for (int w = 0; w < 8; ++w) state_out[(long long)gid * 8 + w] = best_state[w];
}
}
'''


def point_map() -> np.ndarray:
    points = _antipodal_classes(P)
    squares = sorted({value * value % P for value in range(1, P)})
    square_index = {square: index for index, square in enumerate(squares)}
    out = np.full((N, D), -1, dtype=np.int8)
    for point_index, point in enumerate(points):
        for direction, functional in enumerate(_directions(P)):
            value = _evaluate(P, functional, point)
            if value:
                out[point_index, direction] = square_index[value * value % P]
    return out


def verify(state: np.ndarray, mapping: np.ndarray) -> dict[str, object]:
    support = tuple(
        point
        for point in range(N)
        if int(state[point >> 6]) & (1 << (point & 63))
    )
    masks = [0] * D
    for point in support:
        for direction in range(D):
            block = int(mapping[point, direction])
            if block >= 0:
                masks[direction] ^= 1 << block
    silent = tuple(i for i, mask in enumerate(masks) if mask == 0)
    margin = len(silent) - len(support)
    return {
        "support_indices": list(support),
        "silent_direction_indices": list(silent),
        "weight": len(support),
        "silent_group_count": len(silent),
        "best_margin": margin,
        "counterexample": margin > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chains", type=int, default=32768)
    parser.add_argument("--iterations", type=int, default=60000)
    parser.add_argument("--seed", type=int, default=310035)
    args = parser.parse_args()
    if args.chains <= 0 or args.iterations <= 0:
        parser.error("chains and iterations must be positive")

    import cupy as cp

    mapping = point_map()
    device_map = cp.asarray(mapping)
    margins = cp.empty(args.chains, dtype=cp.int32)
    weights = cp.empty(args.chains, dtype=cp.uint16)
    silents = cp.empty(args.chains, dtype=cp.uint16)
    states = cp.empty((args.chains, WORDS), dtype=cp.uint64)
    kernel = cp.RawKernel(CUDA_SOURCE, "hunt_even", options=("-std=c++11",))
    threads = 256
    blocks = (args.chains + threads - 1) // threads
    kernel(
        (blocks,),
        (threads,),
        (
            device_map,
            np.uint32(args.seed),
            np.uint32(args.iterations),
            margins,
            weights,
            silents,
            states,
            np.int32(args.chains),
        ),
    )
    cp.cuda.Stream.null.synchronize()
    host_margins = cp.asnumpy(margins)
    best = int(np.argmax(host_margins))
    device_weight = int(cp.asnumpy(weights[best]))
    device_silent = int(cp.asnumpy(silents[best]))
    state = cp.asnumpy(states[best])
    checked = verify(state, mapping)
    verified = (
        checked["best_margin"] == int(host_margins[best])
        and checked["weight"] == device_weight
        and checked["silent_group_count"] == device_silent
        and checked["weight"] % 2 == 0
        and checked["weight"] > 0
    )
    out = {
        "p": P,
        "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "backend": f"CuPy {cp.__version__}",
        "chains": args.chains,
        "iterations": args.iterations,
        "seed": args.seed,
        "scope": "even support only",
        **checked,
        "device_reported_weight": device_weight,
        "device_reported_silent": device_silent,
        "verified": verified,
        "status": "heuristic search only; absence of a counterexample is not a proof",
    }
    print(json.dumps(out, indent=2))
    raise SystemExit(0 if verified else 2)


if __name__ == "__main__":
    main()
