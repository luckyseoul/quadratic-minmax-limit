#!/usr/bin/env python3
"""V100 meet-in-the-middle sieve for one exceptional p=7 mean leaf.

The sieve uses only explicitly selected rows of the exact mod-3 and mod-7
left-dependency systems.  A zero projected match is therefore a rigorous
infeasibility certificate; a nonzero match is reported as unresolved.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import time
from pathlib import Path

import cupy as cp
import numpy as np


CUDA_SOURCE = r"""
extern "C" {

__device__ __forceinline__ unsigned long long mix64(unsigned long long x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

__device__ __forceinline__ unsigned long long sum_key(
    unsigned int a3, unsigned long long a7,
    unsigned int b3, unsigned long long b7,
    unsigned int c3, unsigned long long c7
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 12; ++i) {
        unsigned int digit = ((a3 >> (2*i)) & 3U)
                           + ((b3 >> (2*i)) & 3U)
                           + ((c3 >> (2*i)) & 3U);
        digit %= 3U;
        key += place * (unsigned long long)digit;
        place *= 3ULL;
    }
    #pragma unroll
    for (int i = 0; i < 13; ++i) {
        unsigned int digit = (unsigned int)((a7 >> (4*i)) & 15ULL)
                           + (unsigned int)((b7 >> (4*i)) & 15ULL)
                           + (unsigned int)((c7 >> (4*i)) & 15ULL);
        digit %= 7U;
        key += place * (unsigned long long)digit;
        place *= 7ULL;
    }
    return key;
}

__device__ __forceinline__ unsigned long long needed_key(
    unsigned int target3, unsigned long long target7,
    unsigned int a3, unsigned long long a7,
    unsigned int b3, unsigned long long b7,
    unsigned int c3, unsigned long long c7
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 12; ++i) {
        unsigned int total = ((a3 >> (2*i)) & 3U)
                           + ((b3 >> (2*i)) & 3U)
                           + ((c3 >> (2*i)) & 3U);
        unsigned int target = (target3 >> (2*i)) & 3U;
        unsigned int digit = (target + 6U - total) % 3U;
        key += place * (unsigned long long)digit;
        place *= 3ULL;
    }
    #pragma unroll
    for (int i = 0; i < 13; ++i) {
        unsigned int total = (unsigned int)((a7 >> (4*i)) & 15ULL)
                           + (unsigned int)((b7 >> (4*i)) & 15ULL)
                           + (unsigned int)((c7 >> (4*i)) & 15ULL);
        unsigned int target = (unsigned int)((target7 >> (4*i)) & 15ULL);
        unsigned int digit = (target + 21U - total) % 7U;
        key += place * (unsigned long long)digit;
        place *= 7ULL;
    }
    return key;
}

__global__ void make_sum_keys3(
    const unsigned int* a3, const unsigned long long* a7, unsigned long long na,
    const unsigned int* b3, const unsigned long long* b7, unsigned long long nb,
    const unsigned int* c3, const unsigned long long* c7, unsigned long long nc,
    unsigned long long offset, unsigned long long count,
    unsigned long long* output
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    output[local] = sum_key(a3[ia], a7[ia], b3[ib], b7[ib], c3[ic], c7[ic]);
}

__global__ void bloom_insert(
    const unsigned long long* keys, unsigned long long count,
    unsigned int* bloom, unsigned long long bit_mask
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    unsigned long long key = keys[index];
    unsigned long long h1 = mix64(key) & bit_mask;
    unsigned long long h2 = mix64(key ^ 0x9e3779b97f4a7c15ULL) & bit_mask;
    atomicOr(&bloom[h1 >> 5], 1U << (h1 & 31ULL));
    atomicOr(&bloom[h2 >> 5], 1U << (h2 & 31ULL));
}

__global__ void probe_needed3(
    const unsigned int* a3, const unsigned long long* a7, unsigned long long na,
    const unsigned int* b3, const unsigned long long* b7, unsigned long long nb,
    const unsigned int* c3, const unsigned long long* c7, unsigned long long nc,
    unsigned int target3, unsigned long long target7,
    unsigned long long offset, unsigned long long count,
    const unsigned int* bloom, unsigned long long bit_mask,
    unsigned long long* candidates, unsigned long long capacity,
    unsigned long long* candidate_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long key = needed_key(
        target3, target7, a3[ia], a7[ia], b3[ib], b7[ib], c3[ic], c7[ic]
    );
    unsigned long long h1 = mix64(key) & bit_mask;
    unsigned long long h2 = mix64(key ^ 0x9e3779b97f4a7c15ULL) & bit_mask;
    unsigned int present1 = (bloom[h1 >> 5] >> (h1 & 31ULL)) & 1U;
    unsigned int present2 = (bloom[h2 >> 5] >> (h2 & 31ULL)) & 1U;
    if (present1 && present2) {
        unsigned long long slot = atomicAdd(candidate_count, 1ULL);
        if (slot < capacity) candidates[slot] = key;
    }
}

__global__ void count_sorted_matches(
    const unsigned long long* sorted_keys, unsigned long long key_count,
    const unsigned long long* probes, unsigned long long probe_count,
    unsigned long long* match_count
) {
    unsigned long long index = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= probe_count) return;
    unsigned long long key = probes[index];
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        if (sorted_keys[middle] < key) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_keys[low] == key) {
        atomicAdd(match_count, 1ULL);
    }
}

__global__ void count_needed_sorted_matches3(
    const unsigned int* a3, const unsigned long long* a7, unsigned long long na,
    const unsigned int* b3, const unsigned long long* b7, unsigned long long nb,
    const unsigned int* c3, const unsigned long long* c7, unsigned long long nc,
    unsigned int target3, unsigned long long target7,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_keys, unsigned long long key_count,
    unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long key = needed_key(
        target3, target7, a3[ia], a7[ia], b3[ib], b7[ib], c3[ic], c7[ic]
    );
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        if (sorted_keys[middle] < key) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_keys[low] == key) {
        atomicAdd(match_count, 1ULL);
    }
}

__device__ __forceinline__ unsigned long long sum_key7_22(
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 11; ++i) {
        unsigned int digit = (unsigned int)((alo >> (4*i)) & 15ULL)
                           + (unsigned int)((blo >> (4*i)) & 15ULL)
                           + (unsigned int)((clo >> (4*i)) & 15ULL);
        key += place * (unsigned long long)(digit % 7U);
        place *= 7ULL;
    }
    #pragma unroll
    for (int i = 0; i < 11; ++i) {
        unsigned int digit = (unsigned int)((ahi >> (4*i)) & 15ULL)
                           + (unsigned int)((bhi >> (4*i)) & 15ULL)
                           + (unsigned int)((chi >> (4*i)) & 15ULL);
        key += place * (unsigned long long)(digit % 7U);
        place *= 7ULL;
    }
    return key;
}

__device__ __forceinline__ unsigned long long needed_key7_22(
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 11; ++i) {
        unsigned int total = (unsigned int)((alo >> (4*i)) & 15ULL)
                           + (unsigned int)((blo >> (4*i)) & 15ULL)
                           + (unsigned int)((clo >> (4*i)) & 15ULL);
        unsigned int target = (unsigned int)((target_lo >> (4*i)) & 15ULL);
        key += place * (unsigned long long)((target + 21U - total) % 7U);
        place *= 7ULL;
    }
    #pragma unroll
    for (int i = 0; i < 11; ++i) {
        unsigned int total = (unsigned int)((ahi >> (4*i)) & 15ULL)
                           + (unsigned int)((bhi >> (4*i)) & 15ULL)
                           + (unsigned int)((chi >> (4*i)) & 15ULL);
        unsigned int target = (unsigned int)((target_hi >> (4*i)) & 15ULL);
        key += place * (unsigned long long)((target + 21U - total) % 7U);
        place *= 7ULL;
    }
    return key;
}

__global__ void make_sum_keys7_22(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long offset, unsigned long long count,
    unsigned long long* output
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    output[local] = sum_key7_22(
        alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]
    );
}

__global__ void count_needed_sorted_matches7_22(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_keys, unsigned long long key_count,
    unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long key = needed_key7_22(
        target_lo, target_hi,
        alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]
    );
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        if (sorted_keys[middle] < key) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_keys[low] == key) {
        atomicAdd(match_count, 1ULL);
    }
}

__device__ __forceinline__ unsigned long long sum_key3_40(
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 20; ++i) {
        unsigned int digit = (unsigned int)((alo >> (2*i)) & 3ULL)
                           + (unsigned int)((blo >> (2*i)) & 3ULL)
                           + (unsigned int)((clo >> (2*i)) & 3ULL);
        key += place * (unsigned long long)(digit % 3U);
        place *= 3ULL;
    }
    #pragma unroll
    for (int i = 0; i < 20; ++i) {
        unsigned int digit = (unsigned int)((ahi >> (2*i)) & 3ULL)
                           + (unsigned int)((bhi >> (2*i)) & 3ULL)
                           + (unsigned int)((chi >> (2*i)) & 3ULL);
        key += place * (unsigned long long)(digit % 3U);
        place *= 3ULL;
    }
    return key;
}

__device__ __forceinline__ unsigned long long needed_key3_40(
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 20; ++i) {
        unsigned int total = (unsigned int)((alo >> (2*i)) & 3ULL)
                           + (unsigned int)((blo >> (2*i)) & 3ULL)
                           + (unsigned int)((clo >> (2*i)) & 3ULL);
        unsigned int target = (unsigned int)((target_lo >> (2*i)) & 3ULL);
        key += place * (unsigned long long)((target + 6U - total) % 3U);
        place *= 3ULL;
    }
    #pragma unroll
    for (int i = 0; i < 20; ++i) {
        unsigned int total = (unsigned int)((ahi >> (2*i)) & 3ULL)
                           + (unsigned int)((bhi >> (2*i)) & 3ULL)
                           + (unsigned int)((chi >> (2*i)) & 3ULL);
        unsigned int target = (unsigned int)((target_hi >> (2*i)) & 3ULL);
        key += place * (unsigned long long)((target + 6U - total) % 3U);
        place *= 3ULL;
    }
    return key;
}

__global__ void make_sum_keys3_40(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long offset, unsigned long long count,
    unsigned long long* output
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    output[local] = sum_key3_40(alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]);
}

__global__ void count_needed_sorted_matches3_40(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_keys, unsigned long long key_count,
    unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long key = needed_key3_40(
        target_lo, target_hi, alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]
    );
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        if (sorted_keys[middle] < key) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_keys[low] == key) atomicAdd(match_count, 1ULL);
}

__device__ __forceinline__ unsigned long long sum_key5_27(
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 13; ++i) {
        unsigned int digit = (unsigned int)((alo >> (3*i)) & 7ULL)
                           + (unsigned int)((blo >> (3*i)) & 7ULL)
                           + (unsigned int)((clo >> (3*i)) & 7ULL);
        key += place * (unsigned long long)(digit % 5U);
        place *= 5ULL;
    }
    #pragma unroll
    for (int i = 0; i < 14; ++i) {
        unsigned int digit = (unsigned int)((ahi >> (3*i)) & 7ULL)
                           + (unsigned int)((bhi >> (3*i)) & 7ULL)
                           + (unsigned int)((chi >> (3*i)) & 7ULL);
        key += place * (unsigned long long)(digit % 5U);
        place *= 5ULL;
    }
    return key;
}

__device__ __forceinline__ unsigned long long needed_key5_27(
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long alo, unsigned long long ahi,
    unsigned long long blo, unsigned long long bhi,
    unsigned long long clo, unsigned long long chi
) {
    unsigned long long key = 0ULL;
    unsigned long long place = 1ULL;
    #pragma unroll
    for (int i = 0; i < 13; ++i) {
        unsigned int total = (unsigned int)((alo >> (3*i)) & 7ULL)
                           + (unsigned int)((blo >> (3*i)) & 7ULL)
                           + (unsigned int)((clo >> (3*i)) & 7ULL);
        unsigned int target = (unsigned int)((target_lo >> (3*i)) & 7ULL);
        key += place * (unsigned long long)((target + 15U - total) % 5U);
        place *= 5ULL;
    }
    #pragma unroll
    for (int i = 0; i < 14; ++i) {
        unsigned int total = (unsigned int)((ahi >> (3*i)) & 7ULL)
                           + (unsigned int)((bhi >> (3*i)) & 7ULL)
                           + (unsigned int)((chi >> (3*i)) & 7ULL);
        unsigned int target = (unsigned int)((target_hi >> (3*i)) & 7ULL);
        key += place * (unsigned long long)((target + 15U - total) % 5U);
        place *= 5ULL;
    }
    return key;
}

__global__ void make_sum_keys5_27(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long offset, unsigned long long count, unsigned long long* output
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    output[local] = sum_key5_27(alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]);
}

__global__ void count_needed_sorted_matches5_27(
    const unsigned long long* alo, const unsigned long long* ahi, unsigned long long na,
    const unsigned long long* blo, const unsigned long long* bhi, unsigned long long nb,
    const unsigned long long* clo, const unsigned long long* chi, unsigned long long nc,
    unsigned long long target_lo, unsigned long long target_hi,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_keys, unsigned long long key_count,
    unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long key = needed_key5_27(
        target_lo, target_hi, alo[ia], ahi[ia], blo[ib], bhi[ib], clo[ic], chi[ic]
    );
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        if (sorted_keys[middle] < key) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_keys[low] == key) atomicAdd(match_count, 1ULL);
}

__global__ void count_needed_sorted_pairs7_22(
    const unsigned long long* a0lo, const unsigned long long* a0hi,
    const unsigned long long* b0lo, const unsigned long long* b0hi, unsigned long long na,
    const unsigned long long* a1lo, const unsigned long long* a1hi,
    const unsigned long long* b1lo, const unsigned long long* b1hi, unsigned long long nb,
    const unsigned long long* a2lo, const unsigned long long* a2hi,
    const unsigned long long* b2lo, const unsigned long long* b2hi, unsigned long long nc,
    unsigned long long target_a_lo, unsigned long long target_a_hi,
    unsigned long long target_b_lo, unsigned long long target_b_hi,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_a, const unsigned long long* sorted_b,
    unsigned long long key_count, unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long ka = needed_key7_22(
        target_a_lo, target_a_hi,
        a0lo[ia], a0hi[ia], a1lo[ib], a1hi[ib], a2lo[ic], a2hi[ic]
    );
    unsigned long long kb = needed_key7_22(
        target_b_lo, target_b_hi,
        b0lo[ia], b0hi[ia], b1lo[ib], b1hi[ib], b2lo[ic], b2hi[ic]
    );
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        unsigned long long ma = sorted_a[middle];
        unsigned long long mb = sorted_b[middle];
        if (ma < ka || (ma == ka && mb < kb)) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_a[low] == ka && sorted_b[low] == kb) {
        atomicAdd(match_count, 1ULL);
    }
}

__global__ void count_needed_sorted_triples7_22(
    const unsigned long long* a0lo, const unsigned long long* a0hi,
    const unsigned long long* b0lo, const unsigned long long* b0hi,
    const unsigned long long* c0lo, const unsigned long long* c0hi, unsigned long long na,
    const unsigned long long* a1lo, const unsigned long long* a1hi,
    const unsigned long long* b1lo, const unsigned long long* b1hi,
    const unsigned long long* c1lo, const unsigned long long* c1hi, unsigned long long nb,
    const unsigned long long* a2lo, const unsigned long long* a2hi,
    const unsigned long long* b2lo, const unsigned long long* b2hi,
    const unsigned long long* c2lo, const unsigned long long* c2hi, unsigned long long nc,
    unsigned long long target_a_lo, unsigned long long target_a_hi,
    unsigned long long target_b_lo, unsigned long long target_b_hi,
    unsigned long long target_c_lo, unsigned long long target_c_hi,
    unsigned long long offset, unsigned long long count,
    const unsigned long long* sorted_a, const unsigned long long* sorted_b,
    const unsigned long long* sorted_c, unsigned long long key_count,
    unsigned long long* match_count
) {
    unsigned long long local = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (local >= count) return;
    unsigned long long index = offset + local;
    unsigned long long ic = index % nc; index /= nc;
    unsigned long long ib = index % nb; index /= nb;
    unsigned long long ia = index;
    unsigned long long ka = needed_key7_22(target_a_lo, target_a_hi,
        a0lo[ia], a0hi[ia], a1lo[ib], a1hi[ib], a2lo[ic], a2hi[ic]);
    unsigned long long kb = needed_key7_22(target_b_lo, target_b_hi,
        b0lo[ia], b0hi[ia], b1lo[ib], b1hi[ib], b2lo[ic], b2hi[ic]);
    unsigned long long kc = needed_key7_22(target_c_lo, target_c_hi,
        c0lo[ia], c0hi[ia], c1lo[ib], c1hi[ib], c2lo[ic], c2hi[ic]);
    unsigned long long low = 0ULL;
    unsigned long long high = key_count;
    while (low < high) {
        unsigned long long middle = low + ((high - low) >> 1);
        unsigned long long ma = sorted_a[middle];
        unsigned long long mb = sorted_b[middle];
        unsigned long long mc = sorted_c[middle];
        if (ma < ka || (ma == ka && (mb < kb || (mb == kb && mc < kc)))) low = middle + 1ULL;
        else high = middle;
    }
    if (low < key_count && sorted_a[low] == ka && sorted_b[low] == kb && sorted_c[low] == kc) {
        atomicAdd(match_count, 1ULL);
    }
}

}
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def add_packed(values: list[int], modulus: int, bits: int, count: int) -> int:
    out = 0
    for index in range(count):
        digit = sum((value >> (bits * index)) & ((1 << bits) - 1) for value in values)
        out |= (digit % modulus) << (bits * index)
    return out


def negate_packed(value: int, modulus: int, bits: int, count: int) -> int:
    out = 0
    mask = (1 << bits) - 1
    for index in range(count):
        digit = (-(value >> (bits * index) & mask)) % modulus
        out |= digit << (bits * index)
    return out


def choose_partition(sizes: tuple[int, ...], max_build: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not sizes:
        return (), ()
    if len(sizes) == 1:
        if sizes[0] > max_build:
            raise ValueError("single catalog exceeds the GPU build cap")
        return (0,), ()
    best = None
    indices = tuple(range(len(sizes)))
    for count in range(1, min(3, len(sizes) - 1) + 1):
        for first in itertools.combinations(indices, count):
            second = tuple(index for index in indices if index not in first)
            if len(second) > 3:
                continue
            build = math.prod(sizes[index] for index in first)
            probe = math.prod(sizes[index] for index in second)
            if build > max_build:
                continue
            score = (probe, build)
            if best is None or score < best[0]:
                best = (score, first, second)
    if best is None:
        raise ValueError("no three-by-three partition fits the GPU build cap")
    return best[1], best[2]


def pad_group(group: list[tuple[cp.ndarray, cp.ndarray]]) -> list[tuple[cp.ndarray, cp.ndarray]]:
    one = (cp.asarray([0], dtype=cp.uint32), cp.asarray([0], dtype=cp.uint64))
    return group + [one] * (3 - len(group))


def run(
    projection_summary_path: Path,
    projection_path: Path,
    mean_path: Path,
    leaf_index: int,
    max_build: int,
    chunk_size: int,
    candidate_capacity: int,
) -> dict:
    started = time.time()
    projection_summary = json.loads(projection_summary_path.read_text())
    if (
        projection_summary.get("experiment") != "p7_exceptional_projected_catalogs"
        or projection_summary.get("projection_group_order") != 3**12 * 7**13
    ):
        raise ValueError("unexpected projected-catalog summary")
    mean_batch = json.loads(mean_path.read_text())
    leaf = mean_batch["leaves"][leaf_index]
    if int(leaf["leaf_index"]) != leaf_index:
        raise AssertionError("mean leaf index mismatch")
    if leaf.get("solver_status") == "INFEASIBLE":
        raise ValueError("the requested leaf is already excluded by CP-SAT")
    means = tuple(int(value) for value in leaf["scaled_means_direction_order"])

    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++14",))
    make_sum = module.get_function("make_sum_keys3")
    bloom_insert = module.get_function("bloom_insert")
    probe_needed = module.get_function("probe_needed3")
    count_matches = module.get_function("count_sorted_matches")
    count_needed_matches = module.get_function("count_needed_sorted_matches3")

    with np.load(projection_path, allow_pickle=False) as source:
        base3 = int(source["base_p3"][0])
        base7 = int(source["base_p7"][0])
        variable = []
        singleton3 = [base3]
        singleton7 = [base7]
        metadata = []
        for direction_index, mean in enumerate(means):
            host3 = source[f"d{direction_index}_m{mean}_p3"]
            host7 = source[f"d{direction_index}_m{mean}_p7"]
            metadata.append(
                {
                    "direction_index": direction_index,
                    "scaled_mean": mean,
                    "catalog_rows": int(len(host3)),
                }
            )
            if len(host3) == 1:
                singleton3.append(int(host3[0]))
                singleton7.append(int(host7[0]))
            else:
                variable.append(
                    (
                        cp.asarray(host3, dtype=cp.uint32),
                        cp.asarray(host7, dtype=cp.uint64),
                    )
                )
    target3 = negate_packed(add_packed(singleton3, 3, 2, 12), 3, 2, 12)
    target7 = negate_packed(add_packed(singleton7, 7, 4, 13), 7, 4, 13)
    sizes = tuple(len(item[0]) for item in variable)
    first_indices, second_indices = choose_partition(sizes, max_build)
    first = pad_group([variable[index] for index in first_indices])
    second = pad_group([variable[index] for index in second_indices])
    build_count = math.prod(len(item[0]) for item in first)
    probe_count = math.prod(len(item[0]) for item in second)

    block = 256
    build_keys = cp.empty(build_count, dtype=cp.uint64)
    make_sum(
        ((build_count + block - 1) // block,),
        (block,),
        (
            first[0][0], first[0][1], np.uint64(len(first[0][0])),
            first[1][0], first[1][1], np.uint64(len(first[1][0])),
            first[2][0], first[2][1], np.uint64(len(first[2][0])),
            np.uint64(0), np.uint64(build_count), build_keys,
        ),
    )
    cp.cuda.get_current_stream().synchronize()

    desired_bits = max(30, math.ceil(math.log2(max(1, build_count * 32))))
    bloom_bits_power = min(33, desired_bits)
    bloom_bits = 1 << bloom_bits_power
    bloom = cp.zeros(bloom_bits // 32, dtype=cp.uint32)
    bloom_insert(
        ((build_count + block - 1) // block,),
        (block,),
        (build_keys, np.uint64(build_count), bloom, np.uint64(bloom_bits - 1)),
    )
    candidates = cp.empty(candidate_capacity, dtype=cp.uint64)
    candidate_count = cp.zeros(1, dtype=cp.uint64)
    for offset in range(0, probe_count, chunk_size):
        count = min(chunk_size, probe_count - offset)
        probe_needed(
            ((count + block - 1) // block,),
            (block,),
            (
                second[0][0], second[0][1], np.uint64(len(second[0][0])),
                second[1][0], second[1][1], np.uint64(len(second[1][0])),
                second[2][0], second[2][1], np.uint64(len(second[2][0])),
                np.uint32(target3), np.uint64(target7),
                np.uint64(offset), np.uint64(count),
                bloom, np.uint64(bloom_bits - 1), candidates,
                np.uint64(candidate_capacity), candidate_count,
            ),
        )
    cp.cuda.get_current_stream().synchronize()
    bloom_candidates = int(candidate_count.get()[0])
    exact_matches = 0
    exact_join_method = "bloom_candidates_binary_search"
    if bloom_candidates > candidate_capacity:
        # The Bloom filter is only an acceleration heuristic.  If its candidate
        # buffer fills, retain exactness by binary-searching every probe key.
        exact_join_method = "full_probe_binary_search_after_bloom_overflow"
        build_keys.sort()
        match_count = cp.zeros(1, dtype=cp.uint64)
        for offset in range(0, probe_count, chunk_size):
            count = min(chunk_size, probe_count - offset)
            count_needed_matches(
                ((count + block - 1) // block,),
                (block,),
                (
                    second[0][0], second[0][1], np.uint64(len(second[0][0])),
                    second[1][0], second[1][1], np.uint64(len(second[1][0])),
                    second[2][0], second[2][1], np.uint64(len(second[2][0])),
                    np.uint32(target3), np.uint64(target7),
                    np.uint64(offset), np.uint64(count),
                    build_keys, np.uint64(build_count), match_count,
                ),
            )
        exact_matches = int(match_count.get()[0])
    elif bloom_candidates:
        build_keys.sort()
        probe = candidates[:bloom_candidates]
        match_count = cp.zeros(1, dtype=cp.uint64)
        count_matches(
            ((bloom_candidates + block - 1) // block,),
            (block,),
            (
                build_keys,
                np.uint64(build_count),
                probe,
                np.uint64(bloom_candidates),
                match_count,
            ),
        )
        exact_matches = int(match_count.get()[0])

    return {
        "experiment": "p7_exceptional_projected_join_gpu",
        "status": "complete_exact_selected_dependency_gpu_join",
        "p": 7,
        "c_H": -1,
        "fixed_boundary": mean_batch["fixed_boundary"],
        "leaf_index": leaf_index,
        "fixed_scaled_means": list(means),
        "projection_summary": str(projection_summary_path),
        "projection_summary_sha256": sha256(projection_summary_path),
        "projection_catalog": str(projection_path),
        "projection_catalog_sha256": sha256(projection_path),
        "projection_group_order": 3**12 * 7**13,
        "catalogs": metadata,
        "variable_catalog_sizes": list(sizes),
        "partition": [list(first_indices), list(second_indices)],
        "build_signatures": build_count,
        "probe_signatures": probe_count,
        "bloom_bits": bloom_bits,
        "bloom_candidates": bloom_candidates,
        "bloom_candidate_capacity": candidate_capacity,
        "exact_join_method": exact_join_method,
        "exact_projected_matches": exact_matches,
        "projected_modularly_infeasible": exact_matches == 0,
        "finite_mean_allocation_exclusion": exact_matches == 0,
        "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "elapsed_seconds": time.time() - started,
    }


def run_tail22(
    projection_summary_path: Path,
    projection_path: Path,
    mean_path: Path,
    leaf_index: int,
    max_build: int,
    chunk_size: int,
) -> dict:
    """Run an injective single-prime projection without a Bloom filter."""
    started = time.time()
    summary = json.loads(projection_summary_path.read_text())
    summary_kind = summary.get("experiment")
    if summary_kind == "p7_exceptional_tail22_catalogs":
        modulus, bits, lo_count, hi_count = 7, 4, 11, 11
        group_order = 7**22
        key_stem = "p7"
        make_kernel = "make_sum_keys7_22"
        count_kernel = "count_needed_sorted_matches7_22"
        projection_mode = "injective_tail22_mod7"
    elif summary_kind == "p7_exceptional_p3_40_catalogs":
        modulus, bits, lo_count, hi_count = 3, 2, 20, 20
        group_order = 3**40
        key_stem = "p3"
        make_kernel = "make_sum_keys3_40"
        count_kernel = "count_needed_sorted_matches3_40"
        projection_mode = "injective_directed40_mod3"
    elif summary_kind == "p7_exceptional_p5_27_catalogs":
        modulus, bits, lo_count, hi_count = 5, 3, 13, 14
        group_order = 5**27
        key_stem = "p5"
        make_kernel = "make_sum_keys5_27"
        count_kernel = "count_needed_sorted_matches5_27"
        projection_mode = "injective_balanced27_mod5"
    else:
        raise ValueError("unexpected injective projected-catalog summary")
    if summary.get("projection_group_order") != group_order:
        raise ValueError("injective projected-catalog group order changed")
    mean_batch = json.loads(mean_path.read_text())
    leaf = mean_batch["leaves"][leaf_index]
    if int(leaf["leaf_index"]) != leaf_index:
        raise AssertionError("mean leaf index mismatch")
    if leaf.get("solver_status") == "INFEASIBLE":
        raise ValueError("the requested leaf is already excluded by CP-SAT")
    means = tuple(int(value) for value in leaf["scaled_means_direction_order"])

    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++14",))
    make_sum = module.get_function(make_kernel)
    count_needed = module.get_function(count_kernel)

    with np.load(projection_path, allow_pickle=False) as source:
        singleton_lo = [int(source[f"base_{key_stem}_lo"][0])]
        singleton_hi = [int(source[f"base_{key_stem}_hi"][0])]
        variable = []
        metadata = []
        for direction_index, mean in enumerate(means):
            host_lo = source[f"d{direction_index}_m{mean}_{key_stem}_lo"]
            host_hi = source[f"d{direction_index}_m{mean}_{key_stem}_hi"]
            metadata.append(
                {
                    "direction_index": direction_index,
                    "scaled_mean": mean,
                    "catalog_rows": int(len(host_lo)),
                }
            )
            if len(host_lo) == 1:
                singleton_lo.append(int(host_lo[0]))
                singleton_hi.append(int(host_hi[0]))
            else:
                variable.append(
                    (
                        cp.asarray(host_lo, dtype=cp.uint64),
                        cp.asarray(host_hi, dtype=cp.uint64),
                    )
                )
    target_lo = negate_packed(
        add_packed(singleton_lo, modulus, bits, lo_count),
        modulus, bits, lo_count,
    )
    target_hi = negate_packed(
        add_packed(singleton_hi, modulus, bits, hi_count),
        modulus, bits, hi_count,
    )
    sizes = tuple(len(item[0]) for item in variable)
    first_indices, second_indices = choose_partition(sizes, max_build)
    zero = (cp.asarray([0], dtype=cp.uint64), cp.asarray([0], dtype=cp.uint64))
    first = [variable[index] for index in first_indices]
    second = [variable[index] for index in second_indices]
    first += [zero] * (3 - len(first))
    second += [zero] * (3 - len(second))
    build_count = math.prod(len(item[0]) for item in first)
    probe_count = math.prod(len(item[0]) for item in second)

    block = 256
    build_keys = cp.empty(build_count, dtype=cp.uint64)
    make_sum(
        ((build_count + block - 1) // block,),
        (block,),
        (
            first[0][0], first[0][1], np.uint64(len(first[0][0])),
            first[1][0], first[1][1], np.uint64(len(first[1][0])),
            first[2][0], first[2][1], np.uint64(len(first[2][0])),
            np.uint64(0), np.uint64(build_count), build_keys,
        ),
    )
    build_keys.sort()
    match_count = cp.zeros(1, dtype=cp.uint64)
    for offset in range(0, probe_count, chunk_size):
        count = min(chunk_size, probe_count - offset)
        count_needed(
            ((count + block - 1) // block,),
            (block,),
            (
                second[0][0], second[0][1], np.uint64(len(second[0][0])),
                second[1][0], second[1][1], np.uint64(len(second[1][0])),
                second[2][0], second[2][1], np.uint64(len(second[2][0])),
                np.uint64(target_lo), np.uint64(target_hi),
                np.uint64(offset), np.uint64(count),
                build_keys, np.uint64(build_count), match_count,
            ),
        )
    exact_matches = int(match_count.get()[0])
    return {
        "experiment": "p7_exceptional_projected_join_gpu",
        "status": "complete_exact_selected_dependency_gpu_join",
        "projection_mode": projection_mode,
        "p": 7,
        "c_H": -1,
        "fixed_boundary": mean_batch["fixed_boundary"],
        "leaf_index": leaf_index,
        "fixed_scaled_means": list(means),
        "projection_summary": str(projection_summary_path),
        "projection_summary_sha256": sha256(projection_summary_path),
        "projection_catalog": str(projection_path),
        "projection_catalog_sha256": sha256(projection_path),
        "projection_group_order": group_order,
        "catalogs": metadata,
        "variable_catalog_sizes": list(sizes),
        "partition": [list(first_indices), list(second_indices)],
        "build_signatures": build_count,
        "probe_signatures": probe_count,
        "exact_join_method": "full_probe_binary_search",
        "exact_projected_matches": exact_matches,
        "projected_modularly_infeasible": exact_matches == 0,
        "finite_mean_allocation_exclusion": exact_matches == 0,
        "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "elapsed_seconds": time.time() - started,
    }


def run_mod7_pair(
    summary_a_path: Path,
    catalog_a_path: Path,
    summary_b_path: Path,
    catalog_b_path: Path,
    mean_path: Path,
    leaf_index: int,
    max_build: int,
    chunk_size: int,
    summary_c_path: Path | None = None,
    catalog_c_path: Path | None = None,
) -> dict:
    """Exact lexicographic join on two or three disjoint mod-7 projections."""
    started = time.time()
    summary_paths = [summary_a_path, summary_b_path]
    catalog_paths = [catalog_a_path, catalog_b_path]
    if summary_c_path is not None and catalog_c_path is not None:
        summary_paths.append(summary_c_path)
        catalog_paths.append(catalog_c_path)
    summaries = [json.loads(path.read_text()) for path in summary_paths]
    summary_kinds = {summary.get("experiment") for summary in summaries}
    if summary_kinds == {"p7_exceptional_tail22_catalogs"}:
        row_key = "selected_dependency_rows_7"
        projection_mode = f"injective_disjoint_mod7_22x{len(summary_paths)}_tuple"
        omitted_direction = None
    elif summary_kinds == {"p7_exceptional_omit_high_catalogs"}:
        if len(summaries) != 3:
            raise ValueError("high-direction omission requires three projections")
        row_key = "selected_conditioned_basis_rows"
        projection_mode = "injective_disjoint_mod7_22x3_omitted_high_direction"
        omitted_values = {int(summary["omitted_direction"]) for summary in summaries}
        family_hashes = {
            summary.get("selected_full_dependency_sha256") for summary in summaries
        }
        full_selections = {
            tuple(summary.get("selected_full_conditioned_basis_rows", []))
            for summary in summaries
        }
        if (
            len(omitted_values) != 1
            or len(family_hashes) != 1
            or None in family_hashes
            or len(full_selections) != 1
            or len(next(iter(full_selections))) != 66
            or any(
                summary.get("conditioned_dependency_block_is_zero") is not True
                or int(summary.get("conditioned_dependency_dimension", -1)) != 112
                or int(summary.get("conditioned_dependency_basis_rank", -1)) != 112
                or summary.get("large_catalog_treatment")
                != "exactly_eliminated_by_zero_dependency_block"
                for summary in summaries
            )
        ):
            raise ValueError("incompatible high-direction omission projections")
        omitted_direction = next(iter(omitted_values))
    else:
        raise ValueError("multi-key join requires one compatible projection family")
    for summary in summaries:
        if summary.get("projection_group_order") != 7**22:
            raise ValueError("multi-key join requires exact 22-row mod-7 projections")
    row_sets = [tuple(summary[row_key]) for summary in summaries]
    if any(len(rows) != 22 for rows in row_sets):
        raise ValueError("each projection must contain 22 dependency rows")
    if any(set(row_sets[i]) & set(row_sets[j]) for i in range(len(row_sets)) for j in range(i)):
        raise ValueError("projection row sets must be pairwise disjoint")

    mean_batch = json.loads(mean_path.read_text())
    leaf = mean_batch["leaves"][leaf_index]
    if int(leaf["leaf_index"]) != leaf_index:
        raise AssertionError("mean leaf index mismatch")
    if leaf.get("solver_status") == "INFEASIBLE":
        raise ValueError("the requested leaf is already excluded by CP-SAT")
    means = tuple(int(value) for value in leaf["scaled_means_direction_order"])
    if omitted_direction is not None:
        high_directions = [
            direction_index
            for direction_index, (row, mean) in enumerate(
                zip(mean_batch["direction_rows"], means)
            )
            if int(row["b"]) == 0
            and int(row["phase"]) == 0
            and mean > 16
        ]
        if high_directions != [omitted_direction]:
            raise ValueError("leaf does not match the omitted high-mean direction")
        boundary = tuple(int(value) for value in mean_batch["fixed_boundary"])
        for summary, catalog_path in zip(summaries, catalog_paths):
            orbit_row = next(
                (
                    row
                    for row in summary.get("orbits", [])
                    if tuple(int(value) for value in row["fixed_boundary"])
                    == boundary
                ),
                None,
            )
            if orbit_row is None or orbit_row.get("sha256") != sha256(catalog_path):
                raise ValueError("omission catalog is not certified by its summary")

    def load_one(path: Path):
        with np.load(path, allow_pickle=False) as source:
            singleton_lo = [int(source["base_p7_lo"][0])]
            singleton_hi = [int(source["base_p7_hi"][0])]
            variable = []
            sizes_local = []
            for direction_index, mean in enumerate(means):
                host_lo = source[f"d{direction_index}_m{mean}_p7_lo"]
                host_hi = source[f"d{direction_index}_m{mean}_p7_hi"]
                sizes_local.append(int(len(host_lo)))
                if len(host_lo) == 1:
                    singleton_lo.append(int(host_lo[0]))
                    singleton_hi.append(int(host_hi[0]))
                else:
                    variable.append((cp.asarray(host_lo), cp.asarray(host_hi)))
        target_lo = negate_packed(add_packed(singleton_lo, 7, 4, 11), 7, 4, 11)
        target_hi = negate_packed(add_packed(singleton_hi, 7, 4, 11), 7, 4, 11)
        return variable, tuple(sizes_local), target_lo, target_hi

    variable_a, all_sizes_a, target_a_lo, target_a_hi = load_one(catalog_a_path)
    variable_b, all_sizes_b, target_b_lo, target_b_hi = load_one(catalog_b_path)
    loaded_c = load_one(catalog_c_path) if catalog_c_path is not None else None
    if all_sizes_a != all_sizes_b or tuple(len(row[0]) for row in variable_a) != tuple(len(row[0]) for row in variable_b):
        raise AssertionError("paired projection catalogs are not aligned")
    if loaded_c is not None and (
        loaded_c[1] != all_sizes_a
        or tuple(len(row[0]) for row in loaded_c[0]) != tuple(len(row[0]) for row in variable_a)
    ):
        raise AssertionError("third projection catalog is not aligned")
    sizes = tuple(len(row[0]) for row in variable_a)
    first_indices, second_indices = choose_partition(sizes, max_build)
    zero = (cp.asarray([0], dtype=cp.uint64), cp.asarray([0], dtype=cp.uint64))

    def groups(variable):
        first = [variable[index] for index in first_indices]
        second = [variable[index] for index in second_indices]
        return first + [zero] * (3 - len(first)), second + [zero] * (3 - len(second))

    first_a, second_a = groups(variable_a)
    first_b, second_b = groups(variable_b)
    if loaded_c is not None:
        variable_c, _all_sizes_c, target_c_lo, target_c_hi = loaded_c
        first_c, second_c = groups(variable_c)
    build_count = math.prod(len(item[0]) for item in first_a)
    probe_count = math.prod(len(item[0]) for item in second_a)
    module = cp.RawModule(code=CUDA_SOURCE, options=("--std=c++14",))
    make_sum = module.get_function("make_sum_keys7_22")
    count_pairs = module.get_function("count_needed_sorted_pairs7_22")
    count_triples = module.get_function("count_needed_sorted_triples7_22")
    block = 256

    def make_build(group):
        keys = cp.empty(build_count, dtype=cp.uint64)
        make_sum(
            ((build_count + block - 1) // block,), (block,),
            (
                group[0][0], group[0][1], np.uint64(len(group[0][0])),
                group[1][0], group[1][1], np.uint64(len(group[1][0])),
                group[2][0], group[2][1], np.uint64(len(group[2][0])),
                np.uint64(0), np.uint64(build_count), keys,
            ),
        )
        return keys

    build_a = make_build(first_a)
    build_b = make_build(first_b)
    build_c = make_build(first_c) if loaded_c is not None else None
    lex_keys = cp.stack((build_c, build_b, build_a)) if build_c is not None else cp.stack((build_b, build_a))
    order = cp.lexsort(lex_keys)
    del lex_keys
    sorted_a = build_a[order]
    sorted_b = build_b[order]
    sorted_c = build_c[order] if build_c is not None else None
    del build_a, build_b, build_c, order

    match_count = cp.zeros(1, dtype=cp.uint64)
    for offset in range(0, probe_count, chunk_size):
        count = min(chunk_size, probe_count - offset)
        if loaded_c is not None:
            count_triples(
                ((count + block - 1) // block,), (block,),
                (
                    second_a[0][0], second_a[0][1], second_b[0][0], second_b[0][1], second_c[0][0], second_c[0][1], np.uint64(len(second_a[0][0])),
                    second_a[1][0], second_a[1][1], second_b[1][0], second_b[1][1], second_c[1][0], second_c[1][1], np.uint64(len(second_a[1][0])),
                    second_a[2][0], second_a[2][1], second_b[2][0], second_b[2][1], second_c[2][0], second_c[2][1], np.uint64(len(second_a[2][0])),
                    np.uint64(target_a_lo), np.uint64(target_a_hi),
                    np.uint64(target_b_lo), np.uint64(target_b_hi),
                    np.uint64(target_c_lo), np.uint64(target_c_hi),
                    np.uint64(offset), np.uint64(count),
                    sorted_a, sorted_b, sorted_c, np.uint64(build_count), match_count,
                ),
            )
        else:
            count_pairs(
                ((count + block - 1) // block,), (block,),
                (
                    second_a[0][0], second_a[0][1], second_b[0][0], second_b[0][1], np.uint64(len(second_a[0][0])),
                    second_a[1][0], second_a[1][1], second_b[1][0], second_b[1][1], np.uint64(len(second_a[1][0])),
                    second_a[2][0], second_a[2][1], second_b[2][0], second_b[2][1], np.uint64(len(second_a[2][0])),
                    np.uint64(target_a_lo), np.uint64(target_a_hi),
                    np.uint64(target_b_lo), np.uint64(target_b_hi),
                    np.uint64(offset), np.uint64(count),
                    sorted_a, sorted_b, np.uint64(build_count), match_count,
                ),
            )
    exact_matches = int(match_count.get()[0])
    return {
        "experiment": "p7_exceptional_projected_join_gpu",
        "status": "complete_exact_selected_dependency_gpu_join",
        "projection_mode": projection_mode,
        "p": 7,
        "c_H": -1,
        "fixed_boundary": mean_batch["fixed_boundary"],
        "leaf_index": leaf_index,
        "fixed_scaled_means": list(means),
        "projection_summaries": [str(path) for path in summary_paths],
        "projection_summary_sha256": [sha256(path) for path in summary_paths],
        "projection_catalogs": [str(path) for path in catalog_paths],
        "projection_catalog_sha256": [sha256(path) for path in catalog_paths],
        row_key: [list(rows) for rows in row_sets],
        "omitted_high_mean_direction": omitted_direction,
        "projection_group_order": (7**22) ** len(summary_paths),
        "direction_catalog_sizes": list(all_sizes_a),
        "variable_catalog_sizes": list(sizes),
        "partition": [list(first_indices), list(second_indices)],
        "build_signatures": build_count,
        "probe_signatures": probe_count,
        "exact_join_method": f"lexicographic_uint64_{len(summary_paths)}tuple_full_probe_binary_search",
        "exact_projected_matches": exact_matches,
        "projected_modularly_infeasible": exact_matches == 0,
        "finite_mean_allocation_exclusion": exact_matches == 0,
        "gpu": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection-summary", type=Path, required=True)
    parser.add_argument("--projection-catalog", type=Path, required=True)
    parser.add_argument("--projection-summary-2", type=Path)
    parser.add_argument("--projection-catalog-2", type=Path)
    parser.add_argument("--projection-summary-3", type=Path)
    parser.add_argument("--projection-catalog-3", type=Path)
    parser.add_argument("--mean-batch", type=Path, required=True)
    parser.add_argument("--leaf-index", type=int, required=True)
    parser.add_argument("--max-build", type=int, default=150_000_000)
    parser.add_argument("--chunk-size", type=int, default=20_000_000)
    parser.add_argument("--candidate-capacity", type=int, default=2_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.projection_summary_2 is None) != (args.projection_catalog_2 is None):
        raise ValueError("both second-projection paths must be supplied together")
    if (args.projection_summary_3 is None) != (args.projection_catalog_3 is None):
        raise ValueError("both third-projection paths must be supplied together")
    if args.projection_summary_3 is not None and args.projection_summary_2 is None:
        raise ValueError("a third projection requires a second projection")
    summary_kind = json.loads(args.projection_summary.read_text()).get("experiment")
    if args.projection_summary_2 is not None:
        out = run_mod7_pair(
            args.projection_summary,
            args.projection_catalog,
            args.projection_summary_2,
            args.projection_catalog_2,
            args.mean_batch,
            args.leaf_index,
            args.max_build,
            args.chunk_size,
            args.projection_summary_3,
            args.projection_catalog_3,
        )
    elif summary_kind in (
        "p7_exceptional_tail22_catalogs",
        "p7_exceptional_p3_40_catalogs",
        "p7_exceptional_p5_27_catalogs",
    ):
        out = run_tail22(
            args.projection_summary,
            args.projection_catalog,
            args.mean_batch,
            args.leaf_index,
            args.max_build,
            args.chunk_size,
        )
    else:
        out = run(
            args.projection_summary,
            args.projection_catalog,
            args.mean_batch,
            args.leaf_index,
            args.max_build,
            args.chunk_size,
            args.candidate_capacity,
        )
    if args.output is not None:
        atomic_json(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
