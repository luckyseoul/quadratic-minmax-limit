"""A380/SYCL emit wrapper. Does not import cupy (safe on jellyfin; not on V100)."""
from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np
from gpu_inner_fast import _vmap_from_av

_LIB = None


def _preload_compiler_sycl():
    """intelpython ships an older libsycl.so.8 that lacks submit_with_event_impl."""
    override = os.environ.get("LIBSYCL_SO")
    cands = []
    if override:
        cands.append(Path(override))
    cands.extend((
        Path("/opt/intel/oneapi/compiler/2025.2/lib/libsycl.so.8"),
        Path("/opt/intel/oneapi/compiler/latest/lib/libsycl.so.8"),
    ))
    for p in cands:
        if p.is_file():
            ctypes.CDLL(str(p), mode=ctypes.RTLD_GLOBAL)
            return str(p)
    return None


def _lib():
    global _LIB
    if _LIB is None:
        _preload_compiler_sycl()
        here = Path(__file__).resolve().parent
        cand = [
            here / "libgpu_gen_sycl.so",
            here / "gpu_gen_sycl.so",
            Path("/tmp/maxplus_enum/libgpu_gen_sycl.so"),
            Path("/tmp/e1work/libgpu_gen_sycl.so"),
            Path("/tmp/e1work/gpu_gen_sycl.so"),
        ]
        env = os.environ.get("GPU_GEN_SYCL_SO")
        if env:
            cand.insert(0, Path(env))
        path = next((p for p in cand if p.is_file()), None)
        if path is None:
            raise FileNotFoundError("gpu_gen_sycl.so not found; compile gpu_gen_sycl.cpp")
        _LIB = ctypes.CDLL(str(path))
        _LIB.sycl_device_name.restype = ctypes.c_char_p
        _LIB.emit_chunk_sycl.restype = ctypes.c_int64
        _LIB.emit_chunk_sycl.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_uint64,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int64,
        ]
    return _LIB


def sycl_device_name() -> str:
    n = _lib().sycl_device_name()
    return n.decode() if n else ""


def gen_candidates_sycl(
    p, k, av, af, an, aull, UU, c0, u_lo, u_hi, probes, cprobes, thi, tlo,
    codes, fsums,
):
    nci = int(u_hi - u_lo)
    if nci <= 0:
        return 0
    ncombo = 1
    for _ in range(k - 1):
        ncombo *= int(p)
    vmap = _vmap_from_av(p, k, av, an)
    av_ = np.ascontiguousarray(av, dtype=np.int32)
    af_ = np.ascontiguousarray(af, dtype=np.int32)
    an_ = np.ascontiguousarray(an, dtype=np.int32)
    aull_ = np.ascontiguousarray(aull, dtype=np.uint8)
    UU_ = np.ascontiguousarray(UU, dtype=np.int32)
    vmap_ = np.ascontiguousarray(vmap, dtype=np.int32)
    probes_ = np.ascontiguousarray(probes, dtype=np.int16)
    cprobes_ = np.ascontiguousarray(cprobes, dtype=np.int16)
    codes_ = np.ascontiguousarray(codes)
    fsums_ = np.ascontiguousarray(fsums)
    n = int(
        _lib().emit_chunk_sycl(
            int(p), int(k), int(u_lo), nci, int(c0), int(thi), int(tlo),
            int(probes_.shape[2]), int(2 * p), ctypes.c_uint64(ncombo),
            av_.ctypes.data, af_.ctypes.data, an_.ctypes.data, aull_.ctypes.data,
            UU_.ctypes.data, vmap_.ctypes.data, probes_.ctypes.data, cprobes_.ctypes.data,
            codes_.ctypes.data, fsums_.ctypes.data, int(codes_.shape[0]),
        )
    )
    if n < 0:
        raise RuntimeError("SYCL emit failed (see stderr)")
    if n > codes.shape[0]:
        raise RuntimeError(f"candidate overflow {n} > {codes.shape[0]}")
    if n and codes_ is not codes:
        codes[:n] = codes_[:n]
        fsums[:n] = fsums_[:n]
    return n
