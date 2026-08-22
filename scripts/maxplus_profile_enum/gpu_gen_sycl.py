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
        if hasattr(_LIB, "sycl_test_load"):
            _LIB.sycl_test_load.restype = ctypes.c_int
            _LIB.sycl_test_load.argtypes = [
                ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                ctypes.c_int, ctypes.c_int,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p,
            ]
            _LIB.sycl_test_batch.restype = ctypes.c_int
            _LIB.sycl_test_batch.argtypes = [
                ctypes.c_int64, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
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


def sycl_has_tester() -> bool:
    return hasattr(_lib(), "sycl_test_load")


def sycl_test_load(p, k, q, thi, tlo, UU, av, af, Tm, cont, cov) -> None:
    UU_ = np.ascontiguousarray(UU, dtype=np.int32)
    av_ = np.ascontiguousarray(av, dtype=np.int32)
    af_ = np.ascontiguousarray(af, dtype=np.int32)
    Tm_ = np.ascontiguousarray(Tm, dtype=np.int32)
    cont_ = np.ascontiguousarray(cont, dtype=np.int16)
    cov_ = np.ascontiguousarray(cov, dtype=np.int8)
    rc = int(
        _lib().sycl_test_load(
            int(p), int(k), int(q), int(thi), int(tlo), int(UU_.shape[0]),
            UU_.ctypes.data, av_.ctypes.data, af_.ctypes.data, Tm_.ctypes.data,
            cont_.ctypes.data, cov_.ctypes.data,
        )
    )
    if rc != 0:
        raise RuntimeError(f"sycl_test_load rc={rc}")


def sycl_test_batch(codes, fsums):
    B = int(len(codes))
    if B <= 0:
        z = np.zeros(0, dtype=np.int32)
        return z, z
    codes_ = np.ascontiguousarray(codes, dtype=np.int64)
    fsums_ = np.ascontiguousarray(fsums, dtype=np.int64)
    f0 = np.empty(B, dtype=np.int32)
    fl = np.empty(B, dtype=np.int32)
    n_f0 = np.zeros(1, dtype=np.int32)
    n_fl = np.zeros(1, dtype=np.int32)
    rc = int(
        _lib().sycl_test_batch(
            ctypes.c_int64(B), codes_.ctypes.data, fsums_.ctypes.data,
            f0.ctypes.data, fl.ctypes.data, n_f0.ctypes.data, n_fl.ctypes.data,
        )
    )
    if rc != 0:
        raise RuntimeError(f"sycl_test_batch rc={rc}")
    # atomic fill is unordered; np.where is sorted
    return np.sort(f0[: int(n_f0[0])]), np.sort(fl[: int(n_fl[0])])
