"""Optional NTL-backed gcd for packed polynomials over F2.

The first call builds a tiny ctypes bridge in /tmp.  If NTL or a C++ compiler
is unavailable, callers can retain their pure-Python fallback.
"""
from __future__ import annotations

import ctypes
import hashlib
import subprocess
import tempfile
from pathlib import Path

import numpy as np


SOURCE = Path(__file__).with_suffix(".cpp")
_LIB = None


def _load():
    global _LIB
    if _LIB is not None:
        return _LIB
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()[:16]
    library = Path(tempfile.gettempdir()) / f"qml_gf2x_ntl_{digest}.so"
    if not library.exists():
        temporary = library.with_suffix(".so.tmp")
        subprocess.run(
            [
                "g++",
                "-O3",
                "-std=c++17",
                "-fPIC",
                "-shared",
                str(SOURCE),
                "-lntl",
                "-lgmp",
                "-o",
                str(temporary),
            ],
            check=True,
        )
        temporary.replace(library)
    loaded = ctypes.CDLL(str(library))
    loaded.qml_gf2x_gcd.argtypes = [
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
    ]
    loaded.qml_gf2x_gcd.restype = ctypes.c_long
    loaded.qml_field_orbits.argtypes = [
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_long,
    ]
    loaded.qml_field_orbits.restype = ctypes.c_int
    _LIB = loaded
    return loaded


def available() -> bool:
    try:
        _load()
    except (OSError, subprocess.SubprocessError):
        return False
    return True


def gcd_bits(left: int, right: int) -> int:
    """Polynomial gcd over F2 for little-endian packed Python integers."""
    if left < 0 or right < 0:
        raise ValueError("packed polynomials must be nonnegative")
    left_bytes = left.to_bytes(max(1, (left.bit_length() + 7) // 8), "little")
    right_bytes = right.to_bytes(max(1, (right.bit_length() + 7) // 8), "little")
    left_buffer = (ctypes.c_uint8 * len(left_bytes)).from_buffer_copy(left_bytes)
    right_buffer = (ctypes.c_uint8 * len(right_bytes)).from_buffer_copy(right_bytes)
    capacity = max(1, min(len(left_bytes), len(right_bytes)))
    output = (ctypes.c_uint8 * capacity)()
    written = _load().qml_gf2x_gcd(
        left_buffer,
        len(left_bytes),
        right_buffer,
        len(right_bytes),
        output,
        capacity,
    )
    if written < 0:
        raise RuntimeError(f"NTL gcd output needs {-written} bytes, had {capacity}")
    return int.from_bytes(bytes(output[:written]), "little")


def field_orbits(
    p: int,
    ia: int,
    ib: int,
    generator: int,
    omega: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate both multiplicative point-orbits and an inverse log table."""
    q = p * p
    orbit_length = (q - 1) // 2
    if orbit_length >= 1 << 31:
        raise ValueError("native uint32 orbit encoding requires p < 65536")
    square = np.empty(orbit_length, dtype=np.uint32)
    nonsquare = np.empty(orbit_length, dtype=np.uint32)
    logarithm = np.empty(q, dtype=np.uint32)
    pointer = ctypes.POINTER(ctypes.c_uint32)
    status = _load().qml_field_orbits(
        p,
        ia,
        ib,
        generator,
        omega,
        square.ctypes.data_as(pointer),
        nonsquare.ctypes.data_as(pointer),
        logarithm.ctypes.data_as(pointer),
        orbit_length,
    )
    if status:
        raise RuntimeError(f"native field orbit generation failed with {status}")
    return square, nonsquare, logarithm


def field_two_orbits(
    p: int,
    ia: int,
    ib: int,
    generator: int,
    omega: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate the two point-orbits without allocating an inverse-log table."""
    q = p * p
    orbit_length = (q - 1) // 2
    if orbit_length >= 1 << 31:
        raise ValueError("native uint32 orbit encoding requires p < 65536")
    square = np.empty(orbit_length, dtype=np.uint32)
    nonsquare = np.empty(orbit_length, dtype=np.uint32)
    pointer = ctypes.POINTER(ctypes.c_uint32)
    status = _load().qml_field_orbits(
        p,
        ia,
        ib,
        generator,
        omega,
        square.ctypes.data_as(pointer),
        nonsquare.ctypes.data_as(pointer),
        None,
        orbit_length,
    )
    if status:
        raise RuntimeError(f"native field orbit generation failed with {status}")
    return square, nonsquare
