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
