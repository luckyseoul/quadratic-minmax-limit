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
                "-fopenmp",
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
    loaded.qml_gf2x_cyclic_product.argtypes = [
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
    ]
    loaded.qml_gf2x_cyclic_product.restype = ctypes.c_long
    loaded.qml_gf2x_cyclic_star_product.argtypes = [
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
    ]
    loaded.qml_gf2x_cyclic_star_product.restype = ctypes.c_long
    loaded.qml_field_primitive.argtypes = [
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
    ]
    loaded.qml_field_primitive.restype = ctypes.c_uint64
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
    loaded.qml_selected_line_bins.argtypes = [
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_long,
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
    ]
    loaded.qml_selected_line_bins.restype = ctypes.c_int
    loaded.qml_selected_line_bins_wide.argtypes = [
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_uint64,
        ctypes.c_uint64,
        ctypes.c_uint64,
        ctypes.c_uint64,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_long,
        ctypes.c_long,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_long,
    ]
    loaded.qml_selected_line_bins_wide.restype = ctypes.c_int
    loaded.qml_selected_line_counts_mod4_wide.argtypes = (
        loaded.qml_selected_line_bins_wide.argtypes
    )
    loaded.qml_selected_line_counts_mod4_wide.restype = ctypes.c_int
    _LIB = loaded
    return loaded


def available() -> bool:
    try:
        _load()
    except (OSError, subprocess.SubprocessError):
        return False
    return True


def field_primitive(p: int, ia: int, ib: int) -> int:
    """Return the least encoded primitive element of F_(p^2)."""
    value = int(_load().qml_field_primitive(p, ia, ib))
    if value == 0:
        raise RuntimeError(f"native primitive-element search failed for p={p}")
    return value


def gcd_bits(left: int, right: int) -> int:
    """Polynomial gcd over F2 for little-endian packed Python integers."""
    if left < 0 or right < 0:
        raise ValueError("packed polynomials must be nonnegative")
    left_bytes = left.to_bytes(max(1, (left.bit_length() + 7) // 8), "little")
    right_bytes = right.to_bytes(max(1, (right.bit_length() + 7) // 8), "little")
    left_buffer = (ctypes.c_uint8 * len(left_bytes)).from_buffer_copy(left_bytes)
    right_buffer = (ctypes.c_uint8 * len(right_bytes)).from_buffer_copy(right_bytes)
    # Usually a gcd fits in the smaller nonzero operand.  If either operand
    # is zero, however, gcd(f, 0) = f, so reserve the larger packed length.
    capacity = max(1, max(len(left_bytes), len(right_bytes)))
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


def cyclic_product_bits(left: int, right: int, order: int) -> int:
    """Multiply packed F2 polynomials modulo X^order+1."""
    if left < 0 or right < 0 or order < 1:
        raise ValueError("polynomials must be nonnegative and order positive")
    left_bytes = left.to_bytes(max(1, (left.bit_length() + 7) // 8), "little")
    right_bytes = right.to_bytes(max(1, (right.bit_length() + 7) // 8), "little")
    left_buffer = (ctypes.c_uint8 * len(left_bytes)).from_buffer_copy(left_bytes)
    right_buffer = (ctypes.c_uint8 * len(right_bytes)).from_buffer_copy(right_bytes)
    capacity = (order + 7) // 8
    output = (ctypes.c_uint8 * capacity)()
    written = _load().qml_gf2x_cyclic_product(
        left_buffer,
        len(left_bytes),
        right_buffer,
        len(right_bytes),
        order,
        output,
        capacity,
    )
    if written < 0:
        raise RuntimeError(f"NTL cyclic product failed with {written}")
    return int.from_bytes(bytes(output[:written]), "little")


def cyclic_star_product_bits(left: int, right: int, order: int) -> int:
    """Return left(X) * right(X^-1) modulo X^order+1 over F2."""
    if left < 0 or right < 0 or order < 1:
        raise ValueError("polynomials must be nonnegative and order positive")
    left_bytes = left.to_bytes(max(1, (left.bit_length() + 7) // 8), "little")
    right_bytes = right.to_bytes(max(1, (right.bit_length() + 7) // 8), "little")
    left_buffer = (ctypes.c_uint8 * len(left_bytes)).from_buffer_copy(left_bytes)
    right_buffer = (ctypes.c_uint8 * len(right_bytes)).from_buffer_copy(right_bytes)
    capacity = (order + 7) // 8
    output = (ctypes.c_uint8 * capacity)()
    written = _load().qml_gf2x_cyclic_star_product(
        left_buffer,
        len(left_bytes),
        right_buffer,
        len(right_bytes),
        order,
        output,
        capacity,
    )
    if written < 0:
        raise RuntimeError(f"NTL cyclic star product failed with {written}")
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


def selected_line_bins(
    p: int,
    ia: int,
    ib: int,
    generator: int,
    omega: int,
    sigma: int,
    pole_t: int,
    levels: list[int],
    orders: list[int],
    force_wide: bool = False,
) -> tuple[np.ndarray, list[int]]:
    """Fold selected affine-line pullbacks without traversing either orbit."""
    if not levels or not orders:
        raise ValueError("levels and orders must be nonempty")
    offsets = np.cumsum([0, *orders[:-1]], dtype=np.uint32)
    total_order = sum(orders)
    level_array = np.ascontiguousarray(levels, dtype=np.uint32)
    order_array = np.ascontiguousarray(orders, dtype=np.uint32)
    output = np.zeros((2, len(levels), total_order), dtype=np.uint8)
    u32_pointer = ctypes.POINTER(ctypes.c_uint32)
    u8_pointer = ctypes.POINTER(ctypes.c_uint8)
    function = (
        _load().qml_selected_line_bins_wide
        if force_wide or p * p >= 1 << 32
        else _load().qml_selected_line_bins
    )
    status = function(
        p,
        ia,
        ib,
        generator,
        omega,
        sigma,
        pole_t,
        level_array.ctypes.data_as(u32_pointer),
        len(levels),
        order_array.ctypes.data_as(u32_pointer),
        offsets.ctypes.data_as(u32_pointer),
        len(orders),
        total_order,
        output.ctypes.data_as(u8_pointer),
        output.nbytes,
    )
    if status:
        raise RuntimeError(f"native selected-line folding failed with {status}")
    return output, [int(value) for value in offsets]


def selected_line_counts_mod4(
    p: int,
    ia: int,
    ib: int,
    generator: int,
    omega: int,
    sigma: int,
    pole_t: int,
    levels: list[int],
    orders: list[int],
) -> tuple[np.ndarray, list[int]]:
    """Fold selected affine lines into exponent-class counts modulo four."""
    if not levels or not orders:
        raise ValueError("levels and orders must be nonempty")
    offsets = np.cumsum([0, *orders[:-1]], dtype=np.uint32)
    total_order = sum(orders)
    level_array = np.ascontiguousarray(levels, dtype=np.uint32)
    order_array = np.ascontiguousarray(orders, dtype=np.uint32)
    output = np.zeros((2, len(levels), total_order), dtype=np.uint8)
    u32_pointer = ctypes.POINTER(ctypes.c_uint32)
    u8_pointer = ctypes.POINTER(ctypes.c_uint8)
    status = _load().qml_selected_line_counts_mod4_wide(
        p,
        ia,
        ib,
        generator,
        omega,
        sigma,
        pole_t,
        level_array.ctypes.data_as(u32_pointer),
        len(levels),
        order_array.ctypes.data_as(u32_pointer),
        offsets.ctypes.data_as(u32_pointer),
        len(orders),
        total_order,
        output.ctypes.data_as(u8_pointer),
        output.nbytes,
    )
    if status:
        raise RuntimeError(
            f"native selected-line mod-four folding failed with {status}"
        )
    return output, [int(value) for value in offsets]
