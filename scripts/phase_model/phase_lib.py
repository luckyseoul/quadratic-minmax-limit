"""Split-torus phase model of the Paley conference form.

Coordinates.  q = 1 (mod 4) prime, g a primitive root, M = (q-1)/2.  A sign
vector on GF(q)^* is written z_s = x(g^s), s in Z/(q-1).  Then

    chi(g^s - g^t) = (-1)^s h(t - s),   h(d) = chi(1 - g^d),

and with J_k = sum_d h(d) e(kd/(q-1)) = J(chi, psi^k) (a Jacobi sum):

    sum_{s,t} chi(g^s - g^t) z_s z_t = (q-1) sum_k J_k zhat(k) zhat(M - k).

All evaluations of Q used for certificates are exact integer arithmetic
(Kronecker substitution); floating point is used only to locate candidates.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "parity_ceiling"))
from paley_lib import chi_table, is_prime, prim_root  # noqa: E402


def split_model(q):
    """Return (g, chi, pw, h, J): primitive root, character table, powers g^s, h(d), Jacobi sums J_k."""
    assert is_prime(q) and q % 4 == 1
    g = prim_root(q)
    chi = np.asarray(chi_table(q), dtype=np.int64)
    pw = np.empty(q - 1, dtype=np.int64)
    x = 1
    for s in range(q - 1):
        pw[s] = x
        x = x * g % q
    h = chi[(1 - pw) % q]
    J = np.conj(np.fft.fft(h.astype(float)))
    return g, chi, pw, h, J


def square_state(q, k, delta):
    """Discrete square-wave pair of frequency k and relative shift delta (torus coordinates)."""
    n1 = q - 1
    M = n1 // 2
    g = math.gcd(k, M)
    M1, k1 = M // g, k // g
    assert M1 % 2 == 0, "M / gcd(k, M) must be even"
    u = np.arange(M)
    z = np.empty(n1, dtype=np.int64)
    z[0::2] = np.where((k1 * u) % M1 < M1 // 2, 1, -1)
    z[1::2] = np.where((k1 * u - delta) % M1 < M1 // 2, 1, -1)
    return z


def to_elements(q, pw, z, x0=1):
    """Torus coordinates -> vector indexed by the elements of GF(q)."""
    x = np.empty(q, dtype=np.int64)
    x[pw] = z
    x[0] = x0
    return x


def exact_fields(q, chi, x_inf, x):
    """Exact (Cx) for the Paley conference matrix of order q+1: returns ((Cx)_inf, (Cx)_i for i in GF(q))."""
    a = (chi + 1).astype(np.uint32)
    b = (x + 1).astype(np.uint32)

    def pack(v):
        buf = np.zeros((len(v), 3), dtype=np.uint8)
        buf[:, 0] = v & 255
        buf[:, 1] = (v >> 8) & 255
        buf[:, 2] = (v >> 16) & 255
        return int.from_bytes(buf.tobytes(), "little")

    assert 4 * q < 2 ** 24
    prod = pack(a) * pack(b)
    raw = np.frombuffer(prod.to_bytes(3 * (2 * q), "little"), dtype=np.uint8).reshape(2 * q, 3).astype(np.int64)
    co = raw[:, 0] + (raw[:, 1] << 8) + (raw[:, 2] << 16)
    conv = co[:q] + co[q:] - int(x.sum()) - q
    return int(x.sum()), x_inf + conv


def exact_Q(q, chi, x_inf, x):
    """Exact Q_C(x) = sum_{i<j} C_ij x_i x_j for the Paley conference matrix of order q+1."""
    y_inf, y = exact_fields(q, chi, x_inf, x)
    two = x_inf * y_inf + int((x * y).sum())
    assert two % 2 == 0
    return two // 2


def rho_star(q, phi, k, delta):
    """Closed formula: 2 Q* / ((q-1) sqrt(q)) for square_state(q, k, delta)."""
    n1 = q - 1
    M = n1 // 2
    g = math.gcd(k, M)
    M1, k1 = M // g, k // g
    m = np.arange(1, M1, 2)
    return float((4 / M1 ** 2) * np.sum(np.sin(phi[(m * k) % n1]) * np.sin(np.pi * m * (2 * delta + k1) / M1)
                                        / np.sin(np.pi * m / M1) ** 2))


def rho_star_all_t(q, phi, k):
    """For gcd(k, M) = 1: rho*(k, delta) for every t = 2 delta + k (index t in [0, q-1); only odd t are states)."""
    n1 = q - 1
    M = n1 // 2
    m = np.arange(1, M, 2)
    c = np.zeros(n1, dtype=complex)
    c[m] = (4 / M ** 2) * np.sin(phi[(m * k) % n1]) / np.sin(np.pi * m / M) ** 2
    return np.imag(np.fft.ifft(c) * n1)


def admissible(q):
    """Odd k with gcd(k, M) = 1."""
    n1 = q - 1
    k = np.arange(1, n1, 2)
    return k[np.gcd(k, n1 // 2) == 1]
