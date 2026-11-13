"""Best square-wave state of the phase model for each prime q, with exact integer certificates.

usage: scan.py q1 q2 ...   (one JSON record per line on stdout)

For q <= 40000 every admissible frequency k and every shift is examined; above that
the frequencies are preselected by the first 16 harmonics and the record is a lower
bound for the best square-wave state.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase_lib import admissible, exact_Q, rho_star_all_t, split_model, square_state, to_elements  # noqa: E402

FULL_SCAN_LIMIT = 40000
WEYL_VECTORS = [(1, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
                (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, 1, 1), (2, 1, 0), (2, -1, 0),
                (1, 2, 0), (1, -2, 0), (2, 0, 1), (2, 1, 1), (3, -1, 0), (1, -2, 1)]
BOXES = [(1, 0.3), (2, 0.3), (3, 0.3), (3, 0.5), (4, 0.5)]


def preselect(q, phi, r=16, nt=64, top=300):
    n1 = q - 1
    M = n1 // 2
    ks = admissible(q)
    m = np.arange(1, 2 * r, 2)
    w = (8 / M ** 2) / np.sin(np.pi * m / M) ** 2
    A = np.stack([w[i] * np.sin(phi[(mm * ks) % n1]) for i, mm in enumerate(m)], axis=1)
    sc = A @ np.sin(np.outer(m, np.linspace(0, 2 * np.pi, nt, endpoint=False)))
    return ks[np.argsort(sc.max(axis=1))[::-1][:top]]


def best_square(q, phi, full):
    n1 = q - 1
    best = (-9.0, None, None)
    for k in (admissible(q) if full else preselect(q, phi)):
        F = rho_star_all_t(q, phi, int(k))
        F[0::2] = -9.0
        t = int(np.argmax(F))
        if F[t] > best[0]:
            best = (float(F[t]), int(k), t)
    rho, k, t = best
    return rho, k, ((t - k) // 2) % (n1 // 2), t


def fft_Q(q, chif, xi, x):
    conv = np.rint(np.fft.ifft(chif * np.fft.fft(x.astype(float))).real).astype(np.int64)
    y = xi + conv
    yi = int(x.sum())
    return (xi * yi + int((x * y).sum())) // 2, yi, y


def polish(q, chi, xi, x, rounds=400):
    """Deterministic batch 1-opt ascent (floating point guides, the result is re-certified exactly)."""
    chif = np.fft.fft(chi.astype(float))
    Q, yi, y = fft_Q(q, chif, xi, x)
    for _ in range(rounds):
        if xi * yi < 0:
            xi = -xi
            Q, yi, y = fft_Q(q, chif, xi, x)
            continue
        hf = x * y
        neg = np.where(hf < 0)[0]
        if len(neg) == 0:
            break
        order = neg[np.argsort(hf[neg], kind="stable")]
        frac, moved = 1.0, False
        while True:
            take = order[:max(1, int(len(order) * frac))]
            x2 = x.copy()
            x2[take] *= -1
            Q2, yi2, y2 = fft_Q(q, chif, xi, x2)
            if Q2 > Q:
                x, Q, yi, y, moved = x2, Q2, yi2, y2, True
                break
            if len(take) == 1:
                break
            frac /= 2
        if not moved:
            break
    return xi, x


def statistics(q, phi):
    n1 = q - 1
    ks = np.arange(1, n1, 2)
    p1, p3, p5 = phi[ks % n1], phi[(3 * ks) % n1], phi[(5 * ks) % n1]
    weyl = max(abs(np.mean(np.exp(1j * (v[0] * p1 + v[1] * p3 + v[2] * p5)))) for v in WEYL_VECTORS)
    boxes = []
    for r, th in BOXES:
        eps = [(-1) ** ((m - 1) // 2) for m in range(1, 2 * r, 2)]
        P = np.stack([eps[i] * np.sin(phi[(m * ks) % n1]) for i, m in enumerate(range(1, 2 * r, 2))], axis=1)
        cnt = int((P >= math.cos(th)).all(axis=1).sum() + (-P >= math.cos(th)).all(axis=1).sum())
        boxes.append({"harmonics": r, "theta": th, "count": cnt, "haar_prediction": 2 * len(ks) * (th / math.pi) ** r})
    return float(weyl), boxes


def record(q):
    t0 = time.time()
    g, chi, pw, h, J = split_model(q)
    phi = np.angle(J)
    n1 = q - 1
    M = n1 // 2
    S = 0.5 * (q + 1) * math.sqrt(q)
    full = q <= FULL_SCAN_LIMIT
    rho, k, delta, t = best_square(q, phi, full)
    x = to_elements(q, pw, square_state(q, k, delta), 1)
    Q_sq = exact_Q(q, chi, 1, x)
    assert abs(Q_sq - (1 + 0.5 * n1 * math.sqrt(q) * rho)) < 0.5
    xi, xp = polish(q, chi, 1, x.copy())
    Q_pol = exact_Q(q, chi, xi, xp)
    weyl, boxes = statistics(q, phi)
    return {
        "q": q, "order": q + 1, "primitive_root": int(g), "full_scan": full, "k": k, "delta": delta,
        "t_over_M": t / M, "Q_square": int(Q_sq), "rho_square": Q_sq / S, "Q_polished": int(Q_pol),
        "rho_polished": Q_pol / S,
        "harmonic_products": [round(float(np.sin(phi[(m * k) % n1]) * np.sin(np.pi * m * t / M)), 4)
                              for m in range(1, 16, 2)],
        "weyl_max": weyl, "weyl_max_times_sqrt_q": weyl * math.sqrt(q), "weyl_vectors": len(WEYL_VECTORS),
        "boxes": boxes, "seconds": round(time.time() - t0, 1),
    }


if __name__ == "__main__":
    for q in [int(a) for a in sys.argv[1:]]:
        print(json.dumps(record(q)), flush=True)
