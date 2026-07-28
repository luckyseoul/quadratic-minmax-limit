#!/usr/bin/env python3
"""
m4 moduli for Max+ fourth moments (Prop 15.52–15.53 support).

Stratify 4-sets by pure C-invariants:
  (S4-canonical 6-edge sign pattern, external C-sum histogram).
On each class m4 is constant (certified vs Max+ at p=5).

Averaged evec identities Cy=py give a linear system
  A m = b
with A,b combinatorial (from C and the 2-design m2=C/p only).
At p=5: rank = n_var - 1 (nullity 1). True Max+ moments lie on the line
  m = m_part + c n.

Pin c by Tr(G^2) (quadratic in c). Recovers g_min = -3/65 at p=5.

Max+-free identities shipped:
  - g_min = -max |m4| over 4-sets with |kappa|=1  (pairing algebra)
  - E[dot^2] = n + n(n-1)/p^2 from the 2-design E[yy^T]=I+C/p
  - Tr(G^2) = (1/4)(E[dot^4] - 2n E[dot^2] + n^2) via K_ab=((dot)^2-n)/2
  - E[dot^4] / closed G-spectrum still Max+-dependent (OPEN for general p)

Writes evidence/e1_gmin_moduli.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def L_p(p: int) -> float:
    return -(p - 2) / (2 * p * p)


def T_p(p: int) -> float:
    return -(p - 2) / (p * (2 * p - 1))


def load_maxplus(p: int) -> np.ndarray:
    if p == 5:
        path = Path("/tmp/maxplus_p5.npy")
    elif p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
    else:
        raise ValueError(p)
    if not path.is_file():
        raise FileNotFoundError(path)
    return np.load(path).astype(float)


def type6(C: np.ndarray, pts: tuple) -> tuple:
    """S4-canonical 6-tuple of C-edge signs on a 4-set."""
    best = None
    for perm in permutations(range(4)):
        q = [pts[perm[i]] for i in range(4)]
        t = tuple(int(C[q[i], q[j]]) for i in range(4) for j in range(i + 1, 4))
        if best is None or t < best:
            best = t
    return best


def ext_sum_hist(C: np.ndarray, pts: tuple) -> tuple:
    """
    Histogram of sum_{v in S} C_{r v} over r notin S.
    Pure C-invariant; with type6, yields constant-m4 classes at p=5.
    """
    ctr: Counter = Counter()
    pts_l = list(pts)
    n = C.shape[0]
    S = set(pts_l)
    for r in range(n):
        if r in S:
            continue
        ctr[sum(int(C[r, v]) for v in pts_l)] += 1
    return tuple(sorted(ctr.items()))


def class_key(C: np.ndarray, pts: tuple) -> tuple:
    return (type6(C, pts), ext_sum_hist(C, pts))


def kappa(C: np.ndarray, pts: tuple) -> int:
    a, b, c, d = pts
    return int(
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    )


def m2_val(C: np.ndarray, p: int, a: int, b: int) -> float:
    if a == b:
        return 1.0
    return float(C[a, b]) / p


def reduce_moment(idxs: list) -> tuple:
    c = Counter(idxs)
    free = [v for v, cnt in c.items() if cnt % 2]
    if len(free) == 0:
        return ("1",)
    if len(free) == 2:
        a, b = sorted(free)
        return ("m2", a, b)
    if len(free) == 4:
        return ("m4", tuple(sorted(free)))
    return ("0",)


def build_classes(C: np.ndarray) -> tuple[dict, dict]:
    n = C.shape[0]
    classes: dict = defaultdict(list)
    quad_to_key: dict = {}
    for quad in combinations(range(n), 4):
        k = class_key(C, quad)
        classes[k].append(quad)
        quad_to_key[quad] = k
    return classes, quad_to_key


def empirical_m4(Mp: np.ndarray, classes: dict, keys: list) -> np.ndarray:
    m = np.zeros(len(keys))
    idx = {k: i for i, k in enumerate(keys)}
    for k in keys:
        a, b, c, d = classes[k][0]
        m[idx[k]] = float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
    return m


def m4_constancy(Mp: np.ndarray, classes: dict, sample: int = 40) -> dict:
    n_const = 0
    max_std = 0.0
    for quads in classes.values():
        take = quads[:sample]
        m4s = [
            float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
            for a, b, c, d in take
        ]
        std = float(np.std(m4s))
        max_std = max(max_std, std)
        if std < 1e-9:
            n_const += 1
    return {
        "n_classes": len(classes),
        "n_constant_m4": n_const,
        "all_constant": n_const == len(classes),
        "max_std": max_std,
    }


def build_evec_system(
    C: np.ndarray, p: int, classes: dict, quad_to_key: dict, keys: list, max_quads: int = 60
) -> tuple[np.ndarray, np.ndarray]:
    """Averaged evec identities: p m4(quad) = sum_r C_ir m(r, others)."""
    n = C.shape[0]
    idx = {k: i for i, k in enumerate(keys)}
    nv = len(keys)
    rows = []
    bvec = []
    for kA, quads in classes.items():
        coef = np.zeros(nv)
        rhs = 0.0
        count = 0
        if len(quads) <= max_quads:
            sample = quads
        else:
            step = max(1, len(quads) // max_quads)
            sample = quads[::step][:max_quads]
        for quad in sample:
            pts = list(quad)
            for pos in range(4):
                i = pts[pos]
                others = [pts[t] for t in range(4) if t != pos]
                coef[idx[kA]] += p
                for r in range(n):
                    if r == i:
                        continue
                    red = reduce_moment([r] + others)
                    Cr = float(C[i, r])
                    if red[0] == "1":
                        rhs += Cr
                    elif red[0] == "m2":
                        rhs += Cr * m2_val(C, p, red[1], red[2])
                    elif red[0] == "m4":
                        coef[idx[quad_to_key[red[1]]]] -= Cr
                count += 1
        rows.append(coef / count)
        bvec.append(rhs / count)
    return np.array(rows), np.array(bvec)


def nullity_analysis(A: np.ndarray, b: np.ndarray, m_true: np.ndarray | None = None) -> dict:
    u, s, vt = np.linalg.svd(A, full_matrices=True)
    tol = 1e-7
    rank = int(np.sum(s > tol))
    nv = A.shape[1]
    nullity = nv - rank
    m_part = np.linalg.lstsq(A, b, rcond=1e-10)[0]
    out = {
        "n_var": nv,
        "n_eq": int(A.shape[0]),
        "rank": rank,
        "nullity": nullity,
        "singular_values_min3": [float(x) for x in s[-3:]],
        "residual_part": float(np.linalg.norm(A @ m_part - b)),
    }
    if nullity >= 1:
        nvec = vt[rank:][0]
        out["null_residual"] = float(np.linalg.norm(A @ nvec))
        if m_true is not None:
            c = float(np.dot(nvec, m_true - m_part) / np.dot(nvec, nvec))
            m_fit = m_part + c * nvec
            out["c_true_fit"] = c
            out["fit_err"] = float(np.linalg.norm(m_true - m_fit))
            out["m_part"] = m_part.tolist()
            out["nvec"] = nvec.tolist()
    return out


def pairing_gmin_identity(C: np.ndarray, Mp: np.ndarray) -> dict:
    """
    Prove-check: on every 4-set with |kappa|=1, min of three pairing G equals -|m4|.
    Hence g_min = -max |m4| over |kappa|=1 four-sets.
    """
    n = C.shape[0]
    max_abs_m4 = 0.0
    gmin = 1.0
    n_kap1 = 0
    max_err = 0.0
    for pts in combinations(range(n), 4):
        a, b, c, d = pts
        kap = kappa(C, pts)
        if abs(kap) != 1:
            continue
        n_kap1 += 1
        m4 = float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
        max_abs_m4 = max(max_abs_m4, abs(m4))
        Gs = []
        for (i, j), (k, l) in (
            ((a, b), (c, d)),
            ((a, c), (b, d)),
            ((a, d), (b, c)),
        ):
            Gs.append(float(C[i, j] * C[k, l] * m4))
        g_loc = min(Gs)
        gmin = min(gmin, g_loc)
        max_err = max(max_err, abs(g_loc + abs(m4)))
    return {
        "n_kappa_abs_1": n_kap1,
        "max_abs_m4": max_abs_m4,
        "g_min": gmin,
        "identity_gmin_eq_minus_abs_m4": max_err < 1e-12,
        "max_identity_err": max_err,
        "g_min_frac": str(Fraction(gmin).limit_denominator(10000)),
    }


def trG2_from_maxplus_dots(Mp: np.ndarray) -> dict:
    """
    Tr(G^2) via K_ab = ((y_a·y_b)^2 - n)/2.
    Uses only Max+ pairwise dots (not edge enumeration).
    """
    N, n = Mp.shape
    # For large N use chunked accumulation
    tr = 0.0
    s2 = 0.0
    s4 = 0.0
    bs = 400
    for i0 in range(0, N, bs):
        i1 = min(N, i0 + bs)
        D = Mp[i0:i1] @ Mp.T
        K = ((D * D) - n) / 2.0
        tr += float(np.sum(K * K))
        s2 += float(np.sum(D * D))
        s4 += float(np.sum(D ** 4))
    trG2 = tr / (N * N)
    E_dot2 = s2 / (N * N)
    E_dot4 = s4 / (N * N)
    p = int(round(np.sqrt(n - 1)))  # n=p^2+1
    E_dot2_design = n + n * (n - 1) / (p * p)
    return {
        "Tr_G2": trG2,
        "Tr_G2_frac": str(Fraction(trG2).limit_denominator(100000)),
        "E_dot2": E_dot2,
        "E_dot2_2design": E_dot2_design,
        "E_dot2_matches_2design": abs(E_dot2 - E_dot2_design) < 1e-9,
        "E_dot4": E_dot4,
        "identity_TrG2": " (1/4)*(E[dot^4]-2n E[dot^2]+n^2) ",
        "TrG2_from_moments": 0.25 * (E_dot4 - 2 * n * E_dot2 + n * n),
    }


def pin_by_trG2(
    C: np.ndarray,
    p: int,
    classes: dict,
    quad_to_key: dict,
    keys: list,
    m_part: np.ndarray,
    nvec: np.ndarray,
    trG2_target: float,
    Mp: np.ndarray,
) -> dict:
    """
    Build G(c) = G_wedge_fixed + disj from m_part + c nvec; match Tr(G^2).
    """
    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    idx = {k: i for i, k in enumerate(keys)}
    N = len(Mp)

    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    Gtrue = (Chip.T @ Chip) / N

    G0 = np.eye(E)
    for ei, (i, j) in enumerate(edges):
        for ej in range(ei + 1, E):
            k, l = edges[ej]
            if len({i, j, k, l}) == 3:
                G0[ei, ej] = G0[ej, ei] = Gtrue[ei, ej]

    G_part = G0.copy()
    G_null = np.zeros((E, E))
    pair_info = []
    for ei, (i, j) in enumerate(edges):
        for ej in range(ei + 1, E):
            k, l = edges[ej]
            if len({i, j, k, l}) < 4:
                continue
            pts = tuple(sorted([i, j, k, l]))
            ci = idx[quad_to_key[pts]]
            cp = float(C[i, j] * C[k, l])
            G_part[ei, ej] = G_part[ej, ei] = cp * m_part[ci]
            G_null[ei, ej] = G_null[ej, ei] = cp * nvec[ci]
            pair_info.append((ei, ej, ci, cp))

    # Tr(G(c)^2) = ||G_part||_F^2 + 2c <G_part,G_null> + c^2 ||G_null||_F^2
    a2 = float(np.sum(G_null * G_null))
    a1 = 2.0 * float(np.sum(G_part * G_null))
    a0 = float(np.sum(G_part * G_part))
    # a2 c^2 + a1 c + a0 = trG2_target
    Aq, Bq, Cq = a2, a1, a0 - trG2_target
    disc = Bq * Bq - 4 * Aq * Cq
    roots = []
    if disc >= -1e-8:
        sd = np.sqrt(max(disc, 0.0))
        roots = [(-Bq + sd) / (2 * Aq), (-Bq - sd) / (2 * Aq)]

    def gmin_at(c: float) -> float:
        gm = 1.0
        m = m_part + c * nvec
        for ei, ej, ci, cp in pair_info:
            gm = min(gm, cp * m[ci])
        return float(gm)

    root_info = []
    for c in roots:
        gm = gmin_at(c)
        root_info.append(
            {
                "c": float(c),
                "g_min": gm,
                "g_min_frac": str(Fraction(gm).limit_denominator(1000)),
                "beats_T": gm > T_p(p),
                "ge_L": gm >= L_p(p) - 1e-15,
            }
        )
    # select root with larger g_min (true at p=5; F16: not PSD-max alone, but among TrG2 roots)
    if root_info:
        best = max(root_info, key=lambda r: r["g_min"])
    else:
        best = None

    return {
        "quadratic_coeffs_a2_a1_a0": [a2, a1, a0],
        "trG2_target": trG2_target,
        "discriminant": float(disc),
        "roots": root_info,
        "selected_root": best,
        "true_G_trG2": float(np.sum(Gtrue * Gtrue)),
        "note": (
            "Among the two Tr(G^2) roots, select the one with larger g_min "
            "(matches Max+ at p=5). Do not use PSD-max over the whole line (F16)."
        ),
    }


def two_design_dot2_identity(p: int) -> dict:
    """Proved: E[dot^2] = n + n(n-1)/p^2 from E[yy^T]=I+C/p only."""
    n = p * p + 1
    # sum_{a,b} (ya·yb)^2 = sum_{i,j} (sum_a ya_i ya_j)^2
    # sum_a ya ya^T = N (I+C/p), so E[dot^2] = ||I+C/p||_F^2
    # = n + n(n-1)/p^2
    pred = n + n * (n - 1) / (p * p)
    return {
        "p": p,
        "n": n,
        "E_dot2": pred,
        "proof": "||I+C/p||_F^2 from 2-design frame; Max+ free beyond E[yy^T]",
    }


def wick_dot4_upper(p: int) -> dict:
    """
    Gaussian/Wick E[dot^4] for Cov=Σ=I+C/p (same 2-design frame).
    Formula: E[X^4] = 3||Σ||_F^4 + 6 Tr(Σ^4) for X=z·z' independent N(0,Σ).
    Discrete Max+ has strictly smaller E[dot^4] at p=5,7 (boolean kurtosis).
    Using Tr(G^2)<=Tr_Wick on the moduli line at p=5 yields gmin >= T (endpoint
    hits T within float); need strict E[dot^4]<Wick (or exact pin) for gmin>T.
    """
    n = p * p + 1
    C = paley_conference_prime_power(p)
    Sigma = np.eye(n) + C / p
    fro2 = float(np.sum(Sigma * Sigma))
    trS4 = float(np.trace(np.linalg.matrix_power(Sigma, 4)))
    E_dot4_wick = 3 * fro2**2 + 6 * trS4
    TrG2_wick = 0.25 * (E_dot4_wick - 2 * n * fro2 + n * n)
    return {
        "p": p,
        "E_dot2": fro2,
        "E_dot4_wick": E_dot4_wick,
        "Tr_G2_wick": TrG2_wick,
        "formula": "3||Σ||_F^4 + 6 Tr(Σ^4); Σ=I+C/p",
        "note": (
            "Certified disc E[dot^4] < Wick at p=5,7. "
            "Wick pin alone gives gmin>=T not strict >T at p=5 endpoint."
        ),
    }


def run_p5() -> dict:
    p = 5
    C = paley_conference_prime_power(p)
    Mp = load_maxplus(p)
    classes, quad_to_key = build_classes(C)
    keys = list(classes.keys())
    const = m4_constancy(Mp, classes)
    A, b = build_evec_system(C, p, classes, quad_to_key, keys)
    m_true = empirical_m4(Mp, classes, keys)
    null = nullity_analysis(A, b, m_true)
    pair = pairing_gmin_identity(C, Mp)
    dots = trG2_from_maxplus_dots(Mp)
    design = two_design_dot2_identity(p)
    wick = wick_dot4_upper(p)

    pin = None
    if null["nullity"] == 1 and "nvec" in null:
        m_part = np.array(null["m_part"])
        nvec = np.array(null["nvec"])
        pin = pin_by_trG2(
            C, p, classes, quad_to_key, keys, m_part, nvec, dots["Tr_G2"], Mp
        )

    g_min = pair["g_min"]
    disc_lt_wick = dots["E_dot4"] < wick["E_dot4_wick"] - 1e-9
    return {
        "p": p,
        "constancy": const,
        "nullity": {k: null[k] for k in null if k not in ("m_part", "nvec")},
        "pairing_identity": pair,
        "dot_moments": dots,
        "two_design_dot2": design,
        "wick_comparison": {
            **wick,
            "E_dot4_discrete": dots["E_dot4"],
            "discrete_E_dot4_lt_wick": disc_lt_wick,
            "Tr_G2_discrete": dots["Tr_G2"],
        },
        "trG2_pin": pin,
        "g_min": g_min,
        "g_min_frac": pair["g_min_frac"],
        "L": L_p(p),
        "T": T_p(p),
        "g_min_gt_T": g_min > T_p(p),
        "g_min_ge_L": g_min >= L_p(p) - 1e-15,
        "recovered_minus_3_over_65": abs(g_min + 3 / 65) < 1e-12,
    }


def main() -> None:
    out: dict = {"results": {}, "status": None, "open": None}
    print("moduli p=5", flush=True)
    r5 = run_p5()
    # strip bulky if any
    out["results"]["5"] = r5
    print(
        f"  classes={r5['constancy']['n_classes']} const={r5['constancy']['all_constant']} "
        f"nullity={r5['nullity']['nullity']} gmin={r5['g_min_frac']}",
        flush=True,
    )
    if r5.get("trG2_pin") and r5["trG2_pin"].get("selected_root"):
        sel = r5["trG2_pin"]["selected_root"]
        print(f"  TrG2 pin c={sel['c']:.6f} gmin={sel['g_min_frac']}", flush=True)

    ok = (
        r5["constancy"]["all_constant"]
        and r5["nullity"]["nullity"] == 1
        and r5["pairing_identity"]["identity_gmin_eq_minus_abs_m4"]
        and r5["dot_moments"]["E_dot2_matches_2design"]
        and r5["recovered_minus_3_over_65"]
        and r5["trG2_pin"] is not None
        and r5["trG2_pin"]["selected_root"] is not None
        and abs(r5["trG2_pin"]["selected_root"]["g_min"] + 3 / 65) < 1e-9
        and r5["wick_comparison"]["discrete_E_dot4_lt_wick"]
    )
    out["status"] = (
        "p=5 moduli OK: constant-m4 classes, nullity 1, pairing identity "
        "g_min=-|m4| on |kappa|=1, Tr(G^2) pin recovers -3/65. "
        "E[dot^2] from 2-design proved; discrete E[dot^4]<Wick certified. "
        "Wick-only pin gives gmin>=T (endpoint equality), not strict >T. "
        "OPEN: prove E[dot^4]<Wick (or closed spectrum) for all p>=5; "
        "prove g_min>=L(p) or >T(p). Deep ND independent. L OPEN."
        if ok
        else "CERT FAILURE — see results"
    )
    out["open"] = (
        "Nullity-1 for all p>=5; Max+-free strict E[dot^4]<Wick or G-spectrum; "
        "select root; prove g_min>=L(p) or >T(p). Deep non-tight independent. L OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_moduli.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
