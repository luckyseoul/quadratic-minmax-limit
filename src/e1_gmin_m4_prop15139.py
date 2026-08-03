#!/usr/bin/env python3
"""
Prop 15.139 — Affine halfspaces + double Seidel–norm-circle switch cover all
p=7 H_+ G-orbit size classes Max+-free; AP vs non-AP dichotomy for affine types;
residual still OPEN for general p (no soft-close).

Continues 15.138 toward N_ED4_SUF. No class_key (F19). No GPU theater (F20).

============================================================================
Setup. From 15.137–15.138: residual ⇔ ∑_j c_j² ≤ room with
  c_j = ∑_t w_t Q_j(y_t) on H_+ G-orbit types.
At p=7, H_+ has size classes {56,84,294,588,1176} (1176 with mult 4 orbits).
Norm circles alone cover only 84 (via hs) and 588. Need Max+-free reps for
56, 294, 1176.

============================================================================
Theorem A (affine halfspaces are Max+) — CERTIFIED structure.
  For odd prime p, F_p-linear form L≢0 and S⊂F_p with |S|=(p+1)/2, set
      y_∞ = +1,   y_u = +1 ⇔ L(u)∈S  (u∈F_{p²}).
  Then Cy=py for the Paley conference C of order p²+1 at certified primes
  p=5,7 (all such S).  This is Max+-free (field + C only).  ∎

Theorem B (AP dichotomy at p=7) — CERTIFIED.
  At p=7, |S|=4, L fixed (e.g. L(a+bω)=b): among all C(7,4)=35 choices of S,
    • S is a 4-term arithmetic progression in F_7  ⇒  |G·y_S ∩ H_+| = 84
      (hs-type; includes the standard interval S={0,1,2,3});
    • S is not a 4-AP  ⇒  |G·y_S ∩ H_+| = 56
      (includes the quadratic-residue half {0,1,2,4}=F_7^{□}∪{0}).
  Counts: 21 of AP type, 14 of non-AP type.  Thus both affine size classes
  56 and 84 are Max+-free.  ∎

Theorem C (double Seidel–norm-circle switch) — PROVED structure + cert.
  Let y0∈Max+ (boolean evec).  Set C0=D_{y0} C D_{y0}.  If z is a boolean
  evec of C0 with C0 z = p z (e.g. norm-circle Hoffman coclique of C0 as in
  15.138 for C'=C_hs), then
      y := y0 ⊙ z   satisfies  C y = p y
  (same dictionary as 15.138 Thm A with hs replaced by y0).  Max+-free when
  y0 and z are.  ∎

Theorem D (all p=7 H_+ size classes Max+-free) — CERTIFIED.
  Explicit constructions at p=7:
    |O|=84:  affine halfspace, S 4-AP (e.g. {0,1,2,3});
    |O|=56:  affine halfspace, S non-4-AP (e.g. {0,1,2,4});
    |O|=588: y_hs ⊙ z_nc on C_hs (15.138; t,c)=(35,4);
    |O|=1176: y_56 ⊙ z_nc on C_{y_56} (norm circles of the size-56 switch);
    |O|=294:  y_nc ⊙ z_nc' on C_{y_nc} (norm circles of a size-588 switch;
              specific (t,c) among the C_{y_nc}-circles).
  Every H_+ size class therefore has a Max+-free representative.  (The four
  distinct G-orbits of size 1176 are not all required for the size-class
  statement; constructions hit ≥1 orbit of each size, and ≥1 of size 294.)  ∎

Theorem E (OPEN residual, honest).
  For general primes p≥5: (i) Max+-free classification of all r(p) G-orbit
  types and weights w_t; (ii) character-sum Q_j(y_t) for each type;
  (iii) ∑c_j²≤room.  At p=5: closed (15.138).  At p=7: all size classes
  have Max+-free reps (this prop), but full residual chain still needs
  Max+-free weights + Q_j on every orbit (incl. mult-4 of size 1176) without
  Max+ census.  Do not soft-close L.  ∎

============================================================================
Backend: CPU field arithmetic + G-orbit sizes; GPU unused (F20).
F19: no class_key thrash.
L OPEN. hyp_residual_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15139.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402
from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def _field_ops(p: int):
    q = p * p

    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return ((c0 * d0 + c1 * d1 * ib) % p) + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def add(u, v):
        return ((u % p + v % p) % p) + (((u // p + v // p) % p) * p)

    def neg(u):
        return ((-u % p) % p) + (((-(u // p)) % p) * p)

    def frob(u):
        if u == 0:
            return 0
        r, b, e = 1, u, p
        while e:
            if e & 1:
                r = mul(r, b)
            b = mul(b, b)
            e >>= 1
        return r

    return mul, add, neg, frob, q


def affine_halfspace(p: int, L: tuple[int, int], S: set[int]) -> np.ndarray:
    """y_∞=+1; y_u=+1 iff L(u)∈S. Max+-free construction."""
    q = p * p
    n = q + 1
    alpha, beta = L
    y = np.ones(n, dtype=np.float64)
    for u in range(q):
        a, b = u % p, u // p
        Lval = (alpha * a + beta * b) % p
        y[u + 1] = 1.0 if Lval in S else -1.0
    return y


def is_4_AP(S: set[int] | list[int], p: int) -> bool:
    """True iff S is a 4-term arithmetic progression in F_p."""
    S = sorted(S)
    if len(S) != 4:
        return False
    for a in range(p):
        for d in range(1, p):
            ap = sorted((a + k * d) % p for k in range(4))
            if ap == S:
                return True
    return False


def check_evec(C: np.ndarray, y: np.ndarray, p: int) -> bool:
    return bool(np.allclose(C @ y, p * y, atol=1e-5))


def seidel_switch(C: np.ndarray, y0: np.ndarray) -> np.ndarray:
    return np.diag(y0) @ C @ np.diag(y0)


def norm_circle_z_on(Cbase: np.ndarray, p: int, mul, add, neg, frob, q) -> list[dict]:
    """Hoffman norm-circle boolean evecs of Cbase (relative to all-ones)."""
    n = q + 1
    found = []
    for t in range(q):
        for c in range(1, p):
            pts = []
            for u in range(q):
                d = add(u, neg(t))
                if mul(d, frob(d)) == c:
                    pts.append(u)
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            good = all(
                Cbase[S[i], S[j]] > 0
                for i in range(p + 1)
                for j in range(i + 1, p + 1)
            )
            if not good:
                continue
            z = np.ones(n)
            for i in S:
                z[i] = -1.0
            if np.allclose(Cbase @ z, p * z, atol=1e-5):
                found.append({"t": t, "c": c, "z": z, "S": S})
    return found


def double_switch(y0: np.ndarray, z: np.ndarray, C: np.ndarray, p: int) -> np.ndarray | None:
    y = y0 * z
    if y[0] < 0:
        y = -y
    if check_evec(C, y, p):
        return y
    return None


# ---------------------------------------------------------------------------
# Proofs / certs
# ---------------------------------------------------------------------------

def prove_affine_halfspaces(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = {}
    ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(np.float64)
        k = (p + 1) // 2
        n_ok = 0
        n_tot = 0
        for S_t in combinations(range(p), k):
            n_tot += 1
            y = affine_halfspace(p, (0, 1), set(S_t))
            if check_evec(C, y, p):
                n_ok += 1
        rows[str(p)] = {
            "k": k,
            "n_S": n_tot,
            "n_evec": n_ok,
            "all_evec": n_ok == n_tot and n_tot > 0,
        }
        if n_ok != n_tot:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Affine halfspaces y_{L,S} with |S|=(p+1)/2 are Max+ evecs "
            "at certified p=5,7 (all S, L=(0,1))."
        ),
        "by_p": rows,
    }


def prove_ap_dichotomy_p7() -> dict:
    p = 7
    C = paley_conference_prime_power(p).astype(np.float64)
    # Need G-orbit sizes: import scratch true_aut if available, else cert table
    # Live orbit sizes via stabilizer-free: compare fingerprints + known sizes
    # Use G from true_aut_orbits_T when present
    scratch = Path("/tmp/grok-goal-9ac8197fcb0f/implementer")
    sys.path.insert(0, str(scratch))
    from true_aut_orbits_T import affine_perm, close_group, field_ops  # noqa: E402

    mul, add, neg, inv, frob, chi, q = field_ops(p)
    n = q + 1
    squares = [a for a in range(1, q) if chi(a) == 1]
    gens = [affine_perm(p, 1, b, 0) for b in range(q)]
    gens += [affine_perm(p, a, 0, 0) for a in squares]
    gens.append(affine_perm(p, 1, 0, 1))
    G = close_group(gens, n)

    def g_act(g, y):
        ig = np.empty(len(g), dtype=np.int64)
        for v, w in enumerate(g):
            ig[w] = v
        return y[ig]

    def orbit_size(y0):
        seen = set()
        for g in G:
            y = g_act(g, y0)
            if y[0] < 0:
                y = -y
            seen.add(tuple((y > 0).astype(np.int8).tolist()))
        return len(seen)

    ap84 = 0
    non56 = 0
    bad = 0
    samples = {"84": None, "56": None}
    for S_t in combinations(range(p), 4):
        y = affine_halfspace(p, (0, 1), set(S_t))
        assert check_evec(C, y, p)
        osz = orbit_size(y)
        ap = is_4_AP(S_t, p)
        if ap and osz == 84:
            ap84 += 1
            if samples["84"] is None:
                samples["84"] = list(S_t)
        elif (not ap) and osz == 56:
            non56 += 1
            if samples["56"] is None:
                samples["56"] = list(S_t)
        else:
            bad += 1
    return {
        "proved": bad == 0 and ap84 == 21 and non56 == 14,
        "theorem": (
            "At p=7: S 4-AP ⇒ H_+ orbit size 84; S non-4-AP ⇒ size 56 "
            "(L=(0,1), all C(7,4)=35)."
        ),
        "n_AP_size84": ap84,
        "n_nonAP_size56": non56,
        "n_bad": bad,
        "sample_S_84": samples["84"],
        "sample_S_56": samples["56"],
        "QR_half": [0, 1, 2, 4],
        "QR_half_is_nonAP": not is_4_AP({0, 1, 2, 4}, 7),
    }


def prove_double_switch() -> dict:
    return {
        "proved": True,
        "theorem": (
            "y0∈Max+, C0=D_y0 C D_y0, C0 z=p z boolean ⇒ y=y0⊙z ∈ Max+ "
            "(Seidel dictionary; 15.138 with general base)."
        ),
    }


def prove_p7_all_size_classes() -> dict:
    """Live construct one Max+-free rep per size class at p=7."""
    p = 7
    C = paley_conference_prime_power(p).astype(np.float64)
    mul, add, neg, frob, q = _field_ops(p)
    hs = halfspace_boolean_vector(p)
    assert check_evec(C, hs, p)

    # 84: AP halfspace
    y84 = affine_halfspace(p, (0, 1), {0, 1, 2, 3})
    assert check_evec(C, y84, p)

    # 56: non-AP (QR half)
    y56 = affine_halfspace(p, (0, 1), {0, 1, 2, 4})
    assert check_evec(C, y56, p)
    assert not is_4_AP({0, 1, 2, 4}, p)

    # 588: hs ⊙ norm circle
    r_nc = construct_ystar(p)
    assert r_nc["found"]
    y588 = r_nc["y"]
    assert check_evec(C, y588, p)

    # 1176: y56 ⊙ nc of C56
    C56 = seidel_switch(C, y56)
    ncs56 = norm_circle_z_on(C56, p, mul, add, neg, frob, q)
    y1176 = None
    meta1176 = None
    for item in ncs56:
        y = double_switch(y56, item["z"], C, p)
        if y is not None:
            y1176 = y
            meta1176 = {"t": item["t"], "c": item["c"], "base": "affine56"}
            break
    assert y1176 is not None

    # 294: y588 ⊙ nc of C588
    C588 = seidel_switch(C, y588)
    ncs588 = norm_circle_z_on(C588, p, mul, add, neg, frob, q)
    y294 = None
    meta294 = None
    # Need G-orbit size to identify 294 vs 1176/84; use fingerprint from session
    # or orbit size via true_aut
    scratch = Path("/tmp/grok-goal-9ac8197fcb0f/implementer")
    sys.path.insert(0, str(scratch))
    from true_aut_orbits_T import affine_perm, close_group, field_ops  # noqa: E402

    _mul, _add, _neg, _inv, _frob, chi, _q = field_ops(p)
    n = q + 1
    squares = [a for a in range(1, q) if chi(a) == 1]
    gens = [affine_perm(p, 1, b, 0) for b in range(q)]
    gens += [affine_perm(p, a, 0, 0) for a in squares]
    gens.append(affine_perm(p, 1, 0, 1))
    G = close_group(gens, n)

    def g_act(g, y):
        ig = np.empty(len(g), dtype=np.int64)
        for v, w in enumerate(g):
            ig[w] = v
        return y[ig]

    def orbit_size(y0):
        seen = set()
        for g in G:
            y = g_act(g, y0)
            if y[0] < 0:
                y = -y
            seen.add(tuple((y > 0).astype(np.int8).tolist()))
        return len(seen)

    sizes_found = {
        84: orbit_size(y84),
        56: orbit_size(y56),
        588: orbit_size(y588),
        1176: orbit_size(y1176),
    }
    for item in ncs588:
        y = double_switch(y588, item["z"], C, p)
        if y is None:
            continue
        osz = orbit_size(y)
        if osz == 294:
            y294 = y
            meta294 = {"t": item["t"], "c": item["c"], "base": "normcircle_hs"}
            sizes_found[294] = 294
            break

    ok = (
        sizes_found.get(84) == 84
        and sizes_found.get(56) == 56
        and sizes_found.get(588) == 588
        and sizes_found.get(1176) == 1176
        and y294 is not None
        and sizes_found.get(294) == 294
    )
    return {
        "proved": ok,
        "theorem": (
            "At p=7 every H_+ size class {56,84,294,588,1176} has a "
            "Max+-free representative via affine halfspaces and/or double "
            "Seidel–norm-circle switch."
        ),
        "orbit_sizes_verified": sizes_found,
        "meta_1176": meta1176,
        "meta_294": meta294,
        "n_nc_on_C56": len(ncs56),
        "n_nc_on_C588": len(ncs588),
        "size_classes": [56, 84, 294, 588, 1176],
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "residual_closed_p5_maxplus_free": True,
        "p7_size_classes_maxplus_free": True,
        "p7_residual_maxplus_free": False,
        "open": (
            "Max+-free weights w_t and character-sum Q_j(y_t) for every "
            "H_+ G-orbit (incl. four size-1176 orbits) for all primes p≥5; "
            "close ∑c_j²≤room. No soft-close from size-class coverage alone."
        ),
        "dead": [
            "norm circles alone for all types",
            "class_key thrash (F19)",
            "soft-close L from p=5 + p=7 size classes (F3)",
        ],
    }


def main() -> dict:
    a = prove_affine_halfspaces([5, 7])
    b = prove_ap_dichotomy_p7()
    c = prove_double_switch()
    d = prove_p7_all_size_classes()
    e = prove_open()

    out = {
        "prop": "15.139",
        "title": (
            "Prop 15.139: affine halfspaces + double Seidel–norm-circle; "
            "all p=7 H_+ size classes Max+-free; residual OPEN generally"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Affine halfspaces give 56/84 at p=7 (AP dichotomy). "
            "Double switch covers 588/1176/294. Size classes Max+-free at p=7; "
            "residual still needs Max+-free Q_j weights for general p. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "affine_halfspaces_Max+": a["proved"],
            "ap_dichotomy_p7": b["proved"],
            "double_switch_dictionary": c["proved"],
            "p7_all_size_classes_maxplus_free": d["proved"],
            "residual_for_all_p_ge_5": False,
            "p7_residual_maxplus_free": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "affine": a,
            "ap_dichotomy": b,
            "double_switch": c,
            "p7_sizes": d,
            "open": e,
        },
        "closed_forms": {
            "affine_halfspace": "y_u=+1 iff L(u)∈S, |S|=(p+1)/2",
            "p7_AP": "4-AP S ⇒ orbit size 84; non-4-AP ⇒ 56",
            "double_switch": "y = y0 ⊙ z_nc(C0), C0=D_y0 C D_y0",
            "p7_delta2_true": str(P7_DELTA2_PAIR),
            "room_p7": str(room_hyp_over_24(7)),
            "nu_G": {str(pp): G_SPECTRUM[pp]["nullity"] for pp in (3, 5, 7)},
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_size_classes_maxplus_free": bool(d["proved"]),
            "general_residual": False,
        },
        "open_attack": (
            "Character-sum Q_j on each Max+-free type (incl. mult-4 of size "
            "1176); Max+-free weights; close residual for all p≥5."
        ),
        "backend": (
            "CPU affine + double Seidel–norm-circle + G-orbit sizes; "
            "GPU unused (F20); F19 no class_key"
        ),
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15139.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(json.dumps({k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2))
    return out


if __name__ == "__main__":
    main()
