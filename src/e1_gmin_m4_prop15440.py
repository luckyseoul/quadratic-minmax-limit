#!/usr/bin/env python3
"""
Prop 15.440 — Lemma D 2-plane is the geometric Fejer triple, not A5.

The locked-edge amplitudes
    Geo(c,d) = ( F(c,s0)F(d,s1),
                 F(c−d,s0)F(−d,s2),
                 F(d−c,s1)F(−c,s2) )
sum to 0 whenever s0+s1+s2 ≡ −2, by the rational identity
    1/((u−1)(v−1))+1/((u−v)(1−v))+1/((v−u)(1−u))=0.
A5-free (c1,c2,c3)=(1,2,3) does **not** lie on the plane.
phi_F / e1 / L not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.439 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (0-polynomial; Max+-free).
  For u≠v, u,v≠1,
      1/((u−1)(v−1))+1/((u−v)(1−v))+1/((v−u)(1−u))=0.
  Clearing (u−1)(v−1)(u−v) gives (u−v)−(u−1)+(v−1)=0.
  Fail: drop the third term (then (u−v)−(u−1)≢0).  ∎

Theorem B — PROVED (A + Fejer; certified p=5,7,11,13,17,19).
  After cancelling the common lock character, Geo(c,d) is a
  nonzero scalar times the reciprocal-product triple of A, hence
  Geo(c,d) lies on {x+y+z=0} for every c≠d in F_p^* with
  c−d≠0 and s0+s1+s2≡−2.  Fail: shift the lock by 2
  (s2 ↦ s2+2 at p=5, (c,d)=(1,2) has |sum|=100).  ∎

Theorem C — PROVED as a negative (certified p=5,7,11).
  A5-free (1,2,3) under the same lock has sum ≠0, so it is not
  a permutation of any Geo(c,d).  A3_PROOF §9's claim that
  convolution puts the A5 3-parameter vector on the plane is
  false as written.  Fail: claim |∑ A5(1,2,3)|<10^{−6}.  ∎

Theorem D — PROVED (B + two Geo; certified p=5,7,11).
  Geo(1,2) and Geo(1,3) are not ℂ-parallel and both lie on
  the plane, so they span {x+y+z=0}.  Fail: claim the 01-minor
  vanishes.  Live locked-triple edge products already sum to 0
  (15.276 certify_construction); that check is now the Geo
  identity, not a census of rank-2.  ∎

============================================================================
Backend: Fraction 0-polynomial + Fejer products.  Writes
evidence/e1_gmin_m4_prop15440.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15276 import (  # noqa: E402
    amplitude_vector_A5,
    fejer_F,
    locked_phases,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15440_catalog.json"


def geo_amplitude_vector(
    p: int, c: int, d: int, s0: int, s1: int, s2: int, lam: int
) -> np.ndarray:
    """Geometric Fejer 3-vector on a locked triple."""
    return np.array(
        [
            fejer_F(p, c, s0, lam) * fejer_F(p, d, s1, lam),
            fejer_F(p, (c - d) % p, s0, lam) * fejer_F(p, (-d) % p, s2, lam),
            fejer_F(p, (d - c) % p, s1, lam) * fejer_F(p, (-c) % p, s2, lam),
        ],
        dtype=np.complex128,
    )


def reciprocal_triple(u: complex, v: complex) -> tuple[complex, complex, complex]:
    return (
        1 / ((u - 1) * (v - 1)),
        1 / ((u - v) * (1 - v)),
        1 / ((v - u) * (1 - u)),
    )


def cleared_numerator(u: int, v: int) -> int:
    """(u−v)−(u−1)+(v−1) over Z, the cleared reciprocal identity."""
    return (u - v) - (u - 1) + (v - 1)


def prove_A() -> dict:
    ok = True
    for u, v in ((2, 3), (3, 5), (4, 7), (-2, 5), (8, -3)):
        if cleared_numerator(u, v) != 0:
            ok = False
        drop3 = (u - v) - (u - 1)
        if drop3 == 0:
            ok = False
    # random integer pairs, avoid u=v
    for u in range(-6, 7):
        for v in range(-6, 7):
            if u == v:
                continue
            if cleared_numerator(u, v) != 0:
                ok = False
    return {
        "proved": bool(ok),
        "drop_third_at_23": (2 - 3) - (2 - 1),
        "theorem": (
            "Reciprocal 3-term is the 0-polynomial after clearing. "
            "Fail: drop the third term."
        ),
    }


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    return all(n % d for d in range(2, int(n**0.5) + 1))


def prove_B() -> dict:
    ok = True
    lock_fail_hit = False
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        if not _is_prime(p):
            continue
        worst = 0.0
        n = 0
        n_ok = 0
        for lam, s0, s1 in ((1, 0, 0), (1, 0, 1), (2, 1, 3)):
            if lam % p == 0:
                continue
            s2 = locked_phases(s0, s1, p)[2]
            for c in range(1, p):
                for d in range(1, p):
                    if c == d or (c - d) % p == 0:
                        continue
                    g = geo_amplitude_vector(p, c, d, s0, s1, s2, lam)
                    n += 1
                    err = float(abs(g.sum()))
                    worst = max(worst, err)
                    if err < 1e-8:
                        n_ok += 1
                    else:
                        ok = False
        rows[str(p)] = {"n_ok": n_ok, "n": n, "worst": worst}
        if n_ok != n:
            ok = False
    # fail-when-wrong: s2 ↦ s2+2 (s2+1 still sums to 0)
    s0, s1 = 0, 1
    s2_lock = locked_phases(s0, s1, 5)[2]
    s2_bad = (s2_lock + 2) % 5
    g_bad = geo_amplitude_vector(5, 1, 2, s0, s1, s2_bad, 1)
    if abs(g_bad.sum()) < 1.0:
        lock_fail_hit = True
        ok = False
    if abs(abs(g_bad.sum()) - 100.0) > 1e-6:
        ok = False
    return {
        "proved": bool(ok) and not lock_fail_hit,
        "lock_drop_sum_abs": float(abs(g_bad.sum())),
        "by_p": rows,
        "theorem": (
            "Geo(c,d) sums to 0 on the lock. Fail: s2 ↦ s2+2 at p=5 "
            "(|sum|=100)."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        s = locked_phases(0, 1, p)
        a5 = amplitude_vector_A5(p, 1, 2, 3, *s, 1)
        sm = float(abs(a5.sum()))
        rows[str(p)] = sm
        if sm < 1e-6:
            ok = False
        # not a permutation of any Geo(c,d)
        matched = False
        for c in range(1, p):
            for d in range(1, p):
                if c == d:
                    continue
                g = geo_amplitude_vector(p, c, d, *s, 1)
                for perm in (
                    (0, 1, 2),
                    (0, 2, 1),
                    (1, 0, 2),
                    (1, 2, 0),
                    (2, 0, 1),
                    (2, 1, 0),
                ):
                    gp = g[list(perm)]
                    if abs(a5 - gp).max() < 1e-6:
                        matched = True
        if matched:
            ok = False
    if rows.get("5", 0.0) < 10.0:
        # live p=5 A5(1,2,3) sum is ~72; a near-zero would be the false close
        ok = False
    return {
        "proved": bool(ok),
        "A5_sum_abs": rows,
        "theorem": (
            "A5-free (1,2,3) is not on the plane and is not a "
            "permutation of any Geo(c,d). Fail: |∑A5|<1e-6."
        ),
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        s = locked_phases(0, 1, p)
        g12 = geo_amplitude_vector(p, 1, 2, *s, 1)
        g13 = geo_amplitude_vector(p, 1, 3, *s, 1)
        if abs(g12.sum()) > 1e-8 or abs(g13.sum()) > 1e-8:
            ok = False
        if np.min(np.abs(g12)) < 1e-10 or np.min(np.abs(g13)) < 1e-10:
            ok = False
        minor = g12[0] * g13[1] - g12[1] * g13[0]
        if abs(minor) < 1e-8:
            ok = False
        rows[str(p)] = {"minor_abs": float(abs(minor))}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Geo(1,2) and Geo(1,3) span the plane. Fail: vanishing minor."
        ),
    }


@lru_cache(maxsize=1)
def plane_identity_proved() -> bool:
    return bool(
        prove_A()["proved"]
        and prove_B()["proved"]
        and prove_C()["proved"]
        and prove_D()["proved"]
    )


def prove_open() -> dict:
    return {
        "proved": False,
        "plane_identity": plane_identity_proved(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Lemma D 2-plane is Geo, not A5. Writeup leftover can "
            "import this unit. phi_F / residual (ii) / Type I unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.440  Lemma D plane is Geo, not A5", flush=True)
    A = prove_A()
    print(f"  A 0-poly: {A['proved']} drop3={A['drop_third_at_23']}", flush=True)
    B = prove_B()
    print(
        f"  B Geo lock: {B['proved']} drop_lock={B['lock_drop_sum_abs']:.3g}",
        flush=True,
    )
    C = prove_C()
    print(f"  C A5 not plane: {C['proved']} |∑A5|={C['A5_sum_abs']}", flush=True)
    D = prove_D()
    print(f"  D span: {D['proved']}", flush=True)
    O = prove_open()
    print(f"  open phi_F={O['phi_F_ge_6']} e1={O['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "D": D["proved"],
                "lock_drop": B["lock_drop_sum_abs"],
                "A5_sum_abs": C["A5_sum_abs"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.440",
        "title": "Lemma D 2-plane is Geo Fejer, not A5-free",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "rational_identity": A["proved"],
            "geo_on_plane": B["proved"],
            "A5_free_not_plane": C["proved"],
            "geo_span_plane": D["proved"],
            "plane_identity": plane_identity_proved(),
            "phi_F_ge_6_proved_general": O["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": O["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "open": O},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15440.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
