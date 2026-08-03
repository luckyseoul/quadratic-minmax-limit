"""Tests for Prop 15.138 — Max+-free y_*; p=5 residual without Max+ census."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np  # used by end-to-end residual test

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15138 import (  # noqa: E402
    YSTAR_CERT,
    construct_ystar,
    main,
    prove_hs_switch,
    prove_norm_circle_construction,
    prove_p5_residual_maxplus_free,
    prove_partial_p7,
)
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)


def test_hs_switch_and_ystar_p5():
    assert prove_hs_switch()["proved"] is True
    r = construct_ystar(5)
    assert r["found"] is True
    assert r["t"] == 0 and r["c"] == 3
    assert r["y_inf"] == 1.0
    assert abs(r["sum"] - 6.0) < 1e-9
    # y is true Max+ evec
    C = paley_conference_prime_power(5).astype(float)
    y = r["y"]
    assert np.allclose(C @ y, 5 * y, atol=1e-5)
    assert np.allclose(np.abs(y), 1.0)
    # not equal to hs (non-hs type)
    hs = halfspace_boolean_vector(5)
    assert not np.allclose(y, hs) and not np.allclose(y, -hs)


def test_norm_circle_certs():
    b = prove_norm_circle_construction([5, 7, 11, 13])
    assert b["proved"] is True
    assert YSTAR_CERT[5]["found"] is True
    assert YSTAR_CERT[7]["t"] == 35 and YSTAR_CERT[7]["c"] == 4
    assert YSTAR_CERT[13]["found"] is False


def test_p5_residual_and_p7_partial():
    c = prove_p5_residual_maxplus_free()
    assert c["proved"] is True
    assert c["match"] is True
    assert room_hyp_over_24(5) == Fraction(1536, 65)
    d = prove_partial_p7()
    assert d["proved"] is True
    assert d["complete"] is False
    assert d["types_from_norm_circles"] == 2
    assert d["types_needed_in_Hplus"] == 5


def test_p5_residual_end_to_end_maxplus_free():
    """Drive shipped construct_ystar + G-quotient to get ∑c_j²=room without Max+ npy."""
    # Local imports keep optional heavy path isolated
    sys.path.insert(0, str(Path("/tmp/grok-goal-9ac8197fcb0f/implementer")))
    from true_aut_orbits_T import (  # noqa: E402
        affine_perm,
        build_orbits,
        close_group,
        field_ops,
    )

    p = 5
    r = construct_ystar(p)
    assert r["found"]
    hs = halfspace_boolean_vector(p)
    ystar = r["y"]
    C = paley_conference_prime_power(p).astype(float)
    mul, add, neg, inv, frob, chi, q = field_ops(p)
    squares = [a for a in range(1, q) if chi(a) == 1]
    gens = [affine_perm(p, 1, b, 0) for b in range(q)]
    gens += [affine_perm(p, a, 0, 0) for a in squares]
    gens.append(affine_perm(p, 1, 0, 1))
    G = close_group(gens, q + 1)
    orbits = sorted(build_orbits(G, q + 1), key=lambda o: (-len(o), min(o)))
    n_orb = len(orbits)
    sizes = np.array([len(o) for o in orbits], dtype=float)
    q2o = {S: i for i, o in enumerate(orbits) for S in o}
    T = np.zeros((n_orb, n_orb))
    n = q + 1
    for i, o in enumerate(orbits):
        acc = np.zeros(n_orb)
        for SS in o:
            Sset = set(SS)
            for v in SS:
                others = [x for x in SS if x != v]
                for rr in range(n):
                    if rr in Sset:
                        continue
                    Sp = tuple(sorted(others + [rr]))
                    acc[q2o[Sp]] += C[v, rr]
        T[i] = acc / len(o)
    Ds = np.sqrt(sizes)
    Tonb = (Ds[:, None] / Ds[None, :]) * T
    Tsym = 0.5 * (Tonb + Tonb.T)
    ew, ev = np.linalg.eigh(Tsym)
    idx4 = [i for i, lam in enumerate(ew) if abs(lam - 4 * p) < 1e-6]
    V = ev[:, idx4]

    def Q(y, vcol):
        f = vcol / Ds
        s = 0.0
        for oi, o in enumerate(orbits):
            if abs(f[oi]) < 1e-15:
                continue
            acc = 0.0
            for SS in o:
                a, b, c, d = SS
                acc += y[a] * y[b] * y[c] * y[d]
            s += f[oi] * acc
        return s

    cvec = np.zeros(len(idx4))
    for j in range(len(idx4)):
        cvec[j] = (3 / 13) * Q(hs, V[:, j]) + (10 / 13) * Q(ystar, V[:, j])
    delta2 = float(np.dot(cvec, cvec))
    assert abs(delta2 - float(room_hyp_over_24(5))) < 1e-8


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["p5_residual_maxplus_free"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["status"]["p5_residual_maxplus_free"] is True
    assert out["status"]["general_residual"] is False
    # No soft-close of L from p=5 alone
    assert "L OPEN" in out["settlement_note"] or out["L_status"] == "OPEN"
    path = ROOT / "evidence" / "e1_gmin_m4_prop15138.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.138"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
