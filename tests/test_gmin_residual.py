"""Tests for g_min residual (Prop 15.48–15.49 arc): CR classify + uniform LB algebra."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_uniform_lb_beats_threshold_algebra():
    """Candidate g_min >= -(p-2)/(2p^2) strictly beats bi-tight thresh for odd p>2."""
    def thresh(p):
        return Fraction(-(p - 2), p * (2 * p - 1))

    def lb(p):
        return Fraction(-(p - 2), 2 * p * p)

    for p in range(3, 100):
        if p % 2 == 0:
            continue
        assert lb(p) > thresh(p), (p, lb(p), thresh(p))


def test_uniform_lb_holds_for_certified_gmin():
    """Certified g_min at p=5,7 satisfy the candidate LB; p=3 does not (expected)."""
    known = {
        3: Fraction(-1, 3),
        5: Fraction(-3, 65),
        7: Fraction(-109, 2863),
    }

    def lb(p):
        return Fraction(-(p - 2), 2 * p * p)

    assert known[5] >= lb(5)
    assert known[7] >= lb(7)
    assert known[3] < lb(3)  # bi-tight exists at p=3


def test_cr_classify_evidence_and_threshold():
    """CR classification JSON: p=5,7 beat bi-tight threshold; structure present."""
    path = ROOT / "evidence" / "e1_gmin_cr_classify.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["results"]}
    assert 5 in by_p and 7 in by_p
    assert by_p[5]["g_min_gt_threshold"] is True
    assert by_p[7]["g_min_gt_threshold"] is True
    assert by_p[3]["g_min_gt_threshold"] is False
    # constant-m4 |kappa|=1 classes exist
    for p in (5, 7):
        classes = by_p[p]["classes"]
        const_k1 = [
            c
            for c in classes
            if c["m4_constant"] and abs(c["kappa"]) == 1
        ]
        assert len(const_k1) >= 1
        # g_min achieved on some constant-m4 class
        gmins = [c["gmin"] for c in const_k1]
        assert min(gmins) <= by_p[p]["g_min"] + 1e-12


def test_handoff_still_open():
    """L remains OPEN until uniform g_min proved for all p>=5 and deep residual closed."""
    text = (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in text
    assert "15.49" in text or "uniform LB" in text or "2p^2" in text


def test_prop_15_49_in_solution():
    """Prop 15.49 present in solution.md with OPEN residual and LB candidate."""
    sol = (ROOT / "solution.md").read_text()
    assert "15.49" in sol
    idx = sol.index("15.49")
    chunk = sol[idx : idx + 4000]
    assert "OPEN" in chunk
    assert "2p^2" in chunk or "2 p^2" in chunk or "2p^{2}" in chunk or "2p^2" in chunk.replace(" ", "")
    # LB candidate algebra referenced
    assert "threshold" in chunk.lower() or "bi-tight" in chunk.lower()


def test_prop_15_50_cond_mean_shipped():
    """Prop 15.50: drive e1_gmin_cond_mean.gaussian_cond_mean vs Max+ (p=5)."""
    import sys

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_cond_mean import gaussian_cond_mean, sigma_frame
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.is_file():
        # fall back to evidence-only if cache wiped
        path = ROOT / "evidence" / "e1_gmin_cond_mean.json"
        assert path.is_file()
        data = json.loads(path.read_text())
        assert data["certs"][0]["ok"] is True
        return
    C = paley_conference_prime_power(p)
    Sigma = sigma_frame(C, p)
    Mp = np.load(cache).astype(float)
    # fixed pair of coordinates
    i, j = 0, 1
    for a in (1.0, -1.0):
        for b in (1.0, -1.0):
            mask = (Mp[:, i] == a) & (Mp[:, j] == b)
            assert mask.sum() > 0
            emp = Mp[mask].mean(axis=0)
            mu = gaussian_cond_mean(Sigma, i, j, a, b)
            assert np.max(np.abs(emp - mu)) < 1e-9


def test_disj_mean_formula():
    """avg disj G = 1/(p^2-2) for odd primes p."""
    for p in (3, 5, 7, 11, 13):
        n = p * p + 1
        D = (n - 2) * (n - 3) // 2
        s = n // 2 - 1
        assert abs(s / D - 1 / (p * p - 2)) < 1e-15


def test_prop_15_50_in_solution_and_evidence():
    """Prop 15.50 writeup + evidence present; L still OPEN."""
    sol = (ROOT / "solution.md").read_text()
    assert "15.50" in sol
    idx = sol.index("15.50")
    chunk = sol[idx : idx + 3500]
    assert "OPEN" in chunk
    assert "conditional" in chunk.lower() or "Conditional" in chunk
    path = ROOT / "evidence" / "e1_gmin_cond_mean.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert all(c["ok"] for c in data["certs"])


def test_prop_15_51_equiv_form_algebra():
    """a=(1+G*p)/(p+1) converts T-threshold to 1/(2p-1) for all odd p>2."""
    from fractions import Fraction

    for p in range(3, 60):
        if p % 2 == 0:
            continue
        T = Fraction(-(p - 2), p * (2 * p - 1))
        # a at G=T
        a_at_T = (1 + T * p) / (p + 1)
        assert a_at_T == Fraction(1, 2 * p - 1), (p, a_at_T)
        # a at G=L
        L = Fraction(-(p - 2), 2 * p * p)
        a_at_L = (1 + L * p) / (p + 1)
        assert a_at_L == Fraction(p + 2, 2 * p * (p + 1)), (p, a_at_L)


def test_prop_15_51_structure_evidence_and_shipped():
    """Drive e1_gmin_structure.equiv_form_check when Max+ cache present."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    path = ROOT / "evidence" / "e1_gmin_structure.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert "OPEN" in data["status"] or "open" in data["status"].lower()
    for row in data["results"]:
        assert row["equiv"]["a_ge_a_T"] is True
        assert row["loewner"]["loewner_R_ge_lmin_PW"] is True
    # Live check p=5 if cache exists
    cache = Path("/tmp/maxplus_p5.npy")
    if cache.is_file():
        import numpy as np

        from e1_gmin_structure import equiv_form_check

        Mp = np.load(cache).astype(float)
        r = equiv_form_check(5, Mp)
        assert r["a_ge_a_T"] is True
        assert r["G_from_a_gt_T"] is True


def test_prop_15_51_in_solution():
    sol = (ROOT / "solution.md").read_text()
    assert "15.51" in sol
    idx = sol.index("15.51")
    chunk = sol[idx : idx + 3500]
    assert "OPEN" in chunk
    assert "2p-1" in chunk or "2p-1" in chunk.replace(" ", "")


def test_prop_15_52_coordinate_sum_maxplus():
    """|1^T y| = p+1 for every Max+ vector (Prop 15.52); drive halfspace + cache."""
    import sys

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power

    for p in (3, 5, 7):
        h = halfspace_boolean_vector(p)
        assert abs(abs(h.sum()) - (p + 1)) < 1e-9
        # algebraic check from Cy=py on halfspace
        C = paley_conference_prime_power(p)
        assert np.allclose(C @ h, p * h, atol=1e-6)
        # C1 structure: (C1)[0]=p^2, (C1)[1:]=1
        C1 = C @ np.ones(C.shape[0])
        assert abs(C1[0] - p * p) < 1e-6
        assert np.allclose(C1[1:], 1.0, atol=1e-6)
        # identity: s = (p+1) y_inf
        s = float(h.sum())
        assert abs(s - (p + 1) * h[0]) < 1e-9

    cache = Path("/tmp/maxplus_p5.npy")
    if cache.is_file():
        Mp = np.load(cache).astype(float)
        assert np.allclose(np.abs(Mp.sum(axis=1)), 6.0, atol=1e-6)


def test_prop_15_52_in_solution_and_moduli_doc():
    sol = (ROOT / "solution.md").read_text()
    assert "15.52" in sol
    assert "OPEN" in sol[sol.index("15.52") : sol.index("15.52") + 2500]
    assert (ROOT / "evidence" / "E1_GMIN_MODULI.md").is_file()
    text = (ROOT / "evidence" / "E1_GMIN_MODULI.md").read_text()
    assert "nullity" in text.lower() or "Nullity" in text
    assert "OPEN" in text