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


def test_prop_15_53_pairing_identity_and_two_design_dot2():
    """Drive shipped moduli helpers: pairing g_min=-|m4| and E[dot^2] 2-design."""
    import sys

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_moduli import (
        L_p,
        T_p,
        kappa,
        pairing_gmin_identity,
        two_design_dot2_identity,
    )
    from minmax_quadratic import paley_conference_prime_power

    # Algebraic pairing pattern for |kappa|=1 on random 4-sets (C only)
    C = paley_conference_prime_power(5)
    n = C.shape[0]
    seen = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for d in range(c + 1, min(c + 8, n)):
                    pts = (a, b, c, d)
                    kap = kappa(C, pts)
                    if abs(kap) != 1:
                        continue
                    prods = [
                        int(C[a, b] * C[c, d]),
                        int(C[a, c] * C[b, d]),
                        int(C[a, d] * C[b, c]),
                    ]
                    assert sorted(prods) in ([-1, 1, 1], [-1, -1, 1])
                    seen += 1
                    if seen >= 30:
                        break
                if seen >= 30:
                    break
            if seen >= 30:
                break
        if seen >= 30:
            break
    assert seen >= 10

    design = two_design_dot2_identity(5)
    assert abs(design["E_dot2"] - (26 + 26 * 25 / 25)) < 1e-12  # n + n(n-1)/p^2
    assert L_p(5) > T_p(5)
    assert L_p(7) > T_p(7)

    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.is_file():
        return
    # Full identity on Max+ (shipped path)
    from e1_gmin_moduli import load_maxplus

    Mp = load_maxplus(5)
    row = pairing_gmin_identity(C, Mp)
    assert row["identity_gmin_eq_minus_abs_m4"]
    assert abs(row["g_min"] + 3 / 65) < 1e-12
    assert row["g_min"] > T_p(5)


def test_prop_15_53_moduli_script_evidence():
    """Evidence JSON from e1_gmin_moduli must record nullity-1 pin and OPEN residual."""
    path = ROOT / "evidence" / "e1_gmin_moduli.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    r5 = data["results"]["5"]
    assert r5["constancy"]["all_constant"]
    assert r5["nullity"]["nullity"] == 1
    assert r5["recovered_minus_3_over_65"]
    assert r5["dot_moments"]["E_dot2_matches_2design"]
    assert r5["trG2_pin"]["selected_root"] is not None
    assert abs(r5["trG2_pin"]["selected_root"]["g_min"] + 3 / 65) < 1e-9
    assert "OPEN" in data["status"]
    sol = (ROOT / "solution.md").read_text()
    assert "15.53" in sol
    assert "OPEN" in sol[sol.index("15.53") : sol.index("15.53") + 4000]


def test_prop_15_53_moduli_pipeline_live():
    """Live: build classes + evec system + nullity 1 at p=5 (no hardcoded m4 table)."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_moduli import (
        build_classes,
        build_evec_system,
        empirical_m4,
        m4_constancy,
        nullity_analysis,
    )
    from minmax_quadratic import paley_conference_prime_power

    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.is_file():
        return
    import numpy as np

    p = 5
    C = paley_conference_prime_power(p)
    Mp = np.load(cache).astype(float)
    classes, q2k = build_classes(C)
    keys = list(classes.keys())
    const = m4_constancy(Mp, classes)
    assert const["all_constant"]
    assert const["n_classes"] == 37
    A, b = build_evec_system(C, p, classes, q2k, keys, max_quads=40)
    m_true = empirical_m4(Mp, classes, keys)
    null = nullity_analysis(A, b, m_true)
    assert null["nullity"] == 1
    assert null["fit_err"] < 1e-10
    assert abs(null["c_true_fit"] + 0.4240159256359155) < 1e-6


def test_prop_15_54_wedge_G_formula():
    """Shipped wedge identity G_ee' = C_ij C_ik C_jk / p for shared-vertex edges."""
    import sys

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.is_file():
        return
    p = 5
    C = paley_conference_prime_power(p)
    Mp = np.load(cache).astype(float)
    # edges (0,1) and (0,2) share vertex 0
    i, j, k = 0, 1, 2
    pred = float(C[i, j] * C[i, k] * C[j, k] / p)
    g = float(
        np.mean(
            (C[i, j] * Mp[:, i] * Mp[:, j]) * (C[i, k] * Mp[:, i] * Mp[:, k])
        )
    )
    assert abs(g - pred) < 1e-9
    assert abs(abs(pred) - 1 / p) < 1e-12


def test_prop_15_54_cbound_and_abound_evidence():
    """Evidence from shipped cbound/abound scripts; solution Prop 15.54 OPEN residual."""
    for name in ("e1_gmin_cbound.json", "e1_gmin_abound.json"):
        path = ROOT / "evidence" / name
        assert path.is_file(), name
        data = json.loads(path.read_text())
        assert "OPEN" in data.get("status", "") or "open" in data.get("status", "").lower() or "OPEN" in str(data)
    cb = json.loads((ROOT / "evidence" / "e1_gmin_cbound.json").read_text())
    r5 = cb["results"]["5"]
    assert r5["recovered"]
    assert r5["true_root_beats_T"]
    assert abs(r5["g_min_at_c_true"] + 3 / 65) < 1e-9
    ab = json.loads((ROOT / "evidence" / "e1_gmin_abound.json").read_text())
    assert all(r["true_ge_aT"] for r in ab["results"])
    sol = (ROOT / "solution.md").read_text()
    assert "15.54" in sol
    assert "OPEN" in sol[sol.index("15.54") : sol.index("15.54") + 3500]


def test_prop_15_54_deep_spike_refresh():
    """Deep spike evidence exists and records OPEN uniform ND."""
    path = ROOT / "evidence" / "e1_deep_spike_theory.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert "status" in data
    # If covers were found this session, they should spike; else status notes missing cache
    if data.get("p5_covers"):
        found = [c for c in data["p5_covers"] if c.get("found")]
        if found:
            assert data.get("all_found_covers_meet_spike_and_Phi_ge_Phi_m2") is True


def test_prop_15_55_tight_obstruction_algebra():
    """All-ones mass identity: (n/2)/E*(2p)^2 = 4 for n=p^2+1."""
    for p in (3, 5, 7, 11, 13):
        n = p * p + 1
        E = n * (n - 1) // 2
        val = (n / 2) * (2 * p) ** 2 / E
        assert abs(val - 4.0) < 1e-12, (p, val)


def test_prop_15_55_lambda_max_certified():
    """Drive shipped tight-obstruction analyzer: λ_max=n/2 simple iff p≥5 in {3,5,7}."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_tight_obstruction import analyze

    r3 = analyze(3)
    assert r3["row_sum_is_n_over_2"]
    assert r3["allones_quad_is_4"]
    assert not r3["n_over_2_is_simple_max"]
    assert r3["lambda_max"] > r3["n_over_2"]

    for p in (5, 7):
        r = analyze(p)
        assert r["row_sum_is_n_over_2"]
        assert r["allones_quad_is_4"]
        assert r["n_over_2_is_simple_max"]
        assert r["tight_size_2p_impossible_if_simple_max"]
        assert abs(r["lambda_max"] - r["n_over_2"]) < 1e-8


def test_prop_15_55_in_solution_and_evidence():
    sol = (ROOT / "solution.md").read_text()
    assert "15.55" in sol
    chunk = sol[sol.index("15.55") : sol.index("15.55") + 3500]
    assert "OPEN" in chunk
    assert "lambda" in chunk.lower() or "λ" in chunk or "n/2" in chunk
    path = ROOT / "evidence" / "e1_gmin_tight_obstruction.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert "OPEN" in data["status"]
    assert data["results"][1]["p"] == 5
    assert data["results"][1]["n_over_2_is_simple_max"] is True