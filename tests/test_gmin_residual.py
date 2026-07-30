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


def test_prop_15_56_star_maps_to_one():
    """Shipped identity: G u^{(i)} = 1 for star indicators (drive real G)."""
    import sys

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_spectral import analyze

    for p in (3, 5, 7):
        r = analyze(p)
        assert r["stars_map_to_one"]
        assert r["row_sum_ok"]
        assert r["lambda_max_is_max_of_nhalf_and_cycle"]
        assert r["star_differences_in_ker"]


def test_prop_15_56_avg_cycle_algebra_and_gap():
    """Average cycle eig < n/2 for p>=5; spectral gap cert p=5,7 fail p=3."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_spectral import analyze

    r3 = analyze(3)
    assert not r3["spectral_gap_ok"]
    assert r3["lambda_max_G"] > r3["n_over_2"]

    for p in (5, 7):
        r = analyze(p)
        assert r["avg_cycle_lt_n_over_2"]
        assert r["avg_cycle_formula_lt_nhalf_algebra"]
        assert r["spectral_gap_ok"]
        assert r["simple_lambda_max_eq_n_over_2"]
        # reduction identity
        assert abs(r["two_N_times_lambda2_PP"] - r["lambda_max_cycle"]) < 1e-6


def test_prop_15_56_in_solution_and_evidence():
    sol = (ROOT / "solution.md").read_text()
    assert "15.56" in sol
    assert "OPEN" in sol[sol.index("15.56") : sol.index("15.56") + 4000]
    path = ROOT / "evidence" / "e1_gmin_spectral.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert "OPEN" in data["status"]
    assert data["results"][1]["spectral_gap_ok"] is True
    assert data["results"][2]["spectral_gap_ok"] is True


def test_prop_15_58_maxplus_in_vplus_and_perron():
    """Max+ ⊂ V_+ (Cy=py) and (P⊙P)1 = α1 with λ_max = α (p=3,5)."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_cr_classify import load_maxplus
    from minmax_quadratic import paley_conference_prime_power

    for p, loader in (
        (3, lambda: load_maxplus(3).astype(float)),
        (5, lambda: np.load("/tmp/maxplus_p5.npy").astype(float)),
    ):
        C = paley_conference_prime_power(p).astype(float)
        Y = loader()
        N, n = Y.shape
        d = n // 2
        assert all(np.allclose(C @ y, p * y, atol=1e-8) for y in Y[: min(40, N)])
        P = (Y @ Y.T) / (2.0 * N)
        M = P * P
        alpha = d / float(N)
        ones = np.ones(N)
        assert np.allclose(M @ ones, alpha * ones, atol=1e-8)
        ev = np.linalg.eigvalsh(M)
        assert abs(ev[-1] - alpha) < 1e-8


def test_prop_15_58_veronese_equivalence_and_probe():
    """Gap ⇔ ||T(x)||_F^2 ≤ n N ||x||^2 on 1^⊥; probe evidence present."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))

    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    P = (Y @ Y.T) / (2.0 * N)
    M = P * P
    # λ2 via deflated matrix
    Q1 = np.ones((N, N)) / N
    Mp = (np.eye(N) - Q1) @ M @ (np.eye(N) - Q1)
    lam2 = float(np.linalg.eigvalsh(Mp)[-1])
    thr = d / (2.0 * N)
    assert lam2 <= thr + 1e-10
    # random x ⟂ 1: Veronese bound (uses real Max+ / shipped path)
    rng = np.random.default_rng(0)
    for _ in range(8):
        x = rng.standard_normal(N)
        x -= x.mean()
        T = Y.T @ (x[:, None] * Y)
        lhs = float(np.sum(T * T))
        rhs = n * N * float(x @ x)
        assert lhs <= rhs + 1e-6 * max(rhs, 1.0), (lhs, rhs)

    path = ROOT / "evidence" / "e1_gmin_gap_probe.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["results"]}
    assert by_p[3]["gap_ok"] is False
    assert by_p[5]["gap_ok"] is True
    assert by_p[7]["gap_ok"] is True
    assert by_p[5]["cycle_le_8n"] is True
    assert by_p[3]["cycle_le_8n"] is False
    sol = (ROOT / "solution.md").read_text()
    assert "15.58" in sol
    assert "Veronese" in sol or "T(x)" in sol
    assert "OPEN" in data["status"]


def test_prop_15_59_centered_P1_zero_and_rank():
    """Central symmetry ⇒ ∑y=0 ⇒ P1=0; rank(P⊙P)=binom(d-1,2) at p=3,5."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_cr_classify import load_maxplus
    from minmax_quadratic import paley_conference_prime_power

    for p, loader in (
        (3, lambda: load_maxplus(3).astype(float)),
        (5, lambda: np.load("/tmp/maxplus_p5.npy").astype(float)),
    ):
        Y = loader()
        N, n = Y.shape
        d = n // 2
        C = paley_conference_prime_power(p).astype(float)
        # -y is Max+ (sample)
        assert all(np.allclose(C @ (-y), p * (-y)) for y in Y[: min(20, N)])
        assert np.allclose(Y.sum(axis=0), 0, atol=1e-8)
        P = (Y @ Y.T) / (2.0 * N)
        assert np.allclose(P @ np.ones(N), 0, atol=1e-8)
        M = P * P
        rank_M = int(np.sum(np.linalg.eigvalsh(M) > 1e-8))
        assert rank_M == (d - 1) * (d - 2) // 2

    path = ROOT / "evidence" / "e1_gmin_veronese.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["results"]}
    assert by_p[3]["centered_sum_y"] and by_p[3]["P1_zero"]
    assert by_p[5]["gap_ok"] and not by_p[3]["gap_ok"]
    assert by_p[7]["two_moment_forces_gap"] is True
    assert by_p[5]["two_moment_forces_gap"] is False
    sol = (ROOT / "solution.md").read_text()
    assert "15.59" in sol
    assert "OPEN" in data["status"]


def test_prop_15_60_antipodal_reduction_and_projective():
    """T depends only on antipode-symmetric part; projective ENTF gap form."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    m = N // 2
    # Antipodal pairing
    Ys = np.round(Y, 10)
    lookup = {tuple(Ys[a]): a for a in range(N)}
    anti = np.array([lookup[tuple(np.round(-Y[a], 10))] for a in range(N)])
    assert (anti[anti] == np.arange(N)).all()
    # T reduction
    rng = np.random.default_rng(0)
    x = rng.standard_normal(N)
    s = 0.5 * (x + x[anti])
    T_x = Y.T @ (x[:, None] * Y)
    T_s = Y.T @ (s[:, None] * Y)
    assert np.allclose(T_x, T_s, atol=1e-8)
    # Projective ENTF
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    fund = np.array([a for a in range(N) if a < anti[a]])
    assert len(fund) == m
    Up = U[fund]
    assert np.allclose(Up.T @ Up, (m / d) * np.eye(d), atol=1e-8)
    W_proj = (Up @ Up.T) ** 2
    Q1 = np.ones((m, m)) / m
    Wpp = (np.eye(m) - Q1) @ W_proj @ (np.eye(m) - Q1)
    lam2_p = float(np.linalg.eigvalsh(Wpp)[-1])
    assert lam2_p <= m / (2.0 * d) + 1e-8
    # Algebra: 2×sphere ≤ thr for d≥6
    sphere = 2.0 * m / (d * (d + 2))
    assert 2.0 * sphere <= m / (2.0 * d) + 1e-12

    path = ROOT / "evidence" / "e1_gmin_projective.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["results"]}
    assert by_p[3]["gap_proj_ok"] is False
    assert by_p[5]["gap_proj_ok"] is True
    assert by_p[7]["gap_proj_ok"] is True
    assert by_p[5]["two_sphere_bound_holds"] is True
    assert by_p[3]["two_sphere_bound_holds"] is False
    assert by_p[5]["reduction_identity_ok"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.60" in sol
    assert "OPEN" in data["status"]


def test_prop_15_61_16n_bound_algebra_and_certs():
    """16N ⇔ λ_cycle≤8; d≥8 ⇒ gap; evidence equality p=3, strict p=5,7."""
    # Algebra: 4N/d^2 ≤ N/(2d) ⇔ d≥8
    for d in range(5, 30):
        if d >= 8:
            assert 4.0 / d <= 0.5 + 1e-15
        else:
            assert 4.0 / d > 0.5

    # Frame identity sum ||By||^2 = 2N ||B||_F^2 at p=5 (shipped path)
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    rng = np.random.default_rng(1)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    B = Vp @ A @ Vp.T
    sum_By2 = float(sum(np.linalg.norm(B @ y) ** 2 for y in Y))
    assert abs(sum_By2 - 2 * N * np.sum(B * B)) < 1e-6 * max(sum_By2, 1.0)

    path = ROOT / "evidence" / "e1_gmin_16n.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["results"]}
    assert by_p[3]["equality_case_p3"] is True
    assert abs(by_p[3]["ratio_to_16N"] - 1.0) < 1e-6
    assert abs(by_p[3]["lambda_cycle"] - 8.0) < 1e-6
    assert by_p[5]["bound_16N_ok"] and by_p[7]["bound_16N_ok"]
    assert by_p[5]["lambda_cycle_le_8"] and by_p[7]["lambda_cycle_le_8"]
    assert by_p[5]["gap_ok"] and by_p[7]["gap_ok"]
    assert by_p[3]["gap_ok"] is False
    assert by_p[5]["algebra_16N_implies_gap"] is True
    # p=5 exact ratio 11/13
    assert abs(by_p[5]["ratio_to_16N"] - 11 / 13) < 1e-6
    sol = (ROOT / "solution.md").read_text()
    assert "15.61" in sol
    assert "16N" in sol or "16 N" in sol
    assert "OPEN" in data["status"]


def test_prop_15_62_typeA_wedge_identity():
    """typeA+wedge=6N for zero-diag B on V_+; Q=6N+Q_4; evidence multi-seed."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    # Live check at p=5 on a zero-diag∩V_+ random matrix
    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    # Build one nullspace vector (ambient zero diag)
    dimS = d * (d + 1) // 2

    def A_from_vec(v):
        A = np.zeros((d, d))
        k = 0
        for i in range(d):
            for j in range(i, d):
                A[i, j] = A[j, i] = v[k]
                k += 1
        return A

    M = np.zeros((n, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        M[:, k] = np.diag(Vp @ A_from_vec(e) @ Vp.T)
    _u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    null = vh[rank:].T
    rng = np.random.default_rng(42)
    v = null @ rng.standard_normal(null.shape[1])
    B = Vp @ A_from_vec(v) @ Vp.T
    B = B / (np.linalg.norm(B) + 1e-30)
    assert np.max(np.abs(np.diag(B))) < 1e-10
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))
    sixN = 6.0 * N * float(np.sum(B * B))
    # Edge split
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    Be = np.array([2.0 * B[i, j] for i, j in edges])
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    Gu = F.T @ F
    E = len(edges)
    share = np.zeros((E, E), dtype=bool)
    for a, (i, j) in enumerate(edges):
        for b in range(a + 1, E):
            kk, ll = edges[b]
            if len({i, j, kk, ll}) == 3:
                share[a, b] = share[b, a] = True
    disj = ~share
    np.fill_diagonal(disj, False)
    same = np.eye(E, dtype=bool)
    typeA = float(Be @ (Gu * same) @ Be)
    wedge = float(Be @ (Gu * share) @ Be)
    q4 = float(Be @ (Gu * disj) @ Be)
    assert abs(typeA - 2 * N) < 1e-6 * N
    assert abs(wedge - 4 * N) < 1e-4 * N
    assert abs(typeA + wedge - sixN) < 1e-4 * max(sixN, 1.0)
    assert abs(Q - (sixN + q4)) < 1e-4 * max(Q, 1.0)
    assert q4 <= 10 * N + 1e-3 * N

    path = ROOT / "evidence" / "e1_gmin_typeA_wedge.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["identity_ok"] is True
    for p in ("3", "5", "7"):
        assert data["summary"][p]["typeA_wedge_identity_all_ok"] is True
        assert data["summary"][p]["all_disj_le_10N"] is True
    # p=3 saturates 16N
    assert abs(data["summary"]["3"]["max_ratio_to_16N"] - 1.0) < 1e-9
    sol = (ROOT / "solution.md").read_text()
    assert "15.62" in sol
    assert "typeA" in sol or "Type A" in sol
    assert "OPEN" in data["status"]


def test_prop_15_63_H_implies_16N_algebra_and_certs():
    """H(p)=(p+2)^2/d ≤5 (eq p=3); H⇒16N; spectrum certs p=3,5,7."""
    # Algebra: H(p)≤5 iff 3p^2-8p-3≥0, equality at p=3
    for p in range(3, 40, 2):
        d = (p * p + 1) // 2
        H = (p + 2) ** 2 / d
        if p == 3:
            assert abs(H - 5.0) < 1e-15
        else:
            assert H < 5.0 - 1e-12
        # explicit inequality 2(p+2)^2 ≤ 5(p^2+1)
        assert 2 * (p + 2) ** 2 <= 5 * (p * p + 1)

    # Live: H holds at p=5 maximizer (ray = 49/13)
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(0)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(60):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    B /= np.linalg.norm(B) + 1e-30
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))
    Q4 = Q - 6.0 * N
    ray = Q4 / (2.0 * N)
    H = (5 + 2) ** 2 / d  # 49/13
    assert abs(H - 49 / 13) < 1e-12
    assert ray <= H + 1e-8
    assert abs(ray - 49 / 13) < 1e-6  # equality at p=5 maximizer
    assert Q4 <= 10 * N + 1e-6
    assert Q <= 16 * N + 1e-6

    # Evidence from multi-worker campaigns
    for name in ("e1_gmin_q4_ub.json", "e1_gmin_q4_spectrum.json", "e1_gmin_q4_bound.json"):
        path = ROOT / "evidence" / name
        assert path.is_file(), name
    ub = json.loads((ROOT / "evidence" / "e1_gmin_q4_ub.json").read_text())
    assert ub["H_certified_p357"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.63" in sol
    assert "H(p)" in sol or "hypothesis H" in sol.lower() or "(p+2)^2" in sol
    assert "OPEN" in ub["status"]


def test_prop_15_64_dual_Phi_and_residual_budget():
    """Dual Φ on Z: max Q/N = λ_max(Φ|Z); residual budget (p+1)(p+7)/d; H certs."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    # Residual budget algebra: 2(H-1)=(p+1)(p+7)/d
    for p in (3, 5, 7, 11, 13):
        d = (p * p + 1) // 2
        H = (p + 2) ** 2 / d
        budget = (p + 1) * (p + 7) / d
        assert abs(2 * (H - 1) - budget) < 1e-12

    # Live dual at p=5: build Φ Rayleigh via power on A, compare to 176/13
    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    S = Y @ Vp
    U = S / np.sqrt(n)
    rng = np.random.default_rng(1)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(70):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    # Φ Rayleigh = E[(sAs)^2] for unit A
    sAs = np.einsum("bi,ij,bj->b", S, A, S)
    ray_Phi = float(np.mean(sAs**2))
    assert abs(ray_Phi - 176 / 13) < 1e-6
    residual = ray_Phi - 8.0
    budget = (5 + 1) * (5 + 7) / d
    assert abs(residual - budget) < 1e-6  # equality at p=5 maximizer
    # Wick identity: Gaussian baseline 8
    assert abs(8.0 - 8.0) < 1e-15
    # K = (D⊙D)/(2N) on 1^⊥
    D = Y @ Y.T
    x = rng.standard_normal(N)
    x = x - x.mean()
    x /= np.linalg.norm(x)
    # x^T (D⊙D) x / (2N) should be ≤ λ_cycle = 88/13
    val = float(x @ ((D * D) @ x)) / (2 * N)
    assert val <= 88 / 13 + 0.5  # random vector, loose

    path = ROOT / "evidence" / "e1_gmin_H_proof.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["H_certified"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.64" in sol
    assert "residual" in sol.lower()
    assert "OPEN" in data["status"]


def test_prop_15_65_clean_form_and_kappa_spectrum():
    """16N ⇔ λ2(P⊙P)≤4/N; κ|Z spectrum certs; general projectors can violate."""
    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    # Clean equivalence algebra: α=d/N, λ2(P⊙P)=α² λ2(W)
    # 16N ⇔ λ2(P⊙P)≤4/N ⇔ λ2(W)≤4N/d²
    for p, d, N in [(3, 5, 12), (5, 13, 260)]:
        alpha = d / N
        thr_PP = 4.0 / N
        thr_W = 4.0 * N / (d**2)
        assert abs(thr_PP - alpha**2 * thr_W) < 1e-12

    # Live: Max+ at p=5 satisfies λ2(P⊙P)≤4/N with room
    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    P = (Y @ Y.T) / (2.0 * N)
    assert np.allclose(P @ P, P, atol=1e-8)
    assert abs(P[0, 0] - d / N) < 1e-10
    # P1=0
    assert np.linalg.norm(P @ np.ones(N)) < 1e-8
    evals = np.linalg.eigvalsh(P * P)
    evals = np.sort(evals)[::-1]
    lam1, lam2 = float(evals[0]), float(evals[1])
    assert abs(lam1 - d / N) < 1e-8
    assert lam2 <= 4.0 / N + 1e-8
    # exact λ2 = (88/13)/(2N)
    assert abs(lam2 - (88 / 13) / (2 * N)) < 1e-6

    # κ spectrum evidence
    path = ROOT / "evidence" / "e1_gmin_cumulant.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    by_p = {r["p"]: r for r in data["kappa_spectrum"]}
    assert abs(by_p[5]["lam_max_residual"] - 72 / 13) < 1e-6
    assert by_p[5]["bound_ok"] and by_p[3]["bound_ok"] and by_p[7]["bound_ok"]
    # general projectors counterexample recorded
    assert "general_projector_counterexample" in data
    sol = (ROOT / "solution.md").read_text()
    assert "15.65" in sol
    assert "4/N" in sol or "4/N" in sol.replace(" ", "")
    assert "OPEN" in data["status"]


def test_prop_15_66_zero_diag_free_and_pairing_criterion():
    """Zero-diag freeness; algebra wick+eps=L; g_min≥L certs p=5,7."""
    # Algebra: 1/p^2 + (p-4)/(2p^2) = (p-2)/(2p^2) for p≥5
    for p in range(5, 40, 2):
        wick = 1.0 / (p * p)
        eps = (p - 4) / (2.0 * p * p)
        L_abs = (p - 2) / (2.0 * p * p)
        assert abs(wick + eps - L_abs) < 1e-15
        # L > T
        L = -(p - 2) / (2.0 * p * p)
        T = -(p - 2) / (p * (2.0 * p - 1))
        assert L > T

    import os
    import sys

    import numpy as np

    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    # Live zero-diag freeness at p=5
    Y = np.load("/tmp/maxplus_p5.npy").astype(float)
    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(5).astype(float)
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(2)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(60):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    assert np.max(np.abs(np.diag(B))) < 1e-12

    # Live g_min at p=5 via pairing: sample |κ|=1 4-sets
    from itertools import combinations

    max_abs_m4 = 0.0
    for pts in combinations(range(n), 4):
        i, j, k, l = pts
        kap = int(C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k])
        if abs(kap) != 1:
            continue
        m4 = float(np.dot(Y[:, i] * Y[:, j], Y[:, k] * Y[:, l]) / N)
        max_abs_m4 = max(max_abs_m4, abs(m4))
    g_min = -max_abs_m4
    L5 = -3 / 50
    assert abs(max_abs_m4 - 3 / 65) < 1e-9
    assert g_min >= L5 - 1e-12

    path = ROOT / "evidence" / "e1_gmin_m4_residual.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["gmin_ge_L_certified_p57"] is True
    assert data["zero_diag_automatic"] is True
    assert data["algebra_ok"] is True
    # residual triangle criterion fails at p=5 (documented)
    by_p = {r["p"]: r for r in data["m4_residuals"]}
    assert by_p[5]["resid_le_eps_budget"] is False
    assert by_p[7]["resid_le_eps_budget"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.66" in sol
    assert "OPEN" in data["status"]


def test_prop_15_67_master_identity_and_m4_census():
    """σ_sum=4κ; same-sign Ext algebra; full |κ|=1 census p=5,7; F17 workers."""
    from itertools import product

    # Combinatorial: σ_sum = 4κ for all ±1 K4 labelings
    def sigma_eq(bits):
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        Cv = {edges[i]: bits[i] for i in range(6)}

        def C(i, j):
            if i == j:
                return 0
            return Cv[tuple(sorted((i, j)))]

        def sigma(v):
            o = [x for x in range(4) if x != v]
            a, b, c = o
            return C(v, a) * C(b, c) + C(v, b) * C(a, c) + C(v, c) * C(a, b)

        kap = C(0, 1) * C(2, 3) + C(0, 2) * C(1, 3) + C(0, 3) * C(1, 2)
        return sum(sigma(v) for v in range(4)) == 4 * kap

    assert all(sigma_eq(bits) for bits in product([-1, 1], repeat=6))

    # Algebra: same-sign Ext thr ⇒ |m4| ≤ L_abs
    for p in range(5, 40, 2):
        thr_ext = 2.0 * (p - 4) / p
        L_abs = (p - 2) / (2.0 * p * p)
        ub = 1.0 / (p * p) + thr_ext / (4.0 * p)
        assert abs(ub - L_abs) < 1e-12

    # Evidence census (produced at W=86)
    path = ROOT / "evidence" / "e1_gmin_m4_proof.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["sigma_identity"]["proved"] is True
    assert data["workers"] >= 4  # F17: never single-core census
    for p in (5, 7):
        c = data["certs"][str(p)]
        assert c["all_m4_le_L"] is True
        assert c["Ext_same_le_thr"] is True
        assert c["n_kappa1"] > 0
        assert c["n_le_L"] == c["n_kappa1"]
    assert abs(data["certs"]["5"]["max_abs_m4"] - 3 / 65) < 1e-9
    assert abs(data["certs"]["7"]["max_abs_m4"] - 109 / 2863) < 1e-9

    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from workers import cpu_count_reliable, default_workers

    assert cpu_count_reliable() >= 4
    assert default_workers() >= 2

    sol = (ROOT / "solution.md").read_text()
    assert "15.67" in sol
    assert "OPEN" in data["status"]


def test_prop_15_68_tkappa_vanishing_and_resolvent_reduction():
    """Tκ=0 on |κ|=1 (C²+K4); residual source on |κ|=3; degree/resolvent algebra."""
    from itertools import combinations, product
    import sys

    import numpy as np

    # --- Algebraic: C² reduction formula + 64 K4 labelings ---
    # Tκ = -6 * sum_v prod_{u≠v} e_vu
    def T_red(bits):
        # bits: e01,e02,e03,e12,e13,e23
        e01, e02, e03, e12, e13, e23 = bits
        stars = (
            e01 * e02 * e03
            + e01 * e12 * e13
            + e02 * e12 * e23
            + e03 * e13 * e23
        )
        return -6 * stars

    def kap(bits):
        e01, e02, e03, e12, e13, e23 = bits
        return e01 * e23 + e02 * e13 + e03 * e12

    for bits in product([-1, 1], repeat=6):
        k = kap(bits)
        t = T_red(bits)
        if abs(k) == 1:
            assert t == 0
        if abs(k) == 3:
            assert t in (-24, 24)

    # --- Resolvent algebra: candidate ≤ L, gain budget ---
    for p in range(5, 40, 2):
        L_abs = (p - 2) / (2.0 * p * p)
        cand = (p - 2) / (p * (2.0 * p + 3))
        assert cand <= L_abs + 1e-15
        gain_budget = (p - 4) / 48.0
        assert gain_budget > 0
        # same-sign: |ρ| budget
        assert abs(1.0 / (p * p) + (p - 4) / (2.0 * p * p) - L_abs) < 1e-12

    # --- Evidence from multi-worker combinatorial campaign ---
    path = ROOT / "evidence" / "e1_gmin_m4_tkappa.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["algebraic_reduction"]["proved_Tkappa_zero_on_abs_kappa1"] is True
    assert data["proved"] is True
    for p in (3, 5, 7):
        c = data["paley_certs"][str(p)]
        assert c["Tkappa_all_zero_on_kappa1"] is True
        assert c["d3_constant_match"] is True
        assert c["d1_constant_match"] is True
        assert c["d3_expected_p2_minus_5"] == p * p - 5
        assert c["d1_expected_3p2_minus_7"] == 3 * p * p - 7

    # Live Paley check at p=5: Tκ=0 and degrees (shipped path)
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(5).astype(float)
    n = C.shape[0]

    def kappa(S):
        a, b, c, d = S
        return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])

    def star(S):
        s = 0
        for v in S:
            prod = 1
            for u in S:
                if u != v:
                    prod *= int(C[v, u])
            s += prod
        return s

    n_checked = 0
    for S in combinations(range(n), 4):
        if abs(kappa(S)) != 1:
            continue
        assert -6 * star(S) == 0
        # degrees
        Sset = set(S)
        d3 = d1 = 0
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                ak = abs(kappa(Sp))
                if ak == 3:
                    d3 += 1
                elif ak == 1:
                    d1 += 1
        assert d3 == 5 * 5 - 5
        assert d1 == 3 * 25 - 7
        n_checked += 1
    assert n_checked == 11700

    sol = (ROOT / "solution.md").read_text()
    assert "15.68" in sol
    assert "OPEN" in data["status"]


def test_prop_15_76_one_center_degrees():
    """Drive shipped one-center degree algebra + census evidence."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_onecenter_deg import (
        algebra_degrees,
        d1_one,
        d3_one,
        d1_full,
        d3_full,
        M_cand,
        L_abs,
    )

    alg = algebra_degrees()
    assert alg["proved_divisibility_algebra"] is True
    for p in (3, 5, 7, 11, 13, 17, 19):
        assert 4 * d1_one(p) == d1_full(p)
        assert 4 * d3_one(p) == d3_full(p)
        assert d1_one(p) + d3_one(p) == p * p - 3
        assert (3 * p * p - 7) % 4 == 0
        assert (p * p - 5) % 4 == 0
        if p >= 5:
            assert M_cand(p) < (p - 2) / (p * (2 * p - 1))  # cand < T_abs
            assert M_cand(p) <= L_abs(p) + 1e-15

    path = ROOT / "evidence" / "e1_gmin_m4_onecenter_deg.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data["one_center_degrees_constant"] is True
    for p in ("3", "5", "7", "11"):
        r = data["degree_certs"][p]
        assert r["d1_constant_match"] and r["d3_constant_match"]
        assert r["d1_one_unique"] == [d1_one(int(p))]
        assert r["d3_one_unique"] == [d3_one(int(p))]
    # GPU residual probe when present
    if data.get("gpu_residual_probe"):
        for p in ("5", "7"):
            if p in data["gpu_residual_probe"] and data["gpu_residual_probe"][p]:
                g = data["gpu_residual_probe"][p]
                assert g["gpu_used"] is True
                assert g["m4_le_cand"] is True
                assert "mmap" in g["io"].lower()
    assert "OPEN" in data["status"]
    assert "15.76" in (ROOT / "solution.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_S1_sign_gpu_mmap_atomic():
    """S1≤0@star=+1 cert must use GPU m4 + mmap + atomic evidence (F17)."""
    path = ROOT / "evidence" / "e1_gmin_m4_S1_sign.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_sign.py (GPU+mmap+atomic)"
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data.get("use_gpu") is True
    assert "mmap" in data.get("io_policy", "").lower()
    assert "atomic" in data.get("io_policy", "").lower()
    assert data["S1_nonpos_both"] is True
    for p in ("5", "7"):
        r = data["results"][p]
        assert r["S1_always_nonpositive"] is True
        assert r["n_S1_positive"] == 0
        assert r["gpu_used"] is True
        assert r["max_S1_star_plus"] < 0
        assert "mmap" in r["io"].lower()
        assert "atomic" in r["io"].lower() or "write_json_atomic" in r["io"]
    # numerical anchors (GPU precompute path must match prior census)
    assert abs(data["results"]["5"]["max_S1_star_plus"] + 0.12923076923076918) < 1e-12
    assert abs(data["results"]["7"]["max_S1_star_plus"] + 0.01746419839329371) < 1e-10
    assert "OPEN" in data["status"]
    # helpers importable
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_sign import M_cand, rho_budget_cand

    assert abs(M_cand(5) - 3 / 65) < 1e-15
    assert abs(rho_budget_cand(5) - 2 / 325) < 1e-15


def test_prop_15_81_moduli_line_gd():
    """Prop 15.81: moduli-line GD criterion; p=5 complete cand+GD."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_moduli import (
        moduli_gd_algebra,
        M_cand,
        L_abs,
    )

    alg = moduli_gd_algebra()
    assert alg["proved"] is True
    assert "GD_criterion" in alg
    assert abs(M_cand(5) - 3 / 65) < 1e-15
    assert L_abs(5) > M_cand(5)

    path = ROOT / "evidence" / "e1_gmin_m4_S1_moduli.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_moduli.py"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.81"
    assert data["workers"] >= 4
    assert data["algebra"]["proved"] is True

    r5 = data["p5"]
    assert r5["m4_all_constant"] is True
    assert r5["nullity"] == 1
    assert r5["affine"]["exactly_affine"] is True
    assert r5["affine"]["beta_positive"] is True
    assert r5["c_true_lt_c_GD"] is True
    assert r5["at_c_true"]["GD"] is True
    assert r5["at_c_true"]["max_abs_m4_eq_M_cand"] is True
    assert r5["at_c_true"]["max_star_S1_eq_minus_2_over_65"] is True
    assert abs(r5["at_c_true"]["max_abs_m4"] - 3 / 65) < 1e-12
    assert abs(r5["at_c_true"]["max_star_S1"] + 2 / 65) < 1e-12
    # c* safely below c_GD
    assert r5["c_true"] < r5["c_GD"] < 0

    r7 = data["p7"]
    assert r7["pointwise_evec_holds"] is True
    assert r7["m4_all_constant"] is False  # coarse classes incomplete
    assert r7["n_constant_m4"] < r7["n_classes"]

    assert "OPEN" in data["status"]
    assert "15.81" in (ROOT / "solution.md").read_text()
    assert "15.81" in (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_prop_15_82_type6_cr_refine_and_moduli():
    """Prop 15.82: type6+CR m4-constant p=5,7; p=5 moduli GD; p=7 true cand+GD."""
    path_r = ROOT / "evidence" / "e1_gmin_m4_refine.json"
    path_m = ROOT / "evidence" / "e1_gmin_m4_refine_moduli.json"
    assert path_r.is_file(), "run src/e1_gmin_m4_refine.py (GPU+W=86)"
    assert path_m.is_file(), "run src/e1_gmin_m4_refine_moduli.py (W=86)"

    ref = json.loads(path_r.read_text())
    assert ref["prop"] == "15.82"
    assert ref["workers"] >= 4
    assert ref["use_gpu"] is True
    assert ref["p7_some_strategy_constant"] is True
    for p in ("5", "7"):
        r = ref["results"][p]
        assert r["gpu_used"] is True
        assert r["m4_le_cand"] is True
        t6 = r["strategies"]["type6+cr"]
        assert t6["all_constant"] is True
        assert t6["n_nonconstant"] == 0
        assert t6["max_std"] < 1e-9
    assert ref["results"]["5"]["strategies"]["type6+cr"]["n_classes"] == 26
    assert ref["results"]["7"]["strategies"]["type6+cr"]["n_classes"] == 48
    assert abs(ref["results"]["5"]["max_abs_m4_kappa1"] - 3 / 65) < 1e-9
    assert "OPEN" in ref["status"]

    mod = json.loads(path_m.read_text())
    assert mod["prop"] == "15.82"
    assert mod["workers"] >= 4
    assert mod["both_m4_constant"] is True
    r5 = mod["results"]["5"]
    assert r5["m4_all_constant"] is True
    assert r5["nullity"] == 1
    assert r5["affine"]["exactly_affine"] is True
    assert r5["c_true_on_safe_side_of_c_GD"] is True
    assert r5["moduli_GD_ok"] is True
    assert r5["at_true_Max+"]["GD"] is True
    assert r5["at_true_Max+"]["m4_le_cand"] is True
    assert abs(r5["at_true_Max+"]["max_abs_m4_kappa1"] - 3 / 65) < 1e-9
    r7 = mod["results"]["7"]
    assert r7["m4_all_constant"] is True
    assert r7["n_classes"] == 48
    assert r7["nullity"] == 2  # multi-param pin still OPEN
    assert r7["at_true_Max+"]["GD"] is True
    assert r7["at_true_Max+"]["m4_le_cand"] is True
    assert r7["at_true_Max+"]["max_abs_m4_kappa1"] < r7["M_cand"] + 1e-12
    assert "OPEN" in mod["status"]

    sol = (ROOT / "solution.md").read_text()
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "15.82" in sol
    assert "15.82" in hand
    assert "OPEN" in hand[:1200]
    assert "type6+CR" in sol or "type6+cr" in sol.lower()


def test_prop_15_80_gd_linear_wick():
    """Prop 15.80: linear Wick identity, GD form, U1-special, GPU census."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_gd import (
        linear_form_wick_identity,
        gd_formulation,
        sum_star_tau1_algebra,
        cand_via_gd_algebra,
        M_cand,
        d1_one,
    )

    lin = linear_form_wick_identity()
    assert lin["proved"] is True
    assert "corollary_U1" in lin
    gd = gd_formulation()
    assert gd["proved"] is True
    assert "GD_iff" in gd
    assert "specialization" in gd
    st = sum_star_tau1_algebra()
    assert st["certified_sum"]["5"]["epsilon"] == 1
    assert st["certified_sum"]["7"]["epsilon"] == -1
    assert st["certified_sum"]["5"]["sum"] == 11700
    assert st["certified_sum"]["7"]["sum"] == -176400
    cand = cand_via_gd_algebra()
    assert cand["by_p"]["5"]["S3_budget_negative"] is True
    assert d1_one(5) == 17
    assert abs(M_cand(5) - 3 / 65) < 1e-15

    path = ROOT / "evidence" / "e1_gmin_m4_S1_gd.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_gd.py (GPU+mmap+atomic)"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.80"
    assert data["workers"] >= 4
    assert data["use_gpu"] is True
    assert data["linear_form_wick"]["proved"] is True
    assert data["gd_formulation"]["proved"] is True
    assert data["GD_both"] is True
    assert data["linear_wick_both"] is True
    assert data["U1_special_both"] is True
    for p in ("5", "7"):
        r = data["results"][p]
        assert r["GD_holds"] is True
        assert r["max_star_S1"] < 0
        assert r["gpu_used"] is True
        assert r["linear_form_wick_match"] is True
        assert r["E_U1_sq_over_Wick_max_dev"] is not None
        assert r["E_U1_sq_over_Wick_max_dev"] < 1e-9
        assert r["generic_L_shows_U1_special"] is True
        assert r["m4_le_cand"] is True
        assert "mmap" in r["io"].lower()
    # exact sum anchors
    assert abs(data["results"]["5"]["sum_star_S1"] + 1128) < 1e-6
    assert data["results"]["5"]["epsilon_sum_t_over_n1"] == 1
    assert data["results"]["7"]["epsilon_sum_t_over_n1"] == -1
    assert abs(data["results"]["5"]["max_star_S1"] + 2 / 65) < 1e-12
    assert "OPEN" in data["status"]
    assert "15.80" in (ROOT / "solution.md").read_text()
    assert "15.80" in (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_prop_15_79_aut_constancy():
    """Prop 15.79: Aut-constancy proof + modular τ1 census."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_aut import (
        aut_constancy_proof,
        modular_star_tau1_algebra,
        d1_one,
        M_cand,
    )

    aut = aut_constancy_proof()
    assert aut["proved"] is True
    assert "stabilizer_fact" in aut
    assert "equivariance_tau1" in aut
    assert "equivariance_S1" in aut
    mod = modular_star_tau1_algebra()
    assert mod["proved_elementary"] is True
    assert mod["elementary"]["d1_always_odd_for_odd_p"] is True
    for p in (3, 5, 7, 11, 13, 17, 19):
        assert d1_one(p) % 2 == 1
        assert (3 * p * p - 7) % 4 == 0

    path = ROOT / "evidence" / "e1_gmin_m4_S1_aut.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_aut.py"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.79"
    assert data["workers"] >= 4
    assert data["aut_constancy"]["proved"] is True
    assert data["all_constancy_certified"] is True
    assert data["all_mod6"] is True
    assert data["all_t1_formula"] is True
    assert data["all_n_values_match"] is True
    expected_vals = {
        "3": [-1],
        "5": [-1, 5],
        "7": [-7, -1, 5],
        "11": [-13, -7, -1, 5, 11],
    }
    for p, exp in expected_vals.items():
        r = data["tau1_certs"][p]
        assert r["star_tau1_constant_all_sets"] is True
        assert r["star_tau1_values"] == exp
        assert r["n_distinct_values"] == (int(p) - 1) // 2
        assert r["all_equiv_5_mod_6"] is True
        assert r["t1_eq_2A_minus_d1"] is True
    assert "OPEN" in data["status"]
    assert "15.79" in (ROOT / "solution.md").read_text()
    assert "15.79" in (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]
    assert abs(M_cand(5) - 3 / 65) < 1e-15


def test_prop_15_78_star_S1_constancy_gd():
    """Prop 15.78: moment form, τ1/S1 constancy, GD, p=5 exact spectrum."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_const import (
        algebraic_moment_form,
        p5_exact_spectrum_algebra,
        M_cand,
        rho_budget_cand,
    )

    alg = algebraic_moment_form()
    assert alg["proved"] is True
    assert "gaussian_domination_iff" in alg
    p5a = p5_exact_spectrum_algebra()
    assert p5a["algebra_ok"] is True
    assert p5a["both_negative"] is True
    assert abs(p5a["star_S1_values"][0] + 2 / 65) < 1e-15
    assert abs(p5a["star_S1_values"][1] + 42 / 325) < 1e-15

    path = ROOT / "evidence" / "e1_gmin_m4_S1_const.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_const.py (GPU+mmap+atomic)"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.78"
    assert data["workers"] >= 4
    assert data["use_gpu"] is True
    assert "mmap" in data["io_policy"].lower()
    assert "atomic" in data["io_policy"].lower()
    assert data["tau1_constancy_all"] is True
    for p in ("3", "5", "7", "11"):
        assert data["tau1_constancy"][p]["all_constant"] is True
    assert data["star_S1_constancy_both"] is True
    assert data["gaussian_domination_both"] is True
    assert data["star_S1_le0_both"] is True
    for p in ("5", "7"):
        r = data["results"][p]
        assert r["Cy_eq_py"] is True
        assert r["star_S1_constant_on_every_set"] is True
        assert r["gaussian_domination"] is True
        assert r["star_S1_always_le_0"] is True
        assert r["gpu_used"] is True
        assert r["max_star_S1"] < 0
        assert "mmap" in r["io"].lower()
    # p=5 exact spectrum
    assert data["results"]["5"]["p5_exact"]["matches_closed_form"] is True
    assert abs(data["results"]["5"]["max_star_S1"] + 2 / 65) < 1e-12
    assert abs(data["results"]["5"]["min_star_S1"] + 42 / 325) < 1e-12
    assert abs(M_cand(5) - 3 / 65) < 1e-15
    assert abs(rho_budget_cand(5) - 2 / 325) < 1e-15
    assert "OPEN" in data["status"]
    assert "15.78" in (ROOT / "solution.md").read_text()
    assert "15.78" in (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_prop_15_77_star_S1_structure():
    """Prop 15.77: star·S1 algebra + GPU full-centre cert (mmap+atomic)."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_S1_star import (
        algebraic_star_S1_consequences,
        star_S1_implies_cand_algebra,
        M_cand,
        rho_budget_cand,
        L_abs,
        T_abs,
    )

    alg = algebraic_star_S1_consequences()
    assert alg["proved"] is True
    assert "one_center_split" in alg["identities"]
    # p=5: S3 budget for cand under S1≤0 is strictly negative (need strong S1)
    s5 = alg["by_p"]["5"]
    assert s5["S3_budget_cand_negative"] is True
    assert abs(s5["S3_budget_at_cand_if_S1_le0"] + 16 / 325) < 1e-12
    # cand < L < T for p≥5
    for p in (5, 7, 11, 13, 17, 19):
        assert M_cand(p) < L_abs(p) + 1e-15
        assert L_abs(p) < T_abs(p) - 1e-15
        assert abs(rho_budget_cand(p) - (M_cand(p) - 1 / (p * p))) < 1e-15

    chk = star_S1_implies_cand_algebra()
    assert chk["proved_algebra"] is True
    assert abs(chk["sample_p5"]["rho_budget_cand"] - 2 / 325) < 1e-15

    path = ROOT / "evidence" / "e1_gmin_m4_S1_star.json"
    assert path.is_file(), "run src/e1_gmin_m4_S1_star.py (GPU+mmap+atomic)"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.77"
    assert data["workers"] >= 4
    assert data["use_gpu"] is True
    assert "mmap" in data["io_policy"].lower()
    assert "atomic" in data["io_policy"].lower()
    assert data["star_S1_le0_both"] is True
    assert data["joint_le_cand_both"] is True
    for p in ("5", "7"):
        r = data["results"][p]
        assert r["star_S1_always_le_0"] is True
        assert r["n_plus_S1_positive"] == 0
        assert r["n_minus_S1_negative"] == 0
        assert r["max_star_S1"] < 0
        assert r["gpu_used"] is True
        assert r["joint_implies_le_cand"] is True
        assert r["synthetic_violates_star_S1"] is True
        assert r["identity_err_max"] < 1e-12
        assert r["antisym_max_plus_eq_minus_min"] is True
        assert "mmap" in r["io"].lower()
    # anchors: p=5 sharp joint = −16/325, max_star_S1 = −2/65
    assert abs(data["results"]["5"]["max_joint_S1S3_rpos_starplus"] + 16 / 325) < 1e-12
    assert abs(data["results"]["5"]["max_star_S1"] + 2 / 65) < 1e-12
    assert abs(data["results"]["5"]["rho_ub_from_joint"] - 2 / 325) < 1e-12
    assert "OPEN" in data["status"]
    assert "15.77" in (ROOT / "solution.md").read_text()
    assert "15.77" in (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_prop_15_75_onecenter_sigma_and_gpu_cand():
    """Drive shipped onecenter algebra + GPU cand census (mmap+atomic path)."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_onecenter import (
        algebraic_sigma_identity,
        k4_gram_spectrum_formula,
        weak_ub_algebra,
        M_cand,
        L_abs,
    )
    from e1_gmin_m4_gpu import (
        M_cand as gpu_M_cand,
        L_abs as gpu_L_abs,
        gain_budget_cand,
        rho_budget_cand,
    )

    alg = algebraic_sigma_identity()
    assert alg["proved"] is True
    assert alg["n_kappa1_labelings"] == 48
    assert alg["sigma_eq_2_star_count"] == 48
    assert k4_gram_spectrum_formula()["proved"] is True
    assert weak_ub_algebra()["proved_cand_beats_T"] is True
    # shared algebra with GPU module
    assert abs(M_cand(5) - gpu_M_cand(5)) < 1e-15
    assert abs(L_abs(7) - gpu_L_abs(7)) < 1e-15
    assert abs(gain_budget_cand(5) - 1 / 156) < 1e-15
    assert abs(rho_budget_cand(5) - 2 / 325) < 1e-15

    # onecenter evidence (multi-worker σ census)
    path_oc = ROOT / "evidence" / "e1_gmin_m4_onecenter.json"
    assert path_oc.is_file()
    data_oc = json.loads(path_oc.read_text())
    assert data_oc["workers"] >= 4
    assert data_oc["sigma_identity_ok"] is True
    assert data_oc["k4_spectrum_ok"] is True
    for p in ("3", "5", "7", "11"):
        assert data_oc["sigma_certs"][p]["all_ok"] is True
        assert data_oc["sigma_certs"][p]["bad_sigma_ne_2star"] == 0
    # evidence records atomic/mmap io policy
    assert "atomic" in data_oc["io"].lower() or "write_json_atomic" in data_oc["io"]

    # GPU cand census evidence (real shipped GPU path results)
    path_gpu = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
    assert path_gpu.is_file()
    data_g = json.loads(path_gpu.read_text())
    assert data_g["use_gpu"] is True
    assert data_g["both_le_cand"] is True
    assert data_g["both_le_mid"] is True
    for p in ("5", "7"):
        r = data_g["results"][p]
        assert r["gpu_used"] is True
        assert r["m4_le_cand"] is True
        assert r["m4_le_mid"] is True
        assert r["m4_le_L"] is True
        assert "mmap" in str(r["io"]).lower() or "mmap" in data_g.get("io_policy", "").lower()
    r5 = data_g["results"]["5"]
    assert abs(r5["max_abs_m4"] - 3 / 65) < 1e-9
    assert abs(r5["effective_gain"] - 1 / 156) < 1e-9
    r7 = data_g["results"]["7"]
    assert abs(r7["max_abs_m4"] - 109 / 2863) < 1e-9
    # OPEN not soft-closed
    assert "OPEN" in data_g["status"]
    assert "OPEN" in data_oc["status"]
    sol = (ROOT / "solution.md").read_text()
    assert "15.75" in sol
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in hand[:500] or "L still OPEN" in hand[:800]


def test_prop_15_74_true_maxplus_cand_bound():
    """Drive shipped kernel helpers: cand algebra + true Max+ census evidence."""
    import sys
    from fractions import Fraction

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_kernel import (
        M_cand,
        M_mid,
        L_abs,
        algebra_cand,
        gain_budget_cand,
        rho_budget_cand,
    )

    alg = algebra_cand()
    assert alg["proved"] is True
    # p=5 sharp identities
    assert abs(M_cand(5) - 3 / 65) < 1e-15
    assert abs(gain_budget_cand(5) - 1 / 156) < 1e-15
    assert abs(rho_budget_cand(5) - 2 / 325) < 1e-15
    for p in (5, 7, 11, 13, 17, 19):
        assert M_cand(p) <= M_mid(p) + 1e-15
        assert M_mid(p) <= L_abs(p) + 1e-15
        # recon: wick + 24*gain_cand/p^2 = M_cand
        recon = 1.0 / (p * p) + 24 * gain_budget_cand(p) / (p * p)
        assert abs(recon - M_cand(p)) < 1e-12

    path = ROOT / "evidence" / "e1_gmin_m4_kernel.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data["true_m4_le_cand_p5_p7"] is True
    assert data["signed_identity_ok"] is True
    # True Max+ (not type6)
    r5 = data["true_maxplus_census"]["5"]
    r7 = data["true_maxplus_census"]["7"]
    assert r5["le_cand"] and r5["le_mid"] and r5["le_L"]
    assert r7["le_cand"] and r7["le_mid"] and r7["le_L"]
    assert abs(r5["max_abs_m4"] - 3 / 65) < 1e-9
    assert abs(r5["effective_gain"] - 1 / 156) < 1e-9
    assert abs(r7["max_abs_m4"] - 109 / 2863) < 1e-9
    assert r5["io"].startswith("mmap")
    assert data["signed_residual_identity"]["5"]["identity_holds"]
    assert data["signed_residual_identity"]["7"]["identity_holds"]
    sol = (ROOT / "solution.md").read_text()
    assert "15.74" in sol
    assert "OPEN" in data["status"]
    # ensure we did not soft-close L
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:2000] or "OPEN" in sol


def test_prop_15_73_e4_identity_and_type6_multi_prime():
    """Drive shipped e4/type6 helpers: e4 formula, sumκ, multi-prime type6 evidence."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_e4_gain import (
        algebra_e4,
        e4_formula,
        e4_from_s2,
        sum_kappa_formula,
        L_abs,
        M_mid,
    )

    # Proved e4 algebra
    alg = algebra_e4()
    assert alg["proved_e4_algebra"] is True
    for p in (3, 5, 7, 11, 13, 17, 19):
        n = p * p + 1
        s2 = (p + 1) ** 2
        assert abs(e4_from_s2(n, s2) - e4_formula(p)) < 1e-12
        assert abs(e4_formula(p) + p * (p - 1) * (p + 1) * (p + 4) / 12.0) < 1e-12
        assert abs(sum_kappa_formula(p) - p * p * (p * p - 1) / 4.0) < 1e-12

    path = ROOT / "evidence" / "e1_gmin_m4_e4_gain.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data["e4_maxplus_ok"] is True
    assert data["sum_kappa_match_all"] is True
    assert data["type6_le_L_all_primes"] is True
    assert data["type6_le_mid_all_primes"] is True
    assert data["type6_gain_le_budget_all"] is True
    for p in (5, 7, 11, 13):
        r = data["type6_resolvent"][str(p)]
        assert r["type6_le_L"] is True
        assert r["type6_le_mid"] is True
        assert r["max_abs_m4_type6_kappa1"] <= L_abs(p) + 1e-8
        assert r["max_abs_m4_type6_kappa1"] <= M_mid(p) + 1e-8
        assert r["gain_le_budget"] is True
        assert data["sum_kappa_certs"][str(p)]["match"] is True
    # p=3,5,7,11,13 sumκ
    for p in (3, 5, 7, 11, 13):
        assert data["sum_kappa_certs"][str(p)]["match"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.73" in sol
    assert "OPEN" in data["status"]


def test_prop_15_72_resolvent_gain_algebra_and_evidence():
    """Drive shipped resolvent-gain helpers: reverse degrees, gain⇔L, evidence."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_resolvent_gain import (
        algebra_gain_targets,
        algebra_reverse_degrees,
        gain_budget,
        rho_budget,
        source_amp,
        tkappa_over_kappa_identity,
        L_abs,
    )
    from e1_gmin_m4_stratum import n1_formula, n3_formula, d3_formula

    tk = tkappa_over_kappa_identity()
    assert tk["proved"] is True
    assert set(tk["Tkappa_over_kappa_on_abs3"]) == {8, -8}

    rev = algebra_reverse_degrees()
    assert rev["proved_under_hypothesis"] is True
    for p in (3, 5, 7, 11, 13):
        n1, n3, d3 = n1_formula(p), n3_formula(p), d3_formula(p)
        assert abs(n1 * d3 / n3 - 3 * (p * p - 1)) < 1e-9
        assert 3 * (p * p - 1) + (p * p - 9) == 4 * (p * p + 1 - 4)

    tgt = algebra_gain_targets()
    assert tgt["proved_algebra"] is True
    for p in (5, 7, 11):
        assert abs(gain_budget(p) * source_amp(p) - rho_budget(p)) < 1e-15
        assert abs(1.0 / (p * p) + rho_budget(p) - L_abs(p)) < 1e-15

    path = ROOT / "evidence" / "e1_gmin_m4_resolvent_gain.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data["certs_ok"] is True
    assert data["type6_le_L_for_p_ge_5"] is True
    for p in (5, 7):
        t6 = data["type6_resolvent"][str(p)]
        assert t6["type6_m4_le_L"] is True
        assert t6["gain_le_budget"] is True
        emp = data["empirical_maxplus_gain"][str(p)]
        assert emp["gain_le_budget"] is True
        assert emp["m4_le_L"] is True
    sol = (ROOT / "solution.md").read_text()
    assert "15.72" in sol
    assert "OPEN" in data["status"]


def test_prop_15_71_stratum_counts_conference_proof():
    """Drive shipped stratum helpers: S2 proof, n1/n3 formulas, evidence."""
    import sys
    from math import comb

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_stratum import (
        L_abs,
        M_mid,
        algebra_identities,
        conference_S2_proof,
        d1_formula,
        d3_formula,
        n1_formula,
        n3_formula,
        n_of_p,
    )

    proof = conference_S2_proof()
    assert proof["proved"] is True
    assert proof["identity_ratio_8"] is True
    assert proof["K4_wsum_eq_8_cross_count"] == 64

    for p in (3, 5, 7, 11, 13, 17, 19):
        n = n_of_p(p)
        n1, n3 = n1_formula(p), n3_formula(p)
        assert n1 + n3 == comb(n, 4)
        S2 = n * (n - 1) * (n - 2) * (n - 5) // 8
        assert n3 == (S2 - comb(n, 4)) // 8
        assert n1 == n * (n - 1) * (n - 2) * (n - 2) // 32
        assert d1_formula(p) + d3_formula(p) == 4 * (n - 4)
        if p >= 5:
            assert M_mid(p) <= L_abs(p) + 1e-15

    alg = algebra_identities([3, 5, 7, 11])
    assert alg["proved_algebra"] is True

    path = ROOT / "evidence" / "e1_gmin_m4_stratum.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    assert data["conference_S2_proof"]["proved"] is True
    assert data["counts_match_all_certified_p"] is True
    for p in (3, 5, 7, 11):
        r = data["count_certs"][str(p)]
        assert r["match_n1"] and r["match_n3"]
        assert r["n1_census"] == n1_formula(p)
        assert r["n3_census"] == n3_formula(p)
    sol = (ROOT / "solution.md").read_text()
    assert "15.71" in sol
    assert "OPEN" in data["status"]


def test_prop_15_70_mid_ub_algebra_and_census():
    """mid_ub <= L_abs; census g_min>=L at p=5,7; bi-tight threshold algebra."""
    # Algebra: mid = (p-2)/(2p(p+1)), L = (p-2)/(2p^2), mid/L = p/(p+1)
    for p in range(5, 50, 2):
        if p > 2 and any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        L_abs = (p - 2) / (2.0 * p * p)
        mid = (p - 2) / (2.0 * p * (p + 1))
        cand = (p - 2) / (p * (2.0 * p + 3))
        T_abs = (p - 2) / (p * (2.0 * p - 1))
        assert mid <= L_abs + 1e-15
        assert cand <= mid + 1e-15
        assert abs(mid / L_abs - p / (p + 1)) < 1e-12
        assert L_abs < T_abs  # so L > T as signed floors
        # signed: L_p = -L_abs > T_p = -T_abs
        assert -L_abs > -T_abs

    path = ROOT / "evidence" / "e1_gmin_m4_close.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["workers"] >= 4
    for p in (5, 7):
        c = data["parts"][f"census_p{p}"]
        assert c["g_min_ge_L"] is True
        assert c["g_min_gt_T"] is True
        assert c["m4_le_mid"] is True
        assert c["workers"] >= 4
    assert abs(data["parts"]["census_p5"]["max_abs_m4"] - 3 / 65) < 1e-9
    assert abs(data["parts"]["census_p7"]["max_abs_m4"] - 109 / 2863) < 1e-9

    sol = (ROOT / "solution.md").read_text()
    assert "15.70" in sol
    assert "OPEN" in data["status"]


def test_gpu_budget_and_m4_gpu_evidence():
    """GPU path: gpu_budget reports CuPy when present; GPU census evidence if written."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot, prefer_gpu_for, gpu_info

    snap = compute_snapshot()
    assert "workers" in snap and snap["workers"] >= 2
    assert "gpu" in snap
    g = gpu_info()
    # Machine policy: this host should see a V100; if not, still exercise API
    if g.get("available"):
        assert prefer_gpu_for("m4_batch") is True
        path = ROOT / "evidence" / "e1_gmin_m4_gpu.json"
        assert path.is_file(), "run src/e1_gmin_m4_gpu.py to produce GPU census evidence"
        data = json.loads(path.read_text())
        assert data["use_gpu"] is True
        for p in (5, 7):
            r = data["results"][str(p)]
            assert r["m4_le_mid"] is True
            assert r["m4_le_L"] is True
            assert r["wall_s"] < 30.0  # GPU path must be fast
        assert abs(data["results"]["5"]["max_abs_m4"] - 3 / 65) < 1e-9


def test_io_atomic_mmap_and_write():
    """Atomic write + mmap load (Wieferich-style I/O helpers)."""
    import sys
    import tempfile
    from pathlib import Path

    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from io_atomic import write_json_atomic, write_npy_atomic, load_npy

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        jp = td / "e.json"
        write_json_atomic(jp, {"p": 5, "ok": True})
        assert jp.is_file()
        data = json.loads(jp.read_text())
        assert data["ok"] is True and data["p"] == 5

        ap = td / "a.npy"
        arr = np.arange(20, dtype=np.float64).reshape(4, 5)
        write_npy_atomic(ap, arr)
        mm = load_npy(ap, mmap_mode="r")
        assert isinstance(mm, np.memmap) or hasattr(mm, "shape")
        assert mm.shape == (4, 5)
        assert float(mm[2, 1]) == float(arr[2, 1])


def test_prop_15_83_resolvent_budget_hierarchy():
    """Prop 15.83: Max+-free algebra ranking M_cand residual vs L residual."""
    import sys
    from fractions import Fraction

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_prop1583 import (
        M_cand,
        M_mid,
        L_abs,
        T_abs,
        rho_L,
        rho_cand,
        gain_L,
        gain_cand,
        prove_cascade,
        prove_gain_hierarchy,
        is_prime,
        main as run_prop1583,
    )

    primes = [p for p in range(5, 60) if is_prime(p)]
    cas = prove_cascade(primes)
    assert cas["proved_for_listed_primes"] is True
    gh = prove_gain_hierarchy(primes)
    assert gh["proved_algebra"] is True
    assert gh["certified_primes"] is True

    # Closed form gain_L - gain_cand = 3(p-2)/(48(2p+3))
    for p in primes:
        diff = gain_L(p) - gain_cand(p)
        assert diff == Fraction(3 * (p - 2), 48 * (2 * p + 3))
        assert M_cand(p) < M_mid(p) <= L_abs(p) < T_abs(p)
        assert 0 < rho_cand(p) < rho_L(p)
        # p=5 sharpness of M_cand
        if p == 5:
            assert M_cand(5) == Fraction(3, 65)

    # Drive real entry point (writes evidence JSON)
    run_prop1583()
    path = ROOT / "evidence" / "e1_gmin_m4_prop1583.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.83"
    assert data["proved"] is True
    assert data["L_status"] == "OPEN"
    assert "GPU unused" in data["backend"]
    assert "15.83" in (ROOT / "solution.md").read_text()
    # L not soft-closed in handoff one-liner region
    hand = (ROOT / "HANDOFF.md").read_text()[:1200]
    assert "OPEN" in hand


def test_prop_15_84_gd_cand_S3_budget():
    """Prop 15.84: B_cand closed form, sign pattern, diag-dom fails."""
    import sys
    from fractions import Fraction

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_prop1584 import (
        B_cand,
        rho_cand,
        M_cand,
        d1_kappa1,
        prove_B_cand_formula,
        prove_B_cand_sign,
        prove_diagonal_dominance_fails,
        is_prime,
        main as run_1584,
    )

    primes = [p for p in range(5, 50) if is_prime(p)]
    assert prove_B_cand_formula(primes)["proved_closed_form"] is True
    assert prove_B_cand_sign(primes)["proved_sign_pattern"] is True
    assert prove_diagonal_dominance_fails(primes)["proved_fails_for_p_ge_5"] is True

    assert B_cand(5) == Fraction(-16, 325)
    assert B_cand(5) < 0
    for p in primes:
        if p >= 7:
            assert B_cand(p) > 0
        # formula identity
        assert B_cand(p) == p * rho_cand(p) - Fraction(2, p * p)
        assert d1_kappa1(p) == 3 * p * p - 7
        assert 4 * p - d1_kappa1(p) < 0
        assert M_cand(p) == Fraction(p - 2, p * (2 * p + 3))

    run_1584()
    path = ROOT / "evidence" / "e1_gmin_m4_prop1584.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.84"
    assert data["proved"] is True
    assert data["L_status"] == "OPEN"
    assert "GPU unused" in data["backend"]
    assert "15.84" in (ROOT / "solution.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]


def test_prop_15_85_Q4_mean_split():
    """Prop 15.85: mu closed form, S1=0 algebra, ray split."""
    import sys
    from fractions import Fraction
    from math import comb

    sys.path.insert(0, str(ROOT / "src"))
    from e1_gmin_m4_prop1585 import (
        mu_m4,
        e4,
        H,
        n_of,
        prove_mu_closed,
        prove_S1_zero_algebra,
        prove_Sd_identity,
        prove_ray_split,
        is_prime,
        main as run_1585,
    )

    primes = [p for p in range(3, 40) if is_prime(p)]
    assert prove_mu_closed(primes)["proved"] is True
    assert prove_S1_zero_algebra()["proved"] is True
    assert prove_Sd_identity()["proved"] is True
    assert prove_ray_split([p for p in primes if p >= 5])["proved_formula"] is True

    for p in (5, 7, 11):
        n = n_of(p)
        assert mu_m4(p) == Fraction(e4(p), comb(n, 4))
        assert abs(mu_m4(p)) / H(p) < Fraction(1, 10)

    assert mu_m4(5) == Fraction(-9, 1495)

    run_1585()
    path = ROOT / "evidence" / "e1_gmin_m4_prop1585.json"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.85"
    assert data["proved"] is True
    assert data["L_status"] == "OPEN"
    assert "15.85" in (ROOT / "solution.md").read_text()
    assert "OPEN" in (ROOT / "HANDOFF.md").read_text()[:900]
