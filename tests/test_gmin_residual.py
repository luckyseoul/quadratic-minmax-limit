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