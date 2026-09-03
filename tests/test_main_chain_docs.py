"""Drive e1_main_chain_status: L OPEN, docs_ok, no soft-close."""
from __future__ import annotations

from e1_main_chain_status import check_docs_L_status, run_main_chain


def test_main_chain_L_open_and_docs_ok():
    out = run_main_chain()
    assert out["L_status"] == "OPEN"
    assert out["L_closed"] is False
    assert out["gsum_disj_lb_proved_general"] is False
    assert out["writeup_L_closed"] is False
    docs = out["docs"]
    assert docs["soft_close_detected"] is False
    assert docs["overclaim_detected"] is False
    assert docs["HANDOFF_shows_L_OPEN"] is True
    assert docs["docs_ok"] is True
    assert docs["scanned_denseness_package_full"] is True
    assert docs["e1_closed_general"] is False
    units = docs["four_e1_units"]
    assert units["bitight_levels_2_3"] is True
    assert units["residual_ii_k_ge_4p"] is False
    assert units["type_I_multilevel"] is True
    assert units["closed"] is False


def test_solution_does_not_assert_limit_theorem():
    from pathlib import Path

    text = (Path(__file__).resolve().parents[1] / "solution.md").read_text(
        encoding="utf-8", errors="replace"
    )
    head = text[:5000]
    assert "**Main Theorem (limit).**" not in head
    assert "E(1) and" in head and "are not complete" in head
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text(
        encoding="utf-8", errors="replace"
    )
    assert "blocked**\nby residual **(i)** only" not in readme
    assert "L=\\lim_n\\alpha_n$ is **OPEN**" in readme or "is **OPEN**" in readme

    long_goal = (
        Path(__file__).resolve().parents[1] / "LONG_HORIZON_GOAL.md"
    ).read_text(encoding="utf-8", errors="replace")
    assert "Residual **(ii)** is **CLOSED**" not in long_goal
    assert "**Residual (ii), even \\(k\\ge4p\\):** OPEN" in long_goal
    assert "required_bitight_levels_empty_all_primes" in long_goal
    assert "no longer acceptance gates" in long_goal
    assert "Existence CLOSED; value unidentified" in long_goal
    assert "Proposition 6.3" in long_goal
    assert "Proposition 6.8" in long_goal
    assert "Proposition 6.9" in long_goal
    assert "(6.42)--(6.43)" in long_goal
    assert "multipliers 2 and 3" in long_goal
    assert "Dini-summable" in long_goal

    paley_goal = (
        Path(__file__).resolve().parents[1] / "GOAL.md"
    ).read_text(encoding="utf-8", errors="replace")
    assert "[x] Close Type I bad case when Max− is multi-level (15.750)" in paley_goal
    assert "attack only the two\nopen units" not in paley_goal

    assert "**Proposition 6.3 (two-ray convergence criterion).**" in text
    assert "H(2n)\\le2H(n)" in text
    assert "H(3n)\\le3H(n)" in text
    assert "2^a3^b" in text
    assert "E(N)=\\sum_{j\\ge0}" in text
    assert "**Proposition 6.4 (exact four-state form of Hadamard doubling).**" in text
    assert "**Proposition 6.5 (equal-endpoint skew reduction).**" in text
    assert "**Proposition 6.5a (exact directed-half-cut reformulation).**" in text
    assert "**Proposition 6.5b (sharpness of the directed-half-cut multiplier).**" in text
    assert "**Proposition 6.5c (opposite-diagonal diamond and hybrid-slice complexification).**" in text
    assert "**Proposition 6.5d (bivector cover gate and low-degree no-go).**" in text
    assert "**Proposition 6.5e (signed-regular arcsine rigidity" in text
    assert "**Proposition 6.5f (finite-anchor signature-cell shielding).**" in text
    assert "**Proposition 6.5g (random approximate mate and generic spectral-bridge" in text
    assert "**Proposition 6.5h (exact outgoing-half random criterion and first-moment" in text
    assert "**Proposition 6.5i (Gaussian saturation and central two-half saddle).**" in text
    assert "**Proposition 6.5j (degree-four squared-row preordering no-go).**" in text
    assert "**Proposition 6.5k (weighted multi-anchor integral rounding).**" in text
    assert "**Proposition 6.5l (growing-degree squared-row preordering no-go).**" in text
    assert "**Proposition 6.5m (conference commuting-mate parity obstruction).**" in text
    assert "**Proposition 6.5n (optimal-scale coherent clique-flip counterfamily).**" in text
    assert "w^2\\ge3m-2" in text
    assert "\\mathbb E_0[g_bg_c]" in text
    assert "50\\sum_{a=1}^k{\\rho_a^2\\over B_a^2}\\le1" in text
    assert "c_*={3-2\\sqrt2\\over25\\pi}" in text
    assert "same resulting \\(R\\)" in text
    assert "raw polynomial degree" in text
    assert "(2D+1)^D" in text
    assert "(AR-RA)_{ii}" in text
    assert "\\|AR-RA\\|_F^2\\ge4n" in text
    assert "L_{\\rm cl}(A)-\\sqrt2\\Phi(A)" in text
    assert "This does not\ndisprove the still-open implication for global minimizers" in text
    assert "\\Phi({\\cal K}(A,C))" in text
    assert "D_{\\to}(A,S)" in text
    assert "b_{xy}={x\\wedge y\\over2}" in text
    assert "L_{\\rm cl}(A)=6>4\\sqrt2" in text
    assert "|O-X^2/(4J)|" in text
    assert "|O+J|+\\sqrt{(O-J)^2+X^2}\\le2M" in text
    assert "noncoherent" in text
    assert "NOTE_2026-09-02_COMPLEXIFICATION_OPPOSITE_DIAGONAL_AUDIT.md" in readme
    assert "NOTE_2026-09-02_BIVECTOR_ENERGY_LAYER_MINIMAX.md" in readme
    assert "NOTE_2026-09-02_SIGNED_REGULAR_ARCSINE_RIGIDITY.md" in readme
    assert "NOTE_2026-09-02_FINITE_ANCHOR_SIGNATURE_TOURNAMENT.md" in readme
    assert "NOTE_2026-09-02_RANDOM_SKEW_MATE_SECOND_MOMENT.md" in readme
    assert "NOTE_2026-09-02_DIRECTED_HALFCUT_RANDOM_ORIENTATION.md" in readme
    assert "NOTE_2026-09-02_GAUSSIAN_SATURATION_CENTRAL_SADDLE.md" in readme
    assert "NOTE_2026-09-02_BIVECTOR_DEGREE4_PREORDERING_NO_GO.md" in readme
    assert "NOTE_2026-09-02_BANASZCZYK_WEIGHTED_ANCHOR_ROUNDING.md" in readme
    assert "NOTE_2026-09-02_BIVECTOR_GROWING_DEGREE_PREORDERING_NO_GO.md" in readme
    assert "NOTE_2026-09-02_CONFERENCE_COMMUTING_MATE_NO_GO.md" in readme
    assert "NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md" in readme
    assert "**Proposition 6.6 (balanced Paley-skew shielding).**" in text
    assert "**Proposition 6.7 (tetrahedral tripling frame and exact diamond).**" in text
    assert "**Proposition 6.8 (bi-balanced Hadamard shield" in text
    assert "**Proposition 6.9 (conference obstruction to every fixed-temperature" in text
    assert "**Proposition 6.10 (critical-pressure gate and graphon no-go).**" in text
    assert "**Proposition 6.10a (product-pressure lower-curve no-go and entropy\nfallback).**" in text
    assert "{s_n(c)\\over c}\\le\\alpha_n" in text
    assert "{n\\over2}\\log\\cosh" in text
    assert "{c^2\\over4}" in text
    assert "NOTE_2026-09-02_THERMODYNAMIC_INTERPOLATION_GATE.md" in readme
    assert "NOTE_2026-09-02_PRESSURE_LOWER_CURVE_NO_GO.md" in readme
    assert "{5u^2-1\\over4}" in text
    assert "{c\\sqrt{1-1/n}\\over\\pi}r^2-I(r)" in text
    assert "{\\cal J}=\\begin{pmatrix}" in text
    assert "\\sqrt{q\\,k_A(x)k_B(y)}" in text
    assert "\\sqrt{q,k_A(x)k_B(y)}" not in text
    assert "three diagonal\nmatching terms" in text
    assert "k_A(x)k_B(y)>{n^2\\over100}" in text
    assert "false for every fixed" in text
    binding_docs = text + readme + long_goal
    assert "c=3` remains viable" not in binding_docs
    assert "c=3 remains viable" not in binding_docs
    assert "the **tetrahedral diamond**---is open" in text
    assert "does not prove the multiplier-three ray" in text
    assert "Every hereditary endpoint inequality (6.10) is automatic" in text
    assert r"{\cal F}(A)\ge" in text
    assert "This does **not**\nprove that random orientations fail" in text
    assert "finding a frame that passes all cuts remains open" not in text
    assert "constant \\(\\ge0.28\\) by bilinear lower bound" not in text
    assert "cannot beat the cross-term barrier" not in text
    assert "superseded by Proposition 15.750" in text
    assert "exact identity" in text
    assert "incorrectly combined" in text
    assert "C\\in\\Bigl[" not in text


def test_prop15757_15761_edge_radon_docs_preserve_open_box_gate():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    solution = (root / "solution.md").read_text(encoding="utf-8", errors="replace")
    for prop in range(15757, 15762):
        label = f"15.{prop - 15000}"
        assert f"## Proposition {label}" in solution

    assert "\\operatorname {rank}_{\\mathbf F_2}R_2" in solution
    assert "{(m-1)(4m^2+7m+6)\\over6}" in solution
    assert "{\\cal A}/R\\mathbf Z^{\\binom V2}" in solution
    assert "(z_0+\\ker_{\\mathbf Z}R)\\cap" in solution
    assert "{\\cal Q}(W,P,T)\\le|H|" in solution
    assert "Audited \\(p=31,t=69\\) arbitrary-compact local gate" in solution
    assert "449 free noncentered orbits" in solution
    assert "no `SAT` orbit exists" in solution
    assert "All-prime branch-C odd--Radon centrality" in solution
    assert "\\(p=4r+3\\)" in solution
    assert "\\(r\\ge7\\)" in solution
    assert "N=3(r+b-1)\\le4r-1" in solution
    assert "3b\\le r+2" in solution
    assert "0\\le\\delta\\le(2r+2)" in solution
    assert "69\\le t\\le99" in solution
    assert "unbalanced allocations" in solution
    assert "Seven-channel algebraic-dominance barrier" in solution
    assert "226534996574208000" in solution
    assert "220242357780480000" in solution
    assert "does **not** provide labels or form coefficients in \\(\\mathbf F_p\\)" in solution
    assert "No common simple graph" in solution or "not a common graph construction" in solution
    assert "Residual (ii), E(1), and the limit remain\nopen" in solution

    binding_names = (
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    binding_parts = [
        (root / name).read_text(encoding="utf-8", errors="replace")
        for name in binding_names
    ]
    binding = "\n".join(binding_parts)
    assert "p31 arbitrary-compact" in binding
    assert "nonzero even global forms" in binding
    assert "product_e {0,tau_e}" in binding
    assert "A/R(E)=(Z/pZ)^S(p)" in binding
    assert "no common simple graph is constructed" in binding.lower()
    assert "original MO limit remain **OPEN**" in binding

    for name, text in zip(binding_names, binding_parts):
        flattened = " ".join(text.split())
        assert "p=4r+3" in text, name
        assert "r>=7" in text, name
        assert (
            "zero odd global forms" in flattened
            or "all odd global forms vanish" in flattened
        ), name
        assert "3b<=r+2" in text, name
        assert "68<=t<=116" in text, name
        assert "69<=t<=99" in text, name
        assert "unbalanced" in text, name
        assert "dominant over the algebraic closure" in flattened, name
        assert "F_p" in text, name
        assert "NOTE_2026-09-03_CROSS_RECTANGLE_FOURIER_STABILITY.md" in text, name
        assert "Gram-only" in text, name
        assert "statewise diagonal-payment inequality" in flattened, name

    compact_note = (
        root / "evidence/NOTE_2026-09-02_COMPACT_RAY_HIGHER_MOMENT_GATE.md"
    ).read_text(encoding="utf-8", errors="replace")
    compact_flat = " ".join(compact_note.split())
    assert "no one $\\mathbf F_p$ labelling is yet known" in compact_flat
    assert "3b\\le r+2" in compact_flat
    assert "69\\le t\\le99" in compact_flat
    assert "does not supply admissible $\\mathbf F_p$ labels or form coefficients" in compact_flat
    assert "Residual (ii), E1" in binding


def test_post15761_geometry_conic_and_boolean_reductions_are_scoped():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    names = (
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "solution.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    texts = {
        name: (root / name).read_text(encoding="utf-8", errors="replace")
        for name in names
    }

    for name, text in texts.items():
        flat = " ".join(text.split())
        normalized = text.replace("\\beta", "beta").replace("{,}", ",")
        assert "Couvreur" in text, name
        assert "boundary cubic" in text.lower(), name
        assert "q^3=1" in text, name
        assert "k^2=-3" in text, name
        assert "p=31,b=7,k=11" in text, name
        assert "230,314,710" in normalized, name
        assert "17,076" in normalized, name
        assert "beta_R(y)=0" in normalized, name
        assert "complete Graver" in text, name
        assert "tests/test_p31_equi_zero68_mitm.py" in text, name
        assert (
            "Residual (ii), E1, `L=1/2`, and the original MO limit remain OPEN."
            in flat
        ), name
        assert "compact-threshold conjecture" not in text.lower(), name

    solution = texts["solution.md"]
    assert "{h(h+1)\\over2}>3h-6" in solution
    assert "3p-15>p+1+8\\sqrt p" in solution
    assert "p\\equiv7\\pmod {12}" in solution
    assert "F_6=(11,19,10)" in solution
    assert "F_8=(12,11,23,6)" in solution
    assert "\\ker_{\\mathbf Z}R/{\\cal K}_{\\rm ridge}" in solution
    assert "\\nu_p=dpm^2+{m(m-1)(4m+1)\\over6}" in solution
    assert (
        "not the complete Graver system and not a signed Boolean lift"
        in " ".join(solution.split())
    )

    for name in names:
        if name == "solution.md":
            continue
        text = texts[name]
        assert "p=7 mod 12" in text, name
        assert "K_ridge" in text, name
        assert "nu_p=d p m^2+m(m-1)(4m+1)/6" in text, name

    binding = "\n".join(texts.values())
    assert "finite-fiber theorem" in binding
    assert "not a residual close" in binding
    assert "no complete graver basis" in binding.lower()
    assert "do not establish a Boolean lift or closure" in binding


def test_frozen_component_and_inversion_upgrade_is_canonical_and_scoped():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    binding_names = (
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    binding = {
        name: (root / name).read_text(encoding="utf-8", errors="replace")
        for name in binding_names
    }
    for name, text in binding.items():
        flat = " ".join(text.split())
        assert "b>=(2r+7)/3" in text, name
        assert "3b<=2r+4" in text, name
        assert "68<=t<=164" in text, name
        assert "p=43,b=9" in text, name
        assert "antisymmetric Boolean half" in flat, name
        assert "coupled symmetric" in flat, name
        assert "2^32*3^26*5^2*7*2161" in text, name
        assert "NOTE_2026-09-03_EQUIANHARMONIC_COMPONENT_PACKING.md" in text, name
        assert "NOTE_2026-09-03_HARD_ROW_COMPACT_ODD_RADON_CENTRALITY.md" in text, name
        assert "NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md" in text, name
        assert "NOTE_2026-09-03_HARD_STAR_ANTISYMMETRIC_SUPPORT.md" in text, name
        assert "NOTE_2026-09-03_EQUIANHARMONIC_THRESHOLD_EVEN_BARRIER.md" in text, name
        assert (
            "Residual (ii), E1, `L=1/2`, and the original MO limit remain OPEN."
            in flat
        ), name

    solution = (root / "solution.md").read_text(
        encoding="utf-8", errors="replace"
    )
    solution_flat = " ".join(solution.split())
    assert "Equianharmonic component-packing upgrade" in solution
    assert "{2r+7\\over3}" in solution
    assert "3b\\le2r+4" in solution
    assert "68\\le t\\le164" in solution
    assert "Hard compact residual and the closed antisymmetric Boolean half" in solution
    assert "{\\cal A}^-/R(E^-)" in solution
    assert "h(h-1)(h+1)/3" in solution
    assert "R_Lz=S_{-j}-S_j" in solution
    assert "s_e\\in\\{0,2\\}" in solution
    assert "4128623683475967290061619200" in solution
    assert "2^{32}3^{26}5^2\\,7\\,2161" in solution
    assert "does not construct a rational finite-field zero" in solution_flat
    assert "Residual (ii), E1, `L=1/2`, and the original MO limit remain" in solution
    for note in (
        "NOTE_2026-09-03_EQUIANHARMONIC_COMPONENT_PACKING.md",
        "NOTE_2026-09-03_HARD_ROW_COMPACT_ODD_RADON_CENTRALITY.md",
        "NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md",
        "NOTE_2026-09-03_HARD_STAR_ANTISYMMETRIC_SUPPORT.md",
        "NOTE_2026-09-03_EQUIANHARMONIC_THRESHOLD_EVEN_BARRIER.md",
    ):
        assert note in solution


def test_frozen_symmetric_lattice_mobius_overlap_and_all_active_bound():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    binding_names = (
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    for name in binding_names:
        text = (root / name).read_text(encoding="utf-8", errors="replace")
        flat = " ".join(text.split())
        assert "(h-1)(2h^2+5h+6)/6" in text, name
        assert "mod-two" in text, name
        assert "tau_t=eta(Q(e1-t^2(e1+e2)))" in text, name
        assert "2(p-1)-2" in text, name
        assert "c>=p-1" in text, name
        assert "central Boolean" in flat, name
        assert "NOTE_2026-09-03_INVERSION_SYMMETRIC_LATTICE.md" in text, name
        assert "NOTE_2026-09-03_MOBIUS_HALF_SYMMETRIC.md" in text, name
        assert "NOTE_2026-09-03_ALL_ACTIVE_PENCIL_SUPPORT.md" in text, name
        assert "remain OPEN" in flat or "remain **OPEN**" in flat, name

    solution = (root / "solution.md").read_text(
        encoding="utf-8", errors="replace"
    )
    solution_flat = " ".join(solution.split())
    assert "The inversion-symmetric lattice and corrected Mobius capacity" in solution
    assert "{dh(p^2+1)\\over2}" in solution
    assert "S_+(p)={(h-1)(2h^2+5h+6)\\over6}" in solution
    assert "P_N(E)=" in solution
    assert "Q(e_1-t^2(e_1+e_2))" in solution
    assert "T_U=Y-Rq_U={Y+IY-RC_U\\over2}" in solution
    assert "2(p-1)-2" in solution
    assert "All-active support sharpening" in solution
    assert "c\\ge p-1" in solution
    assert "all-active hypothesis is indispensable" in solution_flat.lower()
    assert "not a restricted central Boolean lift" in solution_flat
    for note in (
        "NOTE_2026-09-03_INVERSION_SYMMETRIC_LATTICE.md",
        "NOTE_2026-09-03_MOBIUS_HALF_SYMMETRIC.md",
        "NOTE_2026-09-03_ALL_ACTIVE_PENCIL_SUPPORT.md",
    ):
        assert note in solution


def test_fixed_elimination_halved_code_and_rigid_overlap_are_canonical():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    names = (
        "solution.md",
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    texts = {
        name: (root / name).read_text(encoding="utf-8", errors="replace")
        for name in names
    }
    notes = (
        "NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md",
        "NOTE_2026-09-03_SYMMETRIC_HALVED_MOD2.md",
        "NOTE_2026-09-03_SYMMETRIC_HALVED_MOBIUS_COVER.md",
        "NOTE_2026-09-03_SYMMETRIC_UNUSED_SLICE_EXCHANGE.md",
        "NOTE_2026-09-03_MOBIUS_HALF_INTERSECTIONS.md",
    )
    for name, text in texts.items():
        flat = " ".join(text.split())
        assert "R+=[[A,2B],[0,C]]" in text or "R^+=" in text, name
        assert "a_[v]" in text or "a_{[v]}" in text, name
        assert "D=(C,Phi)" in text or "D=(C,\\Phi)" in text, name
        assert "d h(h+1)" in text or "dh(h+1)" in text, name
        assert "p h=|Delta|-h" in text or "ph=|\\Delta|-h" in text, name
        assert "A_(h-1)" in text or "A_{h-1}" in text, name
        assert "q=r=1/2" in text, name
        assert "t_max-t+1" in text or "t_{\\max}-t+1" in text, name
        assert "punctur" in text.lower(), name
        assert "row-code gap" in text or "minimum words" in text, name
        assert "remain OPEN" in flat, name
        for note in notes:
            assert note in text, (name, note)

    solution = texts["solution.md"]
    solution_flat = " ".join(solution.split())
    assert "Fixed-edge elimination and the halved symmetric code" in solution
    assert "MM^{\\mathsf T}=M^{\\mathsf T}M=I" in solution
    assert "\\operatorname {rank}D=dh(h+1)" in solution
    assert "X_{L,\\beta}" in solution
    assert "\\boxed{d_{\\rm row}(D)=ph}" in solution
    assert "there is no weight strictly between" in solution_flat
    assert "structured \\(D_U\\) is onto" in solution
    assert "no free parameter for a greedy multi-pair construction" in solution_flat
    assert "Pairwise overlap counts alone control neither test" in solution_flat

    integrated = (
        "solution.md",
        "STATUS.md",
        "HANDOFF.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    new_notes = (
        "NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE.md",
        "NOTE_2026-09-03_PRESCRIBED_CENTER_COMMON_BLOCK.md",
        "NOTE_2026-09-03_RIGID_PAIR_FIXED_WORD.md",
        "NOTE_2026-09-03_GROUPED_UNCERTAINTY_SQUARE.md",
        "NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE_GAP.md",
        "NOTE_2026-09-03_SYMMETRIC_QUOTA_CARDINALITY_BARRIER.md",
    )
    for name in integrated:
        text = texts[name]
        flat = " ".join(text.split())
        assert "group-support" in text, name
        assert "pseudoforest" in text, name
        assert "178-t" in text or "178-t" in flat, name
        assert "a_compact" in text or "a_{\\rm compact}" in text, name
        assert "residual (ii)" in text and "OPEN" in text, name
        for note in new_notes:
            assert note in text, (name, note)

    assert "canonical remainder vanishes to order at least two" in solution_flat
    assert "actual transverse target" in solution_flat
    assert "conditional rigidity, not existence" in solution_flat
    assert "a saturated equal-square common-block incidence cover" in solution_flat


def test_new_direct_no_gos_preserve_the_exact_open_scopes():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    names = (
        "solution.md",
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    )
    texts = {
        name: (root / name).read_text(encoding="utf-8", errors="replace")
        for name in names
    }
    for name, text in texts.items():
        assert "Proposition 6.5n" in text, name
        assert "Proposition 6.10a" in text, name

    binding = "\n".join(texts.values())
    assert "NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md" in binding
    assert "NOTE_2026-09-02_PRESSURE_LOWER_CURVE_NO_GO.md" in binding
    assert "NOTE_2026-09-03_TWO_HALF_SELF_GLUING_OBSTRUCTION.md" in binding
    assert "Phi(A)=m_n" in binding
    assert "large-temperature slope only `1/pi`" in binding
    assert "The quadratic-minmax limit is still OPEN" in texts["HANDOFF.md"]
    assert "original MO limit remain **OPEN**" in texts["STATUS.md"]


def test_soft_close_detector_flags_bare_L_CLOSED():
    """Unit-level: detector pattern must catch bare '**L CLOSED.**' style."""
    import re

    pat = re.compile(r"\*\*L CLOSED\.\*\*|\*\*L CLOSED\*\*|E\(1\)\s*\(CLOSED|E1_closed_general=true\. L CLOSED", re.I)
    assert pat.search("**E(1) (CLOSED, 15.168–171):** ... **L CLOSED.**")
    assert pat.search("E1_closed_general=true. L CLOSED (via E1∧bi-tight).")
    assert not pat.search("**L OPEN.** Denseness path blocked")


def test_props_15170_171_body_not_soft_closed():
    """solution.md Props 15.170–171 must not claim residual/E1/L CLOSED while hinge open."""
    from pathlib import Path
    import re

    sol = Path(__file__).resolve().parents[1] / "solution.md"
    text = sol.read_text(encoding="utf-8", errors="replace")
    # extract from Prop 15.170 onward
    i = text.find("## Prop 15.170")
    assert i >= 0
    tail = text[i:]
    # forbidden soft-close phrases (post-retraction)
    bad = [
        r"E1_closed_general=true\. L CLOSED",
        r"\*\*E\(1\) closed\*\*",
        r"Closes residual \(i\) of E\(1\) for all primes",
        r"Closes residual \(ii\) of E\(1\) for all primes",
        r"association-scheme min \$-12",
        r"disj: association-scheme min",
    ]
    for pat in bad:
        assert not re.search(pat, tail), f"soft-close residue matched: {pat}"
    assert "L OPEN" in tail or "OPEN for general" in tail or "not complete" in text[:4000]
    assert "gsum_disj_lb_proved_general" in tail or "NOT proved for general" in tail


def test_denseness_package_does_not_soft_close_the_theorem():
    """The stand-alone package must agree with its caveats throughout."""
    from pathlib import Path

    package = (
        Path(__file__).resolve().parents[1]
        / "evidence"
        / "share"
        / "denseness_path_package.md"
    ).read_text(encoding="utf-8", errors="replace")
    assert "does **not** prove E(1)" in package
    assert "## Lemma K (required bi-tight levels)." in package
    assert "| Residual (i) Type I, multi-level Max− | **Proved** (15.750) |" in package
    assert "type_I_multilevel_bad_case_ND_closed` stays False" not in package
    assert "Hence E(1) on the whole Paley family" not in package
    assert "Historical remarks “\\(L\\) OPEN”" not in package


def test_docs_track_all_prime_group_support_without_closing_residual():
    """The all-prime theorem must supersede the finite ladder but not residual."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    status = (root / "STATUS.md").read_text(encoding="utf-8")
    handoff = (root / "HANDOFF.md").read_text(encoding="utf-8")
    solution = (root / "solution.md").read_text(encoding="utf-8")
    dedup = (root / "evidence" / "PROPOSITION_DEDUP_AUDIT_2026-08-30.md").read_text(
        encoding="utf-8"
    )
    for text in (status, handoff, dedup):
        flat = " ".join(text.split())
        assert "group-support" in text
        assert "d_row(D)=p h" in text
        assert "minimum words" in text
        assert "superseded" in text
        assert "prescribed" in flat and "Boolean" in text
    assert "canonical remainder vanishes to order at least two" in " ".join(
        solution.split()
    )
    assert "\\boxed{d_{\\rm row}(D)=ph}" in solution
    assert "SYMMETRIC_HALVED_ROW_CODE_GAP" in solution
    assert "residual (ii)" in solution
