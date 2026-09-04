from e1_gmin_m4_antisymmetric_selector_kernel import (
    antisymmetric_selector_countercircuit_certificate,
    antisymmetric_transverse_countercircuit,
    boundary_selector_signature,
    clique_star_kernel_half,
    full_chain_boundary_invariance_theorem,
    mod_two_boundary,
    unsigned_selector_code_theorem,
    unsigned_unit_selector_kernel_half,
)
from e1_gmin_m4_inversion_antisymmetric_radon import (
    _negative_edge,
    edge_radon_image,
)


def test_full_chain_boundary_is_target_determined_symbolically():
    for p in (3, 5, 7, 11, 31):
        row = full_chain_boundary_invariance_theorem(p)
        assert row["full_chain_boundary_target_determined"]
        assert row["full_chain_selector_signature_target_determined"]
        assert not row["physical_half_signature_target_determined"]
        assert row["proved"]


def test_transverse_circuit_is_exactly_antisymmetric_and_Radon_zero():
    for p in (3, 5, 7, 11, 31):
        source, positive = antisymmetric_transverse_countercircuit(p)
        assert len(positive) == 2 * p
        assert len(source) == 4 * p
        assert edge_radon_image(p, source) == {}
        assert all(
            source[_negative_edge(p, edge)] == -coefficient
            for edge, coefficient in source.items()
        )
        assert mod_two_boundary(source) == frozenset()


def test_positive_half_changes_the_full_signature_and_aggregate():
    for p in (3, 5, 7, 11, 31):
        source, positive = antisymmetric_transverse_countercircuit(p)
        boundary = mod_two_boundary(positive)
        assert boundary == frozenset(
            {(2 % p, y) for y in range(p)}
            | {(-2 % p, y) for y in range(p)}
        )
        assert boundary_selector_signature(p, boundary) == (0,) + (1,) * p
        row = antisymmetric_selector_countercircuit_certificate(p)
        assert row["selector_signature_weight"] == p
        assert row["selector_signature_aggregate_mod_2"] == 1
        assert not row["physical_half_signature_target_determined"]
        assert not row["aggregate_selector_parity_target_determined"]
        assert not row["produces_residual_witness"]
        assert row["proved"]


def test_p31_countercircuit_survives_normalized_Paley_signs():
    row = antisymmetric_selector_countercircuit_certificate(31)
    assert row["positive_half_edge_count"] == 62
    assert row["anti_chain_support"] == 124
    assert row["positive_half_boundary_weight"] == 62
    assert row["full_chain_boundary_weight"] == 0
    assert row["p31_paley_edge_sign_set"] == [1]
    assert row["p31_normalized_residual_countercircuit_unchanged"]


def test_clique_star_is_an_exact_unsigned_anti_kernel_half():
    for p in (3, 7, 11):
        positive = clique_star_kernel_half(p, 0, p)
        source = {
            edge: coefficient
            for edge in positive
            for coefficient in (1,)
        }
        anti_source = dict(source)
        anti_source.update(
            {_negative_edge(p, edge): -1 for edge in positive}
        )
        assert len(positive) == p * (p - 1) // 2
        assert edge_radon_image(p, anti_source) == {}
        signature = boundary_selector_signature(p, mod_two_boundary(positive))
        assert signature[0] == signature[p] == 0
        assert sum(signature) == p - 1


def test_every_unsigned_projective_unit_has_a_ternary_kernel_half():
    for p in (3, 7):
        for wanted in range(p + 1):
            positive = unsigned_unit_selector_kernel_half(p, wanted)
            assert len(positive) == p * (p + 3) // 2
            assert boundary_selector_signature(
                p, mod_two_boundary(positive)
            ) == tuple(int(index == wanted) for index in range(p + 1))
    for p in (11, 31):
        row = unsigned_selector_code_theorem(p)
        assert row["minimum_nonzero_codeword_weight"] == 1
        assert row["induced_linear_code_dimension"] == p + 1
        assert row["every_projective_unit_realized"]
        assert not row["normalized_Paley_code_fullness_proved"]
        assert row["normalized_Paley_minimum_weight"] == "OPEN"
        assert row["proved"]
