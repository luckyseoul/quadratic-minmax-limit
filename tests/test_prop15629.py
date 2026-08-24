from src.e1_gmin_m4_prop15629 import (
    PARI_CERT,
    all_circle_index,
    discriminant_group_invariants,
    dual_discriminant_theorem,
    glue_dimension,
    lattice_det,
    lattice_level,
    line_basis_index,
    m_of,
    pari_certificate,
    quotient_theorem,
)


def test_profile_glue_dimension_and_indices():
    for p in (3, 5, 7, 11, 13, 17):
        m = m_of(p)
        assert glue_dimension(p) == (m - 1) * (m - 2) // 2
        assert all_circle_index(p) == p ** ((m - 1) * (m - 2) // 2)
        assert line_basis_index(p) == p ** (m * (m - 1) // 2)
        assert lattice_det(p) == 2 * p ** (m * m)


def test_quotient_theorem_and_pari_certificates():
    assert quotient_theorem()["proved"]
    assert pari_certificate()["certified"]
    for p, row in PARI_CERT.items():
        assert row["det"] == lattice_det(p)
        assert row["all_circle_index"] == all_circle_index(p)
        assert row["line_index"] == line_basis_index(p)
        assert row["dual_den"] == 2 * p
        assert row["level"] == lattice_level(p)


def test_dual_discriminant_group_and_level():
    assert dual_discriminant_theorem()["proved"]
    for p in (3, 5, 7, 11, 13, 17, 19):
        m2 = m_of(p) ** 2
        invariants = discriminant_group_invariants(p)
        assert invariants == [p] * (m2 - 1) + [2 * p]
        assert len(invariants) == m2
        assert lattice_level(p) == 4 * p
