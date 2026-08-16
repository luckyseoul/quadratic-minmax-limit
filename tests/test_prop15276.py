"""Prop 15.276: Lemma D existence (A3) and 2-plane amplitudes.

Fails if Cy=py is skipped (M3 count only), if A3_PROOF.md is missing,
or if amplitudes are hard-coded as 1+0j while the live DFT disagrees.
Drives shipped constructors, not a reimplementation of Paley C.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15276 import (
    A3_PROOF,
    M3_count,
    a3_proof_present,
    amplitude_vector_A5,
    certify_construction,
    construct_locked_triple,
    construct_sawtooth_z,
    lemma_D_2plane_amplitudes_proved,
    lemma_D_existence_written,
    live_check_prime,
    locked_phases,
    sawtooth_dft_G,
    theorem_A5_characters_distinct,
    theorem_amplitude_formula,
    theorem_amplitudes_not_parallel,
    theorem_fejer_never_zero,
)
from minmax_quadratic import paley_conference_prime_power


def test_a3_proof_missing_fails_existence():
    assert A3_PROOF.is_file()
    assert a3_proof_present() is True
    text = A3_PROOF.read_text()
    for needle in (
        "occupanc",
        "N(x)=1+",
        "s_0+s_1+s_2",
        "Cy=py",
        "majority",
        "three dual",
        "F_{\\lambda,s}",
    ):
        assert needle in text
    assert "prize" not in text.lower()
    src = inspect.getsource(lemma_D_existence_written)
    assert "a3_proof_present" in src
    assert "certify_construction" in src or "live_check" in src


def test_m3_count_is_not_existence():
    """M3 matching enum is a check. Existence requires Cy=py."""
    assert M3_count(5) == 100
    assert M3_count(7) == 1176
    assert M3_count(11) == 24200
    src = inspect.getsource(lemma_D_existence_written)
    assert "M3_count" not in src or "cy_eq_py" in src
    rec = live_check_prime(5)
    assert rec["M3"] == 100
    assert rec["cy_eq_py"] is True
    assert rec["zhat0_is_p"] is True
    assert rec["used_paley_conference_prime_power"] is True
    # existence source must demand the live Cy check, not the count
    assert "cy_eq_py" in src


def test_live_cy_py_uses_shipped_paley():
    src_mod = inspect.getsource(live_check_prime)
    assert "paley_conference_prime_power" in src_mod
    src_file = inspect.getsource(
        inspect.getmodule(live_check_prime)
    )
    assert "from minmax_quadratic import paley_conference_prime_power" in src_file
    # drive the shipped matrix, do not rebuild C
    rec = construct_locked_triple(5, lam=1, s0=0, s1=1)
    C = paley_conference_prime_power(5)
    err = float(np.max(np.abs(C @ rec["y"] - 5 * rec["y"])))
    assert err < 1e-8
    rec7 = live_check_prime(7)
    assert rec7["cy_eq_py"]
    assert rec7["support_in_omega"]
    assert rec7["support_on_three_lines"]


def test_amplitudes_not_hardcoded_one():
    live = sawtooth_dft_G(5, 1, 0, 1)
    form = theorem_amplitude_formula()
    assert form["proved"]
    assert form["not_hardcoded_one"]
    assert abs(live - (1 + 0j)) > 1.0
    cert = certify_construction()
    for p, rec in cert["by_p"].items():
        assert rec["amplitudes_not_one"] is True
        assert rec["fejer_matches_live"] is True
        assert rec["edge_formula_matches_live"] is True
        for a in rec["edge_amps_formula"]:
            # stored as strings of complex values; 1+0j is forbidden
            assert a not in ("(1+0j)", "1+0j", "1", "(1+0j)")
    src = inspect.getsource(sawtooth_dft_G)
    assert "omega" in src
    assert "2 * p" in src or "2*p" in src
    assert "1 + 0j" not in src or "hardcoded" in inspect.getsource(
        theorem_amplitude_formula
    )


def test_amplitude_formula_matches_live_and_spans_plane():
    assert theorem_fejer_never_zero()["proved"] is True
    assert theorem_A5_characters_distinct()["proved"] is True
    par = theorem_amplitudes_not_parallel()
    assert par["proved"] is True
    s = locked_phases(0, 0, 5)
    a = amplitude_vector_A5(5, 1, 2, 3, *s, 1)
    b = amplitude_vector_A5(5, 1, 2, 4, *s, 1)
    assert np.min(np.abs(a)) > 1e-8
    assert abs(a[0] * b[1] - a[1] * b[0]) > 1e-8
    rec = live_check_prime(5)
    assert rec["rank2_across_mu"] is True
    rec7 = live_check_prime(7)
    assert rec7["rank2_across_mu"] is True


def test_flags_require_writeup_live_and_formula():
    assert lemma_D_existence_written() is True
    assert lemma_D_2plane_amplitudes_proved() is True
    src_ex = inspect.getsource(lemma_D_existence_written)
    src_pl = inspect.getsource(lemma_D_2plane_amplitudes_proved)
    assert "a3_proof_present" in src_ex
    assert "theorem_amplitude_formula" in src_ex
    assert "theorem_amplitudes_not_parallel" in src_pl
    # 2-plane is Fejer + formula, not only certified rank at 5,7,11
    assert "theorem_fejer_never_zero" in inspect.getsource(
        theorem_amplitudes_not_parallel
    )


def test_construct_occupancy_lock():
    rec = construct_locked_triple(5, lam=2, s0=3, s1=4)
    assert rec["s"][0] + rec["s"][1] + rec["s"][2] == 3 + 4 + ((-2 - 3 - 4) % 5)
    assert sum(rec["s"]) % 5 == (-2) % 5
    z = rec["z"]
    assert z.shape == (25,)
    assert set(np.unique(z)).issubset({-1.0, 1.0})
    # occupancy constructor used
    src = inspect.getsource(construct_sawtooth_z)
    assert "occupancy" in src
    assert "2 * p + 1" in src or "2*p+1" in src
