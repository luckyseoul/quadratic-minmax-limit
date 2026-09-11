"""Paley two-block decomposition: entrywise identity and block-minimizer checks."""
import sys
from pathlib import Path

import numpy as np

TOOL = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(TOOL))

from paley_two_block_decomposition import char_split, check, phi  # noqa: E402


def test_entrywise_identity_small():
    for q in (5, 9, 13):
        C, T, Tc = check(q, verbose=False)
        A = C[np.ix_(T, T)]
        D = C[np.ix_(Tc, Tc)]
        assert (A == -D).all()


def test_decomposition_holds_for_prime_powers():
    for q in (17, 25, 29):
        check(q, verbose=False)


def test_block_phi_matches_recorded_minima():
    # recorded exact values: m_3=3 (q=5), m_7=9 (q=13), m_9=12 (q=17), m_15=27 (q=29)
    for q, m in ((5, 3), (13, 9), (17, 12)):
        C, T, Tc = check(q, verbose=False)
        assert phi(C[np.ix_(T, T)]) == m


def test_char_split_partitions_are_balanced():
    for q in (5, 9, 13, 17, 25):
        T, Tc, nu = char_split(q)
        n = q + 1
        assert len(T) == len(Tc) == n // 2
        assert sorted(T + Tc) == list(range(n))
