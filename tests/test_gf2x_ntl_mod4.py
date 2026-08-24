"""Regression checks for selected-line counts modulo four."""
from __future__ import annotations

import numpy as np
import pytest

import gf2x_ntl


@pytest.mark.skipif(not gf2x_ntl.available(), reason="NTL bridge unavailable")
def test_selected_line_mod4_counts_reduce_to_parity():
    # F_29[sigma]/(sigma^2-2), with the least primitive element returned by
    # the native search.  These are the twelve d=21 boundary levels.
    p = 29
    levels = [2, 3, 4, 11, 12, 13, 16, 17, 18, 25, 26, 27]
    parity, parity_offsets = gf2x_ntl.selected_line_bins(
        p,
        0,
        2,
        250,
        33,
        29,
        28,
        levels,
        [21],
        force_wide=True,
    )
    counts, count_offsets = gf2x_ntl.selected_line_counts_mod4(
        p, 0, 2, 250, 33, 29, 28, levels, [21]
    )

    assert parity_offsets == count_offsets == [0]
    assert counts.shape == parity.shape == (2, 12, 21)
    assert int(counts.max()) <= 3
    assert np.array_equal(counts & 1, parity)
