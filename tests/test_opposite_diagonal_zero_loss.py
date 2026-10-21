"""A global minimizer can have a coherent clique flip above sqrt(2) Phi."""
from __future__ import annotations

import math

from minmax_quadratic import (
    clique_flip_norm,
    doubling_image,
    exact_m,
    order_four_optimal_signing,
    phi,
)


def test_order_four_minimizer_breaks_the_zero_loss_lemma():
    signing = order_four_optimal_signing()
    height = phi(signing)
    assert height == exact_m(signing.shape[0])
    lifted = clique_flip_norm(signing)
    assert lifted > math.sqrt(2.0) * height
    order = signing.shape[0]
    image = doubling_image(height)
    assert 2.0 * lifted - order < image < 2.0 * lifted + order
    excess = (lifted - math.sqrt(2.0) * height) / (order ** 1.5)
    assert 0.0 < excess < 0.5 - 1.0 / math.pi
