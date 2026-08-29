import math

from e1_gmin_m4_prop15691 import (
    adversarial_upper_rate,
    fractional_moment_barrier,
    sufficient_lower_rate,
    theorem_record,
)


def test_c2_target_is_false_by_linear_margin():
    assert adversarial_upper_rate(2) < sufficient_lower_rate(2) == 0
    assert math.isclose(
        -adversarial_upper_rate(2),
        (1 - math.sqrt(math.log(2))) ** 2,
    )


def test_c3_survives_this_no_go_and_barrier_is_correct():
    assert 2.08 < fractional_moment_barrier() < 2.09
    assert adversarial_upper_rate(3) > sufficient_lower_rate(3) == -0.75
    assert theorem_record()["c3"]["excluded_by_fractional_moment"] is False
