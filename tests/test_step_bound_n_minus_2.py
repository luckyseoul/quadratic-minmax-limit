"""The one-step bound m_{n+1} <= m_n + n - 2 is forced for n >= 5.

The Gaussian floor exceeds the flat second-moment value, so a minimizer
has a vector with |Q| <= m_n - 2. The certified table shows the saving of
2 is sharp.
"""

import math


CERTIFIED = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12, 10: 13}


def test_floor_exceeds_flat_second_moment_for_n_at_least_5():
    for n in range(5, 64):
        floor = n * math.sqrt(n - 1) / math.pi
        flat = math.sqrt(n * (n - 1) / 2)
        assert floor > flat


def test_certified_steps_save_at_least_two_and_two_are_sharp():
    for n in range(3, 10):
        assert CERTIFIED[n + 1] <= CERTIFIED[n] + n - 2
    assert CERTIFIED[4] == CERTIFIED[3] + 3 - 2
    assert CERTIFIED[7] == CERTIFIED[6] + 6 - 2
