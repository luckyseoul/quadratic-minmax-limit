"""Drive shipped moments.py; leftover flags stay False. k=6 table not required."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MOMENTS = ROOT / "scripts" / "maxplus_profile_enum" / "moments.py"
ANALYZE = ROOT / "scripts" / "maxplus_profile_enum" / "analyze_k6_p13.py"
P7 = Path("/mnt/storage/e1work/maxplus_p13/gt_k4_p7.npy")
LOG = Path("/mnt/storage/e1work/maxplus_p13/enum_p13_k6.log")


def test_analyze_refuses_until_k6_complete():
    r = subprocess.run(
        [sys.executable, "-u", str(ANALYZE)],
        capture_output=True,
        text=True,
    )
    text = LOG.read_text() if LOG.exists() else ""
    done = "[17805/17805]" in text and "k=6 gauged GPU TOTAL" in text
    if not done:
        assert r.returncode == 2
        assert "NOT FINISHED" in r.stdout
    else:
        assert r.returncode == 0


@pytest.mark.skipif(not P7.exists(), reason="p=7 k=4 ground-truth npy missing")
def test_moments_py_on_p7_k4_table():
    r = subprocess.run(
        [sys.executable, "-u", str(MOMENTS), "7", str(P7)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "max |mu| over |kappa|=1 four-sets:" in r.stdout
    assert "L=(p-2)/2p^2=" in r.stdout


def test_leftover_flags_untouched():
    from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    assert type_I_multilevel_bad_case_ND_closed() is False
    assert phi_F_ge_6_proved_general() is False
