"""Soft-soft stop flags, independent start, starve vs in-flight locks."""
from __future__ import annotations

import os
import sys
from pathlib import Path

TOOL = Path(__file__).resolve().parents[1] / "scripts" / "maxplus_profile_enum"
sys.path.insert(0, str(TOOL))

from k6_control import (  # noqa: E402
    STARVE_MARK,
    prepare_start,
    request_stop,
    should_stop,
    starve_unclaimed,
    drop_starve_locks,
    stop_flag_state,
    try_claim_orbit,
    clear_stop,
)


def test_should_stop_host_and_all(tmp_path):
    root = tmp_path / "stop"
    assert should_stop("v100", root) is False
    request_stop("v100", root)
    assert should_stop("v100", root) is True
    assert should_stop("nuka", root) is False
    request_stop("ALL", root)
    assert should_stop("nuka", root) is True
    assert should_stop("a380", root) is True
    st = stop_flag_state(root)
    assert st["ALL"] is True
    assert st["v100"] is True
    assert st["dash"] is False
    clear_stop("ALL", root)
    assert should_stop("nuka", root) is False
    assert should_stop("v100", root) is True


def test_prepare_start_does_not_unstop_others(tmp_path):
    root = tmp_path / "stop"
    request_stop("ALL", root)
    prepare_start("v100", root)
    st = stop_flag_state(root)
    assert st["ALL"] is False
    assert st["v100"] is False
    assert should_stop("v100", root) is False
    assert should_stop("nuka", root) is True
    assert should_stop("orin", root) is True
    assert should_stop("a380", root) is True
    assert should_stop("cpu44", root) is True
    prepare_start("nuka", root)
    assert should_stop("nuka", root) is False
    assert should_stop("v100", root) is False
    assert should_stop("orin", root) is True


def test_try_claim_orbit_soft_skip_is_minus_two(tmp_path):
    out = tmp_path / "gauge"
    stop = tmp_path / "stop"
    out.mkdir()
    os.environ["K6_HOST"] = "v100"
    os.environ["K6_STOP_DIR"] = str(stop)
    try:
        assert try_claim_orbit(out, 7, host="v100", root=stop) == 0
        assert (out / ".lock_7").is_dir()
        request_stop("v100", stop)
        assert try_claim_orbit(out, 8, host="v100", root=stop) == -2
        assert not (out / ".lock_8").exists()
        (out / "orb9.npy").write_bytes(b"x")
        assert try_claim_orbit(out, 9, host="nuka", root=stop) == -1
    finally:
        os.environ.pop("K6_HOST", None)
        os.environ.pop("K6_STOP_DIR", None)


def test_starve_skips_inflight_and_drop_only_dummy(tmp_path):
    out = tmp_path / "gauge"
    out.mkdir()
    (out / "orb0.npy").write_bytes(b"x")
    inflight = out / ".lock_1"
    inflight.mkdir()
    rec = starve_unclaimed(out, n_tasks=5)
    assert rec["npy"] == 1
    assert rec["inflight"] == [1]
    assert rec["dummy_locks"] == 3  # 2,3,4
    assert (out / ".lock_2" / STARVE_MARK).is_file()
    assert not (inflight / STARVE_MARK).exists()
    n = drop_starve_locks(out)
    assert n == 3
    assert inflight.is_dir()
    assert (out / "orb0.npy").is_file()
    assert not (out / ".lock_2").exists()


def test_worker_nc_minus_two_matches_claim(tmp_path):
    """run_kgauged.worker returns nc=-2 when try_claim_orbit would skip."""
    stop = tmp_path / "stop"
    out = tmp_path / "gauge"
    out.mkdir()
    request_stop("cpu44", stop)
    assert try_claim_orbit(out, 0, host="cpu44", root=stop) == -2
    src = (TOOL / "run_kgauged.py").read_text()
    assert "try_claim_orbit" in src
    assert "claimed < 0" in src


def test_dashboard_html_is_mobile():
    from k6_control import WORKERS
    from k6_dashboard import render

    nodes = {
        w: {
            "name": w,
            "horse": "hp",
            "running": False,
            "stopping": False,
            "pid": "",
            "log": "<script>x</script>",
            "flag": False,
        }
        for w in WORKERS
    }
    page, _, _ = render(
        dict(
            eta_h=1.0, n_done=1, n_tasks=17805, pct=0.01, n_left=1,
            locks=0, starve_locks=0, rate=1.0, gained=1, hours=1.0,
            updated="now", nodes=nodes,
        )
    )
    assert "width=device-width" in page
    assert "@media (max-width:600px)" in page
    assert "min-height:44px" in page
    assert "&lt;script&gt;" in page
    assert "<script>x</script>" not in page.split("<script>")[0]
    mesh = (TOOL / "k6_mesh.sh").read_text()
    cpu = mesh.split("start_cpu44()")[1].split("start_dash()")[0]
    assert "run_kgauged.py 6 44" in cpu
    assert "GPU_WORKERS=44" in cpu
    assert "OMP_NUM_THREADS=1" in cpu
    assert "NUMBA_NUM_THREADS=1" in cpu
    assert "NUMBA_NUM_THREADS=44" not in cpu
    # Start/stop POST must spawn the mesh script next to the dashboard,
    # never the main-tree path that no longer has k6_mesh.sh.
    assert "if (!busy) location.reload()" in page
    from k6_dashboard import HERE, MESH
    assert MESH == HERE / "k6_mesh.sh"
    assert MESH.is_file()
    assert 'CODE="${K6_CODE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"' in mesh
    assert 'CODE="${K6_CODE:-/home/nick/quadratic-minmax-limit/scripts/maxplus_profile_enum}"' not in mesh
