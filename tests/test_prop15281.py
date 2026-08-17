"""Tests for Prop 15.281 — Jacobi/Hoffman/Gauss Cross; floor still OPEN.

Fail-when-wrong is required.  phi_F_ge_6 / e1 / L / 15.280 stay put.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15280 import prove_open as prove_open_15280
from e1_gmin_m4_prop15281 import (
    main,
    prove_hoffman_delta,
    prove_jacobi_cross,
    prove_open,
    prove_triv_1d,
)


def test_jacobi_cross_and_bad_eta():
    a = prove_jacobi_cross([5, 7])
    assert a["proved"] is True
    assert a["bad_eta_is_identity"] is False
    assert a["n_checks"] > 0
    assert a["by_p"]["5"]["ok"] == a["by_p"]["5"]["n"]
    assert a["by_p"]["7"]["ok"] == a["by_p"]["7"]["n"]


def test_hoffman_delta_and_invert_cy():
    b = prove_hoffman_delta([5, 7])
    assert b["proved"] is True
    assert b["invert_Cy_is_identity"] is False
    assert b["nonhoffman_breaks_B_p5"] is True
    assert b["half_formula_always"] is False
    assert b["n_checks"] > 0


def test_triv_1d_and_yp_ync_fail():
    c = prove_triv_1d([5, 7])
    assert c["proved"] is True
    assert c["yp_ratio_is_identity"] is False
    assert c["ync_triv_formula_is_identity"] is False
    assert c["by_p"]["5"]["pred_E"] == 0.0
    assert c["by_p"]["7"]["pred_E"] == 784.0


def test_open_floor_flags_untouched():
    d = prove_open()
    assert d["floor_closed"] is False
    assert d["Nplus_named"] is False
    assert d["mixture_meets_floor"] is False
    assert d["one_d_plus_hs_meets_floor_p7k8"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    # 15.280 D stays: linear Cross still not a general-type identity
    e = prove_open_15280()
    assert e["cross_closed"] is False
    assert e["floor_closed"] is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["jacobi_cross"] is True
    assert out["proved"]["hoffman_delta"] is True
    assert out["proved"]["triv_1d_energy"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    assert out["proved"]["invert_Cy_is_identity"] is False
    assert out["proved"]["bad_eta_is_identity"] is False
    assert out["status"]["floor"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15281.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.281"
    assert data["status"]["floor"] is False
    assert data["status"]["identities_A_B_C"] is True
