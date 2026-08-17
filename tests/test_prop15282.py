"""Tests for Prop 15.282 — Gauss-fibre / circ / n_T / generic Cy; floor OPEN.

Fail-when-wrong is required.  phi_F_ge_6 / e1 / L / 15.280–281 stay put.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15280 import prove_open as prove_open_15280
from e1_gmin_m4_prop15281 import prove_open as prove_open_15281
from e1_gmin_m4_prop15282 import (
    main,
    one_d_hs_yp_mix_p7_k8,
    prove_circ,
    prove_gauss_fibre,
    prove_generic,
    prove_nT,
    prove_open,
)


def test_gauss_fibre_all_slots_and_invert_cy():
    a = prove_gauss_fibre([5, 7])
    assert a["proved"] is True
    assert a["invert_Cy_is_identity"] is False
    assert a["hoffman_ok"] is True
    assert a["n_checks"] > 0
    for sl in ("triv", "circ", "gen"):
        hit, n = a["by_slot"][sl].split("/")
        assert int(hit) == int(n) > 0
    assert a["by_p"]["5"]["ok"] == a["by_p"]["5"]["n"]
    assert a["by_p"]["7"]["ok"] == a["by_p"]["7"]["n"]


def test_circ_and_triv_on_circ_not_identity_at_p7():
    b = prove_circ([5, 7])
    assert b["proved"] is True
    assert b["eta_fp_is_identity"] is False
    assert b["triv_on_circ_p7_is_identity"] is False
    assert b["Lambda_naive_no_chi_is_identity"] is False
    hits, n = (int(x) for x in b["triv_on_circ"]["7"].split("/"))
    assert n > 0
    assert hits < n
    assert b["by_p"]["5"]["Cy"].split("/")[0] == b["by_p"]["5"]["Cy"].split("/")[1]
    assert b["by_p"]["7"]["Cy"].split("/")[0] == b["by_p"]["7"]["Cy"].split("/")[1]


def test_nT_independent_of_t0():
    c = prove_nT([5, 7])
    assert c["proved"] is True
    assert c["t0_wrong_is_identity"] is False
    assert c["n_checks"] > 0
    for p in ("5", "7"):
        gok, gn = (int(x) for x in c["by_p"][p]["geom"].split("/"))
        assert gok == gn > 0
        iok, inn = (int(x) for x in c["by_p"][p]["t0_indep"].split("/"))
        assert iok == inn > 0


def test_generic_and_hs_cy_on_ync_fails():
    d = prove_generic([5, 7])
    assert d["proved"] is True
    assert d["ync_hs_formula_is_identity"] is False
    assert d["n_checks"] > 0
    assert d["by_p"]["5"]["ok"] == d["by_p"]["5"]["n"]
    assert d["by_p"]["7"]["ok"] == d["by_p"]["7"]["n"]
    hits, n = (int(x) for x in d["ync_hs"].split("/"))
    assert n > 0
    assert hits < n


def test_one_d_hs_yp_mix_p7_k8_below_floor():
    mix = one_d_hs_yp_mix_p7_k8()
    assert mix["floor"] == 441.0
    assert mix["mix"] < 441.0
    assert mix["meets_floor"] is False
    assert mix["den_census"] == 5726
    # live energies: Paley vanishes, AP / hs / yP match 15.140 table
    assert abs(mix["energies"]["Paley"]) < 1e-6
    assert abs(mix["energies"]["AP"] - 4116.0) < 1e-3
    assert abs(mix["energies"]["hs⊙nc"] - 1596.0) < 1e-3
    assert abs(mix["energies"]["yP⊙nc"] - 84.0) < 1e-3


def test_open_floor_flags_untouched():
    e = prove_open()
    assert e["floor_closed"] is False
    assert e["Nplus_named"] is False
    assert e["mixture_meets_floor"] is False
    assert e["one_d_plus_hs_yp_meets_floor_p7k8"] is False
    assert e["one_d_plus_hs_yp_mix_p7k8"] < 441.0
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    f = prove_open_15281()
    assert f["floor_closed"] is False
    assert f["mixture_meets_floor"] is False
    g = prove_open_15280()
    assert g["cross_closed"] is False
    assert g["floor_closed"] is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["gauss_fibre"] is True
    assert out["proved"]["circ_cross"] is True
    assert out["proved"]["nT"] is True
    assert out["proved"]["generic_cross"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    assert out["proved"]["invert_Cy_is_identity"] is False
    assert out["proved"]["triv_on_circ_p7_is_identity"] is False
    assert out["proved"]["ync_hs_formula_is_identity"] is False
    assert out["proved"]["mixture_meets_floor"] is False
    assert out["proved"]["Nplus_named"] is False
    assert out["status"]["floor"] is False
    assert out["status"]["identities_A_B_C_D"] is True
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15282.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.282"
    assert data["status"]["floor"] is False
    assert data["proved"]["phi_F_ge_6_proved_general"] is False
