"""AΓL square-multiplication Aut of Paley conference is reducible on V_+."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


def _load_script():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "paley_Vplus_irrep.py"
    )
    spec = importlib.util.spec_from_file_location("paley_Vplus_irrep", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_affine_squares_commute_but_commutant_dim_two_at_p3_p5():
    mod = _load_script()
    for p in (3, 5):
        rec = mod.commutant_dim(p)
        assert rec["generator_commute_err_max"] < 1e-9
        assert rec["dim_Vplus"] == rec["n"] // 2
        assert rec["commutant_dim"] == 2
        assert rec["irreducible"] is False
