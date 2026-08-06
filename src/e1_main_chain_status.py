#!/usr/bin/env python3
"""
Main / E(1) wiring check after Props 15.167–15.171.

L closed only if bi-tight AND full E(1) (m_n≥Φ−2 all p≥5).
Honest 2026-08-06: E1 open until gsum_disj_lb_proved_general (15.170 hinge).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_bitight_chain import run_bitight_chain  # noqa: E402
from e1_gmin_m4_prop15167 import prove_open as open_15167  # noqa: E402
from e1_gmin_m4_prop15168 import (  # noqa: E402
    e1_closed_general,
    e1_residual_open,
    main as prop15168_main,
    main_L_from_e1,
)
from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general  # noqa: E402


def check_docs_L_status() -> dict:
    """
    Docs OK iff top-of-repo status asserts L OPEN (or true L closed with proved hinge)
    and no bare soft-close 'L CLOSED' while E1 is open.
    """
    handoff = (ROOT / "HANDOFF.md").read_text(encoding="utf-8", errors="replace")
    solution = (ROOT / "solution.md").read_text(encoding="utf-8", errors="replace")
    status = (ROOT / "STATUS.md").read_text(encoding="utf-8", errors="replace")
    head = handoff[:5000] + "\n" + status[:3000] + "\n" + solution[:2500]

    e1 = bool(e1_closed_general())
    # Soft-close: assert L CLOSED while chain open, or bare **L CLOSED** in head when E1 open
    soft_patterns = [
        r"\*\*L CLOSED\.\*\*",
        r"\*\*L CLOSED\*\*",
        r"L CLOSED\.\s*$",
        r"L=\s*\\?tfrac\{1\}\{2\}\s*CLOSED",
        r"L=\s*½\s*CLOSED",
        r"lim\s*α_n\s*=\s*1/2\s*CLOSED",
        r"Main claim.*L=.*1/2.*CLOSED",
        r"superseded.*L=.*CLOSED",
    ]
    soft = False
    if not e1:
        for pat in soft_patterns:
            if re.search(pat, head, re.I | re.M):
                soft = True
                break
        # P0-style false closed
        if re.search(r"E\(1\)\s*\(CLOSED", head, re.I):
            soft = True

    handoff_open = bool(
        re.search(r"L\s*(=\s*lim[^.\n]*)?\s*OPEN|L=\s*\\?lim[^\n]*OPEN|is OPEN", head, re.I)
    ) or ("L OPEN" in head[:4000])

    docs_ok = bool(handoff_open and not soft and not e1) or bool(e1 and gsum_disj_lb_proved_general())
    # When E1 open, require OPEN assertion and no soft-close
    if not e1:
        docs_ok = handoff_open and not soft

    return {
        "HANDOFF_shows_L_OPEN": handoff_open,
        "soft_close_detected": soft,
        "docs_ok": docs_ok,
        "e1_closed_general": e1,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
    }


def run_main_chain() -> dict:
    bt = run_bitight_chain()
    o = open_15167()
    p168 = prop15168_main()
    e1_closed = bool(e1_closed_general())
    residual_e1 = e1_residual_open()
    bi = bool(bt["bi_tight_empty_for_all_p_ge_5"])
    residual = bool(o["residual_closed_general"])
    Lwire = main_L_from_e1(e1_closed, bi)
    L_closed = bool(Lwire["L_closed"])
    docs = check_docs_L_status()
    out = {
        "title": "Main/E(1) chain status after 15.167–15.171 (honest Gsum hinge)",
        "bi_tight_empty_for_all_p_ge_5": bi,
        "residual_closed_general": residual,
        "bitight_bypass_residual": True,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "E1_closed": e1_closed,
        "E1_open_residual": residual_e1["open"],
        "E1_structure_15168": {
            "type_I_fail_k_2p_minus_1_ND": p168["proved"].get(
                "type_I_fail_k_2p_minus_1_ND"
            ),
            "type_I_all_classes_closed": p168["proved"].get(
                "type_I_all_classes_closed"
            ),
            "deep_tight_empty": p168["proved"].get("deep_tight_empty_p_ge_5"),
            "deep_auto_freeness": p168["proved"].get(
                "deep_auto_freeness_k_le_3p_minus_2"
            ),
            "deep_all_ND_closed": p168["proved"].get("deep_all_ND_closed"),
        },
        "Main_closed": L_closed,
        "L_closed": L_closed,
        "L_status": "CLOSED" if L_closed else "OPEN",
        "docs": docs,
        "rule": (
            "L closed iff bi-tight (15.167) ∧ full E(1). "
            "E1 requires gsum_disj_lb_proved_general (15.170 hinge; 15.158). "
            "Soft-close banned (F3)."
        ),
    }
    return out


def main() -> dict:
    out = run_main_chain()
    path = ROOT / "evidence" / "e1_main_chain_status.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("=== main/E(1) chain ===")
    print(f"bi_tight={out['bi_tight_empty_for_all_p_ge_5']}")
    print(f"gsum_disj_lb_proved={out['gsum_disj_lb_proved_general']}")
    print(f"E1_closed={out['E1_closed']}")
    print(f"E1_open={out['E1_open_residual']}")
    print(f"Main_closed={out['Main_closed']}")
    print(f"L_status={out['L_status']}")
    print(f"docs={out['docs']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
