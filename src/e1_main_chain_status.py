#!/usr/bin/env python3
"""
Public-doc honesty check for E(1) / L.

    Writeup may assert L=1/2 only after the corrected GOAL.md units are
    actually imported (not the legacy e1_closed_general old AND, the retired
    Gsum hinge, or the retracted 15.167 spectral arrow). Soft-close banned.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_bitight_chain import run_bitight_chain  # noqa: E402
from e1_gmin_m4_prop15168 import (  # noqa: E402
    e1_closed_general,
    e1_residual_open,
    main as prop15168_main,
    main_L_from_e1,
)
from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general  # noqa: E402


def four_e1_units_closed() -> dict:
    """Corrected GOAL.md acceptance, not the old e1 AND."""
    try:
        from e1_gmin_m4_prop15720 import required_bitight_levels_empty_all_primes
        from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
        from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
        from e1_gmin_m4_prop15276 import (
            lemma_D_existence_written,
            lemma_D_2plane_amplitudes_proved,
        )
    except ImportError as exc:  # pragma: no cover
        return {"closed": False, "import_error": str(exc)}
    bitight = bool(required_bitight_levels_empty_all_primes())
    resii = bool(residual_ii_k_ge_4p_ND_closed())
    type_i = bool(type_I_multilevel_bad_case_ND_closed())
    lem_d = bool(lemma_D_existence_written() and lemma_D_2plane_amplitudes_proved())
    return {
        "bitight_levels_2_3": bitight,
        "residual_ii_k_ge_4p": resii,
        "type_I_multilevel": type_i,
        "lemma_D": lem_d,
        "closed": bool(bitight and resii and type_i and lem_d),
    }


def _props_15167_171_slice(solution: str) -> str:
    """Include Prop 15.167–15.171 writeup sections in soft-close scans (not only solution[:2500])."""
    m = re.search(r"## Prop 15\.167", solution)
    if not m:
        return solution[-8000:]  # fallback: tail often has latest props
    return solution[m.start() :]


def check_docs_L_status() -> dict:
    """
    Docs OK iff status asserts L OPEN (or true L closed with proved hinge)
    and no soft-close 'L CLOSED' / residual closed claims while the four-unit
    acceptance gate is open.
    Scans HANDOFF head, STATUS, solution top, **and** Props 15.167–171 body.
    """
    handoff = (ROOT / "HANDOFF.md").read_text(encoding="utf-8", errors="replace")
    solution = (ROOT / "solution.md").read_text(encoding="utf-8", errors="replace")
    status = (ROOT / "STATUS.md").read_text(encoding="utf-8", errors="replace")
    package = ""
    pkg_path = ROOT / "evidence" / "share" / "denseness_path_package.md"
    if pkg_path.exists():
        package = pkg_path.read_text(encoding="utf-8", errors="replace")

    props_tail = _props_15167_171_slice(solution)
    head = (
        handoff[:5000]
        + "\n"
        + status[:3000]
        + "\n"
        + solution[:2500]
        + "\n"
        + props_tail
        + "\n"
        + package
    )

    units = four_e1_units_closed()
    e1 = bool(e1_closed_general())
    # Fixed-width patterns only (variable-width look-behind is invalid in re).
    soft_patterns = [
        r"\*\*L CLOSED\.\*\*",
        r"\*\*L CLOSED\*\*",
        r"L CLOSED \(via E1",
        r"L CLOSED \(via",
        r"E1_closed_general\s*=\s*true\.?\s*L CLOSED",
        r"E1_closed_general=true\. L CLOSED",
        r"L=\s*\\?tfrac\{1\}\{2\}\s*CLOSED",
        r"L=\s*½\s*CLOSED",
        r"lim\s*α_n\s*=\s*1/2\s*CLOSED",
        r"Main claim.*L=.*1/2.*CLOSED",
        r"superseded.*L=.*CLOSED",
        r"\*\*E\(1\) closed\*\*",
        r"E\(1\) closed\s*\\?Rightarrow",
        r"Closes residual \(i\) of E\(1\) for all primes",
        r"Closes residual \(ii\) of E\(1\) for all primes",
        r"CLOSED for general \$p\$ by Prop 15\.170",
        r"disj: association-scheme min",
        r"association-scheme min \$-12",
    ]
    soft = False
    soft_hit = None
    if not units.get("closed"):
        for pat in soft_patterns:
            m = re.search(pat, head, re.I | re.M)
            if m:
                soft = True
                soft_hit = m.group(0)[:80]
                break
        if re.search(r"E\(1\)\s*\(CLOSED", head, re.I):
            soft = True
            soft_hit = soft_hit or "E(1) (CLOSED"

    handoff_open = bool(
        re.search(
            r"L\s*(=\s*lim[^.\n]*)?\s*OPEN|L=\s*\\?lim[^\n]*OPEN|is OPEN|L OPEN",
            handoff[:5000] + status[:2000] + solution[:2500],
            re.I,
        )
    )

    overclaim_patterns = [
        r"\*\*Main Theorem \(limit\)\.\*\*",
        r"E\(1\); \$L=\\tfrac12\$",
        r"\*\*E\(1\) CLOSED:\*\*",
        r"\*\*Claim:\*\* \*\*\\?\$?L=1/2",
        r"Denseness path is \*\*blocked\*\*\s+by residual \*\*\(i\)\*\* only",
        r"Deep freeness-fail ND \(residual ii\) \| \*\*CLOSED\*\*",
        r"Residual \*\*\(ii\)\*\* ND is \*\*CLOSED\*\*",
        r"Residual \(ii\) full ND closed",
        r"Full residual \(ii\) is \*\*CLOSED\*\*",
    ]
    overclaim = False
    overclaim_hit = None
    if not units.get("closed"):
        for pat in overclaim_patterns:
            m = re.search(pat, head, re.I | re.M)
            if m:
                overclaim = True
                overclaim_hit = m.group(0)[:80]
                break

    if not units.get("closed"):
        docs_ok = handoff_open and not soft and not overclaim
    else:
        docs_ok = not soft

    return {
        "HANDOFF_shows_L_OPEN": handoff_open,
        "soft_close_detected": soft,
        "soft_close_hit": soft_hit,
        "overclaim_detected": overclaim,
        "overclaim_hit": overclaim_hit,
        "docs_ok": docs_ok,
        "e1_closed_general": e1,
        "four_e1_units": units,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "scanned_props_15167_171": True,
        "scanned_denseness_package_full": True,
    }


def run_main_chain() -> dict:
    bt = run_bitight_chain()
    p168 = prop15168_main()
    e1_closed = bool(e1_closed_general())
    residual_e1 = e1_residual_open()
    bi = bool(bt["bi_tight_required_levels_empty_for_all_p_ge_5"])
    residual = False
    Lwire = main_L_from_e1(e1_closed, bi)
    docs = check_docs_L_status()
    L_closed = bool(docs["four_e1_units"]["closed"])
    out = {
        "title": "Main/E(1) chain status (corrected discrete bi-tight gate)",
        "bi_tight_required_levels_empty_for_all_p_ge_5": bi,
        "bi_tight_all_levels_empty": False,
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
        "writeup_L_closed": bool(docs["four_e1_units"]["closed"]),
        "docs": docs,
        "rule": (
            "Public writeup may assert L=1/2 only after the valid gates: "
            "bi-tight levels 2/3, residual_ii k≥4p, multi-level Type I, Lemma D. "
            "QVAR/R1 and the spectral floor are no longer acceptance units. "
            "Live e1_closed_general is a separate wiring fact. Soft-close banned (F3)."
        ),
    }
    return out


def main() -> dict:
    out = run_main_chain()
    path = ROOT / "evidence" / "e1_main_chain_status.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("=== main/E(1) chain ===")
    print(f"bi_tight_required_levels={out['bi_tight_required_levels_empty_for_all_p_ge_5']}")
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
