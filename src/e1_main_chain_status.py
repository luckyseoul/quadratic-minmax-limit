#!/usr/bin/env python3
"""Original-problem status with separately retained Paley-route diagnostics.

Legacy E1 names describe the stronger Paley gap-two proof architecture.
They neither certify nor veto the original MathOverflow conclusion.  Global
status comes only from the explicit reviewed registry in original_mo_status.
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
)
from e1_gmin_m4_prop15170 import gsum_disj_lb_proved_general  # noqa: E402
from original_mo_status import original_mo_status  # noqa: E402


def four_e1_units_closed() -> dict:
    """Historical Paley-route diagnostic, not global proof acceptance."""
    try:
        from e1_gmin_m4_prop15720 import required_bitight_levels_empty_all_primes
        from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
        from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
        from e1_gmin_m4_prop15276 import (
            lemma_D_existence_written,
            lemma_D_2plane_amplitudes_proved,
        )
        from e1_gmin_m4_prop15764 import minimal_gap4_shell_bridge_closed_general
    except ImportError as exc:  # pragma: no cover
        return {"closed": False, "import_error": str(exc)}
    bitight = bool(required_bitight_levels_empty_all_primes())
    resii = bool(residual_ii_k_ge_4p_ND_closed())
    type_i = bool(type_I_multilevel_bad_case_ND_closed())
    lem_d = bool(lemma_D_existence_written() and lemma_D_2plane_amplitudes_proved())
    shell_bridge = bool(minimal_gap4_shell_bridge_closed_general())
    return {
        "bitight_levels_2_3": bitight,
        "residual_ii_k_ge_4p": resii,
        "type_I_multilevel": type_i,
        "lemma_D": lem_d,
        "minimal_gap4_shell_bridge": shell_bridge,
        "closed": bool(bitight and resii and type_i and lem_d and shell_bridge),
    }


def _props_15167_171_slice(solution: str) -> str:
    """Include Prop 15.167–15.171 writeup sections in soft-close scans (not only solution[:2500])."""
    m = re.search(r"## Prop 15\.167", solution)
    if not m:
        return solution[-8000:]  # fallback: tail often has latest props
    return solution[m.start() :]


def check_docs_L_status() -> dict:
    """Check global conclusions and historical Paley claims independently.

    An open optional route cannot invalidate an independently reviewed global
    proof.  Conversely, route flags alone never license a global conclusion.
    This is a bounded claim scanner, not a verifier of mathematical prose.
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
    global_status = original_mo_status()
    # Fixed-width patterns only (variable-width look-behind is invalid in re).
    existence_patterns = [
        r"\*\*L CLOSED\.\*\*",
        r"\*\*L CLOSED\*\*",
        r"L CLOSED \(via",
        r"(?<![\w-])(?<!Non )Existence\s+(?:CLOSED|PROVED)\b",
        r"E1_closed_general\s*=\s*true\.?\s*L CLOSED",
    ]
    half_value_patterns = [
        r"L=\s*1/2\s*CLOSED",
        r"L=\s*\\?tfrac\{1\}\{2\}\s*CLOSED",
        r"L=\s*½\s*CLOSED",
        r"lim\s*α_n\s*=\s*1/2\s*CLOSED",
        r"Main claim.*L=.*1/2.*CLOSED",
        r"superseded.*L=.*CLOSED",
        r"E\(1\); \$L=\\tfrac12\$",
        r"\*\*Claim:\*\* \*\*\\?\$?L=1/2",
    ]
    nonexistence_patterns = [r"Non[- ]?existence\s+(?:CLOSED|PROVED)"]
    route_patterns = [
        r"\*\*E\(1\) closed\*\*",
        r"E\(1\) closed\s*\\?Rightarrow",
        r"E\(1\)\s*\(CLOSED",
        r"\*\*E\(1\) CLOSED:\*\*",
        r"Closes residual \(i\) of E\(1\) for all primes",
        r"Closes residual \(ii\) of E\(1\) for all primes",
        r"Deep freeness-fail ND \(residual ii\) \| \*\*CLOSED\*\*",
        r"Residual \*\*\(ii\)\*\* ND is \*\*CLOSED\*\*",
        r"Residual \(ii\) full ND closed",
        r"Full residual \(ii\) is \*\*CLOSED\*\*",
    ]
    obsolete_patterns = [
        r"CLOSED for general \$p\$ by Prop 15\.170",
        r"disj: association-scheme min",
        r"association-scheme min \$-12",
        r"Denseness path is \*\*blocked\*\*\s+by residual \*\*\(i\)\*\* only",
    ]

    def first_hit(patterns: list[str], text: str) -> str | None:
        for pattern in patterns:
            match = re.search(pattern, text, re.I | re.M)
            if match:
                return match.group(0)[:80]
        return None

    global_hits = []
    if not global_status["existence_proved"]:
        hit = first_hit(existence_patterns, head)
        if hit:
            global_hits.append(hit)
    if not (global_status["value_proved"] and global_status["limit_value"] == "1/2"):
        hit = first_hit(half_value_patterns, head)
        if hit:
            global_hits.append(hit)
    if not global_status["nonexistence_proved"]:
        hit = first_hit(nonexistence_patterns, head)
        if hit:
            global_hits.append(hit)
    if not global_status["problem_settled"]:
        hit = first_hit([r"\*\*Main Theorem \(limit\)\.\*\*"], head)
        if hit:
            global_hits.append(hit)

    # These legacy closure phrases are interpreted only in the historical
    # Paley writeup/package.  They are not prerequisites for global claims.
    route_hit = None
    if not units.get("closed"):
        route_hit = first_hit(route_patterns, props_tail + "\n" + package)
    soft_hit = global_hits[0] if global_hits else route_hit
    soft = soft_hit is not None

    handoff_open = bool(
        re.search(
            r"L\s*(=\s*lim[^.\n]*)?\s*OPEN|L=\s*\\?lim[^\n]*OPEN|is OPEN|L OPEN",
            handoff[:5000] + status[:2000] + solution[:2500],
            re.I,
        )
    )

    overclaim_hit = first_hit(obsolete_patterns, head)
    overclaim = overclaim_hit is not None
    docs_ok = bool(
        global_status["registry_valid"]
        and (handoff_open or global_status["problem_settled"])
        and not soft
        and not overclaim
    )

    return {
        "HANDOFF_shows_L_OPEN": handoff_open,
        "soft_close_detected": soft,
        "soft_close_hit": soft_hit,
        "global_claim_hits": global_hits,
        "paley_route_claim_hit": route_hit,
        "overclaim_detected": overclaim,
        "overclaim_hit": overclaim_hit,
        "docs_ok": docs_ok,
        "e1_closed_general": e1,
        "four_e1_units": units,
        "original_mo": global_status,
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
    docs = check_docs_L_status()
    global_status = docs["original_mo"]
    L_closed = global_status["existence_proved"]
    L_status = global_status["limit_status"]
    out = {
        "title": "Original MO status with optional historical Paley diagnostics",
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
            "type_I_historical_two_level_classes_closed": p168["proved"].get(
                "type_I_historical_two_level_classes_closed"
            ),
            "type_I_all_classes_closed": p168["proved"].get(
                "type_I_all_classes_closed"
            ),
            "type_I_multilevel_bad_case_closed": p168["proved"].get(
                "type_I_multilevel_bad_case_closed"
            ),
            "deep_tight_empty": p168["proved"].get("deep_tight_empty_p_ge_5"),
            "deep_auto_freeness": p168["proved"].get(
                "deep_auto_freeness_k_le_3p_minus_2"
            ),
            "deep_all_ND_closed": p168["proved"].get("deep_all_ND_closed"),
        },
        "Main_closed": global_status["problem_settled"],
        "L_closed": L_closed,
        "L_status": L_status,
        "L_value": global_status["limit_value"],
        "existence_proved": global_status["existence_proved"],
        "nonexistence_proved": global_status["nonexistence_proved"],
        "value_proved": global_status["value_proved"],
        "original_mo": global_status,
        "optional_route_diagnostics": {"paley_gap_two": docs["four_e1_units"]},
        "writeup_L_closed": bool(L_closed and docs["docs_ok"]),
        "docs": docs,
        "rule": (
            "Global completion is determined by the independent reviewed original-MO "
            "proof registry. Legacy E1/Paley flags are optional route diagnostics, "
            "not necessary or sufficient machine acceptance for the original question. "
            "Existence, an identified value, and nonexistence are separate conclusions."
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
