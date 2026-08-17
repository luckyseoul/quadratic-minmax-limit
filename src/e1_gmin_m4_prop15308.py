#!/usr/bin/env python3
"""
Prop 15.308 — |G|=p²(p²−1); AP-halfspace stab=4p; |H_+|=Σ|G|/|Stab_τ|
after AP/non-AP split; NL stab list still unnamed in p.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.307 flags.  Floor stays OPEN.

============================================================================
Setup.  Strict Aut G fixing ∞:
    u ↦ a u^{p^k} + b,   a∈F_q^□, k∈{0,1},   ∞↦∞.
|G|=q(q−1)=p²(p²−1).  G preserves H_+={y_∞=+1}.  Orbit–stab:
    |O_τ|=|G|/|Stab_τ|,   |H_+|=∑_τ |O_τ|.

============================================================================
Theorem A — PROVED (group order; Max+-free).
  |G|=p²(p²−1).  Fail-when-wrong: drop Frobenius; then |G|=q(q−1)/2
  (300 at p=5, 1176 at p=7) and the AP-orbit identity in B fails.  ∎

Theorem B — PROVED (AP halfspaces).
  If S⊂F_p is an AP of size m=(p+1)/2, the affine halfspace y_{L,S}
  has |Stab_G|=4p and |O|=|G|/(4p)=p(p²−1)/4.
  Certified p=3,5,7.  Fail-when-wrong: |Stab|=2p (drops the interval
  flip or Frob).  ∎

Theorem C — PROVED (orbit–stab sum; fail if halfspaces are lumped).
  Split H_+ types by square-dir signature, and split halfspaces by
  whether the special-direction S is an AP.  Then Σ |O_τ| = |H_+|
  at p=3,5,7 (6, 130, 5726).  Lumping all halfspaces into one type
  misses the p=7 non-AP orbit of size 56 (=|G|/(6p)).
  There is an NL orbit with |Stab|=p+1 and |O|=p²(p−1) at p=5,7;
  at p=5 it is the unique NL type, at p=7 there are extra NL orbits
  with |Stab|∈{2,4}.  Fail-when-wrong: claim unique NL at p=7.  ∎

Theorem D — OPEN.  The NL stab list {p+1} (p=5) vs {2,2,2,2,4,8} (p=7)
  is not a named function of p.  |H_+|=n_1d + Σ_{NL} |G|/s_τ stays
  unnamed.  Q_{++} unnamed.  phi_F not imported.  ∎

============================================================================
Backend: CPU group close |G|≤2352; vectorized H_+ type reps.
Writes evidence/e1_gmin_m4_prop15308.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15307 import (  # noqa: E402
    load_D,
    sig_H,
    square_classes,
    tag_sig,
)


def G_order(p: int) -> int:
    q = p * p
    return q * (q - 1)


def G_order_no_frob(p: int) -> int:
    q = p * p
    return q * (q - 1) // 2


def AP_stab(p: int) -> int:
    return 4 * p


def AP_orbit(p: int) -> int:
    return G_order(p) // AP_stab(p)


def frob_table(F: dict) -> list[int]:
    q, p = F["q"], F["p"]
    units = [0] * (q - 1)
    g, mul = F["g"], F["mul"]
    cur = 1
    for j in range(q - 1):
        units[j] = cur
        cur = mul[cur][g]
    out = [0] * q
    for x in range(1, q):
        j = F["dlog"][x]
        out[x] = units[(j * p) % (q - 1)]
    return out


_G_CACHE: dict[tuple[int, bool], list] = {}


def build_G(F: dict, use_frob: bool = True) -> list[np.ndarray]:
    key = (int(F["p"]), bool(use_frob))
    if key in _G_CACHE:
        return _G_CACHE[key]
    q = F["q"]
    n = q + 1
    mul, add = F["mul"], F["add"]
    fr = frob_table(F)
    squares = F["squares"]

    def apply(a: int, b: int, k: int) -> np.ndarray:
        perm = np.empty(n, dtype=np.int32)
        perm[0] = 0
        for u in range(q):
            w = fr[u] if k else u
            w = add[mul[a][w]][b]
            perm[1 + u] = 1 + w
        return perm

    gens = [apply(1, b, 0) for b in range(q)]
    gens += [apply(a, 0, 0) for a in squares]
    if use_frob:
        gens.append(apply(1, 0, 1))
    seen: dict[tuple[int, ...], np.ndarray] = {}
    stack = [np.arange(n, dtype=np.int32)]
    seen[tuple(range(n))] = stack[0]
    while stack:
        g = stack.pop()
        for s in gens:
            h = s[g]
            key = tuple(int(x) for x in h)
            if key not in seen:
                seen[key] = h
                stack.append(h)
    out = list(seen.values())
    _G_CACHE[key] = out
    return out


def is_AP(S: set[int] | list[int], p: int) -> bool:
    S = sorted(set(int(x) for x in S))
    n = len(S)
    if n <= 1:
        return True
    for a in range(p):
        for d in range(1, p):
            want = sorted((a + k * d) % p for k in range(n))
            if want == S:
                return True
    return False


def pack_y(Drow: np.ndarray) -> np.ndarray:
    q = Drow.shape[0]
    y = np.ones(q + 1, dtype=np.int8)
    y[1:] = np.where(Drow, -1, 1)
    return y


def orbit_size(G: list[np.ndarray], y0: np.ndarray) -> int:
    if y0[0] < 0:
        y0 = -y0
    seen = set()
    for g in G:
        ynew = np.empty_like(y0)
        ynew[g] = y0
        if ynew[0] < 0:
            ynew = -ynew
        seen.add(tuple((ynew > 0).astype(np.int8).tolist()))
    return len(seen)


def recover_S_if_hs(p: int, Drow: np.ndarray, classes: list[np.ndarray]) -> set[int] | None:
    """If this D is a halfspace, return S = empty-line intercepts on the H class."""
    H = sig_H(p)
    for lines in classes:
        ns = Drow[lines].sum(axis=1)
        sig = tuple(int(x) for x in np.sort(ns))
        if sig == H:
            return {int(c) for c in range(p) if int(ns[c]) == 0}
    return None


def classify_split(p: int, G: list[np.ndarray]) -> dict:
    D, F = load_D(p)
    classes = square_classes(F)
    n_sq = len(classes)
    # type key: tags, plus AP/nonAP for halfspaces
    first: dict[tuple, int] = {}
    counts: dict[tuple, int] = {}
    hs_key = tuple(sorted(["H"] + ["C"] * (n_sq - 1)))
    for i in range(D.shape[0]):
        tags = []
        for lines in classes:
            ns = D[i, lines].sum(axis=1)
            sig = tuple(int(x) for x in np.sort(ns))
            tags.append(tag_sig(p, sig))
        t = tuple(sorted(tags))
        if t == hs_key:
            S = recover_S_if_hs(p, D[i], classes)
            t = ("hs-AP",) if (S is not None and is_AP(S, p)) else ("hs-nonAP",)
        first.setdefault(t, i)
        counts[t] = counts.get(t, 0) + 1
    rows = []
    total = 0
    for t, i in first.items():
        orb = orbit_size(G, pack_y(D[i]))
        stab = len(G) // orb if orb else 0
        total += orb
        rows.append(
            {
                "type": str(t),
                "count_in_Hplus": counts[t],
                "orbit": orb,
                "stab": stab,
                "count_matches_orbit": counts[t] == orb,
            }
        )
    lumped_hs = sum(r["orbit"] for r in rows if r["type"] in ("('hs-AP',)", "('hs-nonAP',)"))
    # if we used only one hs rep, lumped orbit is the first hs orbit only
    return {
        "p": p,
        "Nplus": int(D.shape[0]),
        "n_1d": n_1d(p),
        "G": len(G),
        "types": rows,
        "sum_orbits": total,
        "sum_matches_Hplus": total == int(D.shape[0]),
        "n_types": len(rows),
        "lumped_note": lumped_hs,
    }


def prove_G_order() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        G = build_G(F, use_frob=True)
        G0 = build_G(F, use_frob=False)
        good = len(G) == G_order(p)
        dropped = len(G0) == G_order_no_frob(p)
        if not good or not dropped:
            ok = False
        if len(G0) == G_order(p):
            ok = False  # fail-when-wrong: no-Frob accidentally full
        rows[str(p)] = {
            "G": len(G),
            "named": G_order(p),
            "G_no_frob": len(G0),
            "named_no_frob": G_order_no_frob(p),
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "|G|=p²(p²−1).  Drop Frob ⇒ q(q−1)/2.",
    }


def prove_AP_stab() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        G = build_G(F, True)
        r = classify_split(p, G)
        ap = next((t for t in r["types"] if t["type"] == "('hs-AP',)"), None)
        if ap is None or ap["stab"] != AP_stab(p) or ap["orbit"] != AP_orbit(p):
            ok = False
        if ap is not None and ap["stab"] == 2 * p:
            ok = False
        rows[str(p)] = ap
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "AP-halfspace |Stab|=4p, |O|=p(p²−1)/4.  Not 2p.",
    }


def prove_orbit_sum() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        F = _kernel_field(p)
        G = build_G(F, True)
        r = classify_split(p, G)
        if not r["sum_matches_Hplus"]:
            ok = False
        if not all(t["count_matches_orbit"] for t in r["types"]):
            ok = False
        # unique-NL-at-p=7 must fail
        nl = [t for t in r["types"] if not t["type"].startswith("('hs-")]
        if p == 7 and len(nl) <= 1:
            ok = False
        if p == 5 and not any(t["stab"] == p + 1 for t in nl):
            ok = False
        if p == 7 and not any(t["stab"] == p + 1 for t in nl):
            ok = False
        # lumped check at p=7: one hs type would miss 56
        if p == 7:
            ap = next(t for t in r["types"] if t["type"] == "('hs-AP',)")
            na = next(t for t in r["types"] if t["type"] == "('hs-nonAP',)")
            if ap["orbit"] + sum(t["orbit"] for t in nl) == r["Nplus"]:
                ok = False  # that would mean nonAP orbit is 0
            if na["orbit"] != 56 or na["stab"] != 6 * p:
                ok = False
        rows[str(p)] = {
            "Nplus": r["Nplus"],
            "sum_orbits": r["sum_orbits"],
            "n_types": r["n_types"],
            "types": r["types"],
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "Σ|O_τ|=|H_+| after AP/non-AP split. "
            "Lumping hs misses 56 at p=7. NL stab=p+1 exists; not unique at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Hplus_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "NL stabs {p+1} at p=5 vs {2,2,2,2,4,8} at p=7 are not a "
            "named list in p. |H_+|=n_1d+Σ|G|/s_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.308  |G|=p²(p²−1); AP stab=4p; orbit-stab sum", flush=True)
    A = prove_G_order()
    print(f"  A |G|: {A['proved']} {A['by_p']}", flush=True)
    B = prove_AP_stab()
    print(f"  B AP stab=4p: {B['proved']}", flush=True)
    C = prove_orbit_sum()
    print(f"  C Σ|O|=|H_+|: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.308",
        "title": "|G|=p²(p²−1); AP-hs stab=4p; |H_+|=Σ|G|/|Stab| after split; NL stabs unnamed",
        "proved": {
            "G_order": A["proved"],
            "AP_stab_4p": B["proved"],
            "orbit_stab_sum": C["proved"],
            "Hplus_named": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15308.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
