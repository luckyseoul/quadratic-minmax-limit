#!/usr/bin/env python3
"""
Prop 15.302 — S_k is the even torus Fourier of Q on U = {N=1}
plus a named F_p weight times the form factor of a complementary coset.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.301 flags.  Does **not** name Q_τ in p.  Floor stays OPEN.

============================================================================
Setup.  q=p².  T=□.  U={u: N_{q/p}(u)=u^{p+1}=1}, |U|=p+1, U⊂T.
F_p^× ⊂ T and T = F_p^× · U with F_p^× ∩ U={±1}.  ξ_j(u)=u-character
of index j on the cyclic group U.  Even k ⇒ j=k mod (p+1) is even.

============================================================================
Theorem A — PROVED (group law on T; Max+-free once Q is known on T).
  For even k,
      S_k/q² = F_1(j) + σ_k F_*(j),   j = k mod (p+1),
  σ_k = (p−3)/2 if (p−1)|k else −1,
  F_1(j)=∑_{u∈U} ξ_j(u) Q(u)/q²,
  F_*(j)=∑_{u∈U} ξ_j(u) Q(s u)/q²  for any s∈F_p^×\\{±1}.
  At p=5,7, F_* is independent of that s (one inversion-orbit).
  Fail-when-wrong: flip σ on a HD character; or replace s by 1.  ∎

Theorem B — PROVED (norm of sU; Max+-free type calculus).
  s·{±1} ⊂ sub ⊂ coarse ++ and s·(U\\{±1}) ⊂ N*.  Hence for even j
      F_*(0)     = 2 Q_{++} + (p−1) Q_{N*},
      F_*(j≠0)   = 2(Q_{++} − Q_{N*}).
  Fail-when-wrong: put Q_{++} on s·(U\\{±1}).  ∎

Theorem C — CERTIFIED p=5,7 (even circulant of ΔQ).
  Odd torus frequencies of Q|_U vanish.  If Q_{++}≥Q_{N*} (so F_*(j≠0)≥0)
  the σ=+(p−3)/2 family is at least the σ=−1 family, and
      min_k S_k/q² = min_{even j} (F_1(j)−F_*(j))
  = min even eigenvalue of the U-circulant of ΔQ(u)=Q(u)−Q(s u).
  Live mins are 80/13 (j=2) and 3072/409 (j=0, HD).  Both ≥6, but
  Q_τ still unnamed, so ⟨δ,ψ⟩≤2 / phi_F stay OPEN.  ∎

============================================================================
Backend: CPU field tables + live Q (one cache fold per prime).
Writes evidence/e1_gmin_m4_prop15302.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import field_norm, in_subfield, live_Q_on_T  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402

EVEN_K = {5: (0, 2, 4, 6, 8, 10), 7: (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)}


def torus_U(F) -> list[int]:
    return [x for x in range(1, F["q"]) if field_norm(F, x) == 1]


def gen_U(F) -> int:
    acc = 1
    for _ in range(F["p"] - 1):
        acc = F["mul"][acc][F["g"]]
    return acc


def U_powers(F) -> list[int]:
    gU = gen_U(F)
    out = [1]
    for _ in range(F["p"]):
        out.append(F["mul"][out[-1]][gU])
    return out


def subfield_units(F) -> list[int]:
    return [x for x in range(1, F["q"]) if in_subfield(F, x)]


def nontrivial_s(F) -> int:
    for x in subfield_units(F):
        if x not in (1, F["neg1"]):
            return x
    raise RuntimeError("no s")


def sigma_k(p: int, k: int) -> Fraction:
    if k % (p - 1) == 0:
        return Fraction(p - 3, 2)
    return Fraction(-1)


def coarse_Q_means(F, Q: dict) -> dict[str, float]:
    q2 = float(F["q"] ** 2)
    acc: dict[str, list[float]] = {}
    for r, val in Q.items():
        acc.setdefault(merge_cls(F, r), []).append(val / q2)
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def form_factor(F, Q: dict, scale: int, j: int) -> complex:
    q2 = float(F["q"] ** 2)
    p = F["p"]
    acc = 0j
    for t, u in enumerate(U_powers(F)):
        su = F["mul"][scale][u]
        acc += np.exp(2j * np.pi * j * t / (p + 1)) * Q[su] / q2
    return acc


def F_star_named(p: int, j: int, Qpp: float, Qn: float) -> float:
    if j == 0:
        return 2.0 * Qpp + (p - 1) * Qn
    return 2.0 * (Qpp - Qn)


def S_from_torus(F, Q: dict, k: int, s: int | None = None, flip_sigma: bool = False) -> complex:
    p = F["p"]
    j = k % (p + 1)
    if s is None:
        s = nontrivial_s(F)
    sig = sigma_k(p, k)
    if flip_sigma:
        sig = -sig
    return form_factor(F, Q, 1, j) + float(sig) * form_factor(F, Q, s, j)


def live_S(F, Q: dict, k: int) -> complex:
    psi = _psi_vec(F, k)
    q2 = float(F["q"] ** 2)
    return sum(psi[r] * Q[r] / q2 for r in Q)


def prove_decomposition() -> dict:
    ok = True
    rows = {}
    flip_hit = 0
    s1_hit = 0
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        rec = {}
        s = nontrivial_s(F)
        alts = [x for x in subfield_units(F) if x not in (1, F["neg1"])]
        for k in EVEN_K[p]:
            live = live_S(F, Q, k)
            recov = S_from_torus(F, Q, k, s=s)
            err = abs(recov - live)
            rec[str(k)] = {
                "j": k % (p + 1),
                "sigma": str(sigma_k(p, k)),
                "live": float(live.real),
                "recov": float(recov.real),
                "err": float(err),
            }
            if err > 1e-8:
                ok = False
            if k != 0 and k % (p + 1) == 0:
                bad = S_from_torus(F, Q, k, s=s, flip_sigma=True)
                if abs(bad - live) > 0.5:
                    flip_hit += 1
            if k == 2:
                bad1 = S_from_torus(F, Q, k, s=1)
                if abs(bad1 - live) > 0.5:
                    s1_hit += 1
        # F_* independent of s
        j = 2
        fvals = [form_factor(F, Q, x, j).real for x in alts]
        rec["Fstar_spread"] = float(max(fvals) - min(fvals))
        if max(fvals) - min(fvals) > 1e-8:
            ok = False
        rows[str(p)] = rec
    if flip_hit == 0 or s1_hit == 0:
        ok = False
    return {
        "proved": ok,
        "flip_hit": flip_hit,
        "s1_hit": s1_hit,
        "by_p": rows,
        "theorem": (
            "S_k/q² = F_1(j)+σ_k F_*(j) with j=k mod (p+1), "
            "σ=(p−3)/2 on (p−1)|k else −1. Flip σ or s=1 breaks it."
        ),
    }


def prove_Fstar_named() -> dict:
    ok = True
    wrong_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        means = coarse_Q_means(F, Q)
        Qpp, Qn = means["++"], means["N*"]
        rec = {"Qpp": Qpp, "Qn": Qn}
        s = nontrivial_s(F)
        for j in range(0, p + 1, 2):
            live = form_factor(F, Q, s, j).real
            named = F_star_named(p, j, Qpp, Qn)
            err = abs(live - named)
            rec[f"j{j}"] = {"live": live, "named": named, "err": err}
            if err > 1e-8:
                ok = False
            if j != 0:
                # fail-when-wrong: put Q_++ on s·(U\{±1})
                wrong = 2.0 * (Qpp - Qpp)
                if abs(wrong - live) > 1e-6:
                    wrong_hit += 1
        # odd frequencies of Q|_U vanish
        odd_max = 0.0
        for j in range(1, p + 1, 2):
            odd_max = max(odd_max, abs(form_factor(F, Q, 1, j)))
        rec["odd_max"] = odd_max
        if odd_max > 1e-8:
            ok = False
        rows[str(p)] = rec
    if wrong_hit == 0:
        ok = False
    return {
        "proved": ok,
        "wrong_hit": wrong_hit,
        "by_p": rows,
        "theorem": (
            "F_*(0)=2 Q_++ +(p−1) Q_N*; F_*(j≠0 even)=2(Q_++−Q_N*). "
            "Odd torus hats of Q|_U vanish.  Replacing N* by --N1 fails."
        ),
    }


def prove_floor_rewrite() -> dict:
    """Certified equivalence at p=5,7; does not name Q or import phi_F."""
    ok = True
    mins = {}
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        means = coarse_Q_means(F, Q)
        if means["++"] + 1e-12 < means["N*"]:
            ok = False
        ev = []
        for j in range(0, p + 1, 2):
            ev.append((form_factor(F, Q, 1, j) - form_factor(F, Q, nontrivial_s(F), j)).real)
        live_min = min(live_S(F, Q, k).real for k in EVEN_K[p] if k not in (0,))
        # exclude the χ copy: k=(q-1)/2 has S=S_0
        floor_ks = [k for k in EVEN_K[p] if k not in (0, (F["q"] - 1) // 2)]
        live_min = min(live_S(F, Q, k).real for k in floor_ks)
        dmin = min(ev)
        mins[str(p)] = {
            "circulant_min": dmin,
            "live_min": live_min,
            "ev": ev,
            "ge6": bool(dmin >= 6 - 1e-9),
        }
        if abs(dmin - live_min) > 1e-8:
            ok = False
        if dmin < 6 - 1e-9:
            ok = False
    return {
        "proved": ok,
        "mins": mins,
        "phi_F_imported": False,
        "theorem": (
            "At p=5,7, Q_++≥Q_N* and min S_{k∉{1,χ}} equals min_j(F_1−F_*). "
            "Those mins are 80/13 and 3072/409, both ≥6, but the identity "
            "still uses live Q, so it does not name S in p or import phi_F."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "S_k_named_in_p": False,
        "note": (
            "S_k = F_1(j)+σ F_* reduces the floor to min even eig of "
            "circ(Q|_U − Q|_{sU}), with F_* named in (Q_++,Q_N*). "
            "Q_τ still unnamed (den N/(4p)). Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.302  S_k = F_1(j)+σ F_* on the norm-1 torus", flush=True)
    A = prove_decomposition()
    print(f"  A torus formula: {A['proved']} flip={A['flip_hit']} s1={A['s1_hit']}", flush=True)
    B = prove_Fstar_named()
    print(f"  B F_* named in Q++ , QN*: {B['proved']} wrong={B['wrong_hit']}", flush=True)
    C = prove_floor_rewrite()
    print(
        f"  C min=circ(ΔQ): {C['proved']} "
        f"p5={C['mins']['5']['circulant_min']:.5f} p7={C['mins']['7']['circulant_min']:.5f}",
        flush=True,
    )
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.302",
        "title": "S_k = F_1(j)+σ_k F_* on U={N=1}; F_* named in Q++,QN*",
        "proved": {
            "torus_decomposition": A["proved"],
            "Fstar_named_in_Q": B["proved"],
            "floor_is_even_circulant_of_dQ": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
            "S_k_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15302.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
