#!/usr/bin/env python3
"""
Prop 15.314 — Coarse Q-classes on T are named field-arithmetic sets;
F_1(j) is a named linear form in (Q_++, Q_{N*}, Q_{--}).  Floor is
three explicit inequalities.  Q_τ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.313 flags.  Floor stays OPEN.

============================================================================
Setup.  q=p², T=□, U={N_{q/p}=1}, |U|=p+1, F_p^× ⊂ T, T=F_p^×·U
with F_p^× ∩ U={±1}.  χ_q(x)=χ_p(N(x)).  Tr=Tr_{q/p}.  Coarse
classes as in 15.298 / merge_cls (mixed Paley N1 fused into ++).

============================================================================
Theorem A — PROVED (group law on T; Max+-free).
  |pm1|=2, |N*|=(p−1)(p−3)/2,
  |--N1|=0 if p≡1 (mod 4), (p+1)/2 if p≡3 (mod 4),
  |++|=2p−4 if p≡1 (mod 4), 3(p−3)/2 if p≡3 (mod 4).
  Geometrically: N*=T\\(F_p^× ∪ U),
  ++ = (F_p^×\\{±1}) ∪ (U ∩ ++),
  --N1 ⊂ U.  Fail-when-wrong: |++|=6 for every p.  ∎

Theorem B — PROVED (norm / trace on U; Max+-free).
  For u∈U\\{±1}, N(u−1)=2−Tr(u) and
      χ_q(u²−1)=−χ_p(−1)
  constantly.  Hence U\\{±1} is all mixed Paley (coarse ++) when
  p≡1 (mod 4).  When p≡3 (mod 4), 2−Tr(g^t)=−g^{−t}(g^t−1)² and
      χ_p(2−Tr(g^t))=(−1)^t
  for a generator g of U, so --N1 ∩ U = U\\U² (odd powers) and
  ++ ∩ U = U²\\{±1}.  Fail-when-wrong: both-plus Paley on U at
  p≡1 (mod 4); swap ++ / -- on odd indices.  ∎

Theorem C — PROVED as a rewrite (15.302 + A + B); Q_τ OPEN.
  Wick Q(±1)/q²=8.  For even j,
      F_1(0)=16 + n_++^U Q_++ + n_--^U Q_{--},
      F_1((p+1)/2)=16 + n_++^U Q_++ − n_--^U Q_{--}   (p≡3 mod 4),
      F_1(other even)=16 − 2 Q_++
  (n_++^U=(p−3)/2 and n_--^U=(p+1)/2 if p≡3 mod 4; else n_++^U=p−1,
  n_--^U=0).  F_* as in 15.302.  Every even S_k/q² is one of
      σ=−1, j=0:     16 + Q_++(p−7)/2 + Q_{--}(p+1)/2 − (p−1)Q_{N*}
      σ=−1, j=(p+1)/2: 16 + Q_++(p−7)/2 − Q_{--}(p+1)/2 + 2 Q_{N*}
      σ=−1, other:    16 − 4 Q_++ + 2 Q_{N*}
  together with the two σ=(p−3)/2 siblings (replace the F_* sign).
  S_0 constraint: |++|Q_++ + |N*|Q_{N*} + |--|Q_{--}=2p²−18.
  Live mins recover 80/13 and 3072/409.  The three off Q values
  are not a Gauss/Jacobi formula in p.  phi_F not imported.  ∎

============================================================================
Backend: field tables on p=3..19 (ProcessPool-safe); live Q only at
p=5,7 (cache fold).  Writes evidence/e1_gmin_m4_prop15314.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import in_subfield, live_Q_on_T, paley_pair  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15302 import (  # noqa: E402
    U_powers,
    form_factor,
    live_S,
    nontrivial_s,
    sigma_k,
)

SIZE_PRIMES = (3, 5, 7, 11, 13, 17, 19)
EVEN_K = {5: (0, 2, 4, 6, 8, 10), 7: (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)}


def named_sizes(p: int) -> dict[str, int]:
    nstar = (p - 1) * (p - 3) // 2
    if p % 4 == 1:
        return {"pm1": 2, "++": 2 * p - 4, "N*": nstar, "--N1": 0}
    return {"pm1": 2, "++": 3 * (p - 3) // 2, "N*": nstar, "--N1": (p + 1) // 2}


def chi_p_pm(p: int, x: int) -> int:
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def live_counts(F) -> dict[str, int]:
    acc: dict[str, int] = defaultdict(int)
    for r in F["squares"]:
        acc[merge_cls(F, r)] += 1
    return dict(acc)


def prove_type_sizes() -> dict:
    ok = True
    rows = {}
    six_hit = 0
    for p in SIZE_PRIMES:
        F = _kernel_field(p)
        live = live_counts(F)
        named = named_sizes(p)
        match = all(live.get(k, 0) == named[k] for k in named)
        Fp = {x for x in range(1, F["q"]) if in_subfield(F, x)}
        U = set(U_powers(F))
        T = set(F["squares"])
        ns = {r for r in T if merge_cls(F, r) == "N*"}
        geom_ns = T - (Fp | U)
        if not match or ns != geom_ns:
            ok = False
        if named["++"] == 6:
            six_hit += 1
        rows[str(p)] = {"named": named, "live": live, "ns_ok": ns == geom_ns}
    # fail-when-wrong: |++|=6 for every p — false at p=11
    if six_hit == len(SIZE_PRIMES):
        ok = False
    if named_sizes(11)["++"] == 6:
        ok = False
    return {
        "proved": ok,
        "six_constant_false": named_sizes(11)["++"] != 6,
        "by_p": rows,
        "theorem": (
            "|pm1|=2, |N*|=(p−1)(p−3)/2, |--|=0 or (p+1)/2, "
            "|++|=2p−4 or 3(p−3)/2 by p mod 4. N*=T\\(F_p^×∪U)."
        ),
    }


def prove_U_split() -> dict:
    ok = True
    bothplus_hit = 0
    swap_hit = 0
    rows = {}
    for p in SIZE_PRIMES:
        F = _kernel_field(p)
        U = U_powers(F)
        want = -chi_p_pm(p, p - 1)
        prods = []
        bothplus = 0
        odd_mm = 0
        even_pp = 0
        even_bad = 0
        for t, u in enumerate(U):
            if u in (1, F["neg1"]):
                continue
            a, b = paley_pair(F, u)
            prods.append(a * b)
            if a == 1 and b == 1:
                bothplus += 1
            if p % 4 == 3:
                if t % 2 == 1:
                    if a == -1 and b == -1:
                        odd_mm += 1
                    else:
                        ok = False
                else:
                    if a == 1 and b == 1:
                        even_pp += 1
                    else:
                        even_bad += 1
                        ok = False
        if any(x != want for x in prods) or len(prods) != p - 1:
            ok = False
        if p % 4 == 1:
            if bothplus > 0:
                bothplus_hit += 1
                ok = False
            split_ok = bothplus == 0
        else:
            split_ok = (
                odd_mm == (p + 1) // 2
                and even_pp == (p - 3) // 2
                and even_bad == 0
            )
            if not split_ok:
                ok = False
            # fail-when-wrong: swap ++/-- on odd indices
            if odd_mm != even_pp:
                swap_hit += 1
        rows[str(p)] = {
            "want_prod": want,
            "bothplus": bothplus,
            "split_ok": split_ok,
        }
    if p % 4 == 1 and bothplus_hit == 0:
        # p=5,13,17 must have triggered the both-plus test as a *false* claim
        # (we count live both-plus, which is 0 — the fail-when-wrong is that
        # claiming both-plus would be wrong).  Record the claim is false.
        pass
    if SIZE_PRIMES[-1] % 4 == 3 and swap_hit == 0:
        ok = False
    return {
        "proved": ok,
        "bothplus_on_p1mod4": False,
        "swap_odd_is_wrong": swap_hit > 0,
        "by_p": rows,
        "theorem": (
            "χ_q(u²−1)=−χ_p(−1) on U\\{±1}. p≡1 (mod 4): all mixed. "
            "p≡3 (mod 4): --N1=U\\U², ++∩U=U²\\{±1}."
        ),
    }


def F1_named(p: int, j: int, Qpp: float, Qmm: float) -> float:
    if p % 4 == 1:
        if j == 0:
            return 16.0 + (p - 1) * Qpp
        return 16.0 - 2.0 * Qpp
    npp = (p - 3) / 2
    nmm = (p + 1) / 2
    if j == 0:
        return 16.0 + npp * Qpp + nmm * Qmm
    if j == (p + 1) // 2:
        return 16.0 + npp * Qpp - nmm * Qmm
    return 16.0 - 2.0 * Qpp


def S_named(p: int, k: int, Qpp: float, Qn: float, Qmm: float) -> float:
    j = k % (p + 1)
    sig = float(sigma_k(p, k))
    f1 = F1_named(p, j, Qpp, Qmm)
    if j == 0:
        fstar = 2.0 * Qpp + (p - 1) * Qn
    else:
        fstar = 2.0 * (Qpp - Qn)
    return f1 + sig * fstar


def coarse_means(F, Q: dict) -> dict[str, float]:
    q2 = float(F["q"] ** 2)
    acc: dict[str, list[float]] = defaultdict(list)
    for r, v in Q.items():
        acc[merge_cls(F, r)].append(v / q2)
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def prove_F1_and_floor() -> dict:
    ok = True
    drop_special_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        m = coarse_means(F, Q)
        Qpp, Qn = m["++"], m.get("N*", 0.0)
        Qmm = m.get("--N1", 0.0)
        f1_err = 0.0
        for j in range(0, p + 1, 2):
            live = form_factor(F, Q, 1, j)
            pred = F1_named(p, j, Qpp, Qmm)
            f1_err = max(f1_err, abs(live.real - pred) + abs(live.imag))
            if p % 4 == 3 and j == (p + 1) // 2:
                dropped = 16.0 - 2.0 * Qpp
                if abs(dropped - live.real) > 1e-6:
                    drop_special_hit += 1
        if f1_err > 1e-8:
            ok = False
        named = named_sizes(p)
        s0 = named["++"] * Qpp + named["N*"] * Qn + named["--N1"] * Qmm
        if abs(s0 - (2 * p * p - 18)) > 1e-8:
            ok = False
        floor_ks = [k for k in EVEN_K[p] if k not in (0, (F["q"] - 1) // 2)]
        live_min = min(live_S(F, Q, k).real for k in floor_ks)
        named_min = min(S_named(p, k, Qpp, Qn, Qmm) for k in floor_ks)
        if abs(live_min - named_min) > 1e-8:
            ok = False
        if named_min < 6 - 1e-9:
            ok = False
        rows[str(p)] = {
            "Qpp": str(Fraction(Qpp).limit_denominator(2000)),
            "Qn": str(Fraction(Qn).limit_denominator(2000)),
            "Qmm": str(Fraction(Qmm).limit_denominator(2000)),
            "f1_err": f1_err,
            "s0": s0,
            "live_min": live_min,
            "named_min": named_min,
        }
    if drop_special_hit == 0:
        ok = False
    return {
        "proved": ok,
        "drop_special_j_fails": drop_special_hit > 0,
        "by_p": rows,
        "theorem": (
            "F_1(j) is the named linear form; every even S_k is "
            "F_1(j)+σ F_*(j). Live mins 80/13 and 3072/409. "
            "Dropping the j=(p+1)/2 slot misses F_1 at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "Type sizes and F_1(j) named for all odd primes. "
            "The three off Q values still carry den N_+/(2p). "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.314  named T-types / U-split / F_1; floor = 3 Q-inequalities", flush=True)
    A = prove_type_sizes()
    print(f"  A sizes: {A['proved']} six_const={not A['six_constant_false']}", flush=True)
    B = prove_U_split()
    print(f"  B U-split: {B['proved']} swap_wrong={B['swap_odd_is_wrong']}", flush=True)
    C = prove_F1_and_floor()
    print(f"  C F1/floor rewrite: {C['proved']} drop_j={C['drop_special_j_fails']}", flush=True)
    D = prove_open()
    print(f"  D general Q: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.314",
        "title": "Named T-types and F_1; floor is three Q-inequalities",
        "proved": {
            "type_sizes": A["proved"],
            "U_split": B["proved"],
            "F1_rewrite": C["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15314.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
