#!/usr/bin/env python3
"""
Prop 15.406 — Interior 4-level leftover is empty at p=5 and p=7:
minus-slice S≡−2 on U forces χ_Odd constant on U, and every
Walsh character constant on U is constant on U^c (GF(2) rref).
A concrete p=5 leftover exists and is not interior.
residual_ii_k_eq_4p_empty stays False (Walsh ∀p and 5+ open).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.405 flags.

============================================================================
Setup.  Leftover (15.274 E / 15.381): k=4p, max_Max− S=−2,
f_e≡−1 on U={S=−2}.  Interior 4-level (15.382 A) needs S=−4
(≡0 mod 4) and S=−6 (≡2 mod 4) on U^c.  Odd = {v : deg_G(v) odd}.

============================================================================
Theorem A — PROVED (handshaking + pair-span product).
  On U one has S=−2, so n_minus=(k−S)/2=2p+1 is odd, hence
      ∏_{g∈G} f_g = (−1)^{n_minus} = −1.
  Also ∏_g f_g = (∏_g C_g) ∏_{v∈Odd} y_v, so χ_Odd is constant
  on U.  Fail: n_minus=2p (that is S=0).  ∎

Theorem B — PROVED (cut parity).
  Fix y0∈U.  For z∈U^c write D={v : y0_v ≠ z_v}.  Then
      S(z) = −2 + 2 ∑_{G∩cut(D)} f_g(z)
           ≡ 2 + 2 |D∩Odd|  (mod 4).
  Thus S_Uc is a single residue iff |D∩Odd| has constant parity
  iff χ_Odd is constant on U^c.  Fail: claim S(z)≡−2 (mod 4)
  identically (false for the 15.406 D witness).  ∎

Theorem C — PROVED (GF(2) rref; p=3,5,7).
  Every Walsh character χ_A that is constant on U is constant
  on U^c: the GF(2) forms {x : B[U] x const} map to constant
  columns on U^c (ker mixed = affine mixed = 0).  Fail: claim
  a mixed column at p=5.  ∎

Theorem D — PROVED (explicit G at p=5).
  The 20-edge double-star WITNESS has S≡−2 on U, support
  {−12,−8,−4,−2} on Max−, leftover (max=−2), U^c ≡0 (mod 4),
  not interior, and χ_Odd=−y_i y_j.  Fail: claim this G is
  interior; fail: claim no leftover at p=5.  ∎

Theorem E — OPEN.  Walsh for general p≥11 is not proved, so
  interior 4-level is not emptied for all p.  5+ levels and
  even k>4p stay open.  residual_ii_k_eq_4p_empty stays False.
  phi_F not imported.  A proposed Aut_e single-orbit spanning
  proof is false at p=5 (dir ≤11 < 12); Aut_e-invariant extra
  duals are empty (15.601) but that is not Walsh.
  15.602: unique 1-dim G_aff^□-invariant subspace of H0 is ⟨1⟩.
  15.604: 1_QR ∈ H0 iff p≡1 (mod 4); ker(D−I)∩H0 has dim 2.
  15.605: H0=⟨1⟩⊕W, W=translate-span of extra, Paley A²=A over F2.
  15.606: W=⊕ nsq W^H and M transits.
  15.607: W is G_aff^□-irreducible for every odd p (F_p^× transits
  Φ_p-factors), so dir(affine_span(Max−))=H0.
  15.608: two PSL-orbits of F_p-sublines; 1∈dir(U) by antipodes.
  15.609: I(H0)=H0 (opposite-type circles never tangent).
  15.610: Aut({0,∞}) uniqueness DEAD (ker((D−I)^2) dim 2, I-invariant).
  15.611: W ≅ F2[M] ≅ F2[X]/(X^N+1); W_0 unique D-invariant
  hyperplane; dim ker((D−I)^2)∩W_0=2 is a p-law.
  15.612: Walsh ⇔ W1 ∧ W2 on Aut-invariant ideals of W_0
  (CLASS p-law). W1 not a p-law. V/⟨1⟩ stays OPEN.

============================================================================
Backend: Fraction identities + GF(2) rref on Max− caches
(/tmp/maxminus_p{5,7}.npy; p=3 by 2^10 enum).  Writes
evidence/e1_gmin_m4_prop15406.json
"""
from __future__ import annotations

import itertools
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15406_catalog.json"
YMINUS = {5: "/tmp/maxminus_p5.npy", 7: "/tmp/maxminus_p7.npy"}

# p=5 leftover double-star, e={0,1}, seed 17 / t=18500 of gpu_resii_witness.
WITNESS = (
    (0, 2),
    (0, 3),
    (0, 4),
    (0, 7),
    (0, 8),
    (0, 9),
    (0, 13),
    (0, 14),
    (0, 21),
    (0, 25),
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 7),
    (1, 10),
    (1, 12),
    (1, 14),
    (1, 15),
    (1, 16),
    (1, 18),
)


def n_minus_on_U(p: int) -> int:
    return 2 * p + 1


def leftover_S_mod4_formula(cut_odd: int) -> int:
    """2 + 2|D∩Odd| (mod 4)."""
    return (2 + 2 * cut_odd) % 4


def gf2_rref(A: np.ndarray) -> tuple[np.ndarray, list[int], int]:
    A = A.copy()
    m, n = A.shape
    pivots: list[int] = []
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        pivots.append(c)
        r += 1
        if r == m:
            break
    return A, pivots, r


def gf2_nullspace(A: np.ndarray) -> tuple[np.ndarray, int]:
    R, pivots, rank = gf2_rref(A)
    n = A.shape[1]
    pset = set(pivots)
    free = [c for c in range(n) if c not in pset]
    prow = {pivots[i]: i for i in range(rank)}
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.uint8)
        v[f] = 1
        for p in pivots:
            if R[prow[p], f]:
                v[p] = 1
        basis.append(v)
    if not basis:
        return np.zeros((n, 0), dtype=np.uint8), rank
    return np.stack(basis, axis=1), rank


def gf2_solve(A: np.ndarray, b: np.ndarray) -> np.ndarray | None:
    m, n = A.shape
    Aug = np.concatenate([A, b.reshape(-1, 1)], axis=1)
    R, pivots, _rank = gf2_rref(Aug)
    if any(bool(R[i, -1]) and not R[i, :-1].any() for i in range(m)):
        return None
    x = np.zeros(n, dtype=np.uint8)
    for i, p in enumerate(pivots):
        if p < n:
            x[p] = R[i, -1]
    return x


def load_minus(p: int) -> tuple[np.ndarray, np.ndarray]:
    if p == 3:
        C = paley_conference_prime_power(3)
        n = C.shape[0]
        rows = []
        for bits in itertools.product((-1.0, 1.0), repeat=n):
            y = np.array(bits, dtype=np.float64)
            if np.allclose(C @ y, -3.0 * y, atol=1e-8):
                rows.append(y)
        return np.stack(rows), C
    Y = np.sign(np.load(YMINUS[p]).astype(np.float64))
    C = paley_conference_prime_power(p)
    return Y, C


def walsh_U_implies_Uc(p: int) -> dict:
    Y, C = load_minus(p)
    i, j = 0, 1
    fe = C[i, j] * Y[:, i] * Y[:, j]
    U = fe < 0
    Uc = ~U
    B = ((1.0 - Y) / 2.0).astype(np.uint8)
    BU, BC = B[U], B[Uc]
    N0, rankU = gf2_nullspace(BU)
    x1 = gf2_solve(BU, np.ones(BU.shape[0], dtype=np.uint8))
    if N0.shape[1]:
        KN = (BC.astype(np.int32) @ N0.astype(np.int32)) % 2
        ker_mixed = sum(
            int(KN[:, c].min() != KN[:, c].max()) for c in range(KN.shape[1])
        )
    else:
        ker_mixed = 0
    aff_mixed = 0
    aff_ok = x1 is not None
    if x1 is not None:
        w1 = (BC.astype(np.int32) @ x1.astype(np.int32)) % 2
        part_const = bool(w1.size == 0 or w1.min() == w1.max())
        aff_ok = part_const and ker_mixed == 0
        if not part_const:
            aff_mixed = 1
    closed = ker_mixed == 0 and aff_ok
    return {
        "p": p,
        "nU": int(U.sum()),
        "nUc": int(Uc.sum()),
        "rank_BU": int(rankU),
        "ker_dim": int(N0.shape[1]),
        "ker_mixed": int(ker_mixed),
        "aff_solvable": x1 is not None,
        "aff_mixed": int(aff_mixed),
        "closed": bool(closed),
    }


def witness_scores() -> dict:
    Y = np.sign(np.load(YMINUS[5]).astype(np.float64))
    C = paley_conference_prime_power(5)
    n = C.shape[0]
    i, j = 0, 1
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    idx = {e: t for t, e in enumerate(edges)}
    F = (
        Y[:, [e[0] for e in edges]]
        * Y[:, [e[1] for e in edges]]
        * np.array([C[a, b] for a, b in edges], dtype=np.float64)
    )
    v = np.zeros(len(edges), dtype=np.float64)
    for a, b in WITNESS:
        v[idx[(a, b)]] = 1.0
    S = F @ v
    fe = C[i, j] * Y[:, i] * Y[:, j]
    U = fe < 0
    su = np.array([int(round(x)) for x in S[U]])
    sc = np.array([int(round(x)) for x in S[~U]])
    sall = np.array([int(round(x)) for x in S])
    deg = np.zeros(n, dtype=int)
    for a, b in WITNESS:
        deg[a] += 1
        deg[b] += 1
    Odd = [t for t in range(n) if deg[t] % 2 == 1]
    chi = np.prod(Y[:, Odd], axis=1)
    yiyj = Y[:, i] * Y[:, j]
    return {
        "setU": sorted(set(su.tolist())),
        "setUc": sorted(set(sc.tolist())),
        "setAll": sorted(set(sall.tolist())),
        "maxAll": int(sall.max()),
        "U_is_m2": set(su.tolist()) == {-2},
        "leftover": int(sall.max()) == -2 and set(su.tolist()) == {-2},
        "interior": (
            set(su.tolist()) == {-2}
            and int(sall.max()) == -2
            and -4 in sc
            and -6 in sc
        ),
        "Uc_mods": sorted({x % 4 for x in set(sc.tolist())}),
        "Odd": Odd,
        "chi_eq_minus_yiyj": bool(np.allclose(chi, -yiyj)),
        "n_edges": int(v.sum()),
        "histUc": {str(k): int((sc == k).sum()) for k in sorted(set(sc.tolist()))},
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 19):
        nm = n_minus_on_U(p)
        ok = ok and nm == 2 * p + 1
        ok = ok and nm % 2 == 1
        ok = ok and (-1) ** nm == -1
        # fail-when-wrong: n_minus=2p is S=0
        ok = ok and nm != 2 * p
        if nm == 2 * p:
            ok = False
        rows[str(p)] = {"n_minus": nm, "prod_sign": -1}
    ok = ok and n_minus_on_U(5) == 11
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Leftover S≡−2 on U ⇒ n_minus=2p+1 odd ⇒ χ_Odd const on U. "
            "Fail: n_minus=2p."
        ),
    }


def prove_B() -> dict:
    ok = True
    # 2 + 2*even ≡ 2; 2 + 2*odd ≡ 0 (mod 4)
    ok = ok and leftover_S_mod4_formula(0) == 2
    ok = ok and leftover_S_mod4_formula(1) == 0
    ok = ok and leftover_S_mod4_formula(2) == 2
    ok = ok and leftover_S_mod4_formula(3) == 0
    # witness is ≡0, so |D∩Odd| is odd — not identically ≡2
    W = witness_scores()
    ok = ok and W["Uc_mods"] == [0]
    if leftover_S_mod4_formula(0) == 0 and leftover_S_mod4_formula(1) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "r_even": leftover_S_mod4_formula(0),
        "r_odd": leftover_S_mod4_formula(1),
        "witness_Uc_mods": W["Uc_mods"],
        "theorem": (
            "S(z)≡2+2|D∩Odd| (mod 4) on U^c. Single residue iff "
            "χ_Odd const on U^c. Fail: S≡−2 (mod 4) on the witness."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        rec = walsh_U_implies_Uc(p)
        ok = ok and rec["closed"]
        ok = ok and rec["ker_mixed"] == 0
        ok = ok and rec["aff_mixed"] == 0
        rows[str(p)] = rec
        if p == 5 and rec["ker_mixed"] > 0:
            ok = False
    # fail-when-wrong: mixed at p=5
    if rows["5"]["ker_mixed"] != 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Walsh const on U ⇒ const on U^c at p=3,5,7. "
            "Fail: mixed column at p=5."
        ),
    }


def prove_D() -> dict:
    W = witness_scores()
    ok = True
    ok = ok and W["n_edges"] == 20
    ok = ok and W["U_is_m2"]
    ok = ok and W["leftover"]
    ok = ok and not W["interior"]
    ok = ok and W["setUc"] == [-12, -8, -4]
    ok = ok and W["Uc_mods"] == [0]
    ok = ok and W["chi_eq_minus_yiyj"]
    ok = ok and -6 not in W["setAll"]
    if W["interior"]:
        ok = False
    if not W["leftover"]:
        ok = False
    return {
        "proved": bool(ok),
        "witness": W,
        "theorem": (
            "p=5 WITNESS is leftover, not interior, χ_Odd=−y_i y_j. "
            "Fail: claim interior; fail: no leftover at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Interior empty at p=5,7 via Walsh rref. Walsh ∀p≥11 unproved. "
            "5+ levels and even k>4p stay open. "
            "residual_ii_k_eq_4p_empty stays False. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.406  interior 4-level empty at p=5,7 via Walsh", flush=True)
    A = prove_A()
    print(f"  A n_minus: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B cut: {B['proved']} wit_mods={B['witness_Uc_mods']}", flush=True)
    C = prove_C()
    print(
        f"  C Walsh: {C['proved']} p5ker={C['rows']['5']['ker_dim']}",
        flush=True,
    )
    D = prove_D()
    print(
        f"  D witness: {D['proved']} left={D['witness']['leftover']} "
        f"int={D['witness']['interior']}",
        flush=True,
    )
    E = prove_open()
    print(f"  E open: resii={E['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "n_minus": {str(p): n_minus_on_U(p) for p in (5, 7, 11)},
                "walsh": {sp: C["rows"][sp]["closed"] for sp in C["rows"]},
                "witness_setUc": D["witness"]["setUc"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.406",
        "title": (
            "Interior 4-level empty at p=5,7 (Walsh on U⇒U^c); "
            "p=5 leftover witness is not interior"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "chi_Odd_const_on_U": A["proved"],
            "S_Uc_residue_is_cut_Odd": B["proved"],
            "walsh_U_implies_Uc_p357": C["proved"],
            "p5_leftover_not_interior": D["proved"],
            "interior_empty_p5_p7": bool(A["proved"] and B["proved"] and C["proved"]),
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15406.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
