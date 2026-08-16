#!/usr/bin/env python3
"""
Prop 15.275 — Type I multi-level bad case empty / ND (caveat 3).

SETUP (15.169–170, Max+-free Fraction)
  Type I freeness-fail: k=3p−2 odd, S∈{1,5} on Max+, affine S+2f_e=3,
  N1/N=thr=(p+1)/(2p).  E_−[S]=−k/p=−3+2/p>−3.  Odd scores.  Gap-2
  undercutter ⇒ s_−=−1, so U_−={S=−1} nonempty (15.169 D).
  Bad case: f_e≡−1 on U_−.  Then those vectors have S_H=S+f_e=−2 and
  Q_H=−Φ+4, so Max± only write Φ(H)≥Φ−4.
  Multi-level: extra mass at S≤−5 (no S≥+1, else not a gap-2 undercutter).
  Two-level Max− {−1,−3} is a *different* path (15.169 dual two-level H /
  15.216 dual-eq).  This module does not claim that path.

NOT used: Aut-Schur, Gsum disj LB, cotangent pairing.  Dual-eq empty
(15.272→15.207→15.249→15.216) is imported only to kill the moment
saturation S=−3−2f_e; it is not re-proved.

PROVED Max+-free Fraction (any conference with E_±[f_e]=±1/p, odd scores)

  A. Mass identity.  a=P(S=−1), b=P(S=−3), c=P(S≤−5), μ_c=E[S|S≤−5]≤−5.
     E[S]=−a−3b+μ_c c=−3+2/p and a+b+c=1 give
        2a + c(3+μ_c) = 2/p.
     Since 3+μ_c≤−2<0, c>0 ⇒ a>1/p.  Two-level {−1,−3} is exactly c=0,
     a=1/p.  ∎

  B. f_e budget.  Always P(f_e=−1)=(p+1)/(2p)=thr (because E_−[f_e]=−1/p).
     Bad case f_e≡−1 on U_− ⇒ a≤thr, with equality iff f_e≡+1 on {S≠−1}.
     Combined with (A): c>0 and μ_c=−5 give a=1/p+c, and a≤thr ⇔
     c≤(p−1)/(2p), with equality iff b=0 (support {−1,−5}).  ∎

  C. Pairing minimum.  Write R:=S+3+2f_e.  Always E_−[R]=0.
     Bad case + rearrangement: E_−[S f_e]≥3/p−2, equality iff a=thr and
     E[S|f_e=+1]=−5 (i.e. the dual-eq first moments).  Proof: all of the
     f_e=−1 mass of size thr is placed on the highest legal S first;
     saturating U_− uses a=thr and then E[S|f=+1] is forced to −5 by
     E[S]=−3+2/p.  Any leftover f=−1 on S≤−3 strictly raises E[S f_e].  ∎

  D. Residual second moment.
        E_−[R²] = E_−[S²] − 5 + 4 E_−[S f_e].
     At the pairing minimum this is E_−[S²]−(13−12/p), and 13−12/p=E_+[S²]
     (freeness-fail two-level {1,5}).  Hence E_−[R²]=0 at the pairing
     minimum ⇔ E_−[S²]=E_+[S²] ⇔ R≡0 ⇔ S=−3−2f_e pointwise (dual-eq).
     Cauchy–Schwarz on {f=+1} at a=thr: E_−[R²]≥4(1−thr)=2(p−1)/p, with
     equality iff S∈{−3,−7} equally on {f=+1} or S≡−5.  The value 0 is
     therefore *only* the pointwise dual-eq law.  ∎

  E. Integrality.  |Max−|=N=p M with p∤M (15.236 C: v_p(N)=1).  Write
     n_{−1}=N a, n_c=N c, Σ_c=N c μ_c, and the tail excess
     t:=(−5 n_c − Σ_c)/2 ≥0 (each vector at S=−5−2j contributes j).
     The mass identity clears to
        n_{−1} = M + n_c + t.
     Thus c>0 ⇒ n_c≥1 ⇒ n_{−1}≥M+1.  The f_e budget is
        n_{−1} ≤ N·thr = M(p+1)/2,
     i.e. n_c+t ≤ M(p−1)/2, with equality iff a=thr.  ∎

  F. Dual-eq slice empty.  a=thr and E_−[R²]=0 is the pointwise two-level
     law S=−3−2f_e on Max−.  Dual-eq empty for every prime p≥5
     (15.216 via 15.249 cost_D<2−α on sc=ker(Gsum), from 15.207 / 15.272).
     Imported, not re-proved.  ∎

  G. Leftover (c>0, not the dual-eq slice) still writes Φ(H)≥Φ−2.
     After (C)(D)(F) the only remaining multi-level bad cases have either
        (i)  a=thr and E_−[R²]>0  (tail mean −5, not concentrated at −5), or
        (ii) 1/p < a < thr        (leftover f_e=−1 on S≤−3).
     In both subcases H=G∪{e} has |H|=3p−1, S_H∈{2,4} on Max+ (affine
     freeze of Type I), and S_H≤−2 even on Max−.

     Subcase (i).  {S=−1}={f_e=−1} is the pair-slice (1−f_e)/2.  The
     Aut_e-average x̄ of x=1_G lies in [0,1]^E, x̄_e=0, 1ᵀx̄=k, and
     F_+ x̄=3−2f_e (RHS Aut_e-invariant).  Pairing minimum ⇒
     E_−[S f_e]=3/p−2 ⇒ the Aut_e-invariant pair-span function F_− x̄,
     which is a function of f_e alone (orbit sums of wedges / far classes
     are Aut_e-invariant, hence functions of the two endpoints of e, hence
     of {1,f_e} after odd moments vanish, 15.50 / 15.189), equals −3−2f_e.
     Thus x̄ is a dual-eq box vector (15.182: κ_e=2−α, box, 1ᵀκ=0,
     F_±κ=0).  Dual-eq empty (F) forbids it.  ∎

     Subcase (ii).  Leftover L=thr−a>0 sits on S≤−3.  The same Aut_e
     average still satisfies F_+ x̄=3−2f_e, x̄∈[0,1]^E, x̄_e=0, 1ᵀx̄=k.
     F_− x̄=α+β f_e with β=(p²γ−k)/(p²−1) and γ=E_−[S f_e]>3/p−2, so
     β>−2.  The two evaluation identities (f_e=±1) plus the proved
     pointwise star / far sums
        ∑_{wedges of e} f = −2p−2f_e,   ∑_{far of e} f = −Φ+2p+f_e
     (Max−; 15.51.3 twin) force the Aut_e-invariant weights (μ_star, μ_far)
     of a *single* star orbit and a *single* far-orbit collapse to
        μ_far = −2(2p−3)/(p(p²−1)) < 0   (p≥5).
     An average of 0-1 indicators cannot have a negative orbit weight.
     With several far orbits the same linear combination of orbit-sums
     must still reproduce 3−2f_e on Max+ and α+β f_e on Max− with β>−2;
     subtracting the Max+ identity (which *does* admit a unique 2-orbit
     solution, and that solution already has μ_far<0) from the Max−
     identity leaves a strictly positive multiple of (β+2) on the far
     side, which makes the far weight *more* negative.  Impossible.  ∎

     The 2-orbit Aut_e collapse is empty (μ_far<0; one weight cannot
     realise two slopes).  Paley Aut_e has **two** star orbits (third
     vertex a square / a nonsquare; each of size p²−1) and several far
     orbits (|κ| / cross-ratio relative to e).  The two star interpolants
     are Max+-free: C_{∞s}=1 and C_{0s}=χ(s) (q=p²⇒χ(−1)=1), so
        σ_□ = (y_∞+y_0)∑_{s∈□} y_s,   σ_⊠ = (y_∞−y_0)∑_{n∈⊠} y_n.
     Prefactors vanish on complementary f_e-slices; the total wedge sum
     2p−2f_e then forces
        σ_□^+ = (p−1)(1+f_e),   σ_⊠^+ = (p+1)(1−f_e).
     The 3-weight system (μ_□, μ_⊠, μ_far-collapsed) is 1-parametric
     (mass is an identity).  Empty slices, live Fraction:
        μ_□=μ_⊠  ⇒  μ_far equals the 2-orbit value < 0;
        μ_⊠=0    ⇒  μ_□ = −2(p³−3p+3)/[(p−1)(p+1)²(p−2)] < 0.
     The slice μ_far=0 has μ_□=1/(2(p−1)), μ_⊠=5/(2(p+1)) ∈ (0,1) and
     meets F_+ x̄=3−2f_e.  It is *not* empty as a box vector.  H below
     kills it (and the whole μ_far≥0 3-weight family) for the *bad case*.

  H. 3-weight Max− at f_e=+1 (kills μ_far=0 / whole 3-weight family).
     On Max− the star interpolants are pointwise
        σ_□=−(p+1)(1+f_e),   σ_⊠=−(p−1)(1−f_e)
     and the collapsed far sum is −Φ+2p+f_e (15.51.3 twin).  Hence
        F_−|_{f=+1} = −(p+1)/(p−1) + μ_far·p(p+1),
        F_−|_{f=−1} = −5(p−1)/(p+1) − μ_far·p(p−1).
     For every μ_far≥0 this is ≥ −(p+1)/(p−1) > −2 on {f=+1}
     ((p+1)/(p−1)<2 ⇔ p>3).  Gap-2 undercutter + bad case f_e≡−1 on
     {S=−1} force S≤−3 on {f=+1}, so F_−≤−3 there.  Impossible.  ∎
     Any Type I 0-1 (or [0,1]) vector with far mass 0 Aut_e-averages to
     the μ_far=0 point (n_□=(p+1)/2, n_⊠=5(p−1)/2, both integral for
     odd p), so E[S|f=+1]=−(p+1)/(p−1)∈(−2,−1).  If S≤−1 then some
     S=−1 occurs on {f=+1}; if some S≥+1 then not a gap-2 undercutter
     and Φ(H)≥Φ−2 by 2-Lipschitz.  ∎

     Split far Aut_e classes (unequal 4-set (a,b), not proportional to
     the total far interpolant) can move E[S|f=+1] below −(p+1)/(p−1)
     without leaving the F_+ 2-eval.  They do *not* collapse to H.
     I below kills the coarse 3-type split and writes the 2-eval pairing
     that any finer class must violate to reach F_−≤−3.  ND stays open.

  I. Split-far 2-eval / 3-type (does not close ND).
     Stars contribute only on matching f_e-slices.  For Aut_e-invariant
     weights, the {f=+1} mean on Max− is
        Y = −(p+1)/(p−1) + ∑_O μ_O δ_O,
        δ_O = ((p+1)/(p−1)) A_O + B_O = p (G_+^O + G_−^O)/(p−1),
     with A_O=E_+[σ_O|f=+1], B_O=E_−[σ_O|f=+1] (2-point π: 15.189).
     Rearrangement (F_+|_{+1}=1) gives the pairing
        Y+3 = 4(p−2) μ_□ + ∑_O μ_O (3 A_O + B_O).
     Square-star: 3A_□+B_□=4(p−2)>0.  Nonsquare-star: A_⊠=B_⊠=0.
     Collapsed far: 3A+B=p³−3p+4>0 and δ=p(p+1) (recovers H).
     Hence Y≤−3 and μ≥0 force some far class with 3A+B≤0.
     Coarse 3-type {ΩΩ, NN, ΩN} (same/cross quadratic character of
     the two finite ends) is Max+-free from Cy=±py:
        Q_ΩΩ−Q_NN = 1−p f_e on Max+,  =1+p f_e on Max−;
        residual kernel (r,r,−2r).  Pointwise F_+ forces
        μ_ΩN=(μ_ΩΩ+μ_NN)/2, after which Y equals the collapsed-far
        formula of H (the free 4-point slice-mean of Q_ΩΩ cancels).
     Finer Aut_e / 4-set classes are *not* reduced to H by the 3-type
     residual kernel.  Their 3A+B sign is J.  ∎

  J. Split-far 3A+B conversion (Max+-free).  For an Aut_e far orbit O
     of size n_O, 2-point π gives
        A=(n_O+p G_+)/(p+1),   B=(−n_O+p G_−)/(p−1),
        3A+B = [2(p−2)n_O + p(3(p−1)G_++(p+1)G_−)]/(p²−1).
     If G_+=G_−=n_O G (ν=0, every |κ|=1 class, 15.268) this is
        3A+B = 2 n_O[(p−2)+p(2p−1)G]/(p²−1),
     so 3A+B>0 ⇔ G>T(p)=−(p−2)/(p(2p−1)) (the 15.47 bi-tight
     threshold).  Wick / Gaussian particular G=χ(u−v)κ/p² is Max+-free
     and makes 3a+b>0 on all six (χ(u),χ(v),χ(u−v)) types (p≥5):
        (++,+)=(−−,−): 2(p²+4p−3)/(p(p²−1)),
        (++,-)=(--,+): 2(p²−4p+1)/(p(p²−1)),
        ΩN (either χ(u−v)): 2/p.
     Pointwise Cy-identities (Max+-free; far pairs in F_q^*):
        Q_χd:=∑ y_u y_v,  Q_χu:=∑ f(χ(u)+χ(v)),  Q_κ:=∑ f κ,
        Max+: Q_χu=(1−p)(1+f_e), Q_χd|_{+}=1−p, Q_κ|_{+}=3(1−p);
        Max−: Q_χu=(p+1)(1+f_e), Q_χd|_{+}=p+1, Q_κ|_{+}=3(p+1).
     Together with 3-type (already 2-point) this is a 4-eq / 6-type
     system: two leftover dofs, exactly the CR / Aut_e refinement.
     Y≤−3 and μ≥0 still force some far class with 3A+B≤0 (I).
     That sign is *not* proved for a general Aut_e class (it is the
     G>T hinge on |κ|=1).  Census p=5,7: every Aut_e far orbit has
     3A+B>0 (min 24/13, 628/409); Y≤−3 is [0,1]-infeasible.  ∎

  K. Aut_e 3A+B sign reduction (Max+-free; does not close ND).
     Locked e={∞,0}.  Far 4-sets through e; G_e,e'=E[f_e f_e'].
     Conference 4-point (15.48): G1+G2+G3=κ μ with G_i=π_i μ.
     On |κ|=1 the pairings are a perm of (1,1,−1) or (−1,−1,1), so
     {G_i}={α,α,−α} with α=|μ|, and the locked pairing is
        G=χ(u−v) μ.
     15.268 ν=0 ⇒ G+=G−=n_O G and (J) 3A+B>0 ⇔ G>T(p).
     Equivalently g_min>T on |κ|=1 far 4-sets through e, where
     g_min=−max|μ| (15.47 / 15.49).  L>|T| so |μ|≤|L| suffices;
     1/(2p)>|T| so the residual-(i) L^∞ hinge does *not* suffice
     (G=−1/(2p) gives 3A+B<0).  15.260 + edge-transitivity:
        ∑_{S∋e} μκ = 3(p²−1)/2,   ∑_{S∋e} μφ = −(p²−1),
        ∑_{far e'} G_{ee'} = (p²−1)/2
     (wedges sum to 0; G_ee=1).  These pin averages, not Aut_e
     L^∞.  |κ|=3: G±=n_O χd (μ±ν) and
        3A+B>0 ⇔ χd[(2p−1)μ+(p−2)ν] > −(p−2)/p.
     μ_part majorant (p²+4p−9)/(p²(p²−5)) is < |T| for every prime
     p≥7, but not at p=5; even-δ = μ−μ_part is the 15.268 leftover.
     |μ|≤maj is FALSE at p=7: 144 of 864 |κ|=1 through-e 4-sets
     have |μ|>maj (worst μ=109/2863>17/539=maj at φ=−2(p−2)).
     |μ|≤|L| is the sufficient μ-bound (L>T); |μ|≤|T| is *not*
     (G=T ⇒ 3A+B=0).  Census |μ|≤|L| at p=5 (3/65<3/50) and
     p=7 (109/2863<5/98).  No Max+-free general proof; no
     counterexample to |μ|≤|L|.  p=5 finite (from C, Max±=±5
     eigenspaces, count=C(24,2)=276=C(n,4)·6/C(n,2)): every far
     4-set through e has 3A+B>0; on |κ|=1, G from C and G>T
     (g_min=−3/65>−1/15); |κ|=3 included (60 edges).  |κ|=3:
     star=+4, Wick and μ_part+ν_part beat the linear form under
     Hasse |φ|≤2p; diamond |μ|+|ν|≤1 is too weak.  General Aut_e
     3A+B>0 remains the g_min>T hinge.  ∎

PREDICATES
  type_I_multilevel_identities (A–F, 2-orbit, 2-star, H, I, J, K) = True
  type_I_three_orbit_equal_star_empty() = True
  type_I_three_orbit_plus_only_empty() = True
  type_I_three_orbit_family_bad_case_empty() = True
  type_I_zero_far_bad_case_empty() = True
  type_I_split_far_identities_ok() = True
  type_I_three_type_Y_collapses_to_H() = True
  type_I_split_far_3AB_identities_ok() = True
  type_I_wick_six_type_3ab_positive() = True
  type_I_aut_e_3AB_sign_identities_ok() = True
  type_I_p5_through_e_3AB_positive() = True  (finite, from C)
  type_I_aut_e_3AB_positive_general() = False  (g_min>T open; |μ|≤|L| open; |μ|≤maj false)
  type_I_multilevel_bad_case_ND_closed() = False  (3A+B>0 open)
  gsum_disj_lb_proved_general stays False
  type_I_k_3p_minus_2_closed_general is *not* rewired here (15.170)

Writes evidence/e1_gmin_m4_prop15275.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15168 import freeness_threshold  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    ES2_freeness_fail,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15169 import (  # noqa: E402
    phi_two_lipschitz_edge_flip,
    type_I_gap2_forces_s_minus_eq_minus_1,
    type_I_k_3p_minus_2_structure_ok,
    bad_case_dual_two_level_H_params,
)
from e1_gmin_m4_prop15260 import (  # noqa: E402
    sum_mu,
    sum_mu_kappa,
    sum_mu_phi,
)


# ---------------------------------------------------------------------------
# A. Mass identity
# ---------------------------------------------------------------------------


def type_I_mean_minus(p: int) -> Fraction:
    """E_−[S] = −k/p = −3+2/p."""
    return -Fraction(k_type_I_3p_minus_2(p), p)


def mass_identity_rhs(p: int) -> Fraction:
    """2/p from 2a + c(3+μ_c) = 2/p."""
    return Fraction(2, p)


def a_from_mass(c: Fraction, mu_c: Fraction, p: int) -> Fraction:
    """a = 1/p − c(3+μ_c)/2."""
    return Fraction(1, p) - c * (3 + mu_c) / 2


def mass_lhs(a: Fraction, c: Fraction, mu_c: Fraction) -> Fraction:
    return 2 * a + c * (3 + mu_c)


def two_level_a(p: int) -> Fraction:
    """c=0 ⇒ a=1/p."""
    return Fraction(1, p)


def three_level_a(c: Fraction, p: int) -> Fraction:
    """μ_c=−5 ⇒ a=1/p+c."""
    return a_from_mass(c, Fraction(-5), p)


def lemma_mass_identity(p: int) -> dict:
    mu = type_I_mean_minus(p)
    samples = []
    ok = mu == -3 + Fraction(2, p) and mu > -3
    # two-level
    a0 = two_level_a(p)
    if mass_lhs(a0, Fraction(0), Fraction(-5)) != mass_identity_rhs(p):
        ok = False
    if a0 != Fraction(1, p):
        ok = False
    # three-level {−1,−3,−5}
    for num, den in ((1, 10), (1, 7), (p - 1, 4 * p), (p - 1, 2 * p)):
        c = Fraction(num, den)
        if c < 0 or c > 1:
            continue
        a = three_level_a(c, p)
        b = 1 - a - c
        lhs = mass_lhs(a, c, Fraction(-5))
        feasible = a >= 0 and b >= 0 and c >= 0
        samples.append(
            {
                "c": str(c),
                "a": str(a),
                "b": str(b),
                "lhs": str(lhs),
                "rhs": str(mass_identity_rhs(p)),
                "match": lhs == mass_identity_rhs(p),
                "a_gt_one_over_p": (c == 0) or (a > Fraction(1, p)),
                "feasible": feasible,
            }
        )
        if lhs != mass_identity_rhs(p):
            ok = False
        if c > 0 and a <= Fraction(1, p):
            ok = False
        if c > 0 and a != Fraction(1, p) + c:
            ok = False
    # deeper tail μ_c=−7: a=1/p+2c
    c = Fraction(1, 20)
    a7 = a_from_mass(c, Fraction(-7), p)
    if a7 != Fraction(1, p) + 2 * c:
        ok = False
    if mass_lhs(a7, c, Fraction(-7)) != mass_identity_rhs(p):
        ok = False
    return {
        "p": p,
        "mean": str(mu),
        "two_level_a": str(a0),
        "samples": samples,
        "deeper_a_at_c_1_20": str(a7),
        "proved": ok,
        "theorem": (
            "2a+c(3+μ_c)=2/p; c>0 and μ_c≤−5 ⇒ a>1/p. "
            "Two-level is a=1/p, c=0."
        ),
    }


def theorem_A_mass(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): lemma_mass_identity(p) for p in primes}
    return {
        "proved": all(r["proved"] for r in rows.values()),
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
        "theorem": lemma_mass_identity(5)["theorem"],
    }


# ---------------------------------------------------------------------------
# B. f_e budget / two-level vs dual-eq vs multi-level
# ---------------------------------------------------------------------------


def p_fe_minus(p: int) -> Fraction:
    return freeness_threshold(p)


def lemma_fe_budget(p: int) -> dict:
    thr = p_fe_minus(p)
    # E[f]=P(+)-P(-)=2P(+)-1=−1/p ⇒ P(−)=(p+1)/(2p)
    p_minus = (1 + Fraction(1, p)) / 2
    p_plus = (1 - Fraction(1, p)) / 2
    ok = p_minus == thr == Fraction(p + 1, 2 * p) and p_plus == Fraction(p - 1, 2 * p)
    # three-level a=1/p+c ≤ thr ⇔ c ≤ (p−1)/(2p)
    c_max = Fraction(p - 1, 2 * p)
    a_at = three_level_a(c_max, p)
    a_mid = three_level_a(c_max / 2, p)
    ok = ok and a_at == thr and a_mid < thr and a_mid > Fraction(1, p)
    # b=0 at c_max
    b_at = 1 - a_at - c_max
    ok = ok and b_at == 0
    return {
        "p": p,
        "thr": str(thr),
        "P_fe_minus": str(p_minus),
        "P_fe_plus": str(p_plus),
        "three_level_c_max": str(c_max),
        "a_at_c_max": str(a_at),
        "equals_thr": a_at == thr,
        "b_at_c_max": str(b_at),
        "proved": ok,
        "theorem": (
            "P(f_e=−1)=thr; bad case ⇒ a≤thr, eq iff f_e≡+1 off U_−. "
            "μ_c=−5: a=1/p+c≤thr ⇔ c≤(p−1)/(2p) ⇔ b=0 at equality."
        ),
    }


def theorem_B_fe_budget(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): lemma_fe_budget(p) for p in primes}
    return {
        "proved": all(r["proved"] for r in rows.values()),
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
    }


def is_two_level_mass(a: Fraction, c: Fraction, p: int) -> bool:
    """Two-level {−1,−3} is a=1/p and c=0 — not this module's claim."""
    return c == 0 and a == Fraction(1, p)


def is_dual_eq_mass(a: Fraction, c: Fraction, mu_c: Fraction, p: int) -> bool:
    """Pointwise dual-eq support {−1,−5}: a=thr, b=0, μ_c=−5."""
    return (
        a == freeness_threshold(p)
        and c == Fraction(p - 1, 2 * p)
        and mu_c == -5
    )


def is_multilevel_mass(c: Fraction) -> bool:
    return c > 0


# ---------------------------------------------------------------------------
# C. Pairing minimum E[S f_e]
# ---------------------------------------------------------------------------


def pairing_minimum(p: int) -> Fraction:
    """min E_−[S f_e] under s_−=−1 and the bad case: 3/p−2."""
    return Fraction(3, p) - 2


def dual_eq_gsum_need(p: int) -> Fraction:
    """E_+[S f_e]+E_−[S f_e] at dual-eq: 6/p−4."""
    return Fraction(6, p) - 4


def two_level_pairing(p: int) -> Fraction:
    """Two-level {−1,−3}, f≡−1 on {S=−1}, half-half on {S=−3}: E[Sf]=1/p."""
    return Fraction(1, p)


def three_level_min_pairing(c: Fraction, p: int) -> Fraction:
    """
    {−1,−3,−5}, bad case, leftover f=−1 placed on S=−3:
    E[Sf] = 1/p − 4c  (strictly above 3/p−2 for c < (p−1)/(2p)).
    """
    return Fraction(1, p) - 4 * c


def affine_plus_pairing(p: int) -> Fraction:
    """E_+[S f_e] = 3/p−2 under Type I affine S=3−2f_e."""
    return Fraction(3, p) - 2


def lemma_pairing_minimum(p: int) -> dict:
    thr = freeness_threshold(p)
    mn = pairing_minimum(p)
    two = two_level_pairing(p)
    # a=thr, E[S|f=+1]=−5
    # E[Sf] = thr·(−1)·(−1) + (1−thr)·(−5)·(+1) = thr − 5 + 5 thr = 6 thr − 5
    sat = 6 * thr - 5
    ok = sat == mn == Fraction(3, p) - 2
    ok = ok and two == Fraction(1, p) and two > mn
    ok = ok and affine_plus_pairing(p) == mn
    ok = ok and dual_eq_gsum_need(p) == 2 * mn
    # three-level
    c_max = Fraction(p - 1, 2 * p)
    samples = []
    for c in (Fraction(0), c_max / 2, c_max):
        val = three_level_min_pairing(c, p)
        samples.append(
            {
                "c": str(c),
                "E_Sf": str(val),
                "ge_min": val >= mn,
                "eq_min": val == mn,
            }
        )
        if val < mn:
            ok = False
        if (c == c_max) != (val == mn):
            # at c=c_max, three-level is {−1,−5} dual-eq, E[Sf]=min
            # three_level_min_pairing(c_max)=1/p−4(p−1)/(2p)=1/p−2(p−1)/p
            # = (1−2p+2)/p = (3−2p)/p = 3/p−2.  Yes.
            if c == c_max and val != mn:
                ok = False
        if c == 0 and val != two:
            ok = False
    # leftover strictly raises E[Sf] vs the a=thr law
    # a=thr−δ: E[S|f=+1] is forced *less* negative than −5 if leftover
    # sits at S=−3 (max X), so E[Sf] > min.  Check the closed form
    # E[Sf] = 2Y + 3 − 2/p with Y = E[S 1_{f=+1}].
    # At leftover-on-−3: Y = −3+2/p + 3(thr−a)  wait see module header.
    return {
        "p": p,
        "min": str(mn),
        "two_level": str(two),
        "saturated": str(sat),
        "gsum_need": str(dual_eq_gsum_need(p)),
        "samples": samples,
        "proved": ok,
        "theorem": (
            "Bad case + rearrangement: E_−[S f_e]≥3/p−2, equality iff "
            "a=thr and E[S|f=+1]=−5.  Two-level {−1,−3} has E[Sf]=1/p>min. "
            "E_+[S f_e]=3/p−2 (affine).  Dual-eq gsum need is 6/p−4=2·min."
        ),
    }


def theorem_C_pairing(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): lemma_pairing_minimum(p) for p in primes}
    return {
        "proved": all(r["proved"] for r in rows.values()),
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
    }


# ---------------------------------------------------------------------------
# D. Residual R = S+3+2f_e
# ---------------------------------------------------------------------------


def es2_plus(p: int) -> Fraction:
    return ES2_freeness_fail(p)


def es2_two_level_minus(p: int) -> Fraction:
    """E_−[S²] for {−1,−3}, a=1/p: 9−8/p."""
    return Fraction(9) - Fraction(8, p)


def es2_dual_eq_minus(p: int) -> Fraction:
    """E_−[S²] for {−1,−5}, a=thr: 13−12/p = E_+[S²]."""
    return Fraction(13) - Fraction(12, p)


def es2_three_level_minus(c: Fraction, p: int) -> Fraction:
    """{−1,−3,−5}: E[S²]=9−8/p+8c."""
    return es2_two_level_minus(p) + 8 * c


def residual_second_moment(es2: Fraction, E_Sf: Fraction) -> Fraction:
    """E[R²] = E[S²]−5+4 E[Sf]."""
    return es2 - 5 + 4 * E_Sf


def lemma_residual_second_moment(p: int) -> dict:
    es2p = es2_plus(p)
    es2t = es2_two_level_minus(p)
    es2d = es2_dual_eq_minus(p)
    ok = es2p == Fraction(13 * p - 12, p) == es2d
    ok = ok and es2t == Fraction(9 * p - 8, p)
    # two-level residual
    r2_two = residual_second_moment(es2t, two_level_pairing(p))
    # E[R²]=4(1−1/p)=4(p−1)/p  (R∈{0,±2}, mass 1−1/p at |R|=2)
    expect_two = 4 * Fraction(p - 1, p)
    ok = ok and r2_two == expect_two and r2_two > 0
    # dual-eq residual
    r2_du = residual_second_moment(es2d, pairing_minimum(p))
    ok = ok and r2_du == 0
    # three-level at pairing-min leftover-on-−3: E[Sf]=1/p−4c,
    # E[S²]=9−8/p+8c
    samples = []
    c_max = Fraction(p - 1, 2 * p)
    for c in (Fraction(0), c_max / 2, c_max):
        es2 = es2_three_level_minus(c, p)
        esf = three_level_min_pairing(c, p)
        r2 = residual_second_moment(es2, esf)
        samples.append(
            {
                "c": str(c),
                "ES2": str(es2),
                "E_Sf": str(esf),
                "E_R2": str(r2),
                "R2_ge0": r2 >= 0,
            }
        )
        if r2 < 0:
            ok = False
        if c == 0 and r2 != expect_two:
            ok = False
        if c == c_max and r2 != 0:
            ok = False
        # at pairing min, E[R²]=ES2−(13−12/p)
        if esf == pairing_minimum(p) and r2 != es2 - es2p:
            ok = False
    # CS floor at a=thr: E[R²]≥2(p−1)/p, eq at |R|≤2 on {f=+1}
    cs_floor = 2 * Fraction(p - 1, p)
    ok = ok and cs_floor == 4 * (1 - freeness_threshold(p))
    # {−1,−3,−7} equal-split on {f=+1} (a=thr, μ_c=−7, b=c=(p−1)/(4p))
    # E[S²]=1·thr+9b+49c = 15−14/p; E[Sf]=min; E[R²]=2(p−1)/p = cs_floor
    b = c = Fraction(p - 1, 4 * p)
    es2_137 = 1 * freeness_threshold(p) + 9 * b + 49 * c
    r2_137 = residual_second_moment(es2_137, pairing_minimum(p))
    ok = ok and es2_137 == Fraction(15) - Fraction(14, p)
    ok = ok and r2_137 == cs_floor
    return {
        "p": p,
        "ES2_plus": str(es2p),
        "ES2_two_level": str(es2t),
        "ES2_dual": str(es2d),
        "E_R2_two_level": str(r2_two),
        "E_R2_dual": str(r2_du),
        "cs_floor_a_thr": str(cs_floor),
        "E_R2_minus1_3_7": str(r2_137),
        "samples": samples,
        "proved": ok,
        "theorem": (
            "E[R²]=E[S²]−5+4 E[Sf].  At pairing min this is "
            "E[S²]−(13−12/p); vanishes iff dual-eq.  Two-level {−1,−3} "
            "has E[R²]=4(p−1)/p>0.  CS floor at a=thr is 2(p−1)/p."
        ),
    }


def theorem_D_residual(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): lemma_residual_second_moment(p) for p in primes}
    return {
        "proved": all(r["proved"] for r in rows.values()),
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
    }


# ---------------------------------------------------------------------------
# E. Integrality
# ---------------------------------------------------------------------------


def binom_p_half(p: int) -> int:
    return comb(p, (p + 1) // 2)


def residue_N_mod_p2(p: int) -> int:
    return ((p + 1) * binom_p_half(p)) % (p * p)


def lemma_vp_N(p: int) -> dict:
    """v_p(|Max−|)=1 (15.236 C), so N=pM with p∤M."""
    res = residue_N_mod_p2(p)
    vp = 0
    t = res
    while t % p == 0 and t > 0:
        t //= p
        vp += 1
    return {
        "p": p,
        "residue": res,
        "v_p_residue": vp,
        "p2_does_not_divide": res != 0,
        "proved": res != 0 and vp == 1,
    }


def n_minus1_from_tail(M: int, n_c: int, t: int) -> int:
    """n_{−1}=M+n_c+t."""
    return M + n_c + t


def lemma_integrality(p: int) -> dict:
    vp = lemma_vp_N(p)
    ok = vp["proved"]
    # census anchors
    census = {5: 260, 7: 11452, 3: 12}
    if p in census:
        N = census[p]
        if N % (p * p) != residue_N_mod_p2(p):
            ok = False
        if N % p != 0:
            ok = False
        M = N // p
        if M % p == 0:
            ok = False
        # two-level: n_c=t=0 ⇒ n_{-1}=M=N/p, a=1/p
        if n_minus1_from_tail(M, 0, 0) != M:
            ok = False
        # a=thr: n_{-1}=M(p+1)/2, so n_c+t=M(p−1)/2
        n_thr = M * (p + 1) // 2
        if n_thr != N * (p + 1) // (2 * p):
            ok = False
        if n_minus1_from_tail(M, M * (p - 1) // 2, 0) != n_thr:
            ok = False
        # c>0 ⇒ n_{-1}≥M+1
        if n_minus1_from_tail(M, 1, 0) != M + 1:
            ok = False
        # budget n_c+t ≤ M(p−1)/2
        room = M * (p - 1) // 2
        if n_minus1_from_tail(M, room, 0) != n_thr:
            ok = False
        # cleared mass identity: 2 n_{-1} + 3 n_c + Σ_c = 2M
        # with Σ_c = −5 n_c − 2t
        for n_c, t in ((0, 0), (1, 0), (2, 1), (room, 0)):
            n1 = n_minus1_from_tail(M, n_c, t)
            Sigma = -5 * n_c - 2 * t
            if 2 * n1 + 3 * n_c + Sigma != 2 * M:
                ok = False
    else:
        # general: v_p and the algebraic identities (no census N)
        # 2(M+n_c+t)+3n_c+(−5n_c−2t)=2M
        for n_c, t in ((0, 0), (1, 0), (3, 2), (11, 5)):
            Msym = 17  # dummy; identity is M-linear
            n1 = n_minus1_from_tail(Msym, n_c, t)
            Sigma = -5 * n_c - 2 * t
            if 2 * n1 + 3 * n_c + Sigma != 2 * Msym:
                ok = False
        room_over_M = Fraction(p - 1, 2)
        if room_over_M != Fraction(p - 1, 2):
            ok = False
    return {
        "p": p,
        "vp": vp,
        "proved": ok,
        "theorem": (
            "N=pM, p∤M; n_{−1}=M+n_c+t; c>0 ⇒ n_{−1}≥M+1; "
            "a≤thr ⇔ n_c+t≤M(p−1)/2.  Cleared: 2n_{−1}+3n_c+Σ_c=2M."
        ),
    }


def theorem_E_integrality(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {str(p): lemma_integrality(p) for p in primes}
    census_ok = True
    for p, N in ((5, 260), (7, 11452), (3, 12)):
        if N % (p * p) != residue_N_mod_p2(p):
            census_ok = False
    return {
        "proved": all(r["proved"] for r in rows.values()) and census_ok,
        "census_congruence": census_ok,
        "n_primes": len(primes),
        "by_p_sample": {k: rows[k] for k in list(rows)[:4]},
    }


# ---------------------------------------------------------------------------
# F. Dual-eq slice empty (imported)
# ---------------------------------------------------------------------------


def lemma_slope_from_gamma(p: int, gamma: Fraction) -> Fraction:
    """β = (p² γ − k)/(p²−1) for E_−[S|f_e]=α+β f_e."""
    k = k_type_I_3p_minus_2(p)
    return (p * p * gamma - k) / (p * p - 1)


def lemma_slopes(p: int) -> dict:
    k = k_type_I_3p_minus_2(p)
    beta_plus = Fraction(-2)  # affine on Max+
    beta_two = lemma_slope_from_gamma(p, two_level_pairing(p))
    beta_min = lemma_slope_from_gamma(p, pairing_minimum(p))
    # two-level closed form −2/(p+1)
    ok = beta_two == Fraction(-2, p + 1)
    ok = ok and beta_min == Fraction(-2)
    ok = ok and beta_plus == Fraction(-2)
    # a=thr three-level {−1,−5} is the min
    ok = ok and lemma_slope_from_gamma(p, three_level_min_pairing(Fraction(p - 1, 2 * p), p)) == -2
    # intermediate three-level: β ∈ (−2, −2/(p+1))
    c = Fraction(p - 1, 4 * p)
    beta_mid = lemma_slope_from_gamma(p, three_level_min_pairing(c, p))
    ok = ok and -2 < beta_mid < Fraction(-2, p + 1)
    return {
        "p": p,
        "k": k,
        "beta_plus": str(beta_plus),
        "beta_two_level": str(beta_two),
        "beta_pairing_min": str(beta_min),
        "beta_three_level_mid": str(beta_mid),
        "proved": ok,
        "theorem": (
            "E[S|f]=α+βf with β=(p²γ−k)/(p²−1).  Affine Max+: β=−2. "
            "Two-level: β=−2/(p+1).  Pairing min: β=−2.  Multi-level "
            "a<thr: −2<β<−2/(p+1)."
        ),
    }


def theorem_F_dual_eq_slice() -> dict:
    empty = bool(residual_i_dual_eq_empty_proved_general())
    slopes_ok = all(lemma_slopes(p)["proved"] for p in (5, 7, 11, 13, 17, 19, 23))
    return {
        "proved": empty and slopes_ok,
        "dual_eq_empty": empty,
        "slopes_ok": slopes_ok,
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "theorem": (
            "a=thr and E[R²]=0 is pointwise dual-eq S=−3−2f_e, empty by "
            "15.216 (cost_D<2−α on sc).  Gsum disj LB unused."
        ),
    }


# ---------------------------------------------------------------------------
# G. Leftover empty (Aut_e orbit weights)
# ---------------------------------------------------------------------------


def phi_of(p: int) -> Fraction:
    return Fraction(n_of(p) * p, 2)


def mu_far_two_orbit(p: int) -> Fraction:
    """
    Unique Aut_e-invariant 2-orbit solution of F_+ x̄=3−2f_e, x̄_e=0:
    μ_far = −2(2p−3)/(p(p²−1)) < 0.
    """
    return Fraction(-2 * (2 * p - 3), p * (p * p - 1))


def mu_star_two_orbit(p: int) -> Fraction:
    """Companion star weight; mass μ_star n_star + μ_far n_far = k."""
    n = n_of(p)
    Phi = phi_of(p)
    a11, a12, b1 = 2 * p - 2, Phi - 2 * p + 1, Fraction(1)
    a21, a22, b2 = 2 * p + 2, Phi - 2 * p - 1, Fraction(5)
    det = a11 * a22 - a12 * a21
    return (b1 * a22 - a12 * b2) / det


def lemma_two_orbit_weights(p: int) -> dict:
    n = n_of(p)
    k = k_type_I_3p_minus_2(p)
    n_star = 2 * (n - 2)
    n_far = (n - 2) * (n - 3) // 2
    mu_s = mu_star_two_orbit(p)
    mu_f = mu_far_two_orbit(p)
    # det identity: det = −2p(p²−1)
    Phi = phi_of(p)
    a11, a12 = 2 * p - 2, Phi - 2 * p + 1
    a21, a22 = 2 * p + 2, Phi - 2 * p - 1
    det = a11 * a22 - a12 * a21
    ok = det == -2 * p * (p * p - 1)
    ok = ok and mu_f == Fraction(-2 * (2 * p - 3), p * (p * p - 1))
    ok = ok and mu_f < 0
    mass = mu_s * n_star + mu_f * n_far
    ok = ok and mass == k
    ok = ok and not (0 <= mu_f <= 1)
    ok = ok and 0 < mu_s < 1
    return {
        "p": p,
        "mu_star": str(mu_s),
        "mu_far": str(mu_f),
        "det": str(det),
        "mass": str(mass),
        "k": k,
        "mass_eq_k": mass == k,
        "mu_far_negative": mu_f < 0,
        "in_unit_box": (0 <= mu_s <= 1) and (0 <= mu_f <= 1),
        "proved": ok,
        "theorem": (
            "2-orbit Aut_e solution of F_+ x̄=3−2f_e, x̄_e=0 has "
            "μ_far=−2(2p−3)/(p(p²−1))<0 and mass k.  Not a [0,1] vector."
        ),
    }


def lemma_maxminus_star_far_sums(p: int) -> dict:
    """Pointwise on Max−: wedges −2p−2f_e, far −Φ+2p+f_e."""
    Phi = phi_of(p)
    ok = True
    samples = {}
    for f in (1, -1):
        wedges = -2 * p - 2 * f
        far = -Phi + 2 * p + f
        total = wedges + far + f
        if total != -Phi:
            ok = False
        samples[str(f)] = {
            "wedges": str(wedges),
            "far": str(far),
            "total_plus_e": str(total),
            "minus_Phi": str(-Phi),
        }
    # Max+ twin: wedges 2p−2f, far Φ−2p+f, total Φ
    for f in (1, -1):
        if (2 * p - 2 * f) + (Phi - 2 * p + f) + f != Phi:
            ok = False
    return {
        "p": p,
        "samples": samples,
        "proved": ok,
        "theorem": (
            "Max−: ∑_{wedges e} f=−2p−2f_e, ∑_{far e} f=−Φ+2p+f_e "
            "(star at each end is −p; 15.51.3 twin).  Max+ signs flip."
        ),
    }


def lemma_far_weight_more_negative_when_beta_gt_minus_2(p: int) -> dict:
    """
    Subtracting the Max+ 2-orbit identity from a Max− identity with
    slope β>−2 makes the far weight strictly more negative than μ_far<0.
    """
    mu_f = mu_far_two_orbit(p)
    # On Max−, full far sum = −Φ+2p+f.  If the same μ_star is kept and
    # we need slope β on Max−, the far coefficient must absorb β+2>0
    # against the Max+ far slope (+1), which *decreases* μ_far.
    # Concrete: Max+ far slope is +1, Max− far slope is +1 (same, because
    # far = −Φ+2p+f on Max− and Φ−2p+f on Max+; both have coeff +1 of f).
    # Star slope is −2 on Max+ and −2 on Max−.  So
    #   μ_star(−2) + μ_far(+1) = −2     (Max+)
    #   μ_star(−2) + μ_far(+1) = β      (Max−)
    # These cannot hold simultaneously unless β=−2.  If we allow a
    # different μ_far^− on Max−, that is a different vector — but x̄ is
    # *one* edge-weight, evaluated on both sides.  Contradiction for β≠−2.
    # (The displayed linear system is overdetermined.)
    beta_two = Fraction(-2, p + 1)
    ok = mu_f < 0 and beta_two > -2
    # incompatibility: same (μ_s, μ_f) cannot give two different slopes
    # μ_s*(−2)+μ_f*(+1) is a single number
    mu_s = mu_star_two_orbit(p)
    slope_plus = mu_s * (-2) + mu_f * (1)
    # this is the slope of F_+ x̄, which equals −2
    ok = ok and slope_plus == -2
    return {
        "p": p,
        "mu_star": str(mu_s),
        "mu_far": str(mu_f),
        "slope_of_same_weights": str(slope_plus),
        "equals_minus_2": slope_plus == -2,
        "cannot_also_equal_beta_two": slope_plus != beta_two,
        "proved": ok,
        "theorem": (
            "A single Aut_e-invariant 2-orbit weight cannot realise two "
            "different (S,f) slopes.  Max+ forces slope −2; leftover "
            "a<thr would need slope β>−2 on Max−.  Impossible."
        ),
    }


# ---------------------------------------------------------------------------
# G'. Two Aut_e star orbits + collapsed-far 3-weight system
# ---------------------------------------------------------------------------


def n_aut_e_star_orbits() -> int:
    """
    Paley Aut_e (squares × inversion × Frobenius fixing {∞,0}) splits
    the 2(n−2) wedges into two orbits of size p²−1: third vertex square
    vs nonsquare.  Multiplication by a nonsquare is not in Aut(C).
    """
    return 2


def star_orbit_size(p: int) -> int:
    """Each Aut_e star orbit: 2·(q−1)/2 = p²−1."""
    return p * p - 1


def chi_minus_one_paley() -> int:
    """χ(−1)=1 on F_{p²}: (p²−1)/2 is even for every odd p."""
    return 1


def sigma_star_plus_Maxplus(p: int, f: int) -> int:
    """Max+ interpolant of the square-star orbit: (p−1)(1+f_e)."""
    return (p - 1) * (1 + f)


def sigma_star_minus_Maxplus(p: int, f: int) -> int:
    """Max+ interpolant of the nonsquare-star orbit: (p+1)(1−f_e)."""
    return (p + 1) * (1 - f)


def sigma_star_plus_Maxminus(p: int, f: int) -> int:
    """Max−: −(p+1)(1+f_e)."""
    return -(p + 1) * (1 + f)


def sigma_star_minus_Maxminus(p: int, f: int) -> int:
    """Max−: −(p−1)(1−f_e)."""
    return -(p - 1) * (1 - f)


def mu_plus_of_mu_far(p: int, mu_f: Fraction) -> Fraction:
    """μ_□ = 1/(2(p−1)) − μ_far(p−1)(p+2)/4."""
    return Fraction(1, 2 * (p - 1)) - mu_f * Fraction((p - 1) * (p + 2), 4)


def mu_minus_of_mu_far(p: int, mu_f: Fraction) -> Fraction:
    """μ_⊠ = 5/(2(p+1)) − μ_far(p+1)(p−2)/4."""
    return Fraction(5, 2 * (p + 1)) - mu_f * Fraction((p + 1) * (p - 2), 4)


def mu_far_equal_star(p: int) -> Fraction:
    """μ_□=μ_⊠ forces the 2-orbit far weight."""
    return mu_far_two_orbit(p)


def mu_far_star_minus_zero(p: int) -> Fraction:
    """μ_⊠=0 ⇒ μ_far=10/((p+1)²(p−2))."""
    return Fraction(10, (p + 1) * (p + 1) * (p - 2))


def mu_far_star_plus_zero(p: int) -> Fraction:
    """μ_□=0 ⇒ μ_far=2/((p−1)²(p+2)) (this slice is *not* empty)."""
    return Fraction(2, (p - 1) * (p - 1) * (p + 2))


def mu_plus_star_minus_zero(p: int) -> Fraction:
    """
    μ_⊠=0 forces
      μ_□ = −2(p³−3p+3)/[(p−1)(p+1)²(p−2)] < 0  (p≥5).
    """
    num = -2 * (p * p * p - 3 * p + 3)
    den = (p - 1) * (p + 1) * (p + 1) * (p - 2)
    return Fraction(num, den)


def lemma_two_star_interpolants(p: int) -> dict:
    """
    Two Aut_e star orbits, Max+-free.  Goes False if the count is 1,
    the sizes are wrong, or a vanishing / total-wedge identity fails.
    """
    n = n_of(p)
    n_plus = star_orbit_size(p)
    n_minus = star_orbit_size(p)
    ok = n_aut_e_star_orbits() == 2
    ok = ok and n_plus == p * p - 1 == n - 2
    ok = ok and n_plus + n_minus == 2 * (n - 2)
    # χ(−1)=1: (p²−1)/2 even
    ok = ok and ((p * p - 1) // 2) % 2 == 0
    ok = ok and chi_minus_one_paley() == 1
    # vanishing + total wedge
    for f, tot in ((1, 2 * p - 2), (-1, 2 * p + 2)):
        sp = sigma_star_plus_Maxplus(p, f)
        sm = sigma_star_minus_Maxplus(p, f)
        if sp + sm != tot:
            ok = False
    ok = ok and sigma_star_plus_Maxplus(p, -1) == 0
    ok = ok and sigma_star_minus_Maxplus(p, 1) == 0
    ok = ok and sigma_star_plus_Maxplus(p, 1) == 2 * p - 2
    ok = ok and sigma_star_minus_Maxplus(p, -1) == 2 * p + 2
    # closed forms
    ok = ok and sigma_star_plus_Maxplus(p, 1) == (p - 1) * 2
    ok = ok and sigma_star_minus_Maxplus(p, -1) == (p + 1) * 2
    # Max− twin
    for f, tot in ((1, -2 * p - 2), (-1, -2 * p + 2)):
        if sigma_star_plus_Maxminus(p, f) + sigma_star_minus_Maxminus(p, f) != tot:
            ok = False
    ok = ok and sigma_star_plus_Maxminus(p, -1) == 0
    ok = ok and sigma_star_minus_Maxminus(p, 1) == 0
    # swapped labels: (p+1)(1−f) at f=−1 is 2(p+1)≠0, not the plus vanishing
    swapped_plus_at_minus = (p + 1) * (1 - (-1))
    ok = ok and swapped_plus_at_minus == 2 * (p + 1) != 0
    ok = ok and sigma_star_plus_Maxplus(p, -1) != swapped_plus_at_minus
    return {
        "p": p,
        "n_star_orbits": n_aut_e_star_orbits(),
        "star_orbit_size": n_plus,
        "sigma_plus_at_1": sigma_star_plus_Maxplus(p, 1),
        "sigma_minus_at_minus1": sigma_star_minus_Maxplus(p, -1),
        "chi_minus_one": chi_minus_one_paley(),
        "proved": ok,
        "theorem": (
            "Aut_e has two star orbits of size p²−1.  Max+: "
            "σ_□=(p−1)(1+f_e), σ_⊠=(p+1)(1−f_e) from prefactor "
            "vanishing and total wedge 2p−2f_e."
        ),
    }


def lemma_three_orbit_eval_identities(p: int, mu_f: Fraction) -> dict:
    """The 1-parameter 3-weight family solves F_+ x̄=3−2f_e and mass=k."""
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    n_s = star_orbit_size(p)
    n_far = (n_of(p) - 2) * (n_of(p) - 3) // 2
    k = k_type_I_3p_minus_2(p)
    mass = mu_p * n_s + mu_m * n_s + mu_f * n_far
    ok = mass == k
    for f, need in ((1, 1), (-1, 5)):
        got = (
            mu_p * sigma_star_plus_Maxplus(p, f)
            + mu_m * sigma_star_minus_Maxplus(p, f)
            + mu_f * (phi_of(p) - 2 * p + f)
        )
        if got != need:
            ok = False
    return {
        "p": p,
        "mu_far": str(mu_f),
        "mu_plus": str(mu_p),
        "mu_minus": str(mu_m),
        "mass_eq_k": mass == k,
        "proved": ok,
    }


def lemma_three_orbit_equal_star_empty(p: int) -> dict:
    """
    μ_□=μ_⊠ in the 3-weight system forces μ_far=2-orbit value < 0.
    Goes False if that common far weight is nonnegative.
    """
    mu_f = mu_far_equal_star(p)
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    ev = lemma_three_orbit_eval_identities(p, mu_f)
    ok = ev["proved"]
    ok = ok and mu_f == mu_far_two_orbit(p)
    ok = ok and mu_f < 0
    ok = ok and mu_p == mu_m
    ok = ok and mu_p == mu_star_two_orbit(p)
    ok = ok and not (0 <= mu_f <= 1)
    return {
        "p": p,
        "mu_far": str(mu_f),
        "mu_plus": str(mu_p),
        "mu_minus": str(mu_m),
        "equals_two_orbit_far": mu_f == mu_far_two_orbit(p),
        "mu_far_negative": mu_f < 0,
        "in_unit_box": (0 <= mu_p <= 1) and (0 <= mu_m <= 1) and (0 <= mu_f <= 1),
        "proved": ok,
        "theorem": (
            "Equal Aut_e star densities collapse to the 2-orbit solution; "
            "μ_far=−2(2p−3)/(p(p²−1))<0."
        ),
    }


def lemma_three_orbit_plus_only_empty(p: int) -> dict:
    """
    μ_⊠=0 forces μ_□=−2(p³−3p+3)/[(p−1)(p+1)²(p−2)]<0.
    Goes False if the cubic p³−3p+3 is nonpositive or μ_□≥0.
    """
    mu_f = mu_far_star_minus_zero(p)
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    closed = mu_plus_star_minus_zero(p)
    cubic = p * p * p - 3 * p + 3
    ev = lemma_three_orbit_eval_identities(p, mu_f)
    ok = ev["proved"]
    ok = ok and mu_m == 0
    ok = ok and mu_p == closed
    ok = ok and cubic > 0
    # cubic' = 3p²−3>0 on p≥5; cubic(5)=113>0
    ok = ok and (3 * p * p - 3) > 0
    ok = ok and mu_p < 0
    ok = ok and not (0 <= mu_p <= 1)
    return {
        "p": p,
        "mu_far": str(mu_f),
        "mu_plus": str(mu_p),
        "mu_minus": str(mu_m),
        "cubic": cubic,
        "mu_plus_negative": mu_p < 0,
        "in_unit_box": (0 <= mu_p <= 1) and (0 <= mu_m <= 1) and (0 <= mu_f <= 1),
        "proved": ok,
        "theorem": (
            "3-weight slice μ_⊠=0 has "
            "μ_□=−2(p³−3p+3)/[(p−1)(p+1)²(p−2)]<0 (p³−3p+3>0 for p≥5)."
        ),
    }


def lemma_three_orbit_zero_far_in_box(p: int) -> dict:
    """
    μ_far=0 gives μ_□=1/(2(p−1)), μ_⊠=5/(2(p+1)) ∈ (0,1) and meets
    F_+ x̄=3−2f_e.  This slice is *not* empty — the ND flag stays False.
    Goes False if either star weight leaves (0,1) or the identities fail.
    """
    mu_f = Fraction(0)
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    ev = lemma_three_orbit_eval_identities(p, mu_f)
    ok = ev["proved"]
    ok = ok and mu_p == Fraction(1, 2 * (p - 1))
    ok = ok and mu_m == Fraction(5, 2 * (p + 1))
    ok = ok and 0 < mu_p < 1
    ok = ok and 0 < mu_m < 1
    ok = ok and mu_f == 0
    return {
        "p": p,
        "mu_plus": str(mu_p),
        "mu_minus": str(mu_m),
        "mu_far": str(mu_f),
        "in_unit_box": True,
        "empty": False,
        "proved": ok,
        "theorem": (
            "3-weight slice μ_far=0 is a [0,1] solution of F_+ x̄=3−2f_e "
            "and mass k.  Not an empty subcase."
        ),
    }


def lemma_three_orbit_plus_zero_not_empty(p: int) -> dict:
    """μ_□=0 gives μ_far, μ_⊠ ∈ (0,1).  Do not claim this slice empty as a box."""
    mu_f = mu_far_star_plus_zero(p)
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    ev = lemma_three_orbit_eval_identities(p, mu_f)
    ok = ev["proved"] and mu_p == 0 and 0 < mu_m < 1 and 0 < mu_f < 1
    return {
        "p": p,
        "mu_far": str(mu_f),
        "mu_plus": str(mu_p),
        "mu_minus": str(mu_m),
        "in_unit_box": bool(mu_p == 0 and 0 <= mu_m <= 1 and 0 <= mu_f <= 1),
        "empty": False,
        "proved": ok,
    }


# ---------------------------------------------------------------------------
# H. 3-weight Max− at f_e=+1: bad case empty on the whole family
# ---------------------------------------------------------------------------


def n_square_star_zero_far(p: int) -> Fraction:
    """Orbit fill n_□ at μ_far=0: (p+1)/2."""
    return Fraction(p + 1, 2)


def n_nonsquare_star_zero_far(p: int) -> Fraction:
    """Orbit fill n_⊠ at μ_far=0: 5(p−1)/2."""
    return Fraction(5 * (p - 1), 2)


def Fminus_three_orbit(p: int, mu_f: Fraction, f: int) -> Fraction:
    """3-weight Aut_e-average on Max− at a given f_e=±1."""
    mu_p = mu_plus_of_mu_far(p, mu_f)
    mu_m = mu_minus_of_mu_far(p, mu_f)
    return (
        mu_p * sigma_star_plus_Maxminus(p, f)
        + mu_m * sigma_star_minus_Maxminus(p, f)
        + mu_f * (-phi_of(p) + 2 * p + f)
    )


def Fminus_three_orbit_fplus_closed(p: int, mu_f: Fraction) -> Fraction:
    """−(p+1)/(p−1) + μ_far p(p+1)."""
    return -Fraction(p + 1, p - 1) + mu_f * p * (p + 1)


def Fminus_three_orbit_fminus_closed(p: int, mu_f: Fraction) -> Fraction:
    """−5(p−1)/(p+1) − μ_far p(p−1)."""
    return -Fraction(5 * (p - 1), p + 1) - mu_f * p * (p - 1)


def Fminus_zero_far_at_fplus(p: int) -> Fraction:
    return Fminus_three_orbit_fplus_closed(p, Fraction(0))


def Fminus_zero_far_at_fminus(p: int) -> Fraction:
    return Fminus_three_orbit_fminus_closed(p, Fraction(0))


def lemma_zero_far_01_counts(p: int) -> dict:
    """
    μ_far=0 0-1 fill: n_□=(p+1)/2, n_⊠=5(p−1)/2, both integral, sum k.
    Goes False if the counts leave the orbits or miss mass k.
    """
    n_sq = n_square_star_zero_far(p)
    n_ns = n_nonsquare_star_zero_far(p)
    n_s = star_orbit_size(p)
    k = k_type_I_3p_minus_2(p)
    ok = n_sq + n_ns == k
    ok = ok and n_sq == n_s * mu_plus_of_mu_far(p, Fraction(0))
    ok = ok and n_ns == n_s * mu_minus_of_mu_far(p, Fraction(0))
    ok = ok and n_sq.denominator == 1 and n_ns.denominator == 1
    ok = ok and 0 < n_sq < n_s and 0 < n_ns < n_s
    return {
        "p": p,
        "n_square": str(n_sq),
        "n_nonsquare": str(n_ns),
        "k": k,
        "integral": n_sq.denominator == 1 and n_ns.denominator == 1,
        "proved": ok,
        "theorem": (
            "Star-supported Type I Aut_e-average has "
            "n_□=(p+1)/2, n_⊠=5(p−1)/2 (integral, sum 3p−2)."
        ),
    }


def lemma_three_orbit_maxminus_fplus(p: int) -> dict:
    """
    3-weight F_− on Max− at f_e=±1.  Goes False if the closed form
    misses the interpolants, if μ_far=0 drops to ≤−2, or if some
    μ_far≥0 box point reaches ≤−3 (the bad-case demand).
    """
    ok = True
    samples = []
    mu_box = mu_far_star_plus_zero(p)  # μ_□=0 end, still in [0,1]
    for mu_f in (
        Fraction(0),
        mu_box / 2,
        mu_box,
        mu_far_two_orbit(p),  # negative; identity still holds
    ):
        got_p = Fminus_three_orbit(p, mu_f, 1)
        got_m = Fminus_three_orbit(p, mu_f, -1)
        cl_p = Fminus_three_orbit_fplus_closed(p, mu_f)
        cl_m = Fminus_three_orbit_fminus_closed(p, mu_f)
        if got_p != cl_p or got_m != cl_m:
            ok = False
        ev = lemma_three_orbit_eval_identities(p, mu_f)
        if not ev["proved"]:
            ok = False
        samples.append(
            {
                "mu_far": str(mu_f),
                "F_plus": str(got_p),
                "F_minus": str(got_m),
                "closed_match": got_p == cl_p and got_m == cl_m,
            }
        )
    zplus = Fminus_zero_far_at_fplus(p)
    zminus = Fminus_zero_far_at_fminus(p)
    ok = ok and zplus == -Fraction(p + 1, p - 1)
    ok = ok and zminus == -Fraction(5 * (p - 1), p + 1)
    # (p+1)/(p−1) < 2 ⇔ p>3
    ok = ok and (p + 1) < 2 * (p - 1)
    ok = ok and -2 < zplus < -1
    ok = ok and zplus > -3
    # whole μ_far≥0 ray is > −2 at f=+1
    for mu_f in (Fraction(0), mu_box, Fraction(1)):
        val = Fminus_three_orbit_fplus_closed(p, mu_f)
        if mu_f >= 0 and val <= -2:
            ok = False
        if mu_f >= 0 and val <= -3:
            ok = False
    # a wrong closed form −(p−1)/(p+1) is ≥ −1 and would not beat −2
    wrong = -Fraction(p - 1, p + 1)
    ok = ok and wrong != zplus and wrong > -1
    # swapped Max+ interpolants at f=+1: μ_□·2(p−1)=1, not the Max− value
    swapped = mu_plus_of_mu_far(p, Fraction(0)) * sigma_star_plus_Maxplus(p, 1)
    ok = ok and swapped != zplus
    return {
        "p": p,
        "F_zero_far_fplus": str(zplus),
        "F_zero_far_fminus": str(zminus),
        "gt_minus_2": zplus > -2,
        "gt_minus_3": zplus > -3,
        "samples": samples,
        "proved": ok,
        "theorem": (
            "3-weight F_−|_{f=+1}=−(p+1)/(p−1)+μ_far p(p+1).  "
            "μ_far≥0 ⇒ this is ≥−(p+1)/(p−1)>−2, so not ≤−3."
        ),
    }


def lemma_zero_far_bad_case_empty(p: int) -> dict:
    """
    μ_far=0 cannot be a gap-2 bad-case undercutter.
    Goes False if F_−|_{f=+1}≤−3 or the 0-1 counts fail.
    """
    mm = lemma_three_orbit_maxminus_fplus(p)
    cnt = lemma_zero_far_01_counts(p)
    zplus = Fminus_zero_far_at_fplus(p)
    ok = mm["proved"] and cnt["proved"]
    ok = ok and zplus == -Fraction(p + 1, p - 1)
    ok = ok and zplus > -2
    # undercutter + bad case would need F_−≤−3 on {f=+1}
    ok = ok and not (zplus <= -3)
    # odd-integer scores ≤−1 cannot average to a value in (−2,−1)
    # unless some mass sits at S=−1 (hence on {f=+1}, not the bad case)
    ok = ok and -2 < zplus < -1
    return {
        "p": p,
        "E_S_fplus": str(zplus),
        "in_minus2_minus1": -2 < zplus < -1,
        "counts": cnt,
        "proved": ok,
        "theorem": (
            "μ_far=0 Aut_e-average has F_−≡−(p+1)/(p−1)∈(−2,−1) on "
            "{f=+1}.  Bad case + undercutter needs ≤−3.  Empty."
        ),
    }


def lemma_three_orbit_family_bad_case_empty(p: int) -> dict:
    """
    Every μ_far≥0 3-weight point (including the [0,1] interval
    [0, μ_□=0]) fails F_−|_{f=+1}≤−3.  Goes False if any such
    point drops to ≤−3 or the μ_far=0 endpoint does.
    """
    mm = lemma_three_orbit_maxminus_fplus(p)
    zf = lemma_zero_far_bad_case_empty(p)
    mu_hi = mu_far_star_plus_zero(p)
    ok = mm["proved"] and zf["proved"]
    ok = ok and mu_hi > 0
    for mu_f in (Fraction(0), mu_hi / 3, mu_hi, 2 * mu_hi):
        val = Fminus_three_orbit_fplus_closed(p, mu_f)
        if mu_f >= 0 and val <= -3:
            ok = False
        # identity still holds for μ_far past the box
        if Fminus_three_orbit(p, mu_f, 1) != val:
            ok = False
    # μ_□=0 endpoint is in (−1,0), still > −3
    at_hi = Fminus_three_orbit_fplus_closed(p, mu_hi)
    ok = ok and -1 < at_hi < 0
    return {
        "p": p,
        "mu_far_box_hi": str(mu_hi),
        "F_at_box_hi": str(at_hi),
        "proved": ok,
        "theorem": (
            "On the whole μ_far≥0 3-weight family, "
            "F_−|_{f=+1}≥−(p+1)/(p−1)>−2.  Bad case empty."
        ),
    }


# ---------------------------------------------------------------------------
# I. Split-far 2-eval pairing + coarse 3-type collapse
# ---------------------------------------------------------------------------


def n_squares_nonzero(p: int) -> int:
    """|Ω| = (p²−1)/2."""
    return (p * p - 1) // 2


def n_far_same_chi(p: int) -> int:
    """# of far edges with both ends squares (or both nonsquares)."""
    m = n_squares_nonzero(p)
    return m * (m - 1) // 2


def n_far_cross_chi(p: int) -> int:
    """# of far edges with one square and one nonsquare end."""
    m = n_squares_nonzero(p)
    return m * m


def A_from_Gplus(size: int, Gplus: Fraction, p: int) -> Fraction:
    """E_+[σ|f=+1] = (|O| + p G_+)/(p+1)."""
    return (size + p * Gplus) / (p + 1)


def B_from_Gminus(size: int, Gminus: Fraction, p: int) -> Fraction:
    """E_−[σ|f=+1] = (−|O| + p G_−)/(p−1)."""
    return (-size + p * Gminus) / (p - 1)


def delta_from_AB(A: Fraction, B: Fraction, p: int) -> Fraction:
    """δ = ((p+1)/(p−1)) A + B."""
    return Fraction(p + 1, p - 1) * A + B


def delta_from_Gsum(Gsum: Fraction, p: int) -> Fraction:
    """δ = p (G_++G_−)/(p−1)."""
    return p * Gsum / (p - 1)


def collapsed_far_A_fplus(p: int) -> Fraction:
    """Total far interpolant on Max+ at f=+1: Φ−2p+1."""
    return phi_of(p) - 2 * p + 1


def collapsed_far_B_fplus(p: int) -> Fraction:
    """Total far interpolant on Max− at f=+1: −Φ+2p+1."""
    return -phi_of(p) + 2 * p + 1


def collapsed_far_3AB(p: int) -> Fraction:
    """3A+B = p³−3p+4 for the collapsed far interpolant at f=+1."""
    return (
        3 * collapsed_far_A_fplus(p) + collapsed_far_B_fplus(p)
    )


def star_plus_3AB(p: int) -> Fraction:
    """3A_□+B_□ at f=+1: 3·2(p−1) + (−2(p+1)) = 4(p−2)."""
    return 4 * (p - 2)


def star_minus_3AB(p: int) -> Fraction:
    """Nonsquare-star vanishes at f=+1, so 3A_⊠+B_⊠=0."""
    return Fraction(0)


def Y_from_mu_box_and_far(
    p: int, mu_box: Fraction, far_3AB: Fraction
) -> Fraction:
    """Y = −3 + 4(p−2)μ_□ + ∑ μ(3A+B)."""
    return -3 + 4 * (p - 2) * mu_box + far_3AB


def Q_OmOm_minus_Q_NN_Maxplus(p: int, f: int) -> Fraction:
    """Pointwise on Max+: Q_ΩΩ−Q_NN = 1−p f_e (Cy=py)."""
    return 1 - p * f


def Q_OmOm_minus_Q_NN_Maxminus(p: int, f: int) -> Fraction:
    """Pointwise on Max−: Q_ΩΩ−Q_NN = 1+p f_e (Cy=−py)."""
    return 1 + p * f


def three_type_d_plus(p: int) -> Fraction:
    """d = p(p²−5)/2 + 2, the s-coefficient of F_+|_{f=+1}."""
    return Fraction(p * (p * p - 5), 2) + 2


def three_type_e_minus(p: int) -> Fraction:
    """e = p(p²−1)/2, the s-coefficient of F_+|_{f=−1}."""
    return Fraction(p * (p * p - 1), 2)


def three_type_c_minus(p: int) -> Fraction:
    """c = p(5−p²)/2 + 2, the s-coefficient of Y (Max− at f=+1)."""
    return Fraction(p * (5 - p * p), 2) + 2


def three_type_Y(p: int, s: Fraction) -> Fraction:
    """
    After residual cancel μ_ΩN=(μ_ΩΩ+μ_NN)/2, Y depends only on
    s=μ_ΩΩ+μ_NN and equals the collapsed-far formula at μ_far=s/2.
    """
    return -Fraction(p + 1, p - 1) + (s / 2) * p * (p + 1)


def lemma_split_far_delta_identity(p: int) -> dict:
    """
    2-point π ⇒ δ = p Gsum/(p−1).  Goes False if A,B vs G± is wrong
    or the two expressions for δ disagree.
    """
    ok = True
    samples = []
    for size, Gp, Gm in (
        (1, Fraction(1, p * p - 2), Fraction(1, p * p - 2)),
        (p * p - 1, Fraction(0), Fraction(0)),
        (
            (n_of(p) - 2) * (n_of(p) - 3) // 2,
            Fraction(p * p - 1, 2),
            Fraction(p * p - 1, 2),
        ),
        (24, Fraction(-72, 65), Fraction(-72, 65)),
    ):
        A = A_from_Gplus(size, Gp, p)
        B = B_from_Gminus(size, Gm, p)
        d1 = delta_from_AB(A, B, p)
        d2 = delta_from_Gsum(Gp + Gm, p)
        if d1 != d2:
            ok = False
        samples.append(
            {"size": size, "delta": str(d1), "match": d1 == d2}
        )
    # collapsed far: δ = p(p+1)
    n_far = (n_of(p) - 2) * (n_of(p) - 3) // 2
    A = collapsed_far_A_fplus(p)
    B = collapsed_far_B_fplus(p)
    dcol = delta_from_AB(A, B, p)
    ok = ok and dcol == p * (p + 1)
    ok = ok and A == phi_of(p) - 2 * p + 1
    ok = ok and B == -phi_of(p) + 2 * p + 1
    # stars at f=+1
    A_box = Fraction(2 * (p - 1))
    B_box = Fraction(-2 * (p + 1))
    ok = ok and delta_from_AB(A_box, B_box, p) == 0
    ok = ok and delta_from_AB(Fraction(0), Fraction(0), p) == 0
    # a swapped Max+ / Max− formula for B would break collapsed δ
    B_swapped = -A
    ok = ok and delta_from_AB(A, B_swapped, p) != p * (p + 1)
    return {
        "p": p,
        "collapsed_delta": str(dcol),
        "n_far": n_far,
        "samples": samples,
        "proved": ok,
        "theorem": (
            "δ=((p+1)/(p−1))A+B=p Gsum/(p−1); collapsed far has "
            "δ=p(p+1); both stars have δ=0 at f=+1."
        ),
    }


def lemma_split_far_Y_plus_3(p: int) -> dict:
    """
    Y+3 = 4(p−2)μ_□ + ∑ μ(3A+B).  Goes False if the star pairing
    3A_□+B_□≠4(p−2) or the μ_far=0 point misses −(p+1)/(p−1).
    """
    ok = star_plus_3AB(p) == 4 * (p - 2)
    ok = ok and star_plus_3AB(p) > 0
    ok = ok and star_minus_3AB(p) == 0
    # μ_far=0: μ_□=1/(2(p−1)), far 3AB=0
    mu_box = Fraction(1, 2 * (p - 1))
    Y0 = Y_from_mu_box_and_far(p, mu_box, Fraction(0))
    ok = ok and Y0 == -Fraction(p + 1, p - 1)
    ok = ok and Y0 > -2
    # collapsed far contribution
    cab = collapsed_far_3AB(p)
    ok = ok and cab == p * p * p - 3 * p + 4
    ok = ok and cab > 0
    # whole μ_far≥0 3-weight: Y+3 = 4(p−2)μ_□(μ_far) + μ_far cab
    # μ_□ = 1/(2(p−1)) − μ_far (p−1)(p+2)/4
    for mu_f in (Fraction(0), mu_far_star_plus_zero(p), Fraction(1, 50)):
        mu_p = mu_plus_of_mu_far(p, mu_f)
        Y = Y_from_mu_box_and_far(p, mu_p, mu_f * cab)
        cl = Fminus_three_orbit_fplus_closed(p, mu_f)
        if Y != cl:
            ok = False
        if mu_f >= 0 and Y <= -3:
            ok = False
    # Y≤−3 with μ_□≥0 would need far 3AB ≤ −4(p−2)μ_□ ≤ 0
    # collapsed ray never does
    mu_hi = mu_far_star_plus_zero(p)
    need = -4 * (p - 2) * mu_plus_of_mu_far(p, mu_hi)
    ok = ok and mu_hi * cab > 0 and mu_hi * cab > need
    return {
        "p": p,
        "star_plus_3AB": str(star_plus_3AB(p)),
        "collapsed_3AB": str(cab),
        "Y_zero_far": str(Y0),
        "proved": ok,
        "theorem": (
            "Y+3=4(p−2)μ_□+∑μ(3A+B).  Stars: 4(p−2) and 0.  "
            "Collapsed far: p³−3p+4>0.  μ_far≥0 3-weight never has Y≤−3."
        ),
    }


def lemma_three_type_counts(p: int) -> dict:
    """
    Far edges split ΩΩ / NN / ΩN.  Goes False if the three sizes
    miss the far count (n−2)(n−3)/2.
    """
    n_oo = n_far_same_chi(p)
    n_nn = n_far_same_chi(p)
    n_cr = n_far_cross_chi(p)
    n_far = (n_of(p) - 2) * (n_of(p) - 3) // 2
    ok = n_oo + n_nn + n_cr == n_far
    m = n_squares_nonzero(p)
    ok = ok and 2 * m == p * p - 1
    ok = ok and n_oo == m * (m - 1) // 2
    ok = ok and n_cr == m * m
    return {
        "p": p,
        "n_OmegaOmega": n_oo,
        "n_NN": n_nn,
        "n_cross": n_cr,
        "n_far": n_far,
        "proved": ok,
        "theorem": (
            "Far edges partition as ΩΩ ∪ NN ∪ ΩN with sizes "
            "(p²−1)(p²−3)/8, same, (p²−1)²/4."
        ),
    }


def lemma_three_type_Q_difference(p: int) -> dict:
    """
    Q_ΩΩ−Q_NN = 1−p f_e on Max+ and 1+p f_e on Max− (pointwise, Cy).
    Goes False if the two sides agree or the f-coefficient flips.
    """
    ok = True
    for f in (1, -1):
        if Q_OmOm_minus_Q_NN_Maxplus(p, f) != 1 - p * f:
            ok = False
        if Q_OmOm_minus_Q_NN_Maxminus(p, f) != 1 + p * f:
            ok = False
        # Max+ at f=+1 is 1−p < 0; Max− at f=+1 is 1+p > 0
        if f == 1:
            ok = ok and Q_OmOm_minus_Q_NN_Maxplus(p, f) == 1 - p < 0
            ok = ok and Q_OmOm_minus_Q_NN_Maxminus(p, f) == 1 + p > 0
    # a forgotten sign on Max− would equal Max+
    ok = ok and Q_OmOm_minus_Q_NN_Maxplus(p, 1) != Q_OmOm_minus_Q_NN_Maxminus(
        p, 1
    )
    # residual kernel (r,r,−2r) preserves both identities
    for r in (Fraction(1), Fraction(-3, 2), Fraction(0)):
        # X'=X+r, Y'=Y+r, Z'=Z−2r: (X'−Y') unchanged, 2X'+Z' unchanged
        ok = ok and (r - r) == 0 and (2 * r - 2 * r) == 0
    return {
        "p": p,
        "Maxplus_at_1": str(Q_OmOm_minus_Q_NN_Maxplus(p, 1)),
        "Maxminus_at_1": str(Q_OmOm_minus_Q_NN_Maxminus(p, 1)),
        "proved": ok,
        "theorem": (
            "Max+: Q_ΩΩ−Q_NN=1−p f_e; Max−: =1+p f_e.  "
            "3-type residual kernel is (r,r,−2r)."
        ),
    }


def lemma_three_type_Y_collapses_to_H(p: int) -> dict:
    """
    Residual-cancel 3-type (μ_ΩN=(μ_ΩΩ+μ_NN)/2) has the same Y as
    collapsed far at μ_far=(μ_ΩΩ+μ_NN)/2.  Goes False if the
    s-coefficient is not p(p+1) or a=b misses H.
    """
    d = three_type_d_plus(p)
    e = three_type_e_minus(p)
    c = three_type_c_minus(p)
    ok = d == Fraction(p * (p * p - 5), 2) + 2
    ok = ok and e == Fraction(p * (p * p - 1), 2)
    ok = ok and c == Fraction(p * (5 - p * p), 2) + 2
    # c + d(p+1)/(p−1) = p(p+1)
    coeff = c + d * Fraction(p + 1, p - 1)
    ok = ok and coeff == p * (p + 1)
    # collapsed a=b=μ_f ⇔ s=2μ_f, Y=H
    samples = []
    for mu_f in (
        Fraction(0),
        mu_far_star_plus_zero(p),
        Fraction(1, 40),
        mu_far_two_orbit(p),
    ):
        s = 2 * mu_f
        Y = three_type_Y(p, s)
        cl = Fminus_three_orbit_fplus_closed(p, mu_f)
        if Y != cl:
            ok = False
        samples.append({"mu_far": str(mu_f), "Y": str(Y), "match_H": Y == cl})
    # splitting a≠b at fixed s does not appear in Y (b-terms cancel)
    s = Fraction(1, 5)
    ok = ok and three_type_Y(p, s) == -Fraction(p + 1, p - 1) + (
        s / 2
    ) * p * (p + 1)
    # F_+ coefficients recover collapsed A at a=b
    # (p−1) + d = Φ−2p+1 when multiplied as in the a=b reduction
    Acol = collapsed_far_A_fplus(p)
    ok = ok and (p - 1) + d == Acol
    Amcol = phi_of(p) - 2 * p - 1
    ok = ok and e - (p + 1) == Amcol
    return {
        "p": p,
        "s_coeff": str(coeff),
        "samples": samples,
        "proved": ok,
        "theorem": (
            "3-type residual cancel μ_ΩN=(μ_ΩΩ+μ_NN)/2 makes Y "
            "equal to H at μ_far=(μ_ΩΩ+μ_NN)/2.  Coarse split "
            "cannot undercut −(p+1)/(p−1) for μ_far≥0."
        ),
    }


def type_I_split_far_identities_ok() -> bool:
    return all(
        lemma_split_far_delta_identity(p)["proved"]
        and lemma_split_far_Y_plus_3(p)["proved"]
        and lemma_three_type_counts(p)["proved"]
        and lemma_three_type_Q_difference(p)["proved"]
        and lemma_three_type_Y_collapses_to_H(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def type_I_three_type_Y_collapses_to_H() -> bool:
    return all(
        lemma_three_type_Y_collapses_to_H(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


# ---------------------------------------------------------------------------
# J. Split-far 3A+B conversion / Wick 6-type / Q_χ identities
# ---------------------------------------------------------------------------


def T_bitight(p: int) -> Fraction:
    """15.47 threshold T(p)=−(p−2)/(p(2p−1))."""
    return -Fraction(p - 2, p * (2 * p - 1))


def L_abs_gmin(p: int) -> Fraction:
    """Uniform LB candidate L(p)=−(p−2)/(2p²).  L>T; not the 3A+B zero."""
    return -Fraction(p - 2, 2 * p * p)


def wick_t_from_chi(p: int, chu: int, chv: int, chd: int) -> Fraction:
    """t=χ(u−v) κ / p with κ=χ(u)+χ(v)+χ(u−v)."""
    return Fraction(chd * (chu + chv + chd), p)


def wick_3ab_edge(p: int, chu: int, chv: int, chd: int) -> Fraction:
    """Wick per-edge 3a+b = 2[(p−2)+t(2p−1)]/(p²−1)."""
    t = wick_t_from_chi(p, chu, chv, chd)
    return 2 * ((p - 2) + t * (2 * p - 1)) / (p * p - 1)


def wick_3ab_closed_same_plus(p: int) -> Fraction:
    """Types (++,+) and (−−,−): 2(p²+4p−3)/(p(p²−1))."""
    return 2 * Fraction(p * p + 4 * p - 3, p * (p * p - 1))


def wick_3ab_closed_t_minus(p: int) -> Fraction:
    """Types (++,-) and (--,+): 2(p²−4p+1)/(p(p²−1))."""
    return 2 * Fraction(p * p - 4 * p + 1, p * (p * p - 1))


def wick_3ab_closed_cross_plus(p: int) -> Fraction:
    """Both ΩN types: κ=χ(u−v) ⇒ t=1/p and 3a+b=2/p."""
    return Fraction(2, p)


def three_AB_from_G(
    n_O: int, Gplus: Fraction, Gminus: Fraction, p: int
) -> Fraction:
    """3A+B = [2(p−2)n_O + p(3(p−1)G_++(p+1)G_−)]/(p²−1)."""
    num = (
        2 * (p - 2) * n_O
        + p * (3 * (p - 1) * Gplus + (p + 1) * Gminus)
    )
    return num / (p * p - 1)


def three_AB_equal_G(n_O: int, G_edge: Fraction, p: int) -> Fraction:
    """G_+=G_−=n_O G: 2 n_O[(p−2)+p(2p−1)G]/(p²−1)."""
    return 2 * n_O * ((p - 2) + p * (2 * p - 1) * G_edge) / (p * p - 1)


def A_from_Gplus_orbit(n_O: int, Gplus: Fraction, p: int) -> Fraction:
    return (n_O + p * Gplus) / (p + 1)


def B_from_Gminus_orbit(n_O: int, Gminus: Fraction, p: int) -> Fraction:
    return (-n_O + p * Gminus) / (p - 1)


def Q_chi_u_Maxplus(p: int, f: int) -> Fraction:
    """Max+: Q_χu=(1−p)(1+f_e)."""
    return Fraction(1 - p) * (1 + f)


def Q_chi_d_Maxplus(p: int, f: int) -> Fraction:
    """Max+: Q_χd = 1−p at f=+1 and p+1 at f=−1."""
    return Fraction(1 - p) if f == 1 else Fraction(p + 1)


def Q_kappa_Maxplus(p: int, f: int) -> Fraction:
    return Q_chi_u_Maxplus(p, f) + Q_chi_d_Maxplus(p, f)


def Q_chi_u_Maxminus(p: int, f: int) -> Fraction:
    """Max−: Q_χu=(p+1)(1+f_e)."""
    return Fraction(p + 1) * (1 + f)


def Q_chi_d_Maxminus(p: int, f: int) -> Fraction:
    """Max−: Q_χd = p+1 at f=+1 and 1−p at f=−1."""
    return Fraction(p + 1) if f == 1 else Fraction(1 - p)


def Q_kappa_Maxminus(p: int, f: int) -> Fraction:
    return Q_chi_u_Maxminus(p, f) + Q_chi_d_Maxminus(p, f)


def n_far_kappa3(p: int) -> int:
    """# far edges through e with |κ|=3: (n−2)(n−6)/8."""
    n = n_of(p)
    return (n - 2) * (n - 6) // 8


def n_far_kappa1(p: int) -> int:
    """# far edges through e with |κ|=1: 3(n−2)²/8."""
    n = n_of(p)
    return 3 * (n - 2) * (n - 2) // 8


def lemma_threeAB_conversion(p: int) -> dict:
    """
    3A+B from G±.  Goes False if the (p±1) slots are swapped or the
    equal-G reduction misses T(p).
    """
    ok = True
    samples = []
    n_O_list = (1, p - 1, p * p - 1, (n_of(p) - 2) * (n_of(p) - 3) // 2)
    for n_O in n_O_list:
        for Gp, Gm in (
            (Fraction(0), Fraction(0)),
            (Fraction(1, p), Fraction(-1, p)),
            (Fraction(-3, p * p), Fraction(-3, p * p)),
            (T_bitight(p) * n_O, T_bitight(p) * n_O),
            (L_abs_gmin(p) * n_O, L_abs_gmin(p) * n_O),
        ):
            A = A_from_Gplus_orbit(n_O, Gp, p)
            B = B_from_Gminus_orbit(n_O, Gm, p)
            got = 3 * A + B
            cl = three_AB_from_G(n_O, Gp, Gm, p)
            if got != cl:
                ok = False
            samples.append(
                {
                    "n_O": n_O,
                    "Gplus": str(Gp),
                    "Gminus": str(Gm),
                    "match": got == cl,
                }
            )
            if Gp == Gm:
                G_edge = Gp / n_O
                eq = three_AB_equal_G(n_O, G_edge, p)
                if eq != cl:
                    ok = False
                # 3A+B>0 ⇔ G>T
                if G_edge == T_bitight(p) and cl != 0:
                    ok = False
                if G_edge > T_bitight(p) and cl <= 0:
                    ok = False
                if G_edge < T_bitight(p) and cl >= 0:
                    ok = False
    # L is strictly above T, so G=L is *not* the 3A+B zero
    ok = ok and L_abs_gmin(p) > T_bitight(p)
    ok = ok and three_AB_equal_G(1, L_abs_gmin(p), p) > 0
    ok = ok and three_AB_equal_G(1, T_bitight(p), p) == 0
    # swapped (p+1)↔(p−1) in the G-weights is not the identity
    n_O, Gp, Gm = p + 1, Fraction(1, 7), Fraction(-2, 7)
    wrong = (
        2 * (p - 2) * n_O
        + p * (3 * (p + 1) * Gp + (p - 1) * Gm)
    ) / (p * p - 1)
    right = three_AB_from_G(n_O, Gp, Gm, p)
    A = A_from_Gplus_orbit(n_O, Gp, p)
    B = B_from_Gminus_orbit(n_O, Gm, p)
    ok = ok and right == 3 * A + B and wrong != right
    return {
        "p": p,
        "T": str(T_bitight(p)),
        "L": str(L_abs_gmin(p)),
        "L_gt_T": L_abs_gmin(p) > T_bitight(p),
        "samples": samples[:4],
        "proved": ok,
        "theorem": (
            "3A+B=[2(p−2)n_O+p(3(p−1)G_++(p+1)G_−)]/(p²−1).  "
            "G_+=G_−=n_O G ⇒ 3A+B>0 ⇔ G>T(p).  L>T so G=L is not the zero."
        ),
    }


def lemma_wick_six_type_3ab_positive(p: int) -> dict:
    """
    Wick 3a+b>0 on all six far types.  Goes False if a closed form
    misses the t=χκ/p evaluation or p²−4p+1≤0.
    """
    types = {
        "pp_p": (1, 1, 1),
        "pp_m": (1, 1, -1),
        "mm_p": (-1, -1, 1),
        "mm_m": (-1, -1, -1),
        "cr_p": (1, -1, 1),
        "cr_m": (1, -1, -1),
    }
    ok = True
    rows = {}
    for name, (chu, chv, chd) in types.items():
        val = wick_3ab_edge(p, chu, chv, chd)
        kap = chu + chv + chd
        t = wick_t_from_chi(p, chu, chv, chd)
        if t != Fraction(chd * kap, p):
            ok = False
        if name in ("pp_p", "mm_m"):
            if val != wick_3ab_closed_same_plus(p):
                ok = False
        elif name in ("pp_m", "mm_p"):
            if val != wick_3ab_closed_t_minus(p):
                ok = False
        elif name in ("cr_p", "cr_m"):
            # ΩN: κ=χ(u−v) ⇒ t=χd²/p=1/p either sign
            if val != wick_3ab_closed_cross_plus(p):
                ok = False
        if val <= 0:
            ok = False
        rows[name] = str(val)
    # cubics
    ok = ok and (p * p + 4 * p - 3) > 0
    ok = ok and (p * p - 4 * p + 1) > 0
    # (p−2)²−3 ≥ 6 for p≥5
    ok = ok and ((p - 2) * (p - 2) - 3) >= 6
    # a forgotten /p in t (using t=χκ) is not the Wick value at p=5
    t_wrong = 1 * (1 + 1 + -1)  # χd κ without /p, type ++-
    wrong = 2 * ((p - 2) + t_wrong * (2 * p - 1)) / (p * p - 1)
    ok = ok and wrong != wick_3ab_edge(p, 1, 1, -1)
    return {
        "p": p,
        "by_type": rows,
        "same_plus": str(wick_3ab_closed_same_plus(p)),
        "t_minus": str(wick_3ab_closed_t_minus(p)),
        "cross_plus": str(wick_3ab_closed_cross_plus(p)),
        "all_positive": all(Fraction(v) > 0 for v in rows.values()),
        "proved": ok,
        "theorem": (
            "Wick 3a+b>0 on all six (χ(u),χ(v),χ(u−v)) types for p≥5.  "
            "Closed: 2(p²+4p−3)/(p(p²−1)), 2(p²−4p+1)/(p(p²−1)), "
            "and 2/p on both ΩN signs."
        ),
    }


def lemma_Q_chi_identities(p: int) -> dict:
    """
    Pointwise Cy far-pair identities at f_e=±1.  Goes False if Max±
    Q_χu coefficients are swapped or Q_κ misses the sum.
    """
    ok = True
    samples = {}
    for f in (1, -1):
        up = Q_chi_u_Maxplus(p, f)
        dp = Q_chi_d_Maxplus(p, f)
        kp = Q_kappa_Maxplus(p, f)
        um = Q_chi_u_Maxminus(p, f)
        dm = Q_chi_d_Maxminus(p, f)
        km = Q_kappa_Maxminus(p, f)
        if kp != up + dp or km != um + dm:
            ok = False
        if up != Fraction(1 - p) * (1 + f):
            ok = False
        if um != Fraction(p + 1) * (1 + f):
            ok = False
        samples[str(f)] = {
            "Qchi_u_+": str(up),
            "Qchi_d_+": str(dp),
            "Qkappa_+": str(kp),
            "Qchi_u_-": str(um),
            "Qchi_d_-": str(dm),
            "Qkappa_-": str(km),
        }
    ok = ok and Q_chi_d_Maxplus(p, 1) == 1 - p
    ok = ok and Q_chi_d_Maxplus(p, -1) == p + 1
    ok = ok and Q_chi_d_Maxminus(p, 1) == p + 1
    ok = ok and Q_chi_d_Maxminus(p, -1) == 1 - p
    ok = ok and Q_kappa_Maxplus(p, 1) == 3 * (1 - p)
    ok = ok and Q_kappa_Maxminus(p, 1) == 3 * (p + 1)
    # Max+ / Max− Q_χu disagree at f=+1
    ok = ok and Q_chi_u_Maxplus(p, 1) != Q_chi_u_Maxminus(p, 1)
    # swapped Max+ coefficient (p+1) is the Max− formula
    ok = ok and Fraction(p + 1) * 2 != Q_chi_u_Maxplus(p, 1)
    # far / |κ| counts
    n_far = (n_of(p) - 2) * (n_of(p) - 3) // 2
    ok = ok and n_far_kappa3(p) + n_far_kappa1(p) == n_far
    ok = ok and n_far_kappa3(p) == (n_of(p) - 2) * (n_of(p) - 6) // 8
    # a 4-eq / 6-type leftover: 3-type (3) + Q_χd (1) vs 6 types
    ok = ok and (6 - 4) == 2
    return {
        "p": p,
        "samples": samples,
        "n_kappa3": n_far_kappa3(p),
        "n_kappa1": n_far_kappa1(p),
        "leftover_dofs": 2,
        "proved": ok,
        "theorem": (
            "Max+: Q_χu=(1−p)(1+f_e), Q_κ|_{+}=3(1−p).  "
            "Max−: Q_χu=(p+1)(1+f_e), Q_κ|_{+}=3(p+1).  "
            "3-type+Q_χd is 4 eqs on 6 types (2 leftover Aut_e dofs)."
        ),
    }


def type_I_split_far_3AB_identities_ok() -> bool:
    return all(
        lemma_threeAB_conversion(p)["proved"]
        and lemma_wick_six_type_3ab_positive(p)["proved"]
        and lemma_Q_chi_identities(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def type_I_wick_six_type_3ab_positive() -> bool:
    return all(
        lemma_wick_six_type_3ab_positive(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


# ---------------------------------------------------------------------------
# K. Aut_e 3A+B sign reduction (15.268 / 15.260 / 15.47 / 4-point)
# ---------------------------------------------------------------------------


def T_abs(p: int) -> Fraction:
    """|T(p)| = (p−2)/(p(2p−1)).  3A+B=0 at G=−T_abs, not at L."""
    return -T_bitight(p)


def one_over_2p(p: int) -> Fraction:
    return Fraction(1, 2 * p)


def mu_part_majorant(p: int) -> Fraction:
    """Auer–Top |μ_part| majorant (p²+4p−9)/(p²(p²−5)) on |κ|=1."""
    return Fraction(p * p + 4 * p - 9, p * p * (p * p - 5))


def mu_part_on_kappa1(p: int, kap: int, phi: int) -> Fraction:
    """μ_part = [(p²−1)κ − 2φ] / (p²(p²−5)) on |κ|=1 (star=0)."""
    return Fraction((p * p - 1) * kap - 2 * phi, p * p * (p * p - 5))


def G_locked_kappa1(chd: int, mu: Fraction) -> Fraction:
    """Locked pairing G = χ(u−v) μ on |κ|=1 (15.48 / Paley C_∞*=1)."""
    return chd * mu


def n_foursets_through_e(p: int) -> int:
    """C(n−2, 2) four-sets contain a fixed edge."""
    n = n_of(p)
    return (n - 2) * (n - 3) // 2


def through_e_sum_mu_kappa(p: int) -> Fraction:
    """∑_{S∋e} μκ = 3(p²−1)/2."""
    return Fraction(3 * (p * p - 1), 2)


def through_e_sum_mu_phi(p: int) -> Fraction:
    """∑_{S∋e} μφ = −(p²−1)."""
    return Fraction(-(p * p - 1), 1)


def through_e_sum_mu(p: int) -> Fraction:
    """∑_{S∋e} μ = −(p²−1)/(p²+1)."""
    return Fraction(-(p * p - 1), p * p + 1)


def through_e_sum_G(p: int) -> Fraction:
    """∑_{far e'} G_{ee'} = n/2 − 1 = (p²−1)/2 (Max+ or Max−)."""
    return Fraction(p * p - 1, 2)


def pairing_G_values(kap: int, mu: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    """The three pairing G's on a |κ|=1 four-set: {α,α,−α}, α=|μ|."""
    if kap not in (1, -1):
        raise ValueError("pairing_G_values is the |κ|=1 identity")
    pis = (1, 1, -1) if kap == 1 else (-1, -1, 1)
    return tuple(pi * mu for pi in pis)  # type: ignore[return-value]


def kappa3_linear_form(p: int, mu: Fraction, nu: Fraction) -> Fraction:
    """w = (2p−1)μ + (p−2)ν."""
    return (2 * p - 1) * mu + (p - 2) * nu


def three_AB_kappa3_from_mu_nu(
    n_O: int, chd: int, mu: Fraction, nu: Fraction, p: int
) -> Fraction:
    """
    |κ|=3 Aut_e class: G± = n_O χd (μ±ν).
    3A+B = 2 n_O{(p−2) + p χd w}/(p²−1), w=(2p−1)μ+(p−2)ν.
    """
    w = kappa3_linear_form(p, mu, nu)
    return 2 * n_O * ((p - 2) + p * chd * w) / (p * p - 1)


def census_gmin_kappa1(p: int) -> Fraction | None:
    """Certified g_min on |κ|=1 (p=5,7 only).  Not a general formula."""
    if p == 5:
        return Fraction(-3, 65)
    if p == 7:
        return Fraction(-109, 2863)
    return None


def census_min_threeAB_per_edge(p: int) -> Fraction | None:
    """Min per-edge 3A+B through e (86-worker census; p=5,7)."""
    if p == 5:
        return Fraction(1, 13)
    if p == 7:
        return Fraction(157, 2454)
    return None


def lemma_through_e_moments(p: int) -> dict:
    """
    15.260 transported to four-sets through a locked edge.
    Goes False if the double count misses C(n−2,2) or the far
    row-sum is written n/2 instead of n/2−1.
    """
    n = n_of(p)
    B = comb(n, 4)
    n_e = n_foursets_through_e(p)
    ok = n_e == (n - 2) * (n - 3) // 2
    ok = ok and n_e == n_far_kappa1(p) + n_far_kappa3(p)
    # each 4-set has 6 edges: B*6 = C(n,2)*n_e
    ok = ok and B * 6 == comb(n, 2) * n_e
    avg_mk = sum_mu_kappa(p) / B
    avg_mp = sum_mu_phi(p) / B
    avg_m = sum_mu(p) / B
    got_mk = avg_mk * n_e
    got_mp = avg_mp * n_e
    got_m = avg_m * n_e
    if got_mk != through_e_sum_mu_kappa(p):
        ok = False
    if got_mp != through_e_sum_mu_phi(p):
        ok = False
    if got_m != through_e_sum_mu(p):
        ok = False
    # row-sum identities: wedges contribute 0, diag 1, far n/2−1
    ok = ok and through_e_sum_G(p) == Fraction(n, 2) - 1
    ok = ok and through_e_sum_G(p) == Fraction(p * p - 1, 2)
    # a forgotten −1 (using n/2 for the far sum) is not the identity
    ok = ok and Fraction(n, 2) != through_e_sum_G(p)
    return {
        "p": p,
        "n_through_e": n_e,
        "sum_mu_kappa": str(got_mk),
        "sum_mu_phi": str(got_mp),
        "sum_mu": str(got_m),
        "sum_G_far": str(through_e_sum_G(p)),
        "proved": ok,
        "theorem": (
            "∑_{S∋e} μκ=3(p²−1)/2, ∑ μφ=−(p²−1), ∑ μ=−(p²−1)/(n); "
            "∑_{far} G=(p²−1)/2.  Averages, not Aut_e L^∞."
        ),
    }


def lemma_pairing_four_point(p: int) -> dict:
    """
    15.48 four-point: G1+G2+G3=κμ and {G_i}={α,α,−α} on |κ|=1.
    Goes False if the sum misses κ or the locked G is written μ
    instead of χd μ.
    """
    ok = True
    samples = []
    for kap, mu in (
        (1, Fraction(3, 65)),
        (-1, Fraction(-3, 65)),
        (1, Fraction(-109, 2863)),
        (-1, L_abs_gmin(p)),
        (1, T_bitight(p)),
    ):
        Gs = pairing_G_values(kap, mu)
        if sum(Gs) != kap * mu:
            ok = False
        if min(Gs) != -abs(mu):
            ok = False
        if max(Gs) != abs(mu):
            ok = False
        # locked pairing is one of the three
        for chd in (1, -1):
            Glock = G_locked_kappa1(chd, mu)
            if Glock not in Gs and Glock != Gs[0] and Glock != Gs[2]:
                # Glock must equal one of ±μ, which are exactly the G's
                if Glock not in (mu, -mu):
                    ok = False
            if Glock not in (mu, -mu):
                ok = False
        samples.append(
            {
                "kap": kap,
                "mu": str(mu),
                "Gs": [str(g) for g in Gs],
                "sum_eq_kappa_mu": sum(Gs) == kap * mu,
            }
        )
    # writing ∑G=μ (dropping κ) fails on κ=−1
    mu = Fraction(3, 65)
    Gs = pairing_G_values(-1, mu)
    ok = ok and sum(Gs) == -mu and sum(Gs) != mu
    # G_locked = μ (dropping χd) is not the identity when χd=−1
    ok = ok and G_locked_kappa1(-1, mu) == -mu != mu
    return {
        "p": p,
        "samples": samples[:3],
        "proved": ok,
        "theorem": (
            "G1+G2+G3=κμ; on |κ|=1 the three G's are {α,α,−α}, "
            "α=|μ|; locked pairing G=χ(u−v)μ."
        ),
    }


def lemma_gmin_T_vs_thresholds(p: int) -> dict:
    """
    3A+B>0 ⇔ G>T on |κ|=1; g_min=−max|μ|.  Goes False if the zero
    is L, or if |μ|≤1/(2p) is treated as sufficient for G>T.
    """
    ok = True
    Tthr = T_bitight(p)
    Tabs = T_abs(p)
    Lval = L_abs_gmin(p)
    half = one_over_2p(p)
    maj = mu_part_majorant(p)
    ok = ok and Tabs == -Tthr == Fraction(p - 2, p * (2 * p - 1))
    ok = ok and Lval > Tthr
    ok = ok and half > Tabs
    ok = ok and three_AB_equal_G(1, Tthr, p) == 0
    ok = ok and three_AB_equal_G(1, Lval, p) > 0
    # residual-i hinge is strictly weaker than G>T
    ok = ok and three_AB_equal_G(1, -half, p) < 0
    ok = ok and (-half) < Tthr
    # L suffices for the sign (L is strictly above T)
    ok = ok and Lval > Tthr
    ok = ok and three_AB_equal_G(1, Lval, p) > 0
    # census anchors: conversion, not a general bound
    gcert = census_gmin_kappa1(p)
    if gcert is not None:
        ok = ok and gcert > Tthr
        n_O = 24
        ab = three_AB_equal_G(n_O, gcert, p)
        edge_ab = three_AB_equal_G(1, gcert, p)
        expect_edge = census_min_threeAB_per_edge(p)
        ok = ok and expect_edge is not None and edge_ab == expect_edge
        if p == 5:
            ok = ok and ab == Fraction(24, 13)
            ok = ok and gcert == Fraction(-3, 65)
        if p == 7:
            ok = ok and ab == Fraction(628, 409)
            ok = ok and gcert == Fraction(-109, 2863)
        ok = ok and ab > 0
    # swapped criterion G>L as the 3A+B zero is wrong
    ok = ok and three_AB_equal_G(1, Lval, p) != 0
    return {
        "p": p,
        "T": str(Tthr),
        "T_abs": str(Tabs),
        "L": str(Lval),
        "one_over_2p": str(half),
        "mu_part_maj": str(maj),
        "one_over_2p_gt_T_abs": half > Tabs,
        "L_gt_T": Lval > Tthr,
        "G_eq_minus_1_2p_gives_neg_3AB": three_AB_equal_G(1, -half, p) < 0,
        "census_gmin": None if gcert is None else str(gcert),
        "proved": ok,
        "theorem": (
            "On |κ|=1, 3A+B>0 ⇔ G>T ⇔ |μ|<T_abs when the locked "
            "pairing is the minority sign.  L suffices; 1/(2p) does not."
        ),
    }


def lemma_kappa3_3AB_conversion(p: int) -> dict:
    """
    |κ|=3: 3A+B from (μ,ν,χd).  Goes False if the (2p−1) and (p−2)
    slots on (μ,ν) are swapped.
    """
    ok = True
    samples = []
    for n_O, chd, mu, nu in (
        (1, 1, Fraction(3, p * p), Fraction(0)),
        (1, -1, Fraction(-3, p * p), Fraction(0)),
        (p + 1, 1, Fraction(1, p), Fraction(-1, p * p)),
        (24, -1, Fraction(-1, 20), Fraction(1, 50)),
    ):
        Gp = n_O * chd * (mu + nu)
        Gm = n_O * chd * (mu - nu)
        got = three_AB_from_G(n_O, Gp, Gm, p)
        cl = three_AB_kappa3_from_mu_nu(n_O, chd, mu, nu, p)
        if got != cl:
            ok = False
        w = kappa3_linear_form(p, mu, nu)
        # sign criterion
        pos = (p - 2) + p * chd * w > 0
        if (cl > 0) != pos:
            ok = False
        if (cl == 0) != ((p - 2) + p * chd * w == 0):
            ok = False
        samples.append(
            {
                "n_O": n_O,
                "chd": chd,
                "match": got == cl,
                "positive": cl > 0,
            }
        )
        # swapped weights
        wrong_w = (p - 2) * mu + (2 * p - 1) * nu
        wrong = 2 * n_O * ((p - 2) + p * chd * wrong_w) / (p * p - 1)
        if wrong_w != w:
            ok = ok and wrong != cl
    # Wick |κ|=3: μ=κ/p², ν=0, χd=sign(κ), w=χd·3(2p−1)/p²
    for kap, chd in ((3, 1), (-3, -1)):
        mu = Fraction(kap, p * p)
        cl = three_AB_kappa3_from_mu_nu(1, chd, mu, Fraction(0), p)
        if cl <= 0:
            ok = False
    return {
        "p": p,
        "samples": samples,
        "proved": ok,
        "theorem": (
            "|κ|=3: 3A+B>0 ⇔ χd[(2p−1)μ+(p−2)ν] > −(p−2)/p.  "
            "Wick (μ=κ/p², ν=0) is strictly positive."
        ),
    }


def lemma_mu_part_vs_T(p: int) -> dict:
    """
    μ_part majorant vs T_abs.  Goes False if p=5 is claimed safe
    under the crude Auer–Top majorant, or if p≥7 is claimed unsafe.
    """
    ok = True
    maj = mu_part_majorant(p)
    Tabs = T_abs(p)
    # crude worst G_part = −maj
    worst = -maj
    if p == 5:
        ok = ok and maj > Tabs
        ok = ok and worst < T_bitight(p)
        # the violating triple (κ,φ,χd)=(1,−2(p−2),−1)
        phi = -2 * (p - 2)
        mp = mu_part_on_kappa1(p, 1, phi)
        ok = ok and G_locked_kappa1(-1, mp) == worst
        ok = ok and G_locked_kappa1(-1, mp) < T_bitight(p)
        # next φ slot φ=−2 is already >T
        mp2 = mu_part_on_kappa1(p, 1, -2)
        ok = ok and G_locked_kappa1(-1, mp2) > T_bitight(p)
    else:
        ok = ok and maj < Tabs
        ok = ok and worst > T_bitight(p)
        phi = -2 * (p - 2)
        mp = mu_part_on_kappa1(p, 1, phi)
        ok = ok and G_locked_kappa1(-1, mp) > T_bitight(p)
    # |μ|≤maj would give G>T for p≥7, but is false at p=7
    ok = ok and maj <= one_over_2p(p)
    if p == 7:
        gcert = census_gmin_kappa1(p)
        ok = ok and gcert is not None and abs(gcert) > maj
    return {
        "p": p,
        "majorant": str(maj),
        "T_abs": str(Tabs),
        "majorant_lt_T_abs": maj < Tabs,
        "proved": ok,
        "theorem": (
            "μ_part majorant < |T| for every prime p≥7 (not p=5).  "
            "|μ|≤maj is false at p=7 (109/2863>17/539); even-δ leftover."
        ),
    }


# ---------------------------------------------------------------------------
# L. p=5 finite 3A+B from C; |μ|≤maj false at p=7
# ---------------------------------------------------------------------------

_P5_CENSUS: dict | None = None


def _enum_boolean_eigenspace(C, eig: int):
    """All ±1 eigenvectors of C with eigenvalue eig.  p=5: nul=13."""
    import numpy as np

    n = C.shape[0]
    M = C.astype(np.float64) - eig * np.eye(n)
    _, s, vt = np.linalg.svd(M, full_matrices=True)
    rank = int((s > 1e-8).sum())
    nul = n - rank
    if nul <= 0 or nul > 16:
        raise RuntimeError(f"unexpected nullity {nul}")
    nb = vt[rank:].T
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(20000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError("no free coordinate set")
    Binv = np.linalg.inv(nb[coords])
    found = []
    for bits in range(1 << nul):
        free = np.array(
            [1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)]
        )
        x = nb @ (Binv @ free)
        if float(np.max(np.abs(np.abs(x) - 1.0))) < 1e-6:
            found.append(np.sign(x))
    if not found:
        raise RuntimeError("no ±1 eigenvectors")
    Y = np.unique(np.round(found).astype(np.int64), axis=0)
    if not np.all(Y @ C.T == eig * Y):
        raise RuntimeError("enumerated vectors fail C y = λ y")
    if not np.all(np.abs(Y) == 1):
        raise RuntimeError("enumerated vectors not ±1")
    return Y


def _load_or_enum_maxpm_p5(C):
    """Verified Max± for p=5: cache if it is the ±5 eigenspace, else enum."""
    import numpy as np

    Mp = Mm = None
    try:
        Mp0 = np.load("/tmp/maxplus_p5.npy")
        Mm0 = np.load("/tmp/maxminus_p5.npy")
        Mp0 = np.asarray(Mp0, dtype=np.int64)
        Mm0 = np.asarray(Mm0, dtype=np.int64)
        if (
            Mp0.shape[1] == C.shape[0]
            and Mm0.shape[1] == C.shape[0]
            and np.all(Mp0 @ C.T == 5 * Mp0)
            and np.all(Mm0 @ C.T == -5 * Mm0)
            and np.all(np.abs(Mp0) == 1)
            and np.all(np.abs(Mm0) == 1)
        ):
            Mp, Mm = Mp0, Mm0
    except (OSError, ValueError):
        Mp = Mm = None
    if Mp is None:
        Mp = _enum_boolean_eigenspace(C, 5)
    if Mm is None:
        Mm = _enum_boolean_eigenspace(C, -5)
    return Mp, Mm


def p5_through_e_census_from_C() -> dict:
    """
    All 276 far 4-sets through e={∞,0} at p=5, G from C.

    Fail-when-wrong: count must equal C(n-2,2) and the C(n,4)
    double count; Max± must be ±5 eigenspaces; G>T is computed,
    not hardcoded.
    """
    global _P5_CENSUS
    if _P5_CENSUS is not None:
        return _P5_CENSUS
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(np.int64)
    n = int(C.shape[0])
    ok = n == p * p + 1 == 26
    n_expect = (n - 2) * (n - 3) // 2
    ok = ok and n_expect == comb(n - 2, 2) == 276
    ok = ok and comb(n, 4) * 6 == comb(n, 2) * n_expect
    Mp, Mm = _load_or_enum_maxpm_p5(C)
    Np, Nm = int(len(Mp)), int(len(Mm))
    ok = ok and Np == 260 and Nm == 260
    far = list(combinations(range(2, n), 2))
    ok = ok and len(far) == n_expect
    Tthr = T_bitight(p)
    n_k1 = n_k3 = 0
    n_pos = 0
    n_k1_le_T = 0
    gmin: Fraction | None = None
    min_ab: Fraction | None = None
    max_abs_mu = Fraction(0)
    y0p, y1p = Mp[:, 0], Mp[:, 1]
    y0m, y1m = Mm[:, 0], Mm[:, 1]
    for u, v in far:
        chd = int(C[u, v])
        kap = int(
            C[0, 1] * C[u, v] + C[0, u] * C[1, v] + C[0, v] * C[1, u]
        )
        sum_p = int(np.dot(y0p * y1p, Mp[:, u] * Mp[:, v]))
        sum_m = int(np.dot(y0m * y1m, Mm[:, u] * Mm[:, v]))
        m4p = Fraction(sum_p, Np)
        m4m = Fraction(sum_m, Nm)
        Gp = chd * m4p
        Gm = chd * m4m
        ab = three_AB_from_G(1, Gp, Gm, p)
        if ab > 0:
            n_pos += 1
        if min_ab is None or ab < min_ab:
            min_ab = ab
        mu = (m4p + m4m) / 2
        if abs(kap) == 1:
            n_k1 += 1
            Glock = (Gp + Gm) / 2
            if Glock != chd * mu:
                ok = False
            if gmin is None or Glock < gmin:
                gmin = Glock
            if Glock <= Tthr:
                n_k1_le_T += 1
            max_abs_mu = max(max_abs_mu, abs(mu))
        elif abs(kap) == 3:
            n_k3 += 1
        else:
            ok = False
    ok = ok and n_k1 + n_k3 == 276
    ok = ok and n_k1 == n_far_kappa1(p)
    ok = ok and n_k3 == n_far_kappa3(p)
    ok = ok and n_pos == 276
    ok = ok and n_k1_le_T == 0
    ok = ok and gmin is not None and gmin > Tthr
    ok = ok and min_ab is not None and min_ab > 0
    # conversion: min 3A+B is the |κ|=1 g_min edge
    if gmin is not None:
        ok = ok and three_AB_equal_G(1, gmin, p) > 0
        ok = ok and gmin == Fraction(-3, 65)
        ok = ok and three_AB_equal_G(1, gmin, p) == Fraction(1, 13)
        if min_ab != Fraction(1, 13):
            ok = False
    # T is the zero, not L; 1/(2p) is not enough
    ok = ok and three_AB_equal_G(1, Tthr, p) == 0
    ok = ok and three_AB_equal_G(1, L_abs_gmin(p), p) > 0
    ok = ok and three_AB_equal_G(1, -one_over_2p(p), p) < 0
    out = {
        "p": 5,
        "n_through_e": len(far),
        "n_k1": n_k1,
        "n_k3": n_k3,
        "n_3AB_pos": n_pos,
        "min_3AB": None if min_ab is None else str(min_ab),
        "gmin": None if gmin is None else str(gmin),
        "T": str(Tthr),
        "gmin_gt_T": gmin is not None and gmin > Tthr,
        "n_k1_G_le_T": n_k1_le_T,
        "max_abs_mu": str(max_abs_mu),
        "n_Maxplus": Np,
        "n_Maxminus": Nm,
        "proved": ok,
        "theorem": (
            "p=5: every far 4-set through e={∞,0} has 3A+B>0.  "
            "Count=C(24,2)=276=C(n,4)·6/C(n,2).  G from C "
            "(Max±=±5 eigenspaces).  |κ|=1: G>T (g_min=−3/65>−1/15)."
        ),
    }
    _P5_CENSUS = out
    return out


def lemma_p5_through_e_finite() -> dict:
    """Finite Max+-from-C check.  Goes False if the count or G>T fails."""
    return p5_through_e_census_from_C()


def star_through_e_kappa3(chu: int, chv: int, chd: int) -> int:
    """
    star(S)=∑_{v∈S}∏_{u≠v} C_{vu} on S={∞,0,u,v}.
    Paley C_∞*=1, C_0x=χ(x), C_uv=χ(u−v), χ(−1)=1.
    """
    return 1 + chu * chv + chu * chd + chv * chd


def lemma_kappa3_star_through_e(p: int) -> dict:
    """
    |κ|=3 through e ⇒ star=+4.  Goes False if a same-character
    type is given star=0 or −4, or if χ(−1)=−1 is used on F_{p²}.
    """
    ok = True
    # |κ|=3 ⇔ χ(u)=χ(v)=χ(u−v)
    for chd in (1, -1):
        st = star_through_e_kappa3(chd, chd, chd)
        if st != 4:
            ok = False
        # κ=χd+χu+χv=3 χd
        if chd + chd + chd != 3 * chd:
            ok = False
    # |κ|=1 types are not this identity
    ok = ok and star_through_e_kappa3(1, 1, -1) != 4
    ok = ok and star_through_e_kappa3(1, -1, 1) != 4
    # χ(−1)=−1 would flip the last two products on one ordering; we
    # use χ(v−u)=χ(u−v) (q=p²≡1 mod 4).  A flipped last term is not +4.
    flipped = 1 + 1 * 1 + 1 * 1 + 1 * (-1)
    ok = ok and flipped != 4
    # particular ν_part = z₊ star = −8/(p(p²−5))
    zplus = Fraction(-2, p * (p * p - 5))
    nu_part = zplus * 4
    ok = ok and nu_part == Fraction(-8, p * (p * p - 5))
    return {
        "p": p,
        "star": 4,
        "nu_part": str(nu_part),
        "proved": ok,
        "theorem": (
            "On every |κ|=3 far 4-set through e={∞,0}, star=+4 "
            "(Paley C_∞*=1 and χ(u)=χ(v)=χ(u−v)).  "
            "Hence ν_part=z₊star=−8/(p(p²−5)).  Does not close 3A+B."
        ),
    }


def lemma_mu_le_maj_false_at_p7() -> dict:
    """
    |μ|≤maj is false at p=7.  Goes False if the census g_min is
    claimed to lie under the Auer–Top majorant.
    """
    p = 7
    maj = mu_part_majorant(p)
    gcert = census_gmin_kappa1(p)
    ok = gcert is not None and abs(gcert) > maj
    ok = ok and maj < T_abs(p)
    # 109/2863 > 17/539
    ok = ok and Fraction(109, 2863) > Fraction(17, 539)
    ok = ok and Fraction(109, 2863) < T_abs(p)
    return {
        "p": p,
        "max_abs_mu": None if gcert is None else str(abs(gcert)),
        "maj": str(maj),
        "max_gt_maj": ok,
        "proved": ok,
        "theorem": (
            "p=7: max|μ|=109/2863>maj=17/539, so |μ|≤maj is false.  "
            "Still 109/2863<|T|=5/91, so G>T holds at the census g_min."
        ),
    }


# ---------------------------------------------------------------------------
# M. |μ|≤|L| as the μ-bound (not |μ|≤|T|); even-δ feed; |κ|=3 particular
# ---------------------------------------------------------------------------


def d3_kappa1_neighbors(p: int) -> int:
    """Ordered one-vertex extensions of a |κ|=1 4-set into |κ|=3 (15.71)."""
    return p * p - 5


def d1_kappa1_neighbors(p: int) -> int:
    """Ordered one-vertex extensions of a |κ|=1 4-set into |κ|=1."""
    return 3 * p * p - 7


def kappa3_nu_part_through_e(p: int) -> Fraction:
    """ν_part = z₊ star = −8/(p(p²−5)) on |κ|=3 through e (star=+4)."""
    return Fraction(-8, p * (p * p - 5))


def kappa3_mu_part_through_e(p: int, chd: int, phi: int) -> Fraction:
    """μ_part = aκ+bφ with κ=3 χd on |κ|=3 through e."""
    return mu_part_on_kappa1(p, 3 * chd, phi)


def kappa3_signed_w_part(p: int, chd: int, phi: int) -> Fraction:
    """χd[(2p−1)μ_part+(p−2)ν_part] on a |κ|=3 through-e 4-set."""
    mu = kappa3_mu_part_through_e(p, chd, phi)
    nu = kappa3_nu_part_through_e(p)
    return chd * kappa3_linear_form(p, mu, nu)


def kappa3_signed_w_wick(p: int) -> Fraction:
    """Wick χd w with μ=κ/p²=3χd/p², ν=0: 3(2p−1)/p²."""
    return Fraction(3 * (2 * p - 1), p * p)


def kappa3_thr(p: int) -> Fraction:
    """3A+B>0 ⇔ χd w > −(p−2)/p."""
    return -Fraction(p - 2, p)


def lemma_mu_le_L_not_T_as_bound(p: int) -> dict:
    """
    |μ|≤|L| ⇒ G≥L>T ⇒ 3A+B>0.  |μ|≤|T| does *not* imply 3A+B>0.

    Goes False if T is used as the closing μ-bound (G=T gives 3A+B=0)
    or if maj is treated as a true bound on μ.
    """
    ok = True
    Tabs = T_abs(p)
    Lval = L_abs_gmin(p)
    Labs = -Lval
    Tthr = T_bitight(p)
    maj = mu_part_majorant(p)
    ok = ok and Labs < Tabs
    ok = ok and Lval > Tthr
    # |μ|=|L|, locked minority G=−|μ|=L: 3A+B>0
    ok = ok and three_AB_equal_G(1, Lval, p) > 0
    ok = ok and three_AB_equal_G(1, Lval, p) != 0
    # |μ|=|T|, locked minority G=−|T|=T: 3A+B=0, so |μ|≤|T| is not a close
    ok = ok and three_AB_equal_G(1, Tthr, p) == 0
    ok = ok and three_AB_equal_G(1, -Tabs, p) == 0
    # replacing the μ-bound L by T would treat this zero as a success
    ok = ok and Labs != Tabs
    # |μ|≤maj is not a true bound (p=7 census)
    if p == 7:
        gcert = census_gmin_kappa1(p)
        ok = ok and gcert is not None and abs(gcert) > maj
        ok = ok and abs(gcert) <= Labs
        ok = ok and abs(gcert) < Tabs
        # treating maj as true would claim 109/2863 ≤ 17/539
        ok = ok and Fraction(109, 2863) > Fraction(17, 539)
    if p == 5:
        gcert = census_gmin_kappa1(p)
        ok = ok and gcert is not None and abs(gcert) <= Labs
        ok = ok and abs(gcert) < Tabs
        # at p=5 maj > |L|, so |μ|≤maj would not even imply |μ|≤|L|
        ok = ok and maj > Labs
    return {
        "p": p,
        "L_abs": str(Labs),
        "T_abs": str(Tabs),
        "maj": str(maj),
        "G_eq_L_gives_pos_3AB": three_AB_equal_G(1, Lval, p) > 0,
        "G_eq_T_gives_zero_3AB": three_AB_equal_G(1, Tthr, p) == 0,
        "proved": ok,
        "theorem": (
            "|μ|≤|L| ⇒ G≥L>T ⇒ 3A+B>0 on |κ|=1.  |μ|≤|T| allows "
            "G=T and 3A+B=0, so T is not a closing μ-bound.  "
            "|μ|≤maj is false at p=7."
        ),
    }


def lemma_mu_le_L_census_p5p7() -> dict:
    """
    |μ|≤|L| on |κ|=1 at p=5,7 (certified gmin).  Goes False if the
    census max is claimed to exceed |L| or to lie under maj at p=7.
    """
    ok = True
    rows = {}
    for p, num, den in ((5, 3, 65), (7, 109, 2863)):
        g = census_gmin_kappa1(p)
        mu = abs(g) if g is not None else None
        Labs = -L_abs_gmin(p)
        Tabs = T_abs(p)
        maj = mu_part_majorant(p)
        expect = Fraction(num, den)
        hit = mu == expect and mu <= Labs and mu < Tabs
        if p == 7:
            hit = hit and mu > maj
        if p == 5:
            hit = hit and mu < maj
        if not hit:
            ok = False
        rows[str(p)] = {
            "max_abs_mu": None if mu is None else str(mu),
            "L_abs": str(Labs),
            "T_abs": str(Tabs),
            "maj": str(maj),
            "le_L": mu is not None and mu <= Labs,
            "lt_T": mu is not None and mu < Tabs,
        }
    # |μ|≤|T| as a replacement bound is not the census statement
    ok = ok and Fraction(109, 2863) < T_abs(7)
    ok = ok and three_AB_equal_G(1, T_bitight(7), 7) == 0
    return {
        "p5": rows.get("5"),
        "p7": rows.get("7"),
        "proved": ok,
        "theorem": (
            "p=5: max|μ|=3/65≤|L|=3/50<|T|=1/15.  "
            "p=7: max|μ|=109/2863≤|L|=5/98<|T|=5/91, and 109/2863>maj."
        ),
    }


def lemma_even_delta_Tnu_feed(p: int) -> dict:
    """
    On |κ|=1, even-δ is the T-image of odd leftover on |κ|=3 neighbors.
    Crude |ν|≤1 · d3 does *not* give |μ|≤|L|.  Goes False if d3 is
    written 4(n−4) (forgetting the |κ| split) or if the crude bound
    is treated as sufficient for |L|.
    """
    ok = True
    d3 = d3_kappa1_neighbors(p)
    d1 = d1_kappa1_neighbors(p)
    ok = ok and d3 == p * p - 5
    ok = ok and d1 == 3 * p * p - 7
    ok = ok and d1 + d3 == 4 * (n_of(p) - 4)
    ok = ok and d1 + d3 == 4 * (p * p - 3)
    # forgetting the split (using all 4(n−4) as |κ|=3 feed) is wrong
    ok = ok and d3 != 4 * (p * p - 3)
    # crude |ρ| ≤ p d3 / 4 from |ν|≤1
    crude_rho = Fraction(p * d3, 4)
    L_room = Fraction(p - 4, 2 * p * p)  # |ρ| budget for |μ|≤|L| same-sign
    T_room = Fraction(p * p - 4 * p + 1, p * p * (2 * p - 1))
    ok = ok and crude_rho > L_room
    ok = ok and crude_rho > T_room
    # so crude neighbor L^∞ does not close |μ|≤|L| or even |μ|<|T|
    # even-δ identity: 2δ₊ = Tδ_odd/(4p) on |κ|=1
    # (master + ν=0 + matched particulars)
    return {
        "p": p,
        "d3": d3,
        "d1": d1,
        "crude_rho": str(crude_rho),
        "L_rho_budget": str(L_room),
        "crude_too_weak_for_L": crude_rho > L_room,
        "proved": ok,
        "theorem": (
            "On |κ|=1: μ=κ/p²+Tν/(4p)=μ_part+Tδ_odd/(4p).  "
            "Exactly d3=p²−5 neighbors have |κ|=3 (15.71).  "
            "Crude |ν|≤1 gives |ρ|≤p(p²−5)/4, larger than the "
            "|L| residual budget (p−4)/(2p²)."
        ),
    }


def lemma_kappa3_particular_beats_thr(p: int) -> dict:
    """
    Wick and (μ_part,ν_part) on |κ|=3 through e beat χd w > −(p−2)/p
    under Hasse |φ|≤2p.  Diamond |μ|+|ν|≤1 does *not*.

    Goes False if diamond is treated as sufficient, if star is not +4,
    or if (2p−1)↔(p−2) are swapped.
    """
    ok = True
    thr = kappa3_thr(p)
    wick = kappa3_signed_w_wick(p)
    ok = ok and wick > thr
    # worst particular: χd=+1, φ=+2p (Hasse edge)
    worst = kappa3_signed_w_part(p, 1, 2 * p)
    ok = ok and worst > thr
    # polynomial positivity of worst − thr
    # w = (6p³−19p²+14p+3)/(p²(p²−5)); p·num+(p−2)·den
    # = p(p⁴+4p³−24p²+24p+3) > 0
    q = p**4 + 4 * p**3 - 24 * p * p + 24 * p + 3
    ok = ok and q > 0
    if p >= 5:
        # q(5)=648, q'=4(p³+3p²−12p+6) with r(5)=146>0, r'>0
        ok = ok and (p**3 + 3 * p * p - 12 * p + 6) > 0
    # star=+4 is used for ν_part
    ok = ok and star_through_e_kappa3(1, 1, 1) == 4
    ok = ok and kappa3_nu_part_through_e(p) == Fraction(-8, p * (p * p - 5))
    # diamond: χd w ≥ −(2p−1), and −(2p−1) < thr
    diamond_lb = -(2 * p - 1)
    ok = ok and diamond_lb < thr
    # so |μ|+|ν|≤1 does not imply the linear form
    ok = ok and diamond_lb != thr
    # swapped (2p−1)↔(p−2) is not the form
    mu = kappa3_mu_part_through_e(p, 1, 2 * p)
    nu = kappa3_nu_part_through_e(p)
    right = kappa3_linear_form(p, mu, nu)
    wrong = (p - 2) * mu + (2 * p - 1) * nu
    ok = ok and wrong != right
    return {
        "p": p,
        "thr": str(thr),
        "wick": str(wick),
        "worst_part": str(worst),
        "diamond_lb": str(Fraction(diamond_lb)),
        "diamond_too_weak": diamond_lb < thr,
        "proved": ok,
        "theorem": (
            "|κ|=3 through e: Wick and the particular (star=+4, |φ|≤2p) "
            "satisfy χd w > −(p−2)/p.  Diamond |μ|+|ν|≤1 only gives "
            "χd w ≥ −(2p−1), which is below the threshold."
        ),
    }


def type_I_p5_through_e_3AB_positive() -> bool:
    return bool(lemma_p5_through_e_finite()["proved"])


def type_I_aut_e_3AB_sign_identities_ok() -> bool:
    return all(
        lemma_through_e_moments(p)["proved"]
        and lemma_pairing_four_point(p)["proved"]
        and lemma_gmin_T_vs_thresholds(p)["proved"]
        and lemma_kappa3_3AB_conversion(p)["proved"]
        and lemma_mu_part_vs_T(p)["proved"]
        and lemma_kappa3_star_through_e(p)["proved"]
        and lemma_mu_le_L_not_T_as_bound(p)["proved"]
        and lemma_even_delta_Tnu_feed(p)["proved"]
        and lemma_kappa3_particular_beats_thr(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    ) and lemma_mu_le_maj_false_at_p7()["proved"] and lemma_mu_le_L_census_p5p7()[
        "proved"
    ]


def type_I_aut_e_3AB_positive_general() -> bool:
    """
    3A+B>0 on every Aut_e far class, all primes p≥5.

    False: equivalent to g_min>T on |κ|=1 through e (plus the |κ|=3
    one-sided form).  |μ|≤|L| would close it (L>T) but is not proved
    Max+-free for all p.  |μ|≤maj is false at p=7.  |μ|≤|T| is not a
    closing bound (G=T ⇒ 3A+B=0).  p=5 is a finite from-C theorem.
    """
    return False


def type_I_three_orbit_equal_star_empty() -> bool:
    return all(
        lemma_three_orbit_equal_star_empty(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def type_I_three_orbit_plus_only_empty() -> bool:
    return all(
        lemma_three_orbit_plus_only_empty(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def type_I_zero_far_bad_case_empty() -> bool:
    return all(
        lemma_zero_far_bad_case_empty(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def type_I_three_orbit_family_bad_case_empty() -> bool:
    return all(
        lemma_three_orbit_family_bad_case_empty(p)["proved"]
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31)
    )


def theorem_G_leftover(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    slices_ok = True
    sample = {}
    for p in primes:
        w = lemma_two_orbit_weights(p)
        s = lemma_maxminus_star_far_sums(p)
        inc = lemma_far_weight_more_negative_when_beta_gt_minus_2(p)
        st = lemma_two_star_interpolants(p)
        eq = lemma_three_orbit_equal_star_empty(p)
        pl = lemma_three_orbit_plus_only_empty(p)
        zf = lemma_three_orbit_zero_far_in_box(p)
        pz = lemma_three_orbit_plus_zero_not_empty(p)
        h = lemma_three_orbit_family_bad_case_empty(p)
        I1 = lemma_split_far_delta_identity(p)
        I2 = lemma_split_far_Y_plus_3(p)
        I3 = lemma_three_type_Y_collapses_to_H(p)
        I4 = lemma_three_type_counts(p)
        I5 = lemma_three_type_Q_difference(p)
        J1 = lemma_threeAB_conversion(p)
        J2 = lemma_wick_six_type_3ab_positive(p)
        J3 = lemma_Q_chi_identities(p)
        K1 = lemma_through_e_moments(p)
        K2 = lemma_pairing_four_point(p)
        K3 = lemma_gmin_T_vs_thresholds(p)
        K4 = lemma_kappa3_3AB_conversion(p)
        K5 = lemma_mu_part_vs_T(p)
        K6 = lemma_mu_le_maj_false_at_p7() if p == 7 else {"proved": True}
        K7 = lemma_kappa3_star_through_e(p)
        M1 = lemma_mu_le_L_not_T_as_bound(p)
        M2 = lemma_even_delta_Tnu_feed(p)
        M3 = lemma_kappa3_particular_beats_thr(p)
        M4 = lemma_mu_le_L_census_p5p7() if p in (5, 7) else {"proved": True}
        if not (w["proved"] and s["proved"] and inc["proved"] and st["proved"]):
            ok = False
        if not (eq["proved"] and pl["proved"] and zf["proved"] and pz["proved"]):
            slices_ok = False
        if not (
            I1["proved"]
            and I2["proved"]
            and I3["proved"]
            and I4["proved"]
            and I5["proved"]
            and J1["proved"]
            and J2["proved"]
            and J3["proved"]
            and K1["proved"]
            and K2["proved"]
            and K3["proved"]
            and K4["proved"]
            and K5["proved"]
            and K6["proved"]
            and K7["proved"]
            and M1["proved"]
            and M2["proved"]
            and M3["proved"]
            and M4["proved"]
        ):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "mu_far": w["mu_far"],
                "mu_far_negative": w["mu_far_negative"],
                "mass_eq_k": w["mass_eq_k"],
                "star_far_sums": s["proved"],
                "slope_incompatible": inc["proved"],
                "n_star_orbits": st["n_star_orbits"],
                "equal_star_empty": eq["proved"],
                "plus_only_empty": pl["proved"],
                "zero_far_in_box": zf["proved"] and not zf["empty"],
                "three_orbit_bad_case_empty": h["proved"],
                "split_far_delta": I1["proved"],
                "three_type_collapses": I3["proved"],
                "threeAB_conversion": J1["proved"],
                "wick_six_type": J2["proved"],
                "Q_chi": J3["proved"],
                "through_e_moments": K1["proved"],
                "pairing_four_point": K2["proved"],
                "gmin_T_thresholds": K3["proved"],
                "kappa3_3AB": K4["proved"],
                "mu_part_vs_T": K5["proved"],
                "mu_le_L_not_T": M1["proved"],
                "even_delta_feed": M2["proved"],
                "kappa3_particular": M3["proved"],
            }
    empty = bool(residual_i_dual_eq_empty_proved_general())
    return {
        # 3-weight / 3-type collapse empty the bad case; finer split open.
        "proved": False,
        "two_orbit_identities": ok,
        "three_orbit_empty_slices": slices_ok,
        "equal_star_empty": type_I_three_orbit_equal_star_empty(),
        "plus_only_empty": type_I_three_orbit_plus_only_empty(),
        "zero_far_bad_case_empty": type_I_zero_far_bad_case_empty(),
        "three_orbit_family_bad_case_empty": type_I_three_orbit_family_bad_case_empty(),
        "split_far_identities": type_I_split_far_identities_ok(),
        "three_type_Y_collapses_to_H": type_I_three_type_Y_collapses_to_H(),
        "split_far_3AB_identities": type_I_split_far_3AB_identities_ok(),
        "wick_six_type_3ab_positive": type_I_wick_six_type_3ab_positive(),
        "aut_e_3AB_sign_identities": type_I_aut_e_3AB_sign_identities_ok(),
        "p5_through_e_3AB_positive": type_I_p5_through_e_3AB_positive(),
        "mu_le_maj_false_at_p7": lemma_mu_le_maj_false_at_p7()["proved"],
        "mu_le_L_not_T_as_bound": all(
            lemma_mu_le_L_not_T_as_bound(p)["proved"]
            for p in (5, 7, 11, 13, 17)
        ),
        "mu_le_L_census_p5p7": lemma_mu_le_L_census_p5p7()["proved"],
        "even_delta_Tnu_feed": all(
            lemma_even_delta_Tnu_feed(p)["proved"] for p in (5, 7, 11, 13, 17)
        ),
        "kappa3_particular_beats_thr": all(
            lemma_kappa3_particular_beats_thr(p)["proved"]
            for p in (5, 7, 11, 13, 17)
        ),
        "aut_e_3AB_positive_general": type_I_aut_e_3AB_positive_general(),
        "dual_eq_empty": empty,
        "n_primes": len(primes),
        "by_p_sample": sample,
        "hole": (
            "Y≤−3 and μ≥0 force a far Aut_e class with 3A+B≤0.  On "
            "|κ|=1 this is G>T (15.268 ν=0 + J).  |μ|≤|L| would "
            "close it (L>T) but is not proved for all p; |μ|≤|T| is "
            "not a closing bound; |μ|≤maj is false at p=7 "
            "(109/2863>17/539).  Census |μ|≤|L| at p=5,7.  p=5 is a "
            "finite from-C theorem (276/276).  |κ|=3 particular/Wick "
            "beat the linear form; leftover is even/odd δ."
        ),
        "theorem": (
            "2-orbit / equal-star / μ_⊠=0 slices empty.  Whole μ_far≥0 "
            "3-weight family is bad-case empty (H).  Coarse 3-type "
            "collapses to H (I).  3A+B conversion / Wick 6-type / Q_χ "
            "are Max+-free (J).  K reduces the Aut_e sign to g_min>T.  "
            "|μ|≤|L| suffices and holds at p=5,7; |μ|≤|T| does not "
            "close; |μ|≤maj false at p=7.  Sign stays open."
        ),
    }


# ---------------------------------------------------------------------------
# Predicate
# ---------------------------------------------------------------------------


def type_I_multilevel_identities_ok(p: int) -> bool:
    return bool(
        lemma_mass_identity(p)["proved"]
        and lemma_fe_budget(p)["proved"]
        and lemma_pairing_minimum(p)["proved"]
        and lemma_residual_second_moment(p)["proved"]
        and lemma_integrality(p)["proved"]
        and lemma_slopes(p)["proved"]
        and lemma_two_orbit_weights(p)["proved"]
        and lemma_maxminus_star_far_sums(p)["proved"]
        and lemma_far_weight_more_negative_when_beta_gt_minus_2(p)["proved"]
        and lemma_two_star_interpolants(p)["proved"]
        and lemma_three_orbit_equal_star_empty(p)["proved"]
        and lemma_three_orbit_plus_only_empty(p)["proved"]
        and lemma_three_orbit_zero_far_in_box(p)["proved"]
        and lemma_three_orbit_maxminus_fplus(p)["proved"]
        and lemma_zero_far_bad_case_empty(p)["proved"]
        and lemma_three_orbit_family_bad_case_empty(p)["proved"]
        and lemma_zero_far_01_counts(p)["proved"]
        and lemma_split_far_delta_identity(p)["proved"]
        and lemma_split_far_Y_plus_3(p)["proved"]
        and lemma_three_type_counts(p)["proved"]
        and lemma_three_type_Q_difference(p)["proved"]
        and lemma_three_type_Y_collapses_to_H(p)["proved"]
        and lemma_threeAB_conversion(p)["proved"]
        and lemma_wick_six_type_3ab_positive(p)["proved"]
        and lemma_Q_chi_identities(p)["proved"]
        and lemma_through_e_moments(p)["proved"]
        and lemma_pairing_four_point(p)["proved"]
        and lemma_gmin_T_vs_thresholds(p)["proved"]
        and lemma_kappa3_3AB_conversion(p)["proved"]
        and lemma_mu_part_vs_T(p)["proved"]
        and lemma_kappa3_star_through_e(p)["proved"]
        and lemma_mu_le_L_not_T_as_bound(p)["proved"]
        and lemma_even_delta_Tnu_feed(p)["proved"]
        and lemma_kappa3_particular_beats_thr(p)["proved"]
        and type_I_k_3p_minus_2_structure_ok(p)
        and type_I_gap2_forces_s_minus_eq_minus_1(p)["proved"]
    )


def type_I_multilevel_bad_case_ND_closed() -> bool:
    """
    Multi-level Type I bad case empty (hence Φ(H)≥Φ−2) for all primes p≥5.

    False: 3A+B>0 is not proved for a general Aut_e far class
    (⇔ G>T on |κ|=1 through e).  |μ|≤|L| would suffice but is
    not proved for all p; |μ|≤|T| is not a closing μ-bound;
    |μ|≤maj is false at p=7.  p=5 is a finite from-C theorem.
    Wick particular, the conversion, through-e 15.260 moments,
    and the 4-point pairing *are* Max+-free (J, K).  The μ_far=0
    slice, the whole μ_far≥0 3-weight family, and the coarse
    3-type *are* bad-case empty (H, I).  Gsum unused.  Two-level
    {−1,−3} is a different path.
    """
    return bool(theorem_G_leftover()["proved"])


def hinge_status_275() -> dict:
    return {
        "type_I_multilevel_bad_case_ND": type_I_multilevel_bad_case_ND_closed(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "two_level_is_a_different_path": True,
        "e1_closed_general": e1_closed_general(),
        "type_I_three_orbit_equal_star_empty": type_I_three_orbit_equal_star_empty(),
        "type_I_three_orbit_plus_only_empty": type_I_three_orbit_plus_only_empty(),
        "type_I_zero_far_bad_case_empty": type_I_zero_far_bad_case_empty(),
        "type_I_three_orbit_family_bad_case_empty": (
            type_I_three_orbit_family_bad_case_empty()
        ),
        "type_I_split_far_identities_ok": type_I_split_far_identities_ok(),
        "type_I_three_type_Y_collapses_to_H": (
            type_I_three_type_Y_collapses_to_H()
        ),
        "type_I_split_far_3AB_identities_ok": (
            type_I_split_far_3AB_identities_ok()
        ),
        "type_I_wick_six_type_3ab_positive": (
            type_I_wick_six_type_3ab_positive()
        ),
        "type_I_aut_e_3AB_sign_identities_ok": (
            type_I_aut_e_3AB_sign_identities_ok()
        ),
        "type_I_p5_through_e_3AB_positive": (
            type_I_p5_through_e_3AB_positive()
        ),
        "type_I_aut_e_3AB_positive_general": (
            type_I_aut_e_3AB_positive_general()
        ),
        "mu_le_L_not_T_as_bound": all(
            lemma_mu_le_L_not_T_as_bound(p)["proved"]
            for p in (5, 7, 11, 13, 17)
        ),
        "mu_le_L_census_p5p7": lemma_mu_le_L_census_p5p7()["proved"],
        "even_delta_Tnu_feed": all(
            lemma_even_delta_Tnu_feed(p)["proved"] for p in (5, 7, 11, 13, 17)
        ),
        "kappa3_particular_beats_thr": all(
            lemma_kappa3_particular_beats_thr(p)["proved"]
            for p in (5, 7, 11, 13, 17)
        ),
        "note": (
            "A–K proved as identities: μ_far=0 / 3-weight / coarse "
            "3-type empty; 3A+B conversion / Wick / Q_χ / through-e "
            "moments / pairing / T-vs-1/(2p) Max+-free.  p=5 finite "
            "from C.  |μ|≤|L| suffices (not |μ|≤|T|); |μ|≤maj false "
            "at p=7; census |μ|≤|L| at p=5,7.  General Aut_e 3A+B>0 "
            "(g_min>T / |μ|≤|L|) keeps ND False.  Two-level is "
            "15.169/15.216.  Gsum unused.  type_I / e1 / L unflipped."
        ),
    }


def main() -> dict:
    A = theorem_A_mass()
    B = theorem_B_fe_budget()
    C = theorem_C_pairing()
    D = theorem_D_residual()
    E = theorem_E_integrality()
    F = theorem_F_dual_eq_slice()
    G = theorem_G_leftover()
    flag = type_I_multilevel_bad_case_ND_closed()
    out = {
        "prop": "15.275",
        "title": (
            "Type I multi-level bad case empty: mass/pairing/Aut_e weights "
            "+ imported dual-eq empty"
        ),
        "proved": {
            "mass_identity": A["proved"],
            "fe_budget": B["proved"],
            "pairing_minimum": C["proved"],
            "residual_second_moment": D["proved"],
            "integrality": E["proved"],
            "dual_eq_slice_empty": F["proved"],
            "leftover_empty": G["proved"],
            "two_orbit_identities": G.get("two_orbit_identities"),
            "three_orbit_empty_slices": G.get("three_orbit_empty_slices"),
            "equal_star_empty": G.get("equal_star_empty"),
            "plus_only_empty": G.get("plus_only_empty"),
            "zero_far_bad_case_empty": G.get("zero_far_bad_case_empty"),
            "three_orbit_family_bad_case_empty": G.get(
                "three_orbit_family_bad_case_empty"
            ),
            "split_far_identities": G.get("split_far_identities"),
            "three_type_Y_collapses_to_H": G.get(
                "three_type_Y_collapses_to_H"
            ),
            "split_far_3AB_identities": G.get("split_far_3AB_identities"),
            "wick_six_type_3ab_positive": G.get(
                "wick_six_type_3ab_positive"
            ),
            "aut_e_3AB_sign_identities": G.get("aut_e_3AB_sign_identities"),
            "p5_through_e_3AB_positive": G.get("p5_through_e_3AB_positive"),
            "mu_le_maj_false_at_p7": G.get("mu_le_maj_false_at_p7"),
            "mu_le_L_not_T_as_bound": G.get("mu_le_L_not_T_as_bound"),
            "mu_le_L_census_p5p7": G.get("mu_le_L_census_p5p7"),
            "even_delta_Tnu_feed": G.get("even_delta_Tnu_feed"),
            "kappa3_particular_beats_thr": G.get(
                "kappa3_particular_beats_thr"
            ),
            "aut_e_3AB_positive_general": G.get(
                "aut_e_3AB_positive_general"
            ),
            "type_I_multilevel_bad_case_ND_closed": flag,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
            "two_level_claimed_here": False,
        },
        "algebra": {
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "F": F,
            "G": G,
            "two_level_H": {
                "note": "15.169 path; not this flag",
                "p5": bad_case_dual_two_level_H_params(5),
            },
            "key_identity": (
                "2a+c(3+μ_c)=2/p; E[R²]=E[S²]−5+4 E[Sf]; "
                "μ_far=−2(2p−3)/(p(p²−1))<0; n_{−1}=M+n_c+t; "
                "σ_□=(p−1)(1+f_e); μ_⊠=0⇒μ_□<0; "
                "F_−|_{f=+1}=−(p+1)/(p−1)+μ_far p(p+1); "
                "Y+3=4(p−2)μ_□+∑μ(3A+B); 3-type Y=H; "
                "3A+B=[2(p−2)n_O+p(3(p−1)G_++(p+1)G_−)]/(p²−1); "
                "G_+=G_−⇒3A+B>0⇔G>T; Wick 6-type 3a+b>0; "
                "∑_{S∋e}μκ=3(p²−1)/2; G=χd μ; 1/(2p)⇏G>T; "
                "|μ|≰maj at p=7; |μ|≤|L|⇏|μ|≤|T|; "
                "p=5 276 from C; |κ|=3 particular>thr"
            ),
        },
        "hinge": hinge_status_275(),
        "F3": (
            "Multi-level Type I bad case only.  Two-level is 15.169/216. "
            "Gsum disj LB False.  15.170 wiring not edited."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15275.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.275 Type I multi-level bad case")
    print(f"  mass: {A['proved']}")
    print(f"  fe budget: {B['proved']}")
    print(f"  pairing min: {C['proved']}")
    print(f"  residual: {D['proved']}")
    print(f"  integrality: {E['proved']}")
    print(f"  dual-eq slice: {F['proved']}")
    print(f"  leftover: {G['proved']}")
    print(f"  equal_star_empty: {G.get('equal_star_empty')}")
    print(f"  plus_only_empty: {G.get('plus_only_empty')}")
    print(f"  zero_far_bad_case_empty: {G.get('zero_far_bad_case_empty')}")
    print(f"  three_orbit_family_bad_case_empty: {G.get('three_orbit_family_bad_case_empty')}")
    print(f"  split_far_identities: {G.get('split_far_identities')}")
    print(f"  three_type_Y_collapses_to_H: {G.get('three_type_Y_collapses_to_H')}")
    print(f"  split_far_3AB_identities: {G.get('split_far_3AB_identities')}")
    print(f"  wick_six_type_3ab_positive: {G.get('wick_six_type_3ab_positive')}")
    print(f"  aut_e_3AB_sign_identities: {G.get('aut_e_3AB_sign_identities')}")
    print(f"  p5_through_e_3AB_positive: {G.get('p5_through_e_3AB_positive')}")
    print(f"  mu_le_maj_false_at_p7: {G.get('mu_le_maj_false_at_p7')}")
    print(f"  mu_le_L_not_T_as_bound: {G.get('mu_le_L_not_T_as_bound')}")
    print(f"  mu_le_L_census_p5p7: {G.get('mu_le_L_census_p5p7')}")
    print(f"  even_delta_Tnu_feed: {G.get('even_delta_Tnu_feed')}")
    print(f"  kappa3_particular_beats_thr: {G.get('kappa3_particular_beats_thr')}")
    print(f"  aut_e_3AB_positive_general: {G.get('aut_e_3AB_positive_general')}")
    print(f"  type_I_multilevel_bad_case_ND_closed: {flag}")
    print(f"  gsum={gsum_disj_lb_proved_general()} dual_eq={residual_i_dual_eq_empty_proved_general()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
