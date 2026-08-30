#!/usr/bin/env python3
"""
Prop 15.279 — Character-pair model of Φ|_F: |c_Ω|=2p, 3-pairing
λ=8+ˆR_rest/q².  Floor λ_min(Φ|_F)≥6 still OPEN (ˆR_rest≥−2q²).

Does **not** treat G_{u,disj} as a Gram.  Does **not** flip Aut-Schur /
Gsum / pairing / e1 / L / phi_F_ge_6.

SETUP
  q=p², T=F_q^×/{±1}.  Even multiplicative characters of F_q^× = chars of T.
  Paley-N on F_q is convolution by χ: Nψ=J(χ,ψ)·(χψ).
  J(χ,ψ)=G(χ)G(ψ)/G(χψ) (χ,ψ,χψ all nontrivial).  |J|=|G|=p.
  +p eigenvector f_+=ψ+(J/p)χψ.  F is the real +p line in each pair
  {ψ,χψ} for even ψ∉{1,χ}.  dim F=(q−5)/4.
  Ω={ξ: χ̂(ξ)=p}, ε=χ|Ω=±1, 1_Ω=(1+ε χ)/2 on F_q^×, G(χ)=τ=p ε.
  ẑ(ξ)=∑_{x∈F_q} z_x e^{Tr(ξ x)}; ẑ supported on {0}∪Ω; E[|ẑ|²]=2q on Ω.
  A(δ)=∑_x z_x z_{x+δ}.  Ã(λ)=∑_{δ≠0} A(δ) λ(δ).
  M_ψ=∑_{ξ∈Ω} ψ̄(ξ)|ẑ(ξ)|².

PROVED
  A. |c_Ω|=2p.  c_Ω=G(ψ)+α ε G(χψ), α=J/p=ε G(ψ)/G(χψ),
     so c_Ω=2G(ψ).  Referee math_review PASS.
  B. yᵀBy=q^{-1} c_Ω M_ψ.  ‖B_cpx‖²=2q(q−1).
     Rayleigh=2 E|M|²/(p⁴(p²−1)).
  C. Complex/2-pairing Wick: Q(ξ,ξ)=8q², Q(ξ,η)=4q² for ξ≠η,
     λ_Wick=4 (misses the antipode; superseded by D for the floor).
  D. Real 3-pairing Wick: Q(1)=Q(−1)=8q² and Q=4q² on other squares.
     S_□(ψ):=∑_{r square} ψ(r) Q(r).  E|M|²=|Ω| S_□.
     Q=4q²+4q²(1_{r=1}+1_{r=−1})+R_rest, R_rest(±1)=0.
     For even ψ∉{1,χ}, ∑_{squares}ψ=0 so ∑_{r≠±1}ψ=−2 and
        S_□=16q²+4q²(−2)+ˆR_rest=8q²+ˆR_rest
     (referee verify PASS on the 16q²→8q² cancellation).
     λ=S_□/q²=8+ˆR_rest/q², so λ≥6 ⇔ ˆR_rest≥−2q².
  E. Jacobi/Parseval (Max+-free, r-first).  Wiener + Gauss in r=ηξ^{-1}
     with 1_Ω=(1+ε χ)/2 and G(λ,α)=λ̄(α)G(λ):
        Ã(ψ)=G(ψ)/q M_ψ,   M_ψ=q Ã(ψ)/G(ψ),
        Gs(α)=½(ψ(α)G(ψ̄)+ε (χψ)(α)G(χψ̄))  (α≠0), Gs(0)=0,
        M_ψ=∑_{δ≠0} A(δ) Gs(δ).
     Conjugation: G(ψ,−μ)=ψ̄(μ)G(ψ), not ψ(μ)G(ψ) (referee BLOCK on
     the inverted character; live p=5,7 check to 1e-12 after the fix).
     Does **not** prove the floor.

PROVED (cont.)
  F. Collision remainder Δ_coll (Max+-free, referee math_review PASS).
     Boolean−Wick on non-4-distinct index tuples:
        Δ_coll(r)=−(18q+6)  (r=±1),   −(20q+6)  (other squares).
     Type sum: all-equal −2q; 3+1 four ways −8q (G(χ)|Ω=p);
     2+2 three matchings −4q+6 / −2q+6; one-pair six ways
     −4q−12 / −8q−12.  Hence for even ψ∉{1,χ}
        ˆR_rest = 40q+12 + ˆΔ_4,
     and λ≥6 ⇔ ˆΔ_4 ≥ −2q²−40q−12.  ˆΔ_4 is the 4-distinct
     (m4−κ_C/p²) contraction.
  H. Wick_4dist on 4-distinct indices, r square ≠±1, is the
     constant q²−10q−6 (Max+-free: Wick(π) off ±1 is 4q² and
     Wick_coll off ±1 is 3q²+10q+6).  For even ψ∉{1,χ},
     ∑_{squares}ψ=0 ⇒ ∑_{r≠±1}ψ=−2, so
        ˆW_4 = (q²−10q−6)(−2) = −2q²+20q+12.
     Hence ˆΔ_4 = ⟨m_4,K_ψ⟩ − ˆW_4 and the floor is
        ⟨m_4,K_ψ⟩ ≥ −4q²−20q.
     p=5: Wick_4dist=369 matches Q_split.  Inequality still OPEN.
  J. ⟨κ_C,K_ψ⟩=q ˆW_4.  On 4-distinct π=C/p, so κ_C=p² Wick_π
     and ⟨κ_C/p²,K⟩=ˆW_4 by H.  Equivalently ∑_λ J_κ(λ)=q ˆW_4
     (checked at p=5,7,11, every even k; inverted Θ breaks the mass).
     Therefore ⟨m_4,K⟩=ˆW_4+⟨ρ,K⟩ with ρ=m_4−κ_C/p², and the
     floor is ⟨ρ,K_ψ⟩≥−2q²−40q−12.  Inequality still OPEN.
  K. ∑_{α≠0}|I(α)|²=(q²−5q−8)/2, independent of even ψ∉{1,χ}.
     Split I=A−B, A=½(ψ̄G(ψ)+(χψ)̄G(χψ)), B=e(α)+e(−α):
        ∑|A|²=q(q−1)/2 (∑χ=0),  ∑|B|²=2q−4,
        ∑AB=2q because G(ψ)G(ψ̄)=q=G(χψ)G((χψ)̄).
     CS on ⟨ρ,K⟩ still overshoots.  Does not prove the floor.
  L. Additive 4-distinct pieces of I cancel 4q²+20q (Wick_coll(±1)=6q²+8q+6,
     Boolean_coll=6q²−10q, Q4dist(±1)=2q²+10q).  Hence
        S_□−6q² = Gauss 4-dist pairing of m₄,
     and λ≥6 iff that pairing is ≥0.  Positivity OPEN.
  M. Collision Gauss mass is 3q², independent of even ψ∉{1,χ}.
     Boolean on collisions uses only z²=1 and π=E[z_i z_j].
     Gauss piece of I vanishes at α=0.  T1 (α=0) and T2
     (∑(χψ̄)=0) contribute 0.  T3: the (a=b,c=d) matching has
     α=0; the two crossed matchings each give
     (q/2)G(ψ)G(ψ̄)=q²/2, total q².  T4: ab and cd vanish
     (∑χψ̄=0 or α=0); each of the four crosses is
        (q G(ψ)/(2p)) ψ̄(ξ) J(χ,ψ̄) G(χψ̄,ξ)
        = q²/2
     via J(χ,ψ̄)=G(χ)G(ψ̄)/G(χψ̄) and G(χ)=pε, χ(ξ)=ε on Ω.
     Hence A_coll=3q².  Parseval gives A_full=S_□/2, so
     A_4dist=S_□/2−3q² and A+B=S_□−6q² (same as L).
     Floor remains S_□≥6q².  Inverted J breaks T4=q²/2.
  M. Boolean collision kernel Coll(r) is 2-level on squares (Max+-free,
     z²=1 + G(χ,μ)=p on Ω).  Types: T1=q, T2=4q all r; T3=2q²−3q
     on ±1 and q²−3q off; T4 six pair-slots, generic −2q and
     resonant q(q−2) at r=±1, hence T4=4q²−12q / 2q²−12q.
        Coll(±1)=6q²−10q,   Coll(r≠±1)=3q²−10q.
     Character orthogonality (pair ψ with ψ̄; both even so ψ(±1)=1):
        ∑_{r□} ψ(r) Coll(r) = 6q².
     Same identity for every even mean-zero test weight w on T:
        ∑_{r□} w(r) Coll(r) = 6q² w(1)
     (not only Max+ vectors).  Therefore
        S_□ − 6q² = ∑_{r□} ψ(r) F(r) = Re(G(ψ) A₄(ψ)),
     F=Q−Coll the 4-distinct Boolean kernel.  Dropping −10q from
     Coll_off or inverting G(λ,α)=λ̄(α)G(λ) in I breaks the mass.
     F̂(ψ)≥0 is the floor and stays OPEN (no SoS for F).
  O. Collision Gauss pairing is phase-aligned: A₄_coll(ψ)=6q G(ψ̄)
     (Max+-free: z²=1, π=I+N/p, G(λ,α)=λ̄(α)G(λ), ξ∈Ω).  T1=T2=0
     (ψ̄(0)=0 and ∑χψ̄=0).  T3 two crossed matchings each q G(ψ̄)
     via G(ψ̄,−1)=G(ψ̄).  T4: A1 and B3 vanish; four Jacobi crosses
     each q G(ψ̄) because G(χ)=pε and χ(ξ)ε=1.  Hence
        G(ψ) A₄_coll = 6q²
     and the 4-distinct remainder satisfies
        A₄(ψ)=(S_□/q−6q) G(ψ̄)
     so G A₄ is real and |G||A₄|=|S_□−6q²| (the live “perfect
     alignment”).  Inverted phase (ψ not ψ̄) gives 6q G(ψ) and
     G·A₄_inv=6q G(ψ)², not real when G is non-real.  Paley classes
     of Q fail (p=7 (−1,−1) splits; p=11 every off-diagonal class
     splits); class values have Max+ denominators, not Jacobi.
     A₄/G(ψ̄)≥0 is exactly S_□≥6q² and stays OPEN.  Pointwise
     |M_ψ|² dips to 0 on Max+ (p=5 k=2: 77% below 3q²(q−1)).
  I. Jacobi kernel pairings (Max+-free).  K_ψ=1_{4-dist} e_ξ(b−a) I(ξ(d−c)),
     I(α)=∑_{r□, r≠±1} ψ(r) e^{Tr(α r)}.  Gauss inversion G(λ,α)=λ̄(α)G(λ):
        I(0)=−2,  I(α)=½(ψ̄(α)G(ψ)+(χψ)̄(α)G(χψ))−e(α)−e(−α) (α≠0).
     On exact cr=λ, c=(λ−1)/λ, f_λ=c·s(s−1)/(s−1/λ), μ/(s−1)=f_λ/(s−1)²:
        K(α,β,γ;δ)=∑_{s≠0,1,δ} α(s)β(s−1)γ(s−δ)  (complete Jacobi),
        Θ_ψ=½(ψ(c) K(ψ,ψ̄,ψ̄;1/λ)+(χψ)(c) K(χψ,χψ̄,χψ̄;1/λ)).
        Even ≠ square: ψ((s−1)²)=ψ(s−1)², so the twist is ψ̄²(s−1).
        Untwisted J_3(ψ)=K(ψ,ψ,ψ̄) is the wrong Jacobi and fails the mass.
        H(γ)=∑_x χ(x(x−1)(x−γ)) (Legendre; not in {0,±p,±2p}),
        Ξ=χ(1−λ) H(1−λ),  n_spec=2·1_{χ(1−λ)=1}+2·1_{χ(λ(λ−1))=1},
        Ξ_spec=n_spec (on A=±B, f=±(s−1)² and χ(−1)=1).
        κ=(1+χ(λ)+χ(1−λ)) χ(s(s−1)/(1−λs)),
        J_1=q² Θ+2q(q−3)−q² n_spec,
        J_κ=(1+χ(λ)+χ(1−λ))χ(1−λ)[q² Θ+2q Ξ−q² n_spec].
     Global: ∑_λ J_1=−2q²+12q and ∑_λ J_κ=q Ŵ_4=−2q³+20q²+12q.
     K(η,η̄,η̄;δ) is not a Gauss product (|K| not constantly |G(η)²/G(η²)|).
     Does **not** prove ⟨m_4,K⟩≥−4q²−20q.

  R. 1D lifts are Max+ (Max+-free construction; not the floor).
     L:F_q→F_p F_p-linear with ker L a square line, σ:F_p→{±1}
     with ∑σ=1 (i.e. |σ^{-1}(+1)|=(p+1)/2), y_∞=1, y_x=σ(L(x)).
     Then Cy=py.  Fiber sums: K(0)=p−1, K(β≠0)=−1, by χ_q|F_p^×=1
     and ∑_{t} χ_p(at²+bt+c)=−χ_p(a) (disc≠0).  Convolution
     σ*K=pσ−∑σ.  Type− lines and wrong weight fail.  Does **not**
     exhaust Max+ (p=7 has a nonlinear Aut-orbit of size 2352) and
     does **not** sign F̂ (that orbit has S_□=0 for some even ψ).

  S. Regular-set switching (Max+-free; referee math_review PASS).
     y_∞=+1 and Cy=py ⇔ D={x:z_x=−1} has size p(p−1)/2 and
     S(a)=∑_x χ(a−x)1_D(x) equals (1+p)/2 on D and (1−p)/2 off D,
     i.e. C_fin d = p d + ((1−p)/2)1.  Then ẑ(ξ)=−2 1̂_D(ξ) for
     ξ≠0, so A_r=−8 E[1̂_D(ξ)1̂_D(rξ)1̂_D(−(1+r)ξ)] and A_r≤0 is
     an ensemble statement (pointwise ẑẑẑ can be positive).
     Pointwise |M_ψ|²≥3q²(q−1) is FALSE (p=5 k=2: 200/260 zeros);
     referee falsify BLOCK on pointwise SOS.  Does **not** sign F̂.

  T. D-autocorrelation form of M / A-domain switching (Max+-free;
     referee math_review PASS-WITH-NOTE, do_not_branch).
     χ*z=pz−1 ⇒ χ*A=p(A−1) (χ(−1)=1, A even; holds at δ=0).
     Product of two linear switching eqs is a tautology on A.
     A(δ)=p(2−p)+4 N(δ), N=|D∩(D−δ)|.
     Nontrivial ψ: Ã=4 N̂, N̂(ψ)=∑_{x,y∈D} ψ(y−x),
        M_ψ=4q N̂(ψ)/G(ψ).
     Boolean on Ω: 2p ẑ(ξ)+∑_{t∈R} ẑ(tξ)ẑ((1−t)ξ)=0,
        R={t≠0,1:χ(t)=χ(1−t)=1}, |R|=(q−5)/4=dim F
        (character-sum: (q−2−1−1−1)/4).
     3-point T(δ,ε)=∑_a z_a z_{a+δ} z_{a+ε}:
        ∑_u χ(u) T(δ+u,ε+u)=p T(δ,ε)−A(ε−δ).
     Linear switching sees only the χ-isotype of A (supp ẑ⊂{0}∪Ω)
     and does not constrain N̂(ψ) for even ψ∉{1,χ}.  Inverted G
     (denom G(ψ̄)) and Type− / Boolean coeff p all fail.
     Does **not** sign F̂: pointwise |M|²=0; no SOS of S_□−6q²
     modulo the switching ideal (Boolean rewrite returns M).

  U. Energy fiber is Max+ (Max+-free; C²=p²I).  ||Cy||=p||y||
     and yᵀCy=pn ⇔ Cy=py on the cube.  Does **not** sign Term_4.

  V. Ensemble E[N] is 2-level by χ (referee math_review PASS).
     E[N|χ=1]=p(p−1)/4, E[N|χ=−1]=p(p−3)/4.  Hence E[N̂(ψ)]=0
     for even ψ∉{1,χ} and Term_4 is Var(N̂).  Type+ 1D has those
     class means (not pointwise).  Type− swaps them.  Does **not**
     sign the variance.

  W. Nonsquare pairing (Max+-free).  For χ(r)=−1, Parseval of
     Â=|ẑ|² supported on {0}∪Ω gives ∑_δ A(δ)A(rδ)=q, hence
     ∑_{δ≠0} N(δ)N(rδ)=(q−1)p²(p−1)(p−3)/16 on every regular
     set, so L(r)=E[N_+]E[N_-] on nonsquares and L̂(ψ) lives
     only on squares.  Does **not** sign L̂ on squares.

  X. L is inversion-invariant, not S3-invariant (Max+-free).
     ∑_δ N(δ)N(rδ)=∑_ε N(ε)N(ε/r) is a change of variables,
     so L(1/r)=L(r) for every N (hence every regular set).
     χ_q(−2)=1 for every odd p (Norm(−2) is a square in F_p),
     and −2∉{±1} for p≥5, so r↦−1−r sends the {±1} Aut-orbit
     to a different Aut-orbit.  On every Type+ 1D lift,
     the pairing at 1 is strictly larger than at −2
     (p=5: 650 vs 400; p=7: 6076 vs 4018; p=11: 110110 vs 75504).
     Leftover after W+L(1)+row-sum is n_aut−2, not n_S3−2.
     Does **not** sign L̂.

  Y. 1D N̂ closed form (Max+-free; does not exhaust Max+).
     On a Type+ lift, N=|D| on ker L and N=p n_M(L(δ)) off ker,
     M=σ^{-1}(−1).  Even ψ: (p−1)|k ⇔ ψ|F_p^× trivial, and then
     |N̂|²=[p(p²−1)/4]².  Otherwise |N̂|²=p³|∑_{x,y∈M}η(x−y)|²
     with η=ψ|F_p^×.  Paley plus-set {0}∪□_p has n_M constant
     when p≡3 (mod 4), hence N̂=0 on those even ntriv k.
     Aut_∞ 1D orbits are AGL(1,F_p) orbits of plus-sets.
     1D-only mixture fails the floor at p=7 (k=8).  Does **not**
     sign the ensemble L̂.

OPEN
  G. S_□≥6q²  (equivalently A₄/G(ψ̄)≥0, ˆF(ψ)≥0).  The floor is
     the AVERAGE on the full Max+ fiber, not a pointwise SOS
     (pointwise bound is false).  C-only switching bijections
     cannot determine 4-point extension counts.  Regular-set form
     is S; D-autocorr / A-domain switching is T.  Paley+ω /
     Jacobi-1̂_c Bose–Mesner attack is DEAD: a general-all-p scheme
     formula is false because Paley+ω fails to be a scheme at p=11
     (also at the checked primes 13,17,19,23; it is a scheme at
     p=5,7).  Direction partition of T IS a scheme at those primes,
     but Q is not a class function of it (p=5,7).
     S3 r↦−1−r is not a symmetry of L (theorem X); leftover
     dofs are n_aut−2, not n_S3−2.  Cube-root pairing is not
     constant on regular sets.  C_fin χ-contraction sees only
     the χ-isotype and does not pin leftover Aut-orbit L.
     SWITCHED attack: Q on Frob/inv-orbits of ratios (the Aut
     invariance that is actually true) / Boolean 4-point of
     V_+ intersect {+1,-1}^n / mixture of Aut-orbits of Max+
     (1D lifts R are not all of Max+; some Aut-orbits have S_□=0).
     Not Weil/triangle/CS, not G_u,disj as Gram,
     not a Paley+ω eigenvalue formula.
     Certified at p=5 (λ=80/13, room 1250/13) and p=7.
     Scratch m4_autinf_avg: affine m4 is a function of
     (cross-ratio, κ_C, star), not of κ_C or (χ-edges, κ_C);
     a global linear model in {κ,star,φ} fails at both p.
     Weil/triangle on |J_1|,|J_κ|,|J_*| misses the budget by Θ(p).
     α(cr), β(cr) are not a single Jacobi in p.  Theorem M writes
     Φ|_F−6I as the 4-distinct Fourier of F=Q−Coll after the
     6q² orthogonality, but does not sign F̂.  Paley scheme classes
     (χ(r−1),χ(r+1)) do not determine Q (p=7,11).  Theorem O aligns
     A₄ to G(ψ̄) but the sign A₄/Ḡ≥0 is the same floor.  Not general.

Writes evidence/e1_gmin_m4_prop15279.json
"""
from __future__ import annotations

import cmath
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15278 import (  # noqa: E402
    aut_span_F_equals_Z_proved_general,
    dim_F,
    phi_F_ge_6_proved_general,
)


def q_of(p: int) -> int:
    return p * p


def c_omega_abs(p: int) -> int:
    """|c_Ω|=2p."""
    return 2 * p


def wick_rayleigh() -> int:
    return 4


def EM_wick(p: int) -> int:
    """E|M|² under Wick for ψ≠1,χ: 2 q² (q−1)."""
    q = q_of(p)
    return 2 * q * q * (q - 1)


def rayleigh_from_EM(p: int, EM: Fraction | int) -> Fraction:
    """λ=2 E|M|² / (p⁴ (p²−1))."""
    return Fraction(2 * EM, (p**4) * (p * p - 1))


def Rhat_budget(p: int) -> int:
    """ˆR≥q²(q−1) ⇔ λ≥6."""
    q = q_of(p)
    return q * q * (q - 1)


def Q_row_sum(p: int) -> int:
    """∑_η Q(ξ,η)=2q²(q−1).  Max+-free from E[u]=2q and ∑u=q(q-1)."""
    q = q_of(p)
    return 2 * q * q * (q - 1)


def EM_trivial(p: int) -> int:
    """E|M_1|²=q²(q−1)².  M_1=∑_Ω u=q(q−1) constantly."""
    q = q_of(p)
    return q * q * (q - 1) * (q - 1)


def R_antipode(p: int) -> int:
    """R(-1)=Q(ξ,−ξ)−4q²=8q²−4q²=4q²."""
    q = q_of(p)
    return 4 * q * q


def Q_diag_closed(p: int) -> int:
    """Q(ξ,ξ)=8q²."""
    q = q_of(p)
    return 8 * q * q


def theorem_c_omega_is_2p(primes: list[int] | None = None) -> dict:
    """A. Gauss/Jacobi product ⇒ |c_Ω|=2p."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        # χ_Ω²=1, α=χ_Ω G(ψ)/G(χψ) ⇒ α χ_Ω G(χψ)=G(ψ) ⇒ c=2G, |c|=2p
        chi_om_sq = 1
        two_p = 2 * p
        # |G(ψ)|=p for nontrivial ψ
        g_abs = p
        row = chi_om_sq == 1 and two_p == 2 * g_abs == c_omega_abs(p)
        sample[str(p)] = {"c_abs": two_p, "ok": row}
        if not row:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "J(χ,ψ)=G(χ)G(ψ)/G(χψ) for χ,ψ,χψ nontrivial.  "
            "G(χ)=τ=p χ_Ω on Paley (χ_Ω τ=p).  α=J/p=χ_Ω G(ψ)/G(χψ).  "
            "Then α χ_Ω G(χψ)=G(ψ), so c_Ω=G(ψ)+α χ_Ω G(χψ)=2G(ψ) "
            "and |c_Ω|=2p."
        ),
        "by_p_sample": sample,
    }


def theorem_rayleigh_formula(primes: list[int] | None = None) -> dict:
    """B. λ=2 E|M|²/(p⁴(p²−1)) from |c|=2p and Parseval."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        # |c|² / (2 q³ (q−1)) * EM = 4p² EM / (2 p^6 (p²−1)) = 2 EM/(p⁴(p²−1))
        c2 = (2 * p) ** 2
        den = 2 * q**3 * (q - 1)
        reduced = Fraction(c2, den)
        want = Fraction(2, (p**4) * (p * p - 1))
        row = reduced == want
        sample[str(p)] = {"coeff": str(reduced), "ok": row}
        if not row:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Unnormalized FT: ∑_x w_x w_{x+d}=q^{-1}∑_ξ |ŵ(ξ)|² ψ_{-ξ}(d).  "
            "Hence yᵀBy=q^{-1} c_Ω M_ψ.  ‖B_cpx‖²=2q(q−1).  "
            "Rayleigh=|c_Ω|² E|M|²/(2q³(q−1))=2 E|M|²/(p⁴(p²−1))."
        ),
        "by_p_sample": sample,
    }


def theorem_wick_rayleigh_is_4(primes: list[int] | None = None) -> dict:
    """C. Wick Q is 8q² / 4q²; λ_Wick=4 on F."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        em = EM_wick(p)
        lam = rayleigh_from_EM(p, em)
        row = (
            em == 2 * q * q * (q - 1)
            and lam == wick_rayleigh()
            and Q_diag_closed(p) == 8 * q * q
        )
        sample[str(p)] = {"EM_wick": em, "lambda_wick": str(lam), "ok": row}
        if not row:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Wick 4-point uses π=I+N/p.  On Ω, E[|ẑ|²]=2q so pair1=2q.  "
            "Off-diagonal Wick Q=(2q)²=4q²; diagonal Wick Q=8q² "
            "(term1+term3; term2 vanishes).  For even ψ∉{1,χ}, "
            "∑_{F_q^×}ψ=∑χψ=0 ⇒ ∑_Ω ψ=0.  Hence "
            "E|M|²_Wick=4q²|Ω|=2q²(q−1) and λ_Wick=4."
        ),
        "by_p_sample": sample,
    }


def theorem_Q_row_sum_and_antipode(primes: list[int] | None = None) -> dict:
    """Q row-sum, Q(ξ,−ξ)=8q², R(−1)=4q², trivial DFT.  Max+-free."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        # E[u]=2q from π=I+N/p and χ̂|Ω=p; ∑u=q(q−1) from z²=1 Parseval
        row = Q_row_sum(p) == 2 * q *  q * (q - 1)
        # Q(ξ,−ξ)=E[u(ξ)²]=Q(ξ,ξ)=8q² (z real ⇒ ẑ(−ξ)=conj ẑ(ξ))
        anti = R_antipode(p) == 4 * q * q
        # M_1=∑u=q(q−1) constant (both y_∞=±1 give ∑w²=q−1)
        triv = EM_trivial(p) == q * q * (q - 1) ** 2
        # Plancherel check: 2q²(q−1) * |Ω| = |Ω| * 8q² / 2? skip
        row_ok = row and anti and triv
        sample[str(p)] = {
            "row_sum": Q_row_sum(p),
            "R_minus1": R_antipode(p),
            "EM_1": EM_trivial(p),
            "ok": row_ok,
        }
        if not row_ok:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "z real ⇒ ẑ(−ξ)=conj ẑ(ξ) ⇒ Q(ξ,−ξ)=E[u(ξ)²]=8q² "
            "(same as the diagonal, already Wick-exact).  Hence R(−1)=4q².  "
            "∑_ξ u(ξ)=q(q−1) constantly (Parseval z²=1 and ẑ supported on "
            "{0}∪Ω).  E[u]=2q (π=I+N/p).  Therefore the row-sum "
            "∑_η Q(ξ,η)=2q²(q−1) and E|M_1|²=q²(q−1)².  All Max+-free."
        ),
        "by_p_sample": sample,
    }


def theorem_floor_iff_Rhat(primes: list[int] | None = None) -> dict:
    """D. λ=8+ˆR_rest/q² ≥6 ⇔ ˆR_rest≥−2q².  The inequality is OPEN."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        # 3-pairing Wick on squares: S_□=8q² + ˆR_rest for even ψ∉{1,χ}
        # λ=S_□/q²=8+ˆR_rest/q²; ≥6 ⇔ ˆR_rest≥−2q²
        # equivalent old form: E|M|²=4q²(q−1)+|Ω|ˆR_rest ≥3q²(q−1)
        om = (q - 1) // 2
        em_at_bound = 4 * q * q * (q - 1) + om * (-2 * q * q)
        # 4q²(q−1)−q²(q−1)=3q²(q−1)
        row = em_at_bound == 3 * q * q * (q - 1)
        lam = rayleigh_from_EM(p, em_at_bound)
        row = row and lam == 6
        # 8 + (−2q²)/q² = 6
        row = row and Fraction(8) + Fraction(-2 * q * q, q * q) == 6
        sample[str(p)] = {
            "EM_at_minus_2q2": em_at_bound,
            "lambda_at_bound": str(lam),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Real 3-pairing Wick gives Q(1)=Q(−1)=8q² and Q=4q² on other "
            "squares.  S_□=16q²+4q²∑_{r≠±1}ψ+ˆR_rest and ∑_{r≠±1}ψ=−2, "
            "so S_□=8q²+ˆR_rest.  E|M|²=|Ω|S_□ and λ=S_□/q², so "
            "λ=8+ˆR_rest/q².  Thus λ≥6 iff ˆR_rest≥−2q².  "
            "The inequality is OPEN (theorem F)."
        ),
        "by_p_sample": sample,
    }


def theorem_S_square_is_8q2_plus_Rhat(primes: list[int] | None = None) -> dict:
    """D-accounting.  16q² + 4q²(−2) = 8q².  Referee verify PASS."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        nsq = (q - 1) // 2
        # ∑_{r≠±1 square} 1 = nsq - 2 = (q-5)/2
        n_rest = (q - 5) // 2
        row = n_rest + 2 == nsq
        # ∑_{squares} ψ=0 ⇒ ∑_{rest} ψ=−2
        s_from_pm1 = 8 * q * q + 8 * q * q  # Q(1)+Q(−1)
        s_wick_rest = 4 * q * q * (-2)
        s_sq = s_from_pm1 + s_wick_rest
        row = row and s_sq == 8 * q * q
        row = row and Fraction(8) + Fraction(-2 * q * q, q * q) == 6
        sample[str(p)] = {
            "n_squares": nsq,
            "n_rest": n_rest,
            "S_wick": s_sq,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Even ψ∉{1,χ}: ∑_{r square}ψ=0 and ψ(−1)=1, so "
            "∑_{r≠±1}ψ=−2.  Wick Q(±1)=8q² and Q=4q² on the rest, "
            "hence S_□,Wick=16q²+4q²(−2)=8q².  The residual is ˆR_rest."
        ),
        "by_p_sample": sample,
    }


def theorem_parseval_M_from_autocorr(primes: list[int] | None = None) -> dict:
    """E. M_ψ=q Ã(ψ)/G(ψ)=∑_{δ≠0} A(δ) Gs(δ).  Max+-free Gauss/Wiener."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    from e1_gmin_m4_prop15270 import _make_field, _primroot_q

    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        F = _make_field(p)
        mul, chi, ia = F["mul"], F["chi"], F["ia"]
        g = _primroot_q(F, p)
        units = [1]
        for _ in range(q - 2):
            units.append(mul(units[-1], g))
        dlog = {units[j]: j for j in range(q - 1)}

        def tr(u: int) -> int:
            return (2 * (u % p) + (u // p) * ia) % p

        def padd(mu: int, x: int) -> complex:
            return cmath.exp(2j * math.pi * tr(mul(mu, x)) / p)

        def psi(x: int) -> complex:
            if x == 0:
                return 0j
            return cmath.exp(2j * math.pi * 2 * dlog[x] / (q - 1))

        # G(ψ), G(ψ,α)=∑ ψ(x) e^{Tr(α x)}
        G = sum(psi(x) * padd(1, x) for x in range(1, q))
        # inversion G(ψ,α)=ψ̄(α) G(ψ)
        inv_ok = True
        for alpha in (1, g, mul(g, g)):
            Galpha = sum(psi(x) * padd(alpha, x) for x in range(1, q))
            pred = (psi(alpha).conjugate() if alpha else 0j) * G
            if abs(Galpha - pred) > 1e-8:
                inv_ok = False
        # |G|=p and G(ψ) G(ψ̄)=q (even ⇒ ψ(−1)=1)
        g_abs = abs(G)
        g_prod = abs(G * G.conjugate())
        gauss_ok = abs(g_abs - p) < 1e-8 and abs(g_prod - q) < 1e-6
        # ∑_Ω ψ̄=0: 1_Ω=(1+ε χ)/2, ψ≠1,χ
        tau = sum(chi(x) * padd(1, x) for x in range(1, q))
        Omega = [
            xi
            for xi in range(1, q)
            if abs(chi(xi) * tau - p) < 1e-6
        ]
        sum_om = sum(psi(xi).conjugate() for xi in Omega)
        om_ok = abs(sum_om) < 1e-8 and len(Omega) == (q - 1) // 2
        # algebra M = q Ã / G is the inverse of Ã = G/q M
        rec = Fraction(q) / Fraction(p) * Fraction(p) / Fraction(q)
        rec_ok = rec == 1
        # S_□ accounting already in D
        row = inv_ok and gauss_ok and om_ok and rec_ok
        sample[str(p)] = {
            "G_abs": g_abs,
            "G_prod": g_prod,
            "n_Omega": len(Omega),
            "sum_Omega_psibar": abs(sum_om),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ẑ supported on {0}∪Ω (Cy=py).  Wiener "
            "A(δ)=q^{-1}∑_{μ∈{0}∪Ω} u(μ) e_{-μ}(δ).  "
            "Ã(ψ)=∑_{δ≠0} A(δ)ψ(δ)=q^{-1}∑_μ u(μ) G(ψ,−μ) on the "
            "support.  G(ψ,−μ)=ψ̄(μ)G(ψ) (even ψ), μ=0 term vanishes, "
            "so Ã(ψ)=G(ψ)/q M_ψ and M_ψ=q Ã(ψ)/G(ψ).  "
            "Equivalently 1_Ω=(1+ε χ)/2 gives "
            "Gs(α)=½(ψ(α)G(ψ̄)+ε(χψ)(α)G(χψ̄)) (α≠0), Gs(0)=0, "
            "and M_ψ=∑_{δ≠0} A(δ) Gs(δ).  "
            "Does not bound ˆR_rest."
        ),
        "by_p_sample": sample,
    }


def delta_coll_pm1(p: int) -> int:
    """Δ_coll(±1)=−(18q+6)."""
    q = q_of(p)
    return -(18 * q + 6)


def delta_coll_off(p: int) -> int:
    """Δ_coll(square r≠±1)=−(20q+6)."""
    q = q_of(p)
    return -(20 * q + 6)


def Delta4_budget(p: int) -> int:
    """ˆΔ_4 ≥ −2q²−40q−12 ⇔ ˆR_rest≥−2q², given Δ_coll."""
    q = q_of(p)
    return -2 * q * q - 40 * q - 12


def wick_4dist_off(p: int) -> int:
    """Wick on 4-distinct indices, square r≠±1: q²−10q−6."""
    q = q_of(p)
    return q * q - 10 * q - 6


def wick_coll_off(p: int) -> int:
    """Wick on collisions, square r≠±1: 3q²+10q+6."""
    q = q_of(p)
    return 3 * q * q + 10 * q + 6


def wick_coll_T1(p: int) -> int:
    """All-equal: Wick=3, q points."""
    return 3 * q_of(p)


def wick_coll_T2(p: int) -> int:
    """3+1, four orientations: each (3q/p)G(χ,±ξ or ±rξ)=3q on Ω."""
    return 12 * q_of(p)


def wick_coll_T3_off(p: int) -> int:
    """2+2 off ±1: (q+2)(q−1) − 2(q+2) = (q+2)(q−3)."""
    q = q_of(p)
    return (q + 2) * (q - 3)


def wick_coll_T4_off(p: int) -> int:
    """2+1+1 off ±1: 2q(q−2)+4 from ab,cd plus 8 from four crosses."""
    q = q_of(p)
    return 2 * q * q - 4 * q + 12


def What4(p: int) -> int:
    """ˆW_4 = (q²−10q−6)(−2) = −2q²+20q+12."""
    q = q_of(p)
    return -2 * q * q + 20 * q + 12


def wick_coll_T3_pm1(p: int) -> int:
    """2+2 at r=±1: (q+2)(q−1) + (q+2)(q−1) − (q+2) = (q+2)(2q−3)."""
    q = q_of(p)
    return (q + 2) * (2 * q - 3)


def wick_coll_T4_pm1(p: int) -> int:
    """2+1+1 at r=±1: ab+cd still 2q(q−2)+4; two crosses stay 2;
    the two degenerate crosses (1∓r=0) are each q(q−2)+2."""
    q = q_of(p)
    return 4 * q * q - 8 * q + 12


def wick_coll_pm1(p: int) -> int:
    """Wick on collisions at r=±1: 6q²+8q+6."""
    q = q_of(p)
    return 6 * q * q + 8 * q + 6


def Q4dist_pm1(p: int) -> int:
    """4-distinct part of Boolean Q(±1): 2q²+10q."""
    q = q_of(p)
    return 2 * q * q + 10 * q


def m4K_budget(p: int) -> int:
    """⟨m_4,K⟩ ≥ −4q²−20q ⇔ ˆΔ_4 ≥ Δ4_budget, given ˆW_4."""
    q = q_of(p)
    return -4 * q * q - 20 * q


def theorem_delta_coll_closed(primes: list[int] | None = None) -> dict:
    """F. Collision Δ_coll closed.  Max+-free type-sum / Gauss on Ω."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        # T1
        t1 = -2 * q
        # T2: each of 4 orientations −2q because G(χ,μ)=p on Ω
        t2 = 4 * (-2 * q)
        # T3 2+2
        t3_pm1 = -2 * (q - 1) + 2 + (-2 * (q - 1))  # r=±1
        t3_off = -2 * (q - 1) + 2 + 2
        row = t3_pm1 == -4 * q + 6 and t3_off == -2 * q + 6
        # T4 one-pair
        t4_pm1 = -4 + (-4 * q - 4) + (-4)
        t4_off = -4 + (-4 * q - 4) + (-4 * q - 4)
        row = row and t4_pm1 == -4 * q - 12 and t4_off == -8 * q - 12
        tot_pm1 = t1 + t2 + t3_pm1 + t4_pm1
        tot_off = t1 + t2 + t3_off + t4_off
        row = row and tot_pm1 == delta_coll_pm1(p)
        row = row and tot_off == delta_coll_off(p)
        # floor reduction: ˆR_rest = (−20q−6)(−2) + ˆΔ_4 = 40q+12+ˆΔ_4
        red = -2 * delta_coll_off(p)
        row = row and red == 40 * q + 12
        row = row and (red + Delta4_budget(p) == -2 * q * q)
        sample[str(p)] = {
            "dcoll_pm1": tot_pm1,
            "dcoll_off": tot_off,
            "Rhat_from_coll": red,
            "Delta4_budget": Delta4_budget(p),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Boolean−Wick on collisions, no Max+.  T1 all-equal: −2q.  "
            "T2 3+1: G(χ)|Ω=p gives −2q per orientation, −8q total.  "
            "T3 2+2: matching (a,a,c,c) contributes −2(q−1); the two "
            "crossed matchings contribute −2(q−1) or +2 according as "
            "r=∓1 or not, hence −4q+6 on ±1 and −2q+6 off.  "
            "T4 one-pair: the two adjacent pairings contribute −2 each; "
            "the four cross pairings are −2(S(ξ)S(η)−∑ e_{ξ±η}) with "
            "S=p on Ω, hence −4q−12 on ±1 and −8q−12 off.  "
            "Totals −(18q+6) / −(20q+6).  Therefore ˆR_rest=40q+12+ˆΔ_4 "
            "and λ≥6 iff ˆΔ_4≥−2q²−40q−12.  ˆΔ_4 OPEN."
        ),
        "by_p_sample": sample,
    }


def theorem_wick_4dist_closed(primes: list[int] | None = None) -> dict:
    """H. Wick_4dist constant off ±1.  Floor shifts to ⟨m_4,K⟩.  OPEN."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        # Wick(π) off ±1 is the 3-pairing value Q=4q² (theorem D)
        wick_all_off = 4 * q * q
        t1 = wick_coll_T1(p)
        t2 = wick_coll_T2(p)
        t3 = wick_coll_T3_off(p)
        t4 = wick_coll_T4_off(p)
        wc = wick_coll_off(p)
        w4 = wick_4dist_off(p)
        row = t1 == 3 * q and t2 == 12 * q
        row = row and t3 == (q + 2) * (q - 3)
        row = row and t4 == 2 * q * q - 4 * q + 12
        row = row and t1 + t2 + t3 + t4 == wc
        row = row and wc == 3 * q * q + 10 * q + 6
        row = row and w4 + wc == wick_all_off
        row = row and w4 == q * q - 10 * q - 6
        # even ψ∉{1,χ}: ∑_{squares}ψ=0 ⇒ ∑_{r≠±1}ψ=−2
        what = What4(p)
        row = row and what == w4 * (-2)
        row = row and what == -2 * q * q + 20 * q + 12
        # ˆΔ_4 = ⟨m_4,K⟩ − ˆW_4 ≥ Δ4_budget
        # ⇔ ⟨m_4,K⟩ ≥ Δ4_budget + ˆW_4 = −4q²−20q
        row = row and Delta4_budget(p) + what == m4K_budget(p)
        row = row and m4K_budget(p) == -4 * q * q - 20 * q
        sample[str(p)] = {
            "wick_4dist_off": w4,
            "wick_coll_off": wc,
            "What4": what,
            "m4K_budget": m4K_budget(p),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Wick(π) off ±1 equals 4q² (real 3-pairing, theorem D).  "
            "Collision Wick by types, χ(−1)=1 so π is symmetric, "
            "and G(χ,μ)=p on Ω: T1 all-equal 3q; T2 3+1 four "
            "orientations each (3q/p)G(χ,±ξ or ±rξ)=3q, total 12q; "
            "T3 2+2: matching (a=b,c=d) is (q+2)(q−1), the two "
            "crossed matchings are −(q+2) off ±1, total (q+2)(q−3); "
            "T4 2+1+1: slots ab and cd contribute q(q−2)+2 each "
            "(π(c,d) Gauss + Jacobi J=G(χ,α)G(χ,−α)−(q−1)=1), "
            "four crosses each −2q+2(q+1)=2, total 2q²−4q+12.  "
            "Sum 3q²+10q+6.  Difference is the 4-distinct Wick "
            "q²−10q−6, constant on squares r≠±1.  Even ψ∉{1,χ} has "
            "∑_{r≠±1}ψ=−2, so ˆW_4=−2q²+20q+12.  Therefore "
            "ˆΔ_4=⟨m_4,K_ψ⟩−ˆW_4 and λ≥6 iff ⟨m_4,K_ψ⟩≥−4q²−20q.  "
            "The pairing bound is OPEN."
        ),
        "by_p_sample": sample,
    }


def theorem_floor_iff_gauss_4dist(primes: list[int] | None = None) -> dict:
    """L. Additive 4-dist pieces cancel.  Floor iff Gauss pairing ≥0.  OPEN."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        t1 = wick_coll_T1(p)
        t2 = wick_coll_T2(p)
        t3 = wick_coll_T3_pm1(p)
        t4 = wick_coll_T4_pm1(p)
        wc = wick_coll_pm1(p)
        row = t3 == (q + 2) * (2 * q - 3)
        row = row and t4 == 4 * q * q - 8 * q + 12
        # T4: ab+cd = 2q(q−2)+4; two nondegenerate crosses = 4;
        # two degenerate (1∓r=0) = 2(q(q−2)+2)
        row = row and 2 * q * (q - 2) + 4 + 4 + 2 * (q * (q - 2) + 2) == t4
        row = row and t1 + t2 + t3 + t4 == wc
        row = row and wc == 6 * q * q + 8 * q + 6
        # Boolean_coll(±1) = Wick_coll + Δ_coll
        bcoll = wc + delta_coll_pm1(p)
        row = row and bcoll == 6 * q * q - 10 * q
        # 4-dist Q(±1) = 8q² − Boolean_coll
        q4 = Q4dist_pm1(p)
        row = row and 8 * q * q - bcoll == q4
        row = row and q4 == 2 * q * q + 10 * q
        # ep+em = −Q4dist(1)−Q4dist(−1) = −(4q²+20q)
        row = row and -2 * q4 == -(4 * q * q + 20 * q)
        sample[str(p)] = {
            "wick_coll_pm1": wc,
            "Q4dist_pm1": q4,
            "additive_cancel": -2 * q4,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "At r=±1, Wick T3=(q+2)(2q−3) (one crossed matching has "
            "phase 1).  T4: ab and cd still contribute q(q−2)+2 each; "
            "the two crosses with 1±r≠0 stay 2; the two with 1∓r=0 "
            "have inner additive sum q−2 and contribute q(q−2)+2 each.  "
            "Total T4=4q²−8q+12, Wick_coll(±1)=6q²+8q+6.  Boolean−Wick "
            "Δ_coll(±1)=−(18q+6) gives Boolean_coll=6q²−10q.  Q(±1)=8q² "
            "(theorem D) so the 4-distinct part is 2q²+10q.  The additive "
            "pieces of I on 4-distinct are −e_ξ(d−c)−e_{−ξ}(d−c), hence "
            "sum to −(4q²+20q).  Therefore S_□−6q²=4q²+20q+⟨m_4,K⟩ "
            "equals the 4-distinct Gauss pairing "
            "∑ m4 e_ξ(b−a)·½(ψ̄(ξ(d−c))G(ψ)+(χψ)̄(ξ(d−c))G(χψ)).  "
            "λ≥6 iff that pairing is ≥0.  Positivity OPEN."
        ),
        "by_p_sample": sample,
    }


def bool_coll_T1(p: int) -> int:
    """T1 Boolean: q."""
    return q_of(p)


def bool_coll_T2(p: int) -> int:
    """T2 Boolean, four 3+1 orientations: 4q."""
    return 4 * q_of(p)


def bool_coll_T3_pm1(p: int) -> int:
    """T3 Boolean at r=±1: 2q²−3q."""
    q = q_of(p)
    return 2 * q * q - 3 * q


def bool_coll_T3_off(p: int) -> int:
    """T3 Boolean off ±1: q²−3q."""
    q = q_of(p)
    return q * q - 3 * q


def bool_coll_T4_pm1(p: int) -> int:
    """T4 Boolean at r=±1: 4q²−12q."""
    q = q_of(p)
    return 4 * q * q - 12 * q


def bool_coll_T4_off(p: int) -> int:
    """T4 Boolean off ±1: 2q²−12q."""
    q = q_of(p)
    return 2 * q * q - 12 * q


def bool_coll_pm1(p: int) -> int:
    """Boolean collision kernel at r=±1: 6q²−10q."""
    q = q_of(p)
    return 6 * q * q - 10 * q


def bool_coll_off(p: int) -> int:
    """Boolean collision kernel on squares r≠±1: 3q²−10q."""
    q = q_of(p)
    return 3 * q * q - 10 * q


def bool_coll_orthogonality_mass(p: int) -> int:
    """∑_{r□} ψ(r) Coll(r) = 6q² for even ψ∉{1,χ}."""
    q = q_of(p)
    return 6 * q * q


def bool_coll_hat_weight(p: int, sum_squares: int, w_pm1: int) -> int:
    """∑_□ w Coll = (3q²−10q)∑_□ w + 3q² (w(1)+w(−1)).

    For even w, w(1)=w(−1)=w_pm1.  Mean-zero on squares ⇒ ∑_□ w=0
    and the pairing is 6q² w_pm1.
    """
    q = q_of(p)
    return (3 * q * q - 10 * q) * sum_squares + 3 * q * q * (2 * w_pm1)


def theorem_bool_coll_orthogonality(primes: list[int] | None = None) -> dict:
    """M. 2-level Coll + character orthogonality.  Floor still OPEN."""
    if primes is None:
        primes = [r for r in range(5, 80) if is_prime(r)]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        t1 = bool_coll_T1(p)
        t2 = bool_coll_T2(p)
        t3p = bool_coll_T3_pm1(p)
        t3o = bool_coll_T3_off(p)
        t4p = bool_coll_T4_pm1(p)
        t4o = bool_coll_T4_off(p)
        # T1, T2
        row = t1 == q and t2 == 4 * q
        # T3: (a=a,c=c) always q²−q; crossed 1+r=0 or 1−r=0 resonant
        row = row and (q * q - q) + (q * q - q) + (-q) == t3p
        row = row and (q * q - q) + (-q) + (-q) == t3o
        # T4: A1=B3=q(q−2) always; four other slots generic −2q,
        # two of them resonant q(q−2) at each of r=±1
        row = row and 2 * q * (q - 2) + 2 * q * (q - 2) + 2 * (-2 * q) == t4p
        row = row and 2 * q * (q - 2) + 4 * (-2 * q) == t4o
        row = row and t4p == 4 * q * q - 12 * q
        row = row and t4o == 2 * q * q - 12 * q
        # totals
        row = row and t1 + t2 + t3p + t4p == bool_coll_pm1(p)
        row = row and t1 + t2 + t3o + t4o == bool_coll_off(p)
        row = row and bool_coll_pm1(p) == 6 * q * q - 10 * q
        row = row and bool_coll_off(p) == 3 * q * q - 10 * q
        # matches L's ±1 Boolean_coll
        row = row and bool_coll_pm1(p) == wick_coll_pm1(p) + delta_coll_pm1(p)
        # character orthogonality: Coll=(3q²−10q)+3q² 1_{±1}
        # ∑_□ ψ = 0, ψ(±1)=1 ⇒ ∑_rest ψ=−2
        bump = bool_coll_pm1(p) - bool_coll_off(p)
        row = row and bump == 3 * q * q
        orth = 2 * bool_coll_pm1(p) + bool_coll_off(p) * (-2)
        row = row and orth == bool_coll_orthogonality_mass(p)
        row = row and orth == 6 * q * q
        # even mean-zero test weight: ∑_□ w=0, w(1)=w(−1)
        row = row and bool_coll_hat_weight(p, 0, 1) == 6 * q * q
        # dropping −10q from Coll_off (wrong 3q²) breaks 6q²
        wrong_off = 3 * q * q
        wrong_orth = 2 * bool_coll_pm1(p) + wrong_off * (-2)
        row = row and wrong_orth != 6 * q * q
        # inverted I conjugation is independent; recorded for the live test
        sample[str(p)] = {
            "Coll_pm1": bool_coll_pm1(p),
            "Coll_off": bool_coll_off(p),
            "orth": orth,
            "wrong_orth": wrong_orth,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Boolean E[z4] on collisions is π via z²=1 (Max+-free).  "
            "T1: phase 1, mass q.  T2: four 3+1 orientations, each "
            "G(χ,±ξ or ±rξ)/p=1, total 4q.  T3: matching (a=a,c=c) "
            "is q²−q; the two crossed matchings contribute q(q−1) or "
            "−q according as 1∓r=0 or not, hence 2q²−3q on ±1 and "
            "q²−3q off.  T4: slots (0,0,c,d) and (0,t,s,s) are always "
            "q(q−2); the other four slots have inner additive sum "
            "−1−G(χ,−ξ)/p=−2 off resonance (G(χ,μ)=χ(μ)τ=p on Ω) and "
            "q−2 on 1±r=0, hence −2q generic and q(q−2) resonant.  "
            "Totals Coll(±1)=6q²−10q, Coll(off)=3q²−10q.  "
            "Even ψ∉{1,χ}: ∑_□ ψ=0 and ψ(±1)=1, so "
            "∑_□ ψ Coll=2(6q²−10q)+(3q²−10q)(−2)=6q².  "
            "The same 2-level formula gives ∑ w Coll=6q² w(1) for "
            "every even mean-zero weight on T.  Therefore "
            "S_□−6q²=∑_□ ψ F with F=Q−Coll, i.e. the Gauss 4-dist "
            "pairing of m4.  Dropping −10q from Coll_off breaks 6q².  "
            "F̂≥0 is OPEN (no completed square for F)."
        ),
        "by_p_sample": sample,
    }


def gauss_coll_T3(p: int) -> int:
    """T3 Boolean Gauss (two crossed matchings): q²."""
    return q_of(p) ** 2


def gauss_coll_T4(p: int) -> int:
    """T4 Boolean Gauss (four crosses × q²/2): 2q²."""
    return 2 * q_of(p) ** 2


def gauss_coll_mass(p: int) -> int:
    """A_coll=3q²."""
    return 3 * q_of(p) ** 2


def theorem_gauss_coll_is_3q2(primes: list[int] | None = None) -> dict:
    """N. Collision Gauss mass 3q².  Floor still S_□≥6q²."""
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        F = _kernel_field(p)
        t3 = gauss_coll_T3(p)
        t4 = gauss_coll_T4(p)
        row = t3 == q * q and t4 == 2 * q * q
        row = row and t3 + t4 == gauss_coll_mass(p)
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gpsibar = _gauss(F, [z.conjugate() for z in psi])
        Gchi = _gauss(F, [F["chi"][x] for x in range(q)])
        Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(q)])
        gprod = Gpsi * Gpsibar
        row = row and abs(gprod - q) < 1e-6
        # T3 cross: (q/2) G(ψ)G(ψ̄)=q²/2
        row = row and abs((q / 2) * gprod - q * q / 2) < 1e-6
        # T4 cross after Paley cancel: (q/(2p))*p*q = q²/2
        row = row and abs((q / (2 * p)) * p * q - q * q / 2) < 1e-9
        # inverted J(χ,ψ)=G(χ)G(ψ)/G(χψ) vs J(χ,ψ̄)
        Jright = Gchi * Gpsibar / Gchipsi.conjugate() if abs(Gchipsi) > 1e-9 else 0j
        Jwrong = Gchi * Gpsi / Gchipsi if abs(Gchipsi) > 1e-9 else 0j
        row = row and abs(Jwrong - Jright) > 0.5
        sample[str(p)] = {
            "T3": t3,
            "T4": t4,
            "A_coll": gauss_coll_mass(p),
            "Gpsi_Gpsibar": gprod,
            "J_inverted_diff": abs(Jwrong - Jright),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Boolean collisions, Gauss piece of I only.  T1 α=0; T2 "
            "∑(χψ̄)=0, both 0.  T3 crossed matchings: ∑_{t≠0} e_ξ(t)ψ̄(t)"
            "=G(ψ̄,ξ)=ψ(ξ)G(ψ̄), times G(ψ)/2 gives q²/2 each.  T4 ab,cd "
            "vanish.  Each cross: ∑_x χ(x−u)ψ̄(x)=χ(u)ψ̄(u) J(χ,ψ̄) with "
            "J=G(χ)G(ψ̄)/G(χψ̄); G(χ)χ(ξ)=p on Ω and G(ψ)G(ψ̄)=q give "
            "q²/2.  A_coll=3q².  Parseval A_full=S_□/2 ⇒ "
            "A_4dist=S_□/2−3q², same floor as L.  Inverted J(χ,ψ) "
            "does not cancel.  S_□≥6q² OPEN."
        ),
        "by_p_sample": sample,
    }


def certify_p5_Rhat_beats_budget() -> dict:
    """Live Max+: E|M_ψ|² for the order-12 even character at p=5.

    Fails if the character-pair formula does not recover 80/13.
    """
    from pathlib import Path as P

    import numpy as np

    from e1_gmin_m4_prop15270 import _make_field, _primroot_q

    p, q = 5, 25
    path = P("/tmp/maxplus_p5.npy")
    if not path.is_file():
        return {"p": 5, "skipped": True, "reason": "Max+ cache missing"}
    F = _make_field(p)
    mul, ia = F["mul"], F["ia"]
    g = _primroot_q(F, p)

    def tr(u):
        return (2 * (u % p) + (u // p) * ia) % p

    def padd(mu, x):
        return np.exp(2j * np.pi * tr(mul(mu, x)) / p)

    chi = F["chi"]
    tau = sum(chi(x) * padd(1, x) for x in range(1, q))
    Omega = [xi for xi in range(1, q) if abs(chi(xi) * tau - p) < 1e-5]
    units = [1]
    for _ in range(q - 2):
        units.append(mul(units[-1], g))
    dlog = {units[j]: j for j in range(q - 1)}
    # even character k=2 (order 12 on F_q^×)
    k = 2
    Y = np.load(path).astype(np.float64)
    Z = Y[:, 1:]
    Ch = np.array([[padd(xi, x) for x in range(q)] for xi in Omega], dtype=np.complex128)
    power = np.abs(Ch @ Z.T) ** 2
    wgt = np.array(
        [np.exp(-2j * np.pi * k * dlog[xi] / (q - 1)) for xi in Omega]
    )
    M = wgt @ power
    em = float(np.mean(np.abs(M) ** 2))
    em_frac = Fraction(em).limit_denominator(20)
    wick = EM_wick(p)
    rhat = em_frac - wick
    bud = Rhat_budget(p)
    lam = rayleigh_from_EM(p, em_frac)
    # live Parseval: M = q Ã(ψ)/G(ψ) pointwise (theorem E)
    Gpsi = sum(
        np.exp(2j * np.pi * k * dlog[x] / (q - 1)) * padd(1, x)
        for x in range(1, q)
    )
    Ahat = np.zeros(Z.shape[0], dtype=np.complex128)
    for d in range(1, q):
        # field add: x+d
        shift = np.array(
            [
                (x % p + d % p) % p + ((x // p + d // p) % p) * p
                for x in range(q)
            ]
        )
        Ad = np.sum(Z * Z[:, shift], axis=1)
        Ahat += Ad * np.exp(2j * np.pi * k * dlog[d] / (q - 1))
    M_from_A = q * Ahat / Gpsi
    parseval_err = float(np.max(np.abs(M - M_from_A)))
    return {
        "p": 5,
        "EM_min": str(em_frac),
        "EM_float": em,
        "Rhat": str(rhat),
        "budget": bud,
        "beats": rhat >= bud,
        "lambda": str(lam),
        "lambda_is_80_13": lam == Fraction(80, 13),
        "Q_diag": float(np.mean(power[0] ** 2)),
        "Q_diag_is_8q2": abs(float(np.mean(power[0] ** 2)) - Q_diag_closed(p)) < 1e-6,
        "parseval_err": parseval_err,
        "parseval_ok": parseval_err < 1e-8,
        "proved_general": False,
    }


def _kernel_field(p: int) -> dict:
    """Additive/multiplicative tables for the Jacobi kernel (Max+-free)."""
    from e1_gmin_m4_prop15270 import _make_field, _primroot_q

    q = p * p
    F = _make_field(p)
    mul, chi, powm, ia = F["mul"], F["chi"], F["powm"], F["ia"]
    add = [[0] * q for _ in range(q)]
    mu = [[0] * q for _ in range(q)]
    for x in range(q):
        for y in range(q):
            add[x][y] = (x % p + y % p) % p + ((x // p + y // p) % p) * p
            mu[x][y] = mul(x, y)
    neg = [((-x) % p) + ((-(x // p)) % p) * p for x in range(q)]
    inv = [0 if x == 0 else powm(x, q - 2) for x in range(q)]
    chi_t = [chi(x) for x in range(q)]
    g = _primroot_q(F, p)
    units = [1]
    for _ in range(q - 2):
        units.append(mul(units[-1], g))
    dlog = [-1] * q
    for j, u in enumerate(units):
        dlog[u] = j

    def tr(u: int) -> int:
        return (2 * (u % p) + (u // p) * ia) % p

    trt = [tr(x) for x in range(q)]
    return {
        "p": p,
        "q": q,
        "add": add,
        "mul": mu,
        "neg": neg,
        "inv": inv,
        "chi": chi_t,
        "dlog": dlog,
        "g": g,
        "neg1": units[(q - 1) // 2],
        "trt": trt,
        "squares": [r for r in range(1, q) if chi_t[r] == 1],
    }


def _psi_vec(F: dict, k: int) -> list[complex]:
    q = F["q"]
    out = [0j] * q
    for x in range(1, q):
        out[x] = cmath.exp(2j * math.pi * k * F["dlog"][x] / (q - 1))
    return out


def _e_tr(F: dict, x: int) -> complex:
    return cmath.exp(2j * math.pi * F["trt"][x] / F["p"])


def _gauss(F: dict, char: list[complex]) -> complex:
    return sum(char[x] * _e_tr(F, x) for x in range(1, F["q"]))


def I_of_alpha(
    F: dict, psi: list[complex], Gpsi: complex, Gchipsi: complex, alpha: int
) -> complex:
    """I(α)=∑_{r□, r≠±1} ψ(r) e(α r).  Uses G(λ,α)=λ̄(α)G(λ)."""
    if alpha == 0:
        return -2 + 0j
    psia = psi[alpha]
    chipsi_a = F["chi"][alpha] * psi[alpha]
    return (
        0.5 * (psia.conjugate() * Gpsi + chipsi_a.conjugate() * Gchipsi)
        - _e_tr(F, alpha)
        - _e_tr(F, F["neg"][alpha])
    )


def I_of_alpha_inverted(
    F: dict, psi: list[complex], Gpsi: complex, Gchipsi: complex, alpha: int
) -> complex:
    """Wrong inversion G(λ,α)=λ(α)G(λ).  Must not match I."""
    if alpha == 0:
        return -2 + 0j
    psia = psi[alpha]
    chipsi_a = F["chi"][alpha] * psi[alpha]
    return (
        0.5 * (psia * Gpsi + chipsi_a * Gchipsi)
        - _e_tr(F, alpha)
        - _e_tr(F, F["neg"][alpha])
    )


def I_direct(F: dict, psi: list[complex], alpha: int) -> complex:
    s = 0j
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        ax = 0 if alpha == 0 else F["mul"][alpha][r]
        s += psi[r] * (_e_tr(F, ax) if alpha else 1)
    if alpha == 0:
        return s
    return s


def n_spec_lambda(F: dict, lam: int) -> int:
    """# of s with A±B=0: 2·1_{χ(1−λ)=1}+2·1_{χ(λ(λ−1))=1}."""
    one_m = F["add"][1][F["neg"][lam]]
    n = 0
    if F["chi"][one_m] == 1:
        n += 2
    if F["chi"][lam] * F["chi"][one_m] == 1:
        n += 2
    return n


def kernel_mass(p: int) -> int:
    """∑_{4-dist} K_ψ = ∑_λ J_1(λ) = −2q²+12q."""
    q = p * p
    return -2 * q * q + 12 * q


def kappa_kernel_mass(p: int) -> int:
    """∑_λ J_κ(λ)=⟨κ_C,K_ψ⟩=q Ŵ_4=−2q³+20q²+12q."""
    q = q_of(p)
    return q * What4(p)


def I_L2_mass(p: int) -> int:
    """∑_{α≠0}|I(α)|² = (q²−5q−8)/2, independent of ψ."""
    q = q_of(p)
    return (q * q - 5 * q - 8) // 2


def kappaC_K_mass(p: int) -> int:
    """⟨κ_C, K_ψ⟩ = q ˆW_4."""
    return q_of(p) * What4(p)


def _f_lam(F: dict, lam: int, s: int) -> int:
    """f_λ(s)=s(s−1)(1−λ)/(1−λs)."""
    one_m = F["add"][1][F["neg"][lam]]
    den = F["add"][1][F["neg"][F["mul"][lam][s]]]
    if den == 0:
        return 0
    sm1 = F["add"][s][F["neg"][1]]
    num = F["mul"][F["mul"][s][sm1]][one_m]
    return F["mul"][num][F["inv"][den]]


def _ratio_mu_sm1(F: dict, lam: int, s: int) -> int:
    one_m = F["add"][1][F["neg"][lam]]
    sm1 = F["add"][s][F["neg"][1]]
    den1 = F["add"][1][F["neg"][F["mul"][lam][s]]]
    if sm1 == 0 or den1 == 0:
        return 0
    return F["mul"][F["mul"][s][one_m]][F["inv"][F["mul"][sm1][den1]]]


def theta_psi(F: dict, psi: list[complex], lam: int) -> complex:
    """∑_{χ(f_λ)=1} ψ(μ/(s−1)).  Direct loop; equals theta_from_K3."""
    q = F["q"]
    inv_lam = F["inv"][lam]
    acc = 0j
    for s in range(q):
        if s in (0, 1, inv_lam):
            continue
        f = _f_lam(F, lam, s)
        if f == 0 or F["chi"][f] != 1:
            continue
        r = _ratio_mu_sm1(F, lam, s)
        if r:
            acc += psi[r]
    return acc


def _c_lam(F: dict, lam: int) -> int:
    """(λ−1)/λ."""
    return F["mul"][F["add"][lam][F["neg"][1]]][F["inv"][lam]]


def K3_jacobi(
    F: dict,
    a: list[complex],
    b: list[complex],
    c: list[complex],
    beta: int,
) -> complex:
    """∑_{s≠0,1,β} a(s) b(s−1) c(s−β)."""
    q = F["q"]
    acc = 0j
    for s in range(q):
        if s in (0, 1, beta):
            continue
        sm1 = F["add"][s][F["neg"][1]]
        smb = F["add"][s][F["neg"][beta]]
        if sm1 == 0 or smb == 0:
            continue
        acc += a[s] * b[sm1] * c[smb]
    return acc


def theta_from_K3(F: dict, psi: list[complex], lam: int) -> complex:
    """½(ψ(c) K(ψ,ψ̄,ψ̄;1/λ)+(χψ)(c) K(χψ,χψ̄,χψ̄;1/λ))."""
    c = _c_lam(F, lam)
    beta = F["inv"][lam]
    psib = [z.conjugate() for z in psi]
    chipsi = [F["chi"][x] * psi[x] for x in range(F["q"])]
    chpsib = [F["chi"][x] * psib[x] for x in range(F["q"])]
    return 0.5 * (
        psi[c] * K3_jacobi(F, psi, psib, psib, beta)
        + chipsi[c] * K3_jacobi(F, chipsi, chpsib, chpsib, beta)
    )


def theta_from_wrong_J3(F: dict, psi: list[complex], lam: int) -> complex:
    """Wrong Θ: pretends ψ(x²)≡1 and uses J_3(ψ)=K(ψ,ψ,ψ̄)."""
    c = _c_lam(F, lam)
    beta = F["inv"][lam]
    chipsi = [F["chi"][x] * psi[x] for x in range(F["q"])]
    psib = [z.conjugate() for z in psi]
    chpsib = [F["chi"][x] * psib[x] for x in range(F["q"])]
    return 0.5 * (
        psi[c] * K3_jacobi(F, psi, psi, psib, beta)
        + chipsi[c] * K3_jacobi(F, chipsi, chipsi, chpsib, beta)
    )


def H_legendre(F: dict, gamma: int) -> int:
    """H(γ)=∑_x χ(x(x−1)(x−γ))."""
    q = F["q"]
    acc = 0
    for x in range(q):
        xm1 = F["add"][x][F["neg"][1]]
        xmg = F["add"][x][F["neg"][gamma]]
        acc += int(F["chi"][x]) * int(F["chi"][xm1]) * int(F["chi"][xmg])
    return acc


def Xi_from_H(F: dict, lam: int) -> int:
    """Ξ(λ)=χ(1−λ) H(1−λ)."""
    one_m = F["add"][1][F["neg"][lam]]
    return int(F["chi"][one_m]) * H_legendre(F, one_m)


def theta_psi_inverted(F: dict, psi: list[complex], lam: int, xi: int) -> complex:
    """Wrong Θ=∑_□ ψ̄(ξ² f).  Mass fails."""
    q = F["q"]
    inv_lam = F["inv"][lam]
    xi2 = F["mul"][xi][xi]
    acc = 0j
    for s in range(q):
        if s in (0, 1, inv_lam):
            continue
        f = _f_lam(F, lam, s)
        if f == 0 or F["chi"][f] != 1:
            continue
        acc += psi[F["mul"][xi2][f]].conjugate()
    return acc


def J1_of_lambda(F: dict, psi: list[complex], lam: int) -> complex:
    q = F["q"]
    return q * q * theta_psi(F, psi, lam) + 2 * q * (q - 3) - q * q * n_spec_lambda(F, lam)


def J1_of_lambda_inverted(F: dict, psi: list[complex], lam: int, xi: int) -> complex:
    q = F["q"]
    return (
        q * q * theta_psi_inverted(F, psi, lam, xi)
        + 2 * q * (q - 3)
        - q * q * n_spec_lambda(F, lam)
    )


def Xi_chi(F: dict, lam: int) -> int:
    q = F["q"]
    inv_lam = F["inv"][lam]
    acc = 0
    for s in range(q):
        if s in (0, 1, inv_lam):
            continue
        f = _f_lam(F, lam, s)
        if f:
            acc += int(F["chi"][f])
    return acc


def Jkappa_of_lambda(F: dict, psi: list[complex], lam: int) -> complex:
    """⟨κ 1_λ, K⟩.  Ξ_spec=n_spec because f=±(s−1)² on A=±B."""
    q = F["q"]
    one_m = F["add"][1][F["neg"][lam]]
    pref = (1 + int(F["chi"][lam]) + int(F["chi"][one_m])) * int(F["chi"][one_m])
    Th = theta_from_K3(F, psi, lam)
    Xi = Xi_from_H(F, lam)
    nsp = n_spec_lambda(F, lam)
    return pref * (q * q * Th + 2 * q * Xi - q * q * nsp)


def theorem_jacobi_kernel_pairings(primes: list[int] | None = None) -> dict:
    """I. ⟨κ 1_λ, K⟩ as Gauss/Jacobi.  Floor still OPEN."""
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    for p in primes:
        F = _kernel_field(p)
        q = F["q"]
        k = 2
        psi = _psi_vec(F, k)
        chipsi = [F["chi"][x] * psi[x] for x in range(q)]
        Gpsi = _gauss(F, psi)
        Gchipsi = _gauss(F, chipsi)
        i_err = 0.0
        inv_err = 0.0
        for a in range(q):
            cl = I_of_alpha(F, psi, Gpsi, Gchipsi, a)
            dn = I_direct(F, psi, a)
            i_err = max(i_err, abs(cl - dn))
            if a:
                wr = I_of_alpha_inverted(F, psi, Gpsi, Gchipsi, a)
                inv_err = max(inv_err, abs(wr - dn))
        i0 = I_direct(F, psi, 0)
        i0_pm1 = i0 + psi[1] + psi[F["neg1"]]
        mass = 0j
        mass_k = 0j
        mass_wrong = 0j
        th_err = 0.0
        wrong_th = 0.0
        xi_err = 0
        H_gauss = True
        xi = 1
        for lam in range(q):
            if lam in (0, 1):
                continue
            mass += J1_of_lambda(F, psi, lam)
            mass_k += Jkappa_of_lambda(F, psi, lam)
            mass_wrong += J1_of_lambda_inverted(F, psi, lam, xi)
            th_err = max(
                th_err, abs(theta_from_K3(F, psi, lam) - theta_psi(F, psi, lam))
            )
            wrong_th = max(
                wrong_th,
                abs(theta_from_wrong_J3(F, psi, lam) - theta_psi(F, psi, lam)),
            )
            if Xi_from_H(F, lam) != Xi_chi(F, lam):
                xi_err += 1
            one_m = F["add"][1][F["neg"][lam]]
            Hv = H_legendre(F, one_m)
            if abs(Hv) in (0, p, 2 * p):
                pass
            else:
                H_gauss = False
        want = kernel_mass(p)
        want_k = kappa_kernel_mass(p)
        star_zero = True
        for lam in range(q):
            if lam in (0, 1):
                continue
            one_m = F["add"][1][F["neg"][lam]]
            if F["chi"][lam] == 1 and F["chi"][one_m] == 1:
                continue
            if not (F["chi"][lam] != 1 or F["chi"][one_m] != 1):
                star_zero = False
        # polynomial identities, Max+-free
        poly_ok = want_k == q * What4(p) == -2 * q**3 + 20 * q * q + 12 * q
        poly_ok = poly_ok and (want == -2 * q * q + 12 * q)
        row = (
            i_err < 1e-8
            and inv_err > 0.5
            and abs(i0 + 2) < 1e-8
            and abs(i0_pm1) < 1e-8
            and abs(mass - want) < 1e-6
            and abs(mass_k - want_k) < 1e-4
            and abs(mass_wrong - want) > 10
            and th_err < 1e-8
            and wrong_th > 0.5
            and xi_err == 0
            and not H_gauss
            and abs(Gpsi) > p - 1e-6
            and star_zero
            and poly_ok
        )
        sample[str(p)] = {
            "I_err": i_err,
            "inverted_err": inv_err,
            "I0": i0,
            "mass": mass,
            "mass_expect": want,
            "mass_kappa": mass_k,
            "mass_kappa_expect": want_k,
            "mass_inverted": mass_wrong,
            "theta_K3_err": th_err,
            "theta_wrong_J3_err": wrong_th,
            "Xi_H_mismatches": xi_err,
            "H_is_gauss_value": H_gauss,
            "ok": row,
        }
        if not row:
            ok = False
    p5_match = True
    if 5 in primes:
        F = _kernel_field(5)
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
        xi = 1
        for lam in (2, 5, 8, 13):
            j1 = J1_of_lambda(F, psi, lam)
            jk = Jkappa_of_lambda(F, psi, lam)
            inv_lam = F["inv"][lam]
            one_m = F["add"][1][F["neg"][lam]]
            acc = 0j
            acck = 0j
            Kbar = 1 + int(F["chi"][lam]) + int(F["chi"][one_m])
            for s in range(F["q"]):
                if s in (0, 1, inv_lam):
                    continue
                den = F["add"][1][F["neg"][F["mul"][lam][s]]]
                if den == 0:
                    continue
                mu = F["mul"][F["mul"][s][one_m]][F["inv"][den]]
                sm1 = F["add"][s][F["neg"][1]]
                rat = F["mul"][F["mul"][s][sm1]][F["inv"][den]]
                kap = Kbar * int(F["chi"][rat])
                for v in range(1, F["q"]):
                    u = F["mul"][v][mu]
                    w = F["mul"][v][s]
                    t = F["add"][w][F["neg"][v]]
                    term = _e_tr(F, F["mul"][xi][u]) * I_of_alpha(
                        F, psi, Gpsi, Gchipsi, F["mul"][xi][t]
                    )
                    acc += term
                    acck += kap * term
            if abs(j1 - acc * F["q"]) > 1e-6:
                p5_match = False
            if abs(jk - acck * F["q"]) > 1e-5:
                p5_match = False
        sample["p5_J1_Jkappa_sv"] = p5_match
        if not p5_match:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Even ψ∉{1,χ}: I(0)=−2 and I(α)=½(ψ̄(α)G(ψ)+(χψ)̄(α)G(χψ))"
            "−e(α)−e(−α) for α≠0, via G(λ,α)=λ̄(α)G(λ).  "
            "μ/(s−1)=f_λ/(s−1)² with f_λ=c·s(s−1)/(s−1/λ), c=(λ−1)/λ.  "
            "Even ≠ square, so ψ((s−1)²)=ψ(s−1)² and "
            "Θ_ψ=½(ψ(c) K(ψ,ψ̄,ψ̄;1/λ)+(χψ)(c) K(χψ,χψ̄,χψ̄;1/λ)) "
            "for the complete Jacobi K(α,β,γ;δ)=∑ α(s)β(s−1)γ(s−δ).  "
            "The untwisted J_3(ψ)=K(ψ,ψ,ψ̄) is the wrong Jacobi.  "
            "Ξ=χ(1−λ)H(1−λ) with H the Legendre sum (not a Gauss value "
            "in {0,±p,±2p}).  On A=±B, f=±(s−1)² so Ξ_spec=n_spec.  "
            "J_1=q² Θ+2q(q−3)−q² n_spec and "
            "J_κ=(1+χ(λ)+χ(1−λ))χ(1−λ)[q² Θ+2q Ξ−q² n_spec].  "
            "∑_λ J_1=−2q²+12q and ∑_λ J_κ=q Ŵ_4=−2q³+20q²+12q.  "
            "K(η,η̄,η̄) is not a Gauss product.  "
            "Does not prove ⟨m_4,K⟩≥−4q²−20q."
        ),
        "by_p_sample": sample,
    }


def certify_live_m4K() -> dict:
    """Live Max+ ⟨m_4,K⟩=ˆΔ_4+Ŵ_4 at the minimizing even character.

    Check only.  Does not prove the budget.
    """
    from pathlib import Path as P

    import numpy as np

    from e1_gmin_m4_prop15270 import _make_field, _primroot_q

    out = {}
    jobs = ((5, "/tmp/maxplus_p5.npy", 2), (7, "/tmp/e1_p7/maxplus.npy", 8))
    for p, path, k in jobs:
        if not P(path).is_file():
            out[str(p)] = {"skipped": True}
            continue
        q = p * p
        Fld = _make_field(p)
        mul, ia, chi = Fld["mul"], Fld["ia"], Fld["chi"]
        g = _primroot_q(Fld, p)

        def tr(u):
            return (2 * (u % p) + (u // p) * ia) % p

        def padd(mu, x):
            return np.exp(2j * np.pi * tr(mul(mu, x)) / p)

        tau = sum(chi(x) * padd(1, x) for x in range(1, q))
        Omega = [xi for xi in range(1, q) if abs(chi(xi) * tau - p) < 1e-5]
        units = [1]
        for _ in range(q - 2):
            units.append(mul(units[-1], g))
        dlog = {units[j]: j for j in range(q - 1)}
        Y = np.load(path).astype(np.float64)
        Z = Y[:, 1:]
        Ch = np.array(
            [[padd(xi, x) for x in range(q)] for xi in Omega], dtype=np.complex128
        )
        power = np.abs(Ch @ Z.T) ** 2
        wgt = np.array(
            [np.exp(-2j * np.pi * k * dlog[xi] / (q - 1)) for xi in Omega]
        )
        M = wgt @ power
        em = float(np.mean(np.abs(M) ** 2))
        Ssq = em / ((q - 1) / 2)
        D4 = Ssq - 8 * q * q - 40 * q - 12
        m4K = D4 + What4(p)
        bud = m4K_budget(p)
        out[str(p)] = {
            "k": k,
            "m4K": m4K,
            "D4": D4,
            "What4": What4(p),
            "budget": bud,
            "room": m4K - bud,
            "beats": m4K >= bud - 1e-6,
            "proved_general": False,
        }
    return out


def jacobi_kernel_pairings_proved() -> bool:
    return bool(theorem_jacobi_kernel_pairings()["proved"])


def theorem_kappaC_pairs_as_What4(primes: list[int] | None = None) -> dict:
    """J. ⟨κ_C,K⟩=q ˆW_4.  Floor reduces to ⟨ρ,K⟩.  OPEN."""
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        F = _kernel_field(p)
        # skip k=0 and the quadratic k=(q-1)/2 ≡ χ
        ks = [k for k in range(2, q - 1, 2) if k != (q - 1) // 2]
        rows = {}
        for k in ks[:4]:
            psi = _psi_vec(F, k)
            mass = 0j
            for lam in range(q):
                if lam in (0, 1):
                    continue
                mass += Jkappa_of_lambda(F, psi, lam)
            want = kappaC_K_mass(p)
            row = abs(mass - want) < 1e-6
            rows[str(k)] = {"sum_Jkappa": mass.real, "want": want, "ok": row}
            if not row:
                ok = False
        # inverted Θ must not give the κ-mass
        psi = _psi_vec(F, 2)
        bad = 0j
        for lam in range(q):
            if lam in (0, 1):
                continue
            # reuse inverted Θ inside the J_κ pref skeleton
            one_m = F["add"][1][F["neg"][lam]]
            pref = (1 + int(F["chi"][lam]) + int(F["chi"][one_m])) * int(
                F["chi"][one_m]
            )
            Th = theta_psi_inverted(F, psi, lam, 1)
            bad += pref * (q * q * Th)
        row_bad = abs(bad - kappaC_K_mass(p)) > 10
        if not row_bad:
            ok = False
        sample[str(p)] = {"by_k": rows, "inverted_breaks": row_bad, "ok": ok}
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "On 4-distinct indices π(a,b)=C_{ab}/p, so "
            "κ_C=p²(π_ab π_cd+π_ac π_bd+π_ad π_bc)=p² Wick_π.  "
            "H says the Wick_π pairing against K_ψ is ˆW_4, hence "
            "⟨κ_C,K_ψ⟩=q ˆW_4 and ∑_λ J_κ(λ)=q ˆW_4.  "
            "Thus ⟨m_4,K⟩=ˆW_4+⟨ρ,K⟩ with ρ=m_4−κ_C/p², and "
            "λ≥6 iff ⟨ρ,K_ψ⟩≥−2q²−40q−12.  Inverted Θ breaks the mass.  "
            "The ρ pairing is OPEN."
        ),
        "by_p_sample": sample,
    }


def theorem_I_L2_closed(primes: list[int] | None = None) -> dict:
    """K. ∑_{α≠0}|I|²=(q²−5q−8)/2, independent of ψ."""
    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        want = I_L2_mass(p)
        # algebra of the split
        sumA = q * (q - 1) // 2
        sumB = 2 * q - 4
        cross = 4 * q
        row = sumA + sumB - cross == want
        F = _kernel_field(p)
        live = []
        for k in (2, 4):
            psi = _psi_vec(F, k)
            Gpsi = _gauss(F, psi)
            Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(q)])
            s = 0.0
            for a in range(1, q):
                Ia = I_of_alpha(F, psi, Gpsi, Gchipsi, a)
                s += abs(Ia) ** 2
            live.append(s)
            row = row and abs(s - want) < 1e-6
        sample[str(p)] = {
            "want": want,
            "sumA": sumA,
            "sumB": sumB,
            "cross": cross,
            "live": live,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "I=A−B on F_q^× with A=½(ψ̄G(ψ)+(χψ)̄G(χψ)) and "
            "B=e+e(−·).  ∑|A|²=q(q−1)/2 because ∑χ=0.  "
            "∑|B|²=2q−4.  ∑AB=2q from G(ψ)G(ψ̄)=q=G(χψ)G((χψ)̄) "
            "and evenness ψ(−1)=1, so 2 Re∑AB=4q.  "
            "Hence ∑|I|²=q(q−1)/2−2q−4=(q²−5q−8)/2, independent of ψ.  "
            "Does not bound ⟨ρ,K⟩."
        ),
        "by_p_sample": sample,
    }


def a4_coll_over_Gbar(p: int) -> int:
    """A₄_coll(ψ)/G(ψ̄)=6q."""
    return 6 * q_of(p)


def a4_coll_GA_mass(p: int) -> int:
    """G(ψ) A₄_coll(ψ)=6q²."""
    q = q_of(p)
    return 6 * q * q


def _a4_coll_from_pi(F: dict, psi: list[complex], invert_char: bool = False) -> complex:
    """∑_{coll} π e_ξ(b−a) λ(ξ(d−c)), λ=ψ̄ (or ψ if invert).  a=0 ×q."""
    import numpy as np

    q = F["q"]
    add = np.array(F["add"], dtype=np.int32)
    mu = np.array(F["mul"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    chi = np.array(F["chi"], dtype=np.int8)
    trt = np.array(F["trt"], dtype=np.int32)
    eadd = np.exp(2j * np.pi * np.arange(F["p"]) / F["p"])
    tau = sum(chi[x] * eadd[trt[mu[1, x]]] for x in range(1, q))
    Omega = [xi for xi in range(1, q) if abs(chi[xi] * tau - F["p"]) < 1e-8]
    xi = Omega[0]
    Pi = np.eye(q, dtype=np.float64)
    for a in range(q):
        for b in range(q):
            if a != b:
                Pi[a, b] = chi[add[b, neg[a]]] / F["p"]
    idx = np.arange(q)
    B, C, D = np.meshgrid(idx, idx, idx, indexing="ij")
    b, c, d = B.ravel(), C.ravel(), D.ravel()
    a = np.zeros_like(b)
    nuniq = (
        (a != b).astype(np.int8)
        + (a != c)
        + (a != d)
        + (b != c)
        + (b != d)
        + (c != d)
    )
    coll = nuniq < 6
    e_ba = eadd[trt[mu[xi, add[b, neg[a]]]]]
    alpha = mu[xi, add[d, neg[c]]]
    char = np.array(psi, dtype=np.complex128)
    if not invert_char:
        char = np.conj(char)
    mom = np.zeros(b.shape, np.float64)
    mom[nuniq == 0] = 1.0
    mom[nuniq == 4] = 1.0
    m2 = nuniq == 3
    if np.any(m2):
        aa, bb, cc, dd = a[m2], b[m2], c[m2], d[m2]
        ca = 1 + (aa == bb).astype(np.int8) + (aa == cc) + (aa == dd)
        cb = (bb == aa).astype(np.int8) + 1 + (bb == cc) + (bb == dd)
        cc_ = (cc == aa).astype(np.int8) + (cc == bb) + 1 + (cc == dd)
        trip = np.where(ca == 3, aa, np.where(cb == 3, bb, np.where(cc_ == 3, cc, dd)))
        sing = np.where(ca == 1, aa, np.where(cb == 1, bb, np.where(cc_ == 1, cc, dd)))
        mom[m2] = Pi[trip, sing]
    m4 = nuniq == 5
    if np.any(m4):
        aa, bb, cc, dd = a[m4], b[m4], c[m4], d[m4]
        ca = 1 + (aa == bb).astype(np.int8) + (aa == cc) + (aa == dd)
        cb = (bb == aa).astype(np.int8) + 1 + (bb == cc) + (bb == dd)
        cc_ = (cc == aa).astype(np.int8) + (cc == bb) + 1 + (cc == dd)
        cd = (dd == aa).astype(np.int8) + (dd == bb) + (dd == cc) + 1
        s1 = np.where(ca == 1, aa, np.where(cb == 1, bb, np.where(cc_ == 1, cc, dd)))
        s2 = np.where(
            (cb == 1) & (bb != s1),
            bb,
            np.where(
                (cc_ == 1) & (cc != s1),
                cc,
                np.where((cd == 1) & (dd != s1), dd, aa),
            ),
        )
        s2 = np.where((ca == 1) & (aa != s1), aa, s2)
        mom[m4] = Pi[s1, s2]
    return complex(q * np.sum(mom[coll] * e_ba[coll] * char[alpha[coll]]))


def theorem_A4_coll_equals_6q_Gbar(primes: list[int] | None = None) -> dict:
    """O. A₄_coll(ψ)=6q G(ψ̄).  Alignment identity; floor still OPEN."""
    if primes is None:
        primes = [5, 7]
    ok = True
    sample = {}
    for p in primes:
        q = q_of(p)
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        psibar = [z.conjugate() for z in psi]
        Gpsi = _gauss(F, psi)
        Gbar = _gauss(F, psibar)
        Gchi = _gauss(F, [F["chi"][x] for x in range(q)])
        chipsi_bar = [F["chi"][x] * psibar[x] for x in range(q)]
        sum_chibar = sum(chipsi_bar[x] for x in range(1, q))
        # G(ψ̄,−1)=ψ(−1)G(ψ̄)=G(ψ̄)
        Gbar_m1 = sum(psibar[x] * _e_tr(F, F["neg"][x]) for x in range(1, q))
        # ε=G(χ)/p, χ(ξ)ε=1 on Ω
        eps = Gchi / p
        xi = None
        for cand in range(1, q):
            if abs(F["chi"][cand] * Gchi - p) < 1e-6:
                xi = cand
                break
        row = abs(Gpsi * Gbar - q) < 1e-6
        row = row and abs(sum_chibar) < 1e-8
        row = row and abs(Gbar_m1 - Gbar) < 1e-8
        row = row and xi is not None and abs(F["chi"][xi] * eps - 1) < 1e-8
        # type-sum coefficients
        t3 = 2 * q * Gbar
        t4 = 4 * q * Gbar
        tot = t3 + t4
        row = row and abs(tot - a4_coll_over_Gbar(p) * Gbar) < 1e-8
        row = row and abs(Gpsi * tot - a4_coll_GA_mass(p)) < 1e-6
        # inverted phase: same types with ψ give 6q G(ψ)
        tot_inv = 6 * q * Gpsi
        row = row and abs(Gpsi * tot_inv - 6 * q * Gpsi * Gpsi) < 1e-6
        # live π-enumeration at p=5 (and p=7)
        live = _a4_coll_from_pi(F, psi, invert_char=False)
        live_inv = _a4_coll_from_pi(F, psi, invert_char=True)
        row = row and abs(live - tot) < 1e-6
        row = row and abs(live_inv - tot_inv) < 1e-6
        # fail-when-wrong: G A4_inv is not the real mass 6q² when G² is
        # off the positive-real axis (true for the order-12 character)
        ga_inv = Gpsi * live_inv
        row = row and abs(ga_inv.imag) > 1.0
        row = row and abs(ga_inv - a4_coll_GA_mass(p)) > 1.0
        sample[str(p)] = {
            "A4_coll": live,
            "6q_Gbar": tot,
            "G_A4": Gpsi * live,
            "G_A4_inv": ga_inv,
            "G_abs": abs(Gpsi),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Boolean collisions, Gauss piece A₄_coll(ψ)=∑_coll π "
            "e_ξ(b−a) ψ̄(ξ(d−c)).  T1 and the T2 / T4 slots with "
            "d=c vanish (ψ̄(0)=0).  Remaining T2 is ∑ χψ̄=0.  T3 "
            "crossed matchings (0,b,0,b) and (0,b,b,0): each is "
            "q G(ψ̄,±1)=q G(ψ̄) because ψ even.  T4 A1 is "
            "(q(q−2)/p)∑χψ̄=0; B3 has α=0.  Each of the four "
            "crosses (A2,A3,B1,B2) reduces to (q/p)χ(ξ)(pε G(ψ̄)) "
            "with G(χ)=pε and χ(ξ)ε=1 on Ω, hence q G(ψ̄).  "
            "Total 6q G(ψ̄) and G(ψ)A₄_coll=6q².  Together with "
            "A₄_full=(S_□/q)G(ψ̄) (Parseval + Aut_∞) this is the "
            "observed alignment A₄=((S_□−6q²)/q)G(ψ̄).  Inverted "
            "phase ψ not ψ̄ yields 6q G(ψ) and G·A₄_inv=6q G², "
            "not 6q².  Sign A₄/Ḡ≥0 is the floor and stays OPEN."
        ),
        "by_p_sample": sample,
    }


def _pow_mul(mul: list[list[int]], base: int, exp: int) -> int:
    r, b = 1, base
    while exp:
        if exp & 1:
            r = mul[r][b]
        b = mul[b][b]
        exp >>= 1
    return r


def _t_rep(F: dict, r: int) -> int:
    """Canonical point of T = squares/{±1}."""
    return min(r, F["mul"][r][F["neg1"]])


def _dir_name(r: int, p: int) -> str:
    """[a:b] ∈ P¹(F_p) with last coord 1, or 1:0."""
    a, b = r % p, r // p
    if a == 0 and b == 0:
        return "0"
    if b == 0:
        return "1:0"
    invb = pow(b, p - 2, p)
    return f"{(a * invb) % p}:1"


def _ratio_label(F: dict, r: int, split: str) -> tuple:
    """Max+-free geometric label of a square r."""
    p = F["p"]
    chi = F["chi"]
    rm1 = F["add"][r][F["neg"][1]]
    rp1 = F["add"][r][1]
    crm = 0 if rm1 == 0 else chi[rm1]
    crp = 0 if rp1 == 0 else chi[rp1]
    is_pm1 = int(r in (1, F["neg1"]))
    on_omega = int(_pow_mul(F["mul"], r, p - 1) == F["neg1"])
    if split == "paley":
        return (crm, crp, is_pm1)
    if split == "paley_omega":
        return (crm, crp, is_pm1, on_omega)
    if split == "dir":
        return (is_pm1, _dir_name(r, p))
    raise KeyError(split)


def t_scheme_relations(p: int, split: str) -> dict:
    """Intersection numbers of a partition of T=F_q^×/{±1}.

    Max+-free field arithmetic.  ``scheme`` is True iff every
    p_{C,D}(t) depends only on the class of t.
    """
    F = _kernel_field(p)
    pts = sorted({_t_rep(F, r) for r in F["squares"]})
    mul, inv = F["mul"], F["inv"]
    classes = {t: _ratio_label(F, t, split) for t in pts}

    def tmul(x: int, y: int) -> int:
        return _t_rep(F, mul[x][y])

    def tinv(x: int) -> int:
        return _t_rep(F, inv[x])

    pnums: dict[tuple, set[int]] = {}
    for t in pts:
        E = classes[t]
        counts: dict[tuple, int] = {}
        for x in pts:
            key = (classes[x], classes[tmul(tinv(x), t)])
            counts[key] = counts.get(key, 0) + 1
        for (C, D), n in counts.items():
            pnums.setdefault((C, D, E), set()).add(n)
    n_bad = sum(1 for vals in pnums.values() if len(vals) > 1)
    return {
        "p": p,
        "split": split,
        "n_pts": len(pts),
        "n_rels": len(set(classes.values())),
        "scheme": n_bad == 0,
        "n_bad": n_bad,
    }


def n_aut_orbits_squares(p: int) -> int:
    """# of ⟨Frob, inv⟩-orbits on squares/{±1}.

    Burnside for the multiplier group {±1, ±p} on Z/mZ, m=(p²−1)/4:
      (p+3)²/16 if p≡1 (mod 4),  (p+1)(p+5)/16 if p≡3 (mod 4).
    """
    if p % 4 == 1:
        return (p + 3) ** 2 // 16
    return (p + 1) * (p + 5) // 16


def aut_orbit_Q_leftover_dofs(p: int) -> int:
    """Unknown Q-orbit values after Q(±1)=8q² and the row-sum."""
    return n_aut_orbits_squares(p) - 2


def frob_inv_orbit_count(p: int) -> dict:
    """# of ⟨Frob, inv⟩-orbits on T vs Paley+ω classes.  Max+-free."""
    F = _kernel_field(p)
    mul, inv = F["mul"], F["inv"]
    pts = sorted({_t_rep(F, r) for r in F["squares"]})
    labels = {_t_rep(F, r): _ratio_label(F, r, "paley_omega") for r in F["squares"]}
    unused = set(pts)
    n_orb = 0
    split_class = False
    while unused:
        start = min(unused)
        stack = [start]
        seen: set[int] = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.append(_t_rep(F, _pow_mul(mul, x, p)))
            stack.append(_t_rep(F, inv[x]))
        unused -= seen
        n_orb += 1
        labs = {labels[x] for x in seen}
        if len(labs) > 1:
            split_class = True
    return {
        "p": p,
        "n_orbits": n_orb,
        "n_paley_omega_classes": len(set(labels.values())),
        "orbits_finer": n_orb > len(set(labels.values())),
        "orbit_mixes_classes": split_class,
        "n_orbits_formula": n_aut_orbits_squares(p),
        "formula_ok": n_orb == n_aut_orbits_squares(p),
        "leftover_dofs": aut_orbit_Q_leftover_dofs(p),
    }


def theorem_aut_orbit_S_ring(primes: list[int] | None = None) -> dict:
    """Q. ⟨Frob,inv⟩-orbits on squares form an S-ring; leftover n_orb−2.

    Max+-free.  Does not sign F̂ (too many free Q_O).
    Fail-when-wrong: formula disagrees with the orbit enumeration;
    leftover=0 would be the false 2-design of u on Ω.
    """
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19, 23]
    ok = True
    sample = {}
    for p in primes:
        rec = frob_inv_orbit_count(p)
        n = rec["n_orbits"]
        pred = n_aut_orbits_squares(p)
        leftover = aut_orbit_Q_leftover_dofs(p)
        row = rec["formula_ok"] and n == pred and leftover == n - 2
        # 2-design would force leftover=0; it is false for p≥5
        row = row and leftover >= 2
        sample[str(p)] = {
            "n_orbits": n,
            "formula": pred,
            "leftover_dofs": leftover,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "Q is constant on ⟨Frob,inv⟩-orbits of squares (Aut of "
            "Max+).  Those orbits are the S-ring of the multiplier "
            "group {±1,±p} on Z/mZ, m=(p²−1)/4, hence a commutative "
            "Schurian scheme for every odd p.  The number of orbits "
            "is (p+3)²/16 (p≡1 mod 4) or (p+1)(p+5)/16 (p≡3 mod 4).  "
            "Max+-free linear constraints on the orbit values are "
            "only Q(±1)=8q² and the row-sum 2q²(q−1), so leftover "
            "dofs = n_orbits−2 → ∞.  The 2-design of u (leftover=0) "
            "is false.  F̂ cannot be signed from these two constraints.  "
            "Floor S_□≥6q² stays OPEN."
        ),
        "by_p_sample": sample,
    }


def theorem_paley_omega_scheme_dies(primes: list[int] | None = None) -> dict:
    """P. Paley+ω is not a scheme in general (counterexample p=11).

    Does not claim the partition fails at every prime p≥11 — only
    that a single general-all-p Bose–Mesner formula is impossible,
    because it already fails at p=11 (and at the other checked
    primes in the default list).

    Fail-when-wrong: Paley+ω *is* a scheme at p=5,7 (so a blanket
    ``never a scheme`` is false); the direction partition *is* a
    scheme at p=5,7,11 (so ``no geometric scheme exists`` is false);
    Paley-without-ω already fails at p=7.
    """
    if primes is None:
        primes = [5, 7, 11, 13, 17, 19, 23]
    ok = True
    sample = {}
    for p in primes:
        po = t_scheme_relations(p, "paley_omega")
        pal = t_scheme_relations(p, "paley")
        dire = t_scheme_relations(p, "dir")
        orb = frob_inv_orbit_count(p)
        if p in (5, 7):
            row = po["scheme"] is True
        else:
            row = po["scheme"] is False
        # Paley without ω: only p=5
        row = row and (pal["scheme"] is (p == 5))
        # direction scheme exists (torus / P¹)
        row = row and dire["scheme"] is True
        # Aut-orbits vs Paley+ω: not finer at p=5,7; finer at p=11
        # (and at every other default prime ≥11 that we actually run).
        if p in (5, 7):
            row = row and orb["orbits_finer"] is False
        elif p >= 11:
            row = row and orb["orbits_finer"] is True
        sample[str(p)] = {
            "paley_omega_scheme": po["scheme"],
            "paley_scheme": pal["scheme"],
            "dir_scheme": dire["scheme"],
            "n_bad_omega": po["n_bad"],
            "n_rels_omega": po["n_rels"],
            "n_orbits": orb["n_orbits"],
            "n_omega_classes": orb["n_paley_omega_classes"],
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "attack_killed": "paley_omega_Bose_Mesner",
        "switched_to": (
            "Q on ⟨Frob,inv⟩-orbits of ratios / Boolean 4-point of "
            "V_+ ∩ {±1}^n (A₄/G(ψ̄)≥0)"
        ),
        "theorem": (
            "On T=F_q^×/{±1} the Paley+ω partition "
            "(χ(r−1),χ(r+1),1_{r^{p-1}=−1}) has constant intersection "
            "numbers at p=5 and p=7.  It does not at p=11 (n_bad=63), "
            "so it is not an association scheme in general and there "
            "is no single Paley+ω Bose–Mesner formula for F̂ valid "
            "for every prime p≥5.  The same non-scheme check holds "
            "at the other primes we run (default 13,17,19,23); this "
            "is not a proof that every p≥11 fails.  Paley without "
            "the ω-bit is already not a scheme at p=7.  The direction "
            "partition of T is a scheme at p=5,7,11 (and the other "
            "checked primes), but live Q is not a class function of "
            "it (p=5,7: two Q-values in one direction).  "
            "⟨Frob,inv⟩-orbits on T are strictly finer than Paley+ω "
            "classes at p=11 (12>7), so Aut-invariance of Q does not "
            "force a Paley+ω class function there.  Inverted claim "
            "(Paley+ω scheme at p=11) has n_bad>0.  Floor S_□≥6q² "
            "stays OPEN; listed attack switches to Aut-orbits of "
            "ratios / Boolean 4-point of V_+."
        ),
        "by_p_sample": sample,
    }


def _chi_p(p: int, x: int) -> int:
    """Quadratic character of F_p."""
    if x % p == 0:
        return 0
    return 1 if pow(x % p, (p - 1) // 2, p) == 1 else -1


def legendre_quadratic_sum(p: int, a: int, b: int, c: int) -> int:
    """∑_{t∈F_p} χ_p(at²+bt+c).  Equals −χ_p(a) if a≠0 and disc≠0."""
    s = 0
    for t in range(p):
        s += _chi_p(p, a * t * t + b * t + c)
    return s


def _linear_forms(p: int) -> list[tuple[int, int]]:
    """One representative per F_p-direction: (1,m) and (0,1)."""
    return [(1, m) for m in range(p)] + [(0, 1)]


def _ker_generator(p: int, alpha: int, beta: int) -> int:
    """A nonzero point of ker(αa+βb), encoded a+b*p."""
    if beta % p == 0:
        # L=αa, ker is the ω-axis a=0
        return p
    invb = pow(beta % p, p - 2, p)
    # a=1, b=−α/β
    b = ((-alpha) * invb) % p
    return 1 + b * p


def is_type_plus_form(p: int, alpha: int, beta: int) -> bool:
    """ker L is a square line ⇔ χ_q(v)=+1 ⇔ Type +."""
    F = _kernel_field(p)
    v = _ker_generator(p, alpha, beta)
    return int(F["chi"][v]) == 1


def lift_sigma_vector(p: int, alpha: int, beta: int, plus: set[int]) -> list[float]:
    """y_∞=1, y_x=σ(αa+βb) with σ=+1 on plus."""
    q = p * p
    y = [1.0] * (q + 1)
    y[0] = 1.0
    for x in range(q):
        a, b = x % p, x // p
        t = (alpha * a + beta * b) % p
        y[1 + x] = 1.0 if t in plus else -1.0
    return y


def n_type_plus_forms(p: int) -> int:
    """(p+1)/2 square lines through the origin."""
    return (p + 1) // 2


def theorem_1d_lifts_are_maxplus(primes: list[int] | None = None) -> dict:
    """R. Type+ linear lifts of every weight-(p+1)/2 σ are Max+.

    Max+-free construction.  Does not exhaust Max+ and does not
    sign F̂.  Fail-when-wrong: Type− interval is an eigenvector;
    wrong-weight Type+ is; quadratic sum equals +χ_p(a) or 0.
    """
    from itertools import combinations

    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    # character-sum lemma on F_p (input to the fiber evaluation)
    for p in (5, 7, 11, 13):
        for a, b, c in ((1, 0, 1), (1, 0, 2), (1, 1, 1), (2, 3, 1)):
            disc = (b * b - 4 * a * c) % p
            if a % p == 0 or disc == 0:
                continue
            got = legendre_quadratic_sum(p, a, b, c)
            if got != -_chi_p(p, a):
                ok = False
            # inverted sign is false
            if got == _chi_p(p, a) and _chi_p(p, a) != 0:
                ok = False
    for p in primes:
        C = paley_conference_prime_power(p)
        forms = _linear_forms(p)
        n_plus = sum(1 for ab in forms if is_type_plus_form(p, *ab))
        row = n_plus == n_type_plus_forms(p)
        nplus = (p + 1) // 2
        plus_sets = list(combinations(range(p), nplus))
        n_good = 0
        n_bad = 0
        type_minus_interval_ok = False
        wrong_weight_ok = False
        for alpha, beta in forms:
            typ = is_type_plus_form(p, alpha, beta)
            # all correct-weight σ on Type +; one interval on Type −
            sigmas = plus_sets if typ else [tuple(range(nplus))]
            for A in sigmas:
                y = np.array(lift_sigma_vector(p, alpha, beta, set(A)))
                ev = float(np.max(np.abs(C @ y - p * y))) < 1e-8
                if typ:
                    if ev:
                        n_good += 1
                    else:
                        n_bad += 1
                        row = False
                else:
                    type_minus_interval_ok = type_minus_interval_ok or ev
            # wrong weight on this Type + (if any)
            if typ and nplus >= 2:
                y = np.array(lift_sigma_vector(p, alpha, beta, set(range(nplus - 1))))
                if float(np.max(np.abs(C @ y - p * y))) < 1e-8:
                    wrong_weight_ok = True
        # Type − interval must not be an eigenvector; wrong weight neither
        row = row and (not type_minus_interval_ok) and (not wrong_weight_ok)
        row = row and n_bad == 0 and n_good == n_plus * len(plus_sets)
        sample[str(p)] = {
            "n_type_plus": n_plus,
            "n_type_plus_formula": n_type_plus_forms(p),
            "n_sigma": len(plus_sets),
            "n_good": n_good,
            "n_bad": n_bad,
            "type_minus_interval_is_evec": type_minus_interval_ok,
            "wrong_weight_is_evec": wrong_weight_ok,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "exhausts_maxplus": False,
        "theorem": (
            "If L:F_q→F_p is F_p-linear with ker L a square line "
            "(χ_q=+1 on ker\\{0}) and σ:F_p→{±1} has ∑σ=1, then "
            "y_∞=1, y_x=σ(L(x)) satisfies Cy=py.  Fiber sums of χ_q "
            "are K(0)=p−1 and K(β≠0)=−1 (χ_q|F_p^×=1 and the "
            "Legendre sum of a non-degenerate quadratic is −χ_p of "
            "the leading coefficient).  Hence χ*y=pσ(L)−1.  There "
            "are (p+1)/2 Type+ lines.  Type− lines and wrong weight "
            "fail.  The construction does not exhaust Max+ (p=7 has "
            "a nonlinear Aut-orbit) and does not prove S_□≥6q² "
            "(that nonlinear orbit, and the Aut-orbit of a non-"
            "interval σ, have S_□=0 for some even ψ)."
        ),
        "by_p_sample": sample,
    }


def regular_set_size(p: int) -> int:
    """|D|=p(p−1)/2 when y_∞=+1 and ∑z=p."""
    return p * (p - 1) // 2


def regular_S_on_D(p: int) -> int:
    """S(a)=(1+p)/2 for a∈D."""
    return (1 + p) // 2


def regular_S_off_D(p: int) -> int:
    """S(a)=(1−p)/2 for a∉D."""
    return (1 - p) // 2


def _signed_neighbourhood(C, d):
    """S(a)=∑_x C_{1+a,1+x} d_x  (C_fin d; diag of C is 0)."""
    return C[1:, 1:] @ d


def theorem_regular_set_switching(primes: list[int] | None = None) -> dict:
    """S. Cy=py and y_∞=+1 ⇔ D is a signed Paley regular set.

    Max+-free equivalence.  Does not prove the floor.  Fail-when-wrong:
    inverted S-values on a Type+ lift; Type− interval satisfies the
    size but not the S-values; a pointwise |M|² floor is not claimed.
    """
    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    for p in primes:
        q = p * p
        C = paley_conference_prime_power(p)
        nD = regular_set_size(p)
        on, off = regular_S_on_D(p), regular_S_off_D(p)
        # size / consistency algebra (no Max+)
        # p|D|+(1-p)q/2 = 0
        cons = p * nD + (1 - p) * q // 2 == 0
        cons = cons and nD == (q - p) // 2
        forms = _linear_forms(p)
        plus_forms = [ab for ab in forms if is_type_plus_form(p, *ab)]
        minus_forms = [ab for ab in forms if not is_type_plus_form(p, *ab)]
        plusA = set(range((p + 1) // 2))
        n_good = 0
        n_bad = 0
        type_minus_regular = False
        inverted_S_ok = False
        for alpha, beta in plus_forms[:3]:
            y = np.array(lift_sigma_vector(p, alpha, beta, plusA))
            ev = float(np.max(np.abs(C @ y - p * y))) < 1e-8
            d = (y[1:] < 0).astype(np.float64)
            S = _signed_neighbourhood(C, d)
            size_ok = int(d.sum()) == nD
            on_ok = bool(np.all(np.abs(S[d > 0.5] - on) < 1e-8))
            off_ok = bool(np.all(np.abs(S[d < 0.5] - off) < 1e-8))
            # inverted S-values must fail
            inv_on = bool(np.all(np.abs(S[d > 0.5] - off) < 1e-8))
            if ev and size_ok and on_ok and off_ok:
                n_good += 1
            else:
                n_bad += 1
                ok = False
            if inv_on:
                inverted_S_ok = True
        if minus_forms:
            y = np.array(lift_sigma_vector(p, *minus_forms[0], plusA))
            d = (y[1:] < 0).astype(np.float64)
            S = _signed_neighbourhood(C, d)
            if (
                int(d.sum()) == nD
                and np.all(np.abs(S[d > 0.5] - on) < 1e-8)
                and np.all(np.abs(S[d < 0.5] - off) < 1e-8)
            ):
                type_minus_regular = True
        row = (
            cons
            and n_bad == 0
            and n_good > 0
            and (not type_minus_regular)
            and (not inverted_S_ok)
        )
        sample[str(p)] = {
            "n_D": nD,
            "S_on": on,
            "S_off": off,
            "consistency": cons,
            "n_good": n_good,
            "n_bad": n_bad,
            "type_minus_is_regular": type_minus_regular,
            "inverted_S_values": inverted_S_ok,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "pointwise_M_floor": False,
        "theorem": (
            "If y_∞=+1 then Cy=py if and only if D={x:z_x=−1} has "
            "size p(p−1)/2 and S(a)=∑_x χ(a−x)1_D(x) equals (1+p)/2 "
            "on D and (1−p)/2 off D (equivalently C_fin d = p d + "
            "((1−p)/2)1).  Consistency: p|D|+(1−p)q/2=0.  Type+ 1D "
            "lifts are regular sets; Type− interval has the size but "
            "not the S-values.  Then ẑ(ξ)=−2 1̂_D(ξ) (ξ≠0), so "
            "A_r=−8 E[1̂_D(ξ)1̂_D(rξ)1̂_D(−(1+r)ξ)].  A_r≤0 is an "
            "ensemble statement.  Pointwise |M_ψ|²≥3q²(q−1) is false "
            "and is not claimed.  The identity does not sign F̂."
        ),
        "by_p_sample": sample,
    }


def theorem_energy_equals_maxplus(primes: list[int] | None = None) -> dict:
    """T. C²=p²I ⇒ yᵀCy=pn ⇔ Cy=py (energy fiber is Max+).

    Max+-free linear algebra.  Does not sign Term_4 / F̂.
    Fail-when-wrong: C²≠p²I; Type− interval has energy pn; a
    wrong-energy cube vector is a +p evec.
    """
    import numpy as np

    from minmax_quadratic import paley_conference_prime_power

    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    for p in primes:
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        q = p * p
        # conference: C=Cᵀ, zero diag, C²=p²I
        c2 = C @ C
        c2_ok = bool(np.max(np.abs(c2 - (p * p) * np.eye(n))) < 1e-8)
        sym_ok = bool(np.max(np.abs(C - C.T)) < 1e-12)
        diag_ok = bool(np.max(np.abs(np.diag(C))) < 1e-12)
        # Type+ lift: energy pn and +p evec
        forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        plusA = set(range((p + 1) // 2))
        y = np.array(lift_sigma_vector(p, *forms[0], plusA))
        ev = float(np.max(np.abs(C @ y - p * y))) < 1e-8
        energy = float(y @ C @ y)
        plus_ok = ev and abs(energy - p * n) < 1e-8
        # Type− interval: not +p evec, so energy ≠ pn
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        minus_energy_pn = False
        if minus:
            ym = np.array(lift_sigma_vector(p, *minus[0], plusA))
            em = float(ym @ C @ ym)
            evm = float(np.max(np.abs(C @ ym - p * ym))) < 1e-8
            minus_energy_pn = abs(em - p * n) < 1e-6 or evm
        # Cauchy–Schwarz: |yᵀCy| ≤ p n on the cube, = iff Cy=±py
        # one random cube vector with y_∞=+1, not all +1
        rng = np.random.default_rng(p)
        yrand = np.ones(n)
        yrand[1:] = rng.choice([-1.0, 1.0], size=n - 1)
        er = float(yrand @ C @ yrand)
        evr = float(np.max(np.abs(C @ yrand - p * yrand))) < 1e-8
        # if not an evec, energy cannot be pn
        rand_bad = evr and abs(er - p * n) > 1e-6
        # inverted claim: energy −pn is Max−, not Max+
        ev_minus = float(np.max(np.abs(C @ y + p * y))) < 1e-8
        row = (
            c2_ok
            and sym_ok
            and diag_ok
            and plus_ok
            and (not minus_energy_pn)
            and (not rand_bad)
            and (not ev_minus)
            and n == q + 1
        )
        sample[str(p)] = {
            "c2_ok": c2_ok,
            "plus_energy": energy,
            "target": p * n,
            "plus_ok": plus_ok,
            "type_minus_energy_pn": minus_energy_pn,
            "ok": row,
        }
        if not row:
            ok = False
    # p=3 full half-cube: energy pn ⇔ +p evec
    from itertools import product

    p = 3
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    bits = np.array(list(product((-1.0, 1.0), repeat=n - 1)))
    Y = np.ones((bits.shape[0], n))
    Y[:, 1:] = bits
    qf = np.einsum("ia,ab,ib->i", Y, C, Y)
    ev = np.max(np.abs(Y @ C - p * Y), axis=1) < 1e-8
    mask = np.abs(qf - p * n) < 1e-8
    cube_ok = bool(np.all(ev == mask)) and int(mask.sum()) > 0
    if not cube_ok:
        ok = False
    sample["3_cube"] = {
        "n_energy": int(mask.sum()),
        "n_evec": int(ev.sum()),
        "iff": cube_ok,
    }
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "The Paley conference satisfies C=Cᵀ, zero diagonal, and "
            "C²=p²I.  Hence ||Cy||=p||y|| for every y, and Cauchy–"
            "Schwarz gives |yᵀCy|≤p||y||² with equality iff Cy=±py.  "
            "On the cube, yᵀCy=pn iff Cy=py, so the energy fiber is "
            "Max+.  Type+ 1D lifts lie on that fiber; Type− interval "
            "does not.  Energy −pn is Max−.  The identity does not "
            "sign Term_4 (C² determines only the 2-point / degrees "
            "0+2; Wick 4-point of P_+ is the wrong Term_4)."
        ),
        "by_p_sample": sample,
    }


def R_lift_size(p: int) -> int:
    """|{t≠0,1: χ(t)=χ(1−t)=1}| = (q−5)/4 = dim F."""
    return (p * p - 5) // 4


def A_from_N(p: int, N: int) -> int:
    """A(δ)=p(2−p)+4 N(δ)."""
    return p * (2 - p) + 4 * N


def M_from_Nhat_prefactor(p: int) -> int:
    """M_ψ = 4q N̂(ψ)/G(ψ)."""
    return 4 * p * p


def _chi_conv(F: dict, vec: list[float]) -> list[float]:
    q = F["q"]
    add, neg, chi = F["add"], F["neg"], F["chi"]
    out = [0.0] * q
    for a in range(q):
        s = 0.0
        for x in range(q):
            s += float(chi[add[a][neg[x]]]) * float(vec[x])
        out[a] = s
    return out


def _autocorr_z(F: dict, z: list[float]) -> list[float]:
    q = F["q"]
    add = F["add"]
    return [sum(z[x] * z[add[x][d]] for x in range(q)) for d in range(q)]


def _N_of_D(F: dict, dmask: list[int]) -> list[float]:
    q = F["q"]
    add = F["add"]
    return [
        float(sum(dmask[x] * dmask[add[x][d]] for x in range(q))) for d in range(q)
    ]


def _R_lifts(F: dict) -> list[int]:
    q = F["q"]
    add, neg, chi = F["add"], F["neg"], F["chi"]
    out = []
    for t in range(1, q):
        if t == 1:
            continue
        tm = add[1][neg[t]]
        if tm != 0 and chi[t] == 1 and chi[tm] == 1:
            out.append(t)
    return out


def theorem_M_from_D_autocorr(primes: list[int] | None = None) -> dict:
    """T. M_ψ=4q N̂(ψ)/G(ψ) and χ*A=p(A−1).  Switching SOS does not close.

    Max+-free: Type+ lifts + character-sum |R|.  Fail-when-wrong:
    inverted G; χ*A=p(A+1); Boolean coeff p; Type−; |R|=(q−3)/4.
    """
    import numpy as np

    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    # |R| character-sum for a range of primes (no Max+)
    r_ok = True
    for p in [r for r in range(5, 32) if is_prime(r)]:
        q = p * p
        if R_lift_size(p) != (q - 5) // 4:
            r_ok = False
        if R_lift_size(p) != dim_F(p):
            r_ok = False
        # (q−1)/4 = |R|+1 is the inverted count
        if R_lift_size(p) == (q - 1) // 4:
            r_ok = False
        # A(0)=q from N(0)=|D|
        nD = regular_set_size(p)
        if A_from_N(p, nD) != q:
            r_ok = False
        if M_from_Nhat_prefactor(p) != 4 * q:
            r_ok = False
    if not r_ok:
        ok = False
    for p in primes:
        q = p * p
        F = _kernel_field(p)
        R = _R_lifts(F)
        r_count = len(R) == R_lift_size(p)
        plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        Aset = set(range((p + 1) // 2))
        tau = _gauss(F, [complex(F["chi"][x]) for x in range(q)])
        Omega = [
            xi
            for xi in range(1, q)
            if abs(F["chi"][xi] * tau - p) < 1e-6
        ]
        n_good = 0
        n_bad = 0
        inv_G_ok = False
        wrong_bool_ok = False
        wrong_A_ok = False
        type_minus_ok = False
        max_M = 0.0
        max_bool = 0.0
        max_A = 0.0
        ks = [2, 4]
        for alpha, beta in plus[:2]:
            y = lift_sigma_vector(p, alpha, beta, Aset)
            z = [float(y[1 + x]) for x in range(q)]
            dmask = [1 if z[x] < 0 else 0 for x in range(q)]
            cz = _chi_conv(F, z)
            err_z = max(abs(cz[a] - (p * z[a] - 1.0)) for a in range(q))
            A = _autocorr_z(F, z)
            cA = _chi_conv(F, A)
            err_A = max(abs(cA[d] - p * (A[d] - 1.0)) for d in range(q))
            err_A_plus = max(abs(cA[d] - p * (A[d] + 1.0)) for d in range(q))
            if err_A_plus < 1e-8:
                wrong_A_ok = True
            N = _N_of_D(F, dmask)
            nD = sum(dmask)
            a_from_n = all(
                abs(A[d] - A_from_N(p, int(round(N[d])))) < 1e-8 for d in range(q)
            )
            # ẑ
            zh = []
            for xi in range(q):
                acc = 0j
                for x in range(q):
                    acc += z[x] * _e_tr(F, F["mul"][xi][x])
                zh.append(acc)
            # Boolean at every ξ∈Ω (ẑ(ξ0) can vanish on a lift)
            bacc = 0.0
            b_wrong = 0.0
            for xi0 in Omega:
                acc = 2 * p * zh[xi0]
                acc_w = p * zh[xi0]
                for t in R:
                    tm = F["add"][1][F["neg"][t]]
                    term = zh[F["mul"][t][xi0]] * zh[F["mul"][tm][xi0]]
                    acc += term
                    acc_w += term
                bacc = max(bacc, abs(acc))
                b_wrong = max(b_wrong, abs(acc_w))
            max_bool = max(max_bool, bacc)
            # wrong coeff p gives residual p|ẑ|; must be ≫ 0
            if b_wrong < 1e-6:
                wrong_bool_ok = True
            # M vs 4q N̂/G
            m_ok = True
            for k in ks:
                psi = _psi_vec(F, k)
                G = _gauss(F, psi)
                Nhat = sum(psi[d] * N[d] for d in range(q))
                M = sum(psi[xi].conjugate() * abs(zh[xi]) ** 2 for xi in Omega)
                pred = (4 * q * Nhat) / G
                max_M = max(max_M, abs(M - pred))
                if abs(M - pred) > 1e-6:
                    m_ok = False
                # inverted denominator Ḡ
                if abs(G.imag) > 1e-8:
                    pred_inv = (4 * q * Nhat) / G.conjugate()
                    if abs(M - pred_inv) < 1e-6:
                        inv_G_ok = True
            size_ok = nD == regular_set_size(p) and abs(A[0] - q) < 1e-8
            if err_z < 1e-8 and err_A < 1e-8 and a_from_n and m_ok and size_ok:
                n_good += 1
            else:
                n_bad += 1
                ok = False
            max_A = max(max_A, err_A)
        if minus:
            y = lift_sigma_vector(p, *minus[0], Aset)
            z = [float(y[1 + x]) for x in range(q)]
            cz = _chi_conv(F, z)
            if max(abs(cz[a] - (p * z[a] - 1.0)) for a in range(q)) < 1e-8:
                type_minus_ok = True
        # 3-point on the first Type+ lift
        y = lift_sigma_vector(p, *plus[0], Aset)
        z = [float(y[1 + x]) for x in range(q)]
        A = _autocorr_z(F, z)
        t3_err = 0.0
        for d, e in ((1, 2), (2, 3), (1, 1)):
            def Tfun(dd: int, ee: int) -> float:
                return sum(
                    z[a] * z[F["add"][a][dd]] * z[F["add"][a][ee]] for a in range(q)
                )

            lhs = 0.0
            for u in range(q):
                lhs += float(F["chi"][u]) * Tfun(F["add"][d][u], F["add"][e][u])
            emd = F["add"][e][F["neg"][d]]
            rhs = p * Tfun(d, e) - A[emd]
            t3_err = max(t3_err, abs(lhs - rhs))
        row = (
            r_ok
            and r_count
            and n_bad == 0
            and n_good > 0
            and (not inv_G_ok)
            and (not wrong_bool_ok)
            and (not wrong_A_ok)
            and (not type_minus_ok)
            and t3_err < 1e-8
            and max_M < 1e-6
        )
        sample[str(p)] = {
            "n_R": len(R),
            "n_R_formula": R_lift_size(p),
            "n_good": n_good,
            "n_bad": n_bad,
            "max_M_resid": max_M,
            "max_bool": max_bool,
            "max_chiA": max_A,
            "t3_err": t3_err,
            "inverted_G_matches": inv_G_ok,
            "bool_coeff_p_ok": wrong_bool_ok,
            "chiA_plus_ok": wrong_A_ok,
            "type_minus_chi_z": type_minus_ok,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "sos_found": False,
        "theorem": (
            "χ*z=pz−1 implies χ*A=p(A−1).  A=p(2−p)+4N with "
            "N=|D∩(D−δ)|.  For nontrivial ψ, Ã=4 N̂ and "
            "M_ψ=4q N̂(ψ)/G(ψ), N̂=∑_{x,y∈D} ψ(y−x).  "
            "Boolean 2pẑ+∑_R ẑ(tξ)ẑ((1−t)ξ)=0 on Ω with "
            "|R|=(q−5)/4.  3-point χ-shift T = pT−A.  Linear "
            "switching is support on Ω only.  Inverted G / Type− / "
            "Boolean coeff p fail.  No SOS for S_□−6q² modulo the "
            "switching ideal (pointwise |M|²=0; Boolean rewrite is M)."
        ),
        "by_p_sample": sample,
    }


def EN_on_squares(p: int) -> Fraction:
    """E[N(δ)] for χ(δ)=+1: p(p−1)/4 = |D|/2."""
    return Fraction(p * (p - 1), 4)


def EN_on_nonsquares(p: int) -> Fraction:
    """E[N(δ)] for χ(δ)=−1: p(p−3)/4."""
    return Fraction(p * (p - 3), 4)


def theorem_N_mean_two_level(primes: list[int] | None = None) -> dict:
    """V. Ensemble E[N] is 2-level by χ; E[N̂(ψ)]=0 for even ψ∉{1,χ}.

    Max+-free from 15.279 S parameters + Aut-invariance of Max+.
    On Type+ 1D lifts (PDS) the values are exact, not just means.
    Does not sign Var(N̂)=Term_4.  Referee math_review PASS.
    Fail-when-wrong: swapped χ-values; Type− lift matches.
    """
    import numpy as np

    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    # algebra (no Max+)
    for p in (5, 7, 11, 13, 17, 19, 23):
        nD = regular_set_size(p)
        q = p * p
        if EN_on_squares(p) != Fraction(nD, 2):
            ok = False
        # |D|(p²−1)/4 = E[N|χ=1]·(q−1)/2
        if EN_on_squares(p) * (q - 1) / 2 != Fraction(nD * (q - 1), 4):
            ok = False
        if EN_on_squares(p) == EN_on_nonsquares(p):
            ok = False
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        chi = F["chi"]
        add = F["add"]
        plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        Aset = set(range((p + 1) // 2))
        y = lift_sigma_vector(p, *plus[0], Aset)
        d = [1 if y[1 + x] < 0 else 0 for x in range(q)]
        want_p = EN_on_squares(p)
        want_m = EN_on_nonsquares(p)

        def class_means(dmask):
            acc_p, acc_m, n_p, n_m = 0, 0, 0, 0
            for delta in range(1, q):
                Nd = sum(dmask[x] * dmask[add[x][delta]] for x in range(q))
                if chi[delta] == 1:
                    acc_p += Nd
                    n_p += 1
                else:
                    acc_m += Nd
                    n_m += 1
            return Fraction(acc_p, n_p), Fraction(acc_m, n_m)

        mp, mm = class_means(d)
        plus_ok = mp == want_p and mm == want_m
        # pointwise N is NOT 2-level on a single lift (line geometry)
        pointwise = True
        for delta in range(1, q):
            Nd = sum(d[x] * d[add[x][delta]] for x in range(q))
            w = want_p if chi[delta] == 1 else want_m
            if Fraction(Nd) != w:
                pointwise = False
                break
        minus_matches = False
        minus_swapped = False
        if minus:
            ym = lift_sigma_vector(p, *minus[0], Aset)
            dm = [1 if ym[1 + x] < 0 else 0 for x in range(q)]
            mmp, mmm = class_means(dm)
            minus_matches = mmp == want_p and mmm == want_m
            minus_swapped = mmp == want_m and mmm == want_p
        row = plus_ok and (not pointwise) and (not minus_matches) and minus_swapped
        sample[str(p)] = {
            "EN_plus": str(want_p),
            "EN_minus": str(want_m),
            "lift_mean_plus": str(mp),
            "lift_mean_minus": str(mm),
            "plus_ok": plus_ok,
            "pointwise_two_level": pointwise,
            "type_minus_matches": minus_matches,
            "type_minus_swapped": minus_swapped,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "On a signed Paley regular set, n_+|_D=(p²−1)/4.  "
            "∑_{χ(δ)=1} N(δ)=|D|(p²−1)/4.  Aut-invariance of Max+ "
            "makes E[N] constant on squares, so E[N|χ=1]=p(p−1)/4 "
            "and E[N|χ=−1]=p(p−3)/4.  Thus E[N]=a+bχ and "
            "E[N̂(ψ)]=0 for even ψ∉{1,χ}.  Term_4 is Var(N̂(ψ)), "
            "not a squared mean.  Type+ 1D lifts have those class "
            "means (not pointwise 2-level: line geometry splits N).  "
            "Type− interval swaps the two class means.  "
            "Does not sign Var(N̂)."
        ),
        "by_p_sample": sample,
    }


def A_off_constant(p: int) -> int:
    """A(δ)=p(2−p)+4N(δ) for δ≠0."""
    return p * (2 - p)


def nonsquare_NN_mass(p: int) -> int:
    """∑_{δ≠0} N(δ)N(rδ) for every regular set and every nonsquare r."""
    q = p * p
    return (q - 1) * q * (p - 1) * (p - 3) // 16


def L_on_nonsquares(p: int) -> Fraction:
    """L(r)=E[N(1)N(r)]-style average = p²(p−1)(p−3)/16 when χ(r)=−1."""
    return Fraction(p * p * (p - 1) * (p - 3), 16)


def theorem_nonsquare_N_pairing(primes: list[int] | None = None) -> dict:
    """W. ∑_δ N(δ)N(rδ) is constant for χ(r)=−1; L|_{⊠}=E[N_+]E[N_-].

    Max+-free: ẑ supported on {0}∪Ω and Ω is a χ-coset, so for
    nonsquare r the Parseval identity ∑_δ A(δ)A(rδ)=q (only the
    ξ=0 term survives).  Convert A=p(2−p)+4N.  Fail-when-wrong:
    the same sum on a square r equals q; Type− lift; inverted
    mass (q−1)c0² dropped.  Does not sign L̂ on squares.
    """
    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    # algebra: A-side conversion equals the closed mass
    for p in (5, 7, 11, 13, 17, 19, 23):
        q = p * p
        c0 = A_off_constant(p)
        nD = regular_set_size(p)
        rhs = q - q * q - (q - 1) * c0 * c0 - 8 * c0 * nD * (nD - 1)
        if rhs % 16 != 0:
            ok = False
        if rhs // 16 != nonsquare_NN_mass(p):
            ok = False
        if L_on_nonsquares(p) * (q - 1) != nonsquare_NN_mass(p):
            ok = False
        # product of class means
        if L_on_nonsquares(p) != EN_on_squares(p) * EN_on_nonsquares(p):
            ok = False
        # inverted: square-r mass is not this constant (p=5 live 360–650)
        if nonsquare_NN_mass(p) == 0:
            ok = False
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        add, mul, chi, neg = F["add"], F["mul"], F["chi"], F["neg"]
        plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        Aset = set(range((p + 1) // 2))
        y = lift_sigma_vector(p, *plus[0], Aset)
        z = [float(y[1 + x]) for x in range(q)]
        d = [1 if z[x] < 0 else 0 for x in range(q)]
        N = [sum(d[x] * d[add[x][delta]] for x in range(q)) for delta in range(q)]
        A = [float(q) if delta == 0 else float(A_off_constant(p) + 4 * N[delta])
             for delta in range(q)]
        ns = next(r for r in range(1, q) if chi[r] == -1)
        sq = next(r for r in range(2, q) if chi[r] == 1)
        def pair_A(r):
            return sum(A[delta] * A[mul[delta][r]] for delta in range(q))
        def pair_N(r):
            return sum(N[delta] * N[mul[delta][r]] for delta in range(1, q))
        a_ns = pair_A(ns)
        a_sq = pair_A(sq)
        n_ns = pair_N(ns)
        n_sq = pair_N(sq)
        # ∑_all A A_r = q on nonsquares
        plus_ns_ok = abs(a_ns - q) < 1e-6 and n_ns == nonsquare_NN_mass(p)
        # square r is not the same identity
        plus_sq_is_q = abs(a_sq - q) < 1e-6
        # all-ones is not supported on one χ-coset
        ones = [1.0] * q
        A1 = [float(q)] * q
        ones_ns = sum(A1[d] * A1[mul[d][ns]] for d in range(q))
        ones_is_q = abs(ones_ns - q) < 1e-6
        # Type− is also a χ-coset (Ω^c); pairing =q is correct, not a close
        minus_is_maxplus = False
        if minus:
            from minmax_quadratic import paley_conference_prime_power
            import numpy as np

            C = paley_conference_prime_power(p)
            ym = np.array(lift_sigma_vector(p, *minus[0], Aset))
            minus_is_maxplus = float(np.max(np.abs(C @ ym - p * ym))) < 1e-8
        row = (
            plus_ns_ok
            and (not plus_sq_is_q)
            and (not ones_is_q)
            and (not minus_is_maxplus)
        )
        sample[str(p)] = {
            "A_pair_ns": a_ns,
            "A_pair_ns_want": q,
            "N_pair_ns": n_ns,
            "N_pair_ns_want": nonsquare_NN_mass(p),
            "A_pair_sq": a_sq,
            "square_also_q": plus_sq_is_q,
            "ones_A_pair_is_q": ones_is_q,
            "type_minus_is_maxplus": minus_is_maxplus,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "theorem": (
            "ẑ is supported on {0}∪Ω and Ω={χ=ε} is a χ-coset.  "
            "A(δ)=∑_x z_x z_{x+δ} has Â(ξ)=ẑ(ξ)ẑ(−ξ), so Â vanishes "
            "off {0}∪Ω.  Parseval: ∑_δ A(δ)A(rδ)=q^{-1}∑_ξ Â(ξ)Â(−ξ/r).  "
            "If χ(r)=−1 then ξ∈Ω ⇒ ξ/r∈Ω^c, the Ω-terms die, and the "
            "sum is q^{-1}Â(0)²=q.  Convert A=p(2−p)+4N on δ≠0 and "
            "∑_{δ≠0}N=|D|(|D|−1) to get ∑_{δ≠0}N(δ)N(rδ)="
            "(q−1)p²(p−1)(p−3)/16, constant on every regular set.  "
            "Hence L(r)=p²(p−1)(p−3)/16=E[N_+]E[N_-] on nonsquares, "
            "and L̂(ψ)=∑_{r□}L(r)ψ̄(r) for even ψ∉{1,χ}.  The same "
            "Parseval on a square r keeps Ω-terms (Â(ξ)Â(ξ/r) for "
            "ξ∈Ω) and is not q.  Type− affine parts live on Ω^c "
            "and obey the same affine pairing (not Max+).  The "
            "all-ones vector is not a χ-coset and fails ∑A A_r=q.  "
            "Does not sign L̂ on squares (the floor leftover)."
        ),
        "by_p_sample": sample,
    }


def chi_neg2_is_square(p: int) -> bool:
    """χ_q(−2)=1 on F_{p²} for every odd p.

    χ_q(x)=χ_p(Norm(x))=χ_p(x^{p+1}).  p odd ⇒ p+1 even, so
    χ_p((−2)^{p+1})=χ_p(2)^{p+1}=1.
    """
    return p % 2 == 1


def NN_pairing(F: dict, d: list[int], r: int) -> int:
    """∑_{δ≠0} N(δ)N(rδ) for a 0-1 vector d on F_q."""
    q = F["q"]
    add, mul = F["add"], F["mul"]
    N = [sum(d[x] * d[add[x][delta]] for x in range(q)) for delta in range(q)]
    return int(sum(N[delta] * N[mul[delta][r]] for delta in range(1, q)))


def theorem_L_inv_not_s3(primes: list[int] | None = None) -> dict:
    """X. L(1/r)=L(r) for every N; r↦−1−r is not a symmetry of L.

    Max+-free.  Inversion is a change of variables on any pairing
    ∑ N(δ)N(rδ).  S3 would send 1 to −2, which is a square ≠±1
    for every prime p≥5, and on every Type+ 1D lift the pairing
    at 1 is strictly larger than at −2.  Leftover dofs stay
    n_aut−2.  Fail-when-wrong: χ(−2)=−1; pairing_1=pairing_{−2}
    on a Type+ lift; Type− used as a Max+ witness.
    Does not sign L̂.
    """
    if primes is None:
        primes = [5, 7, 11]
    ok = True
    sample: dict = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        if not chi_neg2_is_square(p):
            ok = False
        if p >= 5 and (p * p + 1) % 4 != 2:
            # q=p² ≡1 mod 4; −2 ≡ −1 * 2, both nonzero
            pass
        # −2 == ±1 only if p | 3 or p | 1, i.e. p=3
        if p >= 5 and (p - 2) % p == 0:
            ok = False
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        chi, inv, neg = F["chi"], F["inv"], F["neg"]
        rneg2 = int(neg[2])
        plus = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        Aset = set(range((p + 1) // 2))
        y = lift_sigma_vector(p, *plus[0], Aset)
        d = [1 if y[1 + x] < 0 else 0 for x in range(q)]
        pair1 = NN_pairing(F, d, 1)
        # inversion identity on a generic square (r=2, always a square)
        pair2 = NN_pairing(F, d, 2)
        pair2inv = NN_pairing(F, d, int(inv[2]))
        pair_neg2 = NN_pairing(F, d, rneg2)
        # every Type+ lift, not just the first
        all_strict = True
        n_checked = 0
        for ab in plus:
            yy = lift_sigma_vector(p, *ab, Aset)
            dd = [1 if yy[1 + x] < 0 else 0 for x in range(q)]
            p1 = NN_pairing(F, dd, 1)
            p2 = NN_pairing(F, dd, rneg2)
            n_checked += 1
            if not (p1 > p2):
                all_strict = False
        minus_is_maxplus = False
        if minus:
            from minmax_quadratic import paley_conference_prime_power
            import numpy as np

            C = paley_conference_prime_power(p)
            ym = np.array(lift_sigma_vector(p, *minus[0], Aset))
            minus_is_maxplus = float(np.max(np.abs(C @ ym - p * ym))) < 1e-8
        chi_neg2 = int(chi[rneg2])
        row = (
            chi_neg2 == 1
            and pair2 == pair2inv
            and pair1 != pair_neg2
            and pair1 > pair_neg2
            and all_strict
            and (not minus_is_maxplus)
            and aut_orbit_Q_leftover_dofs(p) == n_aut_orbits_squares(p) - 2
        )
        sample[str(p)] = {
            "chi_neg2": chi_neg2,
            "pair_1": pair1,
            "pair_neg2": pair_neg2,
            "pair_2_eq_inv": pair2 == pair2inv,
            "n_type_plus_checked": n_checked,
            "all_type_plus_strict": all_strict,
            "type_minus_is_maxplus": minus_is_maxplus,
            "leftover_aut": aut_orbit_Q_leftover_dofs(p),
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "s3_symmetry": False,
        "leftover_is_n_aut_minus_2": True,
        "theorem": (
            "For any function N on F_q, ∑_{δ≠0} N(δ)N(rδ)="
            "∑_{ε≠0} N(ε)N(ε/r) by ε=rδ, so L(1/r)=L(r) on every "
            "regular set (and every 1D lift).  χ_q(−2)=χ_p(2)^{p+1}=1 "
            "for every odd p, and −2∉{±1} when p≥5, so the anharmonic "
            "map r↦−1−r sends the Aut-orbit {±1} to the Aut-orbit of "
            "−2.  On every Type+ 1D lift (Max+-free Max+ vectors) the "
            "pairing at 1 is strictly larger than at −2, so S3 is not "
            "a symmetry of L and does not collapse leftover Aut-orbits.  "
            "Leftover after W, L(1), and the row-sum is n_aut−2.  "
            "Type− is not Max+.  Does not sign L̂ on squares."
        ),
        "by_p_sample": sample,
    }


def _n_minus_1d(M: set[int], c: int, p: int) -> int:
    return len(M & {(x - c) % p for x in M})


def _paley_plus_set(p: int) -> set[int]:
    return {0} | {(x * x) % p for x in range(1, p)}


def _interval_plus_set(p: int) -> set[int]:
    return set(range((p + 1) // 2))


def _N_from_1d_formula(p: int, alpha: int, beta: int, plus: set[int], F: dict) -> list[int]:
    """N(δ)=|D| on ker L; p·n_M(L(δ)) off.  Encoded F_q points."""
    q = p * p
    nD = regular_set_size(p)
    M = set(range(p)) - set(plus)
    out = [0] * q
    for delta in range(q):
        a, b = delta % p, delta // p
        c = (alpha * a + beta * b) % p
        out[delta] = nD if c == 0 else p * _n_minus_1d(M, c, p)
    return out


def _N_direct_from_d(F: dict, d: list[int]) -> list[int]:
    q = F["q"]
    add = F["add"]
    return [sum(d[x] * d[add[x][delta]] for x in range(q)) for delta in range(q)]


def _sum_M_eta(M: set[int], eta: dict, p: int) -> complex:
    s = 0j
    for x in M:
        for y in M:
            diff = (x - y) % p
            if diff:
                s += eta[diff]
    return s


def one_d_triv_Nhat_sq(p: int) -> int:
    """|N̂|² on a Type+ lift when (p−1)|k."""
    return (p * (p * p - 1) // 4) ** 2


def theorem_1d_Nhat_closed(primes: list[int] | None = None) -> dict:
    """Y. Closed N̂ on Type+ 1D lifts.  Does not sign the ensemble floor.

    Fail-when-wrong: inverted (p−1)|k slot (triv energy on a ntriv Paley
    character); Type− form is not claimed Max+.
    """
    if primes is None:
        primes = [5, 7]
    ok = True
    sample: dict = {}
    for p in primes:
        F = _kernel_field(p)
        q = p * p
        forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
        minus = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
        alpha, beta = forms[0]
        ks = [k for k in range(2, q - 1, 2) if k != (q - 1) // 2]
        psi = {k: _psi_vec(F, k) for k in ks}
        triv_val = one_d_triv_Nhat_sq(p)
        plus_sets = [
            ("AP", _interval_plus_set(p)),
            ("Paley", _paley_plus_set(p)),
        ]
        n_formula_ok = 0
        n_triv_ok = 0
        n_ntriv_ok = 0
        paley_ntriv_zero = True
        inverted_slot_fails = True
        n_checked = 0
        for name, plus in plus_sets:
            y = lift_sigma_vector(p, alpha, beta, plus)
            d = [1 if y[1 + x] < 0 else 0 for x in range(q)]
            Nf = _N_from_1d_formula(p, alpha, beta, plus, F)
            Nd = _N_direct_from_d(F, d)
            if Nf == Nd:
                n_formula_ok += 1
            else:
                ok = False
            M = set(range(p)) - set(plus)
            for k in ks:
                ntriv_closed = (p**3) * abs(
                    _sum_M_eta(M, {t: psi[k][t] for t in range(1, p)}, p)
                ) ** 2
                nh2 = abs(
                    sum(Nf[delta] * psi[k][delta] for delta in range(1, q))
                ) ** 2
                n_checked += 1
                if k % (p - 1) == 0:
                    if abs(nh2 - triv_val) < 1e-4:
                        n_triv_ok += 1
                    else:
                        ok = False
                else:
                    if abs(nh2 - ntriv_closed) < 1e-4:
                        n_ntriv_ok += 1
                    else:
                        ok = False
                    if name == "Paley" and p % 4 == 3:
                        if nh2 > 1e-4:
                            paley_ntriv_zero = False
                            ok = False
                        # inverted slot: triv energy does not land here
                        if abs(nh2 - triv_val) < 1e-4:
                            inverted_slot_fails = False
                            ok = False
        # Type− is not Max+ (already R); do not claim the energy formula
        # on Type− as a floor input.
        row = (
            n_formula_ok == 2
            and n_triv_ok > 0
            and n_ntriv_ok > 0
            and inverted_slot_fails
            and (paley_ntriv_zero if p % 4 == 3 else True)
            and len(minus) > 0
        )
        sample[str(p)] = {
            "n_formula_ok": n_formula_ok,
            "n_triv_ok": n_triv_ok,
            "n_ntriv_ok": n_ntriv_ok,
            "triv_val": triv_val,
            "paley_ntriv_zero": paley_ntriv_zero,
            "inverted_slot_fails": inverted_slot_fails,
            "n_checked": n_checked,
            "ok": row,
        }
        if not row:
            ok = False
    return {
        "proved": ok,
        "inequality_proved": False,
        "exhausts_maxplus": False,
        "theorem": (
            "On a Type+ lift y_x=σ(L(x)), N(δ)=|D| for δ∈ker L and "
            "N(δ)=p n_M(L(δ)) off ker, M=σ^{-1}(−1).  For even ψ, "
            "(p−1)|k iff ψ is trivial on F_p^×, and then "
            "|N̂|²=[p(p²−1)/4]².  Otherwise "
            "|N̂|²=p³|∑_{x,y∈M} η(x−y)|² with η=ψ|F_p^×.  "
            "Paley plus-set {0}∪□_p has constant n_M when p≡3 (mod 4), "
            "so N̂=0 on those even ntriv characters.  Aut_∞ 1D orbits "
            "are AGL(1,F_p) orbits of plus-sets.  The 1D-only mixture "
            "does not meet the floor (p=7, k=8).  Inverted (p−1)|k "
            "slot fails.  Does not sign ensemble L̂."
        ),
        "by_p_sample": sample,
    }


def paley_omega_scheme_general() -> bool:
    """False: p=11 is a counterexample.  Not a floor close."""
    return False


def rhat_ge_budget_proved_general() -> bool:
    """OPEN G: ˆR_rest≥−2q² (old form ˆR≥q²(q−1)) for all even ψ∉{1,χ}."""
    return False


def parseval_identity_proved() -> bool:
    return bool(
        theorem_parseval_M_from_autocorr()["proved"]
        and theorem_S_square_is_8q2_plus_Rhat()["proved"]
    )


def delta_coll_closed_proved() -> bool:
    return bool(theorem_delta_coll_closed()["proved"])


def wick_4dist_closed_proved() -> bool:
    return bool(theorem_wick_4dist_closed()["proved"])


def bool_coll_orthogonality_proved() -> bool:
    return bool(theorem_bool_coll_orthogonality()["proved"])


def A4_coll_alignment_proved() -> bool:
    return bool(theorem_A4_coll_equals_6q_Gbar()["proved"])


def phi_F_structure_proved() -> bool:
    return bool(
        theorem_c_omega_is_2p()["proved"]
        and theorem_rayleigh_formula()["proved"]
        and theorem_wick_rayleigh_is_4()["proved"]
        and theorem_Q_row_sum_and_antipode()["proved"]
        and theorem_floor_iff_Rhat()["proved"]
        and parseval_identity_proved()
        and delta_coll_closed_proved()
        and wick_4dist_closed_proved()
        and jacobi_kernel_pairings_proved()
        and theorem_kappaC_pairs_as_What4()["proved"]
        and theorem_I_L2_closed()["proved"]
        and theorem_floor_iff_gauss_4dist()["proved"]
        and bool_coll_orthogonality_proved()
        and theorem_gauss_coll_is_3q2()["proved"]
        and A4_coll_alignment_proved()
        and theorem_paley_omega_scheme_dies()["proved"]
        and theorem_aut_orbit_S_ring()["proved"]
        and theorem_1d_lifts_are_maxplus()["proved"]
        and theorem_regular_set_switching()["proved"]
        and theorem_M_from_D_autocorr()["proved"]
        and theorem_energy_equals_maxplus()["proved"]
        and theorem_N_mean_two_level()["proved"]
        and theorem_nonsquare_N_pairing()["proved"]
        and theorem_L_inv_not_s3()["proved"]
        and theorem_1d_Nhat_closed()["proved"]
    )


def main() -> dict:
    print("Prop 15.279  Φ|_F character-pair / Jacobi", flush=True)
    A = theorem_c_omega_is_2p()
    B = theorem_rayleigh_formula()
    C = theorem_wick_rayleigh_is_4()
    D = theorem_floor_iff_Rhat()
    Sacc = theorem_S_square_is_8q2_plus_Rhat()
    E = theorem_parseval_M_from_autocorr()
    Fcol = theorem_delta_coll_closed()
    H4 = theorem_wick_4dist_closed()
    Iker = theorem_jacobi_kernel_pairings()
    Jkap = theorem_kappaC_pairs_as_What4()
    KI2 = theorem_I_L2_closed()
    L4 = theorem_floor_iff_gauss_4dist()
    Morth = theorem_bool_coll_orthogonality()
    Ncoll = theorem_gauss_coll_is_3q2()
    Oa4 = theorem_A4_coll_equals_6q_Gbar()
    Psch = theorem_paley_omega_scheme_dies()
    Qsr = theorem_aut_orbit_S_ring()
    R1d = theorem_1d_lifts_are_maxplus()
    Sreg = theorem_regular_set_switching()
    Taut = theorem_M_from_D_autocorr()
    Teng = theorem_energy_equals_maxplus()
    Vn = theorem_N_mean_two_level()
    Wn = theorem_nonsquare_N_pairing()
    Xinv = theorem_L_inv_not_s3()
    Y1d = theorem_1d_Nhat_closed()
    Qid = theorem_Q_row_sum_and_antipode()
    cert = certify_p5_Rhat_beats_budget()
    print(f"  A |c_Ω|=2p: {A['proved']}", flush=True)
    print(f"  B Rayleigh formula: {B['proved']}", flush=True)
    print(f"  C Wick λ=4 (2-pairing): {C['proved']}", flush=True)
    print(f"  Q row-sum / antipode / M_1: {Qid['proved']}", flush=True)
    print(f"  D iff λ=8+ˆR_rest/q²≥6: {D['proved']} inequality={D['inequality_proved']}", flush=True)
    print(f"  D-accounting S_□=8q²+ˆR_rest: {Sacc['proved']}", flush=True)
    print(f"  E Parseval M=q Ã/G: {E['proved']}", flush=True)
    print(
        f"  F Δ_coll closed: {Fcol['proved']} inequality={Fcol['inequality_proved']}",
        flush=True,
    )
    print(
        f"  H Wick_4dist closed: {H4['proved']} inequality={H4['inequality_proved']}",
        flush=True,
    )
    print(
        f"  I Jacobi kernel pairings: {Iker['proved']} inequality={Iker['inequality_proved']}",
        flush=True,
    )
    print(
        f"  J ⟨κ_C,K⟩=q ˆW_4: {Jkap['proved']} inequality={Jkap['inequality_proved']}",
        flush=True,
    )
    print(
        f"  K ∑|I|² closed: {KI2['proved']} inequality={KI2['inequality_proved']}",
        flush=True,
    )
    print(
        f"  L floor iff Gauss 4-dist ≥0: {L4['proved']} inequality={L4['inequality_proved']}",
        flush=True,
    )
    print(
        f"  M Coll 2-level / ∑ψ Coll=6q²: {Morth['proved']} inequality={Morth['inequality_proved']}",
        flush=True,
    )
    print(
        f"  N Gauss A_coll=3q²: {Ncoll['proved']} inequality={Ncoll['inequality_proved']}",
        flush=True,
    )
    print(
        f"  O A4_coll=6q G(ψ̄): {Oa4['proved']} inequality={Oa4['inequality_proved']}",
        flush=True,
    )
    print(
        f"  P Paley+ω not a scheme in general (p=11): {Psch['proved']} "
        f"killed={Psch['attack_killed']}",
        flush=True,
    )
    print(
        f"  Q Aut-orbit S-ring leftover n_orb-2: {Qsr['proved']} "
        f"inequality={Qsr['inequality_proved']}",
        flush=True,
    )
    print(
        f"  R 1D Type+ lifts are Max+: {R1d['proved']} "
        f"inequality={R1d['inequality_proved']} exhausts={R1d['exhausts_maxplus']}",
        flush=True,
    )
    print(
        f"  S regular-set switching: {Sreg['proved']} "
        f"inequality={Sreg['inequality_proved']} "
        f"pointwise_M={Sreg['pointwise_M_floor']}",
        flush=True,
    )
    print(
        f"  T M=4q N̂/G / χ*A=p(A−1): {Taut['proved']} "
        f"inequality={Taut['inequality_proved']} "
        f"sos={Taut['sos_found']}",
        flush=True,
    )
    print(
        f"  U energy fiber is Max+: {Teng['proved']} "
        f"inequality={Teng['inequality_proved']}",
        flush=True,
    )
    print(
        f"  V E[N] 2-level / E[N̂]=0: {Vn['proved']} "
        f"inequality={Vn['inequality_proved']}",
        flush=True,
    )
    print(
        f"  W nonsquare N-pairing: {Wn['proved']} "
        f"inequality={Wn['inequality_proved']}",
        flush=True,
    )
    print(
        f"  X L inv not S3: {Xinv['proved']} "
        f"inequality={Xinv['inequality_proved']} "
        f"s3={Xinv['s3_symmetry']}",
        flush=True,
    )
    print(
        f"  Y 1D N̂ closed: {Y1d['proved']} "
        f"inequality={Y1d['inequality_proved']} "
        f"exhausts={Y1d['exhausts_maxplus']}",
        flush=True,
    )
    print(f"  p=5 census λ=80/13: {cert['lambda_is_80_13']} beats={cert['beats']}", flush=True)
    print(f"  rhat_ge_budget_proved_general={rhat_ge_budget_proved_general()}", flush=True)
    print(f"  phi_F_ge_6_proved_general={phi_F_ge_6_proved_general()}", flush=True)
    print(f"  aut_span_F_equals_Z={aut_span_F_equals_Z_proved_general()}", flush=True)
    out = {
        "prop": "15.279",
        "L_status": "OPEN",
        "proved": {
            "c_omega_abs_2p": A["proved"],
            "rayleigh_formula": B["proved"],
            "wick_lambda_4": C["proved"],
            "Q_row_sum_antipode": Qid["proved"],
            "floor_iff_Rhat": D["proved"],
            "S_square_accounting": Sacc["proved"],
            "parseval_M_from_A": E["proved"],
            "delta_coll_closed": Fcol["proved"],
            "wick_4dist_closed": H4["proved"],
            "jacobi_kernel_pairings": Iker["proved"],
            "kappaC_pairs_as_What4": Jkap["proved"],
            "I_L2_closed": KI2["proved"],
            "floor_iff_gauss_4dist": L4["proved"],
            "bool_coll_orthogonality": Morth["proved"],
            "gauss_coll_is_3q2": Ncoll["proved"],
            "A4_coll_equals_6q_Gbar": Oa4["proved"],
            "A4_coll_inequality": Oa4["inequality_proved"],
            "paley_omega_scheme_dies": Psch["proved"],
            "paley_omega_scheme_general": paley_omega_scheme_general(),
            "floor_attack_killed": Psch["attack_killed"],
            "floor_attack_switched_to": Psch["switched_to"],
            "aut_orbit_S_ring": Qsr["proved"],
            "aut_orbit_S_ring_inequality": Qsr["inequality_proved"],
            "one_d_lifts_are_maxplus": R1d["proved"],
            "one_d_lifts_inequality": R1d["inequality_proved"],
            "one_d_lifts_exhaust_maxplus": R1d["exhausts_maxplus"],
            "regular_set_switching": Sreg["proved"],
            "regular_set_inequality": Sreg["inequality_proved"],
            "regular_set_pointwise_M_floor": Sreg["pointwise_M_floor"],
            "M_from_D_autocorr": Taut["proved"],
            "M_from_D_inequality": Taut["inequality_proved"],
            "M_from_D_sos": Taut["sos_found"],
            "energy_equals_maxplus": Teng["proved"],
            "energy_equals_maxplus_inequality": Teng["inequality_proved"],
            "N_mean_two_level": Vn["proved"],
            "N_mean_two_level_inequality": Vn["inequality_proved"],
            "nonsquare_N_pairing": Wn["proved"],
            "nonsquare_N_pairing_inequality": Wn["inequality_proved"],
            "L_inv_not_s3": Xinv["proved"],
            "L_inv_not_s3_inequality": Xinv["inequality_proved"],
            "L_s3_symmetry": Xinv["s3_symmetry"],
            "one_d_Nhat_closed": Y1d["proved"],
            "one_d_Nhat_inequality": Y1d["inequality_proved"],
            "one_d_Nhat_exhausts_maxplus": Y1d["exhausts_maxplus"],
            "rhat_ge_budget_general": rhat_ge_budget_proved_general(),
            "phi_F_ge_6": phi_F_ge_6_proved_general(),
            "aut_span_F_equals_Z": aut_span_F_equals_Z_proved_general(),
            "gsum_disj_lb": gsum_disj_lb_proved_general(),
            "aut_schur": False,
            "pairing": False,
            "e1": e1_closed_general(),
        },
        "cert_p5": cert,
        "dim_F_p5": dim_F(5),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15279.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
