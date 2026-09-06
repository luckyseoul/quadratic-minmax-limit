# The actual cross radius sharpens the original weighted-shell upper

2026-09-06. Analytic finite theorem, uniform actual-diagonal stability,
and one exact fixed-metric evaluation. No mathematical computation,
parameter scan, signing census, or numerical optimization was run.
This is a separate strengthening of the upper argument, not a change
to the source, covariance, spectral law, or original norm objective.

## 1. Finite actual-source upper with the full valid metric window

Use an actual complete paired signing and a feasible positive diagonal:

    K=[[A,B],[B^T,-A]],       N=2n,
    D-K>=0, D+K>=0,         S=tr D, dbar=S/N.

Set H_B=[[0,B],[B^T,0]], L=D^(-1/2)H_BD^(-1/2), and

    r=||L||op=||D_L^(-1/2) B D_R^(-1/2)||op in (0,1],
    V=L/r,                 ||V||op=1.

The upper r<=1 follows by block-sign conjugation and averaging of
D+-K. Positivity r>0 follows because an actual cross signing is nonzero.
For any actual original/weighted refined cell, let c=x^TBy be its
common ORIGINAL cross energy, theta its representative chosen within
that final cell, and M=M_theta its genuine positive field covariance
from the 381-line weighted-shell theorem listed in Section 6. Thus

    M>=0,       diag M=wn,       tr M=wN^2/2,
    w=k+v,      0<=k<=w<=1.

For every |t|<1 put b=1-|t| and define

    eta=t/r,         P=D-eta H_B,
    E=bD,            F=(I-tV)^(-1),
    T=tr(MP^(-1)),   R=tr(DP^(-1)MP^(-1)),
    v_o=2c/(rS).

Then P>0, 0<=E<=P, and |v_o|<=1. Indeed P=D^(1/2)(I-tV)D^(1/2),
while the cross-only contraction gives |z^TH_Bz|<=r z^TDz=rS.
The exact ORIGINAL shell radius is

    z^TPz=S-2eta c=S(1-tv_o),
    z^TPz-tr E=S(|t|-tv_o).

The Gaussian Boolean ellipsoid remainder, including Cauchy--Schwarz
on its diagonal sum, therefore gives the finite valid field-width upper

    E max_(z in cell) g^Tz <= B_D(theta,c;t)
      :=sqrt(S) {sqrt((|t|-tv_o)(T-bR))
                                  +sqrt(kappa)b sqrt(R)},       (1)
    kappa=2/pi,        g~N(0,M).

All square roots are nonnegative: T-bR is the trace of the positive
matrix product corresponding to P-E. Equivalently, in the old eta
coordinate, the valid range is |eta|r<1, and the first radius factor
is |eta|r-eta(2c/S). Both the enlarged range and the larger diagonal
remainder use the ACTUAL operator radius r, not a bulk spectral edge.

For completeness, the ellipsoid inequality used in (1) is

    E max g^Tz <= sqrt((q-tr E)
          [tr(MP^(-1))-tr(E P^(-1)MP^(-1))])
          +sqrt(kappa) sum_i E_ii sqrt((P^(-1)MP^(-1))_ii),

for z^TPz=q and diagonal 0<=E<=P. Complete the square with P, retain
the Boolean remainders E_ii(|(P^(-1)g)_i|/a-1)^2, take expectation,
and optimize a>0. Substitution E=bD and Cauchy--Schwarz give (1).
This does not require independent field coordinates or commutation.

The change also improves the bound at every old admissible eta.
For fixed P,M,q, its Cauchy--Schwarz upper with E=hD, divided by
sqrt(S), is sqrt((q/S-h)(T-hR))+sqrt(kappa)h sqrt(R).
Where the first product is positive, its derivative is at most
-(1-sqrt(kappa))sqrt(R), by arithmetic--geometric mean. Continuity
covers its boundary. Increasing h from 1-|eta| to 1-|eta|r hence
never worsens this upper. The latter is the largest admissible h
within the family E=hD, since the extreme eigenvalues of L are +-r.

Consequently (1), infimized over |t|<1, replaces the field-width
term in the weighted-shell theorem's (4.6), for EVERY actual cell.
Its original drift, independent cushion, padding, and cell-selection
errors are unchanged. No source optimality is used in this theorem.

## 2. Uniform diagonal stability, even for a varying actual radius

Define delta=S tr(D^(-1))/N^2-1, and let c_D=x^T D_L^(-1/2)
B D_R^(-1/2)y be the representative's weighted cross value. Put

    v_D=c_D/(nr),          bar T=tr(MF), bar R=tr(MF^2),
    B_flat(theta;t)=sqrt(N) {sqrt((|t|-tv_D)(bar T-b bar R))
                                      +sqrt(kappa)b sqrt(bar R)}.

This is a numerical reference functional using the SAME actual M,V.
It need not be an exact ellipsoid constraint on every state of a
weighted bin. For 0<=delta<=1,

    |B_D(theta,c;t)-B_flat(theta;t)|
       <=3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|t|).                (2)

Here the constant is independent of r>0, including if r varies with n.
The important normalized energy comparison is

                         |v_o-v_D|<=2sqrt(delta).              (3)

To prove it, write H_B/r=D^(1/2)V D^(1/2) and ||V||op=1. The elementary
dispersion bound sum_i(sqrt(d_i)-sqrt(dbar))^2<=S delta gives

    Phi(H_B/r-dbar V)<=S sqrt(delta).

This matrix is cross-only, so its Phi is its cross Boolean bilinear
norm. Dividing the resulting bound on |c/r-dbar c_D/r| by n dbar
proves (3). No division of an r-independent error by r is needed.

Here are the remaining trace estimates proving (2). Put J=sqrt(dbar)
D^(-1/2), and E_delta=(wN^2/2)sqrt(delta)(1+sqrt(1+delta)). The constant
diagonal of M gives

    ||J M J-M||_1<=E_delta,                                  (4)

by expanding (J-I)MJ+M(J-I) and applying the two Frobenius-product
inequalities. Specifically, the three squared Frobenius norms of
M^(1/2), J M^(1/2), and (J-I)M^(1/2) are respectively wN^2/2,
(wN^2/2)(1+delta), and at most (wN^2/2)delta.
Moreover

    dbar T=tr(JMJ F),       dbar R=tr(JMJ F^2),
    ||F||op<=1/b,          0<=F-bF^2<=I/(4b).

Hence

    |dbar R-bar R|<=E_delta/b^2,
    |dbar(T-bR)-(bar T-b bar R)|<=E_delta/(4b),
    0<=bar T-b bar R<=wN^2/(8b).                             (5)

Both radius coefficients are in [0,2], and their difference is at
most 2|t|sqrt(delta) by (3). The difference between the products
inside the first square roots, after rewriting B_D with outer
factor sqrt(N), is therefore at most

    E_delta/(2b)+|t|wN^2sqrt(delta)/(4b)
                                  <=wN^2sqrt(delta)/b.

Use |sqrt(a)-sqrt(a')|<=sqrt(|a-a'|) for a,a'>=0. The second-term
difference is at most sqrt(kappa E_delta) before that outer factor.
Their sum is bounded by the right side of (2). The cases w=0 or
delta=0 are immediate. This proves uniform all-cell stability.

Thus for any fixed b_0>0 the field term in the complete actual-cell
upper can be replaced by inf_(|t|<=1-b_0) B_flat(theta;t), adding
3sqrt(w)N^(3/2)delta^(1/4)/sqrt(b_0) to its existing errors. When
delta tends to zero this is o(N^(3/2)). The fixed compact t-window
is essential; no uncontrolled n-dependent endpoint limit is taken.

## 3. Exact pure-cross reference functional

For an actual original zero-source cell p=q_A=0, the 303-line weighted
compatibility theorem compares M_theta with the genuine positive
pure-cross field

    M_0=wnI-k c_D L=wn(I-sV),
    u=c_D/n,         a=k/w,         s=a u r,                  (6)

with o(n^(3/2)) field-width error when delta tends to zero and w>0.
For w=0 both fields vanish. Its finite non-augmented error is at most
2sqrt(k log(2))N^(3/2)delta^(1/4). Formula (2) applies to M_0 directly:
it is positive and has the required constant diagonal.

Let nu_r be the empirical law of the n squared singular values of
W/r, where W=D_L^(-1/2)B D_R^(-1/2). It includes zeros and has support
[0,1]. For 0<t<1 define

    A_s,t(y)=[1+(t^2-2t+s(1-2t))y+s t^2 y^2]/(1-t^2 y)^2,
    B_s,t(y)=[1+(t^2-2st)y]/(1-t^2 y)^2.

Evenness of the spectrum of V and the commutation in (6) give

    B_flat(M_0;t)/(2n^(3/2))
      =sqrt(w) {t sqrt((1-u/r) integral A_s,t dnu_r)
               +(1-t)sqrt(kappa) sqrt(integral B_s,t dnu_r)}.  (7)

These are exact reference traces. For instance, F-bF^2=t(I-V)F^2
for t>0; pair its two eigenvalue signs after multiplying by I-sV.
This yields A_s,t, while the paired trace of (I-sV)F^2 yields B_s,t.
The covariance coefficient is s=a u r, NOT a u/r. The original
cell radius and actual weighted cross value were compared by (2)-(3),
not silently equated at finite n.

## 4. One exact rational metric removes the literal middle-law gap

Take the centered sign parameters w=1, k=kappa. Consider the specified
reference data

    r^2=4/5,       u=4/5,       nu_r=(delta_0+delta_1)/2.

Equivalently, the unscaled squared-singular law is
(delta_0+delta_(4/5))/2 AND its actual operator radius is r=2/sqrt(5).
Then u/r=r and s=kappa r^3. At the single fixed choice t=9/10,

    integral A_s,t dnu_r=(461+100s)/722,
    integral B_s,t dnu_r=(18461-18000s)/722.                   (8)

Indeed A_s,t(0)=B_s,t(0)=1, A_s,t(1)=(1+s)/(1+t)^2,
and B_s,t(1)=(1+t^2-2st)/(1-t^2)^2.
The following entirely analytic rational enclosure suffices.

First r>25/28 and r<9/10 follow by squaring against r^2=4/5.
Use the elementary bounds 7/11<kappa<16/25. Then

    s< (16/25)(4/5)(9/10)=288/625<1/2,
    integral A_s,t dnu_r<511/722<71/100.

The square of the first term in (7) is therefore less than

    (81/100)(3/28)(71/100)=17253/280000<1/16,                 (9)

where the final comparison is 276048<280000. That term is <1/4.
Also

    s>(7/11)(4/5)(25/28)=5/11>4/9,
    integral B_s,t dnu_r<10461/722<15.

The square of the second term in (7) is less than

    (1/100)(16/25)15=12/125<(31/100)^2,                     (10)

since 960<961 after denominator 10000. It follows that

    B_flat(M_0;9/10)/(2n^(3/2))<1/4+31/100=14/25
                                      <2sqrt(2)/5.          (11)

The last squared comparison is 196/625<200/625. No parameter was
optimized or scanned and no arithmetic checker was executed.

## 5. Actual-sequence use and the remaining outlier condition

For any actual sequence of zero-original-source cells with a fixed
S=O(n^(3/2)) cap, delta->0, u->4/5, r->2/sqrt(5), and unscaled empirical law
nu->(delta_0+delta_(4/5))/2, equations (1)-(10), continuity on the
fixed t=9/10 window, and the recorded cell/padding errors prove

    limsup E max_cell X_z/(2n^(3/2))<=14/25<2sqrt(2)/5.       (12)

Here X_z is the base Gaussian cross process of the weighted-shell
theorem. The fixed trace cap controls the recorded padding error. The
unscaled weak-law convergence AND the actual radius convergence imply
nu_r->(delta_0+delta_1)/2; neither one is substituted for the other.
Thus the previously uncovered literal middle-law test is inside this
sharper upper region; no claim that every old r=1 metric fails is needed.
In particular its normalized field upper is strictly below sqrt(2)alpha
at alpha=2/5.

A weak empirical law alone does NOT control r: even one singular
outlier can keep r near one while disappearing from that law. The
finite theorem always uses the actual r and remains valid with such
outliers, but the evaluation (11)-(12) is not then automatic. No
actual source realizing the diagnostic data, automatic control of
outliers by conditional optimality, or bound on all other cells is
claimed. The required conditional paired-norm comparison and original
all-orders convergence remain OPEN.

## 6. Complete prerequisites and collaboration provenance

The author read these four complete sources before this extension:

- `original_mo_diagonal_majorizer_weighted_shell_upper.md`, 381 lines,
  SHA256 `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`;
- `original_mo_diagonal_majorizer_metric_stability.md`, 252 lines,
  SHA256 `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`;
- `original_mo_full_sdp_gap_weighted_compatibility.md`, 303 lines,
  SHA256 `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`;
- `original_mo_small_gap_pure_cross_upper.md`, 312 lines,
  SHA256 `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.

Root identified the enlarged actual-radius metric and the changed
coefficient s=a u r. The proof author derived the finite formula,
the r-uniform normalized compatibility/stability extension, and the
single-metric rational enclosure; root separately checked its
rational inequalities before this source was written. These are
contributing roles, not independent reviews of the combined source.
