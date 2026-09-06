# Independent complete review: actual-cross-radius shell upper

2026-09-06. Verdict: PASS, with no mathematical correction requested.
This is an independent review of the new actual-radius extension, not an
author or contributing review. The reviewer supplied no derivation,
coefficient, metric choice, or rational enclosure to this source.
Root and the proof author have contributing roles disclosed in Section 6.
No mathematical execution, scan, optimization, census, or checker was run.
No canonical repository file, gate, backup, or publication was changed.

## 1. Complete frozen sources read and identified

The entire final source was read directly:

- `/tmp/original_mo_actual_cross_radius_shell_upper.md`, 279 lines,
  SHA256 `44fa3e7361e2142b20dce58d2dde727458db786529690f15e752390b8081725f`.

The following complete prerequisite sources were also read directly,
including their hypotheses, normalization, error terms, and scope.
Their hashes were checked in this review:

- `/tmp/original_mo_diagonal_majorizer_weighted_shell_upper.md`, 381 lines,
  SHA256 `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`;
- `/tmp/original_mo_diagonal_majorizer_metric_stability.md`, 252 lines,
  SHA256 `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`;
- `/tmp/original_mo_full_sdp_gap_weighted_compatibility.md`, 303 lines,
  SHA256 `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`;
- `/tmp/original_mo_small_gap_pure_cross_upper.md`, 312 lines,
  SHA256 `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.

One combined tool output truncated portions of the 252/303-line sources;
both were then reread completely in separate untruncated outputs. The
381/312-line sources were fully present in the original output. No
conclusion rests on a partial read or on prerequisite review summaries.

## 2. Finite original-shell inequality and comparison with the old metric

The same actual feasible diagonal D majorizes the cross-only matrix by
block-sign conjugation and averaging. Its normalized cross matrix L has
extreme eigenvalues +/-r with 0<r<=1; the positive radius follows because
the actual complete cross signing is nonzero. V=L/r has norm one.

For |t|<1, P=D-(t/r)H_B is positive, and E=(1-|t|)D satisfies 0<=E<=P.
The actual ORIGINAL cell cross energy gives exactly
z^TPz=S(1-tv_o), v_o=2c/(rS). The cross-only contraction yields |v_o|<=1,
so q-tr E=S(|t|-tv_o)>=0. Weighted representative energies have not been
substituted into this exact radius constraint.

The retained Boolean square-completion remainder gives the two terms
in (1). For E=bD its diagonal sum is at most b sqrt(SR), by weighted
Cauchy--Schwarz. The first trace is
tr(MP^(-1)(P-E)P^(-1))>=0. No covariance/metric commutation or independent
field-coordinate assumption is required; singular positive M is allowed.

At fixed old admissible eta, P is unchanged and h increases from
1-|eta| to 1-|eta|r. Writing Q=q/S, the derivative of
sqrt((Q-h)(T-hR))+sqrt(kappa)h sqrt(R) is at most
-(1-sqrt(kappa))sqrt(R) where the product is positive: apply AM--GM to
T-hR and R(Q-h). Continuity covers the boundary; R=0 is immediate.
The extreme eigenvalues +/-r make 1-|eta|r the maximal admissible h
in E=hD. The extension never worsens the old bound at an old eta.

The representative still lies inside its final original/weighted cell.
Thus the substitution in the full weighted-shell upper is valid for
every actual cell. Drift, padding, independent cushion, binning, and
cell-selection terms are not removed by this finite improvement.

## 3. Radius-uniform actual-diagonal stability

The key energy comparison applies dispersion to H_B/r directly:
H_B/r=D^(1/2)V D^(1/2), ||V||=1. Hence its discrepancy from dbar V has
Boolean quadratic norm at most S sqrt(delta). Because it is cross-only,
the same bound is its cross bilinear norm. Dividing by n dbar gives
|v_o-v_D|<=2sqrt(delta). This does not divide an r-independent old error
by r and is valid for arbitrarily small positive r.

For J=sqrt(dbar)D^(-1/2), positivity and diag M=wn give the three
Frobenius identities/bounds stated before (4). Factoring JMJ-M through
M^(1/2) proves its nuclear norm bound E_delta. Cyclic traces give
dbar T=tr(JMJ F) and dbar R=tr(JMJ F^2), without commuting J and F.
The scalar function (x-b)/x^2 on x>=b proves
0<=F-bF^2<=I/(4b), while ||F||<=1/b.

Both radius coefficients lie in [0,2]. Their difference is at most
2|t|sqrt(delta). The product error is consequently at most
E_delta/(2b)+|t|wN^2sqrt(delta)/(4b)<=wN^2sqrt(delta)/b
for 0<=delta<=1. Square-root continuity bounds its contribution by
sqrt(w)N delta^(1/4)/sqrt(b), before the common factor sqrt(N).
The second contribution is at most sqrt(kappa E_delta). Their sum is
bounded by the constant 3 in (2), uniformly in r and the actual cell.
The w=0 and delta=0 cases are valid directly.

The reference uses the same actual M,V and only the representative's
weighted value. It is a compared numerical functional, not an exact
weighted-bin ellipsoid constraint. Infimizing on a fixed compact
|t|<=1-b_0 preserves the all-cell error bound. No endpoint exchange with
n, or outlier-insensitive replacement of V, is used.

## 4. Pure-cross reference and exact rational evaluation

Original p=q_A=0 controls the internal weighted covariance contribution
through the already proved delta-only compatibility estimate. This part
does not require trace optimality: that hypothesis in the 303-line
source is needed for obtaining delta from its canonical SDP gap.
The finite positive-field comparison error is exactly the stated
2sqrt(k log(2))N^(3/2)delta^(1/4) for a non-augmented maximum.

M_0=wnI-k c_D L=wn(I-sV) has s=(k/w)ur, not (k/w)u/r.
Since |u|<=r<=1, its positivity and constant diagonal hold. For t>0,
F-bF^2=t(I-V)F^2. Pairing eigenvalues +/-sqrt(y) gives exactly the
displayed A and B numerators. The traces are 2wn^2t integral A and
2wn^2 integral B, respectively. Together with the radius coefficient
t(1-u/r), this proves every factor in normalization (7).

For r^2=4/5, u=4/5, the identities u/r=r and s=kappa r^3 are exact.
At t=9/10 the two half-endpoint means are (461+100s)/722 and
(18461-18000s)/722. These follow directly from A(0)=B(0)=1,
A(1)=(1+s)/(1+t)^2 and B(1)=(1+t^2-2st)/(1-t^2)^2.

Every rational comparison in (9)--(11) was independently checked
algebraically, without an arithmetic execution:

- 25/28<r<9/10 follows by squaring against 4/5;
- s<(16/25)(4/5)(9/10)=288/625<1/2;
- integral A<511/722<71/100;
- the first term squared is <17253/280000<1/16 because
  16*17253=276048<280000;
- s>(7/11)(4/5)(25/28)=5/11>4/9;
- integral B<10461/722<15;
- the second term squared is <12/125=960/10000<961/10000;
- the sum is <1/4+31/100=14/25, whose square 196/625 is below
  (2sqrt(2)/5)^2=200/625.

Only the elementary already established interval 7/11<kappa<16/25 is
used here. No prior Machin calculation or other certificate was rerun.

## 5. Sequence implication and limits of the result

For the centered-sign parameters w=1,k=kappa, the specified actual
zero-original-source sequence, delta->0, and fixed trace cap make all
retained comparison/padding errors negligible after normalization.
At fixed t=9/10 all reference integrands are continuous with denominators
bounded away from zero. The unscaled weak-law convergence together with
ACTUAL r->2/sqrt(5) gives the stated normalized endpoint law. Thus the
limsup upper <=14/25 and its strict comparison with 2sqrt(2)/5 follow.

The operator-radius condition is genuinely separate from a bulk law;
a vanishing proportion of singular outliers can keep r near one. The
finite all-cell bound remains valid with those outliers, but the displayed
diagnostic evaluation is not automatically available. No realizing
source, automatic conditional-optimality implication, all-cell closure,
conditional paired-norm theorem, or original limit theorem is asserted.
The original all-orders convergence question remains OPEN.
