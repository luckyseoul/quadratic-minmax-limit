# Independent frozen-source review: low-rank cross-spike mass upper

2026-09-06. PASS, no correction requested. This reviewer contributed
neither an argument nor a coefficient to the new source. Root, exact,
and proof workers' contributing roles are correctly disclosed in Section 7
and the author receipt; their checks are not independent whole-source
reviews. No mathematical execution or canonical change was performed.

## Complete sources and scope

Read directly and completely, with frozen hashes verified:

- Proof: `/tmp/original_mo_low_rank_cross_spike_mass_upper.md`, 347 lines,
  SHA256 `30347140ecf9fb2458444fb152490c601fe81d8a1733e90f31be692126ecdf1c`.
- Author receipt: `/tmp/original_mo_low_rank_cross_spike_mass_upper_author_receipt.md`,
  63 lines, SHA256 `e7c5cb7d680ec77c08af4c83b012f1e4ad7470f42abffa6c5bbe4347ab8966de`.

The entire 279-line actual-radius prerequisite was independently reviewed
earlier in this same review session and its unchanged hash rechecked:
`44fa3e7361e2142b20dce58d2dde727458db786529690f15e752390b8081725f`.
The full 381-line weighted-shell, 252-line metric-stability, 303-line
compatibility, and 312-line pure-cross prerequisites were also directly
read in this session, with hashes respectively `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`,
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`,
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`,
and `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
The earlier actual-radius theorem remains valid; only its separate
operator-radius convergence premise is removed by the new corollary.

## Finite allocation, projection, and Boolean remainder

For the actual weighted state zeta, orthogonal spectral decomposition
gives |v_z|<=t_z and |u-v_z|<=r0(1-t_z), hence the stated relaxed interval
I_u. No small Boolean spike-mass conclusion is drawn from low rank.
Projection of g alone costs at most sqrt(S tr(Pi Mhat)). The identities
sum q_i=N(1+delta), sum 1/q_i=N, and sum (q_i-1)^2/q_i=N delta give
sum |q_i-1|<=N sqrt(delta(1+delta)); together with M<=2wnI this proves
the precise 2n sqrt(w[d+N sqrt(delta(1+delta))]) loss.

The metric is positive with P>=E=eD. Direct expansion proves the exact
original radius in (2.3); only v_z<=t_z is relaxed. The PSD
Cauchy--Schwarz split uses g0-h=(P-E)P^(-1)g0. Since
P^(-1)g0=D^(-1/2)F D^(-1/2)g, the two expected squared quantities are
T-eR and sum Var(h_i)/d_i=e^2R. No constant projected diagonal or
covariance/resolvent commutation is needed.

For independent standard copies, rotation gives exactly
H(b)^2-b^2=E[(V^2-(U+sqrt(2)b)^2)_+]. Centered intervals maximize
Gaussian mass under translation, proving (3.1), including zero variance.
The constrained Fenchel shift is lambda sqrt(S d_i)a_i; weighted
Cauchy--Schwarz then yields the displayed one-variable bound.
Its infimum is sqrt(kappa S)e sqrt(R)sqrt(1-t), with the t=1 value
obtained by a limit. Correlated, unequal Gaussian marginals are allowed.

## Conditioning, nets, and finite mass bins

The identity h=K0 A0 J0 g and its split against hbar=A0g give (4.1):
the first term uses ||K0-I||_F and E||J0g||^2, the second uses
||x||_1<=sqrt(N)||x||_2 and E||(J0-I)g||^2. A0 is an operator
contraction, so hbar has covariance norm at most 2wn.
The weighted state direction admits the stated xi-net; its inner product
is at least (1-xi) times the actual projected norm. Selecting constrained
subsets uses Lipschitz constant 2n sqrt(w), then two uniform coupling
costs transfer the bound back to h. Empty subsets and d=0 are handled.

Within a nonempty mass bin, (2.4) and Jensen give the first term at t_+;
the constrained remainder uses t_-. Selecting bins uses the ACTUAL g,
whose covariance bound is known, not the transformed field. The identity
S tr(Mhat)=wnN^2(1+delta) bounds all normalized trace factors uniformly.
Thus endpoint and net relaxations cost O(sqrt(h0)) and O(sqrt(xi)).
Taking n first at fixed xi,h0 makes projection, conditioning, and net/bin
selection errors vanish; only then are xi,h0 sent to zero. This justifies
the uniform o(n^(3/2)) claim even when u,T,R vary. The old original drift,
padding, representative choice, and polynomial-cell errors remain intact.

## Weak-law corollary and its remaining hypotheses

At fixed r0>rho=2/sqrt(5), weak convergence of the complete squared law
forces rank(Pi)=o(n), without controlling the largest actual singular
value. M0=n(I-kappa u_D L) is positive, has diagonal n, and norm at most
2n. Its congruence error controls the bounded truncated resolvents.
Traces keep the ORIGINAL dimension denominator; the omitted mass vanishes.
The cutoff has no limiting atom, so weak convergence applies to these
bounded truncated traces. First n tends to infinity, then r0 decreases
to rho with inverse gap e=1/10 fixed throughout.

The exact original radius yields the limiting envelope
F0 sqrt(1+t/rho)+G0 sqrt(1-t), uniformly on [0,1]. Independently checked:
F0/rho<(1/4)(28/25)=7/25; s<1/2 gives IB>9461/722>13, hence
G0^2>91/1100>49/625, the last comparison being 56875>53900.
Its derivative at zero is negative and its derivative decreases, so
every allocation is bounded by F0+G0<14/25<2sqrt(2)/5.

The resulting base Gaussian cell upper is valid with arbitrary low-rank
spectral outliers, under the stated trace cap, delta->0, original
p=q_A=0, u_D->4/5, and weak middle law. It does not establish those
hypotheses for conditional optimizers, assert a realizing source, cover
all other internal-energy cells, or close original MO convergence.
