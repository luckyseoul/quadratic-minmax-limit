## 2026-09-25: field-width potential controls all transformed normals

The non-Boolean coefficient-growth problem in the Banaszczyk recursion is
now resolved by a common slab potential.

For a transformed slab `|c.z|<b`, padded by zeros on eliminated
coordinates, define

`P_A(c,b)=b+(1/2)||A c||_1`.

Every original stable slab satisfies exactly `P_A=F+r`. For the general
compatibility child
`h=d_i c-c_i d`,
`beta=d_i b_c+c_i b_d-2c_i d_i`,
normalize by `g=h/2, b_g=beta/2`. If the parents have
`||c||_infty,||d||_infty<=1`, then so does the child, and

`P_A(g,b_g) <= [d_i P_A(c)+c_i P_A(d)]/2-c_i d_i`.

If both parent potentials are at most `L>=2`, the right side is at most
`L-1` for every nonzero pivot pair `0<c_i,d_i<=1`. Hence every genuine
compatibility generation drops the potential by at least one unit, while
zero-pivot slabs persist unchanged.

So arbitrary non-Boolean normals remain infinity-normalized and field-width
controlled at all depths. The explicit obstruction is now narrower:
this is an upper potential bound and does not yet force every generated
half-width to stay positive.

Proof:
`evidence/NOTE_2026-09-25_FIELD_WIDTH_POTENTIAL_RECURSION.md`.
Commit: `4c2c9bde89cd5b790bf0dc835d9aa2d5c6c32dca`.

Next live target: find a complementary LOWER estimate on generated
half-widths, or bypass scalar positivity through an exact section/density
argument.

## 2026-09-25: opposite-phase children have vertexwise bipartite dominance

The current Banaszczyk recursion had one explicit structural gap: opposite-
phase parent pairs produced child supports with only an aggregate size lower
bound. That branch now has a stronger inherited constraint.

Switch by the +stable parent x. If y is -stable and T is their disagreement
cut, write a_i for the signed degree of vertex i inside its own side and c_i
for its signed degree across T|T^c. Then, vertex by vertex,

`c_i >= |a_i|`.

Thus every signed cross-degree is nonnegative and dominates the magnitude of
the corresponding internal signed degree. If
`C=F-(D_x+D_y)/2` is the parent cross energy, then

`||B[T]1||_1 <= C` and `||B[T^c]1||_1 <= C`.

So opposite-phase children are not arbitrary unsupported principal states:
their support is one side of a vertexwise-dominant signed cut. This is
strictly stronger than the prior aggregate bound `|T|(n-|T|)>=C` and feeds
directly into the one-vertex transform route.

Proof:
`evidence/NOTE_2026-09-25_OPPOSITE_PHASE_BIPARTITE_DOMINANCE.md`.
Commit: `aaaa0da2bb2a2468b5175b1064b7b841b9b4d5b3`.

This does NOT yet close the recursion or convergence. The next live question
is whether vertexwise dominance is preserved strongly enough under the next
coordinate transform to control the resulting non-Boolean normals.

## 2026-09-25: core-signature pairing attacks the one-vertex target directly

The frozen one-vertex route now has a correlation-sensitive certificate that
does not pay the full reciprocal-square mass `Xi_A(r)`.

Choose a core `C` of active stable classes. Coordinates are grouped by their
projective sign signatures on `C`; pair coordinates inside each class so
that every pair cancels every core constraint identically. If `u` coordinates
remain unmatched and `u<r`, Komlos is applied only to the residual pair
differences. For a pair `p=(i,j)` with projective relation
`sigma_j=lambda_p sigma_i`, define

`M_p(r,u)=sum_(s outside C, x_i^s != lambda_p x_j^s) 1/(D_s+r-u)^2`.

If the pairing satisfies

`4 max_p M_p(r,u) <= 1/(18 pi)`,

then there is an ACTUAL incident sign row with extension increment strictly
below `r`. At `r=r_n` this gives `alpha_(n+1)<alpha_n`; the usual
Dini-summable excess version also follows.

This replaces the fatal global `Xi` charge by a bottleneck weighted
DISAGREEMENT charge after annihilating a selected dangerous core. It therefore
uses redundancy among stable constraints instead of treating them as
independent rows.

Proof and exact finite search formulation:
`evidence/NOTE_2026-09-25_CORE_SIGNATURE_PAIRING_EXTENSION.md`.
Commit: `733be9f82108fb3999208a6720a6b1a5b86bebb3`.

A diagnostic on the recorded exact order-15 minimizer explains why this is
needed: the neutral old certificate has 301 active antipodal stable classes
and `Xi_A(r_15)~=14.1414`, versus the required `1/(18 pi)~=0.01768`,
even though the exact best one-vertex extension only raises `Phi` from
27 to 30. The next live task is therefore to find/prove a small-core pairing
with uniformly bounded residual disagreement mass, not to bound `Xi`
itself.

# Handoff: original convergence problem

## 2026-09-24: new Komlós theorem hits the cross-order bottleneck

A literature search for September 2026 mathematics found the Guo--Fang--Lu
resolution of the Komlós vector-balancing conjecture
(arXiv:2609.11189, submitted 2026-09-10). Their dimension-free theorem gives
constant `K=3 sqrt(2 pi)`.

Apply it to the EXACT stable-state one-vertex constraints.  After identifying
antipodes and merging the two orientations, let `D_A([x])` be the stable
deficit and, for desired increment r, define

`Xi_A(r)=sum_{D_A([x])+r<=n} 1/(D_A([x])+r)^2`.

Weight the constraint [x] by `1/(D_A([x])+r)`. Every old-vertex column
then has Euclidean norm exactly `sqrt(Xi_A(r))`. Komlós therefore gives
one ACTUAL sign row balancing every active stable state simultaneously.
The exact theorem is

`Xi_A(r)<=1/(18 pi)  =>  Phi(extension)<Phi(A)+r`.

For an exact minimizer and the neutral increment
`r_n=m_n[(1+1/n)^(3/2)-1]`, this implies

`alpha_(n+1)<alpha_n`.

With an added `rho_n` satisfying
`sum rho_n/(n+1)^(3/2)<infinity`, the same theorem gives convergence.

The unconditional contrapositive is already informative: if
`alpha_(n+1)>=alpha_n`, then EVERY exact minimizer has
`Xi_A(r_n)>1/(18 pi)`. Since no term exceeds `1/r_n^2`, every nondecreasing
step forces more than

`r_n^2/(18 pi)=(alpha_n^2/(8 pi)+o(1))n`

active one-flip-stable antipodal classes within raw deficit at most `n-r_n` of
the edge.  So even a flat normalized step has a proved LINEAR near-edge degeneracy
cost.

This is qualitatively stronger than the old random-row union certificate.
At deficit zero and r=Theta(sqrt(n)), the random union criterion charges
constant probability per state and effectively tolerates only O(1)
constraints. The Komlós certificate charges Theta(1/n) per state and
tolerates a linear near-edge family before it can fail.

Proof:
`evidence/NOTE_2026-09-24_KOMLOS_ONE_VERTEX_EXTENSION.md`.
Main theorem commit: `efcc0257ba77b610709feadac3c04e6a3715b147`.

If alpha_n failed to converge, such nondecreasing steps would have to occur infinitely often, so this degeneracy obstruction would recur infinitely often. This is a genuine cross-order advance, but the reciprocal-square bound is
not yet proved for all minimizers, so convergence remains OPEN.

## 2026-09-24: zero-field ties remove the finite-depth wall

The previous actual-successor proof kept the old spin at zero local field,
which forced the universal count recursion down to `mu'>=p^2 mu`.
That convention was unnecessarily destructive. A zero field is a genuine
best-response tie, so choose the opposite spin instead.

For every original changed-or-tied coordinate i, the successor field is
independent of i's own flip coin. Conditional on all other coins, if
`a=x_i F'_i`, then successor disagreement-with-tie occurs with probability

`(1-p) 1_{a<=0}+p 1_{a>=0} >= p`.

At `a=0` the probability is one. Therefore, FINITELY and for every order,

`mu' >= p mu`.

All weighted identities are unchanged because zero fields carry zero
weight. The shifted-noise proof also survives unchanged, so

`g' >= 2p d_p mu^(3/2)-O_p(n^-1/2)`

and the interior count now propagates as `eta'>=p^3 mu^3-o(1)`, replacing
the old `p^6` loss.

This produces an arbitrary-fixed-depth scalar recursion for actual
non-Gaussian successors. If a state has
`e>=a-b alpha`, count floor M and disagreement floor G, then one update
with probability p gives

`a'=(1-p^2)a+p(1-p)G`,
`b'=(1-p^2)b+p^2`,
`M'=pM`,
`G'=2p d_p M^(3/2)`.

The deficit closure is
`[(1+b)alpha-a](alpha+G/2)>=G^2/8`.

There is also a strict continuation theorem. If x is the positive closure
root and `X=(1+b)x-a`, then choosing
`p_*=G/[2(2x+G-X)]` makes the NEXT mean certificate strictly exceed x.
The new count remains positive, so the process can be repeated.

Applied after the existing B_2 state, this gives

`B_2 < B_tie=x_2 < x_3 < x_4 < ... < B_infty <= 1/2`

and, because every rung is a fixed finite-depth theorem before taking the
scalar supremum,

`liminf alpha_n >= B_infty > B_2`.

Proof:
`evidence/NOTE_2026-09-24_TIE_RESOLVED_INFINITE_LADDER.md`.

This is the first current result showing that the actual repeated-update
lower-bound mechanism has NO finite-depth saturation: every certified finite
stage admits a strictly stronger finite successor stage. It does not yet
connect different matrix orders, so the original convergence problem remains
OPEN.

## 2026-09-24: second actual update gives another strict gain

The first ACTUAL successor already has a positive disagreement floor.
Applying the arbitrary-law mean-update inequality once more to that
non-Gaussian successor gives, for a second fixed update probability v,

`alpha >= T(t,u,v)=D/C`,

with
`C=1+v^2+(1-v^2)u^2` and
`D=(1-v^2)A+v(1-v)G`.

Define
`B_2=sup_{99/100<=t<=1,0<u<=1/2,0<=v<=1/2} T(t,u,v)`.
Then

`B_2 > B_opt > B_noise > B_both > B_int > B_tilt`.

The strict step needs no numerical optimizer. At a pair attaining B_opt,
write `x=B_opt`, `X=(1+u^2)x-A`. The first-step closure gives
`X=G^2/[4(2x+G)]`. The second-step slack at alpha=x is
`X-Gv+(2x+G-X)v^2`; its minimizing admissible v makes the slack strictly
negative. Thus one further actual update excludes alpha=B_opt.

Numerical locator only:
`t~=0.9934102, u~=0.0974342, v~=0.0046322, B_2~=0.3258669876`.
Proof:
`evidence/NOTE_2026-09-24_TWO_STEP_SUCCESSOR_GAIN.md`.

This is fixed-time repeated-update progress. A nondecaying long-time gain
or cross-order convergence mechanism is still OPEN.

## 2026-09-24: optimized successor envelope is strictly stronger

Verification pass on the merged September 24 branch checked the
both-phase residual-square identity, the shifted-noise moment bound, the
phase/state normalizations, the successor energy recursion, and the algebra
leading to `B_noise`. No contradiction was found in those derivations.

The next positive step removes two unnecessary historical parameter freezes.
For each fixed `99/100<=t<=1`, retain the exact both-phase disagreement
angle

`H(t)=min(1/4, atan(sqrt(beta(t))/C(t))/pi)>1/5`.

For an independent actual partial-flip probability `0<u<=1/2`, define

`G(t,u)=2u*2sqrt(u(1-u)/3)*H(t)^(3/2)`,

`A(t,u)=(1-u)^2 e(t)+u(1-u)f(t)`, and `c=1+u^2`.

The actual successor law then satisfies the exact asymptotic envelope

`liminf alpha_n >= L(t,u)`,

where

`L=(A-cG/2+sqrt((A+cG/2)^2+cG^2/2))/(2c)`.

Hence

`B_opt=sup_{99/100<=t<=1,0<u<=1/2} L(t,u)`

is a valid unconditional lower bound. It is STRICTLY stronger than the
merged `B_noise`: at the old parameters the exact angle already improves
the coarse `H>=1/5` replacement, and differentiating the new envelope
shows the old pre-successor update optimizer has positive one-sided
derivative. A numerical locator is approximately
`t=0.9934103, u=0.0974340, B_opt~=0.3258669873`; the theorem is the exact
supremum, not the decimal.

Proof:
`evidence/NOTE_2026-09-24_OPTIMIZED_SUCCESSOR_ENVELOPE.md`.
Original convergence remains OPEN.

## 2026-09-24 follow-up: both phases and actual successor flip noise

Two further analytic results, with author review only:

1. [Both-phase disagreement](evidence/NOTE_2026-09-24_BOTH_PHASE_DISAGREEMENT.md)
   retains the phase previously discarded. A square-completion identity
   gives a residual-variance floor in both phases, followed by convexity
   of the Gaussian disagreement angle. For the initial tilted source,
   `mu_0>=1/5-o_L(1)` and `eta_0>=1/125-o_L(1)`: an eightfold
   improvement of the earlier conservative eta floor. The interior-only
   bound becomes `B_both=[B+sqrt(B^2+p^2/(750*(1+p^2)))]/2`.
2. [Actual successor flip noise](evidence/NOTE_2026-09-24_SUCCESSOR_FLIP_NOISE.md)
   proves finite cap-free estimates after independently flipping each bad
   coordinate with probability `0<p<=1/2`. Writing `mu=average E[k]/n`,
   `g=f-2e`, `d_p=2sqrt(p(1-p)/3)`, and `E_p=d_p+4/3`,
   `mu'>=p^2 mu` and `g'>=2p d_p mu^(3/2)-2p E_p/sqrt(n)`.
   These concern the ACTUAL non-Gaussian successor, using conditional
   flip-coin independence and a shifted fourth-moment estimate.

Keep `t=993/1000`, its existing optimizing p, and `B=B_tilt`. Set
`G=2p d_p/5^(3/2)`, so `G^2=16p^3(1-p)/375`. The source energy
bound and sharp deficit closure applied to the actual successor prove

`liminf alpha_n >= B_noise = [B-G/2+sqrt((B+G/2)^2+G^2/(2*(1+p^2)))]/2`.

An exact rational comparison gives `B_noise>B_both>B_int>B_tilt`.
No new decimal evaluation is claimed. Same-order regularization removes
the initial fixed cap with the established order of limits. The two gap
corrections are compared, NOT added to the same successor energy.

There is now quantitative control at every FIXED number j of these
constant-p updates: `eta_j>=p^(6j)/125-o_L(1)` and
`g_(j+1)>=G p^(3j)-o_L(1)`, alongside the existing energy recursion.
The floors decay exponentially. They do not give a useful long-time
energy limit, state-dependent update control, or cross-order comparison.
Those are the next unresolved implications. Original convergence is OPEN;
the completion registry is unchanged.

New checkers, both UNEXECUTED:
`evidence/both_phase_disagreement_20260924/check.py` and
`evidence/successor_flip_noise_20260924/check.py`. Each is a serial SymPy
scalar check with 60-digit (not interval-certified) decimal evaluation.
Pinned hashes and a fresh-staging remote command are in
`evidence/successor_flip_noise_20260924/README.md`.
The NUKA SSH preflight failed with a socket permission error. No controller
mathematical test, independent review, or new major-milestone backup is
claimed. The user's earlier Orin `0` report covers only the OLD checker:
`evidence/paired_disagreement_floor_20260924/USER_REPORTED_CHECK.md`.
The research is published on
`research/paired-disagreement-floor-20260924` at commit
`bd5f2465ec40615d99da0622c26487b42bc496af`. The checker results remain
pending.

## 2026-09-24: the paired source activates the interior correction

Review of main at `889a4c88c8aa0634693e9b7db54de6db4fa7c18d` and all
eight requested September 23 commits:
`evidence/REVIEW_2026-09-24_CURRENT_GAP.md`.

The initial tilted Gaussian phases at `t=993/1000` have phase-averaged
changed-coordinate fraction at least `1/10-o_L(1)` and therefore
`eta>=1/1000-o_L(1)` for every fixed source cap L, with constants independent
of L. The new step uses the existing distinguished-coordinate joint
Gaussianization lemma and a scalar covariance comparison; it does not
assume a Gaussian post-update law. Keeping the already-optimized weight p
and writing B for the previous tilted lower constant gives

`liminf alpha_n >= B_int = [B+sqrt(B^2+p^2/(6000*(1+p^2)))]/2 > B`.

The exact expression was the lower bound at the first checkpoint; it is
now superseded above. Proof: `evidence/NOTE_2026-09-24_PAIRED_DISAGREEMENT_FLOOR.md`.
Review is analytic author review only. No check was run at that checkpoint:
the NUKA preflight failed to resolve its hostname. A subsequent user-reported
Orin scalar check and its limited scope are recorded in the receipt linked
above. This supporting lemma required no new major-milestone backup.

This supplies eta control for the INITIAL paired source and makes the
interior term strictly positive. Quantitative control of actual repeated
updates, or a neutral-increment extension with summable excess, remains
missing. The original convergence problem remains OPEN.

## 2026-09-23: stable-skeleton Hamming geometry

If x,y are stable in the same orientation and differ on T, `|T|=t`,
then every i in T obeys `w_i(x)<=2(t-1)`, and
`sum_{i in T} w_i(x)<=4 sigma Q_{A[T]}(x_T)`.
Thus the disagreement set is contained in the `2(t-1)`-light coordinates
of either endpoint.  A stability margin h forces Hamming separation
`d_H>=1+ceil(h/2)`, and the local number of stable states within radius r
is bounded by the corresponding binomial volume on the light-coordinate
set.  This gives actual packing structure on the exact state family left by
the one-vertex stable-state reduction.
Proof:
`evidence/NOTE_2026-09-23_STABLE_SKELETON_GEOMETRY.md`.

## 2026-09-23: one-vertex extension reduces exactly to stable states

For any proposed incident row `a`, the extension objective is
`max_{sigma,x}[sigma Q_A(x)+|a.x|]`.  If
`sigma x_i(Ax)_i<0`, flipping coordinate i raises the oriented quadratic
energy by at least 2, while the absolute row correlation can drop by at most
2.  Therefore the extension score never decreases under such an improving
flip.  Iteration terminates at a one-flip-stable state, giving exactly

`Phi(A+a)=max_{sigma,x sigma-stable}[sigma Q_A(x)+|a.x|]`.

Hence the random-row exact-tail criterion needs to sum only over stable
states, not the full cube.  This is an exact cross-order simplification,
not a heuristic pruning.  Proof:
`evidence/NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md`.

## 2026-09-23: sharp pointwise deficit-disagreement envelope

For any oriented Boolean state x, put
`D=alpha*n-Q(x)` and let
`R=sum_{i: sign(Mx)_i != x_i}|(Mx)_i|`.
Optimizing the exact mean-update interpolation parameter state by state gives

`R <= D + sqrt(2 alpha n D)`.

This is equivalent to the whole fixed-p inequality family, but stronger than
the previously recorded averaged adaptive closure because it is pointwise.
In particular an O(sqrt(n))-near-edge state has only O(n^(3/4)) weighted
best-response disagreement.  Averaging and Jensen give exactly
`f-2e <= 2(alpha-e)+2sqrt(2alpha(alpha-e))`, whose solved form is the
optimized mean-update bound already used for the tilted Gaussian theorem.
Proof:
`evidence/NOTE_2026-09-23_POINTWISE_DEFICIT_DISAGREEMENT.md`.

## 2026-09-23: cap-free weighted-disagreement closure

For any actual phase law, let `R` be the total local-field mass on
coordinates whose synchronous best response changes sign.  Independently
flip each such coordinate with the sample-dependent probability
`p=R/(4 alpha n)`.  Principal restriction bounds the internal flipped-set
energy, and polarization gives `p<=3/4`.  The resulting actual Boolean
update gains at least

`R^2/(4 alpha n)`

conditionally.  Averaging the two phases yields the new distribution-free
closure

`alpha >= e + (f-2e)^2/(16 alpha)`,

or `f-2e <= 4 sqrt(alpha(alpha-e))`.  This is valid for Gaussian, Gibbs,
or arbitrary post-update laws and does not use an operator cap or the
changed-coordinate count `eta`.  It therefore advances the repeated-update
bottleneck: weighted disagreement is now forced to vanish whenever the
source energy approaches the Boolean norm.  Proof:
`evidence/NOTE_2026-09-23_ADAPTIVE_BEST_RESPONSE_CLOSURE.md`.

## 2026-09-23: tilted paired phases improve the bound again

The September 17 Gaussian phases were not optimal at their fixed tilt
`t=1`.  Keeping `t` symbolic in
`R_s(t)=(I+s t M)^2/(1+t^2q)` preserves the paired-variance proof and
the same-order cap removal.  After optimizing the existing Boolean update,
the exact bound is

`B(t)=e(t)-f(t)/2+0.5*sqrt(f(t)^2+(f(t)-2e(t))^2)`,

with `e(t)=kappa*t/(1+t^2)` and
`f(t)=sqrt(kappa*(1-kappa*asin(t^2/(1+t^2))+kappa*t^2/(1+t^2)))`.
The envelope derivative satisfies `B'(1)<0`, so this is a strict
analytic improvement, not parameter noise.  At the explicit rational
choice `t=993/1000`,

`liminf alpha_n >= 0.3258530333538241...`

with update weight `p=0.0971801303973218...`.  This gains about
`5.63291603e-6` over the already optimized `t=1` theorem.
Proof:
`evidence/NOTE_2026-09-23_TILTED_PAIRED_FIELD_LOWER.md`.
Scalar regression:
`tests/test_tilted_paired_field_lower.py`.
Original convergence remains OPEN.

## 2026-09-23: optimized paired-field lower bound

The September 17 paired-field theorem proved the cap-independent limiting
mean-update inputs but froze the update weight at `p=1/10`. Optimizing the
same proved inequality exactly gives

`liminf alpha_n >= B_* = 0.3258474004377944...`

at
`p_*=0.097102496812...`, an unconditional increase of
`6.6449822202e-6` over the previous bound. The closed form is

`B_*=e0-f0/2+0.5*sqrt(f0^2+(f0-2e0)^2)`,

with `e0=kappa/2`, `f0=sqrt(kappa*(2/3+kappa/2))`,
`kappa=2/pi`. Same-order regularization removes the fixed cap exactly as
before. No new premise, family assumption, or finite census is used.
Proof: `evidence/NOTE_2026-09-23_OPTIMIZED_PAIRED_FIELD_LOWER.md`.
Original convergence remains OPEN.

## 2026-09-23: finite interior gap and mean-update correction

`evidence/NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md` gives the cap-free finite
bound `Phi(A)-|Q_A(z)| >= epsilon^2 k^3/(96 Phi(A))` for `k>=2`
epsilon-interior coordinates. Its short product proof retains both
one-sided extrema: `max Q_B * max(-Q_B) >= k^3/96`. Principal restriction
and symmetric perturbation inside the cube then give the interior gap.
The layer-cake corollary is
`sum_i(1-|z_i|) <= 1+3[96 Phi(A)(Phi(A)-|Q_A(z)|)]^(1/3)`.

For the existing mean update at `0<=p<=1/2`, let k_sigma count its changed
coordinates (keep the old sign at zero local field), and let
`eta=average_sigma E[k_sigma^3 1_(k_sigma>=2)]/n^3`. The finite inequality is

`(1+p^2)alpha >= (1-p)^2 e+p(1-p)f+p^2 eta/(24 alpha)`.

The added term is new in this update bound. The older asymptotic cubic
one-sided estimate in the September 6 archive is explicitly reused as
prior context, not claimed as a new pairing route. The new direct proof
works at every subset size without a spectral cutoff. It does not supply
the archived missing linear posterior estimate or a one-vertex row that
controls all near-maximizers. No lower-bound constant was retuned.

NUKA ran the new serial algebra check once: 15 checks passed, exit 0.
Full command, hashes, staging path, raw result, and author-review scope
are in `evidence/boolean_interior_gap_20260923/`. There is no live job.
Only local algebra was computationally checked; the analytic theorem is
author-reviewed, without independent human or proof-assistant review.
This supporting result required no new major-milestone backup.

Next unresolved implication: quantitative control of `e`, `f`, and `eta`
under actual repeated updates, or an independent cross-order comparison.
Equation (9) alone can be used for non-Gaussian laws, but gives no field
dynamics. The later September 24 successor note adds a separate finite-time
field estimate; long-time control is still missing. Original convergence
remains OPEN; the canonical lower bound is the latest value recorded above.

## 2026-09-23: coherent excess is pinned at the doubling constant

For every signing, \(m_{2n}\le\Psi(A)+n\), by the Seidel block
\([[A,A+I],[A+I,-A]]\). Combined with \(\liminf\alpha_n>0\), this forces
\[
\limsup\Psi(A_n)/m_n\ge 2\sqrt2
\]
for exact minimizers, and \(\Psi\le 2m_n+Cn\) fails infinitely often for
every fixed \(C\). Claude `deep_review` (`claude-opus-5-5`, effort `max`)
returned PASS-WITH-NOTE, `do_not_branch`. Since \(m_n=\Theta(n^{3/2})\),
\(\limsup\Psi/m_n\le 2\sqrt2\) is the same statement as
\(\Psi\le 2\sqrt2\, m_n+o(n^{3/2})\). That upper bound is sufficient for
the doubling estimate and is not proved. Summability of the error is a
stronger demand and still does not, alone, close existence.

OpenAI `suggest_direction` (`gpt-6-astra`) on the next step: attack the
Paley gaps. Prove \(\gamma_q\le\gamma_p+\varepsilon(p)\) for every prime
\(q\ge Q(p)\), with \(\varepsilon(p)\to0\) independent of \(q\). A liminf
seed then forces the Paley limsup and liminf to agree, hence existence,
without a value. It said not to treat the \(\Psi\) rate or the factor-3
ray alone as sufficient, because a fixed-factor bound with a vanishing
error can still oscillate. Claude's direction call at Opus 5.5 / max
timed out at 680s and returned no verdict.

Restriction gives only the local half of that comparison. For odd primes
\(p<q\), with \(n=p^2+1\) and \(N=q^2+1\),
\[
\gamma_q\le\gamma_p+\tfrac12\Bigl(
\sqrt{1-1/N}-\sqrt{1-1/n}\,(n/N)^{3/2}\Bigr).
\]
A fixed rise \(\delta\) forces \(N/n\ge 1+(4/3)\delta+o(1)\). Consecutive
Paley gaps therefore satisfy \(\limsup(\gamma_{p_{k+1}}-\gamma_{p_k})\le0\).
The same majorant tends to \(1/2\) as \(q/p\to\infty\), so monotonicity
of \(m_n\) does not prove a tail bound \(\varepsilon(p)\to0\) independent
of \(q\). Note: `evidence/NOTE_2026-09-23_PALEY_GAP_MODULUS.md`.
Test: `tests/test_paley_gap_modulus.py`.

OpenAI `deep_review` (`gpt-6-astra`) of the closure strategy returned
BLOCK on one sentence: an increment \(\delta_n=\Theta(\sqrt n)\) does not
by itself produce a non-summable \(\alpha\)-error. The neutral increment
is \(m_n((1+1/n)^{3/2}-1)\). The strategy is the one-vertex extension
bound \(\delta_n\le m_n((1+1/n)^{3/2}-1)+r_n\) with
\(\sum r_n n^{-3/2}<\infty\). That forces \(\alpha_n\), hence \(\gamma_p\),
to converge. S1 and S3 remain sufficient side routes. Doubling alone does
not.

The finite form \(\Psi\le 2m_n+2n\) is already false. The matrix
`evidence/psi_minimizer_excess_20260923/A8_psi40.txt` has \(\Phi=10=m_8\)
and \(\Psi=40\). A census of the \(2^{21}\) first-row-positive switching
representatives found 4200 classes at \(\Phi=10\), none with
\(|Q-2Q_T|>20\). Note:
`evidence/NOTE_2026-09-23_PSI_MINIMIZER_EXCESS.md`.
Test: `tests/test_psi_minimizer_excess.py`. Limit OPEN.

## 2026-09-21: actual-Gibbs covariance floor

The author-reviewed supporting proof is
`evidence/NOTE_2026-09-21_ACTUAL_GIBBS_COVARIANCE_FLOOR.md`.
Its conditional-score identity gives the all-directions matrix bound
`Cov(X) >= D^2/[1+||J||op+2 max_i sum_j J_ij^2]`,
`D=diag(E sech^2(h_i))`. Combining it with actual edge-optimality
controls both half-product phases and the unsigned symmetric Gibbs law.
Whole-row replacement plus a verified density comparison extends the
matrix floor to each actual single-row cavity tilt. The older scalar
row-response result is explicitly reused, not republished as new.

NUKA ran the new serial symbolic check once: twelve identities passed,
exit 0. Full command, input hashes, result, and author-review limits are in
`evidence/actual_gibbs_covariance_floor_20260921/`. No signing census,
unchanged test suite, or old family scan was run. Unfinished pre-existing
untracked work was left untouched. This supporting lemma did not require
a new major-milestone backup; prior checkpoints remain preserved.

Next unresolved implication: a lower bound on the signed internal/cross
pressure-defect integral. The new LOWER covariance bounds do not give an
upper fluctuation bound or control a growing collection of pinned rows.
Keep all operator-norm dependence and optimizer hypotheses. No original
convergence result or improved lower-bound constant is claimed.

## 2026-09-21: correct unsupported asymptotic exclusions

The September 19 operator upper bound does not exclude a sub-half Boolean
norm for plus-I. Its valid conclusion is only that the coarse certificate
`(sqrt(2)||B||op+1)/(2sqrt(2m))` is always at least 1/2 for `m>=2`.
Larger operator norm weakens this certificate, not necessarily the actual
construction. The three construction gaps 2,4,6 do not prove an all-orders
formula or bound the optimal gap from above. CORE §8 and the Paley notes
now distinguish an upper bound on a construction from a lower bound on
the optimum. No closed family scan is restarted.

The September 20 finite shell fixtures cannot exclude an asymptotic
comparison with a Dini error: any finite exceptional set can be covered
by a nonincreasing finite-support remainder. The Monte Carlo threshold
grid also does not certify an infimum over every real threshold. These
scope errors are corrected in the shell note. Its numerical receipts
are preserved, not reclassified as exact asymptotic counterexamples.

These are proof-status corrections, not a new convergence theorem.
The original goal remains convergence, and the paired-field lower bound
is unchanged. Historical entries below are superseded where they make
the corrected exclusions.

NUKA verification: the one changed scalar test
`test_minimal_op_input_gives_certificate_above_half` passed by direct
invocation on its three retained scalar inputs; exit 0. The unchanged
random-signing test and other suites were not rerun. Fresh staging:
`/tmp/mo-scope-check-20260921.Rd9fTlnz`. Staged/local test SHA-256:
`5d61c2a3e8e23c922abc2f96c6c92a3aa000b3974d6a68d35d2edfdc9d7c7b86`.
The general certificate inequality and finite-support Dini argument are
proved in the corrected notes; this small regression does not prove them.
This supporting correction is not a major milestone or destructive edit,
so no new large-drive snapshot was taken. Prior backups remain intact.

## 2026-09-19: direction 2 — discrete-derivative superlevel (started)

2026-09-20 (gpt-6-astra exchange run): numerical evidence challenges the
zero-error energy-shell bound \(W\le 2\sqrt2\Phi\) on recorded exact
minimizers (C5, Paley6, plus-I C5, \(B_{13}\)); the widths were estimated,
not certified exactly. Linearized midpoint \(\mathbb E\mathcal M\)
tracks \(W\) and is also above target. The shifted-threshold samples
approach \(\Psi(A)=\Phi([[A,A],[A,-A]])\) at \(n=5,6\) (12 and 18),
still above \(2\sqrt2\Phi\) on the tested grid. This is not a certified
real-threshold infimum. Φ-drop from a global minimizer is empty.
Identity \(\Psi=\max|2Q-4Q_T|\) stands. The slack
\(|Q-2Q_T|\le\Phi+n\) is exhaustive through order 7 and **false at
\(n=8\)** (plus \(K_4\), rest minus: \(|L|=28>\Phi+n=24\)). On that
signing the centered threshold beats \(2\sqrt2\Phi\) in samples, and
\(h\to\infty\) is worse. Doubling is an \(\inf_h\) at finite \(h\), not
\(\Psi\). `evidence/NOTE_2026-09-20_PSI_COHERENT_PAIRING.md`. Limit OPEN.

OpenAI path-choice: superlevel structure, not Fourier/cut-norm and not
flag algebras. Proved D1–D14. D10 uses a maximizer of \(Q\), not \(|Q|\);
write \(M_+=\max Q\). D14: at \(p=5\), \(m=13\), every Q-maximizer with
\(M_+=22\) has a minus triangle, so \(\Gamma\ge6\ge4\) and
\(\Phi(K)\ge63\ge61=\Phi(C_{26})-4\). Proof: nauty geng enumerated all
163477 unlabeled triangle-free \(\Delta\le6\) graphs with \(e\ge28\);
zero satisfied the Q-max cut inequalities. OpenAI math_review PASS
after a BLOCK on the false \(|Q|\)-maximizer form. Remaining Case B:
\(M_+\le20\) at \(m=13\) (need \(\Gamma\ge8\); the known block attains
equality) and all \(p\ge7\). Limit OPEN.
Note: `evidence/NOTE_2026-09-19_SUPERLEVEL_DISCRETE_DERIVATIVE.md`.
Tests: `tests/test_superlevel_discrete_derivative.py`. Limit OPEN.

## 2026-09-19: exact Phi(K_50)=169 on nuka RX 9070 XT

hipBLAS gfx1201 full projective cube (\(2^{49}\) states, 8311 s):
\(\Phi(K(B_{25}))=169\) exactly. Paley \(C_{50}=175\). Certified undercut
by \(6=p-1\), so \(m_{50}\le169\) and \(\gamma_7\ge0.01697\). Same binary
smoked n=10 \(\Phi=13\) and n=26 \(\Phi=61\). This is still
a finite construction gap, not an asymptotic classification. Receipt:
`evidence/k50_nuka_gpu_20260919/`. Limit OPEN.

## 2026-09-19: Paley-gap attack — spectral scope corrected September 21

The missing inequality is still unproved. The estimate
\(\|K(B)\|_{\mathrm{op}}\le\sqrt2\|B\|_{\mathrm{op}}+1\)
gives \(\Phi(K)/n^{3/2}\le\tfrac12+o(1)\) for near-min-op B;
this upper bound does not exclude an actual smaller limit. A larger
operator norm only weakens the certificate. The three certified undercuts
fit \(p-1\), without proving that formula. Known beaters have 2–3 odd local-field
values and \(A^2\) off-diagonals in \(\{0,\pm4\}\). Matrix Aut (AΓL with
square multiplications) has commutant dimension 2 on \(V_+\), so Schur
does not prove the Max+ frame identity. CORE §4 \(u_{ij}\)-correction is
\(<10^{-3}\) on the n=26 beater. OpenAI `suggest_direction` (after one
timeout): fourth-order non-Gaussian rounding on \(D=A^2-(n-1)I\);
Gibbs interpolation; affine-line amplification (not excluded by the
displayed spectral upper bound). These are suggestions, not new dispatch
authority or a restart of closed families. Note:
`evidence/NOTE_2026-09-19_PALEY_GAP_ATTACK.md`. Tests:
`test_plus_i_opnorm_bound`, `test_beater_magnitude`,
`test_paley_Vplus_commutant`. Limit OPEN.

## 2026-09-19: plus-I lift is the n=10 and n=26 Paley-beater

K(B)=[[B, B+I],[B+I, -B]]. K(5-cycle)=m_10=13. K(m_13-minimizer)=61,
which *is* the n=26 record (not a Paley matching). Paley-F_25 lift to
n=50 has ||K||_op=√61 and **exact** Φ=169 vs Paley 175 (nuka GPU).

## 2026-09-19: why γ_p is not closed

Gaussian CORE §4 is saturated at conference (u_ij=0). Lipschitz plus
small F is false as a route: the n=26 Φ=61 undercutter has switching
distance 122 from Paley, A^2 off-diagonal in {0,±4} with six nonzeros
per row. Max+ can be killed while Φ rises (160 flips, cube max 107).
The remaining object is an almost-conference family in the Paley-op
shell with ρ bounded away from 1. See
`evidence/NOTE_2026-09-19_BEATER_SHELL.md`.

## 2026-09-19: cluster set is Paley n=p^2+1

CORE §8 (corollary of §5 + ρ=1): existence of lim α_n is equivalent to
convergence of γ_p = (Φ(C)-m_n)/n^{3/2} on n=p^2+1. If γ→0 then L=1/2;
if γ→c>0 then L=1/2-c; two cluster points of γ would prove divergence.
Certified: γ_3=0.06325 exact; γ_5≥0.03017 from m_26≤61. No uniform c>0
construction, no o(n^{3/2}) gap proof. Do not reopen gap-2 covers.
Next: a construction with γ≥c>0 on this sequence, or an o(n^{3/2}) bound
on Φ(C)-m_n along it — not another 2n residue census.

## 2026-09-19: n=9 zero-error two-half diamond fails

The directed-half-cut minimax `min_R B(A,R)` on an exact `m_9=12`
witness is **36**. CP-SAT with 86 workers: `B<=32` and `B<=34` are
INFEASIBLE (~3 min each); ILS attains 36. Target `2√2·12≈33.941`, so
the zero-error form of (6.13) fails at n=9. Orders 7 and 8 still pass
(existing two-half geometry, not rerun). Normalized excess sequence
`+0.42,+0.26,−0.19,−0.01,+0.08` at n=5..9 does not stay negative.
Dini `n^{3/2}` still covers the gap of 2. See
`evidence/NOTE_2026-09-19_TWO_HALF_N9_ORIENTATION.md`.

## 2026-09-19: residue (6.20) finite census; honest X cards

Uncommitted work checked on soulkiller and nuka. Soulkiller main had the
two 2026-09-18 cross-transfer handoffs, the reverted-but-kept 2026-09-15
lower-bound note (leave untouched), and `check_real_domain.py`. The
`mo_*.py` sweep scripts were in `/home/nick/`, not the repo; they are now
under `evidence/cross_transfer_20260918/`. Nuka main matched those two
handoffs only. Nuka `quadratic-minmax-limit-leftover-work` still has dirty
k7/p13 Paley-orbit files on `codex/leftover-moment-attack` — E1/residual
route, not reopened.

Finite census of Prop 6.6 residue (6.20): on exact minimizers
`m_5=4`, `m_6=5`, `m_8=10` the residue is 55%/91%/51% of folded pairs and
the Paley-skew `R` misses the diamond by `0.42/0.53/0.34` in units of
`n^{3/2}`. High-α random signings have sparse residue and typically hold.
See `evidence/NOTE_2026-09-19_RESIDUE_620_FINITE_CENSUS.md`. This is not a
doubling proof. Next implication: an `A`-dependent bound of `x^T R y` on
(6.20) for optimizers, or emptiness for large n (already false at n=5,6,8).

X-thread cards in `x-cards/` no longer claim `L=1/2`. Regenerated 2026-09-19
from `scripts/render_x_twopager.py` and `scripts/render_x_lemmas.py`.

OpenAI referee `suggest_direction` on (6.20) timed out once this session
(Codex 180s). Do not treat that miss as a Claude consult.

## 2026-09-17: new cap-independent paired-field lower bound

The new proof `evidence/NOTE_2026-09-17_PAIRED_POLYNOMIAL_FIELD_LOWER.md`
gives `liminf alpha_n >= B = 0.3258407554555742... > 13/40`, with
`B=[81 kappa/2+9 sqrt(kappa(2/3+kappa/2))]/101`, `kappa=2/pi`.

New mathematics: split the actual polynomial-phase sign covariances as
`C_s=C_0+s 2kappa M/(1+q)+E_s`, with `||E_s||_F=O_L(1)`, retain the
PSD even tail, and pair the two local variances using
`(M^3)_ii^2<=q[(M^4)_ii-q^2]`. The resulting field lower is independent
of fixed L. Reusing scalar Gaussianization and the already-completed
Boolean-norm mean update gives the displayed B. Existing same-order
regularization removes L by taking n to infinity before the cap.

Root reviewed the complete analytic chain. NUKA ran the new nine-check
exact scalar/polynomial certificate once, all passed, exit 0. No signing
or family census was run. Proof SHA-256:
`1f4ddcc4eb7f1098ccc480165834fbd784b763ab1ecc7734d0361d7698383a3c`.
Checker SHA-256:
`73fda9a3a7c8dd87831e9eec40cb22e6ce1ca958e7e10c3fe3250cfefb4840e9`.
The receipt and detailed author-review scope are in
`evidence/paired_polynomial_field_20260917/`. This is not independent
human review or formal verification of the entire proof.

Nick clarified that the September 17 revert was requested to remove
duplicate work, not to reject the lower-bound result. Do not interpret
its commit-message paraphrase as a mathematical retraction or permission
to rerun completed work. The prior untracked September 15 note/checker
remain untouched; the old update is an input, not a new contribution.

Convergence remains OPEN. The next unresolved implication is a genuine
cross-order comparison, such as a lower bound on the integrated optimized
pressure defect; this theorem does not supply one. Do not spend a new
campaign retuning the last decimal of B. A separate all-real-domain
checker was drafted for Soulkiller but not run (Z3 absent in default
Python); it is outside this commit and not evidence. No owned mesh
research process was left running.

Verified major-milestone checkpoint (before publication):
`/mnt/storage/backups/codex/quadratic-minmax-limit-2212b99a4c24/20260917T053401Z-ec5a643658a3-paired-field-global-lower.CFqifY`.
It covers HEAD `ec5a643658a3c8cb74174679fd4f31f483d34aa4` plus the full
uncommitted tree, including the new proof and all untracked work. Actual
mount: `/dev/sda2`, ext4, read-write; bundle, archives, source comparisons,
unchanged refs/status, and checksums passed. SHA-256 of `SHA256SUMS`:
`7ffc45fde203f19bd89b438e85d7af907699371c171bc549746ee3c877e89556`.
NUKA result SHA-256:
`73b0b787df8def4686675cd4804563b62c33843545a32b794fe74a2731212c32`.
This receipt paragraph was added after the snapshot. No earlier backup was
overwritten or deleted.

## 2026-09-12: session summary — closures, m_38 <= 109, nuka ops

**Conditional spectral bridge diagnostic (finite only):** for the fixed
order-four `Phi=4` source, exhaustive enumeration of all 65,536 cross
signings gives conditional value 10 and 184 conditional minimizers, all with
`||B||op <= 2 sqrt(2)`.  The controller and Soulkiller JSON receipts agree
(SHA-256 `0b23ae0975aa6063282f097b0221fb95a9feb876703f18d910b7d47ba22e354f`).
This is evidence for, not a proof of, the unproved asymptotic conditional
spectral bridge; see `evidence/conditional_spectral_bridge_probe_20260912/`.
For the `Phi=4` order-five negative-cycle source, an 88-worker Soulkiller
CP-SAT solve proves conditional value 13; one returned minimizer has
`||B||op=3.3722813232690143` and independently replays on the controller.
This is one finite witness only; see `evidence/conditional_spectral_bridge_n5_20260912/`.
For the order-six `Phi=5` source, the complete 36-bit, 4,096-pair
88-worker CP-SAT model proves conditional value 18; its returned witness has
`beta(B)=14` and `||B||op=1+sqrt(5)`, independently replayed on controller.
Again this is one finite witness only; see `evidence/conditional_spectral_bridge_n6_20260912/`.
For the order-seven `Phi=9` source, the complete 49-bit, 16,384-pair
88-worker CP-SAT model proves conditional value 21; its returned witness has
`beta(B)=21` and `||B||op=2 sqrt(3)`, independently replayed on controller.
The n=5--7 normalized one-witness operator ratios are 1.508, 1.321, and
1.309; descriptive only, not an asymptotic bridge. See
`evidence/conditional_spectral_bridge_n7_20260912/`.

**New finite records and closures since 2026-09-11:**

- `m_16 <= 30` (tie intermediate; `m_16 in [27,30]`). The order-16 completion
  "witnesses" claiming 30 were RETRACTED (zero off-diagonal entries; not
  signings); clean-seed completion floor is 32. Two distinct valid `Phi=30`
  order-16 objects exist (tie intermediate; nuka16 source). The S15
  one-vertex gap is closed: exact optimum 30, reproducing the tie
  intermediate entrywise via a different code path.
- `m_19 <= 39` (one-vertex extension of the 18-winner; witness landed and
  cross-verified on both GPUs).
- **`m_38 <= 109`** — explicit Paley conference evaluation `Phi(C38)=109`,
  improving the spectral bound 115.58 and the coherent-lift 121; verified on
  V100 and RX 9070 XT. Conferences: C14=21, C18=33, C26=65 (record 61
  stands), C30=75 (=record), C38=109. Two-block theorem verified at q=25, 37;
  the `block = m_k` pattern breaks beyond k=9 (blocks 30>20, 47>39).
- **Exhaustive diagonal-lift closures**: order 16 min 32 over all 256
  diagonals (12 attain), order 18 min 39 over 512 (46), order 20 min 40 over
  1024 (2 attain — family exactly closed at its record). Orders 28/30/32/36/38
  sampled: `D=I` optimal-so-far everywhere (order-38 anchor 121 = the archived
  coherent-lift value, cross-validated).
- **One-vertex scans 20->24**: 46, 44, 49, 53, 56 — records 21-24 matched
  EXACTLY by extensions of the 20/21/22/23 winners ("tight chain"
  42->44->49->53->56; the better order-20 object (40) extends to 46 —
  extension quality is not monotone in source value). Bank winners are
  DISTINCT objects from the scanned extensions (same values). The 24->25 scan
  (bank 60) is the open shot.
- ILS floors at ALL records 16-21 under the linear-response machinery; the
  wide-GEMM reformulation gives 100x-class passes (order-36 pass: 2-5 s).
  Tools landed: `scripts/hunt.py`, `scripts/lift_opt.py`,
  `scripts/exhaustive_lift.py`, `scripts/comp_sweep.py`,
  `scripts/paley_conference_eval.py`.

**Nuka ops (important):** the `amd_cupy` wheel is ROCm 7.2.4-soname-locked
(cannot use the installed ROCm 10.0 at `/opt/rocm/core-10.0`; no ROCm-10
cupy wheels are published). Full-array cupy reductions (max/sum/argmax) are
100-450x slow on gfx1201 under the 7.2 stack — use two-stage row reductions
(`av.reshape(-1,4096).max(axis=1).max()`) or hipcc/hipBLAS (ROCm 10: wide
sgemm 0.766 ms/call). The V100 remains ~2x faster than the 9070 XT on these
GEMM shapes. Turn wall-clock budget raised to 7200 s in
`~/.codewhale/config.toml` (backup `config.toml.bak-20260912`).

**Status:** the convergence question is OPEN as before (bracket unchanged);
the finite frontier is saturated under all current methods. Best known:
16:30, 17:32, 18:33, 19:39, 20:40, 21:44, 22:49, 23:53, 24:56, 25:60,
26:61, 28:70, 30:75, 32:80, 36:108, 38:109.

**Historical repo-state note:** the statement below reflects the September 12
session only. It is not a current publication instruction; commits on the
active research mission are normally pushed under the standing agreement.

## 2026-09-11: structural Paley lift family closed as an asymptotic route

The [family study](evidence/paley_structural_lift_family_20260911/README.md)
answers the two questions the earlier structural-lift artifact left open.
A field-generic `GF(q)` construction reaches the fourth and last enumerable
order `q=9` (order 20, exact `Phi=40`), and all four exact points fit
`Phi=N(N+12)/16` exactly. The fit is false: `K^2=(2q+1)I+2(C(+)C)` gives the
exact all-orders bound `Phi(K)<=(q+1)sqrt(2q+1+2 sqrt q)`, violated by the
fit at `q=25` — the first order past the four it was read off. Normalized,
that bound decreases to `1/2` from above, so the family never certifies below
`1/2`; mesh witnesses for every prime power `q=1 mod 4` up to `q=197`
(order 396) give `Phi>=0.468 N^(3/2)` with no downward trend, so it cannot
certify below about `0.47` either. Six independent deep shards agree on
`Phi=222` at order 60 against the fit's 270.

Do not re-run this family for a better constant, and do not restore the
"scalable mechanism" framing: the ceiling and the floor are now both known.
The useful transferable object is the `K^2` identity, which applies to any
`B=C+I` cross block over a conference `C`. Soulkiller, jellyfin, NUKA and
orin all ran; no GPU was needed.

## 2026-09-11: cone/correction integration on actual signings

The [finite mesh experiment](evidence/fixed_repair_path_20260911/README.md)
integrates the imported cone solver and correction-cycle driver with exact
Boolean norm acceptance. Cone-ranked integral pairs escape two of ten
single-edge stalls (22 to 20), outperforming two fixed proposal controls;
only one input reaches 18 under every method. Separately, exhaustive
repair scheduling excludes a no-overshoot bridge for the archived fixed
two-repair witness: optimum normalized peak is 17/32. Soulkiller did the
parallel work; NUKA independently replayed the paths and corrections.
No all-orders comparison follows. Do not rerun these same inputs as a new
search; the missing implication is uniform defect/repair control with order.

The separate full order-4 canonical-gap diagnostic in
`evidence/fixed_repair_path_20260911/` finds no numerical zero-gap cross
completion for an exact source minimizer: observed minimum relative gap is
1/42 across all 65,536 blocks. It is an SDP numerical diagnostic, with raw
primal residuals and independent exact structural replay, not a theorem.

The imported cone/correction pattern was then applied directly to that
same finite numerical majorizer residual, not merely to a Boolean-norm toy.
Across all 65,536 cross blocks, cone-ranked integral repairs improved the
relative gap at 65,046 starts (mean reduction 0.0772162), ahead of the
index-pair (64,754; 0.0712843) and seeded random-pair (64,776; 0.0725753)
controls. Soulkiller used 88 workers and NUKA independently replayed all
262,144 reported endpoints and aggregates. The complete compressed output
and exact scope are in `evidence/fixed_repair_path_20260911/majorizer_gap_repair_README.md`.
This is finite floating-SDP method evidence only: no zero-gap block, exact
certificate, optimizer premise, or all-orders implication follows.

Updated 2026-09-06. Start with `CORE.md` and `STATUS.md`.
The original MO limit is OPEN; `L=1/2` is also OPEN.
The reviewed bounds now satisfy `1/pi < liminf alpha_n <= limsup alpha_n <= 1/2`.

## Low-quota stop, 2026-09-06 09:17 UTC

The user reported only 3% quota remaining. New hunts were stopped; preserve
the [growth checkpoint](evidence/original_mo_growth_checkpoint_20260906/README.md)
and its archive before resuming. It records the exact two-extension method,
finite source-dependent traps, the exact Phi64=248 matrix certificate,
reverse-entry exclusions, and the final exact two-old-edge repair witness.
The proposed delayed-repair path DP was NOT launched. No global convergence
argument was obtained. Do not restart a broad maze or treat finite progress
as a demonstrated route to solving the original problem.

## Latest checkpoint: broad compute campaign

Research commit `902539250598fd21e78a1d9a1f6dbf0233ed0f16` was pushed to
`origin/main` and verified remotely. The one major large-drive backup is
`/mnt/storage/backups/codex/quadratic-minmax-limit-2212b99a4c24/20260906T053702Z-902539250598-strict-floor-broad-campaign.bzktfR`.
It covers that exact research HEAD, the strict-1/pi theorem, all four raw
campaign packages, and the full repository/history. Archive comparison,
bundle verification, checksums, and unchanged-source checks passed.
The backup receipt SHA-256 is
`30d03c8dfaeb295218ee5ef1f603c21d801e21ac0ca9b4829d7d20accb265ede`;
see [publication.json](evidence/original_mo_broad_campaign_20260906/publication.json)
for complete provenance. This navigation/receipt update is later metadata,
not an additional research result. The research snapshot's focused status
and documentation gate passed 38 tests in 2.11 seconds.

The [campaign handoff](evidence/original_mo_broad_campaign_20260906/README.md)
and its [hash/count manifest](evidence/original_mo_broad_campaign_20260906/manifest.json)
preserve the final two-GPU/100-worker run, started from
`6bcd7d46f307d0e55141244daa3e58044f32cf55`. At least 28,555,014 exhaustive
candidate-score evaluations were recorded; these are calls, not globally
distinct signings. All owned search jobs completed naturally.

Four self-contained archives include sources, seed snapshots, retained raw
candidate data, profiles, witnesses, exact replays, and explicit provenance
limitations. In particular, recoverable CPU-v1 counts are lower bounds,
heuristic spin maxima are not upper bounds, and rejected kernel diagnostics
are not accepted results. Independent field-update enumeration verified
norms 61 at order 26 and 80 at order 32. The campaign README gives the exact
replay command and explains why the new finite data do not imply convergence.
Unpublished analytic scratch in the independent archive is not a reviewed
theorem artifact; the last reviewed global strengthening remains the strict
1/pi lower at the starting commit. Do not restart unchanged censuses.

The next target is still an all-orders convergence or nonconvergence
argument, through any valid route. In the optional amplification route,
sampled lifting performance and an additive-n finite audit supply no uniform
near-minimizer comparison or Dini error. The existing coherent counterfamily
already excludes an all-source version based only on order-scale bounds.

## Preservation and reset

The previous residual worktree, including all 48 uncommitted files, is
preserved at `archive/2026-09-05-paley-research`:
`ad8c6920412af0b3c23629afe2a9e95060c5471e`.
The separate 22-file dirty main snapshot is preserved at
`archive/2026-09-05-main-local-edits`:
`c2e13218cceb7e1fb36de8f2625bf4c4a7c0a606`.
The active checkout is `/home/nick/quadratic-minmax-limit`, branch `main`.
See `ARTIFACTS.md` for the exact scope and replay of historical documents.

The canonical entry documents no longer enforce Paley residual (ii),
gap-two optimality, a conjectured value, or a particular amplification
construction. The original-problem status is separate from route-local
predicates. Valid local results retain their stated scope and evidence.

## Next mathematical work

The [cap-free source milestone](evidence/original_mo_cap_free_source_global_gap_milestone.json)
preserves two reviewed original-source extensions. Section 8 of the
[all-law tail theorem](evidence/NOTE_2026-09-06_ALL_LAW_SECOND_MOMENT_TAIL_GAIN.md)
proves an unspecified eventual uniform gap above 1/pi for ALL complete
symmetric zero-diagonal signings. It proves neither convergence nor a
numeric gap; F_1(1) bounds only the excluded saturation sequence.
For M=A/sqrt(n), the same note gives `liminf alpha(A)>=2/5+7/55000`
if EVERY fixed C'>5/3 has only o(n) eigenvalues with |lambda(M)|>C'
and `limsup tr|M|/n<=4/5`. No operator cap or full limiting law is needed.
The [cap-free near-flat theorem](evidence/NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_CAP_FREE_GAIN.md)
retains `5kappa/8+16/3125>2/5+3/1100` from the FULL law
`(9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4))` alone, kappa=2/pi.
The gain uses the full original source and Boolean rounding, with fixed
smoothing first and n->infinity before smoothing is removed. Do not infer
a fourth-moment bound, a uniform tail rate, or covariance control for the
unsmoothed update. Preserve prior records; global order comparison remains OPEN.

The new [cross-spike/bulk upper milestone](evidence/original_mo_cross_spike_bulk_upper_milestone.json)
contains an [arbitrary-profile ORIGINAL-target region](evidence/NOTE_2026-09-06_WHOLE_PROFILE_CROSS_BULK_UPPER.md).
Use ACTUAL K=[[A,B],[B^T,-A]], positive diagonal D with D+/-K>=0,
S=tr D=O(n^(3/2)), and `delta=S tr(D^(-1))/(2n)^2-1->0`.
Evaluated cells have ORIGINAL internal energies p=q_A=0, common cross
energy c, and representatives chosen within the final refined cells.
Put W=D_L^(-1/2) B D_R^(-1/2), u_n=2c/S, m_n=tr(W^TW)/n.
If `liminf u_n>=7/8` and, for every fixed R in (9/10,1), only o(n)
singular values of W exceed R, then the centered-sign base process obeys
`limsup [E max_cell X_z-2sqrt(2)Phi(A)]/n^(3/2)`
`<=2/sqrt(5)-2sqrt(2)/pi<0`.
Completeness supplies limsup m_n<=1/2; CORE's actual lower
`Phi(A)>=n sqrt(n-1)/pi` supplies the ORIGINAL target. No alpha=2/5 premise,
full limiting spectral law, diagonal optimality, or top-radius convergence
is assumed. This is a bulk condition: low-rank operator outliers remain.
The [actual-radius metric](evidence/NOTE_2026-09-06_ACTUAL_CROSS_RADIUS_SHELL_UPPER.md)
is valid for every actual cell. The [spike-mass extension](evidence/NOTE_2026-09-06_LOW_RANK_CROSS_SPIKE_MASS_UPPER.md)
retains the state's spike mass inside the Boolean remainder and pays for
Gaussian projection and weighted conditioning; its allocation theorem
also applies to general weighted fields under its low-rank/dispersion premises.
For the zero-source weak middle law (delta_0+delta_(4/5))/2 and x^TWy/n->4/5,
it gives limsup E max_cell X_z/(2n^(3/2))<=14/25 without top-radius control.
The wider liminf u_n>=4/5, limsup m_n<=2/5, bulk<=9/10 region is an alpha=2/5
diagnostic comparison, not a uniform ORIGINAL-target result below that scale.
All recorded cell/padding errors remain. The trace cap, dispersion, bulk,
and high-u premises are not derived for arbitrary conditional optimizers;
other internal-energy cells and original convergence remain OPEN.

The separately classified [actual central-cell boundary](evidence/NOTE_2026-09-06_ACTUAL_CENTRAL_CELL_LINEAR_WIDTH_BOUNDARY.md)
now shows why improving the CENTERED all-cell linear-field bound cannot
finish this route: actual central cells saturate normalized field width
sqrt(kappa*w), for any feasible D, without dispersion or trace-cap premises.
This is NOT a lower on the Gaussian cross process. At source scale alpha,
a shifted-sign field-upper argument necessarily needs w<=2alpha^2/kappa;
centered w=1 fails for alpha<=1/2. Next address biased central/active-cell
competition or sharpen the cross-process comparison, not another centered
linear-field metric. The scoped upper corollaries above remain valid.

The newer [all-law source-gain milestone](evidence/original_mo_all_law_source_gain_milestone.json)
contains the [uniform adaptive gain](evidence/NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md).
For ACTUAL complete symmetric zero-diagonal A, the two caps
`limsup ||A||op/sqrt(n)<=5/3`, `limsup tr|A|/n^(3/2)<=4/5` imply
`liminf Phi(A)/n^(3/2)>=35/88+3/1250=2/5+7/55000`.
No limiting spectral law or diagonal homogeneity is assumed; the update
probabilities adapt to actual local fields on the same Boolean source.
The [two-moment transfer](evidence/NOTE_2026-09-06_TWO_CROSS_MOMENT_SOURCE_NUCLEAR_TRANSFER.md)
uses ACTUAL K=[[A,B],[B^T,-A]], positive diagonal D=diag(D_L,D_R)
with D+/-K>=0, and `delta=tr(D)tr(D^(-1))/(2n)^2-1->0`.
Put W=D_L^(-1/2) B D_R^(-1/2), Y=WW^T, and retain the FULL moments
`m_D=tr(Y)/n`, `Delta_D=tr[Y(I-Y)]/n`.
If EVERY accumulation point (m,Delta) belongs to
`R={9/25<=m<=1/2, 0<=Delta<=m[4sqrt(m)-3sqrt(1-m)]^2/25}`,
a common large original principal source supplies both caps above;
its norm comparison transfers the same strict lower back to A.
No full cross or internal law, extra trace cap, optimal diagonal, or
active-face premise is required. A simple sufficient condition is
`liminf m_D>=2/5`, `limsup Delta_D<=1/1600`, with delta->0 retained.
The [endpoint-law transfer](evidence/NOTE_2026-09-06_CROSS_ENDPOINT_SOURCE_NUCLEAR_TRANSFER.md)
is preserved as the Delta=0 special case. The older near-flat theorem
below keeps its larger gap under its narrower law hypothesis.
The ACTUAL source operator cap and paired diagonal-dispersion/moment
premises are not established for arbitrary optimizers; the paired field
is not replaced. No mathematical run was used. Other source regions,
the all-cell implication, and original convergence remain OPEN.

The latest reviewed actual-source milestone (2026-09-06) is
`evidence/original_mo_original_source_strict_gain_milestone.json`.
Its full proof is
`evidence/NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md`
(SHA-256 `7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`).
For ACTUAL complete symmetric zero-diagonal signings A, assume
`limsup ||A||op/sqrt(n)<=5/3` and the FULL empirical eigenvalue law
of A/sqrt(n) tends to
`(9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4))`. Then, kappa=2/pi,

`liminf Phi(A)/n^(3/2)>=5kappa/8+16/3125>2/5+3/1100`.

This is an ORIGINAL same-order quadratic-norm lower bound. It uses a
robustly normalized positive spectral projector, uniform higher-chaos
mean variance at least 1-kappa-o(1), and trace-of-square first-chaos
alignment. The Gaussianization extension proves a joint limit of ONE
local field and ONE distinguished Gaussian input coordinate, uniformly
over rows; it does not claim a growing-dimensional joint field limit.
An actual independent-coordinate Boolean update with FIXED probability
1/10 improves the positive phase. Its penalty retains the actual 5/3
operator cap, not the limiting nonzero atom 5/4. Weak empirical flatness
does not imply exact finite-order flatness, a large exact kernel, or
absence of spectral outliers. The distribution-free scalar support is
`evidence/NOTE_2026-09-06_ORIGINAL_SOURCE_LOCAL_UPDATE_SCALAR_GAIN.md`
(SHA-256 `7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`).

The separate cap-free transfer is
`evidence/NOTE_2026-09-06_NEAR_SCALAR_INTERNAL_FLAT_LAW_TRANSFER.md`
(SHA-256 `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`).
For the ACTUAL paired K=[[A,B],[B^T,-A]], take ANY positive diagonal
D=diag(D_L,D_R) with D+/-K>=0. If
`delta=tr(D)tr(D^(-1))/(2n)^2-1->0` and the FULL actual law of
H_L=D_L^(-1/2) A D_L^(-1/2) tends to
`chi_*=(9/25)delta_0+(8/25)(delta_(3/4)+delta_(-3/4))`,
one common original principal A_J, q=|J| with q/n->1, has the
law at 0 and plus or minus 5/4 for A_J/sqrt(q),
limsup ||A_J||op/sqrt(q)<=5/3, and satisfies
`Phi(A_J)/q^(3/2)<=Phi(A)/n^(3/2)+o(1)`.
Completeness and the full second moment force dbar/sqrt(q)->5/3,
dbar=tr(D)/(2n); no separate trace cap or trace optimality is needed.
No cross-law or active-state condition is required. The auxiliary
source only lower-bounds the original norm and never replaces the
paired covariance, cross block, or active field.

Composition excludes the specified ACTUAL near-scalar internal-law
regime with Phi(A)/n^(3/2)->2/5 underlying the strengthened formal
profile below. The older formal certificate boundary remains correct
for its explicitly listed relaxation, which did not include this new
original-source entry constraint; no actual signing was supplied there.
The following older checkpoints are retained as history. Neither small
delta nor chi_* is proved for arbitrary candidate extremizers. Other
actual profiles and the remaining all-cell implication still require
work; unchanged trace-only scans do not establish those implications.
No mathematical run was used. Original convergence remains OPEN.

The latest paired milestone is
`evidence/original_mo_weighted_cross_gain_boundary_milestone.json`.
It completes the weighted transfer previously marked unpublished in
the older scalar-gain checkpoint below, which is retained as history.
Its two proof sources are
`evidence/NOTE_2026-09-05_NEAR_SCALAR_CROSS_SPECTRAL_GAIN.md`
(SHA-256 `ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9`)
and `evidence/NOTE_2026-09-05_STRENGTHENED_TRACE_PROFILE_BOUNDARY.md`
(SHA-256 `903ac72c78c60706fbcfef09e50abeda0a18fe05976e3efab89d65becdbfccf1`).

For the ACTUAL paired signing K, ANY positive feasible D is allowed
with the SEPARATE cap S=tr D<=C N^(3/2), N=2n. Put dbar=S/N,
delta=S tr(D^(-1))/N^2-1, and let nu be the FULL actual squared-
singular-value law of W_D, including zeros. With m=integral y dnu,
v_2=integral y^2 dnu and g_kappa=sqrt(kappa)-kappa, kappa=2/pi,
uniformly for 0<=delta<=1/512 the actual original cross norm satisfies

`beta(B)/(n dbar)>=kappa v_2/m+g_kappa m`
`-[25kappa C^2+6g_kappa]delta^(1/3)-R_C(n)`, R_C(n)->0.

Here m>=1/(2C^2). No trace optimality, small canonical gap, maximum
diagonal bound or global unweighted operator cap on B is assumed.
A balanced complete submatrix is only an auxiliary ORIGINAL-norm
lower bound; interlacing and congruence return its second and fourth
singular powers to the full original nu. It never replaces W_D,
the covariance, source, or active cell. The uniform marginal CLT
prerequisite and its tail envelope are retained; no finite-n rate is claimed.
For the SEPARATE actual active conditions p=q_A=0 and c=Phi(K),
the same lower holds for u_D=c_D/n with an additional 2sqrt(delta)
loss. Thus near-flat full weighted laws exclude u_D=kappa+o(1)
in this fixed-cap, delta-to-zero branch, without the earlier operator cap.

The paired FORMAL diagnostic then tests a changed profile:
alpha=2/5, f=4/3, u=4/5, m=9/25 and
nu=(16/25)delta_0+(9/25)delta_1. Its specified full/source/cross
conditions, INCLUDING the new entry gain, all pass. Nevertheless the
same reference functional has U_s(t)>71/125>2sqrt(2)/5.
Adding the full ORIGINAL drift z f/2 keeps the certificate above
the target for EVERY signed metric and shifted Gaussian threshold,
including endpoints. One global supporting-line inequality supplies
the squared target margin 41/15625; the kappa enclosure was reused.

This is a lower bound on the FORMAL UPPER certificate, not on actual
Gaussian width or a Boolean norm. No complete signing or actual active
optimizer realizing the profile is supplied. The transfer theorem is
not retracted, and the earlier formal obstruction retains its scope.
A next attempt must add source-entry, Boolean-active-state, frame or
conditional-optimality information, or change the upper argument;
unchanged trace-only threshold/metric scans do not supply that step.
Small delta for optimizers and the all-cell original upper remain open.
No new mathematical run was used. Original convergence remains OPEN.

The latest actual-entry restriction is
`evidence/original_mo_complete_cross_flat_spectral_gain_milestone.json`
and `evidence/NOTE_2026-09-05_COMPLETE_CROSS_FLAT_SPECTRAL_GAIN.md`
(SHA-256 `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`).
For an ACTUAL n by n complete sign matrix B, choose the SEPARATE
unweighted operator bound d>=||B||op and put
`m=n/d^2>=m_0>0`, `epsilon=1-tr[(B^T B)^2]/(n^2 d^2)`.
Uniformly in these data as n tends to infinity, with kappa=2/pi,

`beta(B)/(nd)>=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-o_(m_0)(1)`.

Here 0<=epsilon<=1-m. Exact flatness of all nonzero singular values
at d is equivalent to epsilon=0; exact or asymptotic flatness therefore
forces a strictly positive leading-order gain over kappa for fixed m_0.
The proof retains the actual complete-entry higher-Hermite variance gain
and supplies its necessary Gaussianization bridge: mixed contractions,
a Gamma fluctuation bound and a characteristic-function equation prove
a uniform scalar absolute-moment lemma. Uniform second moments and an
L2 Hermite-tail estimate justify the limit passage. This is not an
absolute-moment inference from variance alone or a growing-dimensional
joint column CLT. No finite-n error rate or mathematical run is claimed.

For the SEPARATE actual active conditions p=q_A=0 and c=Phi(K),
one has c=beta(B); at scalar scale d the same lower applies to u=c/(nd).
Consequently the earlier FORMAL flat cross endpoint with u=kappa is
excluded in this actual scalar, bounded-operator setting. The formal
trace-relaxation certificate obstruction remains valid on its own stated
relaxation; it is not an actual-signing counterexample.

The unweighted d>=||B||op and m>=m_0 hypotheses are NOT obtained here
from small dispersion delta or trace control of a diagonal majorizer.
Transfer to the actual near-scalar weighted cross law is the next
unpublished implication, being treated separately. This note neither
replaces W_D nor evaluates every active cell, and original convergence
remains OPEN. Preserve the older scoped proofs and their exact premises.

The latest actual coupling and formal certificate boundary are recorded in
`evidence/original_mo_source_cross_trace_boundary_milestone.json` and
`evidence/NOTE_2026-09-05_SOURCE_CROSS_NUCLEAR_TRACE_BOUNDARY.md`
(SHA-256 `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`).
For ANY positive feasible D, let N=2n, dbar=tr(D)/N,
delta=tr(D)tr(D^(-1))/N^2-1, alpha=Phi(A)/n^(3/2), and let nu be
the ACTUAL squared-singular-value law of W_D, with mean m. The finite
source/cross inequality is

`integral sqrt(1-y)dnu >= [kappa sqrt(m)/(2alpha)](1-1/n)/(1+delta)`
`-sqrt[(2N/dbar^2)(2delta+delta^2)]`, kappa=2/pi.

This needs no trace cap, optimal diagonal, small canonical gap, or active
cell. Its O(sqrt(delta)+1/n) consequence keeps the original matrices;
the direct nuclear transfer requires no auxiliary trimming.

The SAME note also identifies a FORMAL limitation of the listed trace/block
relaxation. At alpha=2/5, f=4/3, u=kappa and
nu=(1-m)delta_0+m delta_1, m=9kappa^2/16, the retained full/source/cross
moment inequalities and block contraction all pass. Nevertheless EVERY
shifted Gaussian sign threshold and signed ellipsoid metric, including
endpoint limits, has formal drift-plus-certificate value above the target
sqrt(2)alpha. The full ORIGINAL drift z f/2 is retained, with
z=|2Phi_Gauss(h)-1|; the certificate's noise term is bounded below by
sqrt(1-z^2)L_0, L_0^2=40501/125000>8/25.

These formal parameters are not realized by complete signings or actual
active Boolean states in this note. A lower bound on this UPPER certificate
is not a lower bound on actual Gaussian width or the original norm.
The next step within this route is an additional actual entry/active-state
constraint, or a different upper argument, beyond the listed trace data.
Do not repeat threshold/metric scans on the unchanged relaxation. The pi
enclosure was reused analytically; no new mathematical run was needed.
Original convergence remains OPEN; this is not a signing counterexample
or an impossibility theorem for all methods.

The latest delta-only normalization is
`evidence/original_mo_delta_normalization_milestone.json` and its proof
`evidence/NOTE_2026-09-05_NEAR_SCALAR_DIAGONAL_SPECTRAL_NORMALIZATION.md`
(SHA-256 `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`).
It permits ANY positive feasible D with the SEPARATE cap
`S=tr D<=C N^(3/2)`: neither trace optimality nor small canonical gap
is required. With dbar=S/N, r=(N-1)/dbar^2 and mu_j the j-th absolute
moment of the ACTUAL full weighted contraction T, it proves

`2Phi(K)/S>=max{kappa(1+r)mu_3/(2r),kappa r/mu_1}`
`-O_C(delta^(1/3)+N^(-1/2))`, with `mu_2=r+O(delta^(1/3))`.

The proof's principal signing is only an auxiliary ORIGINAL-norm lower
bound. Interlacing and congruence transfer its moments back to the FULL
actual T, including exceptional coordinates. It never replaces the
source, covariance or W_D. The finite phase and nuclear inequalities
are (4.3) and (5.2); no mathematical execution was needed.
For EXACT scalar D, mu_3=r(1-gamma), so the phase term remains useful
at positive gamma where the previous masked-gap estimate was vacuous.
That identity is not silently imposed on near-scalar D with outliers.

The same lower transfers to u_D=c_D/n only on the separately active
original face p=q_A=0, c=Phi(K)>=0, at error 2sqrt(delta). Small delta
is still a hypothesis, not established by source or conditional
near-minimality. The next use must couple actual full and cross spectra
at the actual norm scale, rather than enlarge scalar D and recover the
old attenuation loss. Smaller-normalization and all-cell width estimates
remain open, as does the original convergence problem.

The latest evaluated small-gap package is
`evidence/original_mo_small_gap_evaluation_milestone.json`:

- `NOTE_2026-09-05_FULL_SDP_GAP_ORIGINAL_PHASE_BOUND.md`
  (SHA-256 `1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`)
  constructs two ACTUAL unit-diagonal PSD Gaussian phases with the same
  coordinate normalization. Subtracting their ORIGINAL quadratic energies
  gives the coefficient kappa/2 directly, not through beta(K)<=4Phi(K).
  Here D is trace-optimal for the FULL SDP, S=tr D=tau(K), q=N-1,
  and gamma=(S-tr|K|^3/q)/S. Under a fixed original norm cap,
  `Phi(K)>=kappa S/2-O(N^(3/2)sqrt(gamma)+N^(5/4))`.
  Its finite bound retains the actual weighted mask loss. It is vacuous
  for gamma>=1/4; do not reuse it as a complementary large-gap gain.
  The conclusion `u=c_D/n>=kappa-o(1)` additionally needs the ACTUAL
  active pure-cross conditions p=q_A=0 and c=Phi(K), with N=2n.
- `NOTE_2026-09-05_SMALL_GAP_PURE_CROSS_UPPER.md`
  (SHA-256 `035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`)
  evaluates the actual
  squared-singular-value measure of W_D, including zeros. The exact
  sign-square identity and compatibility give
  `m=integral y dnu(y)=u^2/f_n^2+o(1)`, where c=f_n n^(3/2).
  For standard centered signs, at fixed t=3/5, concavity of A and B's convexity bound
  EVERY actual measure by the algebraic expression (3.5); no Dirac law
  or common extremizing measure is assumed. For
  f_n tending to sqrt(2), gamma tending to zero and the active conditions
  above, the resulting cell upper is
  `limsup E max X_z/(2n^(3/2))<=17677/25000<1/sqrt(2)`.
  The monotonicity argument covers all u>=kappa up to vanishing errors.

Both proofs and their complete review provenance are in the manifest.
Exactly eleven new fixed Fraction comparisons passed one soulkiller run;
result SHA-256
`0ea064435322e698b8e33a4d9bce8ab29156e3cfe013c9885f1f35e205156e41`.
The squared strict margin is `23671/625000000`. The earlier pi enclosure
was reused, not recomputed. The replay artifact is
`python3 evidence/original_mo_small_gap_pure_cross_rational_certificate.py`;
do not rerun unchanged arithmetic merely for another receipt.

This evaluates a previously unevaluated ACTUAL-law diagnostic face, not
the supremum over all coupled original/weighted cells. The general
formula retains f_n; the desired bound is still `F<=2sqrt(2)Phi(A)`.
If Phi(A)/n^(3/2) is below 1/2, the f=sqrt(2) result alone does not meet
that smaller target. The next implications are a bound at the actual
smaller normalization, control of nonzero original internal energies,
and a genuinely complementary positive-gap argument. A norm cap or
optimizer label does not supply small canonical gap. Original all-orders
convergence remains OPEN; none of these route-specific premises is a
necessary condition imposed on every possible convergence proof.

The latest quantified compatibility package is
`evidence/original_mo_gap_compatibility_milestone.json`:

- `NOTE_2026-09-05_FULL_SDP_GAP_WEIGHTED_COMPATIBILITY.md`
  (SHA-256 `3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`)
  uses the canonical primal of the LITERAL complete K and ANY attained
  trace-optimal same-diagonal majorizer D. With S=tr D, q=N-1 and
  g=S-tr|K|^3/q, weighted residual squares are at most 4qg. The exact
  squared Frobenius norm of the inverse-weighted commutator is
  `2(S tr(D^(-1))-N^2)`, giving
  `delta<=4Sg/(qN^2)`. Uniform cube rescaling then proves
  `Phi(K-(S/N)T)<=S sqrt(delta)` for T=D^(-1/2)KD^(-1/2).
  The ORIGINAL source energy errors are at most 2N sqrt(delta), and
  the cross error at most N sqrt(delta), after scaling by S/N.
  No maximum-diagonal bound, nonsingular K, or unique optimum is needed.
  Its original-zero-source corollary compares two individually PSD
  Gaussian fields at O(N^(3/2)delta^(1/4)) cost. The pure-cross field
  retains actual W_D,c_D; its width is NOT evaluated by this corollary.
- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_METRIC_STABILITY.md`
  (SHA-256 `ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`)
  applies to ANY actual majorizer D, not only an optimum. Constant field
  diagonal and a nuclear-norm congruence estimate control both exact
  resolvent traces and their cancellation. For 0<=delta<=1,
  `|B_D-B_flat|<=3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|eta|)`.
  Its (5.1) transfers the all-ACTUAL-cell upper on each fixed compact
  eta window, with the original drift and old bin/selection errors.
  B_flat is a NUMERICAL reference, retaining actual PSD M_theta and
  actual contraction L_D. The representative's weighted c_D need not
  be constant through a bin, so this is not an exact scalar-I shell
  constraint throughout that bin. Every representative is chosen within
  its FINAL refined original/weighted cell.

Both paths above are under `evidence/`. Complete author, root and
independent reviews are recorded in the manifest. These two analytic
proofs needed no mathematical execution, signing census or optimization;
the previous arithmetic certificates were reused without rerunning them.
No source signing, src module, test or global predicate was changed.

For bounded `S/N^(3/2)`, a relative canonical gap g/S tending to zero
gives delta tending to zero. This is a CONDITIONAL actual regime, not a
property established for every exact original or conditional minimizer.
Do not replace K by a purported contraction K/(S/N), remove its rare
diagonal outliers without accounting for them, or invoke an indefinite
scalar tensor covariance. The metric window must be fixed before its
asymptotic limit; the error is not uniform at |eta|=1.

The live next implications are an evaluated actual weighted trace upper
in the small-gap range and a correctly normalized full-K ORIGINAL-norm
argument for the complementary range. Rectangular beta bounds can lose
an essential factor through beta(K)<=4Phi(K); they cannot silently be
read as a quadratic-norm rounding gain. Neither these route-specific
targets nor the new compatibility bounds prove original convergence.

The latest same-source package is
`evidence/original_mo_weighted_covariance_milestone.json`.
Its four analytic results are:

- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_CROSS_COVARIANCE.md`
  (SHA-256 `0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`)
  chooses D with D+-K positive for the LITERAL complete block signing.
  The exact cross correlation is `R_D=I+Q(A tensor A-S_B+I)Q`,
  `q_ij=1/sqrt(d_i d_(n+j))`; its operator norm is below three.
  A norm cap supplies `tr D=O(N^(3/2))` and local correlations O(1/N)
  without trimming or changing A,-A,B. The complete weighted Hermite
  decomposition gives an O(n) retained-profile cost; the resulting
  conditional ORIGINAL-norm Gaussian floor has O(n^(16/11)) error.
  The separated even series is used only when epsilon<=1/2; the
  remaining bounded orders are treated separately. Dropping only the
  independent Gaussian variance padding costs O(n).
- `NOTE_2026-09-05_DIAGONAL_MAJORIZER_WEIGHTED_SHELL_UPPER.md`
  (SHA-256 `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`)
  proves exact positive weighted linear fields. Actual representatives
  and width-1/n bins handle real weighted energies with O(sqrt(n))
  comparison error; selecting all original/weighted cells costs
  O(n sqrt(log n)). The same D majorizes the ORIGINAL cross form H_B,
  so D-eta H_B has exact shell radius `tr D-2eta c` using unweighted c.
  Equation (4.6) is the full all-actual-cell upper, with exact two-trace
  field width (4.5). Its weighted feedback c_D and internal p_D,q_D
  remain distinct from the original drift `(p-q)/2+s c`.
- `NOTE_2026-09-05_SCALAR_TEMPLATE_CUBIC_ALIGNMENT.md`
  (SHA-256 `60037f67234fbca8c17ee90bf52c7f4346b24e5f18eb5f2c922ebbd2d9382c2a`)
  proves `j_3>=-1+2s^2/mu` for ACTUAL matched optimal frames.
  For finite `tau(C)=p q`, this strengthens the Gamma UPPER-certificate
  exclusion to `q>=12/5 => Gamma(C)>283/200>sqrt(2)`.
  It crosses q=1+sqrt(2), the particular weak-Dirac diagnostic barrier.
  Section 6 separately stipulates leading energy sqrt(2), derives its
  variable-u expression, and takes n to infinity at fixed eta before
  eta tends to one. It does not infer actual Boolean saturation from
  a Gamma cap or prove arbitrary nonsymmetric attainability.
- `NOTE_2026-09-05_TENSOR_DEFLATION_FIXED_CAP_RATE.md`
  (SHA-256 `22febfa722afb3e18878f23f8e140895da90a3eb41fe0179356b08232d44f27a`)
  constructs actual complete sources under each FIXED norm cap C>1/2.
  Their tensor positive-part Gaussian repair has a uniform-rate lower
  of order n^(3/2)/sqrt(K), for both tensor signs and for the symmetric
  zero-diagonal quadratic norm. The amplitude is fixed before K, then
  n grows. No assertion is made for C=1/2, exact source minimizers,
  adaptively coupled slack and K, or the full A tensor A-S_B+I repair.

All four note paths are under `evidence/`. Each complete proof has an
independent review, with author/root collaboration disclosed in the
milestone manifest. The cubic note adds exactly five rational comparisons,
run once on soulkiller; the previous 28 clipping comparisons are reused.
The replay artifact is
`python3 evidence/original_mo_scalar_template_cubic_rational_certificate.py`.
Do not rerun unchanged arithmetic for another receipt. There were no
signing censuses, numerical optimizations, src-module or predicate changes.

The live unresolved implication is a sharp upper evaluation of the
weighted-shell note's (4.6) on ACTUAL coupled cells, using source and
conditional optimality. In particular its (5.3)-(5.6) retain the explicit
weighted/unweighted Delta_B and internal Delta_A discrepancies; setting
them to zero is not justified by the trace cap. The conditional Gaussian
floor does not supply this upper evaluation. The new covariance avoids
scalar attenuation and the generic deflation loss, but original all-orders
convergence and the proposed sharp dyadic inequality remain OPEN.

The new actual-sign realization package is
`evidence/original_mo_hadamard_template_milestone.json`. Its three complete
proofs passed independent reviews:

- `NOTE_2026-09-05_HADAMARD_SPARSE_FLIP_TEMPLATE.md`
  (SHA-256 `0d2355f94734b4c1e950c1e05c6df75df38b5ce181ba7fce550a4245e11328ed`)
  constructs actual cross sign matrices with flat singular bulk, finite
  outliers and asymptotically scalar-optimal SDP duals. The finite-template
  completion Gamma is an UPPER on their Boolean norm, not an attained value.
- `NOTE_2026-09-05_SCALAR_TEMPLATE_GAMMA_BOUND.md`
  (SHA-256 `bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`)
  uses matched optimal frames and exact quadratic Hermite cancellation
  to prove Gamma>283/200>sqrt(2) when q>=5/2 AND tau(C)=p q for
  the finite template. This excludes that upper certificate, not the
  actual matrices' Boolean cap.
- `NOTE_2026-09-05_HADAMARD_BOOLEAN_ALGEBRA_LOWER.md`
  (SHA-256 `68ce3f2f2a8fa2280208a9e145f508b6c2b2520d81e83185f579aaac89838a5d`)
  modifies the background to fix an entire Boolean algebra. Dense exact
  block-mean grids give a genuine actual lower Lambda_I. For PSD templates
  Lambda_I=Gamma and the actual normalized Boolean norm converges to it.
  With a symmetric POSITIVE top frame, actual liminf>1.524049912 for
  q>=5/2. Bipartite scalar SDP optimality alone is not that hypothesis;
  symmetric dilation changes the template and the actual sign family.

All three note paths above are under `evidence/`. The published rectangular
Bernstein theorem is stated and applied explicitly. One bounded soulkiller
run verified 28 exact rational comparisons for the clipping constants;
the later actual lower reuses that result without another run. No matrices
were sampled or enumerated. The live missing implication remains actual
low-norm nonsymmetric attainability and source/conditional-shell compatibility,
or a general upper using those constraints. These constructions are not
conditional optimizers and do not establish original convergence.

The preceding evaluated moment frontier is
`evidence/NOTE_2026-09-05_SCALAR_MOMENT_FEEDBACK_DIAGNOSTIC.md`
(SHA-256 `cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`).
At the old scalar endpoint, strongest source feedback plus rank-four
repaired positivity gives normalized squared upper `<9/20<1/2`.
The explicit fixed-metric repair-trace bounds require bounded scalar q,
not a conference-scale cap on A. Weak feedback still leaves a uniformly
positive gap for a formal Dirac moment law, even with the entire mixed
rounding family and its exact fourth-moment refinement. This is only a
counterexample to sufficiency of that MOMENT RELAXATION. The law does not
supply the top singular value and optimal Gram of an actual scalar dual,
an actual complete signing, or compatible source and joint-shell data.
Do not turn failure of an upper certificate into an actual width lower.

The literal old Krivine endpoint is already excluded for actual matrices
by [Braverman--Makarychev--Makarychev--Naor, Theorem 1.1](https://web.math.princeton.edu/~naor/homepage%20files/GroKri.pdf).
The explicitly defined `K_G=pi/(2 asinh(1))` in these diagnostic notes is
the elementary Krivine bound, not the exact real Grothendieck constant.
The retained-moment insufficiency theorem is still valid; it is not a
claim of insufficiency after imposing every published constraint. Actual
realizability work must move off that excluded exact endpoint and examine
admissible nearby ratios or a uniform range, not repeat its realization.

The new actual-matrix constraints supporting this evaluation are:

- `NOTE_2026-09-05_CROSS_SINGULAR_MOMENT_ROUNDING.md`: filtered cubic
  rounding with an exact quartic error and a cap-free clipped variant.
- `NOTE_2026-09-05_CROSS_SDP_COMPLEMENTARITY.md`: every optimal diagonal
  is block-balanced; its canonical-primal gap controls weighted residuals.
  Zero gap is equivalent to equal NONZERO singular values and forces a
  scalar optimum. A scalar optimum alone is not the converse.
- `NOTE_2026-09-05_ORIGINAL_PHASE_SPECTRAL_MOMENT.md` and
  `NOTE_2026-09-05_ORIGINAL_PHASE_MOMENT_BOOTSTRAP.md`: actual positive
  and negative source phases yield cubic and nuclear constraints. The
  bootstrap retains `2(n-1)+osc diag(A|A|)` with `O_C(n^(5/4))` ORIGINAL
  norm error under only `Phi(A)<=C n^(3/2)`, without replacing A.
- `NOTE_2026-09-05_CROSS_TENSOR_MIXTURE_SIGN_DEFECT.md`: actual canonical
  negative sign mass is at most half its SDP gap, or one quarter for a
  scalar optimal dual. Tensor mixing yields an explicit limiting curve
  strictly stronger than the cubic-only constraint below beta/tau=2/pi.
  Its `O_(C,t)(n^(5/4))` error requires t fixed before the n limit.

All paths in that list are under `evidence/`. Each complete proof has an
independent review; exact hashes, source aliases and reviewer roles are in
`evidence/original_mo_spectral_rounding_milestone.json`. No numerical
mathematics or signing census was used. The live next step is to exploit
actual optimal-Gram/coordinate and source-cross compatibility, or obtain
a sharper upper using it; optimizing the same weak-feedback moment-only
functional again does not address its now-proved insufficiency. The
original all-orders convergence question remains OPEN.

The underlying Gaussian upper is
`evidence/NOTE_2026-09-05_BOOLEAN_ELLIPSOID_SHELL_UPPER.md`
(SHA-256 `ede1b62a26a636179d918ba84a48d122ab013c38175bdb9cd164bcfd8bfeb9aa`).
For actual PSD C, positive P and nonnegative diagonal E<=P on a Boolean
shell z^T P z=q, it proves the exact completion-square remainder (3),
not merely a uniform sphere-to-cube multiplier. The diagonal-metric
limit is the exact cube width. Equations (12)-(17) completely evaluate
the weaker diagonal-affine specialization, including negative parameters
and singular-metric limiting infima. The stronger two-trace expression
(18) retains a trace with BOTH signs; do not substitute separate upper
bounds into it without checking the combined expression. The actual
cushioned scalar diagnostic improves the old bound but still exceeds
the desired leading constant. An indefinite reference is not a covariance.

The accompanying fixed-internal-block tool is
`evidence/NOTE_2026-09-05_CROSS_ONLY_OPERATOR_REGULARIZATION.md`
(SHA-256 `27d9ab77768e8b7afa2d48d041cf3fe6bf3b66e8b16e481ca12abcf906a28d4f`).
Its exact loss is the two selected A-cut norms plus `2n sqrt(s)`, where
s counts exceptional cross rows and columns. For `||A||<=K_A sqrt(n)`,
the normalized loss is at most `(2+sqrt(2)K_A)sqrt(Lambda C/K)` and
the new cross operator cap is `(K+8)sqrt(n)`. The source A,-A is unchanged.
Near-source selection and subsequent cross regularization have separate
slacks and cross cap `O(epsilon_A^(-4)epsilon_B^(-2))`. Do not transfer
exact optimizer properties to either regularized object or overlook
competition between this cap loss and the evaluated Gaussian gain.

Both complete proofs have independent reviews. Their provenance is in
`evidence/original_mo_boolean_ellipsoid_milestone.json`; the two elementary
scalar-check offloading exceptions are recorded there, not relabelled
as remote checks. The actual all-shell leading comparison remains open.

The current joint-shell package has three independently reviewed proofs:

- `evidence/NOTE_2026-09-05_CONDITIONAL_CROSS_JOINT_SHELL_UPPER.md`
  (SHA-256 `64d68bb2feaa59a8049d6bcc42f3ab94c845249c3088fa618916522412d0a68a`)
  proves the exchange-preserving cushioned field upper and the separate
  masked-cross conditional floor, with raw error `O(n^(16/11))`.
- `evidence/NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`
  (SHA-256 `e4919c8e16461c35efdf2963eaf9fdc1b45c07ccfba33ae1549a07e904f7ac8a`)
  uses the intrinsic cross operator `H=A tensor A-S_B+I`,
  `mu=max(2,||H||)`, and `R_mu=I+H/mu`. The entire threshold covariance
  correction is controlled; `||K||^2<=8 Phi(K)` and an elementary
  conditional norm cap give the same `O(n^(16/11))` Gaussian floor for
  actual cross optimizers over ANY original exact minimizer A.
- `evidence/NOTE_2026-09-05_INTRINSIC_CROSS_JOINT_SHELL_REPAIR.md`
  (SHA-256 `1dcd9b1e76b00887e406e505113c854b80f0661bb3bd69283f6486fb59fa2d53`)
  repairs the intrinsic linear fields by a rank-four PSD correction.
  It proves the genuine upper and retains the leading joint-mismatch
  formula with error `O(n^(5/4))`, uniformly even in vanishing-noise tails.

Use actual attainable `(x^TAx,y^TAy,x^TB_*y)` shells and conditional
optimality to sharpen/evaluate this upper. Do not drop the independent
cushion or mixed exchange term, declare the unrepaired intrinsic field
PSD, or replace conditional optimality by full-order optimality. The
needed leading comparison `F_A^*<=2sqrt(2) Phi(A)+o(n^(3/2))` is not proved.
Even proving a little-o dyadic inequality would not by itself settle
the original all-orders problem. These are optional analytic tools,
not newly mandatory proof architecture. Provenance and backup coverage:
`evidence/original_mo_conditional_joint_shell_milestone.json`.

The new whole-source variational tool is
`evidence/NOTE_2026-09-05_WHOLE_EDGE_SOURCE_PRESERVING_GAUSSIAN_REDUCTION.md`
(SHA-256 `6b22fb3ab1cc878b08fe79b5b57e0e661eaaa792dfc67f850d35db9f1b68bead`).
It uses all UNORDERED original edges and the normalized symmetric
compression `R=(L^2 I-T)/(L^2-1)`, with `T(X)=KXK` compressed to that
edge basis. For n>=3, `0<=R<=3I` for EVERY complete source; at n=2 the
compression is zero and the theorem uses an explicit independent fallback.
The entire even-Hermite correction is handled by a rank-one term and a
four-cycle operator bound. The resulting expected ORIGINAL whole-order
norm error is absolute `O(n^(16/11))`, uniform over deterministic h.

Its full symmetric Gaussian lift removes the diagonal at expected cost
at most `sqrt(n/pi)` and retains the exact augmented replica matrix
`Gamma=<sigma xx^T>`, not a positive-semidefinite substitute. Equations
(25)-(29) prove a negative-current-energy-square variance upper and an
integrated constraint for actual ORIGINAL norm minimizers, with error
`O(n^(16/11))` at `beta=n^(-5/11)`. Both independent complete reads passed.
These are SAME-order constraints. A valid mapping to the required order
upper remains unproved; do not reverse the lower inequality, assume
opposite diagonal blocks for a full optimizer, or treat a shifted
disorder-dependent posterior as another minimizer.

The threshold-optimized extension is now
`evidence/NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md`
(SHA-256 `a3ed6d9c3ee73b863c91d069e75baf9973911318a8efe9156ca61e30f55d7e25`)
and `evidence/NOTE_2026-09-05_SHIFTED_THRESHOLD_COVARIANCE_REDUCTION.md`
(SHA-256 `74457650912a515eaf6a209b184e5c1404a13fc48c68464068871ebd61236680`).
The mean-preserving OU proof is uniform in all real thresholds and keeps
every Hermite order and the actual posterior. The even covariance term
is not discarded: its PSD low-rank part has actual Gaussian Boolean-norm
cost `O(Phi(A) sqrt(log(2n)/n))`; the remaining operator error is `O(1/n)`.
Thus ANY exact original minimizer A satisfies the proved one-sided bound
`m_(2n)<=inf_h E Phi([[A,Z_h],[Z_h^T,-A]])+D n^(16/11)`,
where `Z_h=s_h A+2 phi(h)G+sqrt(1-s_h^2-4 phi(h)^2)W`,
`s_h=2 Phi_Gauss(h)-1`, and G has the universal midpoint covariance.
The threshold is fixed before drawing disorder, not chosen adaptively.
Both complete proofs have independent full-read PASS receipts.
The remaining target is an actual evaluated Gaussian upper bound, not
another covariance identity or an unsupported derivative sign.

The underlying zero-threshold rounding reduction is
`evidence/NOTE_2026-09-05_UNIVERSAL_SPECTRAL_MIDPOINT_GAUSSIAN_REDUCTION.md`
(SHA-256 `1fc6f5bbb69038b6ac4ed845d26e0724a0ceb0b5a9d96d01b4554a8e37e6f968`).
For EVERY complete source A with extreme eigenvalues `a,-b`, freely choose
`alpha=(a-b)/2`, `mu=(a^2+b^2)/2`. The exact covariance has operator
norm `(a+b)^2/(a^2+b^2)<=2`, and its arcsine-linearization remainder
is bounded by `(1-2/pi)(2/(n-1)^2+1/(n-1))`. The generic quenched
theorem therefore gives an absolute `D n^(16/11)` expected ORIGINAL
paired-norm error for all sources and all fixed internal energies.
No source regularization is needed: for ANY exact original minimizer A,
`m_(2n)<=E Phi([[A,Z],[Z^T,-A]])+D n^(16/11)`.
The full alpha domain, operator-optimal midpoint and actual-pressure
derivative are proved; operator optimality is not pressure optimality.

The immediate Gaussian upper-bound tools are
`evidence/NOTE_2026-09-05_GAUSSIAN_ENERGY_SHELL_UPPER.md`
(SHA-256 `8bd3507b722d13077cdb47e8eaa47024b8e95144900226ae4e38272795c5c728`)
and `evidence/NOTE_2026-09-05_ONE_PHASE_GAUSSIAN_VARIANCE_UPPER.md`
(SHA-256 `1646f57b060db7fdaf15c2cc8a8766806d2f00297c6749e236d8e814e467bae0`).
The first retains exact source-energy shells and gives a quantitative
one-block width deficit; its central-shell two-field comparison is still
too weak. The second removes the augmented phase at subleading expected
cost and retains the coupled posterior in the actual variance derivative.
The live target is an evaluated Gaussian order upper bound on actual
original minima. Neither a favorable integral sign nor a sufficient
cross-order inequality has been proved. An unspecified little-o dyadic
inequality alone would still not establish convergence.

The following Gibbs-generated-law results remain valid separately. They
are not prerequisites for the freely chosen universal midpoint law.

The new comparison is
`evidence/NOTE_2026-09-05_CORRELATED_SIGN_GAUSSIAN_FREE_ENERGY.md`
(SHA-256 `2e6537d0b1e2c4d8a72cc920e3fee50600d82be32417ba77c733aaedabc141c7`).
Root and both independent complete proof reads passed. For any bounded
Gaussian covariance operator, `n^2` correlated Gaussian signs have the
same quenched critical pressure as their matched Gaussian, up to
`O(n^(17/18))`. The entire posterior and singular endpoint are retained.
The key third central moment is contracted BEFORE taking absolute
values; Gaussian Holder and sign smoothing control the nonsmooth limit.
Covariance matching alone was not the proof.

For canonical sources with `||A||op<=K sqrt(n)`,
`evidence/NOTE_2026-09-05_CANONICAL_COVARIANCE_GAUSSIAN_LINEARIZATION.md`
(SHA-256 `44188dde396587f1d148e01857365b44d1bddbe83d81dc5b085ccee0cdff9854`)
has an exact disjoint-support tensor identity. Its Gaussian covariance
remainder has operator norm `O_K(1/n)`, hence pressure cost `O_c,K(1)`.
The resulting simpler covariance is `(2/pi)Sigma+(1-2/pi)I`.
The remaining live implication is an UPPER comparison of this actual
Gaussian paired pressure against the appropriate optimized smaller-order
endpoint. No such bound has been proved. The new reduction also applies
to sufficiently slowly growing caps (fixed-c error
`O(K_n^4 n^(17/18))`), hence to leading ORIGINAL norm near-minimizers
provided by same-order regularization. It does not assert that every
unregularized or quartic-penalized minimizer has bounded operator norm.

The direct original-norm consequence is
`evidence/NOTE_2026-09-05_EXPECTED_PAIRED_NORM_GAUSSIAN_EQUIVALENCE.md`
(SHA-256 `bff778718c0f357598c035edba4598f2ed67b1c49359c668958afe1c39207df3`).
It compares EXPECTED maximum absolute paired energies with normalized
error `O_K(n^(-1/22))`, with constant `O(1+K^4)`. The auxiliary choices
`c=n^(1/22)` and `epsilon=n^(-1/11)` use the full explicit bounds;
the source covariance-generating temperature is held fixed separately.
Regularizing actual original-norm minima at threshold `n^(1/99)` gives
both objective loss and Gaussian-reduction error `O(n^(-1/198))`.
The next target is therefore a Gaussian doubled-norm upper comparison
on these genuine near-minimizers, not an identification of pressure
surrogates or pointwise closeness of individual cross outcomes.

Work directly with the global optima `m_n`. A genuine advance must compare
orders or otherwise control their normalized oscillation. A construction
checked only on a selected low-norm example need not extend to actual
minimizers. A theorem for every signing under a proved norm cap does apply
to minimizers; neither distinction may be hidden in a hypothesis.

The new same-order reduction is
`evidence/NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md`.
Its SHA-256 is `8a52b7e4f171cc2089a00a6fd288e041d52605f820e49ace419ddd5fe850bec8`;
root and both independent complete proof reads passed.
For every complete signing with `Phi(A)<=C N^(3/2)`, a diagonal SDP
majorizer, vertex trimming and one jointly good random recompletion
give a complete signing at the SAME order with
`||A'||_op<=(K+8)sqrt(N)` and normalized norm increase at most
`2sqrt(Gamma C/K)`, where `Gamma=4pi/log(1+sqrt(2))`.
This applies directly to ORIGINAL norm minima. Bounded-operator
constrained minima therefore approximate the actual normalized minimum
uniformly as the bound increases; an arbitrarily slowly diverging bound
admits leading norm near-minimizers. This is one-sided objective control,
not small `Phi(A'-A)`. The missing implication is a useful order
comparison in this controlled class. Neither bounded operator norm nor
typical restriction has been proved to supply one.

The next regularized comparison is now explicit in
`evidence/NOTE_2026-09-05_QUARTIC_PENALIZED_PROFILE_IDENTITY.md`
(SHA-256 `ad393709abb35ed760986b102e1b86ab4d23c80261efec04f35d03104c821013`).
For the actual minimum of `F_c(M)+lambda tr(M^4)` on the balanced
profile, all edge flips are admissible. Every row obeys
`E_i+8lambda(M^4)_ii+8lambda sum_j M_ij^4<=c^2 d+16lambda d^2`.
Thus the fourth diagonal moments are uniformly bounded. Tensor rounding
and the Boolean norm cap give `sum|Gamma_e|=O_c(N^(3/2))`, uniformly
for `0<lambda<=1`; the diagonal SDP majorizer also bounds `tr|M|^5`.
The exact identity is
`G_N(1)-G_N(0)=c^2/4+lambda(5-9/N)-integral J_N^lambda+O_c(sqrt(N))`.
The error is uniform over `0<lambda<=1`. Each actual penalized flip
gap is nonnegative, has bounded row sums and is `O_(c,lambda)(N^(-1/2))`,
but the mixed weighted gap integral still has no proved favorable sign.
Do not substitute a permutation average for the selected envelope
derivative. Its zero-cross endpoint is twice the penalized HALF-PRODUCT
minimum and is at most twice the penalized symmetric minimum; equality
with the latter is not needed or claimed.

The pressure approximation uses one and the same recompletion in
`evidence/NOTE_2026-09-05_SPECTRAL_REGULARIZATION_PRESSURE_PROFILES.md`
(SHA-256 `2f9f63f603fcae42a952fbae53a2301eaa6b95bbe7bac2e35bcab8997d28d7d7`).
It controls both actual phases for all c in a prescribed compact interval,
with normalized cost `Gamma C c^2/(2K)+O(log(N)/N)`, while retaining
the operator and norm bounds. At fixed c, quartic penalized minima
therefore approximate original symmetric minima within
`O_c(lambda^(1/3))`. Vanishing regularized oscillation would suffice;
it has not been established. The whole-row and multi-edge variational
constraints are now proved in `NOTE_2026-09-05_QUARTIC_PROFILE_ROW_RESET.md`.
The weighted signed force kernel is controlled in
`NOTE_2026-09-05_QUARTIC_FORCE_KERNEL_BOUNDS.md`; the weighted row-tilt
identity does not assert a fourth moment at its endpoint. Independent
coefficient refills in `NOTE_2026-09-05_QUENCHED_BIASED_COEFFICIENT_REFILL.md`
retain the full quenched posterior and exact quartic correction, with
`O(sqrt(N))` replacement error even over all edges. The separate actual
canonical cross law adds at most `(41+88C^2)lambda t n` to the paired
quartic penalty. None of these same-order finite variations supplies
the missing Gaussian endpoint inequality; no new signing census is needed.

The actual-Gibbs structural proofs are
`evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md`
(SHA-256 `46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`),
`evidence/NOTE_2026-09-05_EXACT_HALFPRODUCT_SUBCRITICAL_SPECTRAL.md`
(SHA-256 `10dfe02b63aa3c4aa987ce48d4a3e660e90509b43e6a50a1a002ba9ecc1cc522`),
and its strengthening
`evidence/NOTE_2026-09-05_HALFPRODUCT_NEARMINIMIZER_STRUCTURE.md`
(SHA-256 `dccc256d3b7119c666102e54cffe3a2026d31edc1bcd0c4366a15ce92c762f0f`).
A Boolean energy cap gives a positive extensive response to any field
with a positive density of moderate nonzero coordinates, even with
unbounded outside coordinates; the ACTUAL Gaussian posterior is retained.
For EVERY leading half-product near-minimizer at fixed `c/sqrt(N)`,
approximate optimality, eigenvector truncation and sparse pinning prove
`||A||_op=o(N^(3/4))`. Deleting ANY `o(N)` vertices also changes the
full Boolean energy uniformly by `o(N^(3/2))`. Every gap is retained.
These are not exact-minimizer-only properties. Half-product pressure
approaches half the energy WIDTH, not necessarily the absolute norm;
no original-norm transfer or fixed-fraction comparison is inferred.
Complete root and independent reads passed; see
`evidence/original_mo_spectral_structural_root_review.md`.

The singular full-strength criterion is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_SPECTRAL_DEFICIT.md`:
`V_r=tr[-rI-H/mu]_+=o(N)` implies the stated quenched mean failure
and vanishing success probability, not an exponential original-law tail.
The construction-cap example is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_CONSTRUCTION_CAP_HOSTS.md`.
The new subcritical base strengthens this to actual leading HALF-PRODUCT
near-minimizers in
`evidence/NOTE_2026-09-05_FULL_STRENGTH_HALFPRODUCT_NEARMINIMIZERS.md`
(SHA-256 `ad83095163cf8e969e542a6626382dabaa5adb5e2ffce4bfffea274c813b53e4`).
One reused sparse module costs `o(N)` pressure and gives exactly
`V_r=2(1-r)` eventually. This is not an exact-minimum, original-norm
near-minimum or unrestricted selected-outcome exclusion. The separate
`evidence/NOTE_2026-09-05_NUCLEAR_SPECTRAL_BUDGET.md` gives
`Phi(A)>=N^2(N-1)/(pi tr|A|)` and a linear nuclear effective rank under
the relevant objective caps, not spectral flatness. No new mathematical
computation was used in these results.

The actual covariance corollary
`evidence/NOTE_2026-09-05_ACTUAL_GIBBS_COVARIANCE_FLOOR.md` gives a
positive diagonal component of trace at least `chi_c N` in each actual
phase covariance. Thus rank-`o(N)` truncations cannot have `o(N)` tails,
and every integral full cross block has `qbar>=chi_c^2 N^2` under a
fixed Boolean cap. This retires the conditional sublinear-rank/tail
escape at fixed c, not the radial upper comparison, which can be
quadratic as well. The constant is not uniform at zero temperature.

The new unconditional selected-restriction theorem is
`evidence/NOTE_2026-09-05_SELECTED_HALF_RESTRICTION.md`
(SHA-256 `c8a9aa0b8c44fb14f444955fbe3eec8cba8e7f19c01fb8eeb2596418d3416c02`).
Its complete independent root read passed. An explicit odd cycle of
disjoint subsets gives a half-norm restriction at order `2n+1`; a
complementary-phase exchange argument gives boundary error `(n-1)/2`
at order `2n`. Applied to actual global norm minima, these prove
`m_(2n+1)>=2m_n` and `m_(2n)>=2m_n-(n-1)`.
This improves the old fixed-partition/monotonicity estimate, but the
normalized comparison still has factor `sqrt(2)`, not `1+o(1)`.
The missing leading factor is the issue, not the linear boundary error.
Do not confuse selected half-norm restrictions with typical tiny-n
restrictions at the sharper source-normalized threshold.

The exact full-strength identity is
`evidence/NOTE_2026-09-05_FULL_STRENGTH_BOUNDARY_LIKELIHOOD.md`
(SHA-256 `8703433f6118f00dd589d711e9541f558489caa3d13059f8e71405333401fdb2`).
Root and independent complete reads passed. The derivative of the actual
planted log likelihood is a weighted sum of mixed finite differences
under PAIR-DEPENDENT Gaussian boundary laws, and its integral is valid
through singular `rho=1`. An exact actual order-three minimizer refutes
the coordinatewise sign premise; its counterexample context disappears
from the endpoint support. This does not refute the weighted average or
the full-strength finite-step comparison. Retain the boundary support,
full posterior and possible cancellations; the negative prior trace
alone does not control this integral. No computation was used.
The complete root review for both notes is
`evidence/original_mo_boundary_restriction_root_review.md`.

The fixed-strength strengthening is
`evidence/NOTE_2026-09-05_FIXED_WEAK_GAUSSIAN_CUSHION.md`
(SHA-256 `5df7258c4cf99deac09eaeb4a206e1270ffa7add1e49e176b70a4a232eb54d12`).
Root and independent complete reads passed. For ANY latent Gaussian
correlation matrix `S`, conditional independent-sign replacement and
convexity give the actual quenched floor
`[c sqrt(2t) K0-2log2-c^2 t arcsin(rho)/(2pi)]n-o(n)` for covariance
`(1-rho)I+rho S`. Thus a NONEMPTY interval of fixed positive strengths
is excluded at suitable fixed `c,t`, even though its information is
not `o(n)`. Precisely, the gap
`Delta_rho=c(sqrt(2t)K0-1)-2log2-c^2 t arcsin(rho)/(2pi)` must be positive.
The actual centered latent law has `||S||_op<=4n-3` for EVERY generating
host, so the heat-martingale and conditional bounded-difference bounds
make success exponentially rare. Even `exp(o(n))` proposals with these
marginals, including legitimate pre-draw host mixtures, fail; proposals
need not be independent and the internal host may be selected afterward.
This is not a mixture entropy lower bound or an exclusion of `rho=1`.
The complete independent review is
`evidence/original_mo_fixed_weak_gaussian_cushion_exact_review.md`.
No new computation is used. Do not repeat weak fixed-strength sampling
at these parameters or extend this result outside its explicit gap.

The preceding dependent-rounding information theorem is
`evidence/NOTE_2026-09-05_GAUSSIAN_SIGN_INFORMATION_SCALE.md`
(SHA-256 `5846e981204f03230bbfd415443824d1a320840d56b6163267e37ee1b8e5e566`).
Complete proof reads passed. Every sign law satisfies
`D(Q||iid)>=||C-I||_F^2/(4||C||_op)`, with the SECOND-MOMENT matrix
`C=E bb^T`, not a centered covariance absent a mean-zero hypothesis.
For any Gaussian correlation matrix, including singular ones, arcsine
and the Schur product give
`D(sign N(0,Sigma)||iid)>=||Sigma-I||_F^2/(pi^2||Sigma||_op)`.
For the actual centered tensor `H`, `mu=-lambda_min(H)` and
`Sigma_rho=I+rho H/mu`, this implies `Omega(n)` discrete information
at every fixed `rho>0` on norm-capped hosts, including singular `rho=1`.
The proof uses `||A||_op^2<=16Phi(A)` and retains `mu` in the operator
denominator before combining the ratio. No Gaussian determinant upper
bound is used at the singular endpoint. Thus full-strength canonical
rounding is OUTSIDE the low-information exclusion, not proved successful.
On actual half-product minimizers, `rho=o(n^(-1/2))` instead gives
`o(n)` information and is excluded in mean and with substantial success
probability by the following quenched theorem. Strengths outside the new
cushion criterion, the actual Gram--Schmidt law and unrestricted selected
cross blocks remain open. A conditional-law entropy lower
bound must not be extended to arbitrary mixtures over hosts.
No computation is used; the independent general entropy audit is
`evidence/original_mo_entropy_covariance_review.md`.

The preceding iid all-orders cross-block result is
`evidence/NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md`
(SHA-256 `97e1aeb3ac25c2570072d9f0ebdb0c4387f739ed3c005ec7b43d30409dd7ade4`).
Root and independent full reads passed. An explicit Gaussian martingale
control in the sourced zero-temperature Parisi formula proves
`P_SK>=K0=4/(3sqrt(pi))>1/sqrt(2)`. A host-free pure-cross pressure lower
bound, Gaussian covariance interpolation, and direct fixed-temperature
Bernoulli replacement give `E F>=(c sqrt(2t) K0-2log2)n-o(n)`.
Against `2R_n<=cn+o(n)`, the gap is positive when
`Delta=c(sqrt(2t)K0-1)-2log2>0`. Bounded differences prove an exponentially
small iid probability of ANY good internal host at such a cross block,
so `exp(o(n))` proposals with iid matrix marginals cannot succeed even
when dependent across proposals. A successful law must have relative
entropy at least `(Delta^2/(c^2 t)+o(1))n` from iid signs.
Thus iid QUENCHED and `o(n)`-information selection are excluded at those
fixed parameters; arbitrary dependent selection is not. The exact
planted-channel identity retains the reverse relative entropy and the
full actual Gibbs prior. No numerical SK constant or new computation is
used. Do not rerun iid samples or confuse this theorem with the earlier
annealed Gaussian-sign obstruction.

The independently reviewed all-orders coefficient results are
`evidence/NOTE_2026-09-05_POSITIVE_CONE_TRUNCATION.md`
(SHA-256 `632adeb92932db37ba1ac218621eb3f7d1b8bd24e8461273abf74a379d79d304`),
`evidence/NOTE_2026-09-05_EXTENSIVE_COEFFICIENT_MOMENTS.md`
(SHA-256 `b07772332265dea635c59a7d293562feedb5c57cb7b66d7850f77c1ffbd4107e`),
and `evidence/NOTE_2026-09-05_POSITIVE_DEGREE_SELECTOR.md`
(SHA-256 `20dae4c37ece2f5c5808595c54941de1b10a241d03b63c4431c76dc373849875`).
The exact central-factorial coefficients are nonnegative. A cutoff at
`k<=K_N=o(N)` loses extensive pressure for fixed `c>pi log2`, even on
actual norm or symmetric-pressure minimizers; the latter require their
separately justified norm cap. In any fixed positive `k/N` band, the
coefficients are within dimension-uniform multiplicative constants of
`E|Q_A|^(2k)/(2^k(2k)!)` for EVERY complete signing. Convergence of the
separately optimized coefficient rates at unbounded fixed `k/N` values
would imply convergence of `alpha_N`, but that transport is still open.
One selected extensive degree per signing already captures log pressure
to `O(log N)`; do not misstate the cutoff theorem as excluding sparse
degree selection. Mixed minimax is legitimate with its quantifiers;
moving a pure minimum through the coefficient sum is not. No new
census or simulation is used by these proofs. The complete independent
review and correction record is `evidence/original_mo_coefficient_quenched_review.md`.

The new fixed-order analytic theorem is
`evidence/NOTE_2026-09-05_EXACT_OPTIMIZED_ORDER_SIX_PROFILE.md`
(SHA-256 `a1469b34118da1bf971c7d53ad0fb8c50525f588a42bf7c4f2dda9b132966fd4`).
Root and independent complete proof reads passed. For ALL `u>=v>=0`,
the minimum of `E cosh(u I+v C)` over complete order-six signings is
`cosh(v)*(3X^2+3Y^2+2Y-4)/4`, where `X=cosh(2u),Y=cosh(2v)`.
For `v>0`, minimizers are exactly `A^2=5I`; there are twelve after
first-row-positive switching normalization. The proof compares every
coefficient on `X=1+p+q,Y=1+q`, not a finite grid. Its success uses the
candidate's exactly quadratic polynomial: low-moment comparisons do NOT
control higher positive candidate coefficients at larger orders.
Along `u=c sqrt((2-t)/6),v=c sqrt(t/6)`, the optimized endpoints cross
exactly once for positive `c`, and
`f6(c,2/17)-f6(c,0)=((sqrt(17)-4)/sqrt(3))*c-log(2)+o(1)` as `c->infinity`.
Thus no temperature-uniform bounded interior excursion holds even at an
actual global optimum. This does not refute a fixed-`c` small-oh order
comparison or convergence. The finite-temperature maximum is not asserted
to occur at `2/17`; the exact left derivative at `t=1` is positive.
Do not rerun the same catalog/grid or treat fourth moments alone as an
all-orders extension. A new argument must control the actual higher
coefficients, selected finite-step pressure, or another order comparison.

The latest finite-step theorem is
`evidence/NOTE_2026-09-05_FINITE_STEP_ROUNDING_ANNEALING.md`
(SHA-256 `058cdd3e17972be45a664b21e720fafd44c744194c2e8f3bcb37e818d474ee0a`).
The full proof passed root and independent reviews. Its Gram--Schmidt
bound retains `log E_nu exp((gamma^2/2) v^T G^-1 v)`; replacing this
log moment generating function by the average quadratic form is invalid.
For ALL `G>0,diag G<=1`, the resulting proxy has floor `c^2 t n/4`.
A separate Gaussian entropy-tilt proof gives the ACTUAL canonical
Gaussian-sign annealed floor `c^2 t n/(2pi)-o(n)`, uniformly over hosts
and admissible centering, with no covariance operator-norm hypothesis.
Here `c,t` are fixed as `n` grows. Since `2R_n<=cn+o(n)`, the respective
annealed certificates cannot give the needed small-oh finite-step
comparison when `c>4/t` or `c>2pi/t`. The Gaussian claim concerns
`log E_B exp F_B`, not `E_B F_B` or `min_B F_B`. The actual Gram--Schmidt
law is not excluded. The preceding local second-moment theorem survives.
Do not optimize the same quadratic proxy or substitute annealing for
selected-outcome control again; shrinking steps and different laws remain
outside the stated obstruction. No finite sample proves this theorem.

The preceding integral construction is
`evidence/NOTE_2026-09-05_INTEGRAL_CROSS_BLOCK_COVARIANCE_ROUNDING.md`
(SHA-256 `c02bcc4d73ca58ba701b80a1fd73fa1c54f928effd5a62fe77daa9925c7d5c01`).
Root and two independent agents checked the complete proof. With actual
opposite-temperature covariances `U,V`, define
`qbar(B)=(tr(B^T U B V)+tr(B^T V B U))/2`, allowing all `n^2` entries
of `B` to be independent choices of signs, including its diagonal.
The complete-sign Gaussian construction proves
`min_B qbar(B)<=n^2-8(a_A')^2/(pi ||A||_op^2)`.
Its sharper form uses the exact negative spectral edge of the centered
tensor matrix; scalar entrywise arcsine is justified by disjoint entry
types, not matrix functional calculus. The host is fixed during rounding.
General Gram--Schmidt covariance rounding additionally gives integral
spectral-tail bounds and retains fixed coordinate squares through diagonal
shifts. The conference-form scalar-shift optimum is not a limitation
theorem for all diagonal shifts or all rounding methods.
The general comparison `min_B qbar(B)<=2a_A'/beta+o(n^2)` remains open;
the sufficient low-effective-rank case is not established for minimizers.
Even that endpoint bound would not control the integrated balanced path.
No new computation is used by these analytic theorems.

The preceding analytic method check is
`evidence/NOTE_2026-09-05_NEAR_MINIMIZER_OPPOSITE_PHASE_COUNTERFAMILY.md`
(SHA-256 `8130ca8c0af67d9976f71f086a79607d0b7b640b1e5c35ba6eb08d87e81324f7`).
Its complete proof passed root and independent review. For fixed `c>0`
and every sufficiently large `N`, paired modules can be planted into an
arbitrary old signing while changing every spin energy by at most
`O_c(N^(11/8))`. Conditional entropy costs only the new vertices;
Rademacher averaging retains the full Gibbs feedback. A simultaneous
thermal/operator event and exact even-module covariance decomposition give
`tr(A U A V)=Omega_c(N^(9/4))`. Choosing an old norm minimizer gives norm
excess `O_c(N^(11/8))`; separately, choosing a half-product minimizer at
the SAME raw `beta=c/sqrt(N)` gives pressure excess `O_c(N^(7/8))`.
Both are leading-order near-minimizers, not merely correct-scale hosts.
They need not be the same family and are not edge-local half-product
minima. Do not extend the counterexample to exact minimizers, or treat
this moment condition as necessary for convergence. No computation is
used in the all-orders proof.

The preceding analytic method check is
`evidence/NOTE_2026-09-05_FULL_ROW_CAVITY_COUNTEREXAMPLE.md`
(SHA-256 `54de76afacf34c7443ece9f5a34c42ef32d741b6fa381a7f5b9412675a1b331f`).
Root and independent full reads passed; no computation was used. A positive,
exchangeable, even arbitrary cavity has a strict minimizing sign row up to
global reversal, yet its actual row-tilted second and fourth moments grow
at least as `sqrt(d)` and `d^(3/2)` at fixed critical row scale. Every subset
replacement and the complete row-noise hierarchy holds. Thus these
inequalities alone cannot prove the desired bounded tilted moments.
No actual quadratic-host realization was supplied. Do not relabel this
as an Ising counterexample or make bounded moments mandatory for convergence.

An independently derived local endpoint calculation retains actual phases.
Let `A` minimize `a_A(beta)=(log Z_+(A)+log Z_-(A))/2`, where
`beta=c/sqrt(n)`, and let `U,V` be its two phase covariance matrices.
The earlier coherent admissible choice in the paired family `A,-A`,
with cross block `B=A+D` and a fair independent signing of its diagonal, gives
`f_(2n)'(0+)<=-beta a_A'(beta)/2+(beta^2/4)(tr(A U A V)+n)`.
The exact derivative minimizes over ALL active block pairs and cross
signings; pairing with the negative is only an admissible upper bound.
The integral construction above allows other cross blocks; this coherent
trace condition is not required for that enlarged choice.
The trace comparison that would bound this derivative above by `o(n)` is unproved.
Even such an endpoint derivative bound would not by itself control the
whole interpolation. The formula is analytic; the fixed-order numerical
check below neither proves nor refutes an asymptotic small-oh comparison.

The preceding two analytic results are
`evidence/NOTE_2026-09-05_ADAPTIVE_PERTURBATION_CORRELATIONS.md`
(SHA-256 `054063ac00e2fda45b676fc9a257cb901f43627e83cf31ce0e061b7c8816bb5f`)
and `evidence/NOTE_2026-09-05_OPTIMIZED_GAUSSIAN_SWITCH_MEASURE.md`
(SHA-256 `9c5090ddf4e1b43222716182ce5de3c51216ad87cb6f216d4bcd3ea70571fa0a`).
Both complete proofs passed root and independent review; they have no
new finite-check coverage. For arbitrary fixed additive edge noise `E`,
edge-local sign optimality gives `sum |Gamma_e|<=4k tanh(beta)+||E||_2^2`,
with rowwise and balanced-profile versions. The physical edge flip is
`-2A_e`, not `-2(A_e+E_e)`. The sharper alternative is
`sum |Gamma_e|<=2k tanh(beta)+Phi(E)`, also groupwise. Thus Gaussian noise
`E=epsilon G` preserves signed Frobenius diffuseness whenever `beta->0`
and `epsilon=o(sqrt(N))`, including bounded noise and arbitrary edge-local
adaptive selections. At fixed critical `c` and bounded `epsilon`, the
balanced-path squared-correlation error is `O_(c,epsilon)(sqrt(N))`.
No corresponding bounded-noise row-operator improvement is claimed.

For `Psi(G)=min_A log E_(sigma,x) exp(sum u_e(A_e+epsilon G_e) sigma x_i x_j)`,
the weak Hessian is the selected smooth Gibbs covariance minus a PSD
switch measure `M`. Its standard-Gaussian-weighted trace is at most
`epsilon^2 ||u||_2^2+epsilon sqrt(2(N+1)log 2)||u||_2`, hence `o(N)`
for critical profiles and vanishing `epsilon`. The proof keeps optimizer
adaptation throughout; an exact order-two cusp shows why switches cannot
simply be dropped. Mixed variance paths need a bounded relative velocity
for this trace estimate. Changing deterministic weights still produces
`sum u'_e A_e Gamma_e`; a heat identity does not remove that term.
The physical-noise-coordinate switch bound is divided by `epsilon^2`.
Do not claim that bounded covariance transport supplies different
deterministic block endpoints or a cross-order comparison.

The coordinatewise switch identity also gives, for GLOBAL Gaussian
envelope minima, `E L_g<=2K_g tanh(u_g)+epsilon^2 u_g K_g`. This is
distinct from the Boolean-energy bound for arbitrary edge-local choices.
The exact changing-profile derivative is
`sum u'_e E[A_e Gamma_e]+epsilon^2 sum u_e u'_e E(1-Gamma_e^2)`
minus `sum (u'_e/u_e) m_e`, where `m_e` is the Gaussian-weighted
diagonal switch mass. The first term remains uncontrolled, even when
monotone variances make switching favorable. In the noiseless case the
signed normalized edge-flip-gap sum differs from the previous `D_N(t)`
by at most `c^4`. A martingale sign-flip generator therefore reproduces
that defect; it does not supply its sign or another independent obstacle.

The preceding analytic result is
`evidence/NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md`
(SHA-256 `96d2675bd0cf1ee48e962b2974a2b8649afc487454a3912044ee1e737c53a9a5`).
Its full integration passed two independent reads. It proves uniform
sparse near-flat rounding, a common diffuse near-maximizing ensemble for
global norm minimizers, and actual signed Gibbs diffuseness for pressure
minima. Along the specified balanced two-block path, the exact formula is
`f_N(1)-f_N(0)=c^2/4-integral D_N+E_N`, with
`|E_N|<=c^3 sqrt(N)/2+c^4/6`. The shared-sign endpoint is exactly twice
the minimum half log-product of the one-sided partition functions, hence
`f_N(0)<=2P_(N/2)(c)`. A lower bound `integral D_N>=-o(N)` would give
dyadic pressure subadditivity; it has not been proved. Do not require
endpoint equality for that direction, and do not treat a bare dyadic
small-oh estimate as all-orders convergence. No pressure minimizer is
silently identified with a norm minimizer.

The added local corollaries give a rowwise signed operator bound and a
bounded cavity exponential normalizer. Slow cooling of pressure minima
produces asymptotically norm-optimal sources with actual near-maximizing
Gibbs ensembles. The missing unsigned star-fluctuation estimate does not
follow from these signed bounds or from the cavity normalizer alone.
These corollaries are analytically reviewed, not new finite-check coverage.

The previous fresh analytic result is
`evidence/NOTE_2026-09-05_INDUCED_OPTIMIZER_RESTRICTIONS.md`
(SHA-256 `ab65d46bb48627170344219850131aa77ed9cbe9d152e143346a7fec71d42409`).
Root and independent review checked the full proof. For `n -> infinity`,
`n^2=o(log N)`, the complete induced signing law has explicit total-variation
control. Every smaller signing occurs, but uniform restrictions have
typical normalized norm at least `(2/3)*sqrt(2/pi)>1/2`.
The failure result extends to `exp(o(n))` samples with uniform marginals,
even when dependent. Do not omit the growth/scale or marginal hypotheses.
Do not confuse existence of an `m_n`-optimal restriction with existence of
one matching the source constant; the latter comparison remains open.

Do not automatically resume residual (ii), the old equation (33), a skew
ansatz, or a finite-prime census. They are optional archived avenues, not
the definition of progress. Before revisiting one, name the changed
premise and the implication for the original question.

## Verification

The order-six optimized-profile runs are recorded in
`evidence/original_mo_optimized_profile_mesh.json`. A single NUKA exact
integer/Fraction run covered all 1,024 switching-normalized signings and
all 64 spin states, produced 23 joint signatures and 22 nonnegative
polynomial difference certificates, and passed 3,397 checks in `0.184`
seconds. An independent V100 run produced all joint histograms and the
prescribed 455 floating-point profiles in `2.938` seconds. A separate
NUKA comparison of the already stored outputs, without re-enumeration or
pressure replay, matched all 20,480 histogram entries exactly and passed
3,285 checks. Floating-point near-minimizer tolerances do not classify
exact ties or derivatives. All three runs exited 0 once; worker absence
was verified. Full result JSONs, the GPU array archive, reviewed sources,
exact commands and raw/preserved hashes are retained. The analytic theorem
does not depend on any of these computations and concerns only order six.

The new finite-step mesh run and independent replay are recorded in
`evidence/original_mo_finite_cross_mesh.json`. Soulkiller's V100 evaluated
8,192 canonical Gaussian-sign cross blocks, 8,192 independent blocks,
and the two coherent references on ONE fixed order-six conference host,
using all 4,096 spin pairs at 20 prescribed `(c,t)` profiles. This is
neither a new host-minimizer census nor an exhaustive cross-block search.
NUKA independently enumerated all 16 order-two cross blocks and the
four-point Gaussian support, then replayed the selected order-six GPU
pressures with full direct spin sums. Sample minima are only upper
bounds; sample log means need not approximate rare-event annealed values.
The V100 run completed once in `4.249` seconds, exit 0. NUKA passed
1,105 order-two formula checks and 480 checks on 160 selected order-six
cases; maximum CPU/GPU pressure difference was `7.11e-15`. This replay
checks the pressure and endpoint, not the GPU's `qbar` values. Both
workers exited normally and absence was verified. At all 20 order-six
profiles the best Gaussian sample is index 1067; exact signed-permutation
algebra identifies it with the known `A-I` construction for every
temperature and step; the complete identity is in
`evidence/NOTE_2026-09-05_SAMPLED_CROSS_BLOCK_ORBIT.md`
(SHA-256 `98923ba2cf14f71b71511b7896734028a48d1c029866fcc88592c40d820da1aa`).
Thus this sample found no new noncoherent winner,
not an exhaustive proof of optimality. Do not enlarge the same sample
without a changed mathematical premise.
The finite results are not all-orders evidence, and no larger sample
or unchanged successful run is required for a cleaner receipt.

The new fixed-order opposite-phase probe ran on soulkiller's V100
(`2.434` seconds, exit 0) and independently on NUKA CPU (`0.603` seconds,
22,528 checks, exit 0). It examined 1,024 switching-normalized order-six
signings at exactly `c=0.5,1,2,4,8`; no larger census was run. CPU used
64 spin states and covariance traces, while CUDA used 32 antipodal representatives
and direct squared bilinear moments. Their candidate values agree within
`2.85e-14` in `T`; all five profiles have the same 12 numerical minimizing
signings. Floating-point comparisons are not a rigorous optimizer
classification. The positive finite virial gaps at the tested `c=1,2,4,8` do not refute
an asymptotic small-oh allowance, and this check is not evidence for the
all-orders planted theorem. Both workers exited normally and absence was
verified afterward. Exact commands, source hashes, results, tolerances,
and cleanup receipts are in `evidence/original_mo_opposite_phase_n6_mesh.json`.

One new soulkiller run of
`scripts/original_mo_weighted_pressure_n4_check.py` passed 7,110 formula
checks, exit 0, with one CPU worker. Its scope was exactly order four,
64 signings, 16 spin states, and six prescribed weight/temperature
profiles. This is a finite regression of the new pressure identities,
not a larger-order census or theorem certificate. The exact command,
reviewed proof input, input hashes, full log, and live preflight are in
`evidence/original_mo_weighted_pressure_regression.json`. There was no
rerun of unchanged mathematics. The all-orders claims rest on the proofs.

The reviewed technical reset replay on soulkiller passed 40 tests in
44.83 seconds: the new global registry, both independence directions,
legacy route aliases, and three existing wrapper regressions. It ran in
`/tmp/original-mo-reset-replay.W2zk3m` with explicit files and one worker.
The separate documentation replay checks the final entry documents and
retained proof scopes. Its result, both exact commands, input manifests,
and log hashes are recorded in
`evidence/original_mo_route_reset_regression.json`.

The earlier diagonal work's receipt is
`evidence/original_mo_diagonal_regression.json`: 65 technical tests and
17 documentation tests passed across two runs after three missing staging
inputs were supplied. It is not a verification of the present reset.

No convergence claim may be accepted merely by toggling a Boolean or by
closing optional-route checkboxes. A complete reviewed proof is required.

## 2026-09-12 late: duplication audit (user-flagged) and the live route

The user flagged that the day's completion-discrepancy work duplicated
existing material. Audit confirms it: `NOTE_2026-09-01_RG2_EQUAL_ENDPOINT_
PALEY_SHIELD.md` (and its successors on the multiplier-two ray) already
contains the cross-term floor, the independent-budget obstruction, the exact
residual reduction (its eq. (14)), an equal-endpoint doubling frame, a
balanced near-conference skew construction with a PROVED Dini-summable error,
and three geometric shields. The 2026-09-12 framework note now carries a
prior-art notice marking its weaker status. The live route per the program's
own documents is the RG2 equal-endpoint diamond on the residual (14) at
multiplier two (multiplier three remains necessary beyond), with the
Paley principal embedding and degree-balancing reversals as the degrees of
freedom to choose. Read `ARTIFACTS.md`, the RG2 note, and the multiplier-ray
notes BEFORE any new computation on this problem.

## 2026-09-12 night: NS-technique port executed at n=26 — exact negatives (do not repeat)

User-directed: the imported Navier--Stokes techniques were executed at full strength
(stress set = complete tied-maximizer band; exact scoring; correction cycle) on the
verified order-26 record signing (norm 61). Exact outcomes: every single flip -> 63
(all 325); every pair -> 65 (all 52,650); best triple -> 67 (all 5,668,650 evaluated);
CP-SAT exact: no k-flip repair reaches <= 59 for k = 1..8 (k=1,2 cross-validated by
brute force; k=9 UNKNOWN at cap). The record is repair-rigid to depth 8. The transfer
is bookkeeping (stress set / correction cycle), not an engine, and has no path to the
asymptotic question. Do not re-run these techniques as a route to E(1) or L.
Artifacts: `evidence/ns_port_n26_undercut_20260912/` (README, receipts.txt, results.json,
scripts, band data).

Also landed this session: `evidence/NOTE_2026-09-12_paley_fourier_identity.md`
(exact Fourier identity for the Paley cube-max, one-class rho=1 criterion, explicit
Cx = p x cylinder family as a cross-check of the k=1 stratum; no prop number, no
status change).

## 2026-09-12 late: outer-layer pair-cover defect

`NOTE_2026-09-12_OUTER_LAYER_PAIR_COVER_DEFECT.md` adds all-orders
identities for the actual one-edge witness layers.  An opposite-sign outer
pair leaves at least `b/2+M/2-(eps_plus+eps_minus)/4` edges uncovered on its
Hamming cut, where `b=d_H(x,y)(n-d_H(x,y))`; in particular outer-two
witnesses have cut size at least `M-2`.  A same-sign pair leaves the analogous
defect `a/2+M/2-(eps_x+eps_z)/4` on its agreement side, where
`a=binom(n,2)-b`.  Thus every pair type has a certified uncovered set, and
no one positive/one negative pair covers every edge when `M>2`.  The exact
sign identities were exhaustively regressed through order four on the
controller and independently replayed on NUKA.

This is a structural strengthening of the edge-flip witness statement, not
the required outer-layer multiplicity theorem: a larger diffuse family may
still cover all edges.  It does not supply the Banaszczyk cover condition,
the RG2 residual diamond, multiplier three, or convergence.

## 2026-09-12 late: homogeneous three-witness extension

The same outer-layer note now gives a genuine multi-state obstruction.  If
three positive states cover every edge by their one-edge witness sets, then
their errors obey
`eps_x+eps_z+eps_w >= 3M-3n/2`; the identical assertion holds for three
negative states.  Hence three same-sign outer-two witnesses cannot cover
every edge once `M>n/2+2`.  The proof is an exact three-energy expansion on
the four signature cells of the two relative sign vectors, followed by the
convex lower bound on the common-agreement edges.  Controller exhaustive
tests through order four and an independent NUKA replay passed.

At this point the homogeneous calculation alone did not address mixed-sign
triples or the RG2 witness-family / rounding gap.

## 2026-09-12 late: mixed triple closure

The same note now closes the apparent three-state mixed-sign exception.
For one positive and two negative outer states, the edge region on which both
relative products are negative is forced to have the positive-gauged sign
`-1` by a putative cover.  The exact four-state Walsh identity with
`w=x circ y circ z` then gives
`eps_x+eps_y+eps_z >= 2M+4N_{--}`.  The sign-reversed statement covers two
positive and one negative state.  Hence no mixed outer-two triple covers all
edges if `M>3`.  Combined with the homogeneous result, no triple of
outer-two witnesses can cover all edges once `M>n/2+2`.

Controller exhaustive algebraic regression through order four and an
independent NUKA replay passed.  This is a cardinality-four lower bound for
any such cover, not the required growing-family theorem; arbitrary larger
mixed families remain the live gap.

The identities retain quantitative mass: a homogeneous outer-two triple
leaves at least `M/2-n/4-1` jointly uncovered edges, and a mixed outer-two
triple at least `M/4-3/4` (with an additional positive term from its
double-negative signature region).  This is the first route toward a
growing-cover lower bound, but no capacity theorem prevents a fourth witness
from covering those sets.

For an exactly four-state cover, the corresponding eight-state Walsh identity
forces its prescribed relative-signature class to contain at most one edge at
the outer-two scale.  This is a signature-collapse condition, not a
four-witness exclusion: a vertex-cell support can avoid one matching of the
signature cube.  Controller exhaustive order-three algebra and independent
NUKA replay passed.  The formerly proposed transport/capacity continuation
was retired on 2026-09-13: it is not a convergence criterion and must not be
treated as a live research queue item.

## 2026-09-13: user-stopped repetition; no Hadamard restart

Nick flagged repeated Hadamard work and directed that the route stay closed.
The detour rebuilt the Sylvester `H-diag(H)` family already present in the
September 6 campaign archive. The archive member is
`original-mo-campaign-root.Z0Um0rAh/original_mo_symmetric_hadamard_trades.py`
inside `evidence/original_mo_broad_campaign_20260906/independent_and_hadamard.tar.gz`,
SHA-256 `f8f038f7e047f58959901eaac7f567ec2909633ab76ddd18734126d901aaaed0`.
The existing amplification losses are in `solution.md`, Proposition 7.2
and Section 13. The one-edge witness limitation is in
`NOTE_2026-09-02_BEST_RESPONSE_EDGE_FLIP_NO_GO.md`, Section 3.

The newly created construction module, test module, standalone proof note,
independent verifier, and result receipt were removed from the active tree.
The earlier committed pair/triple identities are retained, but their proposed
growing-cover continuation is not a convergence criterion. Do not reinterpret
the preceding chronological entries as authorization for another Hadamard
family scan or small-cover detour. No new global milestone was achieved.

All five removed files remain recoverable in the verified pre-removal backup:
`/mnt/storage/backups/codex/quadratic-minmax-limit-2212b99a4c24/20260913T003015Z-a5337eec5fd5-four-cover-route-closure.V0ZWS9`.
Covered HEAD: `a5337eec5fd58a8cb8dfae4a861c6e45c2f80fe2`, plus the full
then-uncommitted tree. SHA-256 of `SHA256SUMS`:
`17a307b3399dbb768255c2e3a9774dfb925c6fa619fc6dbf14577705971e70a6`.
No second backup is needed for this correction. The actual unresolved task
remains convergence of `m_n/n^(3/2)`; a next attack needs a specific
unproved implication checked against linked trees and archived results first.

Infrastructure correction from Nick: Horus has replaced Lucky as the
DNS-only node; Horus is excluded from compute dispatch. ARTIFACTS and the
mesh manual now reflect this.

## 2026-09-13: bounded counting entropy is a direct all-orders reduction

For `E_n=binom(n,2)` and `C_n(t)` the number of labeled signings with
`Phi(A)/n^(3/2)<=t`, the new note
`evidence/NOTE_2026-09-13_BERNOULLI_COUNTING_ENTROPY.md` proves that for every
fixed `0<delta<=1`, all sufficiently large `n` have at least
`exp(n^2 h(delta^2/64)/8)` signings at norm at most
`m_n+delta n^(3/2)`. This uses independent Bernoulli edge flips around an
arbitrary signing, the existing Bernstein bound, and a Shannon entropy
argument. It is not a finite signing census and makes no family assumption.

Thus `s_n(t)=log(1+C_n(t))/E_n` satisfies
`liminf alpha_n=inf{t:limsup s_n(t)>0}` and
`limsup alpha_n=inf{t:liminf s_n(t)>0}`. Convergence of `s_n(t)` on a dense
set of thresholds is therefore sufficient for the original convergence
question. The remaining implication on this direct route is a cross-order
comparison or convergence theorem for these bounded disorder entropies; no
such result is claimed. This route does not reopen Hadamard, finite census,
or the retired outer-witness continuation.

NUKA ran the two real-algebra cutoff counterexample queries once with
Z3 5.1.0: both are `unsat`, recorded in
`evidence/bernoulli_counting_entropy_20260913/result.json`. That receipt
checks only the explicit scalar constants, not the probability or entropy
proof. The historical September 2 rate criterion is updated to remove its
unnecessary unique-threshold premise.

## 2026-09-13: positive-temperature correction to the gauge precision claim

Section 4.1 of `evidence/NOTE_2026-09-02_UPSTREAM_RELATIVE_GAUGE_BRIDGE.md`
corrects the earlier conversational assertion that the remaining Navier
transfer necessarily requires absolute integer-level occupancy cancellation.
It combines the existing pressure/norm sandwich with the exact partition
normalization `E_g Z_g=Z_A Z_B Z_C`, using all three blocks at the same raw
inverse temperature and including the block sign `tau` in the gauge.

If a selected gauge has model value `ell_g` with
`log rho_g<=ell_g+N epsilon` and
`sum_i log Z_i+ell_g<=beta T`, then
`Phi(Y_g)<=T+N(log 2+epsilon)/beta`. For both optional multiplier rays,
bounded `epsilon` and `beta=[log(n+1)]^2/sqrt(N)` give the explicit
`O(1/log(n+1)^2)` Dini-summable error budget. The leading model inequality
and its remainder bound at these growing scaled temperatures are still
unproved. A fixed-temperature estimate with uncontrolled temperature
dependence does not meet this criterion. No all-orders completion or
unconditional composition improvement is claimed.

NUKA ran one serial, exact polynomial normalization check on one fixed
`3+3` source signing, covering its 32 gauges. The identity held as a Laurent
polynomial; omitting `tau` broke it. Exact rational norm-sandwich checks
also passed. This did not repeat the finite connected-layer or repair
campaigns. Evidence: `evidence/relative_gauge_temperature_20260913/`;
result SHA-256
`70a35b1ff14cb26105ff571cc56d557b0e0fbf82b0a87ac3c8bd964615e98d50`.
Root reviewed the analytic derivation; the finite check is not independent
human review. This supporting correction did not require a fresh major
milestone backup. Current mathematical status remains OPEN; the next
unresolved implication on this route is the pair of inequalities (11).

## 2026-09-13: reverse-KL gauge selection reduction

The same normalized gauge partition density yields a second exact
all-orders sufficient condition, now in Section 4.2 of
`evidence/NOTE_2026-09-02_UPSTREAM_RELATIVE_GAUGE_BRIDGE.md`. Let
`P_beta=rho_g U` be the gauge-tilted law. Then
`I_beta=D(U||P_beta)=-E_U log rho_g`, so a gauge has
`log rho_g<=-I_beta`. Consequently, if
`I_beta>=log Z_A+log Z_B+log Z_C-beta T`, the associated block has
`Phi<=T+N log(2)/beta`.

At `N=r n`, `r=2,3`, and `beta=log(n+1)^2/sqrt(r n)`, this has the explicit
Dini-summable normalized error `r^(3/2)log(2)/log(n+1)^2`; therefore the
two Section 4 ray instances would prove convergence. This is a reduction,
not a bound: `E rho=1` yields only `I_beta>=0`, whereas the target may
require order `N log(n+1)^2` reverse KL. Existing finite connected-layer
checks and covariance estimates do not prove that scale. No finite campaign
was run or restarted for this result.

## 2026-09-13: uniform reverse-KL tripling criterion is impossible

`evidence/NOTE_2026-09-13_UNIFORM_GAUGE_TRIPLING_OBSTRUCTION.md` supersedes
the preceding entry's designation of the uniform reverse-KL bound as an
open target. For a free order-k block and pinned order-m block, it proves
`E_g Phi(Y_g)>=Phi(B)+s a_m+(k-s)a_(m+s)` for every `0<=s<=k`, where
`a_j=E|sum_(i=1)^j epsilon_i|`. Pin a maximizing state, align the first
`s` free spins with the random cross field, then correct the others. The
quadratic interaction within the corrected set cancels exactly under
the global block-sign average; no Gaussian post-update closure is assumed.

At `k=n,m=2n,s=floor(n/2)`, the added cost tends to
`(2+sqrt(5))/(2sqrt(pi))` after division by `n^(3/2)`. CORE's existing
`limsup alpha_n<=1/2` bounds the target allowance by
`(3sqrt(3)-2sqrt(2))/2`. Their difference exceeds `1/100`, with an exact
rational certificate. Uniformly for all cross seeds and optimal internal
blocks, the average norm therefore exceeds the tripling target by at
least `n^(3/2)/200` for all sufficiently large n.

The reverse-KL hypothesis is exactly the uniform average-pressure bound;
at the proposed temperature it fails by order `n log(n+1)^2`. Fixed
`O(N)` slack and the stated positive fractional-moment bounds cannot fix
it. Retire this criterion; do not restart its covariance or positive-moment
variants. Selected good gauges, negative moments, and nonuniform laws are
not excluded by this theorem and are not claimed solved or newly justified.

NUKA checked one new fixed `4+4` strategy fixture: 64 exact conditional
cancellations (24 nonvacuous at fixed block sign), 832 valid strategy
outcomes, and mean score `43/4`. It also checked the rational constant
certificate. Receipt: `evidence/uniform_gauge_tripling_obstruction_20260913/`;
result SHA-256
`850c0cfe49ed64614a7b586dc4b8f4d9939f1d774286737527a2fd3638dc84ab`.
Root reviewed the all-orders proof; the finite computation is corroboration,
not independent human review. The original limit remains OPEN. The next
unresolved implication in the retained gauge route is the selected-gauge
bound in Section 4.1, not the retired uniform reverse-KL condition.
