# Root contributing review: full-weighted actual-cross gain

2026-09-06. PASS for the precise conditional theorem. This is a
contributing review, not an independent review of the entire derivation.

I read the complete final 364-line source end to end at SHA256
ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9,
including a full reread before integration. The complete 411-line
Gaussianization prerequisite had already been read and checked in this
same task sequence at its frozen hash b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20.
The contextual normalization and compatibility arguments are rederived
in the new source, and I checked those derivations directly.

## Normalization and transfer

The row-square sum gives dbar^2>=N-1, and the separate trace cap gives
m_0>=1/(2C^2). Literal cross sign squares factor m=m_0 ell h.
Half-wise Cauchy--Schwarz and t_L+t_R=2 give ell h>=1; their
fixed sum gives ell h<=(1+delta)^2. These directions are correct.

At most b_0 N bad coordinates occur over both halves. Balancing the
good-coordinate counts gives theta<=2b_0, q>=n/2 and
m'=q/[(1+eta)dbar]^2>=1/(9C^2). The retained complete sign matrix
alone has operator norm at most (1+eta)dbar. Unbiased extension proves
beta(B)>=beta(B_J); no operator cap on the original B is assumed.

The symmetric dilation loses 2(n-q) coordinates. Positive eigenvalue
power sums and principal interlacing give normalized loss at most
2theta. Congruence perturbation costs 4k a eta, for k=1,2.
Both actual full moments, including exceptional coordinates, survive.

The imported norm lower carries the exact scale factor a(1+eta).
Its fourth-moment bracket is truncated at zero BEFORE its coefficient
is reduced using m'<=m. Thus a negative error is never multiplied by
an enlarged uncontrolled coefficient. The gain contribution is exactly
g_kappa a^2 m_0/(1+eta). With eta=delta^(1/3)<=1/8, the losses
are at most 25kappa C^2 eta and 6g_kappa eta. The tail supremum
over integer q>=ceil(n/2) gives a finite uniform error tending to zero.
The delta=0 case is treated directly, without a singular trimming limit.

For every feasible D the exact Boolean rescaling argument gives
|u_D-c/(n dbar)|<=2sqrt(delta). The substitution c=beta(B)=Phi(K)
separately requires the actual original pure-cross active state.
This does not follow from spectral moments alone.

## Roles, completed receipts, and scope

I proposed balanced good-coordinate trimming and returning its first
and fourth moments to the original W. The proof worker derived and
authored the finite factors and error bounds. The exact worker checked
the candidate formulas and final source; it had earlier prerequisite
contributions. The docs worker independently reviewed this entire new
transfer. I read all three final receipts completely and verified their
source hashes. No correction remains requested.

No mathematical run, test, scan, checker, solver, or numerical spectral
evaluation was used. The proof retains the full actual weighted law;
the auxiliary matrix is only an original norm lower bound, never a
replacement source, covariance, field, measure, or active cell.

The theorem removes the extra global unweighted operator hypothesis
inside the fixed-cap, near-scalar branch. It does not establish small
dispersion for all conditional optimizers, pure-cross activity for every
cell, the remaining all-law upper comparison, or original convergence.
All these original-problem conclusions remain OPEN.
