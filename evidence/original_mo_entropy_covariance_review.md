# Independent entropy/correlation audit

2026-09-05. All deductions below were checked analytically, with no job,
census, or repository edit. They concern second-moment matrices, not
centered covariances unless a zero-mean hypothesis is supplied.

## 1. Every sign law

Let Q be any probability law on d signs, let P_0 be the uniform product
law, and set C=E_Q[bb^T]. Define F_(i-1)=sigma(b_1,...,b_(i-1)) and
m_i=E_Q[b_i|F_(i-1)]. The entropy chain rule and the binary inequality

    ((1+m)/2)log(1+m)+((1-m)/2)log(1-m)>=m^2/2

give D(Q||P_0)>=(1/2)sum_i E m_i^2, including m=+-1 by continuity.
The binary inequality follows by differentiating twice: the left side
has second derivative 1/(1-m^2)>=1 and vanishing value/derivative at 0.

For c_i=(C_ij)_(j<i) and every real unit vector a,

    a^T c_i=E[m_i a^T b_(<i)],
    |a^T c_i|^2<=E m_i^2 (a^T C_(<i,<i) a)
                <=||C||op E m_i^2.

Taking the supremum and summing gives

    sum_i ||c_i||^2=(1/2)||C-I||F^2,
    D(Q||P_0)>=||C-I||F^2/(4||C||op).                    (1)

No zero means, full support of Q, or invertibility of C is required.
Since C is positive semidefinite with diagonal one, ||C||op>=1.

## 2. Gaussian signs, including singular covariance

Let Sigma be any positive-semidefinite matrix with diagonal one and
b=sign Z for Z~N(0,Sigma). Every marginal is a nondegenerate standard
normal, so its sign is defined almost surely even when Sigma is singular.
The Gaussian sign identity, including perfect-correlation endpoints,
gives C=(2/pi)arcsin[entrywise](Sigma). Consequently

    ||C-I||F^2 >= (4/pi^2)||Sigma-I||F^2.                 (2)

Write (2/pi)arcsin x=sum_(k>=0) w_k x^(2k+1), where w_k>=0 and
sum_k w_k=1. For each k, Sigma^[entrywise 2k] is positive semidefinite
with diagonal one by the Schur product theorem; for k=0 it is the
all-ones matrix. Schur multiplication by any such matrix S is a
positive unital linear map on symmetric matrices. Therefore

    -||X||op I<=X<=||X||op I
      implies -||X||op I<=S circ X<=||X||op I,

so it contracts operator norm. Taking X=Sigma gives
||Sigma^[entrywise (2k+1)]||op<=||Sigma||op. The absolutely convergent
entrywise series is also convergent in operator norm in finite dimension,
and hence ||C||op<=||Sigma||op. Combining with (1)--(2) proves

    D(sign N(0,Sigma)||P_0)
      >=||Sigma-I||F^2/(pi^2||Sigma||op).                 (3)

This proof never invokes a Gaussian density KL formula. Thus it remains
valid at singular Sigma, where the Gaussian KL to N(0,I) is infinite.

## 3. The actual centered tensor matrix

For a complete order-n signing A let L=||A||op. With actual phase
energies p=tr(AU)>=0 and q=tr(AV)<=0, set alpha=(p+q)/(2n) and

    H=A tensor A-alpha(A tensor I+I tensor A),
    mu=-lambda_min(H)>0.

The complete proof in Section 3.1 of the integral covariance-rounding
note supplies mu<=L^2 and ||H||op<=2L^2. Explicitly, if the extreme
eigenvalues of A are a>0 and -b<0, then -b/2<=alpha<=a/2 and
mu=ab+alpha(a-b). Its minimum over the allowed alpha interval is
min(a,b)(a+b)/2>0, and its maximum is L(a+b)/2<=L^2.

The tensor terms in H have disjoint entry supports. Equivalently,
their Frobenius cross inner products vanish because tr A=0. Thus

    ||H||F^2=n^2(n-1)^2+2alpha^2 n^2(n-1)
             >=n^2(n-1)^2.                              (4)

For 0<rho<=1, Sigma_rho=I+rho H/mu has diagonal one and minimum
eigenvalue 1-rho>=0. It is therefore admissible in (3), including rho=1.
Keeping the mu dependence in the operator denominator gives

    D(sign N(0,Sigma_rho)||P_0)
      >=rho^2||H||F^2/[pi^2(mu^2+2rho mu L^2)]
      >=rho^2 n^2(n-1)^2/[pi^2(1+2rho)L^4].              (5)

It would be incorrect to assert ||Sigma_rho||op<=1+2rho in general;
the replacement mu<=L^2 is made only after combining the ratio.

## 4. A safe spectral-to-Boolean constant

Set Phi(A)=max_(x signs)|Q_A(x)| with Q_A(x)=x^T A x/2.
The zero diagonal makes Q_A multiaffine, so |Q_A(v)|<=Phi(A) whenever
v lies in the real unit cube. Polarization gives, for real cube vectors
x,y,

    x^T A y=2[Q_A((x+y)/2)-Q_A((x-y)/2)],

and thus the real l-infinity-to-l-one norm is at most 4Phi(A).
Complexifying both vectors and expanding into their real/imaginary
parts bounds the complex l-infinity-to-l-one norm by 16Phi(A).
The complex l-one-to-l-infinity norm is max_ij|A_ij|=1.
Riesz--Thorin interpolation at one half therefore gives

    ||A||op^2<=16Phi(A).                                  (6)

This deliberately conservative constant avoids any real/complex
interpolation convention ambiguity.

If Phi(A)<=C n^(3/2), equations (5)--(6) imply the explicit bound

    D(sign N(0,Sigma_rho)||P_0)
      >=rho^2/[256pi^2(1+2rho)C^2] * n(1-1/n)^2.         (7)

Thus for every fixed rho>0, including rho=1, these canonical Gaussian
sign laws have at least extensive information relative to iid signs.
This rules out applying an o(n)-KL exclusion to them. It does not
prove their selected or mean pressure comparison, nor rule it out.
