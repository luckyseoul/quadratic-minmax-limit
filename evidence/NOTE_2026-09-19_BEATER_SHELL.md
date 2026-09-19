# Paley-beaters live in a spectral shell; the n=26 undercutter is not a matching

**Status:** two proved lemmas and one certified finite structure. Does not
close \(\gamma_p\to0\). Original limit OPEN.

## 1. Gaussian saturation at conference

Let \(n=p^2+1\) and let \(C\) be the Paley conference matrix, so
\(C^2=(n-1)I\), \(\|C\|_{\mathrm{op}}=\sqrt{n-1}\), and
\(\Phi(C)=\tfrac12 n\sqrt{n-1}\).

**Lemma.** In the CORE §4 Gaussian comparison, conference forces every
pair quantity \(u_{ij}=0\). The arcsine difference is exactly
\(2\arcsin(1/\sqrt{n-1})\), and the resulting lower bound is
\[
\Phi(C)\ge\frac{n\sqrt{n-1}}\pi\cdot\frac{\arcsin v}{v},\qquad
v=\frac1{\sqrt{n-1}}.
\]
The extra factor tends to \(1\). Conference is the worst case of that
estimate, not a place it improves toward \(\tfrac12\). The same Gaussian
cannot prove \(\gamma_p\to0\).

Proposition 15.5 still gives \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\) for
every Seidel matrix, with equality only at conference. Proposition 15.7
still says a beater must pay in \(\rho\). Replacing \(\sqrt{n-1}\) by
\(\|A\|_{\mathrm{op}}\) in CORE §4 is *not* claimed: the coordinate
variances of \((I\pm A/\lambda)g\) change when \(\lambda\ne\sqrt{n-1}\).

## 2. Gaussian saturation at conference

For a conference matrix the off-diagonal of \(A^2\) vanishes, so the
CORE pair quantities \(u_{ij}\) are identically zero. The arcsine
difference is then exactly \(2\arcsin(1/\sqrt{n-1})\), and
\[
\Phi(C)\ge\frac{n\sqrt{n-1}}\pi\cdot\frac{\arcsin v}{v},\qquad v=\frac1{\sqrt{n-1}}.
\]
As \(n\to\infty\) the extra factor tends to \(1\). Conference is the
*worst* case of the Gaussian lower bound, not a place it can be
improved. This method cannot prove \(\gamma_p\to0\).

## 3. The n=26 undercutter is a distant almost-conference

The verified order-26 record \(A\) with \(\Phi(A)=61<65=\Phi(C)\)
(`evidence/ns_port_n26_undercut_20260912/stage01.npz`) satisfies:

- best switching Hamming distance to Paley \(C_5\): \(122\) of \(325\)
  edges (\(\approx 0.92\,n^{3/2}\), not \(O(n)\));
- \(A^2=25I+\Delta\) with \(\Delta_{ij}\in\{0,\pm4\}\) and exactly six
  nonzero off-diagonal entries per row;
- \(\|A\|_{\mathrm{op}}\approx5.849>5\), \(\rho(A)\approx0.802\);
- product \(\rho\|A\|_{\mathrm{op}}\approx4.69<5\).

It is not a matching flip of Paley. Lipschitz therefore does not put
this undercut in the \(o(n^{3/2})\) class. The remaining problem on
\(n=p^2+1\) is to rule out or construct an infinite family of matrices
in the conference-op shell whose Boolean alignment stays a definite
factor below \(1\).

Max+-only discrepancy cannot do this: at \(n=26\), CP-SAT can drop
every Boolean \(+p\) eigenvector by \(64\) with \(160\) flips, after
which the cube maximum *rises* to \(107\).

OpenAI referee `math_review` of the Gaussian-saturation lemma: PASS.
