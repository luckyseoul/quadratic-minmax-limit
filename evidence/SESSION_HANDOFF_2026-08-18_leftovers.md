# Leftover status (2026-08-18)

E(1) for every prime \(p\ge5\), then \(L=\lim\alpha_n=1/2\) by denseness. Not done. Sandwich and Paley \(\rho=1\) are proved. `e1_closed_general` is True only by the old wiring. Public \(L\) is OPEN.

Repo: https://github.com/luckyseoul/quadratic-minmax-limit · `main` · HEAD 15.585 (`0f50da5`).

## Leftovers

| Item | Flag | Status |
|---|---|---|
| \(\lambda_{\min}(\Phi\rvert_F)\ge6\Leftrightarrow\langle\delta,\psi\rangle\le2\) | `phi_F_ge_6_proved_general=False` | Open. Name ensemble \(Q_\tau\) or \(D=\lvert H_+\rvert/(2p)\), or prove the pairing bound without that name. |
| Residual (ii), even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | Open. Affine + even \(k\le4p-2\) closed. Official class is leftover+\(s_+\ge2\). |
| Type I, Max− not two-level | `type_I_multilevel_bad_case_ND_closed=False` | Open. Remainder \(A_{\mathrm{full}}\). |
| Lemma D | True | Closed. |

Aut-Schur, Gsum disj LB, and the cotangent pairing stay False.

## Named identities (do not re-derive)

| Unit | Identity |
|---|---|
| 15.550 | \(S(\lambda)=\mathrm{Kl}(1,\lambda^2/4)\) on every odd \(q\) |
| 15.564 | \(F=-2(3p^2+2)\); \(Q_{3,02}=-4N(2p^2+1)/p\) |
| 15.561 | \(n_{k=3}=\binom{m}{3}(p-1)q\); \(A_{k=3}=0\) on 1-line triples |
| 15.562 | \(\sum\omega^{2t}/(\omega^t-1)^3=-(p^2-1)/24\); \(A_{k=3,n3}=-16p^3/((p-1)(p-3))\) |
| 15.573 | Exclusive 1D / \(k=3\) / full mix reconstructs live \(Q\) at \(p=5,7\) |
| 15.574 | \(\mu_{k=3}=96p^4 P(r)/(p^2-1)\) |
| 15.575 | \(\mu_{1d}=2p^4(p^2-3p-2)/(p-2)\) |
| 15.578 | \(2\chi\) fourth moment \(4p/(p-1)\); \(p=7\) occupancy mix is not a \(p\)-law |
| 15.581 | \(p=5\): \(Q_{++}/q^2=48/13<26/7\) from named \(\mu_{1d}\), \(\mu_{k=3}\), \(n_{1d}\), \(n_{k=3}\) |
| 15.582 | 1D 4-point vanishes off \(\mathbb F_p\); \(p=5\) \(Q_{N^*}/q^2=32/13\) |
| 15.507 | \(p\equiv1\): \(J_{N^*}=2\) pairing \(\le2\) iff \(Q_{++}/q^2\le Q_{\mathrm{pp,ub}}(p)\) (\(26/7\) at \(p=5\)); automatic for \(p\ge13\) if \(M\ge0\) i.e. \(Q_{++}\le4q^2\), after two-type constancy |

At \(p=5\) (\(m=3\)) exclusive \(n_{\mathrm{full}}=0\). At \(p=7\), \(n_{\mathrm{full}}=90q\) and \(\mu_{\mathrm{full}}\) is a mix, not a formula in \(p\). Do not interpolate \((p-5)/15\).

## Facts that block short proofs of the floor

- Pointwise \(Q_y^{++}\le4q^2\) fails on a positive fraction of Max+ at \(p=5,7\). Identities true for every conference eigenvector (\(z_i^2=1\), \(Cy=py\)) cannot force the ensemble bound.
- Paley×norm types split into many Frob-inv orbits at \(p\ge11\). Two-type \(Q\)-constancy is certified only for \(p\le7\).
- Frob orbit masses on \(H_+\) are \(1/\lvert H_+\rvert\) or \(2/\lvert H_+\rvert\). That names \(D\in\{13,409\}\), not a polynomial in \(p\).
- The exclusive mix cannot prove \(Q_{++}\le4q^2\): 1D+\(k=3\)-only is \(4.68>4\) at \(p=7\); \(\mu_{k=3}/q^2>4\) at \(p\ge11\).

## Residual (ii)

- \(p=5\), \(k=20\): leftover+\(s_+\) empty for all \(n_F\) (15.528). leftover-only \(n_F=8\) exists.
- \(p=5\), \(k=22\): leftover+\(s_+\) empty for \(n_F=0,3\)–\(9,11\)–\(14\); \(n_F=10\) unfinished.
- 15.566: two-value leftover+\(s_+\) at \(k=4p+2\) only \(\{2,6\}\) and \(p=5\) \(\{2,8\}\).
- 15.585: leftover+\(s_+\) at \(k=4p\) forces \(\min_+=2\); \(\{2,4,6\}\) is not a plus pair-slice.
- No general emptiness at even \(k\ge4p\).

## Type I (dead as a multi-level kill)

| Unit | Why dead |
|---|---|
| 15.559 Aut\(_e\) | Does not preserve full-\(\Omega\) |
| 15.565 Max± of \(C\) | Writes at most \(\Phi-4\) |
| 15.577 Type+ 1D Johnson | \(p=5\) still has \(\lvert G\rvert=13\) shadows |
| 15.580 Galois \(\Omega_-\) plus \(F\) | Names the 02-sum, not \(f_e\) on \(U_-\) |
| \(\lvert\mu\rvert\le\lvert L\rvert\) on \(\lvert\kappa\rvert=1\) | Unsigned \(\lvert\nu_{\mathrm{part}}\rvert\) exceeds \(\lvert L\rvert\) |

## Do not reopen

Occupancy, Aut-involution counts for \(D\), half-net census as a \(p\)-law, Aut\(_e\) as \(A_{\mathrm{full}}\), \((p-5)/15\), \(10p-46\), \(16(p-4)/D\), exclusive mix as \(Q_\tau\), pointwise Wick as the floor.

Do not commit untracked 15.495 catalogs, 15.496, or 15.530.

## Next

A Gauss/Jacobi (or equivalent) formula for ensemble \(Q_\tau\) or for \(D\), failing when the formula is wrong; or a pairing proof of \(\langle\delta,\psi\rangle\le2\) that does not use those names. Import `phi_F_ge_6` only then.
