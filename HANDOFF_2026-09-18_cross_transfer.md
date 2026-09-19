# Handoff 2026-09-18: cross-transfer / C_cross sweep (side investigation)

STATUS: PRELIMINARY. Measurements, not proofs. Off the main frontier.
Do not cite as progress toward the original limit unless re-anchored.

## 1. Where the REAL frontier is (not this work)
From evidence/NOTE_2026-09-01_ORIGINAL_LIMIT_TWO_RAY.md (2026-09-01) and
solution.md (through 15.774, 2026-09-04):
  * Two-ray criterion PROVED (Prop 6.5-6.6). Not open.
  * Multiplier two: narrowed to explicit RESIDUE (6.20), NOT closed.
    Residue = pairs (x,y) with h_x,h_y>rho*n, d_H(x,y)(n-d_H)>rho*n^2,
      h_x h_y > (rho/4) n^2, |Qx+Qy| > 2*sqrt(2)*M - n^(3/2),
    rho = a^2, a=(sqrt(2)-1)/pi. If M/n^(3/2) < 1/(2*sqrt(2))=0.3535 the
    last term is vacuous (low-alpha regime).
  * Multiplier three / 1:2 split: OPEN. Prop 6.7 (tripling diamond), Prop 6.8
    (1:2 reduction, residual 6.42-6.43).
  * "New live gate": close residue (6.20), then prove the 1:2 split.

THIS WORK WAS NOT ON THAT FRONTIER. I was re-measuring the cross term.
Flagged by the user as possible duplication. Reassess before continuing.

## 2. What was done here (side investigation)
Direct measurement of the dyadic transfer constant, NOT via the residue
(6.20) construction:
  2n object: max_{x,y} | x^T Z y + Qx - Qy |   (real Gaussian Z, (alpha,rho)
  3n object: analogous 3-way (script written, NOT completed).
  Z = sqrt(k)*G + sqrt(1-k)*W, k=2/pi, G~N(0,Sig), Sig=I+rho*H_alpha/mu,
  H_alpha = A^otimes A - alpha(A^otimes I + I^otimes A), mu=ab+alpha(a-b).
  A = conference matrix (order q, q prime 1 mod 4). alpha=a=lambda_max(A).

## 3. Findings (soulkiller, empirical)
  A-term at cross-maximizing pair: |Qx*-Qy*| = O(n).   (mo_aterm_verify.py)
    n=13: |Qx*-Qy*|/n p95=1.23 max=1.46.  n=17: p95=1.34 max=1.59.
  C_cross (2n cross floor, sign field ~1.0; real Gaussian field ~1.4):
    (mo_cross_struct.py) real field, conference A:
      n=13 alpha=0 rho=0.5: 0.96-0.98;  alpha=a rho=0.99: ~0.92.
      (C_real can be pushed <0.92 at (alpha=a, rho=0.995)).
  2n FULL transfer at (alpha=a, rho=0.995), n=13, 20 realizations
    (mo_fullmax.py): max_{x,y}|xZy+Qx-Qy|/n^1.5 = mean 0.694, p5 0.579,
    MAX 0.898 < 0.9216 = 2*sqrt(2)*B, B=0.32584.  all< threshold.
  3n transfer: NOT completed (mo_3n.py timed out at n=9).

## 4. Why it might NOT be duplicate
The 2n transfer closing at (a,0.995) (max 0.898 < 0.9216) would close
multiplier two DIRECTLY, bypassing residue (6.20). That is potentially NEW,
but only at n=13, empirical, real-Gaussian Z. Needs:
  (a) verification at larger n (n>13 is intractable by full cube; need a
      different bound or accept n=13 as a signal only),
  (b) conversion to a proof,
  (c) the matching 3n transfer (or 1:2 split) for the full two-ray criterion.

## 5. Concerns / risks
  * Duplication: residue (6.20) already addresses multiplier two. This sweep
    may re-measure the same object the construction handles. Check
    solution.md §10 (bilinear floor) and the shield note before continuing.
  * n=13 is small. 0.898 vs 0.9216 margin is thin (2.5%). Could close at
    larger n.
  * These are measurements, not proofs.

## 6. Files (soulkiller /home/nick/)
  mo_xterm_test.py      A-term at cross-max pair.
  mo_cross_const.py     C_cross (sign field) vs n.
  mo_cross_struct.py    C_cross (structured field) alpha/rho.
  mo_cross_sweep.py     C_cross sweep alpha/rho.
  mo_cross_fine.py      fine alpha/rho sweep.
  mo_3n.py              3n transfer (timed out).
  mo_fullmax.py         full 2n max.
  mo_aterm_verify.py    A-term verify at (a,0.995).

## 7. Next (pick ONE)
  A) Verify 2n transfer at n>13 (new method, not full cube) + 3n transfer.
  B) Stop this; return to frontier: close residue (6.20) + 1:2 split.
  C) Reassess duplication against solution.md §10 and shield note first.
