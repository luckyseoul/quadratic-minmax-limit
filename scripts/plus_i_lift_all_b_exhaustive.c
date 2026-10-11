/*
 * Exhaustive minimum of Phi over all plus-I lifts of order 2m.
 *
 * For a symmetric zero-diagonal +-1 matrix B of order m put
 *     K(B) = [[ B , B+I ], [ B+I , -B ]],
 * a symmetric zero-diagonal +-1 matrix of order n = 2m (the diagonal of the off-diagonal
 * blocks is +1).  Since K(DBD) = (D (+) D) K(B) (D (+) D) for diagonal +-1 matrices D, Phi(K(B)) is
 * invariant under switching B, and B may be taken with first row and column +1.  The program
 * runs over all 2^((m-1)(m-2)/2) such B and counts those with Phi(K(B)) <= T, where
 * Phi(A) = max_x |sum_{i<j} A_ij x_i x_j| is evaluated exactly over all 2^(n-1) classes
 * {x,-x} by Gray-code enumeration with early abort once |Q| > T.
 *
 * Phi(K(B)) is congruent to C(n,2) mod 2.  Runs used in the note:
 *   m=7  (n=14, C(14,2) odd):  T=23 gives 0 matrices; T=25 gives a positive count  => min 25.
 *   m=8  (n=16, C(16,2) even): T=30 gives 0;           T=32 gives a positive count  => min 32.
 *   m=9  (n=18, C(18,2) odd):  T=31 gives 0;           T=33 gives a positive count  => min 33.
 * (With the parity, count 0 at T = v-2 and a positive count at T = v give the minimum v.)
 *
 * Build:  gcc -O2 -fopenmp -o plus_i_lift_all_b_exhaustive plus_i_lift_all_b_exhaustive.c
 * Run:    ./plus_i_lift_all_b_exhaustive m T
 *         seconds for m=7,8; about 40 s (T=31) and 100 s (T=33) for m=9 on 8 cores.
 * Only the C standard library is used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: %s m T\n", argv[0]); return 1; }
  int m = atoi(argv[1]);
  long T = atol(argv[2]);
  int n = 2 * m;
  if (m < 3 || m > 10) { fprintf(stderr, "3 <= m <= 10\n"); return 1; }
  int mm = m - 1, E = mm * (mm - 1) / 2, ei[64], ej[64], e = 0;
  for (int i = 1; i < m; i++)
    for (int j = i + 1; j < m; j++) { ei[e] = i; ej[e] = j; e++; }
  uint64_t total = 1ULL << E;
  long cnt = 0, minphi = 1000000;
#pragma omp parallel for schedule(dynamic, 2048) reduction(+ : cnt)
  for (uint64_t mask = 0; mask < total; mask++) {
    int B[20][20], A[40][40];
    for (int i = 0; i < m; i++) B[i][i] = 0;
    for (int j = 1; j < m; j++) B[0][j] = B[j][0] = 1;
    for (int k = 0; k < E; k++) {
      int s = ((mask >> k) & 1) ? -1 : 1;
      B[ei[k]][ej[k]] = B[ej[k]][ei[k]] = s;
    }
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++) {
        A[i][j] = B[i][j];
        A[i][j + m] = B[i][j] + (i == j);
        A[i + m][j] = B[i][j] + (i == j);
        A[i + m][j + m] = -B[i][j];
      }
    int x[40], h[40];
    long Q = 0;
    for (int i = 0; i < n; i++) {
      x[i] = 1;
      h[i] = 0;
      for (int j = 0; j < n; j++) h[i] += A[i][j];
    }
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++) Q += A[i][j];
    long mx = labs(Q);
    int ok = (mx <= T);
    uint64_t tot = 1ULL << (n - 1);
    for (uint64_t g = 1; g < tot && ok; g++) {
      int i = __builtin_ctzll(g);
      Q -= 2L * x[i] * h[i];
      int xo = x[i];
      x[i] = -xo;
      for (int j = 0; j < n; j++)
        if (j != i) h[j] -= 2 * xo * A[j][i];
      long a = labs(Q);
      if (a > mx) mx = a;
      if (mx > T) ok = 0;
    }
    if (ok) {
      cnt++;
#pragma omp critical
      {
        if (mx < minphi) minphi = mx;
      }
    }
  }
  printf("m=%d (n=%d) T=%ld count(Phi(K(B))<=T)=%ld min Phi(K(B)) among them=%ld\n", m, n, T, cnt,
         cnt ? minphi : -1L);
  return 0;
}
