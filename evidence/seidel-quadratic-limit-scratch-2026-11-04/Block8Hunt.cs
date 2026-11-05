// Hunt Phi([[A,B],[B^T,-A]]) <= 28 with A an order-8 minimizer (Phi=10).
using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;

public static class Block8Hunt {
    const int N = 8;
    const int Proj7 = 1 << 7; // 128 for x with x0=1
    const int Full8 = 1 << 8; // 256

    static sbyte[,] XA; // [128,8]
    static sbyte[,] YA; // [256,8]

    static void Init() {
        XA = new sbyte[Proj7, N];
        for (int m = 0; m < Proj7; m++) {
            XA[m, 0] = 1;
            for (int i = 1; i < N; i++)
                XA[m, i] = ((m >> (i - 1)) & 1) == 0 ? (sbyte)1 : (sbyte)(-1);
        }
        YA = new sbyte[Full8, N];
        for (int m = 0; m < Full8; m++)
            for (int i = 0; i < N; i++)
                YA[m, i] = ((m >> i) & 1) == 0 ? (sbyte)1 : (sbyte)(-1);
    }

    static int Phi8(sbyte[,] a) {
        int best = 0;
        for (int m = 0; m < Proj7; m++) {
            int q = 0;
            for (int i = 0; i < N; i++)
                for (int j = i + 1; j < N; j++)
                    q += a[i, j] * XA[m, i] * XA[m, j];
            if (q < 0) q = -q;
            if (q > best) best = q;
        }
        return best;
    }

    // Exact Phi of M = [[A,B],[B^T,-A]] over {±1}^16 with first coord +1.
    // Q(x,y) = Q_A(x) - Q_A(y) + x^T B y
    // where x,y in {±1}^8. Using projective: fix (x,y) first nonzero...
    // Full: max over all x,y of |QA(x)-QA(y)+xBy|. Antipodal ( -x,-y) same.
    // We enumerate all 256 x and 256 y (fast).
    static double[] QA_cache;
    static sbyte[,] A_fixed;

    static void PrepA(sbyte[,] A) {
        A_fixed = A;
        QA_cache = new double[Full8];
        for (int m = 0; m < Full8; m++) {
            int q = 0;
            for (int i = 0; i < N; i++)
                for (int j = i + 1; j < N; j++)
                    q += A[i, j] * YA[m, i] * YA[m, j];
            QA_cache[m] = q;
        }
    }

    static int PhiBlock(sbyte[,] B) {
        int best = 0;
        // x index px, y index py
        for (int px = 0; px < Full8; px++) {
            // precompute (B^T y-like) = B row dots: sx[j] = sum_i B[i,j]*x[i] wait x^T B = row vector
            int[] xB = new int[N];
            for (int j = 0; j < N; j++) {
                int s = 0;
                for (int i = 0; i < N; i++)
                    s += B[i, j] * YA[px, i];
                xB[j] = s;
            }
            double qx = QA_cache[px];
            for (int py = 0; py < Full8; py++) {
                int cross = 0;
                for (int j = 0; j < N; j++)
                    cross += xB[j] * YA[py, j];
                int val = (int)(qx - QA_cache[py] + cross);
                if (val < 0) val = -val;
                if (val > best) best = val;
            }
        }
        return best;
    }

    // Find an order-8 minimizer Phi=10 by SA/greedy
    static sbyte[,] FindA(Random rng) {
        sbyte[,] bestA = new sbyte[N, N];
        int best = int.MaxValue;
        for (int trial = 0; trial < 200; trial++) {
            sbyte[,] a = new sbyte[N, N];
            for (int i = 0; i < N; i++)
                for (int j = i + 1; j < N; j++) {
                    sbyte s = rng.Next(2) == 0 ? (sbyte)(-1) : (sbyte)1;
                    a[i, j] = s; a[j, i] = s;
                }
            // greedy edge flips
            bool imp = true;
            int cur = Phi8(a);
            while (imp) {
                imp = false;
                for (int i = 0; i < N; i++)
                    for (int j = i + 1; j < N; j++) {
                        a[i, j] = (sbyte)(-a[i, j]); a[j, i] = a[i, j];
                        int v = Phi8(a);
                        if (v < cur) { cur = v; imp = true; }
                        else { a[i, j] = (sbyte)(-a[i, j]); a[j, i] = a[i, j]; }
                    }
            }
            if (cur < best) {
                best = cur;
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++)
                        bestA[i, j] = a[i, j];
                Console.WriteLine("A-candidate Phi8=" + best);
                if (best <= 10) return bestA;
            }
        }
        Console.WriteLine("Best A Phi8=" + best + " (want 10)");
        return bestA;
    }

    public static int Main(string[] args) {
        Init();
        int seconds = args.Length > 0 ? int.Parse(args[0]) : 300;
        int seed = args.Length > 1 ? int.Parse(args[1]) : 2;
        Random rng = new Random(seed);

        Console.WriteLine("Finding order-8 source A with Phi<=10...");
        sbyte[,] A = FindA(rng);
        int pa = Phi8(A);
        Console.WriteLine("Using A with Phi8=" + pa);
        PrepA(A);

        sbyte[,] B = new sbyte[N, N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                B[i, j] = rng.Next(2) == 0 ? (sbyte)(-1) : (sbyte)1;

        int cur = PhiBlock(B);
        int best = cur;
        sbyte[,] bestB = (sbyte[,])B.Clone();
        Console.WriteLine("Initial block Phi=" + cur);

        double T = 5.0;
        Stopwatch sw = Stopwatch.StartNew();
        long evals = 1;
        while (sw.Elapsed.TotalSeconds < seconds) {
            int i = rng.Next(N), j = rng.Next(N);
            B[i, j] = (sbyte)(-B[i, j]);
            int v = PhiBlock(B);
            evals++;
            int d = v - cur;
            if (d <= 0 || rng.NextDouble() < Math.Exp(-d / Math.Max(T, 1e-9))) {
                cur = v;
                if (v < best) {
                    best = v;
                    bestB = (sbyte[,])B.Clone();
                    Console.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "[BEST] block Phi={0} evals={1} t={2:F1}s", best, evals, sw.Elapsed.TotalSeconds));
                    if (best <= 28) {
                        Console.WriteLine("TARGET_REACHED");
                        break;
                    }
                }
            } else {
                B[i, j] = (sbyte)(-B[i, j]);
            }
            T *= 0.9995;
            if (T < 0.05) {
                T = 4.0;
                // reheat from best
                B = (sbyte[,])bestB.Clone();
                cur = best;
            }
        }
        Console.WriteLine(string.Format(CultureInfo.InvariantCulture,
            "DONE best={0} evals={1} t={2:F1}s", best, evals, sw.Elapsed.TotalSeconds));

        string outDir = @"C:\Users\nicho\seidel-quadratic-limit\results";
        Directory.CreateDirectory(outDir);
        using (var w = new StreamWriter(Path.Combine(outDir, "experiment5_block8_hunt.txt"))) {
            w.WriteLine("best_block_Phi=" + best);
            w.WriteLine("A_Phi8=" + pa);
            w.WriteLine("target_reached=" + (best <= 28));
            w.WriteLine("evals=" + evals);
            w.WriteLine("seconds=" + seconds);
        }
        if (best <= 28) {
            // assemble 16x16
            using (var w = new StreamWriter(Path.Combine(outDir, "experiment5_f16_eq28_witness.json"))) {
                w.Write("{\"Phi\":" + best + ",\"form\":\"[[A,B],[B^T,-A]]\",\"A\":[");
                for (int i = 0; i < N; i++) {
                    if (i > 0) w.Write(",");
                    w.Write("[");
                    for (int j = 0; j < N; j++) { if (j > 0) w.Write(","); w.Write(A[i, j]); }
                    w.Write("]");
                }
                w.Write("],\"B\":[");
                for (int i = 0; i < N; i++) {
                    if (i > 0) w.Write(",");
                    w.Write("[");
                    for (int j = 0; j < N; j++) { if (j > 0) w.Write(","); w.Write(bestB[i, j]); }
                    w.Write("]");
                }
                w.WriteLine("],\"claim\":\"F(16)<=28 hence F(16)=28\"}");
            }
        }
        return best <= 28 ? 0 : 1;
    }
}
