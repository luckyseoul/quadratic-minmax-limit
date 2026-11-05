// Exact Phi evaluator + incremental SA hunt for F(16) <= 28.
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;

public static class Order16Hunt {
    const int N = 16;
    const int NEdges = N * (N - 1) / 2; // 120
    const int Proj = 1 << (N - 1);      // 32768

    // Precomputed x_i[mask] for i=0..15 with x_0=+1
    static sbyte[,] X; // [Proj, N]
    static int[] EdgeU = new int[NEdges];
    static int[] EdgeV = new int[NEdges];

    static void InitSpins() {
        X = new sbyte[Proj, N];
        for (int mask = 0; mask < Proj; mask++) {
            X[mask, 0] = 1;
            for (int i = 1; i < N; i++)
                X[mask, i] = ((mask >> (i - 1)) & 1) == 0 ? (sbyte)1 : (sbyte)(-1);
        }
        int k = 0;
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) {
                EdgeU[k] = i;
                EdgeV[k] = j;
                k++;
            }
    }

    static void ScoresFromEdges(sbyte[] edges, int[] scores) {
        for (int mask = 0; mask < Proj; mask++) {
            int q = 0;
            for (int e = 0; e < NEdges; e++)
                q += edges[e] * X[mask, EdgeU[e]] * X[mask, EdgeV[e]];
            scores[mask] = q;
        }
    }

    static int PhiFromScores(int[] scores) {
        int best = 0;
        for (int mask = 0; mask < Proj; mask++) {
            int aq = scores[mask];
            if (aq < 0) aq = -aq;
            if (aq > best) best = aq;
        }
        return best;
    }

    // Flip edge e in place; update scores by delta = -2*old*x_u*x_v
    static void FlipEdge(sbyte[] edges, int e, int[] scores) {
        sbyte old = edges[e];
        edges[e] = (sbyte)(-old);
        int u = EdgeU[e], v = EdgeV[e];
        int twoOld = 2 * old;
        for (int mask = 0; mask < Proj; mask++)
            scores[mask] -= twoOld * X[mask, u] * X[mask, v];
    }

    static sbyte[,] EdgesToMatrix(sbyte[] edges) {
        sbyte[,] a = new sbyte[N, N];
        for (int e = 0; e < NEdges; e++) {
            int u = EdgeU[e], v = EdgeV[e];
            a[u, v] = edges[e];
            a[v, u] = edges[e];
        }
        return a;
    }

    static sbyte[] MatrixToEdges(sbyte[,] a) {
        sbyte[] e = new sbyte[NEdges];
        for (int k = 0; k < NEdges; k++)
            e[k] = a[EdgeU[k], EdgeV[k]];
        return e;
    }

    static sbyte[,] ParseWitnessJson(string path) {
        string text = File.ReadAllText(path);
        int start = text.IndexOf("[[");
        int end = text.IndexOf("]]", start) + 2;
        string body = text.Substring(start, end - start);
        var nums = new List<sbyte>();
        int i = 0;
        while (i < body.Length) {
            char c = body[i];
            if (c == '-' || (c >= '0' && c <= '9')) {
                int j = i + 1;
                while (j < body.Length && body[j] >= '0' && body[j] <= '9') j++;
                nums.Add(sbyte.Parse(body.Substring(i, j - i), CultureInfo.InvariantCulture));
                i = j;
            } else i++;
        }
        if (nums.Count != N * N)
            throw new Exception("expected " + (N * N) + " entries, got " + nums.Count);
        sbyte[,] a = new sbyte[N, N];
        int k = 0;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                a[r, c] = nums[k++];
        return a;
    }

    static void RandomEdges(sbyte[] edges, Random rng) {
        for (int i = 0; i < edges.Length; i++)
            edges[i] = rng.Next(2) == 0 ? (sbyte)(-1) : (sbyte)1;
    }

    static int Greedy(sbyte[] edges, int[] scores) {
        ScoresFromEdges(edges, scores);
        int cur = PhiFromScores(scores);
        bool improved = true;
        while (improved) {
            improved = false;
            for (int e = 0; e < NEdges; e++) {
                FlipEdge(edges, e, scores);
                int v = PhiFromScores(scores);
                if (v < cur) {
                    cur = v;
                    improved = true;
                    if (cur <= 28) return cur;
                } else {
                    FlipEdge(edges, e, scores); // revert
                }
            }
        }
        return cur;
    }

    static int RunSA(sbyte[] bestEdges, int bestPhi, Random rng, int seconds, int[] scores) {
        sbyte[] cur = (sbyte[])bestEdges.Clone();
        ScoresFromEdges(cur, scores);
        int curPhi = PhiFromScores(scores);
        bestPhi = curPhi;
        Array.Copy(cur, bestEdges, cur.Length);

        double T = 4.0;
        Stopwatch sw = Stopwatch.StartNew();
        long evals = 0;
        while (sw.Elapsed.TotalSeconds < seconds) {
            int e = rng.Next(NEdges);
            FlipEdge(cur, e, scores);
            int v = PhiFromScores(scores);
            evals++;
            int d = v - curPhi;
            if (d <= 0 || rng.NextDouble() < Math.Exp(-d / Math.Max(T, 1e-9))) {
                curPhi = v;
                if (v < bestPhi) {
                    bestPhi = v;
                    Array.Copy(cur, bestEdges, cur.Length);
                    Console.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "[BEST] Phi={0} evals={1} t={2:F1}s", bestPhi, evals, sw.Elapsed.TotalSeconds));
                    if (bestPhi <= 28) {
                        Console.WriteLine("TARGET_REACHED");
                        return bestPhi;
                    }
                }
            } else {
                FlipEdge(cur, e, scores); // revert
            }
            T *= 0.99998;
            if (T < 0.04) T = 3.5;
            if ((evals & 2047) == 0 && rng.NextDouble() < 0.03) {
                Array.Copy(bestEdges, cur, cur.Length);
                ScoresFromEdges(cur, scores);
                curPhi = bestPhi;
                T = 1.5;
            }
        }
        Console.WriteLine(string.Format(CultureInfo.InvariantCulture,
            "SA done evals={0} best={1} t={2:F1}s", evals, bestPhi, sw.Elapsed.TotalSeconds));
        return bestPhi;
    }

    public static int Main(string[] args) {
        InitSpins();
        string witness = args.Length > 0 ? args[0]
            : @"C:\Users\nicho\seidel-quadratic-limit\_ref\quadratic-minmax-limit\evidence\m16_phi30_witness_20260911.json";
        int seconds = args.Length > 1 ? int.Parse(args[1]) : 180;
        int seed = args.Length > 2 ? int.Parse(args[2]) : 1;

        int[] scores = new int[Proj];
        Console.WriteLine("Verifying Phi=30 witness...");
        sbyte[,] W = ParseWitnessJson(witness);
        sbyte[] we = MatrixToEdges(W);
        ScoresFromEdges(we, scores);
        int wphi = PhiFromScores(scores);
        Console.WriteLine("Witness Phi=" + wphi + (wphi == 30 ? " OK" : " MISMATCH"));
        if (wphi != 30) return 2;

        Console.WriteLine("Scanning 1-edge neighborhood...");
        int neighBest = wphi;
        for (int e = 0; e < NEdges; e++) {
            FlipEdge(we, e, scores);
            int v = PhiFromScores(scores);
            if (v < neighBest) neighBest = v;
            FlipEdge(we, e, scores);
        }
        Console.WriteLine("Best 1-flip neighbor Phi=" + neighBest);

        Random rng = new Random(seed);
        sbyte[] bestEdges = (sbyte[])we.Clone();
        int best = wphi;

        Console.WriteLine("Multi-start greedy (80 trials)...");
        for (int trial = 0; trial < 80; trial++) {
            sbyte[] edges = new sbyte[NEdges];
            if (trial == 0) Array.Copy(we, edges, NEdges);
            else if (trial < 10) {
                Array.Copy(we, edges, NEdges);
                // random k flips from witness
                int kflips = 3 + rng.Next(8);
                for (int t = 0; t < kflips; t++) {
                    int ee = rng.Next(NEdges);
                    edges[ee] = (sbyte)(-edges[ee]);
                }
            } else RandomEdges(edges, rng);
            int g = Greedy(edges, scores);
            if (g < best) {
                best = g;
                Array.Copy(edges, bestEdges, NEdges);
                Console.WriteLine("[GREEDY BEST] Phi=" + best + " trial=" + trial);
                if (best <= 28) break;
            }
        }

        if (best > 28) {
            Console.WriteLine("SA for " + seconds + "s from best...");
            best = RunSA(bestEdges, best, rng, seconds, scores);
        }

        string outDir = @"C:\Users\nicho\seidel-quadratic-limit\results";
        Directory.CreateDirectory(outDir);
        using (var sw = new StreamWriter(Path.Combine(outDir, "experiment5_order16_hunt.txt"))) {
            sw.WriteLine("best_Phi=" + best);
            sw.WriteLine("witness_Phi=" + wphi);
            sw.WriteLine("one_flip_best=" + neighBest);
            sw.WriteLine("target_reached=" + (best <= 28));
            sw.WriteLine("seconds_budget=" + seconds);
            sw.WriteLine("seed=" + seed);
            sw.Write("edges=");
            for (int i = 0; i < bestEdges.Length; i++) {
                if (i > 0) sw.Write(",");
                sw.Write(bestEdges[i]);
            }
            sw.WriteLine();
        }
        // Also dump matrix JSON if target hit
        if (best <= 28) {
            var a = EdgesToMatrix(bestEdges);
            using (var sw = new StreamWriter(Path.Combine(outDir, "experiment5_f16_eq28_witness.json"))) {
                sw.Write("{\"Phi\":" + best + ",\"A\":[");
                for (int i = 0; i < N; i++) {
                    if (i > 0) sw.Write(",");
                    sw.Write("[");
                    for (int j = 0; j < N; j++) {
                        if (j > 0) sw.Write(",");
                        sw.Write(a[i, j]);
                    }
                    sw.Write("]");
                }
                sw.WriteLine("],\"claim\":\"F(16)<=28, hence F(16)=28 by parity+monotonicity\"}");
            }
        }
        Console.WriteLine("FINAL best_Phi=" + best + " target_reached=" + (best <= 28));
        return best <= 28 ? 0 : 1;
    }
}
