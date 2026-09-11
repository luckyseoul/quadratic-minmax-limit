"""Independent table replay for the finite cone-ranked majorizer repair."""
import gzip
import hashlib
import json
import sys
import numpy as np

table_path, result_path = sys.argv[1:]
with gzip.open(table_path, "rt") as handle:
    scan = json.load(handle)
table = np.empty((1 << 16, 2))
for row in scan["rows"]:
    diagonal = np.asarray(row["diagonal"], dtype=float)
    table[row["mask"]] = (row["relative_gap"], diagonal.sum() * (1 / diagonal).sum() / 64 - 1)
with gzip.open(result_path, "rt") as handle:
    result = json.load(handle)
assert result["input_sha256"] == hashlib.sha256(open(table_path, "rb").read()).hexdigest()
assert result["starts"] == 65536 and len(result["rows"]) == 262144
tol = 1e-11
by_mode = {}
for row in result["rows"]:
    assert np.allclose(row["initial_residual"], table[row["initial"]], atol=2e-15, rtol=0)
    assert np.allclose(row["final_residual"], table[row["final"]], atol=2e-15, rtol=0)
    before, after = np.asarray(row["initial_residual"]), np.asarray(row["final_residual"])
    assert after[0] < before[0] - tol or (abs(after[0] - before[0]) <= tol and after[1] <= before[1] + tol)
    by_mode.setdefault(row["mode"], []).append(row)
for mode, rows in by_mode.items():
    before = np.asarray([r["initial_residual"][0] for r in rows])
    after = np.asarray([r["final_residual"][0] for r in rows])
    claimed = result["summary"][mode]
    assert claimed["strict_gap_improvements"] == int(np.sum(after < before - tol))
    assert abs(claimed["mean_gap_reduction"] - np.mean(before - after)) < 1e-14
    assert abs(claimed["best_final_relative_gap"] - np.min(after)) < 1e-14
print(json.dumps({"status": "PASS", "rows": len(result["rows"]), "modes": sorted(by_mode)}))
