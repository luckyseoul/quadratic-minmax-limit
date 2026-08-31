#!/usr/bin/env python3
"""f(D)(z+Dz)≠0 on g-factors for many odd p. ProcessPool."""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15616 import _g_factors, _w_zDz, _apply_poly  # noqa: E402


def one(p: int) -> dict:
    wfn, q, mul, gen, eigen, inU = _w_zDz(p)
    g, facs = _g_factors(p)
    hits = []
    all_nz = True
    for f in facs:
        r = _apply_poly(wfn, f, mul, gen, q)
        nz = bool(r.max())
        hits.append({"deg": len(f) - 1, "nz": nz})
        all_nz = all_nz and nz
    rec = {
        "p": p,
        "mod4": p % 4,
        "n_fac": len(facs),
        "all_nz": all_nz and len(facs) > 0,
        "hits": hits,
        "eigen": eigen,
        "inU": inU,
    }
    print(f"p={p}≡{p%4} fac={len(facs)} all_nz={rec['all_nz']}", flush=True)
    return rec


def main():
    primes = (3, 5, 7, 11, 13, 17, 19, 23)
    rows = {}
    with ProcessPoolExecutor(max_workers=8) as ex:
        for rec in ex.map(one, primes):
            rows[str(rec["p"])] = rec
    dest = ROOT / "evidence" / "w2_fd_many_p.json"
    dest.write_text(json.dumps({"rows": rows}, indent=2) + "\n")
    print("wrote", dest, flush=True)
    bad = [p for p, r in rows.items() if not r["all_nz"] and r["n_fac"] > 0]
    print("failures", bad, flush=True)


if __name__ == "__main__":
    main()
