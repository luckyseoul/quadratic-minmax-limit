#!/usr/bin/env python3
"""Cruise-control orchestrator for the quadratic-minmax research program.

The controller does scheduling/state management only. Mathematical claims are
never accepted from agent prose; results remain evidence until independently
reviewed. Configure the existing worker/agent launcher with
AUTOPILOT_WORKER_CMD. It is invoked with TASK_JSON and RESULT_JSON environment
variables set and the task JSON path as argv[1].
"""
from __future__ import annotations
import hashlib, json, os, shlex, sqlite3, subprocess, sys, time, uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "autopilot" / "config.json"
DEFAULT_ROUTES = {
    "scaling": [
        "Find an unconditional asymptotic comparison for H(n)=m_n^(2/3) that can force convergence; prioritize composition/restriction inequalities and quantify every error.",
        "Attack H(2n)-2H(n) and H(3n)-3H(n) using the actual minimax definition; either prove a useful bound or give a rigorous obstruction.",
    ],
    "restriction_lifting": [
        "Determine whether selected-half restriction can be sharpened beyond the sqrt(2) normalized loss; seek an all-orders statement, not a finite pattern.",
        "Use retained optimal/near-optimal matrices to discover a deterministic lifting or restriction invariant, then formulate the weakest plausible universal lemma.",
    ],
    "structure_mining": [
        "Mine retained good matrices for an invariant shared across orders and independent construction families; translate any stable pattern into a quantified lemma.",
        "Analyze principal restrictions, switching classes, spectra, row sums, cuts, A^2, and active Boolean maximizers of retained winners; look for a universal structural consequence of near-optimality.",
    ],
    "spectral_rounding": [
        "Find an unconditional spectral/Boolean-rounding inequality stronger than the current lower bounds, without assuming a limiting spectral law or optimizer family.",
        "Stress-test actual-source Gaussian/Hermite rounding arguments and identify the weakest missing hypothesis that would make them all-orders.",
    ],
    "counterexample": [
        "Try to construct two infinite families with provably separated normalized minimax values; a finite numerical oscillation is not sufficient.",
        "Search for a mechanism that makes alpha_n fail to converge despite monotonicity of m_n and vanishing adjacent alpha gaps.",
    ],
    "proof_synthesis": [
        "Read all current route-neutral results and synthesize a genuinely unconditional convergence lemma from combinations that are logically compatible; expose every unproved bridge.",
        "Attempt to replace a conditional route hypothesis with a theorem derived from the original minimax definition.",
    ],
    "adversarial": [
        "Adversarially audit the strongest recent candidate claim for quantifier errors, finite-to-asymptotic leakage, hidden optimality assumptions, and exceptional subsequences.",
    ],
    "exploration": [
        "Develop a substantially different route to convergence or nonconvergence using the original definition; do not recycle a rejected obstruction under a new name.",
    ],
}


def load_cfg():
    return json.loads(CFG.read_text())


def db_init(db):
    db.executescript("""
    CREATE TABLE IF NOT EXISTS tasks(
      id TEXT PRIMARY KEY, route TEXT NOT NULL, prompt TEXT NOT NULL,
      parent TEXT, status TEXT NOT NULL, priority REAL NOT NULL,
      attempts INTEGER NOT NULL DEFAULT 0, created REAL NOT NULL,
      started REAL, finished REAL, result_path TEXT, host TEXT);
    CREATE TABLE IF NOT EXISTS results(
      id TEXT PRIMARY KEY, task_id TEXT NOT NULL, status TEXT, claim TEXT,
      implication TEXT, scope TEXT, raw_path TEXT, digest TEXT, created REAL);
    CREATE TABLE IF NOT EXISTS routes(
      route TEXT PRIMARY KEY, reward REAL NOT NULL DEFAULT 0,
      attempts INTEGER NOT NULL DEFAULT 0, successes INTEGER NOT NULL DEFAULT 0,
      failures INTEGER NOT NULL DEFAULT 0, last_change REAL);
    CREATE TABLE IF NOT EXISTS meta(k TEXT PRIMARY KEY, v TEXT NOT NULL);
    """)


def digest_file(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def add_task(db, route, prompt, parent=None, priority=1.0):
    canonical = json.dumps({"route": route, "prompt": prompt, "parent": parent}, sort_keys=True)
    tid = hashlib.sha256(canonical.encode()).hexdigest()[:20]
    if db.execute("SELECT 1 FROM tasks WHERE id=?", (tid,)).fetchone():
        return tid
    db.execute("INSERT INTO tasks(id,route,prompt,parent,status,priority,created) VALUES(?,?,?,?,?,?,?)",
               (tid, route, prompt, parent, "queued", priority, time.time()))
    return tid


def seed(db):
    for route, prompts in DEFAULT_ROUTES.items():
        db.execute("INSERT OR IGNORE INTO routes(route) VALUES(?)", (route,))
        for p in prompts:
            add_task(db, route, p, priority=1.0)
    db.commit()


def route_score(db, route):
    r = db.execute("SELECT reward,attempts FROM routes WHERE route=?", (route,)).fetchone()
    reward, attempts = r
    # UCB-ish allocation: reward per attempt plus uncertainty bonus.
    return (reward / max(1, attempts)) + 0.75 / (1 + attempts) ** 0.5


def choose(db, tags):
    rows = db.execute("SELECT id,route,prompt,priority FROM tasks WHERE status='queued'").fetchall()
    if not rows:
        return None
    # Keep a little exploration even when one route dominates.
    rows.sort(key=lambda x: x[3] * route_score(db, x[1]), reverse=True)
    for row in rows:
        route = row[1]
        if route == "counterexample" or route == "exploration" or tags is None:
            return row
        if any(t in tags for t in ("cpu", "cuda", "rocm")):
            return row
    return rows[0]


def make_context():
    files = ["CORE.md", "STATUS.md", "HANDOFF.md", "ARTIFACTS.md"]
    ctx = []
    for f in files:
        p = ROOT / f
        if p.exists():
            # Workers can read the repository themselves; this gives them a small provenance pointer.
            ctx.append(str(p.relative_to(ROOT)))
    return ctx


def build_task(cfg, row, result_path):
    tid, route, prompt, priority = row
    return {
        "schema": "quadratic_minmax_autopilot_task_v1",
        "task_id": tid,
        "route": route,
        "objective": prompt,
        "repository": str(ROOT),
        "read_first": make_context(),
        "rules": [
            "The actual goal is convergence or nonconvergence of alpha_n.",
            "Finite computation is evidence, never an all-orders proof.",
            "Do not assume L=1/2, Paley optimality, optimizer classification, or a route-specific necessity claim.",
            "Reuse existing certificates and do not rerun unchanged searches.",
            "State the exact implication your work would establish.",
            "If a claim is promising, identify assumptions a destroyer must attack.",
        ],
        "required_result": str(result_path),
        "time_budget_seconds": cfg["task_timeout_seconds"],
    }


def validate_result(p):
    if not p.exists() or p.stat().st_size == 0:
        raise ValueError("worker produced no result JSON")
    x = json.loads(p.read_text())
    required = ["status", "claim", "implication", "scope"]
    missing = [k for k in required if k not in x]
    if missing:
        raise ValueError("missing result fields: " + ",".join(missing))
    if x["status"] not in {"success", "counterexample", "partial", "dead_end", "error"}:
        raise ValueError("invalid result status")
    return x


def ingest(db, task_id, result_path, x):
    rid = hashlib.sha256((task_id + digest_file(result_path)).encode()).hexdigest()[:24]
    db.execute("INSERT OR REPLACE INTO results VALUES(?,?,?,?,?,?,?,?,?)",
               (rid, task_id, x["status"], str(x["claim"]), str(x["implication"]),
                str(x["scope"]), str(result_path), digest_file(result_path), time.time()))
    positive = x["status"] in {"success", "counterexample"}
    reward = 2.0 if positive else 0.25 if x["status"] == "partial" else -0.75
    route = db.execute("SELECT route FROM tasks WHERE id=?", (task_id,)).fetchone()[0]
    db.execute("UPDATE routes SET reward=reward+?,attempts=attempts+1,successes=successes+?,failures=failures+?,last_change=? WHERE route=?",
               (reward, int(positive), int(not positive), time.time(), route))
    # Positive claims automatically generate a destroyer; partials generate a focused continuation.
    if positive:
        add_task(db, "adversarial", f"Destroy or validate this candidate claim from task {task_id}: {x['claim']}. Exact implication: {x['implication']}. Attack every dependency and quantifier; use independent computation where useful.", parent=task_id, priority=2.0)
    if x.get("next_tasks"):
        for nt in x["next_tasks"][:6]:
            if isinstance(nt, str) and nt.strip():
                add_task(db, route, nt, parent=task_id, priority=1.25 if positive else 0.9)
    db.execute("UPDATE tasks SET status='done',finished=?,result_path=? WHERE id=?", (time.time(), str(result_path), task_id))
    db.commit()


def run_one(cfg, db, row):
    tid, route, prompt, priority = row
    host = cfg["hosts"][0]
    run_dir = ROOT / cfg["artifact_dir"] / tid
    run_dir.mkdir(parents=True, exist_ok=True)
    task_path = run_dir / "task.json"
    result_path = run_dir / "result.json"
    task_path.write_text(json.dumps(build_task(cfg, row, result_path), indent=2))
    cmd = os.environ.get("AUTOPILOT_WORKER_CMD")
    if not cmd:
        raise RuntimeError("AUTOPILOT_WORKER_CMD is not set; no worker was launched")
    env = os.environ.copy()
    env.update(cfg.get("worker_env", {}))
    env["AUTOPILOT_TASK_JSON"] = str(task_path)
    env["AUTOPILOT_RESULT_JSON"] = str(result_path)
    env["AUTOPILOT_ROUTE"] = route
    env["AUTOPILOT_REPO"] = str(ROOT)
    db.execute("UPDATE tasks SET status='running',started=?,attempts=attempts+1,host=? WHERE id=?", (time.time(), host["name"], tid)); db.commit()
    try:
        subprocess.run(["bash", "-lc", cmd + " " + shlex.quote(str(task_path))], cwd=ROOT, env=env,
                       timeout=cfg["task_timeout_seconds"], check=False)
        x = validate_result(result_path)
        ingest(db, tid, result_path, x)
        return x
    except Exception as exc:
        db.execute("UPDATE tasks SET status='error',finished=? WHERE id=?", (time.time(), tid))
        db.execute("UPDATE routes SET reward=reward-0.5,attempts=attempts+1,failures=failures+1,last_change=? WHERE route=?", (time.time(), route)); db.commit()
        (run_dir / "controller_error.txt").write_text(repr(exc))
        return {"status": "error", "claim": str(exc), "implication": "none", "scope": "none"}


def report(db):
    rows = db.execute("SELECT route,attempts,successes,failures,reward FROM routes ORDER BY reward DESC").fetchall()
    done = db.execute("SELECT count(*) FROM tasks WHERE status='done'").fetchone()[0]
    queued = db.execute("SELECT count(*) FROM tasks WHERE status='queued'").fetchone()[0]
    latest = db.execute("SELECT claim,implication,scope FROM results ORDER BY created DESC LIMIT 5").fetchall()
    return {"completed": done, "queued": queued, "routes": [dict(zip(["route","attempts","successes","failures","reward"], r)) for r in rows], "latest": [dict(zip(["claim","implication","scope"], r)) for r in latest]}


def main():
    cfg = load_cfg()
    state = ROOT / cfg["state_dir"]
    state.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(state / "autopilot.sqlite", timeout=30)
    db_init(db); seed(db)
    dry = "--dry-run" in sys.argv
    while True:
        if dry:
            print(json.dumps(report(db), indent=2)); return
        active = db.execute("SELECT count(*) FROM tasks WHERE status='running'").fetchone()[0]
        if active >= cfg["max_parallel"]:
            time.sleep(cfg["poll_seconds"]); continue
        row = choose(db, ["cpu", "cuda", "rocm"])
        if row is None:
            # Refill only from novel directions supplied by successful workers; otherwise idle.
            time.sleep(cfg["poll_seconds"]); continue
        run_one(cfg, db, row)
        if db.execute("SELECT count(*) FROM results").fetchone()[0] % cfg["milestone_every_results"] == 0:
            (state / "latest_report.json").write_text(json.dumps(report(db), indent=2))


if __name__ == "__main__":
    main()
