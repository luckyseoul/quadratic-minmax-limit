#!/usr/bin/env python3
"""Cruise-control orchestrator for the quadratic-minmax research program.

The controller owns scheduling/state management. Mathematical claims are never
accepted from agent prose; they remain evidence until independently reviewed.
Set AUTOPILOT_WORKER_CMD to the existing runner/dispatcher. It receives the task
JSON path as argv[1] and must write the JSON result path from AUTOPILOT_RESULT_JSON.
"""
from __future__ import annotations
import hashlib, json, os, shlex, sqlite3, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, Future
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "autopilot" / "config.json"
ROUTES = {
 "scaling":["Find an unconditional asymptotic comparison for H(n)=m_n^(2/3) that can force convergence; prioritize composition/restriction inequalities and quantify every error.","Attack H(2n)-2H(n) and H(3n)-3H(n) from the actual minimax definition; prove a useful bound or give a rigorous obstruction."],
 "restriction_lifting":["Sharpen selected-half restriction beyond the sqrt(2) normalized loss; seek an all-orders statement, not a finite pattern.","Use retained optimal/near-optimal matrices to discover a deterministic lifting or restriction invariant, then formulate the weakest universal lemma."],
 "structure_mining":["Mine retained good matrices for invariants shared across orders and independent families; turn stable patterns into quantified lemmas.","Analyze principal restrictions, switching classes, spectra, row sums, cuts, A^2, and active Boolean maximizers; seek a universal consequence of near-optimality."],
 "spectral_rounding":["Find an unconditional spectral/Boolean-rounding inequality stronger than current lower bounds, without a limiting spectral law or optimizer-family assumption.","Stress-test actual-source Gaussian/Hermite rounding and identify the weakest missing hypothesis needed for an all-orders theorem."],
 "counterexample":["Construct two infinite families with provably separated normalized minimax values; finite numerical oscillation is insufficient.","Find a mechanism for nonconvergence despite monotonicity of m_n and vanishing adjacent alpha gaps."],
 "proof_synthesis":["Synthesize an unconditional convergence lemma from compatible route-neutral results; expose every unproved bridge.","Replace a conditional route hypothesis with a theorem derived from the original minimax definition."],
 "adversarial":["Destroy or validate a promising candidate claim: attack quantifiers, finite-to-asymptotic leakage, hidden optimality assumptions, and exceptional subsequences."],
 "exploration":["Develop a substantially different route to convergence or nonconvergence from the original definition; do not recycle a rejected obstruction under a new name."]}
VALID = {"success","counterexample","partial","dead_end","error"}

def cfg(): return json.loads(CFG.read_text())

def connect(path):
    db=sqlite3.connect(path,timeout=60); db.execute("PRAGMA journal_mode=WAL"); return db

def init(db):
    db.executescript("""
    CREATE TABLE IF NOT EXISTS tasks(id TEXT PRIMARY KEY,route TEXT,prompt TEXT,parent TEXT,status TEXT,priority REAL,attempts INTEGER DEFAULT 0,created REAL,started REAL,finished REAL,result_path TEXT,host TEXT);
    CREATE TABLE IF NOT EXISTS results(id TEXT PRIMARY KEY,task_id TEXT,status TEXT,claim TEXT,implication TEXT,scope TEXT,raw_path TEXT,digest TEXT,created REAL);
    CREATE TABLE IF NOT EXISTS routes(route TEXT PRIMARY KEY,reward REAL DEFAULT 0,attempts INTEGER DEFAULT 0,successes INTEGER DEFAULT 0,failures INTEGER DEFAULT 0,last_change REAL);
    CREATE TABLE IF NOT EXISTS meta(k TEXT PRIMARY KEY,v TEXT);
    """)

def add(db,route,prompt,parent=None,priority=1.0):
    raw=json.dumps([route,prompt,parent],sort_keys=True); tid=hashlib.sha256(raw.encode()).hexdigest()[:20]
    db.execute("INSERT OR IGNORE INTO tasks VALUES(?,?,?,?,?,?,0,?,NULL,NULL,NULL,NULL)",(tid,route,prompt,parent,"queued",priority,time.time())); return tid

def seed(db):
    for route,ps in ROUTES.items():
        db.execute("INSERT OR IGNORE INTO routes(route) VALUES(?)",(route,))
        for p in ps:add(db,route,p)
    db.commit()

def score(db,route):
    reward,attempts=db.execute("SELECT reward,attempts FROM routes WHERE route=?",(route,)).fetchone()
    return reward/max(1,attempts)+0.75/(1+attempts)**0.5

def claim(db):
    rows=db.execute("SELECT id,route,prompt,priority FROM tasks WHERE status='queued'").fetchall()
    if not rows:return None
    rows.sort(key=lambda r:r[3]*score(db,r[1]),reverse=True)
    r=rows[0]; now=time.time()
    cur=db.execute("UPDATE tasks SET status='running',started=?,attempts=attempts+1 WHERE id=? AND status='queued'",(now,r[0]))
    db.commit(); return r if cur.rowcount else None

def context(): return [f for f in ("CORE.md","STATUS.md","HANDOFF.md","ARTIFACTS.md") if (ROOT/f).exists()]

def run_task(c,r):
    tid,route,prompt,_=r; d=ROOT/c["artifact_dir"]/tid; d.mkdir(parents=True,exist_ok=True)
    tp=d/"task.json"; rp=d/"result.json"
    task={"schema":"quadratic_minmax_autopilot_task_v1","task_id":tid,"route":route,"objective":prompt,"repository":str(ROOT),"read_first":context(),"rules":["Goal: convergence or nonconvergence of alpha_n.","Finite computation is evidence, never an all-orders proof.","Do not assume L=1/2, Paley optimality, optimizer classification, or route necessity.","Reuse existing certificates; do not rerun unchanged searches.","State the exact implication established.","List assumptions for adversarial attack."],"required_result":str(rp),"time_budget_seconds":c["task_timeout_seconds"]}
    tp.write_text(json.dumps(task,indent=2)); cmd=os.environ.get("AUTOPILOT_WORKER_CMD")
    if not cmd: raise RuntimeError("AUTOPILOT_WORKER_CMD is not set")
    env=os.environ.copy(); env.update(c.get("worker_env",{})); env.update({"AUTOPILOT_TASK_JSON":str(tp),"AUTOPILOT_RESULT_JSON":str(rp),"AUTOPILOT_ROUTE":route,"AUTOPILOT_REPO":str(ROOT)})
    subprocess.run(["bash","-lc",cmd+" "+shlex.quote(str(tp))],cwd=ROOT,env=env,timeout=c["task_timeout_seconds"],check=False)
    if not rp.exists(): raise RuntimeError("worker produced no result JSON")
    x=json.loads(rp.read_text())
    if x.get("status") not in VALID or not all(k in x for k in ("claim","implication","scope")): raise RuntimeError("malformed worker result")
    db=connect(ROOT/c["state_dir"]/"autopilot.sqlite"); digest=hashlib.sha256(rp.read_bytes()).hexdigest(); rid=hashlib.sha256((tid+digest).encode()).hexdigest()[:24]
    db.execute("INSERT OR REPLACE INTO results VALUES(?,?,?,?,?,?,?,?,?)",(rid,tid,x["status"],str(x["claim"]),str(x["implication"]),str(x["scope"]),str(rp),digest,time.time()))
    positive=x["status"] in {"success","counterexample"}; reward=2 if positive else .25 if x["status"]=="partial" else -.75
    db.execute("UPDATE routes SET reward=reward+?,attempts=attempts+1,successes=successes+?,failures=failures+?,last_change=? WHERE route=?",(reward,int(positive),int(not positive),time.time(),route))
    if positive: add(db,"adversarial",f"Independently destroy or validate task {tid}. Candidate claim: {x['claim']}. Implication: {x['implication']}. Attack all assumptions and quantifiers.",tid,2.0)
    for nt in x.get("next_tasks",[])[:6]:
        if isinstance(nt,str) and nt.strip(): add(db,route,nt,tid,1.25 if positive else .9)
    db.execute("UPDATE tasks SET status='done',finished=?,result_path=? WHERE id=?",(time.time(),str(rp),tid)); db.commit(); return x

def recover(db,timeout):
    cutoff=time.time()-timeout*1.25
    db.execute("UPDATE tasks SET status='queued' WHERE status='running' AND started<?",(cutoff,)); db.commit()

def report(db):
    q=lambda s:db.execute(s).fetchone()[0]
    routes=db.execute("SELECT route,attempts,successes,failures,reward FROM routes ORDER BY reward DESC").fetchall()
    latest=db.execute("SELECT claim,implication,scope FROM results ORDER BY created DESC LIMIT 8").fetchall()
    return {"completed":q("SELECT count(*) FROM tasks WHERE status='done'"),"queued":q("SELECT count(*) FROM tasks WHERE status='queued'"),"running":q("SELECT count(*) FROM tasks WHERE status='running'"),"routes":[dict(zip(("route","attempts","successes","failures","reward"),r)) for r in routes],"latest":[dict(zip(("claim","implication","scope"),r)) for r in latest]}

def main():
    c=cfg(); state=ROOT/c["state_dir"]; state.mkdir(parents=True,exist_ok=True); db=connect(state/"autopilot.sqlite"); init(db); seed(db); recover(db,c["task_timeout_seconds"])
    if "--dry-run" in sys.argv: print(json.dumps(report(db),indent=2)); return
    workers=min(c["max_parallel"],sum(h["slots"] for h in c["hosts"]))
    futures={}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        while True:
            for f in list(futures):
                if f.done():
                    try:f.result()
                    except Exception as e: print(f"worker failure: {e}",file=sys.stderr)
                    futures.pop(f)
            while len(futures)<workers:
                r=claim(db)
                if not r: break
                futures[pool.submit(run_task,c,r)]=r[0]
            if len(futures)==0 and db.execute("SELECT count(*) FROM tasks WHERE status='queued'").fetchone()[0]==0:
                time.sleep(c["poll_seconds"])
            else: time.sleep(2)
            if db.execute("SELECT count(*) FROM results").fetchone()[0] % c["milestone_every_results"]==0 and db.execute("SELECT count(*) FROM results").fetchone()[0]>0:
                (state/"latest_report.json").write_text(json.dumps(report(db),indent=2))

if __name__=="__main__": main()
