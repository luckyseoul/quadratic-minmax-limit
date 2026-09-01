#!/usr/bin/env python3
"""k6 dashboard + JSON/POST API for per-worker soft start/stop."""
from __future__ import annotations

import argparse
import glob
import html
import json
import os
import subprocess
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from k6_control import WORKERS, stop_flag_state

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("K6_ROOT", "/mnt/storage/e1work/maxplus_p13"))
# Always the copy next to this file. K6_CODE can be stale (dash started from
# main after mesh scripts moved to mesh/k6-p13-enum) and would 500-crash POST.
CODE = HERE
OUT = ROOT / "k6_gpu_out"
RESUME = ROOT / "k6_resume.json"
MARK = ROOT / "k6_mesh_mark.json"
HTML = ROOT / "k6_dashboard.html"
PIDS = ROOT / "mesh_pids"
MESH = HERE / "k6_mesh.sh"
N_TASKS = 17805
SSH = ["ssh", "-F", "/home/nick/.ssh/mesh.config", "-o", "ConnectTimeout=2",
       "-o", "BatchMode=yes"]

HORSE = {
    "v100": "soulkiller Tesla V100 16 GiB · 2 CUDA workers · GEN_CAP 40e6",
    "nuka": "RX 9070 XT HIP · 1 worker · GEN_CAP 40e6",
    "orin": "Orin Ampere sm_87 · 1 CUDA worker · GEN_CAP 8e6 · 6c / 7.5 GiB",
    "a380": "jellyfin Arc A380 SYCL emit+test · GEN_CAP 20e6",
    "cpu44": "soulkiller CPU role (legacy label) · 11×4 threads · no CUDA",
    "dash": "HTTP UI :8765",
}
LOGS = {
    "v100": ROOT / "enum_p13_k6.log",
    "nuka": ROOT / "enum_p13_k6_nuka.log",
    "orin": ROOT / "enum_p13_k6_orin.log",
    "a380": ROOT / "enum_p13_k6_a380.log",
    "cpu44": ROOT / "enum_p13_k6_cpu44.log",
}
REMOTE = {"nuka": "nuka", "orin": "orin", "a380": "jellyfin"}


def _tail(path: Path, n=5) -> str:
    if not path.is_file():
        return "(no log)"
    try:
        data = path.read_bytes()
        if len(data) > 8000:
            data = data[-8000:]
        lines = data.decode("utf-8", "replace").replace("\r", "\n").strip().split("\n")
        return "\n".join(lines[-n:])
    except OSError as e:
        return str(e)


def _pid(name: str) -> str:
    pf = PIDS / f"{name}.pid"
    if not pf.is_file():
        return ""
    return pf.read_text().strip()


def _pgrep_local(pattern: str) -> bool:
    r = subprocess.run(
        ["pgrep", "-f", pattern],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return r.returncode == 0


def _alive(name: str) -> bool:
    if name == "v100":
        return _pgrep_local("run_kgauged.py 6 2")
    if name == "cpu44":
        pid = _pid("cpu44")
        if pid:
            try:
                os.kill(int(pid), 0)
                return True
            except (OSError, ValueError):
                pass
        return False
    if name == "dash":
        pid = _pid("dash")
        if not pid:
            return False
        try:
            os.kill(int(pid), 0)
            return True
        except (OSError, ValueError):
            return False
    pid = _pid(name)
    host = REMOTE.get(name)
    if not host:
        return False
    if pid:
        r = subprocess.run(
            SSH + [host, f"kill -0 {pid}"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if r.returncode == 0:
            return True
    r = subprocess.run(
        SSH + [host, "pgrep -f '[r]un_kgauged.py'"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return r.returncode == 0


def snapshot():
    npy = glob.glob(str(OUT / "orb*.npy"))
    n_done = len(npy)
    n_left = max(N_TASKS - n_done, 0)
    rec = {}
    if RESUME.is_file():
        try:
            rec = json.loads(RESUME.read_text())
        except json.JSONDecodeError:
            rec = {}
    mark = {}
    if MARK.is_file():
        try:
            mark = json.loads(MARK.read_text())
        except json.JSONDecodeError:
            mark = {}
    now = time.time()
    t0 = float(mark.get("t0", now))
    n0 = int(mark.get("n_done", n_done))
    hours = max((now - t0) / 3600.0, 1e-6)
    gained = max(n_done - n0, 0)
    rate = gained / hours
    eta_h = (n_left / rate) if (rate > 0.5 and gained >= 8 and hours >= 0.15) else None
    pct = 100.0 * n_done / N_TASKS
    flags = stop_flag_state()
    nodes = {}
    for w in WORKERS:
        running = _alive(w)
        stopping = bool(flags.get(w)) and running
        nodes[w] = {
            "name": w,
            "horse": HORSE[w],
            "running": running,
            "stopping": stopping,
            "pid": _pid(w),
            "log": _tail(LOGS[w]) if w in LOGS else "",
            "flag": bool(flags.get(w)),
        }
    inflight = 0
    starve = 0
    for d in glob.glob(str(OUT / ".lock_*")):
        if os.path.isfile(os.path.join(d, "STARVE")):
            starve += 1
        else:
            inflight += 1
    return dict(
        n_done=n_done, n_left=n_left, n_tasks=N_TASKS, pct=pct,
        n_solutions=rec.get("n_solutions"),
        rate=rate, gained=gained, hours=hours, eta_h=eta_h,
        updated=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        nodes=nodes, locks=inflight, starve_locks=starve,
        all_flag=bool(flags.get("ALL")),
    )


def render(s):
    if s["eta_h"] is None:
        eta = "warming (need a few new orbs)"
    elif s["eta_h"] >= 48:
        eta = f"{s['eta_h']/24:.1f} days"
    else:
        eta = f"{s['eta_h']:.1f} hours"
    rate = f"{s['rate']:.1f}/h" if s["rate"] > 0 else "—"
    cards = []
    for w in WORKERS:
        n = s["nodes"][w]
        if n["stopping"]:
            st = "SOFT-STOPPING (finish current orbit)"
            cls = "stop"
        elif n["running"]:
            st = f"RUNNING pid {n['pid']}"
            cls = "on"
        else:
            st = "stopped"
            cls = "off"
        log = (
            f"<pre>{html.escape(n['log'])}</pre>" if n["log"] else ""
        )
        cards.append(
            f"<section class='card {cls}' data-w='{w}'>"
            f"<h2>{w}</h2><p class=hp>{html.escape(n['horse'])}</p>"
            f"<p class=st>{html.escape(st)}</p>"
            f"<p class=btns>"
            f"<button type=button onclick=\"act('start','{w}')\">Start</button> "
            f"<button type=button onclick=\"act('stop','{w}')\">Soft stop</button>"
            f"</p>{log}</section>"
        )
    page = f"""<!doctype html>
<html lang=en>
<meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1,viewport-fit=cover">
<meta name=theme-color content="#111111">
<title>k=6 p=13 mesh</title>
<style>
:root{{color-scheme:dark}}
*{{box-sizing:border-box}}
body{{font:16px/1.45 ui-monospace,monospace;background:#111;color:#ddd;
 margin:0 auto;padding:1rem;padding:max(1rem,env(safe-area-inset-top))
 max(1rem,env(safe-area-inset-right)) max(1rem,env(safe-area-inset-bottom))
 max(1rem,env(safe-area-inset-left));max-width:64rem}}
h1{{font-size:1.1rem;color:#fff;margin:.2rem 0 .6rem}}
.big{{font-size:clamp(1.4rem,6vw,2.1rem);color:#9f6;margin:.3rem 0;overflow-wrap:anywhere}}
.bar{{background:#333;height:12px;border-radius:6px;overflow:hidden}}
.bar>i{{display:block;height:100%;width:{s['pct']:.2f}%;background:#6c6}}
.stats{{display:flex;flex-wrap:wrap;gap:.4rem .9rem;margin:.7rem 0}}
.stats span{{background:#1c1c1c;border:1px solid #333;border-radius:6px;padding:.25rem .55rem}}
pre{{background:#1c1c1c;padding:10px;overflow:auto;font-size:11px;color:#bbb;max-height:9rem;
 -webkit-overflow-scrolling:touch;white-space:pre-wrap;word-break:break-word}}
.k{{color:#888;font-size:13px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}}
.card{{background:#1a1a1a;padding:12px 14px;border-radius:8px;border:1px solid #333}}
.card.on{{border-color:#3a6}}
.card.stop{{border-color:#a83}}
.card.off{{border-color:#444}}
.hp{{color:#9ab;font-size:12px;overflow-wrap:anywhere}}
.st{{margin:.3rem 0}}
.btns{{display:flex;flex-wrap:wrap;gap:8px;margin:.5rem 0}}
button{{font:15px ui-monospace,monospace;padding:.6rem 1rem;cursor:pointer;
 background:#222;color:#eee;border:1px solid #555;border-radius:6px;
 min-height:44px;touch-action:manipulation;-webkit-tap-highlight-color:transparent}}
button:hover,button:focus{{background:#333}}
.row{{display:flex;flex-wrap:wrap;gap:8px;margin:.6rem 0 1rem}}
.msg{{color:#fc6;min-height:1.2em}}
@media (max-width:600px){{
  body{{padding:.75rem;padding-bottom:max(.75rem,env(safe-area-inset-bottom))}}
  .grid{{grid-template-columns:1fr}}
  .btns,.row{{flex-direction:column}}
  button{{width:100%}}
  pre{{font-size:12px;max-height:8rem}}
}}
</style>
<h1>p=13 k=6 mesh</h1>
<div class=big>{s['n_done']} / {s['n_tasks']} &nbsp; {s['pct']:.2f}%</div>
<div class=bar><i></i></div>
<div class=stats>
<span>left <b>{s['n_left']}</b></span>
<span>in-flight <b>{s['locks']}</b></span>
<span>starve <b>{s.get('starve_locks', 0)}</b></span>
<span>rate <b>{html.escape(rate)}</b></span>
<span>+{s['gained']} in {s['hours']:.2f} h</span>
<span>ETA <b>{html.escape(eta)}</b></span>
</div>
<p class=k>{html.escape(s['updated'])} · One card never stops the rest.
Soft stop finishes the current orbit.</p>
<p class=msg id=msg></p>
<p class=row>
<button type=button onclick="act('start','all')">Start all</button>
<button type=button onclick="act('stop','all')">Soft stop all</button>
</p>
<div class=grid>
{''.join(cards)}
</div>
<script>
let busy = false;
async function act(op, name) {{
  const m = document.getElementById('msg');
  if (busy) {{ m.textContent = 'wait — request in flight'; return; }}
  busy = true;
  m.textContent = op + ' ' + name + ' …';
  try {{
    const r = await fetch('/api/' + op + '/' + name, {{method: 'POST'}});
    const t = await r.text();
    m.textContent = r.ok ? t : (r.status + ' ' + t);
  }} catch (e) {{
    m.textContent = String(e);
  }} finally {{
    busy = false;
  }}
}}
setInterval(function() {{ if (!busy) location.reload(); }}, 15000);
</script>
"""
    return page, eta, rate


def write_html():
    s = snapshot()
    html, eta, rate = render(s)
    tmp = HTML.with_name(f"k6_dashboard.{os.getpid()}.{time.time_ns()}.tmp")
    tmp.write_text(html)
    tmp.replace(HTML)
    line = (
        f"k6 {s['n_done']}/{s['n_tasks']} {s['pct']:.2f}% left={s['n_left']} "
        f"rate={rate} eta={eta} locks={s['locks']}"
    )
    for w in WORKERS:
        n = s["nodes"][w]
        bit = "R" if n["running"] else "."
        if n["stopping"]:
            bit = "S"
        line += f" {w}:{bit}"
    print(line, flush=True)
    return s


def _spawn_mesh(args: list[str]) -> str:
    if not MESH.is_file():
        raise FileNotFoundError(f"mesh script missing: {MESH}")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(HERE) + os.pathsep + env.get("PYTHONPATH", "")
    env["K6_ROOT"] = str(ROOT)
    env["K6_CODE"] = str(HERE)
    log = ROOT / "k6_mesh_cli.log"
    with log.open("a") as fh:
        fh.write(f"\n===== {' '.join(args)} {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
        fh.flush()
        subprocess.Popen(
            [str(MESH), *args],
            cwd=str(HERE),
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    return f"queued {' '.join(args)} (soft stop finishes current orbit)"


class Handler(SimpleHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def do_GET(self):
        u = urlparse(self.path)
        if u.path in ("/", "/k6_dashboard.html"):
            write_html()
            self.path = "/k6_dashboard.html"
        if u.path == "/api/status":
            body = json.dumps(snapshot()).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        return SimpleHTTPRequestHandler.do_GET(self)

    def _read_body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n:
            self.rfile.read(n)

    def _text(self, code: int, msg: str):
        body = msg.encode()
        self.send_response(code)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        u = urlparse(self.path)
        parts = u.path.strip("/").split("/")
        ok_names = set(WORKERS) | {"all"}
        self._read_body()
        if len(parts) == 3 and parts[0] == "api" and parts[1] in ("start", "stop") and parts[2] in ok_names:
            op, name = parts[1], parts[2]
            try:
                msg = _spawn_mesh([op, name])
            except Exception as e:
                self._text(500, f"{type(e).__name__}: {e}")
                return
            self._text(202, msg)
            return
        self.send_error(404)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--serve", type=int, default=0)
    ap.add_argument("--loop", type=float, default=8.0)
    args = ap.parse_args()
    if args.serve:
        os.chdir(ROOT)

        httpd = ThreadingHTTPServer(("0.0.0.0", args.serve), Handler)

        def loop():
            while True:
                try:
                    write_html()
                except Exception as e:
                    print("dash", e, flush=True)
                time.sleep(args.loop)

        import threading
        threading.Thread(target=loop, daemon=True).start()
        print(f"dashboard http://0.0.0.0:{args.serve}/k6_dashboard.html", flush=True)
        httpd.serve_forever()
    else:
        write_html()


if __name__ == "__main__":
    main()
